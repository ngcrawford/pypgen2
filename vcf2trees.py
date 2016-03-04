#!/usr/bin/env python
# encoding: utf-8
"""
vcf2phylip.py

Created by Nick Crawford on 2011-11-18.
Copyright (c) 2011

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses

The author may be contacted at ngcrawford@gmail.com
"""

import os
import re
import sys
import gzip
import glob
import pysam
import shlex
import random
import argparse
import tempfile
import subprocess
import itertools
import numpy as np
from subprocess import Popen, PIPE
from collections import namedtuple
from ngs_parsers import VCF


def get_args():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser()

    parser.add_argument('-n', '-N', '--cores',
                        required=True,
                        type=int,
                        help="Number of processor cores to use.")

    parser.add_argument('-r', '-R', '--region-file',
                        required=True,
                        type=argparse.FileType('r'),
                        help="File with regions to make trees from in \
                              this format 'Chrm:start-stop'")

    parser.add_argument('-s', '-S', '--samples-2-keep',
                        required=True,
                        type=argparse.FileType('r'),
                        help="File with sample IDs to retain. One line \
                              per sample.")

    parser.add_argument('-o', '--output',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='Path to output. (default is STOUT)')

    parser.add_argument('input',
                        nargs=1,
                        help='bgzipped and indexed VCF file')

    args = parser.parse_args()
    return args


def makeDataTuple(vcf):
    """Setup a labeled tuple to store the data."""

    chrm_data = {}

    for count, line in enumerate(vcf):

        if "##contig" in line:
            contigline = line.split("<")[1]
            contigline = contigline.strip(">\n").split(",")
            contigline = [item.split("=")[1] for item in contigline]
            chrm_data[contigline[0]] = int(contigline[1])

        if line.startswith("#CHROM"):
            field_labels = line.strip().strip("#").split("\t")
            field_labels = [item.strip().split('.')[0].replace("-", "_") for item in field_labels]
            break

    position_data = namedtuple('base', field_labels)
    return (position_data, chrm_data)


def array2OnelinerAlignment(info, taxa, bases):
    """Convert array of array of taxa and an array of bases to one-liner."""

    oneliner = info
    for count, seq in enumerate(bases):
        oneliner += taxa[count] + "," + ''.join(itertools.chain(bases[count])) + ","
    oneliner = oneliner[:-1] + ";"
    return oneliner


def process_snp_call(snp_call, ref, alt, IUPAC_ambiguities=False):
    """Process VCF genotype fields.
        The current version is very basic and
        doesn't directly take into account the
        quality of the call or call hets with
        IUPAC ambiguity codes."""

    # IUPAC ambiguity codes
    IUPAC_dict = {('A', 'C'): 'M',
                  ('A', 'G'): 'R',
                  ('A', 'T'): 'W',
                  ('C', 'G'): 'S',
                  ('C', 'T'): 'Y',
                  ('G', 'T'): 'K',
                  ('A', 'C', 'G'): 'V',
                  ('A', 'C', 'T'): 'H',
                  ('A', 'G', 'T'): 'D',
                  ('C', 'G', 'T'): 'B'}

    #called_base = ""
    snp_call = snp_call.split(":")

    # process blanks
    if snp_call[0] == "./.":
        called_base = "-"

    else:
        allele1, allele2 = snp_call[0].split("/")

        # process "0/0"
        if allele1 == '0' and allele2 == '0':
            called_base = ref

        if allele1 == '1' and allele2 == '1':
            called_base = alt

        # process "0/N"
        if allele1 == '0' and allele2 != '0':

            if IUPAC_ambiguities == False:
                called_base = 'N'

            else:
                call = [ref] + [alt.split(',')[int(allele2) - 1]]
                call.sort()
                call = tuple(call)
                called_base = IUPAC_dict[call]

        # process "2/2, 1/2, etc."
        if int(allele1) >= 1 and int(allele2) > 1:

            # deal with homozygotes
            if allele1 == allele2:
                called_base = alt.split(',')[int(allele1) - 1]

            # deal with heterozygotes
            else:

                if IUPAC_ambiguities == False:
                    called_base = 'N'

                else:
                    ref = alt.split(',')[int(allele1) - 1]
                    alt = alt.split(',')[int(allele2) - 1]
                    call = [ref, alt]
                    call.sort()
                    call = tuple(call)
                    called_base = IUPAC_dict[call]

    return called_base


def callSNPs(current_base, numb_of_seqs):
    """Call the SNPs. Duh!"""

    blanks =  np.zeros(numb_of_seqs, np.string0)

    if current_base.FILTER == 'LowQual':
        blanks.fill("-")

    if current_base.FORMAT == 'GT':
        blanks.fill("-")

    for count, snp_call in enumerate(current_base[9:]):

        base = process_snp_call(snp_call, current_base.REF, current_base.ALT)
        blanks[count] = base

    return blanks


def count_informative_sites(alignment_array):
    """Informative Sites must have two different SNPs"""
    informative_sites = 0
    for site in alignment_array:
        unique_sites = set(site)
        if len(unique_sites) >= 3:
            informative_sites += 1
    return informative_sites




def parse_window_vcf(vcf, start, stop, window_size, chrm, fout):
    # SETUP NAMED TUPLE TO STORE INFO FROM A SINGLE BASE
    field_labels = []
    position_data, chrm_data = makeDataTuple(vcf)

    # SETUP MULTIPLE ALIGNMENT ARRAY
    numb_of_seqs = len(position_data._fields[9:])
    alignment = np.zeros((window_size,numb_of_seqs), np.string0)

    # SETUP COUNTERS
    current_base = None
    current_window = 1
    line_count = 0
    windows = range(0, chrm_data[chrm], window_size)
    current_data = []
    informative_sites = []

    # PARSE VCF FIlE
    snp_count = 0

    for line in vcf:
        # SKIP HEADER
        if line.startswith("#CHROM"):
            line_count = 0

        # START PROCESSING ALIGNED BASES
        if line.startswith(chrm):
            current_base = position_data._make(line.strip().split("\t"))
            base_calls = callSNPs(current_base, numb_of_seqs)
            current_data.append(base_calls.copy())

    alignment = np.array(current_data)
    inform_sites = count_informative_sites(alignment)
    if current_base == None:
        return 'error'
    else:
        taxa = current_base._fields[9:]
        info = 'chrm={0},start={1},stop={2},inform_sites={3}'.format(current_base.CHROM, start, stop, inform_sites)
        oneliner = array2OnelinerAlignment(info, taxa, alignment.transpose())

    if ":" in oneliner and oneliner[-1] == ';': # this prevents bad alignments from getting printed
        return oneliner
    else:
        return 'error'

def header_slices(vcf, window_size=5000):
    vcf = pysam.Tabixfile(vcf)

    slices = {}
    for line in vcf.header:

        if line.startswith('##contig'):

            line_dict = dict([item.split("=") for item in line.split("<")[1].strip('>').split(",")])
            length = int(line_dict['length'])

            if length < window_size: continue

            start = (length % window_size)/2
            stop = ((length/window_size) * window_size) + start

            s = xrange(start,stop,window_size)
            slices[line_dict['ID']] = s

    return slices


def slice_vcf(vcf, chrm, start, stop):
    cli = """tabix {0} {1}:{2}-{3}""".format(vcf, chrm, start, stop)
    cli_parts = shlex.split(cli)
    vcf = Popen(cli_parts, stdin=PIPE, stderr=PIPE, stdout=PIPE).communicate()[0]
    return vcf

def process_vcf_slice(tabix_file, chrm, start, stop, position_data):

    tbx = pysam.Tabixfile(tabix_file)
    tbx_lines = tbx.fetch(chrm, start, stop)

    numb_of_seqs = len(position_data._fields[9:])
    alignment = np.zeros((stop-start, numb_of_seqs), np.string0)

    # This 'error handling' needs to be rewritten.
    current_data = []
    if tbx_lines == None:
        return 'error'

    for line in tbx_lines:
        current_base = position_data._make(line.strip().split("\t"))
        base_calls = callSNPs(current_base, numb_of_seqs)
        current_data.append(base_calls.copy())

    alignment = np.array(current_data)
    inform_sites = count_informative_sites(alignment)

    if current_base == None:
        return 'error'
    else:
        taxa = current_base._fields[9:]
        info = "tree 'chrm={0},start={1},stop={2},inform_sites={3}':".format(current_base.CHROM, start, stop, inform_sites)
        oneliner = array2OnelinerAlignment(info, taxa, alignment.transpose())

    if ":" in oneliner and oneliner[-1] == ';':  # this prevents bad alignments from getting printed
        return oneliner
    else:
        return 'error'

def oneliner2phylip(line):
    """Convert one-liner to phylip format."""
    seqs = line.strip(";").split(':')[-1].split(',')
    label_seqs = zip(seqs[:-1:2],seqs[1::2])
    taxa_count = len(label_seqs)
    seq_length = len(label_seqs[0][1])
    alignment = "%s %s\n" % (taxa_count, seq_length) # add header
    for taxa_name, seq in label_seqs:
        taxa_name = taxa_name.strip()
        alignment += '%-10s%s\n' % (taxa_name, seq)
    return alignment


def oneliner2phylip(line, samples_2_keep):
    """Convert one-liner to phylip format."""

    seqs = line.strip(";").split(':')[-1].split(',')

    sample_ids = seqs[:-1:2]

    seqs = seqs[1::2]

    label_seqs = zip(sample_ids, seqs)

    final_seqs = []
    for i in label_seqs:
        if i[0] in samples_2_keep:
            final_seqs.append(i)

    label_seqs = final_seqs

    taxa_count = len(label_seqs)
    seq_length = len(label_seqs[0][1])

    alignment = "%s %s\n" % (taxa_count, seq_length)  # add header
    for taxa_name, seq in label_seqs:
        taxa_name = taxa_name.strip()
        alignment += '%-10s%s\n' % (taxa_name, seq)
    return alignment




def process_region(chrm, start, stop, position_data, samples_2_keep, args):

    start, stop = int(start), int(stop)



    oneliner = process_vcf_slice(args.input[0], chrm,
                                 start, stop, position_data)

    phylip = oneliner2phylip(oneliner, samples_2_keep)

    #tempfile.tempdir = '/Users/testudines/Code/pypgen2'
    #phylip_file = tempfile.NamedTemporaryFile()
    prefix = '{}_{}_{}'.format(chrm, start, stop)
    phylip_file = open('/Users/testudines/Code/pypgen2/{}.phylip'.format(chrm), 'w')
    phylip_file.write(phylip)
    phylip_file.close()

    cli = "raxmlHPC-PTHREADS-SSE3 \
        -T {} \
        -s {} \
        -m GTRGAMMA \
        -p 12345 \
        -n {}".format(args.cores, phylip_file.name, prefix)

    raxml_args = shlex.split(cli)
    sbp = subprocess.Popen(raxml_args,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    sbp = sbp.communicate()


    # read in bestTree
    tree = open('RAxML_bestTree.{}'.format(prefix)).readline()

    for i in glob.glob('*.{}'.format(prefix)):
        os.remove(i)


    os.remove('RAxML_flagCheck')
    os.remove('Hmel221018.phylip')
    os.remove('Hmel221018.phylip.reduced')

    tree_id =  oneliner.strip(";").split(':')[0]
    tree = '{} {}'.format(tree_id, tree)
    return tree

def main():


    # SETUP ARGS
    args = get_args()

    # OPEN VCF
    vcf_file = gzip.open(args.input[0],'rb')
    position_data, chrm_data = makeDataTuple(vcf_file)
    vcf_file.close()

    samples_2_keep = [i.strip() for i in args.samples_2_keep]

    fout = open(args.region_file.name.strip('.txt')+'_out.trees', 'w')

    for c, r in enumerate(args.region_file):

        if c > 3:
            break

        r = r.strip()

        chrm, start, stop = re.split(r':|-', r)

        #tree =  process_region(chrm, start, stop, position_data, samples_2_keep, args)

        try:
            tree = process_region(chrm, start, stop, position_data, samples_2_keep, args)

        except:
            tree = ('bad_tree = {}:{}-{}\n'.format(chrm, start, stop))

        fout.write(tree)

if __name__ == '__main__':
    main()







