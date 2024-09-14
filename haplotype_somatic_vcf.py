import argparse
import os
import sys
import pysam
from cyvcf2 import VCF, Writer

# this script was adapted with reference to:
# https://github.com/shlokanegi/denovo_smallvars/blob/master/docker/validate_denovos/dnvval.py

def get_AH(chrom, pos, alt):
    print(chrom, pos, alt)
    # only consider SNVs, filter out indel cases
    if len(alt) == 1:
        for pileupcolumn in pysamfile.pileup(chrom, pos - 1, pos, truncate=True, min_base_quality=7):
            base_dict = {
                'A': [],
                'C': [],
                'G': [],
                'T': []
            }
            unphased_read_count = 0
            indel_read_count = 0
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    try:
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        HP = pileupread.alignment.get_tag("HP")
                        base_dict[base] += [HP]
                    # keep track of number of unphased reads
                    except KeyError:
                        unphased_read_count += 1
                else:
                    indel_read_count += 1

                # if next position contains an indel, keep track of it to filter out this position
                if pileupread.indel != 0:
                    indel_read_count += 1

            # number of haplotypes for the alternate variant
            AH = len(set(base_dict[alt]))

            # exclude cases with 3+ variants
            number_variant_type = 0
            A = len(base_dict['A'])
            C = len(base_dict['C'])
            G = len(base_dict['G'])
            T = len(base_dict['T'])
            if A > 0 :
                number_variant_type += 1
            if C > 0:
                number_variant_type += 1
            if G > 0:
                number_variant_type += 1
            if T > 0:
                number_variant_type += 1
            # if only 2 types of nucleotide bases
            if number_variant_type <= 2:
                # filter out positions that have too many indels or unphased reads
                if indel_read_count < indel_read_threshold and unphased_read_count < unphased_read_threshold:
                        return AH


parser = argparse.ArgumentParser()
parser.add_argument('-bam', help='BAM file (indexed)', required=True)
parser.add_argument('-v', help='VCF file (can be bgzipped)', required=True)
parser.add_argument('-u', type=int, help='unphased read threshold (filter out positions with unphased reads above this threshold', default=3)
parser.add_argument('-i', type=int, help='indel read threshold (filter out positions with indel reads above this threshold', default=3)
parser.add_argument('-o', default='out.vcf', help='output (annotated) VCF (will be bgzipped if ending in .gz)')
parser.add_argument('-t', default=2, help='number of threads used by bcftools mpileup')
args = parser.parse_args()

phased_bam = args.bam
unphased_read_threshold = args.u
indel_read_threshold = args.i
pysamfile = pysam.AlignmentFile(phased_bam, "rb")

vcf = VCF(args.v)
vcf.add_info_to_header({'ID': 'AH',
                        'Description': 'Alternate Haplotyping: Indicates if alternate variants are present on 1 or 2 haplotypes',
                        'Type': 'Integer', 'Number': '1'})
vcf_o = Writer(args.o, vcf)

# Read VCF and get haplotype information at each position in VCF
for variant in vcf:
    chrom = variant.CHROM
    pos = variant.POS
    alt = variant.ALT[0]
    ref = variant.REF

    # get_AH() gets haplotype information and filters out low quality positions
    AH = get_AH(chrom, pos, alt)
    if AH == None:
        continue
    else:
        variant.INFO['AH'] = AH
        vcf_o.write_record(variant)

vcf_o.close()
vcf.close()