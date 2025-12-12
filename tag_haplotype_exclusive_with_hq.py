import argparse
import pysam
import pysam.samtools
from cyvcf2 import VCF, Writer
from typing import Optional

parser = argparse.ArgumentParser()
parser.add_argument('-bam', help='BAM file (indexed and phased)', required=False)
parser.add_argument('-v', help='somatic VCF file (can be bgzipped)', required=True)
parser.add_argument('-g', help='germline VCF file (can be bgzipped)', required=True)
parser.add_argument('-u', type=int, help='unphased read threshold (filter out positions with unphased reads above this threshold', default=3)
parser.add_argument('-i', type=int, help='indel read threshold (filter out positions with indel reads above this threshold', default=3)
parser.add_argument('-k', type=int, help='number of germline variants on a read, set as limit for haplotype quality filtering)', default=2)
parser.add_argument('-o', default='out.vcf', help='output (annotated) VCF (will be bgzipped if ending in .gz)')


# optional args
parser.add_argument('-hqtaggedbam', help='haplotag quality tagged BAM file (indexed)', required=False)



args = parser.parse_args()

phased_bam = args.bam
hqtagged_bam = args.hqtaggedbam if args.hqtaggedbam else None
unphased_read_threshold = args.u
indel_read_threshold = args.i
pysamfile = pysam.AlignmentFile(phased_bam, "rb")
pysamfile_outfile = pysam.AlignmentFile(phased_bam.replace('.bam', '.haplotype_quality_tagged.bam'), "wb", template=pysamfile)

vcf = VCF(args.v)
vcf.add_info_to_header({'ID': 'AH',
                        'Description': 'Alternate Haplotyping: Indicates if alternate variants are present on 1 or 2 haplotypes',
                        'Type': 'Integer', 'Number': '1'})
# vcf.add_info_to_header({'ID': 'HQ',
#                         'Description': 'Haplotype Quality: Indicates whether the read have sufficient germline variants in agreement to be considered accurately phased. 1 = sufficient, 0 = insufficient',
#                         'Type': 'Integer', 'Number': '1'})
vcf_o = Writer(args.o, vcf)

germline_vcf = VCF(args.g)


def get_query_read_base(chrom, pos, query_name):

    # get pileup of reads aligned to this position "pos"
    for pileup_column in pysamfile.pileup(chrom, pos-1, pos, truncate=True):

        # check if this is the position we are interested in
        if pileup_column.reference_pos != pos-1:
            continue
            
        #print(f"Reference {chrom}:{pos}")
        for pileup_read in pileup_column.pileups:
            read = pileup_read.alignment

            if pileup_read.is_del:
                #print(read.query_name, "→ deletion (no base)")
                continue
            elif pileup_read.is_refskip:
                #print(read.query_name, "→ reference skip (N)")
                continue
            else:
                base = read.query_sequence[pileup_read.query_position]
                if read.query_name == query_name:
                    #print(read.query_name, "→", base)
                    return base
# source: https://bioinformatics.stackexchange.com/questions/7401/access-base-aligned-to-particular-reference-position
def find_base_in_alignment(alignment: pysam.AlignedSegment,
                           pos: int, 
                           bam_stores_revcomp: bool = False)-> Optional[str]:
    idx_q = 0
    idx_r = pos - alignment.reference_start
    if bam_stores_revcomp:
        seq = alignment.query_sequence
    else:
        seq = alignment.get_forward_sequence()
    
    if seq is None:
        return None
    
    for op, l in alignment.cigartuples:
        ref_consumed = op in {0, 2, 3, 7, 8}
        query_consumed = op in {0, 1, 4, 7, 8}
        
        if ref_consumed:
            idx_r -= l
        if query_consumed:
            idx_q += l
        
        if idx_r < 0:
            if query_consumed:
                # base is in query between idx_q-l , idx_q
                base = seq[idx_q + idx_r - 1]
                return base
            else:
                # position has been deleted
                return None
               
def tag_haplotype_quality():
    read_count = 0
    unphased_read_count = 0
    phased_read_count = 0

    lq_read_count = 0
    hq_read_count = 0
    indel_read_count = 0
    hp_1_count = 0
    hp_2_count = 0
    unexpected_hap_count = 0
    # go through all reads in BAM file
    for alignedsegment in pysamfile.fetch():
        read_count += 1
        print('beginning for loop for read number:', read_count)

        read_start_pos = alignedsegment.reference_start
        read_end_pos = alignedsegment.reference_end - 1

        #read_end_pos = alignedsegment.reference_start + alignedsegment.query_length - 1
        chrom = pysamfile.get_reference_name(alignedsegment.reference_id)
        region = chrom + ":" + str(read_start_pos) + "-" + str(read_end_pos)
        haplotype = alignedsegment.get_tag("HP") if alignedsegment.has_tag("HP") else None

        if haplotype == None:
            #print("read is unphased, skipping...")
            pysamfile_outfile.write(alignedsegment)
            unphased_read_count += 1
            continue
        else:
            phased_read_count += 1
            germline_variants_in_read = [x for x in germline_vcf(region)]

            hp1_gv = []
            hp2_gv = []
            for gv in germline_variants_in_read:
                #print("germline var:", gv.CHROM, gv.POS, "REF:", gv.REF, "ALT:", gv.ALT, "GT:", gv.genotypes)

                # check if it is a phased variant
                if gv.genotypes[0][2] == True:
                    # collect lists of germline variants on each haplotype
                    # add to corresponding haplotype list
                    if gv.genotypes[0][0] == 1:
                        hp1_gv.append(gv)
                    elif gv.genotypes[0][1] == 1:
                        hp2_gv.append(gv)

            if haplotype == 1:
                hp_1_count += 1
                agreeing_gv_count = 0
                non_agreeing_gv_count = 0
                for gv in hp1_gv:
                    ref_base = gv.REF
                    alt_base = gv.ALT[0]

                    #print("HP=1 read, GV at pos", gv.POS, "REF:", ref_base, "ALT:", alt_base, "read base:", read_base)
                    if len(ref_base) > 1 or len(alt_base) > 1:
                        indel_read_count += 1
                        continue
                    else:
                        read_base = find_base_in_alignment(alignedsegment, gv.POS, bam_stores_revcomp=True)
                        if read_base == alt_base:
                            #print("agrees:", gv.POS, " ref:", ref_base, " alt:", alt_base, " read_base:", read_base)
                            agreeing_gv_count += 1
                        else:
                            #print("DISAGREES:", "name", alignedsegment.query_name, gv.POS, " ref:", ref_base, " alt:", alt_base, " read_base:", read_base)
                            non_agreeing_gv_count += 1

            elif haplotype == 2:
                hp_2_count += 1
                agreeing_gv_count = 0
                non_agreeing_gv_count = 0
                for gv in hp2_gv:
                    ref_base = gv.REF
                    alt_base = gv.ALT[0]

                    if len(ref_base) > 1 or len(alt_base) > 1:
                        indel_read_count += 1
                        continue
                    else:
                        read_base = find_base_in_alignment(alignedsegment, gv.POS, bam_stores_revcomp=True)
                        if read_base == alt_base:
                            #print("agrees:", gv.POS, " ref:", ref_base, " alt:", alt_base, " read_base:", read_base)
                            agreeing_gv_count += 1
                        else:
                            #print("DISAGREES:", "name", alignedsegment.query_name, gv.POS, " ref:", ref_base, " alt:", alt_base, " read_base:", read_base)
                            non_agreeing_gv_count += 1

            print("agreeing gv count:", agreeing_gv_count)
            print("non-agreeing gv count:", non_agreeing_gv_count)

            if agreeing_gv_count >= args.k and agreeing_gv_count > non_agreeing_gv_count:
                print("Haplotype Quality: sufficient (1)")
                hq_read_count += 1
                alignedsegment.set_tag("HQ", 1)
                pysamfile_outfile.write(alignedsegment)
            else:
                print("Haplotype Quality: insufficient (0)")
                lq_read_count += 1
                alignedsegment.set_tag("HQ", 0)
                pysamfile_outfile.write(alignedsegment)

    pysamfile.close()
    pysamfile_outfile.close()
    print("Total number of reads processed for haplotype quality tagging:", read_count)

    print("Number of unphased reads:", unphased_read_count)
    print("Number of phased reads:", phased_read_count)
    print("Number of high-quality phased reads:", hq_read_count)
    print("Number of low-quality phased reads:", lq_read_count)
    print("Number of HP=1 reads:", hp_1_count)
    print("Number of HP=2 reads:", hp_2_count)
    print("Number of unexpected haplotype reads:", unexpected_hap_count)
    print("Indel read count:", indel_read_count)

    

# function to calculate AH value, only accounting for alt alleles exclusive to 1 haplotype
def get_AH(chrom, pos, alt):
    # only consider SNVs, filter out indel cases
    if len(alt) == 1:
        # use the outfile pysamfile to consider haplotype quality tag

        #for pileupcolumn in pysamfile_outfile.pileup(chrom, pos - 1, pos, truncate=True, min_base_quality=7):
        for pileupcolumn in hqtagged_bam_temp.pileup(chrom, pos - 1, pos, truncate=True, min_base_quality=7):
            base_dict = {
                'A': [],
                'C': [],
                'G': [],
                'T': []
            }
            unphased_read_count = 0
            indel_read_count = 0
            for pileupread in pileupcolumn.pileups:
                print("haplotype quality tag (HQ):", pileupread.alignment.get_tag("HQ") if pileupread.alignment.has_tag("HQ") else "not present")
                if not pileupread.is_del and not pileupread.is_refskip and not pileupread.alignment.get_tag("HQ") == 0:
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

# function to calculate AH value, accounting for number of ref and alt alleles on the alt haplotype
def assign_AH(chrom, pos, alt, ref):
    
    # only consider SNVs, filter out indel cases
    if len(alt) == 1:
        for pileupcolumn in hqtagged_bam_temp.pileup(chrom, pos - 1, pos, truncate=True, min_base_quality=2):
            haplotype_dict = {
                '1': [],
                '2': []
            }
            unphased_read_count = 0
            print("pileup: ", pileupcolumn.get_num_aligned())
            for pileupread in pileupcolumn.pileups:
                print("haplotype quality tag (HQ):", pileupread.alignment.get_tag("HQ") if pileupread.alignment.has_tag("HQ") else "not present")
                print("HP:", pileupread.alignment.get_tag("HP") if pileupread.alignment.has_tag("HP") else "not present")
                if not pileupread.is_del and not pileupread.is_refskip:
                    try:
                        # skip reads with insufficient haplotype quality
                        if pileupread.alignment.get_tag("HQ") == 1:
                            base = pileupread.alignment.query_sequence[pileupread.query_position]
                            HP = str(pileupread.alignment.get_tag("HP"))
                            #print("HP:", HP)
                            haplotype_dict[HP] += [base]
                    # keep track of number of unphased reads
                    except KeyError:
                        unphased_read_count += 1

            hap1_ref_count = haplotype_dict['1'].count(ref)
            hap1_alt_count = haplotype_dict['1'].count(alt)
            hap2_ref_count = haplotype_dict['2'].count(ref)
            hap2_alt_count = haplotype_dict['2'].count(alt)

            print("hap1_ref_count:", hap1_ref_count, "hap1_alt_count:", hap1_alt_count)
            print("hap2_ref_count:", hap2_ref_count, "hap2_alt_count:", hap2_alt_count)
            

            # conservatively only select positions where one haplotype is 100% ref
            if hap1_alt_count == 0:
                # check if hap2 contains ref and alt alleles above threshold
                if hap2_ref_count > 2 and hap2_alt_count > 2:
                    AH = 1
                else:
                    AH = 2

            elif hap2_alt_count == 0:
                # check if hap1 contains ref and alt alleles above threshold
                if hap1_ref_count > 2 and hap1_alt_count > 2:
                    AH = 1
                else:
                    AH = 2

            else:
                AH = 2
            
            print("AH:", AH)
            print("unphased_read_count:", unphased_read_count)
            print("------------------------------")
            if unphased_read_count < unphased_read_threshold:
                return AH


#hqtagged_bam_temp = pysam.AlignmentFile(hqtagged_bam, "rb") if hqtagged_bam else pysam.AlignmentFile(phased_bam.replace('.bam', '.haplotype_quality_tagged.bam'), "rb")

# check if haplotype quality tagged BAM file is provided
if not hqtagged_bam:
    print("Tagging haplotype quality (HQ) in BAM file...")
    # call the function to tag haplotype quality
    # uses the global pysamfile object
    tag_haplotype_quality()
    print("Finished tagging haplotype quality (HQ) in BAM file.")
    pysam.sort("-o", phased_bam.replace('.bam', '.haplotype_quality_tagged_sort.bam'), phased_bam.replace('.bam', '.haplotype_quality_tagged.bam'))
    pysam.index(phased_bam.replace('.bam', '.haplotype_quality_tagged_sort.bam'))
    print("Finished indexing haplotype quality tagged BAM file.")

    hqtagged_bam_temp = pysam.AlignmentFile(phased_bam.replace('.bam', '.haplotype_quality_tagged_sort.bam'), "rb")
else:
    print("Using provided haplotype quality tagged BAM file.")
    hqtagged_bam_temp = pysam.AlignmentFile(hqtagged_bam, "rb")


print("Reading somatic VCF and annotating AH...")
# Read VCF and evaluate each denovo SNP
for variant in vcf:
    #somatic_variant = {} # dictionary
    
    chrom = variant.CHROM
    pos = variant.POS
    alt = variant.ALT[0]
    ref = variant.REF
    print("Processing variant:", chrom, pos, "REF:", ref, "ALT:", alt)
    

    AH = assign_AH(chrom, pos, alt, ref)
    if AH == None:
        continue
    else:
        variant.INFO['AH'] = AH
        vcf_o.write_record(variant)

vcf_o.close()
vcf.close()
print("Finished writing output VCF with AH annotations.")


# vcf=/private/groups/patenlab/jimin/data/VCF/HBCC_CARD/HBCC_81992_FTX/v1.9.0_ONT/HBCC_81992_FTX_merged.vcf.gz
# vcf=/private/groups/patenlab/jimin/data/VCF/HBCC_CARD/HBCC_81992_FTX/v1.9.0_ONT/HBCC_81992_FTX_deepsomatic_only_GQ20_DP10_segdup.chr22.vcf.gz
# germline_vcf=/private/groups/patenlab/jimin/data/VCF/HBCC_CARD/HBCC_81992_FTX/v1.9.0_ONT/HBCC_81992_FTX_deepvariant_pass_only.vcf.gz
# bam=/private/groups/patenlab/jimin/data/bams/ont_bams/CARD/HBCC_81992_FTX/HBCC_81992_FTX.haplotagged.bam
# bam=/private/groups/patenlab/jimin/data/bams/ont_bams/CARD/HBCC_81992_FTX/HBCC_81992_FTX.chr1.haplotagged.bam
# bam=/private/groups/patenlab/jimin/data/bams/ont_bams/CARD/HBCC_81992_FTX/HBCC_81992_FTX.chr22.haplotagged.bam

# hqtaggedbam=/private/groups/patenlab/jimin/data/bams/ont_bams/CARD/HBCC_81992_FTX/HBCC_81992_FTX.chr22.haplotagged.haplotype_quality_tagged_sort.bam
# python3 tag_haplotype_quality.py -bam $bam -hqtaggedbam $hqtaggedbam -v $vcf -g $germline_vcf -k 3 -o test.hq.vcf.gz

# time python3 tag_haplotype_quality.py -bam $bam -hqtaggedbam $hqtaggedbam -v $vcf -g $germline_vcf -k 3 -o test.hq.vcf.gz
# time python3 tag_haplotype_quality.py -bam $bam -hqtaggedbam $hqtaggedbam -v $vcf -g $germline_vcf -k 3 -o test.hq.thu.vcf.gz

# germline_vcf=/private/groups/patenlab/jimin/data/VCF/HBCC_CARD/HBCC_81992_FTX/v1.9.0_ONT/HBCC_81992_FTX_deepvariant_pass_only.vcf.gz
# vcf=/private/groups/patenlab/jimin/data/VCF/HBCC_CARD/HBCC_81992_FTX/v1.9.0_ONT/HBCC_81992_FTX_deepsomatic_only_GQ20_DP10_segdup.chr22.vcf.gz
# bam=/private/groups/patenlab/jimin/data/bams/ont_bams/CARD/HBCC_81992_FTX/HBCC_81992_FTX.chr22.haplotagged.bam
# time python3 tag_haplotype_quality.py -bam $bam -v $vcf -g $germline_vcf -k 3 -o test.hq.thu.vcf.gz
