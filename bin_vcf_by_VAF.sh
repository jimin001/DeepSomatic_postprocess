while getopts b:p:o:v: flag
do
    case "${flag}" in
        b) vaf_tagged_bed=${OPTARG};;
        p) output_prefix=${OPTARG};;
        o) output_directory=${OPTARG};;
        v) vcf=${OPTARG};;
    esac
done

# script to generate VCFs binned by VAF


echo "vaf_tagged_bed: $vaf_tagged_bed";
echo "output_prefix: $output_prefix";
echo "output_directory: $output_directory";
echo "vcf: $vcf"

set -o pipefail
set -e
set -u

# filter by binned VAF

VAF_range=0_0.1
TARGETS_FILE=${output_directory}/${output_prefix}_${VAF_range}.bed

awk '$4>0 && $4<=0.1' $vaf_tagged_bed > ${TARGETS_FILE}

VAF_range=0.1_0.2
TARGETS_FILE=${output_directory}/${output_prefix}_${VAF_range}.bed

awk '$4>0.1 && $4<=0.2' $vaf_tagged_bed > ${TARGETS_FILE}

VAF_range=0.2_0.3
TARGETS_FILE=${output_directory}/${output_prefix}_${VAF_range}.bed

awk '$4>0.2 && $4<=0.3' $vaf_tagged_bed > ${TARGETS_FILE}

VAF_range=0.3_0.4
TARGETS_FILE=${output_directory}/${output_prefix}_${VAF_range}.bed

awk '$4>0.3 && $4<=0.4' $vaf_tagged_bed > ${TARGETS_FILE}

VAF_range=0.4_0.5
TARGETS_FILE=${output_directory}/${output_prefix}_${VAF_range}.bed

awk '$4>0.4 && $4<=0.5' $vaf_tagged_bed > ${TARGETS_FILE}

VAF_range=0.5_0.6
TARGETS_FILE=${output_directory}/${output_prefix}_${VAF_range}.bed

awk '$4>0.5 && $4<=0.6' $vaf_tagged_bed > ${TARGETS_FILE}

VAF_range=0.6_0.7
TARGETS_FILE=${output_directory}/${output_prefix}_${VAF_range}.bed

awk '$4>0.6 && $4<=0.7' $vaf_tagged_bed > ${TARGETS_FILE}

VAF_range=0.7_0.8
TARGETS_FILE=${output_directory}/${output_prefix}_${VAF_range}.bed

awk '$4>0.7 && $4<=0.8' $vaf_tagged_bed > ${TARGETS_FILE}

VAF_range=0.8_1
TARGETS_FILE=${output_directory}/${output_prefix}_${VAF_range}.bed

awk '$4>0.8 && $4<=1' $vaf_tagged_bed > ${TARGETS_FILE}


# 4. use targets file and bcftools filter to filter VCF for only variants in the targets file

for VAF_range in 0_0.1 0.1_0.2 0.2_0.3 0.3_0.4 0.4_0.5 0.5_0.6 0.6_0.7 0.7_0.8 0.8_1
do
	TARGETS_FILE=${output_directory}/${output_prefix}_${VAF_range}.bed
	FILTERED_VCF=${output_directory}/${output_prefix}_${VAF_range}.vcf.gz
	bcftools filter -T ${TARGETS_FILE} ${vcf} | bgzip > ${FILTERED_VCF}
	bcftools index -t ${FILTERED_VCF}
done
