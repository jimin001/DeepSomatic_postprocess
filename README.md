# DeepSomatic postprocess
Filtering scripts for DeepSomatic VCF outputs

## haplotype_somatic_vcf.py
This script was specifically made to post-process non-tumor samples.
This script uses VCF and aligned, phased BAM information to tag INFO field for each position in VCF with 
- AH=1 if alt allele is exclusive to one haplotype
- AH=2 is alt allele is present on both haplotypes

run locally:
```
Usage:
python3 haplotype_somatic_vcf.py -v ${VCF} -bam ${BAM} -o ${OUTPUT} -i ${indel_read_threshold} -u ${unphased_read_threshold}
```
## deepsomatic_tumor_only_postfiltering_script.sh
This bash script contains commands to filter non-tumor samples for high confidence somatic variants.

Steps:
1. subtract germline variants from somatic variant set
2. filter out variants GQ<20, DP<10
3. subtract segdup regions
4. tag with AH=2: cases where alternate variants are on both haplotypes, AH=1: cases where alternate variants are on only 1 haplotype
5. tag with gnomAD annotations for total population allele frequency
6. filter for gnomAD annotation AF < 0.001 and AH=1

**need to update paths to `GRCh38_segdups.bed` and `haplotypeonly_somatic_vcf.py` in script**

run locally:
```
Usage:
sample="sample name"
echo "FILTERING SAMPLE: $sample"
BAM=${sample}.haplotagged.bam
germline_VCF=${sample}_deepvariant.phased.vcf.gz
somatic_VCF=${sample}_wg_pon_tumor-only.vcf.gz

OUTPUT_DIR=/path/to/output/directory

./deepsomatic_tumor_only_postfiltering_script.sh -g ${germline_VCF} -s ${somatic_VCF} -b ${BAM} -n ${sample} -o ${OUTPUT_DIR}
```

## filter_bam_cal_VAF.sh
This script calculates variant allele frequencies (VAF) of somatic variants using read depth from bam files.
This script consists of 4 steps (which step to start on can be chosen using '-w' flag):
1. tag depth information for each variant in VCF using information from tumor and normal bams
2. calculate vaf using depth information
3. use awk to filter out variants that do not pass given vaf threshold
4. use targets file and bcftools filter to filter VCF for only variants in the targets file

run locally:
```
Usage:
sample="sample name"
TUMOR_BAM="path to tumor bam"
NORMAL_BAM="path to normal bam"
vcf="path to VCF"
filter_vaf="minimum VAF threshold"
output_prefix="desired output prefix"
output_directory="path to output directory"

/private/groups/patenlab/jimin/scripts/deepsomatic/filter_bam_cal_VAF.sh -t ${TUMOR_BAM} -n ${NORMAL_BAM} -v ${vcf} -f ${filter_vaf} -s ${sample} -p ${output_prefix} -o ${output_directory} -w 1

```



