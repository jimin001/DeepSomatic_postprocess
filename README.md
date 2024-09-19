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
