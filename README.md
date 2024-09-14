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

