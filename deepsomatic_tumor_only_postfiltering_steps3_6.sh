while getopts g:s:b:n:o: flag
do
    case "${flag}" in
        g) germline_VCF=${OPTARG};;
        s) somatic_VCF=${OPTARG};;
        b) BAM=${OPTARG};;
        n) sample=${OPTARG};;
        o) OUTPUT_DIR=${OPTARG};;
    esac
done

# Set the exit code of a pipeline to that of the rightmost command
# to exit with a non-zero status, or zero if all commands of the pipeline exit
set -o pipefail
# cause a bash script to exit immediately when a command fails
set -e
# cause the bash shell to treat unset variables as an error and exit immediately
set -u

mkdir -p ${OUTPUT_DIR}

# 3. subtract segdup regions
cd ${OUTPUT_DIR}/${sample}
DS_VCF=${OUTPUT_DIR}/${sample}/${sample}_deepsomatic_only_GQ20_DP10.vcf.gz
SEGDUP_BED=/private/groups/patenlab/jimin/data/BED/GIAB/GRCh38_segdups.bed

bcftools view -h ${DS_VCF} > ${sample}_header_temp.vcf
bedtools subtract -a ${DS_VCF} -b ${SEGDUP_BED} > ${sample}_deepsomatic_only_GQ20_DP10_segdup_temp.vcf
cat ${sample}_header_temp.vcf ${sample}_deepsomatic_only_GQ20_DP10_segdup_temp.vcf | bgzip > ${sample}_deepsomatic_only_GQ20_DP10_segdup.vcf.gz
bcftools index -t ${sample}_deepsomatic_only_GQ20_DP10_segdup.vcf.gz
rm *temp*

# 4. tag with AH=2: cases where alternate variants are on both haplotypes, AH=1: cases where alternate variants are on only 1 haplotype
TOOL=/private/groups/patenlab/jimin/scripts/deepsomatic/haplotypeonly_somatic_vcf.py
VCF=${sample}_deepsomatic_only_GQ20_DP10_segdup.vcf.gz

OUTPUT_VCF=${sample}_deepsomatic_only_GQ20_DP10_segdup_tagAH.vcf

time python3 ${TOOL} -v ${VCF} -bam ${BAM} -o ${OUTPUT_VCF} -u 10

bgzip ${OUTPUT_VCF}
bcftools index -t ${OUTPUT_VCF}.gz

# 5. tag with gnomAD annotations for total population allele frequency
TOOL=/private/home/jpark621/software/snpEff/SnpSift.jar
GNOMAD_VCF=/private/groups/patenlab/jimin/data/VCF/gnomad/gnomad.genomes.r4.0.sites.small.vcf.bgz

VCF=${sample}_deepsomatic_only_GQ20_DP10_segdup_tagAH.vcf.gz
OUTPUT_VCF=${sample}_deepsomatic_only_GQ20_DP10_segdup_tagAH_gnomad.vcf.gz
java -jar ${TOOL} annotate -noId -v ${GNOMAD_VCF} ${VCF} \
| bgzip > ${OUTPUT_VCF}

# 6. filter for gnomAD annotation AF < 0.001 and AH=1
VCF=${sample}_deepsomatic_only_GQ20_DP10_segdup_tagAH_gnomad.vcf.gz
OUTPUT_VCF=${sample}_deepsomatic_only_GQ20_DP10_segdup_tagAH_gnomadrare_AH1.vcf.gz

bcftools view -i 'INFO/AF<0.001 | INFO/AF="."' ${VCF} | bcftools view -i 'INFO/AH=1' | bgzip  > ${OUTPUT_VCF}

# count number of variants at each filter point

echo "Number of variants at each filtering step: " 2>&1 | tee ${sample}_variant_count.log

bcftools view ${sample}_deepsomatic_only.vcf.gz | grep -v "#" | wc -l | tee -a ${sample}_variant_count.log
bcftools view ${sample}_deepsomatic_only_GQ20_DP10.vcf.gz | grep -v "#" | wc -l | tee -a ${sample}_variant_count.log
bcftools view ${sample}_deepsomatic_only_GQ20_DP10_segdup.vcf.gz | grep -v "#" | wc -l | tee -a ${sample}_variant_count.log
bcftools view ${sample}_deepsomatic_only_GQ20_DP10_segdup_tagAH.vcf.gz | grep -v "#" | wc -l | tee -a ${sample}_variant_count.log
bcftools view ${sample}_deepsomatic_only_GQ20_DP10_segdup_tagAH.vcf.gz | grep -v "#" | grep "AH=1" | wc -l | tee -a ${sample}_variant_count.log
bcftools view ${sample}_deepsomatic_only_GQ20_DP10_segdup_tagAH_gnomadrare_AH1.vcf.gz | grep -v "#" | wc -l | tee -a ${sample}_variant_count.log

