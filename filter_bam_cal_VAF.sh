while getopts t:n:v:f:s:p:o:w: flag
do
    case "${flag}" in
        t) tumor_bam=${OPTARG};;
        n) normal_bam=${OPTARG};;
        v) vcf=${OPTARG};;
        f) filter_vaf=${OPTARG};;
        s) sample=${OPTARG};;
        p) output_prefix=${OPTARG};;
        o) output_directory=${OPTARG};;
        w) start_step=${OPTARG};;
    esac
done

echo "tumor bam: $tumor_bam";
echo "normal bam: $normal_bam";
echo "vcf to filter: $vcf";
echo "minimum vaf threshold to filter: $filter_vaf"
echo "sample: $sample";
echo "output_prefix: $output_prefix";
echo "output_directory: $output_directory";
echo "step to start analysis from: $start_step"

set -o pipefail
set -e
set -u

mkdir -p ${output_directory}
export PYTHONPATH="${PYTHONPATH}:/private/groups/patenlab/jimin/scripts/ClairS"

if [[ $start_step -eq 1 ]]
then
    # 1. tag depth information for each variant in VCF using information from tumor and normal bams
    python3 /private/groups/patenlab/jimin/scripts/ClairS/src/cal_af_distribution.py --tumor_bam_fn $tumor_bam --normal_bam_fn $normal_bam --truth_vcf_fn $vcf \
        --output_path ${output_directory}/${output_prefix}_tag_depth.vcf --threads 16

    # 2. calculate vaf using depth information
    DEPTH_FILE=${output_directory}/${output_prefix}_tag_depth.vcf
    VAF_PATH=${output_directory}/${output_prefix}_tag_vaf.bed
    python3 /private/groups/patenlab/jimin/scripts/deepsomatic/calculate_vaf.py --vaf_file ${DEPTH_FILE} --outfile_path ${VAF_PATH}

    # 3. use awk to filter out variants that do not pass given vaf threshold
    TARGETS_FILE=${output_directory}/${output_prefix}_filter_minVAF${filter_vaf}.bed
    awk -v var="$filter_vaf" '$4>var' ${VAF_PATH} > ${TARGETS_FILE}


    # 4. use targets file and bcftools filter to filter VCF for only variants in the targets file
    FILTERED_VCF=${output_directory}/${output_prefix}_tag_vaf_filter_minVAF${filter_vaf}.vcf.gz
    bcftools filter -T ${TARGETS_FILE} ${vcf} | bgzip > ${FILTERED_VCF}
    bcftools index -t ${FILTERED_VCF}
elif [[ $start_step -eq 2 ]]
then
    # 2. calculate vaf using depth information
    DEPTH_FILE=${output_directory}/${output_prefix}_tag_depth.vcf
    VAF_PATH=${output_directory}/${output_prefix}_tag_vaf.bed
    python3 /private/groups/patenlab/jimin/scripts/deepsomatic/calculate_vaf.py --vaf_file ${DEPTH_FILE} --outfile_path ${VAF_PATH}

    # 3. use awk to filter out variants that do not pass given vaf threshold
    TARGETS_FILE=${output_directory}/${output_prefix}_filter_minVAF${filter_vaf}.bed
    awk -v var="$filter_vaf" '$4>var' ${VAF_PATH} > ${TARGETS_FILE}


    # 4. use targets file and bcftools filter to filter VCF for only variants in the targets file
    FILTERED_VCF=${output_directory}/${output_prefix}_tag_vaf_filter_minVAF${filter_vaf}.vcf.gz
    bcftools filter -T ${TARGETS_FILE} ${vcf} | bgzip > ${FILTERED_VCF}
    bcftools index -t ${FILTERED_VCF}
elif [[ $start_step -eq 3 ]]
then
    # 3. use awk to filter out variants that do not pass given vaf threshold
    VAF_PATH=${output_directory}/${output_prefix}_tag_vaf.bed
    TARGETS_FILE=${output_directory}/${output_prefix}_filter_minVAF${filter_vaf}.bed
    awk -v var="$filter_vaf" '$4>var' ${VAF_PATH} > ${TARGETS_FILE}


    # 4. use targets file and bcftools filter to filter VCF for only variants in the targets file
    FILTERED_VCF=${output_directory}/${output_prefix}_tag_vaf_filter_minVAF${filter_vaf}.vcf.gz
    bcftools filter -T ${TARGETS_FILE} ${vcf} | bgzip > ${FILTERED_VCF}
    bcftools index -t ${FILTERED_VCF}
elif [[ $start_step -eq 4 ]]
then
    # 4. use targets file and bcftools filter to filter VCF for only variants in the targets file
    TARGETS_FILE=${output_directory}/${output_prefix}_filter_minVAF${filter_vaf}.bed
    FILTERED_VCF=${output_directory}/${output_prefix}_tag_vaf_filter_minVAF${filter_vaf}.vcf.gz
    bcftools filter -T ${TARGETS_FILE} ${vcf} | bgzip > ${FILTERED_VCF}
    bcftools index -t ${FILTERED_VCF}
fi
