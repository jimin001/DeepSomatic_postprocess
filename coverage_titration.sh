Help()
{
   # Display Help
   echo "This script performs coverage titrations on tumor-normal bam files."
   echo "Creates coverage titration bams at : 30x, 60x, 90x for tumor bams"
   echo "and creates coverage titration bams at : 20x, 25x, 30x for normal bams."
   echo
   echo "Usage: ./coverage_titration.sh \
   -t <tumor_bam> \
   -n <normal_bam> \
   -c <tumor_coverage> \
   -q <normal_coverage> \
   -p <platform> \
   -s <sample> \
   -o <output_directory> "
   echo
   echo "required:"
   echo "-t  path to input tumor bam"
   echo "-n  path to input normal bam"
   echo "-c  coverage of input tumor bam"
   echo "-q  coverage of input normal bam"
   echo "-p  sequencing platform of data (Illumina, HiFi, ONT)"
   echo "-s  sample name"
   echo "-o  path to output directory"
   echo
}

while getopts :h option
do
   case "${option}" in
      h) Help
         exit;;
   esac
done

while getopts t:n:c:q:p:s:o: flag
do
    case "${flag}" in
        t) tumor_bam=${OPTARG};;
        n) normal_bam=${OPTARG};;
        c) tumor_coverage=${OPTARG};;
        q) normal_coverage=${OPTARG};;
        p) platform=${OPTARG};;
        s) sample=${OPTARG};;
        o) output_directory=${OPTARG};;
    esac
done

echo "Tumor bam: $tumor_bam";
echo "Normal bam: $normal_bam";
echo "Tumor coverage: $tumor_coverage";
echo "Normal coverage: $normal_coverage";
echo "Sequencing platform: $platform";
echo "Cell line sample: $sample";
echo "Output directory: $output_directory";

# Set the exit code of a pipeline to that of the rightmost command
# to exit with a non-zero status, or zero if all commands of the pipeline exit
set -o pipefail
# cause a bash script to exit immediately when a command fails
set -e
# cause the bash shell to treat unset variables as an error and exit immediately
set -u

echo "---------------------------------------------------------------------"
echo "Creating coverage titration tumor bams: 30x, 60x, 90x"
echo "-------------------------------------"

# downsample tumor bam to: 30x, 60x, 90x coverage
tumor_goal_list="30 60 90"

echo "Downsampling tumor bams..."
for tumor_goal_coverage in $tumor_goal_list
do
    echo "Downsample to: ${tumor_goal_coverage}x"

    if [ $(($tumor_coverage - $tumor_goal_coverage)) -lt 3 ]
    then
        echo "Downsampling unnecessary"
    else
        echo "Downsampling with samtools..."
        tumor_fraction=$( bc <<< "scale=3;$tumor_goal_coverage/$tumor_coverage")
        echo "Tumor fraction: $tumor_fraction"

        samtools view -@30 -b -s ${tumor_fraction} ${tumor_bam} > ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${tumor_goal_coverage}x.bam
        samtools index -@30 ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${tumor_goal_coverage}x.bam

        echo "Downsampled bam name: ${sample}_${platform}_Tumor.GRCh38.${tumor_goal_coverage}x.bam"

        echo "Verifying downsampling with samtools depth..."
        samtools depth ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${tumor_goal_coverage}x.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'
        echo "-------------------------------------"
    fi
done

echo "---------------------------------------------------------------------"
echo "Creating coverage titration normal bams: 20x, 25x, 30x"
echo "-------------------------------------"
echo "Downsampling normal bams..."
# downsample normal bam to: 20x, 25x, 30x coverage
normal_goal_list="20 25 30"

for normal_goal_coverage in $normal_goal_list
do
    echo "Downsample to: ${normal_goal_coverage}x"

    if [ $(($normal_coverage - $normal_goal_coverage)) -lt 2 ]
    then
        echo "Downsampling unnecessary"
    else
        echo "Downsampling with samtools..."
        normal_fraction=$( bc <<< "scale=3;$normal_goal_coverage/$normal_coverage")
        echo "normal fraction: $normal_fraction"

        samtools view -@30 -b -s ${normal_fraction} ${normal_bam} > ${output_directory}/${sample}_${platform}_Normal.GRCh38.${normal_goal_coverage}x.bam
        samtools index -@30 ${output_directory}/${sample}_${platform}_Normal.GRCh38.${normal_goal_coverage}x.bam

        echo "Downsampled bam name: ${sample}_${platform}_Normal.GRCh38.${normal_goal_coverage}x.bam"

        echo "Verifying downsampling with samtools depth..."
        samtools depth ${output_directory}/${sample}_${platform}_Normal.GRCh38.${normal_goal_coverage}x.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'
        echo "-------------------------------------"
    fi
done
