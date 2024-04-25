# DeepSomatic
Keep track of scripts used for DeepSomatic analysis

## coverage_titration.sh
This script takes in two bam files (`normal.bam`, `tumor.bam`) and downsamples the `tumor.bam` to 30x, 60x, 90x coverage and downsamples the `normal.bam` to 20x, 25x, 30x coverage.
This script uses samtools with 30 threads.

run locally:
```
Usage: ./coverage_titration.sh \
   -t <tumor_bam> \
   -n <normal_bam> \
   -c <tumor_coverage> \
   -q <normal_coverage> \
   -p <platform> \
   -s <sample> \
   -o <output_directory> 
```

## purity_titration.sh
This script takes in two bam files  (`normal.bam`, `tumor.bam`) and creates purity titration bams.
Creates purity titration bams at : 60%, 70%, 80%, 90% tumor purity and 90%, 95% normal purity.
This script uses samtools with 30 threads.

run locally:
```
Usage: ./purity_titration.sh \
   -t <tumor_bam> \
   -n <normal_bam> \
   -c <tumor_coverage> \
   -q <normal_coverage> \
   -g <tumor_goal_total_coverage> \
   -x <normal_goal_total_coverage> \
   -p <platform> \
   -s <sample> \
   -o <output_directory>
```
