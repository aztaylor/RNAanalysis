#!/bin/bash
for fn in FASTQ_Generation_2023-09-06_06_48_25Z-691456767/*;
do
samp=`basename ${fn}`
echo ${samp}
echo "Processing sample ${samp}"
salmon quant -i MG1655_index -l A \
         -r ${fn}/${samp:11}_R1_001.fastq.gz \
         -p 8 --validateMappings -o quants/${samp}_quant --gcBias
done