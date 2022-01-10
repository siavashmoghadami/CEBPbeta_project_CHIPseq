# Last Modified: 7 Jan 2022
##################################################################################################################################################################################################
# used shell
# !/bin/zsh
##################################################################################################################################################################################################
# required directories
project_dir=/path/to/project/dir
threads={number of available threads} # number of threads
fastq_in=$project_dir/fastqFiles
fastq_in_pseudoreplicates=$project_dir/fastqFiles/pseudoreplicates
fastqc_out=$project_dir/analysis/qc
genome_dir=$project_dir/analysis/genomeRef/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
genome_homer=$project_dir/analysis/genomeRef/Mus_musculus.GRCm39.dna.primary_assembly.fa
genome_ref=$project_dir/analysis/genomeRef/Mus_musculus.GRCm39.dna.primary_assembly
align_out=$project_dir/analysis/alignment
bam_out=$project_dir/analysis/samToBamConversion
sorted_bam=$project_dir/analysis/sortingBAM
filtered_sorted_bam=$project_dir/analysis/filteringBAM
multiqc_in=$project_dir
multiqc_out=$project_dir/analysis/multiqc
macs3_out=$project_dir/analysis/peakcalling
macs3_diff_out=$project_dir/analysis/differentialPeakCalling
bamCoverage_out=$project_dir/analysis/visualizations/bamCoverage
bamCompare_out=$project_dir/analysis/visualizations/bamCompare
computeMatrix_out=$project_dir/analysis/visualizations/computeMatrix
multiBamSummary_out=$project_dir/analysis/ChIPQCreport/multiBamSummary
diffBind=$project_dir/analysis/diffBind
macs3_bdgcmp=$project_dir/analysis/MACS3_bdgcmp
macs3_bdgdiff=$project_dir/analysis/MACS3_bdgdiff
homer_out=$project_dir/analysis/homer
##################################################################################################################################################################################################
# making three pseudoreplicates from each biological sample
for file in $fastq_in/*.fastq.gz
do
	sample_name=`basename $file .fastq.gz`
	fastqsplitter -t $threads -i $fastq_in/$sample_name".fastq.gz" -o $fastq_in_pseudoreplicates/$sample_name"_PR1.fastq.gz" -o $fastq_in_pseudoreplicates/$sample_name"_PR2.fastq.gz" -o $fastq_in_pseudoreplicates/$sample_name"_PR3.fastq.gz"
done
##################################################################################################################################################################################################
# quality assessment of all fastq files
fastqc $fastq_in/*.fastq.gz -t $threads --outdir=$fastqc_out
fastqc $fastq_in_pseudoreplicates/*.fastq.gz -t $threads --outdir=$fastqc_out
##################################################################################################################################################################################################
# make Bowtie2 Index from Mus_musculus.GRCm39.dna.primary_assembly.fa.gz (Ensemble)
bowtie2-build --threads $threads $genome_dir $genome_ref
##################################################################################################################################################################################################
for file in $fastq_in_pseudoreplicates/*.fastq.gz
do
	sample_name=`basename $file .fastq.gz`

	# alignment to Mus_musculus.GRCm39.dna.primary_assembly.fa.gz (Ensemble) using Bowtie2
	bowtie2 -p $threads -q --local -x $genome_ref -U $fastq_in/$sample_name".fastq.gz" -S $align_out/$sample_name.sam &> $align_out/$sample_name.log

	# SAM to BAM conversion using samtools
	samtools view -h -S -b -@ $threads -o $bam_out/$sample_name.bam $align_out/$sample_name.sam

	# sorting BAM files by genomics coordinates using sambamba
	sambamba sort -t $threads -o $sorted_bam/$sample_name"_sorted.bam" $bam_out/$sample_name.bam

	# filtering uniquely mapping reads using sambamba
	sambamba view -h -t $threads -f bam -F "[XS] == null and not unmapped  and not duplicate" $sorted_bam/$sample_name"_sorted.bam" > $filtered_sorted_bam/$sample_name"_sorted_filtered.bam"

	# make indices for filtered/sorted bam files using samtools
	samtools index -@ $threads $filtered_sorted_bam/$sample_name"_sorted_filtered.bam"

	# make bigWig files for visualization with bamCoverage
	bamCoverage -b $filtered_sorted_bam/$sample_name"_sorted_filtered.bam" -o $bamCoverage_out/$sample_name"_sorted_filtered.bw" --binSize 5 --normalizeUsing BPM --smoothLength 30 --extendReads 200 -p $threads --centerReads 2> $bamCoverage_out/$sample_name"_sorted_filtered.log"

	# make computeMatrix with deepTools
	computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 -R $project_dir/analysis/genomeRef/Mus_musculus.GRCm39.105.gtf.gz -S $bamCoverage_out/$sample_name"_sorted_filtered.bw" --skipZeros -o $computeMatrix_out/$sample_name"_sorted_filtered.gz" -p $threads --outFileSortedRegions $computeMatrix_out/$sample_name"_sorted_filtered.bed"

	# visualize the results with deepTools
	plotProfile -m $computeMatrix_out/$sample_name"_sorted_filtered.gz" -out $computeMatrix_out/$sample_name"_sorted_filtered_line.png" --perGroup --colors red --plotTitle "" --samplesLabel "$sample_name" --refPointLabel "TSS" -T "$sample_name" -z ""
	plotHeatmap -m $computeMatrix_out/$sample_name"_sorted_filtered.gz" -out $computeMatrix_out/$sample_name"_sorted_filtered_heatmap.png" --colorMap RdBu --whatToShow 'heatmap and colorbar' --zMin -4 --zMax 4

done
##################################################################################################################################################################################################
# calculate correlation among pseudoreplicates with deepTools
multiBamSummary bins --bamfiles $filtered_sorted_bam/*.bam --smartLabels --binSize 1000 -p $threads --extendReads 200 --centerReads -out $multiBamSummary_out/results.npz 2> $multiBamSummary_out/results.log
plotCorrelation -in $multiBamSummary_out/results.npz --corMethod spearman --skipZeros --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o $multiBamSummary_out/results"_SpearmanCorr_readCounts.png" --outFileCorMatrix $multiBamSummary_out/results"_SpearmanCorr_readCounts.tab"
##################################################################################################################################################################################################
for file in $filtered_sorted_bam/*.bam
do
	sample_name=`basename $file _sorted_filtered.bam`
	
	# make bigWig files for visualization with bamCompare
	bamCompare -b1 $filtered_sorted_bam/$sample_name"_sorted_filtered.bam" -b2 $filtered_sorted_bam/Pooled_Input_sorted_filtered.bam -o $bamCompare_out/$sample_name".bw" --binSize 10 --normalizeUsing BPM --scaleFactorsMethod None --smoothLength 30 --extendReads 200 -p $threads --centerReads 2> $bamCompare_out/$sample_name".log"
	
	# peak calling with MACS3
	macs3 callpeak -t $filtered_sorted_bam/$sample_name"_sorted_filtered.bam" -c $filtered_sorted_bam/Pooled_Input_sorted_filtered.bam -f BAM -g mm -n $sample_name --bdg --outdir $macs3_out 2> $macs3_out/$sample_name".log"

	# fold-enrichment analysis using MACS3 bdgcmp
	macs3 bdgcmp -t $macs3_out/$sample_name"_treat_pileup.bdg" -c $macs3_out/$sample_name"_control_lambda.bdg" -o $macs3_bdgcmp/$sample_name"_FE.bdg" -m FE
done
##################################################################################################################################################################################################
# compare KI Ischemia/KI Sham with WT Ischemia/WT Sham
bamCompare -b1 $filtered_sorted_bam/KI_Ischemia_sorted_filtered.bam -b2 $filtered_sorted_bam/WT_Ischemia_sorted_filtered.bam -o $bamCompare_out/KI_Isch_vs_WT_Isch.bw --binSize 10 --normalizeUsing BPM --scaleFactorsMethod None --smoothLength 30 --extendReads 200 -p $threads --centerReads 2> $bamCompare_out/KI_Isch_vs_WT_Isch.log
bamCompare -b1 $filtered_sorted_bam/KI_Sham_sorted_filtered.bam -b2 $filtered_sorted_bam/WT_Sham_sorted_filtered.bam -o $bamCompare_out/KI_Sham_vs_WT_Sham.bw --binSize 10 --normalizeUsing BPM --scaleFactorsMethod None --smoothLength 30 --extendReads 200 -p $threads --centerReads 2> $bamCompare_out/KI_Sham_vs_WT_Sham.log
macs3 bdgdiff --t1 $macs3_out/KI_Ischemia_treat_pileup.bdg --c1 $macs3_out/KI_Ischemia_control_lambda.bdg --t2 $macs3_out/WT_Ischemia_treat_pileup.bdg --c2 $macs3_out/WT_Ischemia_control_lambda.bdg --d1 7285137 --d2 13719244 --outdir $macs3_bdgdiff --o-prefix diff_KI_Isch_WT_Isch
macs3 bdgdiff --t1 $macs3_out/KI_Sham_treat_pileup.bdg --c1 $macs3_out/KI_Sham_control_lambda.bdg --t2 $macs3_out/WT_Sham_treat_pileup.bdg --c2 $macs3_out/WT_Sham_control_lambda.bdg --d1 8960107 --d2 7282373 --outdir $macs3_bdgdiff --o-prefix diff_KI_Sham_WT_Sham
##################################################################################################################################################################################################
for file in $macs3_out/*.bed
do
	sample_name=`basename $file _summits.bed`
	
	# motif analysis with Homer
	/Applications/homer/bin/findMotifsGenome.pl $macs3_out/$sample_name"_summits.bed" $genome_homer $homer_out/$sample_name -size 200 -preparse
done
##################################################################################################################################################################################################
# make a report
multiqc -d $multiqc_in -o $multiqc_out
