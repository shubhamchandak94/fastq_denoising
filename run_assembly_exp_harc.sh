#!/bin/bash

### Inputs
#prefix=$1
#ref_fasta=$2
#denoising_stats_log=$log_folder/overall_denoising_stats.log_file

##### TOOLS to use
spades_tool=../SPAdes-3.11.0-Linux/bin/spades.py
quast_tool=../quast-4.5/quast.py
reckoner_tool=../RECKONER/bin/reckoner
karect_tool=../karect/karect

dirname=assembly_data/
for dataset in SRR823377
do
	
	prefix=$dirname"/"$dataset
	#### Log folders for the experiment
	log_folder=$prefix.log_folder
	mkdir $log_folder
	assembly_output_folder=$prefix.assembly_output
	mkdir $assembly_output_folder
	
	ref_fasta=$dirname"/"$dataset".genome.fasta"

#	for t in 8
#	do
#		for n in 5
#		do
#		#for m in 500 #maxsearch
#		#do
#			suffix="fastq_denoising_t"$t"_n"$n"_alpha0.9_maxnum_50"
#			/usr/bin/time -v ./denoise -d $prefix".fastq" -t 32 -1 $t -2 $n -3 5 -4 50 | tee $log_folder/$suffix
#			python util/combine_into_fastq.py $prefix $suffix 
#			tmp_folder=$assembly_output_folder/$suffix
#			/usr/bin/time -v python $spades_tool --only-assembler -t 32 -1 $prefix"."$suffix"_1.fastq" -2 $prefix"."$suffix"_2.fastq" -o $tmp_folder | tee $log_folder/$suffix".spades_log"
#			quast_output=$log_folder/$suffix".quast_output"
#			python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels $suffix
#		#done
#		done
#	done
#
	for t in 16
	do
		for n in 4
		do
		for param3 in 0.75
		do 
		for param4 in 1000
		do
			suffix="fastq_denoising_new_new"			
#			suffix="fastq_denoising_t"$t"_n"$n"alpha_"$param3"_maxnum"$param4"maxsearch1000"
#			/usr/bin/time -v ./denoise_maxsearch1000 -d $prefix".fastq" -t 32 -1 $t -2 $n -3 $param3 -4 $param4 | tee $log_folder/$suffix
#			python util/combine_into_fastq.py $prefix $suffix 
			tmp_folder=$assembly_output_folder/$suffix
			/usr/bin/time -v python $spades_tool --only-assembler -t 32 -1 $prefix"."$suffix"_1.fastq" -2 $prefix"."$suffix"_2.fastq" -o $tmp_folder | tee $log_folder/$suffix".spades_log"
			quast_output=$log_folder/$suffix".quast_output"
			python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels $suffix
		done
		done
		done
	done
#	for t in 8 16
#	do
#		for n in 5
#		do
#		for param3 in 0.75
#		do 
#		for param4 in 1000
#		do
#			suffix="fastq_denoising_t"$t"_n"$n"alpha_"$param3"_maxnum"$param4"maxsearch20_numdict4_20"
#			/usr/bin/time -v ./denoise_4dict_dictlen20 -d $prefix".fastq" -t 32 -1 $t -2 $n -3 $param3 -4 $param4 | tee $log_folder/$suffix
#			python util/combine_into_fastq.py $prefix $suffix 
#			tmp_folder=$assembly_output_folder/$suffix
#			/usr/bin/time -v python $spades_tool --only-assembler -t 32 -1 $prefix"."$suffix"_1.fastq" -2 $prefix"."$suffix"_2.fastq" -o $tmp_folder | tee $log_folder/$suffix".spades_log"
#			quast_output=$log_folder/$suffix".quast_output"
#			python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels $suffix
#		done
#		done
#		done
#	done
#	for t in 12
#	do
#		for n in 4 5
#		do
#		for param3 in 0.85
#		do 
#		for param4 in 50
#		do
#			suffix="fastq_denoising_t"$t"_n"$n"alpha_"$param3"_maxnum"$param4"maxsearch1000"
#			/usr/bin/time -v ./denoise_maxsearch1000 -d $prefix".fastq" -t 32 -1 $t -2 $n -3 $param3 -4 $param4 | tee $log_folder/$suffix
#			python util/combine_into_fastq.py $prefix $suffix 
#			tmp_folder=$assembly_output_folder/$suffix
#			/usr/bin/time -v python $spades_tool --only-assembler -t 32 -1 $prefix"."$suffix"_1.fastq" -2 $prefix"."$suffix"_2.fastq" -o $tmp_folder | tee $log_folder/$suffix".spades_log"
#			quast_output=$log_folder/$suffix".quast_output"
#			python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels $suffix
#		done
#		done
#		done
#	done
#	for t in 8
#	do
#		for n in 4
#		do
#		for param3 in 0.85
#		do 
#		for param4 in 50
#		do
#			suffix="fastq_denoising_t"$t"_n"$n"alpha_"$param3"_maxnum"$param4"maxsearch20_dictlen15"
#			/usr/bin/time -v ./denoise_dictlen15 -d $prefix".fastq" -t 32 -1 $t -2 $n -3 $param3 -4 $param4 | tee $log_folder/$suffix
#			python util/combine_into_fastq.py $prefix $suffix 
#			tmp_folder=$assembly_output_folder/$suffix
#			/usr/bin/time -v python $spades_tool --only-assembler -t 32 -1 $prefix"."$suffix"_1.fastq" -2 $prefix"."$suffix"_2.fastq" -o $tmp_folder | tee $log_folder/$suffix".spades_log"
#			quast_output=$log_folder/$suffix".quast_output"
#			python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels $suffix
#		done
#		done
#		done
#	done
#	for t in 8
#	do
#		for n in 4
#		do
#		for param3 in 0.85
#		do 
#		for param4 in 50
#		do
#			suffix="fastq_denoising_t"$t"_n"$n"alpha_"$param3"_maxnum"$param4"maxsearch20_4dict_15"
#			/usr/bin/time -v ./denoise_4dict_dictlen15 -d $prefix".fastq" -t 32 -1 $t -2 $n -3 $param3 -4 $param4 | tee $log_folder/$suffix
#			python util/combine_into_fastq.py $prefix $suffix 
#			tmp_folder=$assembly_output_folder/$suffix
#			/usr/bin/time -v python $spades_tool --only-assembler -t 32 -1 $prefix"."$suffix"_1.fastq" -2 $prefix"."$suffix"_2.fastq" -o $tmp_folder | tee $log_folder/$suffix".spades_log"
#			quast_output=$log_folder/$suffix".quast_output"
#			python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels $suffix
#		done
#		done
#		done
#	done
done
