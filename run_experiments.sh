#!/bin/bash

mkdir $1.log_folder
#
denoising_stats_log="log_1" #$1.log_folder/"overall_denoising_karect_style.log_file"
karect_tool=/data/shubham/assembly/karect/karect

#/usr/bin/time -v $karect_tool -correct -inputfile=$1"_1.fastq" -inputfile=$1"_2.fastq" -resultdir=assembly_data -celltype=diploid -matchtype=hamming -threads=32 | tee $1".karect_log"
#$karect_tool -align -threads=32 -matchtype=hamming -inputfile=$1"_1.fastq" -inputfile=$1"_2.fastq" -refgenomefile=assembly_data/c_elegans.PRJNA13758.WS241.uppercase.genomic.fa -alignfile=$1"_align.txt"          

#$karect_tool -eval -threads=32 -matchtype=hamming -inputfile=$1"_1.fastq" -inputfile=$1"_2.fastq" -resultfile="assembly_data/karect_SRR823377_1.fastq" -resultfile="assembly_data/karect_SRR823377_2.fastq" -refgenomefile=assembly_data/c_elegans.PRJNA13758.WS241.uppercase.genomic.fa -alignfile=$1"_align.txt" -evalfile=$1"_karect.eval.txt" | tee $1.log_folder/"karect_summary.txt"          


for t in 4 8 12
do
	for n in 4 8 12 16 
        do
		for param3 in 0.95 0.9 0.85
		do
			for param4 in 20 30 50 70 
	    		do	
				./denoise -d $1".fastq" -t 32 -1 $t -2 $n -3 $param3 -4 $param4
	    #util/combine_into_fastq.out $1
	    #/data/shubham/assembly/karect/karect -eval -threads=32 -matchtype=hamming -inputfile=$1"_1.fastq" -inputfile=$1"_2.fastq" -resultfile=$1"_1.harc.fastq" -resultfile=$1"_2.harc.fastq" -refgenomefile=assembly_data/c_elegans.PRJNA13758.WS241.uppercase.genomic.fa -alignfile=$1"_align.txt" -evalfile=$1"_t"$t"_n"$n"_eval.txt" | tee $1.log_folder/"t"$t"_n"$n"alpha_0.8_maxnum50_summary.txt"          
				./util/get_denoising_stats.out $1 "t"$t."n"$n".alpha"$param3".maxlimit"$param4 >> $denoising_stats_log
			done
		done
        done 
done

#denoise_type="counts_and_quality"
#denoising_stats_log=$1.log_folder/overall_denoising.$denoise_type".log_file"
#
#for n in 0.1 0.2
#do
#    for s in 25 30
#    do	
#	for r in 2 4
#    	do
#        	for a in 8 
#       		do
#            		./harc -c $1".fastq" -p -u $denoise_type -n $n -s $s -r $r -a $a -t 16 | tee $1.log_folder/$denoise_type"_n"$n."_s"$s._r$r._a$a.log_file
#            		./harc -d $1.harc -p
#		        ./util/get_denoising_stats.out $1 $1.log_folder/$denoise_type"_n"$n."_s"$s._r$r._a$a >> $denoising_stats_log
#		done
#        done 
#    done
#done
#
#denoise_type="quality_threshold"
#denoising_stats_log=$1.log_folder/overall_denoising.$denoise_type".log_file"
#
#for n in 25 30
#do
#    for r in 2 4
#    do
#        for a in 8
#        do
#            ./harc -c $1".fastq" -p -u $denoise_type -n $n -r $r -a $a -t 16 | tee $1.log_folder/$denoise_type"_n"$n._r$r._a$a.log_file
#            ./harc -d $1.harc -p
#	   ./util/get_denoising_stats.out $1 $1.log_folder/$denoise_type"_n"$n._r$r._a$a >> $denoising_stats_log
#        done 
#    done
#done
#
#denoise_type="only_using_counts"
#denoising_stats_log=$1.log_folder/overall_denoising.$denoise_type".log_file"
#
#for n in 0.1 0.2
#do
#    for r in 2 4
#    do
#        for a in 8
#        do
#            ./harc -c $1".fastq" -p -u $denoise_type -n $n -r $r -a $a -t 16 | tee $1.log_folder/$denoise_type"_n"$n._r$r._a$a.log_file
#            ./harc -d $1.harc -p
#	   ./util/get_denoising_stats.out $1 $1.log_folder/$denoise_type"_n"$n._r$r._a$a >> $denoising_stats_log
#        done 
#    done
#done
