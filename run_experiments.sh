#!/bin/bash

mkdir $1.log_folder
#
denoising_stats_log="log" #$1.log_folder/"overall_denoising_karect_style.log_file"
#
for t in 8
do
for n in 6 7 8 9 10
        do
	    ./denoise -d $1".fastq" -t 64 -1 $t -2 $n
            ./util/get_denoising_stats.out $1 "t"$t."n"$n >> $denoising_stats_log
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
