#!/bin/bash

if [ $# -lt 2 ];then

	echo "usage: "$0" outdir threshold";

	exit 0;

else

	########## Directories and CONSTANTS  here  ################
	outdir=$1;
	RMSthreshold=$2;
	
	stats_dir=${outdir}/stats;
	QMov_dir=${outdir}/QMov;

	############ RMS ################

	mkdir -p $stats_dir;
	mkdir -p $QMov_dir;
	echo $RMSthreshold > $QMov_dir/RMSthreshold.par

	echo "RMS relative displacement";

	# Add first volumen to RMS relative displacements:
	echo '0' | cat - ${outdir}/prefiltered_func_data_mcf_rel.rms > ${QMov_dir}/RMSdisp.txt;

	awk -v var="$RMSthreshold" '{if ($1 > var ) print 1; else print 0}' ${QMov_dir}/RMSdisp.txt > ${QMov_dir}/RMS_peaks.txt;
	cat ${QMov_dir}/RMS_peaks.txt | grep 1 | wc -l > ${QMov_dir}/NVols_RMSpeaks.txt;

fi
