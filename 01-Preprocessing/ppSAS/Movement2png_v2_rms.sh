#!/bin/bash

if [ $# -eq 0 ];then

	echo "usage: "$0" outdir";

	exit 0;

else
	if [ -e $1 ];then

		########## Directories and CONSTANTS  here  ################
		outdir=$1;
		stats_dir=${outdir}/stats;
		QMov_dir=${outdir}/QMov;
		id=$(basename $1 .rsfMRIv2_nihpd);	

		RMSthreshold=$(cat ${QMov_dir}/RMSthreshold.par);
		
		#############################################
	
		mkdir -p ${QMov_dir};

		awk -v var="$RMSthreshold" '{print $1*0+var}' ${QMov_dir}/RMSdisp.txt > ${QMov_dir}/thresh.txt;
		paste ${QMov_dir}/RMSdisp.txt ${QMov_dir}/thresh.txt > ${QMov_dir}/RMSdisp_thresh.txt;
		rm -f ${QMov_dir}/thresh.txt;

		${FSLDIR}/bin/fsl_tsplot -i ${outdir}/prefiltered_func_data_mcf.par -t 'MCFLIRT estimated rotations (radians)' -u 1 --start=1 --finish=3 -a x,y,z -w 640 -h 144 -o ${QMov_dir}/rot.png;

		${FSLDIR}/bin/fsl_tsplot -i ${outdir}/prefiltered_func_data_mcf.par -t 'MCFLIRT estimated translations (mm)' -u 1 --start=4 --finish=6 -a x,y,z -w 640 -h 144 -o ${QMov_dir}/trans.png;

		${FSLDIR}/bin/fsl_tsplot -i ${QMov_dir}/RMSdisp_thresh.txt -t 'Relative Displacement RMS (mm)' -u 1 --start=1 --finish=2 -a RMS,threshold -w 640 -h 144 -o ${QMov_dir}/relRMSd.png;

		montage ${QMov_dir}/relRMSd.png ${QMov_dir}/rot.png ${QMov_dir}/trans.png -tile x3 -geometry 640x144+2+2 ${QMov_dir}/montage_{$id}.png;

		##### mean and stdev
		echo -e "\t\tMean\t\tStDev" > ${QMov_dir}/relRMSd_mean_std.txt;
		awk '{ if ($1 > 0) print $1}' ${QMov_dir}/RMSdisp_thresh.txt | awk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { print "RMS \t\t"avg"\t\t"sqrt(mean2 / NR); }' >> ${QMov_dir}/relRMSd_mean_std.txt;
	
	
	else
		echo "ERROR: directory doesn't exist";
	fi
fi
