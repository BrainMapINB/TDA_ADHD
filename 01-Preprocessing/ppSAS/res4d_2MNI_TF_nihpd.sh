#!/bin/bash

if [ $# -lt 6 ];then

	echo "usage: "$0" basedir outdir outputbasename HPFc(Hz) LPFc(Hz) KeepIntermMNIfiles(0|1)";
	echo "";

	exit 0;

else

	########## Directories and CONSTANTS  here  ################
	nbasedir=$1
	outdir=$2
	outbn=$3;
	hpf=$4;
	lpf=$5;
	keep=$6;
	
	outdir_bn=$(basename $2 .rsfMRIv2_nihpd);
	stats_dir=${outdir}/stats;
	ppBold_dir=${outdir}/ppBold;
	mkdir -p $ppBold_dir
	TR=$(cat ${outdir}/TR);

	hp_sigma=`echo "scale=2 ;(1/${hpf})/2.35/${TR}" | bc`; # In volumes for fslmaths
	lp_sigma=`echo "scale=2 ;(1/${lpf})/2.35/${TR}" | bc`; # In volumes for fslmaths


	#############################################

	echo "Projecting to Standard space";

	for input in ${stats_dir}/res4d ${stats_dir}/res4d_woGSR; do
		${FSLDIR}/bin/applywarp --ref=/misc/claustrum/zgtabuenca/Tesis/Atlas/nihpd_asym_04.5-18.5_nifti/resample/t1w_brain_2mm.nii.gz --in=${input} --warp=${outdir}/reg/highres2standard_warp.nii.gz --out=${input}_NIHPD2mm --premat=${outdir}/reg/meanfunc2highres.mat --interp=trilinear;
	done

	if [ -e ${outdir}/reg/sym_highres2standard_warp.nii.gz ];then
		for input in ${stats_dir}/res4d ${stats_dir}/res4d_woGSR; do
			${FSLDIR}/bin/applywarp --ref=${nbasedir}/MNI152_T1_4mm_brain.nii --in=${input} --warp=${outdir}/reg/sym_highres2standard_warp.nii.gz --out=${input}_symMNI4mm --premat=${outdir}/reg/meanfunc2highres.mat --interp=trilinear;
		done
	fi

	echo "Bandpass Temporal Filtering to MNI files";
	# Para todos los archivos finales 
	for i in $(ls ${stats_dir}/*NIHPD*gz);do 
	
		en=$(basename $i .nii.gz | cut -d _ -f 2-);

		fslhd -x $i > ${stats_dir}/tmphdr.txt;
		sed -n "s/dt =.*/dt = \'${TR}\'/" ${stats_dir}/tmphdr.txt;
		fslcreatehd ${stats_dir}/tmphdr.txt $i;

		${FSLDIR}/bin/fslmaths $i -bptf $hp_sigma $lp_sigma ${ppBold_dir}/${outbn}_${en}_${outdir_bn};
		
	done

	rm ${stats_dir}/tmphdr.txt;

	case $keep in 

		0)
			echo "Removing Intermediate MNI files (non time-filtering)";
			rm ${stats_dir}/*NIHPD*.*;
			;;
		1)
			echo "Keeping Intermediate MNI files (non time-filtering) in" $stats_dir;
			;;
		*)
			echo "No option selected for keeping Intermediate MNI files, they will be removed !!"
			rm ${stats_dir}/*NIHPD*.*;			
			;;
	esac
fi
