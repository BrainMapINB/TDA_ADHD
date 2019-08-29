#!/bin/bash

if [ $# -lt 2 ];then

	echo "usage: "$0" mprage outdir";
	echo "";

	exit 0;

else
	if [ -e $1 ];then
		########## Directories and CONSTANTS  here  ################
		mprage=$1
		mprage_bn=$(basename $1 .nii.gz);
		outdir=$2;

		cp $mprage ${outdir}/reg/${mprage_bn}_head.nii.gz
		rm $mprage

		highres_head=${outdir}/reg/${mprage_bn}_head.nii.gz;
		highres_brain=${outdir}/reg/${mprage_bn}_brain.nii.gz;

		meanfunc=${outdir}/meanfunc.nii.gz;
	
		nihpd_brain=${basedir}/template/t1w_brain_2mm.nii.gz;
		nihpd_head=${basedir}/template/t1w_2mm.nii.gz;	
		standard_mask=${FSLDIR}/data/standard/MNI152_T1_2mm_brain_mask_dil;


		##############################################################

		echo "functional to highres";

		# flirt meanfunc to highres
		${FSLDIR}/bin/flirt -ref ${highres_brain} -in ${meanfunc} -out ${outdir}/reg/meanfunc2highres -omat ${outdir}/reg/meanfunc2highres.mat -cost corratio -dof 6 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear;

		# Invert, highres to meanfunc
		${FSLDIR}/bin/convert_xfm -inverse -omat ${outdir}/reg/highres2meanfunc.mat ${outdir}/reg/meanfunc2highres.mat;

		echo "highres to standard, affine";

		# flirt highres to standard
		${FSLDIR}/bin/flirt -ref ${nihpd_brain} -in ${highres_brain} -out ${outdir}/reg/highres2standard -omat ${outdir}/reg/highres2standard.mat -cost corratio -dof 12 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear;
	
		# inverse, standard to highres	
		${FSLDIR}/bin/convert_xfm -inverse -omat ${outdir}/reg/standard2highres.mat ${outdir}/reg/highres2standard.mat

		# meanfunc to standard
		${FSLDIR}/bin/convert_xfm -omat ${outdir}/reg/meanfunc2standard.mat -concat ${outdir}/reg/highres2standard.mat ${outdir}/reg/meanfunc2highres.mat

		# inv, standard to meanfunc
		${FSLDIR}/bin/convert_xfm -omat ${outdir}/reg/standard2meanfunc.mat -inverse ${outdir}/reg/meanfunc2standard.mat

		echo "highres to standard, non-linear";

		# fnirt highres to standard			
		${FSLDIR}/bin/fnirt --in=${highres_head} --aff=${outdir}/reg/highres2standard.mat --cout=${outdir}/reg/highres2standard_warp --iout=${outdir}/reg/highres2standard_fnirt --jout=${outdir}/reg/highres2standard_jac --config=T1_2_MNI152_2mm --ref=${nihpd_head} --refmask=${standard_mask} --warpres=10,10,10;

		#Apply transformations to meanfunc
		${FSLDIR}/bin/applywarp --ref=${nihpd_brain} --in=${meanfunc} --out=${outdir}/reg/meanfunc2standard_w --warp=${outdir}/reg/highres2standard_warp --premat=${outdir}/reg/meanfunc2highres.mat;

		#Apply transformations to highres
		${FSLDIR}/bin/applywarp --ref=${nihpd_brain} --in=${highres_brain} --out=${outdir}/reg/highres2standard_w --warp=${outdir}/reg/highres2standard_warp;


	else
		echo "ERROR: file "$1" doesn't exist";
		exit 0;
	fi
fi
