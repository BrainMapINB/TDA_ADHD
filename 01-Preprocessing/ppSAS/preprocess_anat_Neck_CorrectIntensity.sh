#!/bin/bash

if [ $# -lt 3 ];then

	echo "usage: "$0" mprage f_bet outdir";

	exit 0;

else
	# Define mprage
	if [ -e $1 ];then

		########## Directories and CONSTANTS  here  ################
		mprage=$1;
		mprage_bn=$(basename $mprage .nii.gz);
		f_bet=$2;
		outdir=$3      
		reg_dir=${outdir}/reg;   
		#############################################
		mkdir -p ${reg_dir};

		echo "Bet anatomical"
		${FSLDIR}/bin/bet ${mprage} ${reg_dir}/${mprage_bn}_brain -B -f ${f_bet} -g 0;

		${FSLDIR}/bin/fast -t 1 -n 3 -H 0.1 -I 4 -l 20.0 -B -g -o ${reg_dir}/${mprage_bn}_brain ${reg_dir}/${mprage_bn}_brain;
		
		mv ${reg_dir}/${mprage_bn}_brain_restore.nii.gz ${reg_dir}/${mprage_bn}_brain.nii.gz;

	else
		echo "ERROR: file "$1" doesn't exist";
		exit 0;
	fi
fi
