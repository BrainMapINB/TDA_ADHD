#!/bin/bash

if [ $# -eq 0 ];then

	echo "usage: "$0" 4D_bold mprage";

	exit 0;

else

#####################################################
	bold=$1
	pre_mprage=$2

	basedir=$(dirname ${0});
	
	bold_bn=$(basename $1 .nii.gz);
	bold_bn=$(basename $bold_bn .nii);

	outdir=$(dirname $bold)/${bold_bn}.rsfMRIv2_nihpd;

	outputbasename=ppBoldv2; #output name for final data in standard space

######################################################
mkdir -p ${outdir};
echo "Inputs: " >> ${outdir}/inputs.txt;
echo " " >> ${outdir}/inputs.txt;
echo "Bold: " ${bold} >> ${outdir}/inputs.txt;
echo "MPRage: " ${pmprage} >> ${outdir}/inputs.txt;


	echo "##################################";
	echo "Preprocessing 4D_bold";
	echo "##################################";
	#use: ${basedir}/preprocess_v2.sh bold dropTRs AcqOrd(1 seq-up,2 seq-down,3 int, 4 int-Siemens) refVol f_bet outdir;
	${basedir}/preprocess_v2_rms.sh $bold 4 3 0 0.35 $outdir;

	echo "##################################";
	echo "Preprocessing mprage";
	echo "##################################";
	${FSLDIR}/bin/fslreorient2std $pre_mprage ${outdir}/Structural.nii.gz;
	mprage=${outdir}/Structural.nii.gz;
	#use: ${basedir}/preprocess_anat_Neck_CorrectIntensity.sh mprage f_beat outdir;
	${basedir}/preprocess_anat_Neck_CorrectIntensity.sh $mprage 0.35 $outdir;

	echo "##################################";
	echo "Registration";
	echo "##################################";
	#use: ${basedir}/regist_anat_func.sh mprage outdir
	${basedir}/regist_anat_func_nihpd.sh $mprage $outdir 0;
	${basedir}/Check_registration_nihpd.sh $outdir $mprage;

	echo "################################################";
	echo "Relative Displacement RMS estimation";
	echo "################################################";
	${basedir}/RMSdisplacement.sh $outdir 0.25;		
	${basedir}/Movement2png_v2_rms.sh $outdir;

	echo "############################################";
	echo "Extracting confounding variables from data";
	echo "############################################";
	${basedir}/conf_var_reg_1_v2.sh $mprage $outdir;
	${basedir}/conf_var_reg_2_v2_rms.sh $outdir $basedir;	 #<<<<<<< With and without GSR 

	echo "################################################";
	echo "Apply Transformations and Temporal Filtering";
	echo "################################################";
	#use: ${basedir}/res4d_2MNI_TF.sh basedir outdir outputbasename HPFc LPFc KeepIntermMNIfiles(0|1)
	#${basedir}/res4d_2MNI_TF_nihpd.sh $basedir $outdir $outputbasename 0.01 0.08 0;

fi




