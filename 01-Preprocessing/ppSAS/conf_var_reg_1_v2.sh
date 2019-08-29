#!/bin/bash

if [ $# -eq 0 ];then

	echo "usage: "$0" mprage outdir";

	exit 0;

else

	########## Directories and CONSTANTS  here  ################
	mprage=$1;
	mprage_bn=$(basename $1 .nii.gz);
	outdir=$2;

	basedir=$(dirname $0);
	Conf_dir=${outdir}/conf;
	stats_dir=${outdir}/stats;

	#############################################

	mkdir -p ${Conf_dir};
	mkdir -p ${stats_dir};	

	rm -f ${Conf_dir}/*_func*

	## Masks for WM and CSF, based on Fast results
	echo "Masking WM and CSF based on subject's mprage";

	${FSLDIR}/bin/fslmaths ${outdir}/reg/${mprage_bn}_brain_seg_0 -thr 0.99 -bin -ero ${Conf_dir}/sub_CSF;
	${FSLDIR}/bin/fslmaths ${outdir}/reg/${mprage_bn}_brain_seg_2 -thr 0.99 -bin -ero ${Conf_dir}/sub_WM;
	${FSLDIR}/bin/fslmaths ${Conf_dir}/sub_CSF -add ${Conf_dir}/sub_WM ${Conf_dir}/sub_CSFWM;

	# Taking WM and CSF masks to functional space
	for i in $(ls ${Conf_dir}/sub_*);do
		j=$(basename $i .nii.gz);

		${FSLDIR}/bin/flirt -in ${i} -applyxfm -init ${outdir}/reg/highres2meanfunc.mat -out ${Conf_dir}/${j}_func -paddingsize 0.0 -interp nearestneighbour -ref ${outdir}/meanfunc;	
		
	done

	# Extracting Global signal
	${FSLDIR}/bin/fslmeants -i ${outdir}/filtered_func_data -o ${Conf_dir}/Global.txt -m ${outdir}/mask;

	# Extracting CSF signal
	${FSLDIR}/bin/fslmeants -i ${outdir}/filtered_func_data -o ${Conf_dir}/CSF.txt -m ${Conf_dir}/sub_CSF_func;

	# Extracting WM signal
	${FSLDIR}/bin/fslmeants -i ${outdir}/filtered_func_data -o ${Conf_dir}/WM.txt -m ${Conf_dir}/sub_WM_func; 

	# aCompCor
	${FSLDIR}/bin/fslmeants -i ${outdir}/filtered_func_data -o ${Conf_dir}/aCompCor.txt -m ${Conf_dir}/sub_CSFWM_func --eig --order=5; 

fi
