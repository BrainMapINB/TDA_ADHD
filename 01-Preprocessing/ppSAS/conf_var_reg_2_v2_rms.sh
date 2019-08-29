#!/bin/bash

if [ $# -lt 2 ];then

	echo "usage: "$0" outdir basedir";

	exit 0;

else

	########## Directories and CONSTANTS  here  ################
	outdir=$1;
	basedir=$2;

	Conf_dir=${outdir}/conf;
	stats_dir=${outdir}/stats;
	QMov_dir=${outdir}/QMov;
	#############################################
	mkdir -p ${stats_dir};
	rm -f ${outdir}/design.mat;

	ntpoints=$(wc -l < ${Conf_dir}/WM.txt);

	#Paste columns
	echo "Concatenate columns";
	paste -d " " ${Conf_dir}/Global.txt ${Conf_dir}/CSF.txt ${Conf_dir}/WM.txt ${outdir}/prefiltered_func_data_mcf.par > ${Conf_dir}/design_prev.mat;

		# No GSR
		paste -d " " ${Conf_dir}/CSF.txt ${Conf_dir}/WM.txt ${outdir}/prefiltered_func_data_mcf.par > ${Conf_dir}/design_prev_woGSR.mat; 

	# GENERATE derivatives and quadratic trends
	Rscript --vanilla ${basedir}/dp2.R ${Conf_dir}/design_prev.mat ${Conf_dir}/dp2_design_prev.mat
	Rscript --vanilla ${basedir}/dp2.R ${Conf_dir}/design_prev_woGSR.mat ${Conf_dir}/dp2_design_prev_woGSR.mat


	#Create design.mat
	npeaks=$(cat ${QMov_dir}/NVols_RMSpeaks.txt);
	if [ $npeaks -gt 0 ];then

		# GENERATE spike columns
		Rscript --vanilla ${basedir}/colpeaks.R ${QMov_dir}/RMS_peaks.txt ${QMov_dir}/RMS_confvars.txt

		nwaves=`echo "41+$npeaks" | bc`;
		echo "/NumWaves $nwaves" > ${Conf_dir}/design.mat;
		echo "/NumPoints $ntpoints" >> ${Conf_dir}/design.mat;
		echo "" >> ${Conf_dir}/design.mat;
		paste -d "\t" ${Conf_dir}/dp2_design_prev.mat ${Conf_dir}/aCompCor.txt ${QMov_dir}/RMS_confvars.txt >> ${Conf_dir}/design.mat;
		#noGSR
		echo "/NumWaves 37+$npeaks" > ${Conf_dir}/design_woGSR.mat;
		echo "/NumPoints $ntpoints" >> ${Conf_dir}/design_woGSR.mat;
		echo "" >> ${Conf_dir}/design_woGSR.mat;
		paste -d "\t" ${Conf_dir}/dp2_design_prev_woGSR.mat ${Conf_dir}/aCompCor.txt ${QMov_dir}/RMS_confvars.txt >> ${Conf_dir}/design_woGSR.mat;	
	
	else
		echo "/NumWaves 41" > ${Conf_dir}/design.mat;
		echo "/NumPoints $ntpoints" >> ${Conf_dir}/design.mat;
		echo "" >> ${Conf_dir}/design.mat;
		paste -d "\t" ${Conf_dir}/dp2_design_prev.mat ${Conf_dir}/aCompCor.txt >> ${Conf_dir}/design.mat;
		#noGSR
		echo "/NumWaves 37" > ${Conf_dir}/design_woGSR.mat;
		echo "/NumPoints $ntpoints" >> ${Conf_dir}/design_woGSR.mat;
		echo "" >> ${Conf_dir}/design_woGSR.mat;
		paste -d "\t" ${Conf_dir}/dp2_design_prev_woGSR.mat ${Conf_dir}/aCompCor.txt >> ${Conf_dir}/design_woGSR.mat;	

	fi

echo "fsl_glm";

	# fsl_glm
	${FSLDIR}/bin/fsl_glm -i ${outdir}/filtered_func_data.nii.gz -m ${outdir}/mask.nii.gz -d ${Conf_dir}/design.mat --demean --out_res=${stats_dir}/res4d.nii.gz;
	## correct for fsl_glm bug that sets pixdim4 to 1
	fslmodhd ${stats_dir}/res4d.nii.gz pixdim4 $(cat ${outdir}/TR)


	mean_bold=$(${FSLDIR}/bin/fslstats ${outdir}/filtered_func_data -k ${outdir}/mask -M);
	echo $mean_bold >> ${stats_dir}/mean_filteredfunc;
	##${FSLDIR}/bin/fslmaths ${stats_dir}/res4d -add ${mean_bold} ${stats_dir}/res4d_pmean;

	# fsl_glm_woGSR
	${FSLDIR}/bin/fsl_glm -i ${outdir}/filtered_func_data.nii.gz -m ${outdir}/mask.nii.gz -d ${Conf_dir}/design_woGSR.mat --demean --out_res=${stats_dir}/res4d_woGSR.nii.gz;
	fslmodhd ${stats_dir}/res4d_woGSR.nii.gz pixdim4 $(cat ${outdir}/TR)

	##${FSLDIR}/bin/fslmaths ${stats_dir}/res4d_woGSR -add ${mean_bold} ${stats_dir}/res4d_woGSR_pmean;

	# ----- remove temp files

##	rm ${Conf_dir}/*design_prev*.mat

fi
