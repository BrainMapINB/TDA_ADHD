#!/bin/bash

if [ $# -lt 2 ];then

	echo "usage: "$0" outdir mprage";

	exit 0;

else
	#echo $1 $2;


		########## Directories and CONSTANTS  here  ################
		outdir=$1;
		bold_bn=$(basename $outdir .rsfMRIv2_nihpd);
		mprage=$2
		mprage_bn=$(basename $2 .nii.gz);

		##############################################################

	if [ -e ${outdir}/reg/${mprage_bn}_head.nii.gz ];then

		########## Vars  here ######################################

		highres_head=${outdir}/reg/${mprage_bn}_head.nii.gz;
		highres_brain=${outdir}/reg/${mprage_bn}_brain.nii.gz;
	
		nihpd_brain=${basedir}/template/t1w_brain_2mm.nii.gz;
		nihpd_head=${basedir}/template/t1w_2mm.nii.gz;	
		
		standard_mask=${FSLDIR}/data/standard/MNI152_T1_2mm_brain_mask_dil;

		font=/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf

		##############################################################

		echo "Creating checkreg figures"

		${FSLDIR}/bin/slicer ${outdir}/reg/meanfunc2highres.nii.gz ${highres_brain} -a ${outdir}/reg/func_highres.png;
		${FSLDIR}/bin/slicer -e 0.1 ${outdir}/reg/highres2standard_w.nii.gz ${nihpd_brain} -a ${outdir}/reg/highres_standard.png;
		${FSLDIR}/bin/slicer -e 0.1 ${outdir}/reg/meanfunc2standard_w ${nihpd_brain} -a ${outdir}/reg/func_standard.png;

		montage -quality 60 -font ${font} -fill white -label '%f' ${outdir}/reg/func_highres.png ${outdir}/reg/highres_standard.png ${outdir}/reg/func_standard.png -tile 1x3 -background '#000000' -geometry '480' ${outdir}/reg/checkreg_${bold_bn}.jpg
		
		
		rm ${outdir}/reg/*.png
	
	else
		echo "file "$2" doesn't exist";
		exit 0;
	fi
fi
