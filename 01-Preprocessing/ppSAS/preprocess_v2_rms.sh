#!/bin/bash

if [ $# -lt 6 ];then

	echo "usage: "$0"  4D_bold_file  dropTRs  AcqOrder  refVol(starts in 0)  f_bet  outdir";
	echo "";
	echo "AcqOrder: 1-sequential up ";
	echo "		2-sequential down"
	echo "		3-interOther";
	echo "		4-interSiemens";

	exit 0;

else
	echo $1;
	if [ -e $1 ];then

	########## Directories and CONSTANTS  here  ################
	inFile=$1;
	inFile_bn=$(basename $inFile .nii.gz);

	basedir=$(dirname $0);


	dropTRs=$2;
	AcqOrder=$3;
	refVol=$4;
	f_bet=$5;
	outdir=$6;

ST_customFile=${outdir}/ST_customFile_even.txt;
	TR=$(fslval ${inFile} pixdim4);
	echo TR= $TR;
	
	#############################################
	mkdir -p ${outdir};

	echo $TR > ${outdir}/TR

	NVols=$(fslnvols $inFile);
	newNVols=`echo "${NVols}-${dropTRs}" | bc`;

	## 1. Drop first volumes
	echo "Dropping first "${dropTRs}" volumes";
	${FSLDIR}/bin/fslroi $inFile ${outdir}/prefiltered_func_data_tcrop $dropTRs $newNVols;

	## examplefunc
	${FSLDIR}/bin/fslroi ${outdir}/prefiltered_func_data_tcrop ${outdir}/${inFile_bn}_example_func ${refVol} 1;

	## 2. Slice timing correction
	echo "Slice timing correction"
	nslices=$(fslval ${inFile} dim3);

	case "$AcqOrder" in
	1) echo "Acquisition order: sequential bottom-up";
	${FSLDIR}/bin/slicetimer -i ${outdir}/prefiltered_func_data_tcrop --out=${outdir}/prefiltered_func_data_st -r ${TR};
	;;
	2) echo "Acquisition order: sequential up-bottom";
	${FSLDIR}/bin/slicetimer -i ${outdir}/prefiltered_func_data_tcrop --out=${outdir}/prefiltered_func_data_st -r ${TR} --down;
	;;
	3) echo "Acquisition order: interleaved-odd";
	${FSLDIR}/bin/slicetimer -i ${outdir}/prefiltered_func_data_tcrop --out=${outdir}/prefiltered_func_data_st -r ${TR} --odd;
	;;
	4) echo "Acquisition order: interleaved-Siemens";
	out=$(( $nslices % 2 ))
		   if [ $out -eq 0 ];then
			echo $nslices "slices, even number, writing customfile";
	
			rm -f ${ST_customFile};
			seq 2 2 $nslices >> ${ST_customFile};
			seq 1 2 $nslices >> ${ST_customFile};
			
			${FSLDIR}/bin/slicetimer -i ${outdir}/prefiltered_func_data_tcrop --out=${outdir}/prefiltered_func_data_st -r ${TR} --ocustom=${ST_customFile};
		   else
			echo $nslices "slices, ODD number, using slicetimer with --odd option";
			${FSLDIR}/bin/slicetimer -i ${outdir}/prefiltered_func_data_tcrop --out=${outdir}/prefiltered_func_data_st -r ${TR} --odd;
		   fi
	;;
	*) echo "ERROR: not valid entry for acquisition order";

	esac

	## 3. Motion correction using mcflirt, ref volume= $refVol
	echo "Motion correction using mcflirt, ref volume=" ${refVol};
	${FSLDIR}/bin/mcflirt -in ${outdir}/prefiltered_func_data_st -out ${outdir}/prefiltered_func_data_mcf -plots -refvol ${refVol} -rmsrel -verbose 0;

	## 4. Bet
	echo "Applying Bet"
	${FSLDIR}/bin/fslmaths ${outdir}/prefiltered_func_data_mcf -Tmean ${outdir}/pmeanfunc;
	${FSLDIR}/bin/bet2 ${outdir}/pmeanfunc ${outdir}/mask -f ${f_bet} -n -m;
	${FSLDIR}/bin/immv ${outdir}/mask_mask ${outdir}/mask;
	${FSLDIR}/bin/fslmaths ${outdir}/prefiltered_func_data_mcf -mas ${outdir}/mask ${outdir}/prefiltered_func_data_bet;

	## FILTERED_FUNC_DATA; scale global mean 10000
	${FSLDIR}/bin/fslmaths ${outdir}/prefiltered_func_data_bet -ing 10000 ${outdir}/filtered_func_data;
	${FSLDIR}/bin/fslmaths ${outdir}/filtered_func_data -Tmean ${outdir}/meanfunc;

	# Remove prefiltered data
	rm -rf ${outdir}/prefiltered_func_data_*.nii* ${outdir}/pmeanfunc.*;
	rm -rf ${outdir}/prefiltered_func_data_mcf.mat;


	else
		echo "ERROR: file "$1" doesn't exist";
		exit 0;
	fi

fi
