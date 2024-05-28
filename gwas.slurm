#!/bin/bash

tr=6; sn=${SLURM_ARRAY_TASK_ID}


if [ $sn -gt 2088 ]; then cob=C; REPLACE_model=fastGWA-mlm; else cob=B; REPLACE_model=fastGWA-mlm-binary; fi

filePath_F=both.final.filePath; geno_D=GENOTYPE; GRM_D=GRM; GWAS_D=GWAS; ori_covar_F=UKB_All_covariates.covar

phe=$(cat ${filePath_F} | sed -n "${sn}p" | awk -F'/' '{print $8}' | awk -F'.pheno' '{print $1}'); pheno_F=$(cat ${filePath_F} | sed -n "${sn}p")

for pop in EUR ; do

	geno_mbfile_F=${geno_D}/infos/zFILEPATH/Merged.${pop}.maf0.01.hwe-6.chrAuto.bfile
	
	# specify sparse GRM
	if [ $pop == 'EUR' ]; then spGRM_F=${GRM_D}/GRM_EUR/Sparse_UKB_GRM_HM3_QC_R2_09_MAF01; else spGRM_F=${GRM_D}/GRM_POP4/Sparse.Merged.${pop}.maf0.01.hwe-6.chrAuto; fi

	# run fastGWA
	mkdir -p ${GWAS_D}/fastGWA_ori_bench/${pop}/${phe}
	
	assoc_ori=${GWAS_D}/fastGWA_ori_bench/${pop}/${phe}/${phe}.Merged.${pop}.maf0.01.hwe-6.chrAuto

	${gcta64} --mbfile ${geno_mbfile_F} --qcovar ${ori_covar_F} --grm-sparse ${spGRM_F} --${REPLACE_model} --pheno ${pheno_F} --maf 0.01 --thread-num ${tr} --out ${assoc_ori}

	cojo_fmt_TWO=${assoc_ori}.fastGWA.cojoFMT.TWO
	cojo_fmt_SV=${assoc_ori}.fastGWA.cojoFMT.SV
	cojo_fmt_SGV=${assoc_ori}.fastGWA.cojoFMT.SGV

	# ori GWAS bibary header: CHR	SNP POS A1  A2  N   AF1 T   SE_T	P_noSPA BETA	SE  P   CONVERGE
	# ori GWAS contiuos header: CHR  SNP POS A1  A2  N   AF1 BETA	SE  P

	echo -e "SNP\tA1\tA2\tfreq\tb\tse\tp\tN" > ${cojo_fmt_TWO}
	echo -e "SNP\tA1\tA2\tfreq\tb\tse\tp\tN" > ${cojo_fmt_SV}
	echo -e "SNP\tA1\tA2\tfreq\tb\tse\tp\tN" > ${cojo_fmt_SGV}

	if [ $cob == 'B' ];then
		cat ${assoc_ori}.fastGWA | sed '1d' |awk '{print $2,$4,$5,$7,$11,$12,$13,$6}' >> ${cojo_fmt_TWO}
		cat ${assoc_ori}.fastGWA | sed '1d' |grep HPRC|awk '{print $2,$4,$5,$7,$11,$12,$13,$6}' >> ${cojo_fmt_SV}
		cat ${assoc_ori}.fastGWA | sed '1d' |grep -v HPRC|awk '{print $2,$4,$5,$7,$11,$12,$13,$6}' >> ${cojo_fmt_SGV}
	else
		cat ${assoc_ori}.fastGWA | sed '1d' |awk '{print $2,$4,$5,$7,$8,$9,$10,$6}' >> ${cojo_fmt_TWO}
		cat ${assoc_ori}.fastGWA | sed '1d' |grep HPRC|awk '{print $2,$4,$5,$7,$8,$9,$10,$6}' >> ${cojo_fmt_SV}
		cat ${assoc_ori}.fastGWA | sed '1d' |grep -v HPRC|awk '{print $2,$4,$5,$7,$8,$9,$10,$6}' >> ${cojo_fmt_SGV}
	fi

	# get cojo select var
	mkdir -p ${GWAS_D}/cojo_select_bench/${pop}/${phe}/association
	mkdir -p ${GWAS_D}/cojo_select_bench/${pop}/${phe}/covar

	for vTYPE in TWO SV SGV; do 
		cojo_sclt_varid=${GWAS_D}/cojo_select_bench/${pop}/${phe}/association/chrAuto.${vTYPE}.cojo-slct.vid
		rm -rf ${cojo_sclt_varid}

		# COJO chr by chr
		for c in $(seq 1 22); do
			ref_LD_bfile=${geno_D}/LD_ref/${pop}/Merged.${pop}.chr${c}.maf0.01.hwe-6
			cojo_sclt_var=${GWAS_D}/cojo_select_bench/${pop}/${phe}/association/chr${c}.${vTYPE}.cojo-slct

			${gcta64} --bfile ${ref_LD_bfile} --cojo-file ${assoc_ori}.fastGWA.cojoFMT.${vTYPE} --chr ${c} --cojo-slct --cojo-p 5e-8 --thread-num ${tr} --out ${cojo_sclt_var}

			if [ -s "${cojo_sclt_var}.jma.cojo" ]; then cat ${cojo_sclt_var}.jma.cojo | sed '1d' | awk '{print $2}'  >> ${cojo_sclt_varid}; fi
		done

	# condi if there are significant cojo var 
		if [ -s "${cojo_sclt_varid}" ]; then

			# prepare condi files
			cojo_slct_012=${GWAS_D}/cojo_select_bench/${pop}/${phe}/covar/chrAuto.${vTYPE}.cojo-slct.012

			# get cojo-slct var genotype 0 1 2
			${gcta64} --mbfile ${geno_mbfile_F} --extract ${cojo_sclt_varid} --recode-nomiss --thread-num ${tr} --out ${cojo_slct_012}
			gunzip ${cojo_slct_012}.xmat.gz; cat ${cojo_slct_012}.xmat | awk '{$1="";print $0}' | sed 's/^ //' | sed '2d' > ${cojo_slct_012}.xmat.F && rm ${cojo_slct_012}.xmat

			## make new covar with cojo-slct var genotype 0 1 2
			generate_new_covar_R_script=1_make_new_covar.R
			Rscript --slave --no-save ${generate_new_covar_R_script} ${cojo_slct_012}.xmat.F ${ori_covar_F} ${cojo_slct_012}.covar

			# rerun fastGWA
			mkdir -p ${GWAS_D}/fastGWA_condi_bench/${pop}/${phe}
			assoc_condi=${GWAS_D}/fastGWA_condi_bench/${pop}/${phe}/${phe}.Merged.${pop}.maf0.01.hwe-6.chrAuto.condi_${vTYPE}.cojo-slct
			${gcta64} --mbfile ${geno_mbfile_F} --qcovar ${cojo_slct_012}.covar --grm-sparse ${spGRM_F} --${REPLACE_model} --pheno ${pheno_F} --maf 0.01 --thread-num ${tr} --out ${assoc_condi}
		fi
	done
done





pop=EUR; clumpDir=GWAS/clump_bench

# 2p(1-p)b^2 / var(y)
phe=$(cat ${phenoDir}/infos/Continuous.final.fieldID.652|sed -n "${sn}p" | awk '{print $1}')

pheno_F_1=${phenoDir}/Continuous_final_652/${phe}.pheno
pheno_F_2=${phenoDir}/Continuous_final_652/${phe}.residual.pheno

# get sig ori fastgwa, no header
cat ${fastGWA_ori}/${pop}/${phe}/${phe}.Merged.${pop}.maf0.01.hwe-6.chrAuto.fastGWA|awk '{if($10<5e-8) print $0}' > ${fastGWA_ori}/${pop}/${phe}.Merged.${pop}.maf0.01.hwe-6.chrAuto.fastGWA.5e-8

for c in $(seq 1 22);do

	# prepare clump file
	echo -e "SNP P" > ${fastGWA_ori}/${pop}/${phe}/zChr/${phe}.Merged.${pop}.maf0.01.hwe-6.chr${c}.fastGWA.SV.P.5e-8
	cat ${fastGWA_ori}/${pop}/${phe}/${phe}.Merged.${pop}.maf0.01.hwe-6.chrAuto.fastGWA|grep -w ^${c}|awk '{if($10<5e-8) print $2,$10}'|grep SV >> ${fastGWA_ori}/${pop}/${phe}/zChr/${phe}.Merged.${pop}.maf0.01.hwe-6.chr${c}.fastGWA.SV.P.5e-8

	echo -e "SNP P" > ${fastGWA_ori}/${pop}/${phe}/zChr/${phe}.Merged.${pop}.maf0.01.hwe-6.chr${c}.fastGWA.SGV.P.5e-8
	cat ${fastGWA_ori}/${pop}/${phe}/${phe}.Merged.${pop}.maf0.01.hwe-6.chrAuto.fastGWA|grep -w ^${c}|awk '{if($10<5e-8) print $2,$10}'|grep -v SV >> ${fastGWA_ori}/${pop}/${phe}/zChr/${phe}.Merged.${pop}.maf0.01.hwe-6.chr${c}.fastGWA.SGV.P.5e-8

	echo -e "SNP P" > ${fastGWA_ori}/${pop}/${phe}/zChr/${phe}.Merged.${pop}.maf0.01.hwe-6.chr${c}.fastGWA.TWO.P.5e-8
	cat ${fastGWA_ori}/${pop}/${phe}/${phe}.Merged.${pop}.maf0.01.hwe-6.chrAuto.fastGWA|grep -w ^${c}|awk '{if($10<5e-8) print $2,$10}' >> ${fastGWA_ori}/${pop}/${phe}/zChr/${phe}.Merged.${pop}.maf0.01.hwe-6.chr${c}.fastGWA.TWO.P.5e-8


	LD_bfile=LD_ref/${pop}/Merged.${pop}.chr${c}.maf0.01.hwe-6

	for vT in SV SGV TWO; do

		VAR_P_file=${fastGWA_ori}/${pop}/${phe}/zChr/${phe}.Merged.${pop}.maf0.01.hwe-6.chr${c}.fastGWA.${vT}.P.5e-8
		plink --bfile ${LD_bfile} --clump ${VAR_P_file} --clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.01 --clump-kb 500 --out ${clumpDir}/${pop}/${phe}/clump.0.01.chr${c}.${vT}

		if [ -s "${clumpDir}/${pop}/${phe}/clump.0.01.chr${c}.${vT}.clumped" ]; then
			cat ${clumpDir}/${pop}/${phe}/clump.0.01.chr${c}.${vT}.clumped|sed 1d|awk '{if($3!="") print $3}' >> ${clumpDir}/${pop}/${phe}/clump.0.01.chrAuto.${vT}.vid
		fi
	
	done


done

if [ -s "${clumpDir}/${pop}/${phe}/clump.0.01.chrAuto.TWO.vid" ]; then
	cat ${clumpDir}/${pop}/${phe}/clump.0.01.chrAuto.TWO.vid |grep SV > ${clumpDir}/${pop}/${phe}/clump.0.01.chrAuto.TWO_SV.vid
	cat ${clumpDir}/${pop}/${phe}/clump.0.01.chrAuto.TWO.vid |grep -v SV > ${clumpDir}/${pop}/${phe}/clump.0.01.chrAuto.TWO_SGV.vid
fi


for vT in SV SGV TWO TWO_SV TWO_SGV;do
	
	if [ -s "${clumpDir}/${pop}/${phe}/clump.0.01.chrAuto.${vT}.vid" ]; then

		rm -rf ${indi_hsq_D}/${pop}/${phe}/clump.0.01.${vT}.af.beta
		
		for vid in $(cat ${clumpDir}/${pop}/${phe}/clump.0.01.chrAuto.${vT}.vid);do

			cat ${fastGWA_ori}/${pop}/${phe}.Merged.${pop}.maf0.01.hwe-6.chrAuto.fastGWA.5e-8|grep -w ${vid} | awk -v OFS='\t' '{print $7,$8}' >> ${indi_hsq_D}/${pop}/${phe}/clump.0.01.${vT}.af.beta

		done

		Rscript --slave --no-save 11_calculate_individual_hsq.R \
		${phe} ${vT} ${pheno_F_1} ${pheno_F_2} ${indi_hsq_D}/${pop}/${phe}/clump.0.01.${vT}.af.beta ${indi_hsq_D}/${pop}/${phe}/clump.0.01.${vT}.af.beta.indi-hsq

	fi
done


