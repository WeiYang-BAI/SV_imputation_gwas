#!/bin/bash

sn=${SLURM_ARRAY_TASK_ID} # from 2 since to skip header

M=M24sl5F; GP=HGSVC2_HPRCY1


IMPU_D=IMPUTATION
VNTR_D=VNTR

TR_ref_F=VNTR/${GP}/${M}_len7_InDel_SV_match_HG38/${GP}.HG38.${M}.VNTR.chrAuto.summary
loci_D=VNTR/${GP}/${M}_len7_InDel_SV_match_HG38/zDense_VNTR_loci_with_TR_allele

loci_id=$(cat ${TR_ref_F}|sed -n "${sn}p"|awk '{print $1}')
c=$(cat ${TR_ref_F}|sed -n "${sn}p"|awk '{print $3}'|sed 's/chr//g')

pop=EUR

aim_D=${VNTR_D}/prepare_GENO/${pop}/chr${c}

mkdir -p ${aim_D}/ale_vid
mkdir -p ${aim_D}/ale_vcf
mkdir -p ${aim_D}/ale_gt
mkdir -p ${aim_D}/ale_cn

# get VNTR TR ID
cat ${loci_D}/${loci_id}.tsv|sed '1d'|awk '{print $18}' > ${aim_D}/ale_vid/${loci_id}.vid

# extract imputed TR by ID, and filter estimated MAF < 0.01
bcftools view -Oz --include ID=@${aim_D}/ale_vid/${loci_id}.vid ${IMPU_D}/${GP}_VNTR/${pop}/zTR/${GP}.${M}.VNTR.${pop}.chr${c}.impu.dose.TR.reN.sort.vcf.gz |\
bcftools view -Oz --exclude 'INFO/MAF<0.01' |\
bcftools annotate -Ov -x FILTER,INFO,FORMAT/GP --output ${aim_D}/ale_vcf/${loci_id}.dos-maf01.vcf

# extract imputed TR by ID, and filter real MAF < 0.01
bcftools view -Oz --include ID=@${aim_D}/ale_vid/${loci_id}.vid ${IMPU_D}/${GP}_VNTR/${pop}/zTR/${GP}.${M}.VNTR.${pop}.chr${c}.impu.dose.TR.reN.sort.vcf.gz |\
bcftools view -Oz --min-af 0.01:minor |\
bcftools annotate -Ov -x FILTER,INFO,FORMAT/GP --output ${aim_D}/ale_vcf/${loci_id}.int-maf01.vcf

# deal genotype, note that this mtp is for MAF only, i.e. dos-maf and int-maf, nothing to to with genotype
for mtp in dos int; do

	# check num of TR in the loci, only loci with >= 2 allele were retained for analyais
	TR_n=$(cat ${aim_D}/ale_vcf/${loci_id}.${mtp}-maf01.vcf|grep -ac ^c)

	if [ ${TR_n} -gt 1 ]; then

		# get genotype int:dos for each loci
		cat ${aim_D}/ale_vcf/${loci_id}.${mtp}-maf01.vcf | grep -v '^##' |\
		sed 's/0|0/0/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/2/g' | awk -v OFS='-' '{$1="";$2="";$4="";$5="";$6="";$7="";$8="";$9=""; print $0}' |\
		sed 's/--*/-/g' | sed 's/^-//g' > ${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt
		sed -i 's/-/ /g' ${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt

		# get cn for each ale, ## NOTE THAT THE tr_id CANNOT TOO LONG WHEN USE grep, SO HERE USE -f flag !!!##
		for tr_id in $(cat ${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt|grep -a ^c|awk '{print $1}');do
			echo ${tr_id} > ${aim_D}/ale_cn/${loci_id}.${mtp}-maf01.vid-cn.tmp_trid
			tr_id_cn=$(cat ${loci_D}/${loci_id}.tsv|grep -w -f ${aim_D}/ale_cn/${loci_id}.${mtp}-maf01.vid-cn.tmp_trid|awk '{print $18,$3,$9}'|sed 's/INS //'|sed 's/DEL /-/')
			echo -e ${tr_id_cn} >> ${aim_D}/ale_cn/${loci_id}.${mtp}-maf01.vid-cn
		done

		
		# get gt * cn, step 1, int and dos in one file
		python get_VNTR_loci-geno-cn_step-1.py \
		${aim_D}/ale_cn/${loci_id}.${mtp}-maf01.vid-cn \
		${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt \
		${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt-cn

		# get gt * cn, step 2, int and dos in two files
		python get_VNTR_loci-geno-cn_step-2.py \
		${loci_id} \
		${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt-cn \
		${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt-cn-int \
		${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt-cn-dos

		# trans gt format
		for gtp in dos int;do

			# remove ' ' at the tail of each line
			sed -i 's/ $//g' ${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt-cn-${gtp}
			sed -i 's/ $//g' ${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt-cn-${gtp}

			# transpos
			python get_VNTR_loci-geno-cn_step-3.py \
			${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt-cn-${gtp} \
			${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt-cn-${gtp}.T

			cat ${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt-cn-${gtp}.T | awk '{print $2}' > ${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt-cn-${gtp}.T.tmp

		done
	fi
done


for pop in EUR AFR AMR EAS SAS; do
for c in $(seq 1 22);do
bd=06.VNTR/prepare_GENO/${pop}/chr${c}/ale_gt; cd ${bd}
for mtp in dos int; do
for gtp in dos int; do
cat z.${pop}.chr${c}.N2.${mtp}-maf01.vid-gt-cn-${gtp}.T|head -n1|tr ' ' '\n'|sed '1d' > z.${pop}.chr${c}.N2.${mtp}-maf01.vid-gt-cn-${gtp}.vntr_id
done
done
done
done

for c in $(seq 1 22);do
	for mtp in dos int; do
		for gtp in dos int; do
			PT=$(echo ${phe} | cut -b 1)
			if [ $PT == "D" ]; then
				for pop in EUR ; do
					regression_D=${VNTR_D}/regression_unRel05/${pop}/${mtp}_maf01_${gtp}_geno_rsq00/${phe}
					mkdir -p ${regression_D}

					Rscript --slave --no-save 2_regression_unRel.R \
					${VNTR_D}/prepare_GENO/${pop}/chr${c}/ale_gt/z.${pop}.chr${c}.N2.${mtp}-maf01.vid-gt-cn-${gtp}.vntr_id \
					${VNTR_D}/prepare_GENO/${pop}/chr${c}/ale_gt/z.${pop}.chr${c}.N2.${mtp}-maf01.vid-gt-cn-${gtp}.T \
					${phenoDir}/Main_74/${phe}.pheno \
					${phenoDir}/covar/UKB_All_covariates.covar \
					${regression_D}/${phe}.${pop}.chr${c}.${mtp}.MAF01.${gtp}.Rsq00.unRel.lr \
					${GENO_D}/infos/zUnRel_ID/UnRel.0.05.Merged.${pop}.maf0.01.hwe-6.chrAuto.grm.id \
					bina
				done
			else
				for pop in EUR AFR AMR EAS SAS; do
					regression_D=${VNTR_D}/regression_unRel05/${pop}/${mtp}_maf01_${gtp}_geno_rsq00/${phe}
					mkdir -p ${regression_D}

					Rscript --slave --no-save 2_regression_unRel.R \
					${VNTR_D}/prepare_GENO/${pop}/chr${c}/ale_gt/z.${pop}.chr${c}.N2.${mtp}-maf01.vid-gt-cn-${gtp}.vntr_id \
					${VNTR_D}/prepare_GENO/${pop}/chr${c}/ale_gt/z.${pop}.chr${c}.N2.${mtp}-maf01.vid-gt-cn-${gtp}.T \
					${phenoDir}/Main_74/${phe}.pheno \
					${phenoDir}/covar/UKB_All_covariates.covar \
					${regression_D}/${phe}.${pop}.chr${c}.${mtp}.MAF01.${gtp}.Rsq00.unRel.lr \
					${GENO_D}/infos/zUnRel_ID/UnRel.0.05.Merged.${pop}.maf0.01.hwe-6.chrAuto.grm.id \
					conti
				done
			fi
		done
	done
done


for c in $(seq 1 22);do
	if [ $sn -gt 2088 ]; then
		Rscript --slave --no-save 2_regression_unRel.R \
		${VNTR_D}/prepare_GENO/${pop}/chr${c}/ale_gt/z.${pop}.chr${c}.N2.${mtp}-maf01.vid-gt-cn-${gtp}.vntr_id \
		${VNTR_D}/prepare_GENO/${pop}/chr${c}/ale_gt/z.${pop}.chr${c}.N2.${mtp}-maf01.vid-gt-cn-${gtp}.T \
		${phe_F} \
		${phenoDir}/covar/UKB_All_covariates.covar \
		${regression_D}/${phe}.${pop}.chr${c}.${mtp}.MAF01.${gtp}.Rsq00.unRel.lr \
		${GENO_D}/infos/zUnRel_ID/UnRel.0.05.Merged.${pop}.maf0.01.hwe-6.chrAuto.grm.id \
		conti
	else
		Rscript --slave --no-save 2_regression_unRel.R \
		${VNTR_D}/prepare_GENO/${pop}/chr${c}/ale_gt/z.${pop}.chr${c}.N2.${mtp}-maf01.vid-gt-cn-${gtp}.vntr_id \
		${VNTR_D}/prepare_GENO/${pop}/chr${c}/ale_gt/z.${pop}.chr${c}.N2.${mtp}-maf01.vid-gt-cn-${gtp}.T \
		${phe_F} \
		${phenoDir}/covar/UKB_All_covariates.covar \
		${regression_D}/${phe}.${pop}.chr${c}.${mtp}.MAF01.${gtp}.Rsq00.unRel.lr \
		${GENO_D}/infos/zUnRel_ID/UnRel.0.05.Merged.${pop}.maf0.01.hwe-6.chrAuto.grm.id \
		bina
	fi
done


xmat_D=GWAS/cojo_select_bench/${pop}/${phe}/covar
residual_dir=VNTR/bench_condi_cojo_residual/${pop}
mkdir -p ${residual_dir}

if [ $sn -gt 2088 ]; then boc=conti; else boc=bina; fi


for vt in SGV SV;do


if [ -s "${xmat_D}/chrAuto.${vt}.cojo-slct.012.xmat.F" ]; then

	python3 4_generate_pheno_geno.py \
	${phe_F} \
	${xmat_D}/chrAuto.${vt}.cojo-slct.012.xmat.F \
	${residual_dir}/${phe}.${pop}.cojo-slct.chrAuto.${vt}.xmat.pheno

	sed -i '/NoPheno/d' ${residual_dir}/${phe}.${pop}.cojo-slct.chrAuto.${vt}.xmat.pheno
	sed -i 's/[0-9]*:/chrN/g' ${residual_dir}/${phe}.${pop}.cojo-slct.chrAuto.${vt}.xmat.pheno


	## R-mlr with cojo-slct-var.chrAuto.geno.pheno
	conda deactivate
	conda activate R-3.6.3

	# set R-mlr model
	mlrModel=$(cat ${residual_dir}/${phe}.${pop}.cojo-slct.chrAuto.${vt}.xmat.pheno | head -n 1 | sed 's/IID pheno /pheno~/' | sed 's/ /+/g')

	# run R-mlr
	Rscript --slave --no-save 4_generate_cojo-var_residual.R \
	${residual_dir}/${phe}.${pop}.cojo-slct.chrAuto.${vt}.xmat.pheno \
	${residual_dir}/${phe}.${pop}.cojo-slct.chrAuto.${vt}.xmat.pheno.mlr.residual.raw \
	${mlrModel} \
	${boc}

	# preprocess R-mlr residual
	cat ${residual_dir}/${phe}.${pop}.cojo-slct.chrAuto.${vt}.xmat.pheno | awk '{print $1,$1}' >\
	${residual_dir}/${phe}.${pop}.cojo-slct.chrAuto.${vt}.xmat.pheno.tmp.sample

	cat ${residual_dir}/${phe}.${pop}.cojo-slct.chrAuto.${vt}.xmat.pheno.mlr.residual.raw | awk -F',' '{print $2}' >\
	${residual_dir}/${phe}.${pop}.cojo-slct.chrAuto.${vt}.xmat.pheno.tmp.residual

	paste ${residual_dir}/${phe}.${pop}.cojo-slct.chrAuto.${vt}.xmat.pheno.tmp.sample ${residual_dir}/${phe}.${pop}.cojo-slct.chrAuto.${vt}.xmat.pheno.tmp.residual >\
	${residual_dir}/${phe}.${pop}.cojo-slct.chrAuto.${vt}.xmat.pheno.mlr.residual && \
	rm ${residual_dir}/${phe}.${pop}.cojo-slct.chrAuto.${vt}.xmat.pheno.tmp.sample ${residual_dir}/${phe}.${pop}.cojo-slct.chrAuto.${vt}.xmat.pheno.tmp.residual

	sed -i 's/\t/ /g' ${residual_dir}/${phe}.${pop}.cojo-slct.chrAuto.${vt}.xmat.pheno.mlr.residual
	sed -i '1d' ${residual_dir}/${phe}.${pop}.cojo-slct.chrAuto.${vt}.xmat.pheno.mlr.residual


fi

done


boc=conti

if [ -s "${residual_dir}/${phe}.${pop}.cojo-slct.chrAuto.SGV.xmat.pheno.mlr.residual" ]; then

for c in $(seq 1 22);do
	Rscript --slave --no-save 2_regression_unRel.R \
	${VNTR_D}/prepare_GENO/${pop}/chr${c}/ale_gt/z.${pop}.chr${c}.N2.${mtp}-maf01.vid-gt-cn-${gtp}.vntr_id \
	${VNTR_D}/prepare_GENO/${pop}/chr${c}/ale_gt/z.${pop}.chr${c}.N2.${mtp}-maf01.vid-gt-cn-${gtp}.T \
	${residual_dir}/${phe}.${pop}.cojo-slct.chrAuto.SGV.xmat.pheno.mlr.residual \
	${phenoDir}/covar/UKB_All_covariates.covar \
	${VNTR_D}/condi_bench_SGVs_unRel05/${pop}/${phe}/${phe}.${pop}.chr${c}.${mtp}.MAF01.${gtp}.Rsq00.unRel.lr \
	${GENO_D}/infos/zUnRel_ID/UnRel.0.05.Merged.${pop}.maf0.01.hwe-6.chrAuto.grm.id \
	${boc}
done

fi



