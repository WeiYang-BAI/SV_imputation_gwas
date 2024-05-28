#!/bin/bash


pN=${SLURM_ARRAY_TASK_ID}

geno_D=02.GENOTYPE
GRM_D=03.GRM; GRM_spe_D=GRM_EUR/GRM_spe_phe_chunks_cont_652
GREML_D=04.GREML


pop=EUR

# run 1-459 with one chunk#
phe=$(cat GRM_EUR/GRM_spe_phe_chunks_cont_652/infos/sample.size.sort.fid|sed -n "${pN}p"|awk '{print $2}')
sample_size=$(cat GRM_EUR/GRM_spe_phe_chunks_cont_652/infos/sample.size.sort.fid|sed -n "${pN}p"|awk '{print $1}')

####### sealed

spe_baseDir=final

prSNV=$(cat ${spe_baseDir}/HE_bench/EUR/compare_to_GREML/GREML_prior_NA0.1_neg0_gpos1.hsq|grep -w "${phe}"|awk '{print $2}')
prSV=$(cat ${spe_baseDir}/HE_bench/EUR/compare_to_GREML/GREML_prior_NA0.1_neg0_gpos1.hsq|grep -w "${phe}"|awk '{print $3}')
prVNTR=0.001

sprSNV=$(cat ${spe_baseDir}/HE_bench/EUR/continuous/EUR.${phe}.maf0.01.hwe-6.chrAuto.SNV.HEreg|sed -n '4p'|awk -v OFS='\t' '{print $2}')
sprSV=$(cat ${spe_baseDir}/HE_bench/EUR/continuous/EUR.${phe}.maf0.01.hwe-6.chrAuto.SV.HEreg|sed -n '4p'|awk -v OFS='\t' '{print $2}')
sprVNTR=0.002


	echo -e "SGV.unRel_sample_size_${sample_size}" > ${GRM_spe_D}/GRM/one_ck/zFilePath/${phe}.mgrm.2P.path
	echo -e "${GRM_spe_D}/GRM/one_ck/SV.unRel_sample_size_${sample_size}" >> ${GRM_spe_D}/GRM/one_ck/zFilePath/${phe}.mgrm.2P.path

	gcta64 \
	--reml \
	--reml-no-lrt \
	--reml-no-constrain \
	--reml-priors ${prSNV} ${prSV} \
	--mgrm ${GRM_spe_D}/GRM/one_ck/zFilePath/${phe}.mgrm.2P.path \
	--qcovar ${spe_baseDir}/PHENOTYPE/covar/UKB_All_covariates.covar \
	--pheno ${spe_baseDir}/PHENOTYPE/Continuous_final_652/${phe}.pheno \
	--thread-num 10 \
	--out ${GREML_D}/${pop}_bench_meta/cont_mix/zSpe_one_ck/res/${pop}.${phe}.maf0.01.hwe-6.chrAuto.2P

	gcta64 \
	--reml \
	--reml-no-lrt \
	--reml-no-constrain \
	--reml-priors ${sprSNV} \
	--grm SGV.unRel_sample_size_${sample_size} \
	--qcovar ${spe_baseDir}/PHENOTYPE/covar/UKB_All_covariates.covar \
	--pheno ${spe_baseDir}/PHENOTYPE/Continuous_final_652/${phe}.pheno \
	--thread-num 10 \
	--out ${GREML_D}/${pop}_bench_meta/cont_mix/zSpe_one_ck/res/${pop}.${phe}.maf0.01.hwe-6.chrAuto.SGV

	gcta64 \
	--reml \
	--reml-no-lrt \
	--reml-no-constrain \
	--reml-priors ${sprSV} \
	--grm ${GRM_spe_D}/GRM/one_ck/SV.unRel_sample_size_${sample_size} \
	--qcovar ${spe_baseDir}/PHENOTYPE/covar/UKB_All_covariates.covar \
	--pheno ${spe_baseDir}/PHENOTYPE/Continuous_final_652/${phe}.pheno \
	--thread-num 10 \
	--out ${GREML_D}/${pop}_bench_meta/cont_mix/zSpe_one_ck/res/${pop}.${phe}.maf0.01.hwe-6.chrAuto.SV

for phe in $(cat ${GREML_D}/EUR_bench_meta/cont_mix/zSpe_one_ck/res.phe.id.456); do
	
	beta_SGV_J=$(cat ${GREML_D}/EUR_bench_meta/cont_mix/zSpe_one_ck/res/EUR.${phe}.maf0.01.hwe-6.chrAuto.2P.hsq |grep -w 'V(G1)/Vp'|awk -F'\t' '{print $2}')
	se_SGV_J=$(cat ${GREML_D}/EUR_bench_meta/cont_mix/zSpe_one_ck/res/EUR.${phe}.maf0.01.hwe-6.chrAuto.2P.hsq |grep -w 'V(G1)/Vp'|awk -F'\t' '{print $3}')

	beta_SV_J=$(cat ${GREML_D}/EUR_bench_meta/cont_mix/zSpe_one_ck/res/EUR.${phe}.maf0.01.hwe-6.chrAuto.2P.hsq |grep -w 'V(G2)/Vp'|awk -F'\t' '{print $2}')
	se_SV_J=$(cat ${GREML_D}/EUR_bench_meta/cont_mix/zSpe_one_ck/res/EUR.${phe}.maf0.01.hwe-6.chrAuto.2P.hsq |grep -w 'V(G2)/Vp'|awk -F'\t' '{print $3}')

	beta_SUM_J=$(cat ${GREML_D}/EUR_bench_meta/cont_mix/zSpe_one_ck/res/EUR.${phe}.maf0.01.hwe-6.chrAuto.2P.hsq |grep -w 'Sum of V(G)/Vp'|awk -F'\t' '{print $2}')
	se_SUM_J=$(cat ${GREML_D}/EUR_bench_meta/cont_mix/zSpe_one_ck/res/EUR.${phe}.maf0.01.hwe-6.chrAuto.2P.hsq |grep -w 'Sum of V(G)/Vp'|awk -F'\t' '{print $3}')

	beta_SGV_M=$(cat ${GREML_D_SGV_m}/EUR_bench_meta/cont_mix/zSpe_one_ck/res/EUR.${phe}.maf0.01.hwe-6.chrAuto.SGV.hsq |grep -w 'V(G)/Vp'|awk -F'\t' '{print $2}')
	se_SGV_M=$(cat ${GREML_D_SGV_m}/EUR_bench_meta/cont_mix/zSpe_one_ck/res/EUR.${phe}.maf0.01.hwe-6.chrAuto.SGV.hsq |grep -w 'V(G)/Vp'|awk -F'\t' '{print $3}')

	beta_SV_M=$(cat ${GREML_D}/EUR_bench_meta/cont_mix/zSpe_one_ck/res/EUR.${phe}.maf0.01.hwe-6.chrAuto.SV.hsq |grep -w 'V(G)/Vp'|awk -F'\t' '{print $2}')
	se_SV_M=$(cat ${GREML_D}/EUR_bench_meta/cont_mix/zSpe_one_ck/res/EUR.${phe}.maf0.01.hwe-6.chrAuto.SV.hsq |grep -w 'V(G)/Vp'|awk -F'\t' '{print $3}')

	echo -e "cont\t${phe}\tck${i}\t${beta_SGV_J}\t${se_SGV_J}\t${beta_SV_J}\t${se_SV_J}\t${beta_SUM_J}\t${se_SUM_J}\t${beta_SGV_M}\t${se_SGV_M}\t${beta_SV_M}\t${se_SV_M}" >\
	${GREML_D}/EUR_bench_meta/cont_mix/zSpe_one_ck/inte/${phe}.res

	Rscript --slave --no-save meta_GREML_ck_res_529-652.R \
	${phe} cont \
	${GREML_D}/EUR_bench_meta/cont_mix/zSpe_one_ck/inte/${phe}.res \
	${GREML_D}/EUR_bench_meta/cont_mix/zSpe_one_ck/inte/${phe}.res.withP

done

# run 460-528 with two chunk#
phe=$(cat sample.size.sort.fid|sed -n "${pN}p"|awk '{print $2}')
sample_size=$(cat sample.size.sort.fid|sed -n "${pN}p"|awk '{print $1}')

for ck in 1 2;do

	echo -e "analysis/GRM_EUR/GRM_spe_phe_chunks_cont_652/GRM/two_ck/SGV.unRel_sample_size_${sample_size}_sck${ck}" > ${GRM_spe_D}/GRM/two_ck/zFilePath/${phe}.mgrm.2P.path_sck${ck}
	echo -e "${GRM_spe_D}/GRM/two_ck/SV.unRel_sample_size_${sample_size}_sck${ck}" >> ${GRM_spe_D}/GRM/two_ck/zFilePath/${phe}.mgrm.2P.path_sck${ck}

	# two componants
	gcta64 \
	--reml \
	--reml-no-lrt \
	--reml-no-constrain \
	--reml-priors ${prSNV} ${prSV} \
	--mgrm ${GRM_spe_D}/GRM/two_ck/zFilePath/${phe}.mgrm.2P.path_sck${ck} \
	--qcovar ${spe_baseDir}/PHENOTYPE/covar/UKB_All_covariates.covar \
	--pheno ${spe_baseDir}/PHENOTYPE/Continuous_final_652/${phe}.pheno \
	--thread-num 10 \
	--out ${GREML_D}/${pop}_bench_meta/cont_mix/zSpe_two_ck/res_ck${ck}/${pop}.${phe}.maf0.01.hwe-6.chrAuto.2P

	# single SGV
	gcta64 \
	--reml \
	--reml-no-lrt \
	--reml-no-constrain \
	--reml-priors ${sprSNV} \
	--grm analysis/GRM_EUR/GRM_spe_phe_chunks_cont_652/GRM/two_ck/SGV.unRel_sample_size_${sample_size}_sck${ck} \
	--qcovar ${spe_baseDir}/PHENOTYPE/covar/UKB_All_covariates.covar \
	--pheno ${spe_baseDir}/PHENOTYPE/Continuous_final_652/${phe}.pheno \
	--thread-num 10 \
	--out ${GREML_D}/${pop}_bench_meta/cont_mix/zSpe_two_ck/res_ck${ck}/${pop}.${phe}.maf0.01.hwe-6.chrAuto.SGV

	# single SV
	gcta64 \
	--reml \
	--reml-no-lrt \
	--reml-no-constrain \
	--reml-priors ${sprSV} \
	--grm ${GRM_spe_D}/GRM/two_ck/SV.unRel_sample_size_${sample_size}_sck${ck} \
	--qcovar ${spe_baseDir}/PHENOTYPE/covar/UKB_All_covariates.covar \
	--pheno ${spe_baseDir}/PHENOTYPE/Continuous_final_652/${phe}.pheno \
	--thread-num 10 \
	--out ${GREML_D}/${pop}_bench_meta/cont_mix/zSpe_two_ck/res_ck${ck}/${pop}.${phe}.maf0.01.hwe-6.chrAuto.SV

done

for ck in $(seq 1 7);do 

	echo -e "GRM_EUR/GRM_chunks/EUR.maf0.01.hwe-6.chrAuto.SGV.sck${ck}" > ${GRM_D}/GRM_${pop}/GRM_chunks/zFilePath/${pop}.${phe}.maf0.01.hwe-6.chrAuto.sck${ck}.2P.path
	echo -e "${GRM_D}/GRM_${pop}/GRM_chunks/${pop}.maf0.01.hwe-6.chrAuto.SV.sck${ck}" >> ${GRM_D}/GRM_${pop}/GRM_chunks/zFilePath/${pop}.${phe}.maf0.01.hwe-6.chrAuto.sck${ck}.2P.path

	mkdir -p ${GREML_D}/${pop}_bench_meta/cont_mix/ck${ck}

	# two componants
	gcta64 \
	--reml \
	--reml-no-lrt \
	--reml-no-constrain \
	--reml-priors ${prSNV} ${prSV}  \
	--mgrm ${GRM_D}/GRM_${pop}/GRM_chunks/zFilePath/${pop}.${phe}.maf0.01.hwe-6.chrAuto.sck${ck}.2P.path \
	--qcovar PHENOTYPE/covar/UKB_All_covariates.covar \
	--pheno PHENOTYPE/Continuous_final_652/${phe}.pheno \
	--thread-num 10 \
	--out ${GREML_D}/${pop}_bench_meta/cont_mix/ck${ck}/${pop}.${phe}.maf0.01.hwe-6.chrAuto.2P

	# single SGV
	gcta64 \
	--reml \
	--reml-no-lrt \
	--reml-no-constrain \
	--reml-priors ${prSNV} \
	--grm GRM_EUR/GRM_chunks/EUR.maf0.01.hwe-6.chrAuto.SGV.sck${ck} \
	--qcovar PHENOTYPE/covar/UKB_All_covariates.covar \
	--pheno PHENOTYPE/Continuous_final_652/${phe}.pheno \
	--thread-num 10 \
	--out ${GREML_D}/${pop}_bench_meta/cont_mix/ck${ck}/${pop}.${phe}.maf0.01.hwe-6.chrAuto.SGV

	# single SV
	gcta64 \
	--reml \
	--reml-no-lrt \
	--reml-no-constrain \
	--reml-priors ${prSV} \
	--grm ${GRM_D}/GRM_${pop}/GRM_chunks/${pop}.maf0.01.hwe-6.chrAuto.SV.sck${ck} \
	--qcovar PHENOTYPE/covar/UKB_All_covariates.covar \
	--pheno PHENOTYPE/Continuous_final_652/${phe}.pheno \
	--thread-num 10 \
	--out ${GREML_D}/${pop}_bench_meta/cont_mix/ck${ck}/${pop}.${phe}.maf0.01.hwe-6.chrAuto.SV
done


