#!/bin/bash

c=${SLURM_ARRAY_TASK_ID}

GP=HGSVC2_HPRCY1
M=M24sl5F

GTEx_ori_TR_impu_D=GTEx_77/zTR
eQTL_D=eQTL

TR_ref_F=${GP}/${M}_len7_InDel_SV_match_HG38/${GP}.HG38.${M}.VNTR.chrAuto.summary
loci_D=${GP}/${M}_len7_InDel_SV_match_HG38/zDense_VNTR_loci_with_TR_allele

pop=EUR


# get EUR imputed TR 
mkdir -p ${eQTL_D}/GTEx_v8/VNTR/imputed/EUR

bcftools view -Oz --samples-file ${eQTL_D}/GTEx_v8/sample/EUR.sample.oneCol ${GTEx_ori_TR_impu_D}/QC.GRCh38.WGS_chr${c}.GTEx_v8.${GP}.${M}.VNTR.impu.dose.TR.vcf.gz |\
bcftools annotate -Oz --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' --output ${eQTL_D}/GTEx_v8/VNTR/imputed/EUR/QC.GRCh38.WGS_chr${c}.GTEx_v8.${GP}.${M}.VNTR.impu.dose.TR.EUR.vcf.gz

tabix -p vcf ${eQTL_D}/GTEx_v8/VNTR/imputed/EUR/QC.GRCh38.WGS_chr${c}.GTEx_v8.${GP}.${M}.VNTR.impu.dose.TR.EUR.vcf.gz


aim_D=${eQTL_D}/GTEx_v8/VNTR/prepare_GENO/${pop}/chr${c}

mkdir -p ${aim_D}/ale_vid
mkdir -p ${aim_D}/ale_vcf
mkdir -p ${aim_D}/ale_gt
mkdir -p ${aim_D}/ale_cn


for loci_id in $(cat ${TR_ref_F}|grep "_chr${c}_"|awk '{print $1}'); do

	# get VNTR TR ID
	cat ${loci_D}/${loci_id}.tsv|sed '1d'|awk '{print $18}' > ${aim_D}/ale_vid/${loci_id}.vid

	# extract imputed TR by ID, and filter estimated MAF < 0.01
	bcftools view -Oz --include ID=@${aim_D}/ale_vid/${loci_id}.vid ${eQTL_D}/GTEx_v8/VNTR/imputed/EUR/QC.GRCh38.WGS_chr${c}.GTEx_v8.${GP}.${M}.VNTR.impu.dose.TR.EUR.vcf.gz |\
	bcftools view -Oz --exclude 'INFO/MAF<0.01' |\
	bcftools annotate -Ov -x FILTER,INFO,FORMAT/GP --output ${aim_D}/ale_vcf/${loci_id}.dos-maf01.vcf

	# extract imputed TR by ID, and filter real MAF < 0.01
	bcftools view -Oz --include ID=@${aim_D}/ale_vid/${loci_id}.vid ${eQTL_D}/GTEx_v8/VNTR/imputed/EUR/QC.GRCh38.WGS_chr${c}.GTEx_v8.${GP}.${M}.VNTR.impu.dose.TR.EUR.vcf.gz |\
	bcftools view -Oz --min-af 0.01:minor |\
	bcftools annotate -Ov -x FILTER,INFO,FORMAT/GP --output ${aim_D}/ale_vcf/${loci_id}.int-maf01.vcf

	# deal genotype, note that this mtp is for MAF only, i.e. dos-maf and int-maf, nothing to to with genotype
	for mtp in dos int; do

		# check num of TR in the loci, only loci with >= 2 allele were retained for analyais
		TR_n=$(cat ${aim_D}/ale_vcf/${loci_id}.${mtp}-maf01.vcf|grep -ac ^c)

		if [ ${TR_n} -gt 1 ]; then

			# get genotype int:dos for each loci
			cat ${aim_D}/ale_vcf/${loci_id}.${mtp}-maf01.vcf | grep -v '^##' |\
			sed 's/0|0/0/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/2/g' | awk -v OFS='_' '{$1="";$2="";$4="";$5="";$6="";$7="";$8="";$9=""; print $0}' |\
			sed 's/__*/_/g' | sed 's/^_//g' > ${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt
			sed -i 's/_/ /g' ${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt

			# get cn for each ale
			for tr_id in $(cat ${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt|grep -a ^c|awk '{print $1}');do
				echo ${tr_id} > ${aim_D}/ale_cn/${loci_id}.${mtp}-maf01.vid-cn.tmp_trid
				tr_id_cn=$(cat ${loci_D}/${loci_id}.tsv|grep -w -f ${aim_D}/ale_cn/${loci_id}.${mtp}-maf01.vid-cn.tmp_trid|awk '{print "'${tr_id}'",$3,$9}'|sed 's/INS //'|sed 's/DEL /-/')
				echo -e ${tr_id_cn} >> ${aim_D}/ale_cn/${loci_id}.${mtp}-maf01.vid-cn
			done

			# get gt * cn, step 1, int and dos in one file
			python /home/yangjianLab/baiweiyang/scripts/SV_Imputation_project_analysis/06.VNTR/1_make_VNTR_genotype/1_get_VNTR_loci-geno-cn_step-1.py \
			${aim_D}/ale_cn/${loci_id}.${mtp}-maf01.vid-cn \
			${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt \
			${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt-cn

			# get gt * cn, step 2, int and dos in two files
			python /home/yangjianLab/baiweiyang/scripts/SV_Imputation_project_analysis/06.VNTR/1_make_VNTR_genotype/1_get_VNTR_loci-geno-cn_step-2.py \
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
				python /home/yangjianLab/baiweiyang/scripts/SV_Imputation_project_analysis/06.VNTR/1_make_VNTR_genotype/1_get_VNTR_loci-geno-cn_step-3.py \
				${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt-cn-${gtp} \
				${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt-cn-${gtp}.T

				cat ${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt-cn-${gtp}.T | awk '{print $2}' > ${aim_D}/ale_gt/${loci_id}.${mtp}-maf01.vid-gt-cn-${gtp}.T.tmp

			done
		fi
	done
done

# get togther
for mtp in dos int; do
	for gtp in dos int; do

		# get sample ID from anyone loci_id with N2.
		cd ${aim_D}/ale_gt
		cat $(ls | grep ${mtp}-maf01.vid-gt-cn-${gtp}.T|head -n1) | awk '{print $1}' > ${aim_D}/ale_gt/z.${mtp}-maf01.vid-gt-cn-${gtp}.T.colID

		# group all loci
		paste -d' ' ${aim_D}/ale_gt/*.${mtp}-maf01.vid-gt-cn-${gtp}.T.tmp > ${aim_D}/ale_gt/z.${mtp}-maf01.vid-gt-cn-${gtp}.T.colV

		# paste ID and geno
		paste -d' ' ${aim_D}/ale_gt/z.${mtp}-maf01.vid-gt-cn-${gtp}.T.colID ${aim_D}/ale_gt/z.${mtp}-maf01.vid-gt-cn-${gtp}.T.colV >\
		${aim_D}/ale_gt/z.${pop}.chr${c}.N2.${mtp}-maf01.vid-gt-cn-${gtp}.T

	done
done

# clean gloab
rm -rf ${aim_D}/ale_vid
rm -rf ${aim_D}/ale_vcf
rm -rf ${aim_D}/ale_cn

rm -rf ${aim_D}/ale_gt/*.vid-gt
rm -rf ${aim_D}/ale_gt/*.vid-gt-cn
rm -rf ${aim_D}/ale_gt/*.vid-gt-cn-dos
rm -rf ${aim_D}/ale_gt/*.vid-gt-cn-int

rm -rf ${aim_D}/ale_gt/*.T.tmp
rm -rf ${aim_D}/ale_gt/*.T.colID
rm -rf ${aim_D}/ale_gt/*.T.colV



tissue=$(cat ${eQTL_D}/GTEx_v8/info/tissue_chr.for_vntr_lm|sed -n "${n}p"|awk '{print $1}')
c=$(cat ${eQTL_D}/GTEx_v8/info/tissue_chr.for_vntr_lm|sed -n "${n}p"|awk '{print $2}')

exp_F=${eQTL_D}/GTEx_v8/expression/${pop}/chr${c}/${tissue}.chr${c}.v8.${pop}.normalized_expression.bed.T
cov_F=${eQTL_D}/GTEx_v8/covar/${pop}/${tissue}.v8.${pop}.covariates.txt.T


for mtp in dos int; do
	for gtp in dos int; do

		geno_F=${eQTL_D}/GTEx_v8/VNTR/prepare_GENO/${pop}/chr${c}/ale_gt/z.${pop}.chr${c}.N2.${mtp}-maf01.vid-gt-cn-${gtp}.T

		regression_D=${eQTL_D}/GTEx_v8/VNTR/regression/${pop}/${mtp}_maf01_${gtp}_geno_rsq00/${tissue}
		mkdir -p ${regression_D}

		res_F=${regression_D}/GTEx_v8.VNTR.${pop}.chr${c}.${tissue}.lr

		Rscript --slave --no-save /home/yangjianLab/baiweiyang/scripts/SV_Imputation_project_analysis/08.eQTL/2_regression_GETx_v8_VNTR_EUR/2_regression_GETx_v8_VNTR_EUR.R \
		${geno_F} ${exp_F} ${cov_F} ${res_F} 

	done
done


for mtp in dos int; do
	for gtp in dos int; do
		regression_summ_D=${eQTL_D}/GTEx_v8/VNTR/regression_summ/${pop}/${mtp}_maf01_${gtp}_geno_rsq00
		regression_summ_D=${eQTL_D}/GTEx_v8/zSumm_VNTR
		mkdir -p ${regression_summ_D}

		regression_summ_F=${regression_summ_D}/GTEx_v8.VNTR.${tissue}.${pop}.${mtp}.${gtp}.1e-5.summ

		echo -e "chrom\tgene_start\tgene_end\tstrand\tgene_id\tgene_type\tgene_symbol\tgene_midpos\tgene_id2\tvntr_id\tbeta\tse\tp\tsample_N\tVNTR_ID\tALE_N\tCHR\tREF_STR\tREF_END\tREF_SoP\tREF_SoC\tREF_CN\tREF_PoM\tREF_ACGT\tREF_Seq" >\
		${regression_summ_F}

		for c in $(seq 1 22);do

			ori_res_F=${eQTL_D}/GTEx_v8/VNTR/regression/${pop}/${mtp}_maf01_${gtp}_geno_rsq00/${tissue}/GTEx_v8.VNTR.${pop}.chr${c}.${tissue}.lr

			cat ${ori_res_F} | awk '{if($5<1e-5) print $0}' > ${ori_res_F}.1e-5

			N=$(cat ${eQTL_D}/GTEx_v8/expression/${pop}/chr${c}/${tissue}.chr${c}.v8.${pop}.normalized_expression.bed.T | sed '1d' | wc -l)

			total_sig_line=$(cat ${ori_res_F}.1e-5 | wc -l)

			for i in $(seq 1 ${total_sig_line});do
				entry=$(cat ${ori_res_F}.1e-5 | sed -n "${i}p" | sed 's/ /\t/g')

				ensg_id=$(echo ${entry} | awk -F'.' '{print $1}')
				vntr_id=$(echo ${entry} | awk '{print $2}')

				vntr_info=$(cat ${meta_VNTR_info_F} | grep -w ${vntr_id} | head -n1)
				gene_info_YoN=$(cat ${gencode_F} | grep -w ${ensg_id} | head -n1 | awk '{print $1}')
				gene_info=$(cat ${gencode_F} | grep -w ${ensg_id} | head -n1 | sed 's/ /\t/g')

				if [ ${gene_info_YoN} ]; then
					echo -e "${gene_info}\t${entry}\t${N}\t${vntr_info}" >> ${regression_summ_F}
				fi
			done
		done
	done
done



c=${SLURM_ARRAY_TASK_ID}

GP=HGSVC2_HPRCY1
M=M23s5F

SV_VCF_F=QC.GRCh38.WGS_chr${c}.GTEx_v8.${GP}.${M}.impu.dose.SV.vcf.gz
SGV_VCF_F=QC.GRCh38.WGS_chr${c}.GTEx_v8.vcf.gz

sample_F=EUR.sample.twoCol
GENO_out_D=zGENOTYPE


module load plink/2

plink2 --vcf ${SGV_VCF_F} --double-id --keep ${sample_F} --maf 0.01 --make-bed --out ${GENO_out_D}/EUR.SGV.WGS.maf01.chr${c}

plink2 --vcf ${SV_VCF_F} --double-id --keep ${sample_F} --maf 0.01 --make-bed --out ${GENO_out_D}/EUR.SV.${M}.maf01.chr${c}.tmp


python /home/yangjianLab/baiweiyang/scripts/SV_Imputation_project_analysis/02.GENOTYPE/1_prepare_geno_SGV_SV/1_prepare_geno_SGV_SV.py \
${GENO_out_D}/EUR.SV.${M}.maf01.chr${c}.tmp.bim \
/storage/yangjianLab/baiweiyang/SV_Imputation_Project_panel/Joint_HPRC_ccs47_HGSVC2_ccs9clr21_CKCG_ccs55_132/panels/${GP}/SV_ID/${M}/${GP}.chr${c}.${M}.SV.chr_pos_svid_type_ref_alt.tsv \
${GENO_out_D}/EUR.SV.${M}.maf01.chr${c}.tmp.bim_new

mv ${GENO_out_D}/EUR.SV.${M}.maf01.chr${c}.tmp.bim_new ${GENO_out_D}/EUR.SV.${M}.maf01.chr${c}.tmp.bim

cat ${GENO_out_D}/EUR.SV.${M}.maf01.chr${c}.tmp.bim |awk '{print $2}'|grep HGSVC > ${GENO_out_D}/EUR.SV.${M}.maf01.chr${c}.tmp.id

plink2 --bfile ${GENO_out_D}/EUR.SV.${M}.maf01.chr${c}.tmp --extract ${GENO_out_D}/EUR.SV.${M}.maf01.chr${c}.tmp.id --make-bed --out ${GENO_out_D}/EUR.SV.${M}.maf01.chr${c}

rm ${GENO_out_D}/EUR.SV.${M}.maf01.chr${c}.tmp*


sample_F=EUR.sample.twoCol

c_covar_F=${eQTL_D}/GTEx_v8/Stand_alone_SGV_SV/EUR/${tissue}/covar/${tissue}_covariates.txt.covar
q_covar_F=${eQTL_D}/GTEx_v8/Stand_alone_SGV_SV/EUR/${tissue}/covar/${tissue}_covariates.txt.qcovar

bod_F=${eQTL_D}/GTEx_v8/Stand_alone_SGV_SV/EUR/${tissue}/expression/${tissue}_counts_TMM_RINT

SV_geno_F=${eQTL_D}/GTEx_v8/Stand_alone_SGV_SV/EUR/zGENOTYPE/EUR.SV.${M}.maf01.chr${c}

out_whole_SV_prefix=${eQTL_D}/GTEx_v8/Stand_alone_SGV_SV/EUR/${tissue}/eQTL/SV.${M}/chr/${tissue}.chr${c}.eQTL.GRCh38.maf01.SV.${M}
out_query_SV_prefix=${eQTL_D}/GTEx_v8/Stand_alone_SGV_SV/EUR/${tissue}/eQTL/SV.${M}/query/${tissue}.chr${c}.eQTL.GRCh38.maf01.SV.${M}

# SV
${osca} --eqtl --bfile ${SV_geno_F} --befile ${bod_F} --covar ${c_covar_F} --qcovar ${q_covar_F} --thread-num ${tr} --chr ${c} --to-smr --out ${out_whole_SV_prefix}

${smr} --beqtl-summary ${out_whole_SV_prefix} --descriptive-cis --peqtl-cis 5e-8 --out ${out_query_SV_prefix}
${smr} --beqtl-summary ${out_whole_SV_prefix} --descriptive-trans --peqtl-trans 5e-8 --out ${out_query_SV_prefix}

${smr} --beqtl-summary ${out_whole_SV_prefix} --query 1 --out ${out_query_SV_prefix} && gzip -f ${out_query_SV_prefix}.txt


# SGV
SGV_geno_F=${eQTL_D}/GTEx_v8/Stand_alone_SGV_SV/EUR/zGENOTYPE/EUR.SGV.WGS.maf01.chr${c}

out_whole_SGV_prefix=${eQTL_D}/GTEx_v8/Stand_alone_SGV_SV/EUR/${tissue}/eQTL/SGV.WGS/chr/${tissue}.chr${c}.eQTL.GRCh38.maf01.SGV.WGS

#out_query_SGV_prefix=${eQTL_D}/GTEx_v8/Stand_alone_SGV_SV/EUR/${tissue}/eQTL/SGV.WGS/query/${tissue}.chr${c}.eQTL.GRCh38.maf01.SGV.WGS

${osca} --eqtl --bfile ${SGV_geno_F} --befile ${bod_F} --covar ${c_covar_F} --qcovar ${q_covar_F} --thread-num ${tr} --cis --cis-wind 2000 --chr ${c} --to-smr --out ${out_whole_SGV_prefix}

## not query since it could be huge, if need, please specific (e.g. based on SVs results)


