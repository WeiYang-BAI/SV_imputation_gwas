#!/bin/bash

#########################################################################################
## prepare: stage 2 variants for each sample
# 1) SNP: assembly-based called, QCed and matched (pav + dipcall)
# 2) InDel: assembly-based called, QCed and matched (pav + dipcall)
# 3) SV: assembly-based called, QCed and matched (pav main; dipcall and svim support)
# 4) rename, clear, and tabix all: {SID}.stage2.{vTYPE}.chr{c}.vcf.gz

## prepare: 1KG variants for each sample (if exists)
# 1) SNP: high-depth short read called by 1KG
# 2) rename, clear and tabix all: {SID}.1KG.SNP.chr{c}.vcf.gz
#########################################################################################


sn=${SLURM_ARRAY_TASK_ID}

base_D=panel

infos_132_F=${base_D}/infos/SID_gender_seqTYPE_superPOP_subPOP_inOKG_SOURCE.132

SID=$(cat ${infos_132_F} | sed -n "${sn}p" | awk '{print $1}')
inOKG=$(cat ${infos_132_F} | sed -n "${sn}p" | awk '{print $6}')


# check phasing and make stg 3 snp

stage2_geno_SNP_F=${base_D}/stage_2_calls/SNP/${SID}/${SID}.stage2.SNP.chrAuto.FMT.vcf.gz
stage2_geno_SV_F=${base_D}/stage_2_calls/SV/${SID}/${SID}.stage2.SV.chrAuto.FMT.vcf.gz # SV for exclude SNps

mkdir -p ${base_D}/stage_3_calls/SNP/${SID}

stage3_geno_SNP_F=${base_D}/stage_3_calls/SNP/${SID}/${SID}.stage3.SNP.chrAuto.FMT.vcf.gz


if [ $inOKG == 'yes' ]; then

	# check phasing with whatshap
	hash -r
	source ~/.bashrc
	conda deactivate
	conda activate whatshap

	okg_geno_F=${base_D}/1KG_calls/SNP/${SID}/${SID}.1KG.SNP.chrAuto.FMT.vcf.gz
	out_whatshap_F=${base_D}/check_phasing_whatshap/SNP/${SID}.stage2_vs_okg.SNP.chrAuto.tsv

	whatshap compare --names stage2,okg --tsv-pairwise ${out_whatshap_F} ${stage2_geno_SNP_F} ${okg_geno_F}

	# manually check
	mkdir -p ${base_D}/check_phasing_refine_manually/SNP/${SID}
	isec_D=${base_D}/check_phasing_refine_manually/SNP/${SID}

	bcftools isec -Oz -p ${isec_D} ${stage2_geno_SNP_F} ${okg_geno_F}

	bcftools sort -Oz --temp-dir ${isec_D} -o ${isec_D}/sorted.0002.vcf.gz ${isec_D}/0002.vcf.gz && tabix -p vcf ${isec_D}/sorted.0002.vcf.gz
	bcftools sort -Oz --temp-dir ${isec_D} -o ${isec_D}/sorted.0003.vcf.gz ${isec_D}/0003.vcf.gz && tabix -p vcf ${isec_D}/sorted.0003.vcf.gz

	# 0002: HPRC
	# 0003: OKG

	mkdir -p ${base_D}/check_phasing_refine_manually/SNP/${SID}/zMake
	make_D=${base_D}/check_phasing_refine_manually/SNP/${SID}/zMake

	# please check there is no missing gt in the data (e.g., 1|., 0|., etc.)
	zcat ${isec_D}/sorted.0002.vcf.gz|grep -v '^#'|awk -v OFS='\t' '{print $1,$2,$4,$5,$10}'|awk -F':' '{print $1}'|sed 's/|/\t/g' > ${make_D}/stage2-0002-haps_infos.txt
	zcat ${isec_D}/sorted.0003.vcf.gz|grep -v '^#'|awk -v OFS='\t' '{print $1,$2,$4,$5,$10}'|awk -F':' '{print $1}'|sed 's/|/\t/g' > ${make_D}/okg-0003-haps_infos.txt 

	# put 0002 and 0003 in one file
	paste -d'\t' ${make_D}/stage2-0002-haps_infos.txt ${make_D}/okg-0003-haps_infos.txt > ${make_D}/both-haps_infos.txt

	# make shared genotype (diff=0 use stage2 , diff=1 use okg)
	cat ${make_D}/both-haps_infos.txt|awk '{print $5,$6,$11,$12,($5-$11)+($6-$12)}'|awk -v OFS='|' '{if($5==0) {print $1,$2} else {print $3,$4}}' > ${make_D}/shared.geno

	# get header and info
	zcat ${isec_D}/sorted.0002.vcf.gz|grep '^#' > ${make_D}/shared.header
	zcat ${isec_D}/sorted.0002.vcf.gz|grep -v '^#'|awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9}'> ${make_D}/shared.info

	# gather header info gt
	paste -d'\t' ${make_D}/shared.info ${make_D}/shared.geno > ${make_D}/shared.body && cat ${make_D}/shared.header ${make_D}/shared.body > ${make_D}/shared.vcf
	bgzip ${make_D}/shared.vcf && tabix -p vcf ${make_D}/shared.vcf.gz

	# concat shared and stage2 uniq geno
	bcftools concat -Oz -a ${isec_D}/0000.vcf.gz ${make_D}/shared.vcf.gz -o ${isec_D}/refined.vcf.gz && tabix -p vcf ${isec_D}/refined.vcf.gz && rm -rf ${make_D}

	# exclude SV-overlap SNPs and sort, form stage3 SNPs
	bedtools intersect -a ${isec_D}/refined.vcf.gz -b ${stage2_geno_SV_F} | awk -v OFS='\t' '{print $1,$2}'|sort|uniq > ${base_D}/stage_3_calls/SNP/${SID}/${SID}.SNP.excl

	if [ -s ${base_D}/stage_3_calls/SNP/${SID}/${SID}.SNP.excl ]; then
		bcftools view -Oz --targets-file ^${base_D}/stage_3_calls/SNP/${SID}/${SID}.SNP.excl ${isec_D}/refined.vcf.gz |	bcftools sort -Oz --temp-dir ${base_D}/stage_3_calls/SNP/${SID}/tmp -o ${stage3_geno_SNP_F} && tabix -p vcf ${stage3_geno_SNP_F}
	else
		bcftools sort -Oz --temp-dir ${base_D}/stage_3_calls/SNP/${SID}/tmp -o ${stage3_geno_SNP_F} ${isec_D}/refined.vcf.gz && tabix -p vcf ${stage3_geno_SNP_F}
	fi

	# check stage2 vs okg, mismatch
	hash -r
	source ~/.bashrc
	conda deactivate
	conda activate pyvcf

	zcat ${isec_D}/sorted.0002.vcf.gz|grep -v '^#'|awk '{print $10}'|awk -F':' '{print $1}' > ${isec_D}/stage2_0002_haps.txt
	zcat ${isec_D}/sorted.0003.vcf.gz|grep -v '^#'|awk '{print $10}'|awk -F':' '{print $1}' > ${isec_D}/okg_0003_haps.txt

	python3 check_err.py \
	${isec_D}/stage2_0002_haps.txt \
	${isec_D}/okg_0003_haps.txt \
	${isec_D}/stage2_okg_haps_diff_stat.txt \
	${SID}

else
	# exclude SV-overlap SNPs and sort, form stage3 SNPs
	bedtools intersect -a ${stage2_geno_SNP_F} -b ${stage2_geno_SV_F} | awk -v OFS='\t' '{print $1,$2}'|sort|uniq > ${base_D}/stage_3_calls/SNP/${SID}/${SID}.SNP.excl

	if [ -s ${base_D}/stage_3_calls/SNP/${SID}/${SID}.SNP.excl ]; then
		bcftools view -Oz --targets-file ^${base_D}/stage_3_calls/SNP/${SID}/${SID}.SNP.excl ${stage2_geno_SNP_F} |	bcftools sort -Oz --temp-dir ${base_D}/stage_3_calls/SNP/${SID}/tmp -o ${stage3_geno_SNP_F} && tabix -p vcf ${stage3_geno_SNP_F}
	else
		bcftools sort -Oz --temp-dir ${base_D}/stage_3_calls/SNP/${SID}/tmp -o ${stage3_geno_SNP_F} ${stage2_geno_SNP_F} && tabix -p vcf ${stage3_geno_SNP_F}
	fi
fi

# split into chrs
for c in $(seq 1 22); do
	stage3_geno_chr_SNP_F=${base_D}/stage_3_calls/SNP/${SID}/${SID}.stage3.SNP.chr${c}.FMT.vcf.gz
	bcftools view -Oz --targets chr${c} -o ${stage3_geno_chr_SNP_F} ${stage3_geno_SNP_F} && tabix -p vcf ${stage3_geno_chr_SNP_F}
done



# please make sure chrAuto file is clear first
for vTYPE in SNP InDel;do
	stage2_vcf_chrAuto_F=${base_D}/stage_2_calls/${vTYPE}/${SID}/${SID}.stage2.${vTYPE}.chrAuto.vcf.gz
	stage2_vcf_chrAuto_fmt_F=${base_D}/stage_2_calls/${vTYPE}/${SID}/${SID}.stage2.${vTYPE}.chrAuto.FMT.vcf.gz
	zcat ${stage2_vcf_chrAuto_F} | grep ^# > ${stage2_vcf_chrAuto_F}.header
	zcat ${stage2_vcf_chrAuto_F} | grep -v ^# > ${stage2_vcf_chrAuto_F}.body
	cat ${stage2_vcf_chrAuto_F}.header ${stage2_vcf_chrAuto_F}.body > ${stage2_vcf_chrAuto_F}.vcf && rm ${stage2_vcf_chrAuto_F}.header ${stage2_vcf_chrAuto_F}.body
	bgzip ${stage2_vcf_chrAuto_F}.vcf && mv ${stage2_vcf_chrAuto_F}.vcf.gz ${stage2_vcf_chrAuto_fmt_F} && tabix -p vcf ${stage2_vcf_chrAuto_fmt_F}
	for c in $(seq 1 22);do
		stage2_vcf_chrN_fmt_F=${base_D}/stage_2_calls/${vTYPE}/${SID}/${SID}.stage2.${vTYPE}.chr${c}.FMT.vcf.gz
		bcftools view -Oz -t chr${c} -o ${stage2_vcf_chrN_fmt_F} ${stage2_vcf_chrAuto_fmt_F} && tabix -p vcf ${stage2_vcf_chrN_fmt_F}
	done
done

# NOT lowercase to uppercase for stage2 SV and InDel, please make sure chrAuto file is clear first
for vTYPE in SV;do
	stage2_vcf_chrAuto_F=${base_D}/stage_2_calls/${vTYPE}/${SID}/${SID}.stage2.${vTYPE}.chrAuto.vcf.gz
	for c in $(seq 1 22);do
		stage2_vcf_chrN_fmt_F=${base_D}/stage_2_calls/${vTYPE}/${SID}/${SID}.stage2.${vTYPE}.chr${c}.FMT.vcf.gz
		bcftools view -Oz -t chr${c} -o ${stage2_vcf_chrN_fmt_F} ${stage2_vcf_chrAuto_F} && tabix -p vcf ${stage2_vcf_chrN_fmt_F}
	done
done


# subset sample and clear 1KG variants
if [ $inOKG == 'yes' ]; then
	for c in $(seq 1 22);do
		okg_vcf_all=/storage/yangjianLab/sharedata/LD_reference/1000G/1KGP3_3202_GRCh38/1kGP_high_coverage_Illumina.chr${c}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
		# SNP
		vTYPE=SNP
		mkdir -p ${base_D}/1KG_calls/${vTYPE}/${SID}
		out_F=${base_D}/1KG_calls/${vTYPE}/${SID}/${SID}.1KG.${vTYPE}.chr${c}.FMT.vcf.gz
		bcftools view -Oz -s ${SID} ${okg_vcf_all}|bcftools view -Oz --min-ac 1:nref|bcftools view -Oz -v snps|bcftools annotate -Oz -x ID,INFO,QUAL -o ${out_F} && tabix -f -p vcf ${out_F}
		echo "${out_F}" >> ${base_D}/1KG_calls/${vTYPE}/${SID}/${SID}.1KG.${vTYPE}.chrAuto.FMT.VCFFILE

		# InDel
		vTYPE=InDel
		mkdir -p ${base_D}/1KG_calls/${vTYPE}/${SID}
		out_F=${base_D}/1KG_calls/${vTYPE}/${SID}/${SID}.1KG.${vTYPE}.chr${c}.FMT.vcf.gz
		bcftools view -Oz -s ${SID} ${okg_vcf_all}|bcftools view -Oz --min-ac 1:nref|bcftools view -Oz --include 'ILEN >= -49 && ILEN <= 49'|bcftools annotate -Oz -x ID,INFO,QUAL -o ${out_F} && tabix -f -p vcf ${out_F}
		echo "${out_F}" >> ${base_D}/1KG_calls/${vTYPE}/${SID}/${SID}.1KG.${vTYPE}.chrAuto.FMT.VCFFILE

		# SV
		vTYPE=SV
		mkdir -p ${base_D}/1KG_calls/${vTYPE}/${SID}
		out_F=${base_D}/1KG_calls/${vTYPE}/${SID}/${SID}.1KG.${vTYPE}.chr${c}.FMT.vcf.gz
		bcftools view -Oz -s ${SID} ${okg_vcf_all}|bcftools view -Oz --min-ac 1:nref|bcftools view -Oz --include 'ILEN <= -50 || ILEN >= 50'|bcftools annotate -Oz -x ID,INFO,QUAL -o ${out_F} && tabix -f -p vcf ${out_F}
		echo "${out_F}" >> ${base_D}/1KG_calls/${vTYPE}/${SID}/${SID}.1KG.${vTYPE}.chrAuto.FMT.VCFFILE
	done

	for vTYPE in SNP InDel SV;do
		okg_vcf_chrAuto_fmt_F=${base_D}/1KG_calls/${vTYPE}/${SID}/${SID}.1KG.${vTYPE}.chrAuto.FMT.vcf.gz
		bcftools concat -Oz --file-list ${base_D}/1KG_calls/${vTYPE}/${SID}/${SID}.1KG.${vTYPE}.chrAuto.FMT.VCFFILE -o ${okg_vcf_chrAuto_fmt_F} && tabix -p vcf ${okg_vcf_chrAuto_fmt_F}
	done
fi



infos_D=${base_D}/infos
stage3_SNP_D=${base_D}/stage_3_calls/SNP
stage2_SNP_D=${base_D}/stage_2_calls/SNP
stage2_InDel_D=${base_D}/stage_2_calls/InDel

for vTYPE in SNP InDel;do
	for GP in PUB56 CN62 PUB1KCP111;do
		group_V_D=${base_D}/merged_sample/${vTYPE}/${GP}/group
		mkdir -p ${group_V_D}

		if [ $vTYPE == 'SNP' ];then
			
			# create merge sample list for SNP		
			for SID in $(cat ${infos_D}/SID_gender_seqTYPE_superPOP_subPOP_inOKG_SOURCE.CCS.${GP}|awk '{print $1}');do
				echo "${stage3_SNP_D}/${SID}/${SID}.stage3.${vTYPE}.chrAuto.FMT.vcf.gz" >> "${group_V_D}/${GP}.stage3.${vTYPE}.chrAuto.VCFFILE"
				echo "${stage2_SNP_D}/${SID}/${SID}.stage2.${vTYPE}.chrAuto.FMT.vcf.gz" >> "${group_V_D}/${GP}.stage2.${vTYPE}.chrAuto.VCFFILE"
			done
			
			# merge stage2 and 3 SNP
			for stg in 2 3; do

				group_prefix=${group_V_D}/${GP}.stage${stg}.${vTYPE}.chrAuto.zero.addID

				bcftools merge -Oz --file-list ${group_V_D}/${GP}.stage${stg}.${vTYPE}.chrAuto.VCFFILE --merge none --missing-to-ref --threads ${tr} |\
				bcftools annotate -Oz --set-id '%CHROM\:%POS\:%REF\:%ALT' |\
				bcftools sort -Ov --temp-dir ${group_V_D} -o ${group_prefix}.vcf

				cat ${group_prefix}.vcf | grep '^#' > ${group_prefix}.header
				cat ${group_prefix}.vcf | grep -v '^#' | sed 's/\//|/g' > ${group_prefix}.body
				cat ${group_prefix}.header ${group_prefix}.body > ${group_prefix}.phas.vcf

				bgzip ${group_prefix}.phas.vcf && tabix -f -p vcf ${group_prefix}.phas.vcf.gz && rm ${group_prefix}.header ${group_prefix}.body

				mkdir -p ${group_V_D}/zChr
				
				for c in $(seq 1 22); do
					bcftools view -Oz -t chr${c} --threads ${tr} -o ${group_V_D}/zChr/${GP}.stage${stg}.${vTYPE}.chr${c}.zero.addID.phas.vcf.gz ${group_prefix}.phas.vcf.gz &&\
					tabix -p vcf ${group_V_D}/zChr/${GP}.stage${stg}.${vTYPE}.chr${c}.zero.addID.phas.vcf.gz
				done
			done

		else
			# create merge sample list for InDel
			for SID in $(cat ${infos_D}/SID_gender_seqTYPE_superPOP_subPOP_inOKG_SOURCE.CCS.${GP}|awk '{print $1}');do
				echo "${stage2_InDel_D}/${SID}/${SID}.stage2.${vTYPE}.chrAuto.vcf.gz" >> "${group_V_D}/${GP}.stage2.${vTYPE}.chrAuto.VCFFILE"
			done
				
			group_prefix=${group_V_D}/${GP}.stage2.${vTYPE}.chrAuto.zero.addID

			bcftools merge -Oz --file-list ${group_V_D}/${GP}.stage2.${vTYPE}.chrAuto.VCFFILE --merge none --missing-to-ref --threads ${tr} |\
			bcftools annotate -Oz --set-id '%CHROM\:%POS\:%REF\:%ALT' |\
			bcftools sort -Ov --temp-dir ${group_V_D} -o ${group_prefix}.vcf

			cat ${group_prefix}.vcf | grep '^#' > ${group_prefix}.header
			cat ${group_prefix}.vcf | grep -v '^#' | sed 's/\//|/g' > ${group_prefix}.body
			cat ${group_prefix}.header ${group_prefix}.body > ${group_prefix}.phas.vcf

			bgzip ${group_prefix}.phas.vcf && tabix -f -p vcf ${group_prefix}.phas.vcf.gz && rm ${group_prefix}.header ${group_prefix}.body

			mkdir -p ${group_V_D}/zChr
				
			for c in $(seq 1 22); do
				bcftools view -Oz -t chr${c} --threads ${tr} -o ${group_V_D}/zChr/${GP}.stage2.${vTYPE}.chr${c}.zero.addID.phas.vcf.gz ${group_prefix}.phas.vcf.gz &&\
				tabix -p vcf ${group_V_D}/zChr/${GP}.stage2.${vTYPE}.chr${c}.zero.addID.phas.vcf.gz
			done
		fi
	done
done



for GP in PUB56 CN62 PUB1KCP111;do
	group_SV_D=${base_D}/merged_sample/SV/${GP}/group
	merge_SV_D=${base_D}/merged_sample/SV/${GP}/merge

	mkdir -p ${merge_SV_D}/SVanalyzer_distance
	mkdir -p ${merge_SV_D}/merge_s1/zAsign_GT
	mkdir -p ${merge_SV_D}/merge_s2/zAsign_GT
	mkdir -p ${merge_SV_D}/merge_s3/zAsign_GT

	group_prefix=${group_SV_D}/zChr/${GP}.stage2.SV.chr${c}.dot.addID.phas
	distance_prefix=${merge_SV_D}/SVanalyzer_distance/${GP}.stage2.SV.chr${c}.dot.addID.phas.merge
	## get SVanalyzer distance

	hash -r
	source ~/.bashrc
	conda deactivate
	conda activate svanalyzer

	svanalyzer merge \
	--ref ${ref_genome} \
	--variants ${group_prefix}.vcf.gz \
	--maxdist 2000 \
	--prefix ${distance_prefix}

	for s in 2 1 3;do

		asGT_prefix=${merge_SV_D}/merge_s${s}/zAsign_GT/${GP}.stage2.SV.chr${c}.dot.addID.phas.SV_s${s}
		merge_prefix=${merge_SV_D}/merge_s${s}/${GP}.stage2.SV.chr${c}.dot.addID.phas.SV_s${s}

		# get clustered SV id
		cat ${merge_prefix}.clustered.vcf | grep -v '^#'|awk '{print $8}'|awk -F'ClusterIDs=' '{print $2}'|awk -F';' '{print $1}' > ${asGT_prefix}.ClusterID
		cat ${merge_prefix}.clustered.vcf | grep '^#' > ${asGT_prefix}.header
		cat ${merge_prefix}.clustered.vcf | grep -v '^#'|awk -v OFS='\t' '{print $1,$2,$3,$4,$5}' > ${asGT_prefix}.siteInfo1
		cat ${merge_prefix}.clustered.vcf | grep -v '^#'|awk -v OFS='\t' '{print $6,$7,$8}' > ${asGT_prefix}.siteInfo2

		# get raw body for asign GT
		zcat ${group_prefix}.vcf.gz | grep ^c > ${group_prefix}.vcf.body

		for entry in $(cat ${asGT_prefix}.ClusterID);do
			python GRSB_asign_47_sample_gt.py \
			${entry} \
			${group_prefix}.vcf.body \
			${asGT_prefix}.gt
		done

		paste -d '\t' ${asGT_prefix}.siteInfo1 ${asGT_prefix}.siteInfo2 ${asGT_prefix}.gt > ${asGT_prefix}.body && cat ${asGT_prefix}.header ${asGT_prefix}.body > ${asGT_prefix}.asigned-GT.vcf && rm ${group_prefix}.vcf.body

		bcftools sort -Oz --temp-dir ${asGT_prefix} -o ${asGT_prefix}.asigned-GT.sort.vcf.gz ${asGT_prefix}.asigned-GT.vcf && tabix -f -p vcf ${asGT_prefix}.asigned-GT.sort.vcf.gz

	done
done




dipBED_D=dipcal_high_conficant_bed

for GP in PUB56 CN62 PUB1KCP111;do

	SNP_stg2_F=${base_D}/SNP/${GP}/group/zChr/${GP}.stage2.SNP.chr${c}.zero.addID.phas.vcf.gz
	SNP_stg3_F=${base_D}/SNP/${GP}/group/zChr/${GP}.stage3.SNP.chr${c}.zero.addID.phas.vcf.gz
	InDel_stg2_F=${base_D}/InDel/${GP}/group/zChr/${GP}.stage2.InDel.chr${c}.zero.addID.phas.vcf.gz
	SV_stg2_F=${base_D}/SV/${GP}/merge/merge_s2/zAsign_GT/${GP}.stage2.SV.chr${c}.dot.addID.phas.SV_s2.asigned-GT.sort.vcf.gz
	
	SNP_stg2_scatter_D=${base_D}/SNP/${GP}/dotZero/zDenseScatter
	SNP_stg2_neat_D=${base_D}/SNP/${GP}/dotZero

	SNP_stg3_scatter_D=${base_D}/SNP/${GP}/dotZero/zDenseScatter
	SNP_stg3_neat_D=${base_D}/SNP/${GP}/dotZero

	InDel_stg2_scatter_D=${base_D}/InDel/${GP}/dotZero/zDenseScatter
	InDel_stg2_neat_D=${base_D}/InDel/${GP}/dotZero

	SV_stg2_scatter_D=${base_D}/SV/${GP}/dotZero/zDenseScatter
	SV_stg2_neat_D=${base_D}/SV/${GP}/dotZero

	mkdir -p ${SNP_stg2_scatter_D}
	mkdir -p ${SNP_stg3_scatter_D}
	mkdir -p ${InDel_stg2_scatter_D}
	mkdir -p ${SV_stg2_scatter_D}

	SNP_stg2_neat_prefix=${SNP_stg2_neat_D}/${GP}.chr${c}.stage2.SNP
	SNP_stg3_neat_prefix=${SNP_stg3_neat_D}/${GP}.chr${c}.stage3.SNP
	InDel_stg2_neat_prefix=${InDel_stg2_neat_D}/${GP}.chr${c}.stage2.InDel
	SV_stg2_neat_prefix=${SV_stg2_neat_D}/${GP}.chr${c}.stage2.SV

	# subset and bedding each sample with non-var site and the dipcall high confident BED file
	for SID in $(cat ${infos_D}/SID_gender_seqTYPE_superPOP_subPOP_inOKG_SOURCE.CCS.${GP}|awk '{print $1}');do

		dipBED_F=${dipBED_D}/${SID}.dip.bed

		SNP_stg2_scatter_F=${SNP_stg2_scatter_D}/${SID}.${GP}.chr${c}.stage2.SNP.vcf.gz
		SNP_stg3_scatter_F=${SNP_stg3_scatter_D}/${SID}.${GP}.chr${c}.stage3.SNP.vcf.gz
		InDel_stg2_scatter_F=${InDel_stg2_scatter_D}/${SID}.${GP}.chr${c}.stage2.InDel.vcf.gz
		SV_stg2_scatter_F=${SV_stg2_scatter_D}/${SID}.${GP}.chr${c}.stage2.SV.vcf.gz

		bcftools view -Oz -s ${SID} ${SNP_stg2_F} |bcftools view -Oz -T ${dipBED_F} -o ${SNP_stg2_scatter_F} && tabix -p vcf ${SNP_stg2_scatter_F}
		bcftools view -Oz -s ${SID} ${SNP_stg3_F} |bcftools view -Oz -T ${dipBED_F} -o ${SNP_stg3_scatter_F} && tabix -p vcf ${SNP_stg3_scatter_F}
		bcftools view -Oz -s ${SID} ${InDel_stg2_F} |bcftools view -Oz -T ${dipBED_F} -o ${InDel_stg2_scatter_F} && tabix -p vcf ${InDel_stg2_scatter_F}
		bcftools view -Oz -s ${SID} ${SV_stg2_F} |bcftools view -Oz -T ${dipBED_F} -o ${SV_stg2_scatter_F} && tabix -p vcf ${SV_stg2_scatter_F}

		echo -e "${SNP_stg2_scatter_F}" >> ${SNP_stg2_neat_prefix}.VCFFILE
		echo -e "${SNP_stg3_scatter_F}" >> ${SNP_stg3_neat_prefix}.VCFFILE
		echo -e "${InDel_stg2_scatter_F}" >> ${InDel_stg2_neat_prefix}.VCFFILE
		echo -e "${SV_stg2_scatter_F}" >> ${SV_stg2_neat_prefix}.VCFFILE

	done
	
	# group back all samples with missing sites as dot, exclude sites missing > 10%
	bcftools merge -Oz -l ${SNP_stg2_neat_prefix}.VCFFILE -m none | bcftools annotate -Oz -x INFO | bcftools view -Ov -i 'F_MISSING<0.5' -o ${SNP_stg2_neat_prefix}.cr50.vcf
	bcftools merge -Oz -l ${SNP_stg3_neat_prefix}.VCFFILE -m none | bcftools annotate -Oz -x INFO | bcftools view -Ov -i 'F_MISSING<0.5' -o ${SNP_stg3_neat_prefix}.cr50.vcf
	bcftools merge -Oz -l ${InDel_stg2_neat_prefix}.VCFFILE -m none | bcftools annotate -Oz -x INFO | bcftools view -Ov -i 'F_MISSING<0.5' -o ${InDel_stg2_neat_prefix}.cr50.vcf
	bcftools merge -Oz -l ${SV_stg2_neat_prefix}.VCFFILE -m none | bcftools annotate -Oz -x INFO | bcftools view -Ov -i 'F_MISSING<0.5' -o ${SV_stg2_neat_prefix}.cr50.vcf

	# bgzip sort tabix
	bgzip ${SNP_stg2_neat_prefix}.cr50.zero.vcf
	bgzip ${SNP_stg3_neat_prefix}.cr50.zero.vcf
	bgzip ${InDel_stg2_neat_prefix}.cr50.zero.vcf
	bgzip ${SV_stg2_neat_prefix}.cr50.zero.vcf

	bcftools sort -Oz --temp-dir ${SNP_stg2_neat_D} -o ${SNP_stg2_neat_prefix}.cr50.zero.sort.vcf.gz ${SNP_stg2_neat_prefix}.cr50.zero.vcf.gz
	bcftools sort -Oz --temp-dir ${SNP_stg3_neat_D} -o ${SNP_stg3_neat_prefix}.cr50.zero.sort.vcf.gz ${SNP_stg3_neat_prefix}.cr50.zero.vcf.gz
	bcftools sort -Oz --temp-dir ${InDel_stg2_neat_D} -o ${InDel_stg2_neat_prefix}.cr50.zero.sort.vcf.gz ${InDel_stg2_neat_prefix}.cr50.zero.vcf.gz
	bcftools sort -Oz --temp-dir ${SV_stg2_neat_D} -o ${SV_stg2_neat_prefix}.cr50.zero.sort.vcf.gz ${SV_stg2_neat_prefix}.cr50.zero.vcf.gz

	tabix -f -p vcf ${SNP_stg2_neat_prefix}.cr50.zero.sort.vcf.gz && rm ${SNP_stg2_neat_prefix}.cr50.vcf ${SNP_stg2_neat_prefix}.cr50.zero.vcf.gz
	tabix -f -p vcf ${SNP_stg3_neat_prefix}.cr50.zero.sort.vcf.gz && rm ${SNP_stg3_neat_prefix}.cr50.vcf ${SNP_stg3_neat_prefix}.cr50.zero.vcf.gz
	tabix -f -p vcf ${InDel_stg2_neat_prefix}.cr50.zero.sort.vcf.gz && rm ${InDel_stg2_neat_prefix}.cr50.vcf ${InDel_stg2_neat_prefix}.cr50.zero.vcf.gz
	tabix -f -p vcf ${SV_stg2_neat_prefix}.cr50.zero.sort.vcf.gz && rm ${SV_stg2_neat_prefix}.cr50.vcf ${SV_stg2_neat_prefix}.cr50.zero.vcf.gz

	# addtional process for stage 3 snp
	target_F=${SNP_stg3_neat_prefix}.cr50.zero.sort.vcf.gz
	tmp_F=${SNP_stg3_neat_prefix}.cr50.zero.sort.tmp.vcf.gz
	OKG3202_F=/storage/yangjianLab/sharedata/LD_reference/1000G/1KGP3_3202_GRCh38/1kGP_high_coverage_Illumina.chr${c}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
	
	bcftools isec -Oz -p ${SNP_stg3_neat_D}/${GP}_chr${c} ${target_F} ${OKG3202_F}
	
	zcat ${SNP_stg3_neat_D}/${GP}_chr${c}/0000.vcf.gz | grep ^c | awk -v OFS='\t' '{print $1,$2}' | sort | uniq > ${SNP_stg3_neat_D}/${GP}_chr${c}/site.excl
	
	bcftools view -Oz -T ^${SNP_stg3_neat_D}/${GP}_chr${c}/site.excl ${target_F} -o ${tmp_F} && rm ${target_F} ${target_F}.tbi

	mv ${tmp_F} ${target_F} && tabix -p vcf ${target_F} && rm -rf ${SNP_stg3_neat_D}/${GP}_chr${c}

	################################################
	# addtional for s1 and s3 merge strategy
	for s in 1 3;do
		SV_stg2_F=${base_D}/SV/${GP}/merge/merge_s${s}/zAsign_GT/${GP}.stage2.SV.chr${c}.dot.addID.phas.SV_s${s}.asigned-GT.sort.vcf.gz
		
		SV_stg2_scatter_D=${base_D}/SV/${GP}/dotZero/zDenseScatter
		SV_stg2_neat_D=${base_D}/SV/${GP}/dotZero

		mkdir -p ${SV_stg2_scatter_D}

		SV_stg2_neat_prefix=${SV_stg2_neat_D}/${GP}.chr${c}.stage2.SV.s${s}
		
		# subset and bedding each sample with non-var site and the dipcall high confident BED file
		for SID in $(cat ${infos_D}/SID_gender_seqTYPE_superPOP_subPOP_inOKG_SOURCE.CCS.${GP}|awk '{print $1}');do

			dipBED_F=${dipBED_D}/${SID}.dip.bed

			SV_stg2_scatter_F=${SV_stg2_scatter_D}/${SID}.${GP}.chr${c}.stage2.SV.s${s}.vcf.gz

			bcftools view -Oz -s ${SID} ${SV_stg2_F} |bcftools view -Oz -T ${dipBED_F} -o ${SV_stg2_scatter_F} && tabix -f -p vcf ${SV_stg2_scatter_F}

			echo -e "${SV_stg2_scatter_F}" >> ${SV_stg2_neat_prefix}.VCFFILE

		done

		# group back all samples with missing sites as dot, exclude sites missing > 10%
		bcftools merge -Oz -l ${SV_stg2_neat_prefix}.VCFFILE -m none | bcftools annotate -Oz -x INFO | bcftools view -Ov -i 'F_MISSING<0.5' -o ${SV_stg2_neat_prefix}.cr50.vcf

		# make phase and set dot to zero
		cat ${SV_stg2_neat_prefix}.cr50.vcf | sed 's/\//|/g' | sed 's/\.|/0|/g' | sed 's/|\./|0/g' > ${SV_stg2_neat_prefix}.cr50.zero.vcf

		# bgzip sort tabix
		bgzip ${SV_stg2_neat_prefix}.cr50.zero.vcf

		bcftools sort -Oz --temp-dir ${SV_stg2_neat_D} -o ${SV_stg2_neat_prefix}.cr50.zero.sort.vcf.gz ${SV_stg2_neat_prefix}.cr50.zero.vcf.gz

		tabix -f -p vcf ${SV_stg2_neat_prefix}.cr50.zero.sort.vcf.gz && rm ${SV_stg2_neat_prefix}.cr50.vcf ${SV_stg2_neat_prefix}.cr50.zero.vcf.gz

	done
done



for GP in HGSVC2_HPRCY1 ;do
	# concat
	mkdir -p ${base_D}/panels/${GP}/vcf/M23s5F

	M23s5F=${base_D}/panels/${GP}/vcf/M23s5F/${GP}.chr${c}.M23s5F.vcf.gz

	SNPo=${base_D}/merged_sample/SNP/${GP}/dotZero/${GP}.chr${c}.stage2.SNP.cr50.zero.sort.vcf.gz
	SNPr=${base_D}/merged_sample/SNP/${GP}/dotZero/${GP}.chr${c}.stage3.SNP.cr50.zero.sort.vcf.gz
	InDel=${base_D}/merged_sample/InDel/${GP}/dotZero/${GP}.chr${c}.stage2.InDel.cr50.zero.sort.vcf.gz
	SVs5=${base_D}/merged_sample/SV/${GP}/dotZero/${GP}.chr${c}.stage2.SV.s5.cr50.zero.sort.vcf.gz
	SVsl5=${base_D}/merged_sample/SV/${GP}/dotZero/${GP}.chr${c}.stage2.SV.sl5.cr50.zero.sort.vcf.gz

	bcftools concat -Oz -a -d exact ${SNPr} ${SVs5} | bcftools sort -Oz --temp-dir ${base_D}/panels/${GP}/vcf/${GP}_chr${c} -o ${M23s5F} && tabix -p vcf ${M23s5F}
done

# These two parts (above and below) can be merged into one for-do-done.
for GP in HGSVC2_HPRCY1 ;do
	mkdir -p ${base_D}/panels/${GP}/m3vcf/M23s5F
	M23s5F_m3=${base_D}/panels/${GP}/m3vcf/M23s5F/${GP}.chr${c}.M23s5F
	Minimac3-omp --refHaps ${M23s5F} --processReference --chr chr${c} --cpus ${tr} --prefix ${M23s5F_m3}
done




