#!/bin/bash


c=${SLURM_ARRAY_TASK_ID}
ori_D=raw/GRCh38

plink2 \
--pfile ${ori_D}/ukb22828_c${c}_b0_v3_exFailed_b38 \
--maf 0.0001 \
--geno 0.01 \
--rm-dup exclude-all \
--snps-only \
--hwe 1e-6 \
--memory 80000 \
--threads 12 \
--output-chr chr26 \
--make-pgen \
--out ukb22828_c${c}_b0_v3_exFailed_b38.POP5.maf1e-4.geno01.hwe1e-6


for pop in EUR AFR AMR EAS SAS;do
out_D=${pop}
mkdir -p ${out_D}

plink2 \
--pfile ukb22828_c${c}_b0_v3_exFailed_b38.POP5.maf1e-4.geno01.hwe1e-6 \
--keep ${pop}.id \
--geno 0.01 \
--mac 3 \
--rm-dup exclude-all \
--snps-only \
--hwe 1e-6 \
--memory 80000 \
--threads 12 \
--output-chr chr26 \
--make-pgen \
--out ${out_D}/ukb22828_c${c}_b0_v3_exFailed_b38.${pop}.maf1e-4.geno01.hwe1e-6

done

for pop in EUR AFR AMR EAS SAS;do
cd ${pop}
mkdir -p ${pop}/zScatter/rsid_chunks
for c in $(seq 1 22);do
AL=$(wc -l ukb22828_c${c}_b0_v3_exFailed_b38.${pop}.maf1e-4.geno01.hwe1e-6.pvar|awk '{print $1}')
srt=1
nxt=15000
ct=1
while [ $nxt -le $AL ]; do
sed -n "${srt},${nxt}p" ukb22828_c${c}_b0_v3_exFailed_b38.${pop}.maf1e-4.geno01.hwe1e-6.pvar > zScatter/rsid_chunks/${c}.${ct}.rsid
let srt=(nxt-2500)
let nxt=(srt+15000)
let ct=(ct+1)
done
done
cd zScatter/rsid_chunks
ls *rsid > zFileName.txt
done


out_D=EUR

c=$(cat ${out_D}/zScatter/rsid_chunks/zFileName.txt|sed -n "${sn}p"|awk -F'.' '{print $1}')
vck=$(cat ${out_D}/zScatter/rsid_chunks/zFileName.txt|sed -n "${sn}p"|awk -F'.' '{print $2}')


mkdir -p ${out_D}/zScatter/geno_unphase

plink2 \
--pfile ${out_D}/ukb22828_c${c}_b0_v3_exFailed_b38.EUR.maf1e-4.geno01.hwe1e-6 \
--extract ${out_D}/zScatter/rsid_chunks/${c}.${vck}.rsid \
--memory 50000 \
--threads ${tr} \
--output-chr chr26 \
--export vcf \
--out ${out_D}/zScatter/geno_unphase/UKB.IMPUv3.EUR.chr${c}.QC.vck${vck} \
&& bgzip ${out_D}/zScatter/geno_unphase/UKB.IMPUv3.EUR.chr${c}.QC.vck${vck}.vcf \
&& tabix -p vcf ${out_D}/zScatter/geno_unphase/UKB.IMPUv3.EUR.chr${c}.QC.vck${vck}.vcf.gz

mkdir -p ${out_D}/zScatter/geno_phase

shapeit4.2 \
--input ${out_D}/zScatter/geno_unphase/UKB.IMPUv3.EUR.chr${c}.QC.vck${vck}.vcf.gz \
--map human_genetic_map/chr${c}.b38.gmap.gz \
--region chr${c} \
--thread ${tr} \
--pbwt-modulo 0.005 \
--output ${out_D}/zScatter/geno_phase/UKB.IMPUv3.EUR.chr${c}.QC.vck${vck}.shapeit4.vcf.gz

bcftools sort \
	-Oz \
	--temp-dir ${out_D}/zScatter/geno_phase/tmp-chr${c}-vck${vck} \
	--output ${out_D}/zScatter/geno_phase/UKB.IMPUv3.EUR.chr${c}.QC.vck${vck}.shapeit4.sort.vcf.gz \
	${out_D}/zScatter/geno_phase/UKB.IMPUv3.EUR.chr${c}.QC.vck${vck}.shapeit4.vcf.gz && \
tabix -p vcf ${out_D}/zScatter/geno_phase/UKB.IMPUv3.EUR.chr${c}.QC.vck${vck}.shapeit4.sort.vcf.gz && \
rm ${out_D}/zScatter/geno_phase/UKB.IMPUv3.EUR.chr${c}.QC.vck${vck}.shapeit4.vcf.gz ${out_D}/zScatter/geno_unphase/UKB.IMPUv3.EUR.chr${c}.QC.vck${vck}.vcf.gz*





for pop in EUR AFR AMR EAS SAS;do

bd=${pop}/zLigate_phased

mkdir -p ${bd}/zFilePath

vck_n=$(cat ${pop}/zScatter/rsid_chunks/zFileName.txt|grep -w ^${c}|wc -l)

for vck in $(seq 1 ${vck_n});do
echo -e "${pop}/zScatter/geno_phase/UKB.IMPUv3.${pop}.chr${c}.QC.vck${vck}.shapeit4.sort.vcf.gz" >>\
${bd}/zFilePath/chr${c}.vcks.filepath
done

bcftools concat \
	-Oz \
	--ligate \
	--compact-PS \
	--threads 12 \
	--file-list ${bd}/zFilePath/chr${c}.vcks.filepath \
	--output ${bd}/UKB.IMPUv3.${pop}.chr${c}.QC.shapeit4.sort.vcf.gz && \
tabix -p vcf ${bd}/UKB.IMPUv3.${pop}.chr${c}.QC.shapeit4.sort.vcf.gz

done


c=$(cat zChr_sck.txt|sed -n "${sn}p"|awk -F':' '{print $1}')
sk=$(cat zChr_sck.txt|sed -n "${sn}p"|awk -F':' '{print $2}')


# load envs
module load samtools/1.14
module load bcftools/1.14
module load htslib/1.14

mkdir -p ${bd}/zSChunks

bcftools view \
	-Oz \
	--samples-file ${bd}/zSample_Chunks_10K/zChunks_id/sChunks_${sk}.double \
	--force-samples \
	--output ${bd}/zSChunks/UKB.IMPUv3.EUR.chr${c}.QC.shapeit4.sort.Rsq03.ck${sk}.vcf.gz \
	${bd}/UKB.IMPUv3.EUR.chr${c}.QC.shapeit4.sort.Rsq03.vcf.gz && \
tabix -f -p vcf ${bd}/zSChunks/UKB.IMPUv3.EUR.chr${c}.QC.shapeit4.sort.Rsq03.ck${sk}.vcf.gz



k=${SLURM_ARRAY_TASK_ID}

M=M23s5F; GP=HGSVC2_HPRCY1

m3vcfDir=${GP}/m3vcf/${M}
outDir=${GP}

ARR_D=ARRAY

minimac4 \
	--refHaps ${m3vcfDir}/${GP}.chr${c}.${M}.m3vcf.gz \
	--haps ${ARR_D}/EUR/zLigate_phased/zSChunks/UKB.IMPUv3.EUR.chr${c}.QC.shapeit4.sort.Rsq03.ck${k}.vcf.gz \
	--prefix ${outDir}/EUR/${GP}.${M}.EUR.chr${c}.ck${k}.impu \
	--cpus 4 \
	--ignoreDuplicates \
	--vcfBuffer 155 \
	--allTypedSites \
	--ChunkLengthMb 20 \
	--ChunkOverlapMb 4 \
	--format GT,DS,GP \
	--log

tabix -p vcf ${outDir}/EUR/${GP}.${M}.EUR.chr${c}.ck${k}.impu.dose.vcf.gz

# extract SV
bcftools view \
	-Oz \
	--threads 4 \
	--include 'ILEN <= -50 || ILEN >= 50' \
	${outDir}/EUR/${GP}.${M}.EUR.chr${c}.ck${k}.impu.dose.vcf.gz |\
bcftools sort \
	-Oz \
	--temp-dir ${outDir}/EUR/tmp-EUR.chr${c}-${k} \
	--output ${outDir}/EUR/zSV_ck/${GP}.${M}.EUR.chr${c}.ck${k}.impu.dose.SV.sort.vcf.gz
tabix -p vcf ${outDir}/EUR/zSV_ck/${GP}.${M}.EUR.chr${c}.ck${k}.impu.dose.SV.sort.vcf.gz && rm ${outDir}/EUR/${GP}.${M}.EUR.chr${c}.ck${k}.impu.dose.vcf.gz*

cat ${outDir}/EUR/${GP}.${M}.EUR.chr${c}.ck${k}.impu.info | awk '{if(length($2)>49 || length($3)>49) print $0}' > ${outDir}/EUR/${GP}.${M}.EUR.chr${c}.ck${k}.impu.info_SV &&\
rm ${outDir}/EUR/${GP}.${M}.EUR.chr${c}.ck${k}.impu.info


minimac4 \
	--refHaps ${m3vcfDir}/${GP}.chr${c}.${M}.m3vcf.gz \
	--haps ${ARR_D}/${pop}/zLigate_phased/UKB.IMPUv3.${pop}.chr${c}.QC.shapeit4.sort.Rsq03.vcf.gz \
	--prefix ${outDir}/${pop}/${GP}.${M}.${pop}.chr${c}.impu \
	--cpus 4 \
	--ignoreDuplicates \
	--vcfBuffer 155 \
	--allTypedSites \
	--ChunkLengthMb 20 \
	--ChunkOverlapMb 4 \
	--format GT,DS,GP \
	--log

tabix -p vcf ${outDir}/${pop}/${GP}.${M}.${pop}.chr${c}.impu.dose.vcf.gz

# extract SV
bcftools view \
	-Oz \
	--threads 4 \
	--include 'ILEN <= -50 || ILEN >= 50' \
	${outDir}/${pop}/${GP}.${M}.${pop}.chr${c}.impu.dose.vcf.gz |\
bcftools annotate \
	-Oz \
	--set-id '%CHROM\_%POS\_SV' |\
bcftools sort \
	-Oz \
	--temp-dir ${outDir}/${pop}/tmp-${pop}.chr${c} \
	--output ${outDir}/${pop}/zSV/${GP}.${M}.${pop}.chr${c}.impu.dose.SV.reN.sort.vcf.gz
tabix -p vcf ${outDir}/${pop}/zSV/${GP}.${M}.${pop}.chr${c}.impu.dose.SV.reN.sort.vcf.gz && rm ${outDir}/${pop}/${GP}.${M}.${pop}.chr${c}.impu.dose.vcf.gz*

cat ${outDir}/${pop}/${GP}.${M}.${pop}.chr${c}.impu.info | awk '{if(length($2)>49 || length($3)>49) print $0}' > ${outDir}/${pop}/${GP}.${M}.${pop}.chr${c}.impu.info_SV &&\
rm ${outDir}/${pop}/${GP}.${M}.${pop}.chr${c}.impu.info


M=M23s5F

for GP in HGSVC2_HPRCY1; do

	for k in $(seq 0 45); do
		echo "${impu_D}/${GP}/EUR/zSV_ck/${GP}.${M}.EUR.chr${c}.ck${k}.impu.dose.SV.sort.vcf.gz" >> ${impu_D}/${GP}/EUR/zSV/${GP}.${M}.EUR.chr${c}.ckAll.impu.dose.SV.VCFFILE
	done

	# merge back SV
	bcftools merge \
		-Oz \
		--file-list ${impu_D}/${GP}/EUR/zSV/${GP}.${M}.EUR.chr${c}.ckAll.impu.dose.SV.VCFFILE \
		--info-rules AF:avg,MAF:avg,R2:avg \
		--threads ${tr} \
		--merge none |\
	bcftools annotate \
		-Oz \
		--threads ${tr} \
		--set-id '%CHROM\_%POS\_SV' |\
	bcftools sort \
		-Oz \
		--temp-dir ${impu_D}/${GP}/EUR/tmp-EUR.chr${c} \
		--output ${impu_D}/${GP}/EUR/zSV/${GP}.${M}.EUR.chr${c}.impu.dose.SV.reN.sort.vcf.gz
	tabix -p vcf ${impu_D}/${GP}/EUR/zSV/${GP}.${M}.EUR.chr${c}.impu.dose.SV.reN.sort.vcf.gz && rm ${impu_D}/${GP}/EUR/zSV_ck/${GP}.${M}.EUR.chr${c}.ck*.impu.dose.SV.sort.vcf.gz*

done



GP=HGSVC2_HPRCY1; M=M23s5F

IMPU_D=IMPUTATION; GENO_D=GENOTYPE

# get ukb own imputed data
for pop in EUR AFR AMR EAS SAS; do

inDir=/storage/yangjianLab/sharedata/uk_biobank/geno_imp_GRCh38
outDir=/storage/yangjianLab/baiweiyang/SV_Imputation_Project_final/GENOTYPE

plink2 \
--pfile ${inDir}/ukb22828_c${c}_b0_v3_exFailed_b38 \
--keep /storage/yangjianLab/sharedata/uk_biobank/UKB_Ancestry/${pop}.id \
--make-bed \
--memory 90000 \
--hwe 1e-6 \
--maf 0.01 \
--threads 2 \
--out ${outDir}/UOI_hg38/${pop}/UOI.${pop}.maf0.01.hwe-6.chr${c}

# rename UOI rsid to avoid dup id
mv ${outDir}/UOI_hg38/${pop}/UOI.${pop}.maf0.01.hwe-6.chr${c}.bim ${outDir}/UOI_hg38/${pop}/zORI_BIM/UOI.${pop}.maf0.01.hwe-6.chr${c}.bim

python /home/yangjianLab/baiweiyang/scripts/SV_Imputation_project_final/DEAL_GENOPHENO/prepare_genos/rename_UOI_bim.py \
${outDir}/UOI_hg38/${pop}/zORI_BIM/UOI.${pop}.maf0.01.hwe-6.chr${c}.bim \
${outDir}/UOI_hg38/${pop}/UOI.${pop}.maf0.01.hwe-6.chr${c}.bim \
${c}


# convert vcf to bfile
mkdir -p ${GENO_D}/${GP}/${pop}

plink2 \
--vcf ${IMPU_D}/${GP}/${pop}/zSV/${GP}.${M}.${pop}.chr${c}.impu.dose.SV.reN.sort.vcf.gz \
--memory 90000 \
--hwe 1e-6 \
--mac 1 \
--make-bed \
--threads 2 \
--out ${GENO_D}/${GP}/${pop}/${GP}.${M}.${pop}.chr${c}.SV.mac1.hwe-6

# rename OWN rsid to avoid dup id
python /home/yangjianLab/baiweiyang/scripts/SV_Imputation_project_analysis/02.GENOTYPE/1_prepare_geno_SGV_SV/1_prepare_geno_SGV_SV.py \
${GENO_D}/${GP}/${pop}/${GP}.${M}.${pop}.chr${c}.SV.mac1.hwe-6.bim \
/storage/yangjianLab/baiweiyang/SV_Imputation_Project_panel/Joint_HPRC_ccs47_HGSVC2_ccs9clr21_CKCG_ccs55_132/panels/${GP}/SV_ID/${M}/${GP}.chr${c}.${M}.SV.chr_pos_svid_type_ref_alt.tsv \
${GENO_D}/${GP}/${pop}/${GP}.${M}.${pop}.chr${c}.SV.mac1.hwe-6.bim_new

mv ${GENO_D}/${GP}/${pop}/${GP}.${M}.${pop}.chr${c}.SV.mac1.hwe-6.bim_new ${GENO_D}/${GP}/${pop}/${GP}.${M}.${pop}.chr${c}.SV.mac1.hwe-6.bim


# deal sample ID from SID_SID to SID\tSID
sed -i 's/^0\t//g' ${GENO_D}/${GP}/${pop}/${GP}.${M}.${pop}.chr${c}.SV.mac1.hwe-6.fam
sed -i 's/_/\t/g' ${GENO_D}/${GP}/${pop}/${GP}.${M}.${pop}.chr${c}.SV.mac1.hwe-6.fam


# get maf 0.01 sv for merge
plink2 \
--bfile ${GENO_D}/${GP}/${pop}/${GP}.${M}.${pop}.chr${c}.SV.mac1.hwe-6 \
--memory 90000 \
--maf 0.01 \
--make-bed \
--threads 2 \
--out ${GENO_D}/${GP}/${pop}/${GP}.${M}.${pop}.chr${c}.SV.maf0.01.hwe-6


# merge uoi and 0.01 sv
mkdir -p ${GENO_D}/${GP}_HRC_v3_hg38_merged/${pop}

plink \
--bfile ${GENO_D}/HRC_v3_hg38/${pop}/UOI.${pop}.maf0.01.hwe-6.chr${c} \
--bmerge ${GENO_D}/${GP}/${pop}/${GP}.${M}.${pop}.chr${c}.SV.maf0.01.hwe-6 \
--make-bed \
--memory 90000 \
--threads 2 \
--out ${GENO_D}/${GP}_HRC_v3_hg38_merged/${pop}/Merged.${pop}.chr${c}.maf0.01.hwe-6
done


# get EUR rd50k for LD ref
pop=EUR
mkdir -p ${GENO_D}/LD_ref/${pop}

plink2 \
--bfile ${GENO_D}/${GP}_HRC_v3_hg38_merged/EUR/Merged.${pop}.chr${c}.maf0.01.hwe-6 \
--keep ${GENO_D}/infos/unRel_0.05_UKB_GRM_HM3_QC_R2_09_MAF01_EUR.rd50K.grm.id \
--make-bed \
--memory 90000 \
--threads 2 \
--out ${GENO_D}/LD_ref/${pop}/Merged.${pop}.chr${c}.maf0.01.hwe-6



pop=$(cat pop4|sed -n "${pN}p"|awk '{print $1}')

rm -rf ${GRM_D}/GRM_POP4/zSV_SGV/zFILEPATH/${pop}.mgrm.2P.path

for vt in SV SGV; do

gcta64 \
--mbfile ${geno_D}/infos/zFILEPATH/Merged.${pop}.maf0.01.hwe-6.chrAuto.bfile \
--keep ${geno_D}/infos/zUnRel_ID/UnRel.0.05.Merged.${pop}.maf0.01.hwe-6.chrAuto.grm.id \
--extract ${geno_D}/infos/zVAR_ID/Merged.${pop}.chrAuto.maf0.01.hwe-6.${vt}.vid \
--maf 0.01 \
--make-grm \
--thread-num 2 \
--out ${GRM_D}/GRM_POP4/zSV_SGV/${vt}.unRel0.05.${pop}.maf0.01.hwe-6.chrAuto

echo -e "${GRM_D}/GRM_POP4/zSV_SGV/${vt}.unRel0.05.${pop}.maf0.01.hwe-6.chrAuto" >> ${GRM_D}/GRM_POP4/zSV_SGV/zFILEPATH/${pop}.mgrm.2P.path

done

for vt in SV SV_sole; do
gcta64 \
--mbfile ${geno_D}/infos/zFILEPATH/Merged.${pop}.maf0.01.hwe-6.chrAuto.bfile \
--keep ${GRM_D}/GRM_${pop}/GRM_chunks/zSample/${pop}_sck${sk} \
--extract ${geno_D}/infos/zVAR_ID/Merged.${pop}.chrAuto.maf0.01.hwe-6.${vt}.vid \
--maf 0.01 \
--make-grm \
--thread-num 12 \
--out ${GRM_D}/GRM_${pop}/GRM_chunks/${pop}.maf0.01.hwe-6.chrAuto.${vt}.sck${sk}
done

