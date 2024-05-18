#!/bin/bash

c=${SLURM_ARRAY_TASK_ID}

M=M23s5F

# covar: type, length, maf, impuRsq, distance to TSS

vid_D=GENOTYPE/HGSVC2_HPRCY1/EUR
anno_D=ANNOTATE/HGSVC2_HPRCY1/SV_base_AnnotSV
gwa_D=GWAS/fastGWA_ori_T74/EUR/50
rsq_D=IMPUTATION/HGSVC2_HPRCY1
vep_D=ANNOTATE/HGSVC2_HPRCY1/SV_base_vep

em_D=ENRICHMENT

rm ${em_D}/covar_maf0.01_SV/chr${c}.vid.type.len.frq.ImpuRsq.dtss.dss.covar

for vid in $(cat ${vid_D}/HGSVC2_HPRCY1.${M}.EUR.chr${c}.SV.maf0.01.hwe-6.bim|awk '{print $2}');do

    pos=$(cat ${vid_D}/HGSVC2_HPRCY1.${M}.EUR.chr${c}.SV.maf0.01.hwe-6.bim|grep -w ${vid}|awk '{print $4}')
    ref=$(cat ${anno_D}/HGSVC2_HPRCY1.chr${c}.stage2.SV.s5.cr50.AnnotSV.tsv|grep -w ${vid}|head -n1|awk -F'\t' '{print $9}')
    if [ $ref ]; then ## few SVs of chr11 has no AnnotSV ## 
        alt=$(cat ${anno_D}/HGSVC2_HPRCY1.chr${c}.stage2.SV.s5.cr50.AnnotSV.tsv|grep -w ${vid}|head -n1|awk -F'\t' '{print $10}')
        len=$(cat ${anno_D}/HGSVC2_HPRCY1.chr${c}.stage2.SV.s5.cr50.AnnotSV.tsv|grep -w ${vid}|head -n1|awk -F'\t' '{print $5}')
        type=$(cat ${anno_D}/HGSVC2_HPRCY1.chr${c}.stage2.SV.s5.cr50.AnnotSV.tsv|grep -w ${vid}|head -n1|awk -F'\t' '{print $6}')
        freq=$(cat ${gwa_D}/50.Merged.EUR.maf0.01.hwe-6.chrAuto.fastGWA.cojoFMT.SV|grep -w ${vid}|awk '{print $4}')
        cpra=chr${c}:${pos}:${ref}:${alt}
        rsq=$(cat ${rsq_D}/EUR/HGSVC2_HPRCY1.${M}.EUR.chr${c}.ck9.impu.info_SV|grep -w ^${cpra}|awk '{print $7}')
        dtss=$(cat ${vep_D}/chr${c}.vep.anno.txt|grep -w ${vid}|grep TSSDistance|awk -F'TSSDistance=' '{print $2}'|sort -g|head -n1)
        if [ ! $dtss ];then dtss=-9;fi
        dss=$(cat ${anno_D}/HGSVC2_HPRCY1.chr${c}.stage2.SV.s5.cr50.AnnotSV.tsv|grep -w ${vid}|grep -w split|awk -F'\t' '{print $106}'|sort -g|head -n1)
        if [ ! $dss ];then dss=-9;fi
    fi
    echo -e "${vid}\t${vid}\t${type}\t${len}\t${freq}\t${rsq}\t${dtss}\t${dss}" >> ${em_D}/covar_maf0.01_SV/chr${c}.vid.type.len.frq.ImpuRsq.dtss.dss.covar
done


meta_D=meta_info/chromatin_state_EpiMap

cs_D=chromatin_state_processed_bed

sv_bed_F=HGSVC2_HPRCY1.chrAuto.stage2.SV.s5.cr50.AnnotSV.bed
ov_D=chromatin_state_ov


cell=$(cat ${meta_D}/MNEMONIC.sort.uniq_has_data.txt|sed -n "${sn}p"|awk '{print $1}')

bed_files=$(ls /storage/yangjianLab/sharedata/functional_annotation/EpiMap/chromatin_state/persample/${cell}_*)

for anot in $(cat ${meta_D}/chromatin_states.txt);do

zcat ${bed_files}|sed 's/\//_/g'|grep -w ${anot}|awk -v OFS='\t' '{print $1,$2,$3}'|bedtools sort|bedtools merge|awk -v OFS='\t' '{print $1,$2,$3,"'${anot}'"}' >\
${cs_D}/${cell}__${anot}.bed

# get if ov
bedtools annotate -i ${sv_bed_F} -files ${cs_D}/${cell}__${anot}.bed|awk -v OFS='\t' '{if($5>0) print $4,$4,1; else print $4,$4,0}' > ${ov_D}/${cell}__${anot}.chromatin_state_EpiMap.SV.ov

done


meta_D=meta_info/rbp_encode
rbp_D=RBP_eCLIP_processed_bed
ov_D=RBP_eCLIP_ov

for sn in $(seq 1 2);do

if [ $sn -eq 1 ]; then cell=HepG2; else cell=K562; fi

bed_files=$(cat ${meta_D}/rbp_bed_encode.info|grep -w ${cell}|awk -F'\t' '{print $1".bed.gz"}')


# get in ori Bed file dir
cd /storage/yangjianLab/baiweiyang/Reference_data/ENCODE_RBP_eCLIP_GRCh38

# get beds file with fc > 2 and -logP >5, and sort merge
zcat ${bed_files}|awk -v OFS='\t' '{if($7>2 && $8>5) print $1,$2,$3}'|bedtools sort|bedtools merge > ${rbp_D}/${cell}.allRBP.fc2.1e-5.sort.merge.bed

# get if ov
bedtools annotate -i ${sv_bed_F} -files ${rbp_D}/${cell}.allRBP.fc2.1e-5.sort.merge.bed|awk -v OFS='\t' '{if($5>0) print $4,$4,1; else print $4,$4,0}' > ${ov_D}/${cell}.allRBP.fc2.1e-5.sort.merge.SV.ov
done

ov_D=TAD_boundary_ov

ori_D=3D_Genome/hg38/zDoTADs

for sn in $(seq 1 40);do

cell=$(cat ${meta_D}/TAD_boundary_file.txt|sed -n "${sn}p"|awk '{print $1}')

# get if ov
bedtools annotate -i ${sv_bed_F} -files ${ori_D}/${cell}|awk -v OFS='\t' '{if($5>0) print $4,$4,1; else print $4,$4,0}' > ${ov_D}/${cell}.sort.merge.SV.ov
done



LD_D=LD_ref/EUR

cat ${LD_D}/Merged.EUR.chr${c}.maf0.01.hwe-6.bim|grep HGSV|awk '{print $2}' > ${LD_D}/Merged.EUR.chr${c}.maf0.01.hwe-6.svid

plink2 --bfile ${LD_D}/Merged.EUR.chr${c}.maf0.01.hwe-6 \
--extract ${LD_D}/Merged.EUR.chr${c}.maf0.01.hwe-6.svid \
--threads 4 --memory 20000 --indep-pairwise 500kb 1 0.2 \
--out ${LD_D}/zPrune/SV.EUR.chr${c}.maf0.01.hwe-6.500kb-1-0.2

#### R scripts
library(data.table)

args = commandArgs(trailingOnly = T)

ov_folder_path = as.character(args[1])
sigon_file = as.character(args[2])
covar_file = as.character(args[3])
output_file = as.character(args[4])
unrel_file = as.character(args[5])

signo = fread(sigon_file, col.names = c('ID','signo'), header = F)
covar = fread(covar_file, col.names = c('ID','IDc','type','len','frq','ImpuRsq','dtss','dss'), header = F)
unrel = fread(unrel_file, col.names = c('ID','IDu'), header = F)

file_list <- list.files(ov_folder_path)

summary_data=data.frame()

for (file_name in file_list) {
  ovon_file=paste0(ov_folder_path,'/',file_name)

  ovon = fread(ovon_file, col.names = c('ID','IDp','ovon'), header = F)

  combi = Reduce(function(x,y) merge(x,y,by="ID",all.x=F), list(ovon,covar,signo,unrel))

    mlr_Model = 'ovon ~ signo + type + len + frq + ImpuRsq + dtss + dss'
    fit<-glm(mlr_Model, data=combi, family = "binomial")
    summ_v=summary(fit)$coefficients
    beta=summ_v[2,1]
    se=summ_v[2,2]
    OR=exp(beta)
    OR_se=exp(beta)*se
    pvalue=summ_v[2,4]
    CI_1 = exp(beta + 1.96 * se)
    CI_2 = exp(beta - 1.96 * se)
    summ_v_stat=data.frame(expriment=file_name,beta=beta,se=se,P=pvalue,OR=OR,OR_se=OR_se,CI_1=CI_1,CI_2=CI_2)
    
  summary_data=rbind(summary_data,summ_v_stat)

}


write.table(summary_data, output_file, quote = F, sep = " ", row.names = F, col.names = T)

###
phenoDir=PHENOTYPE

for i in $(ls|grep _ov$);do

Rscript --slave --no-save /home/yangjianLab/baiweiyang/scripts/SV_Imputation_project_analysis/14.ENRICHMENT/4_logistic_regression/4_regression.R \
${i} \
sig_or_not_X/SV_sig-or-not.txt \
covar_maf0.01_SV/chrAuto.vid.type.len.frq.ImpuRsq.dtss.dss.covar \
${i}_regression/${i}.SV.logistic.lr \
covar_maf0.01_SV/zIndep.SV.id.two.col

done


