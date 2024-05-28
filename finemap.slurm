#!/bin/bash

tr=2; sn=${SLURM_ARRAY_TASK_ID}

# specify trait type
if [ $sn -gt 2088 ]; then cob=C; else cob=B; fi
filePath_F=/storage/yangjianLab/baiweiyang/SV_Imputation_Project_final/PHENOTYPE/infos/both.final.filePath
phe=$(cat ${filePath_F} | sed -n "${sn}p" | awk -F'/' '{print $8}' | awk -F'.pheno' '{print $1}')

FINEMAP_D=FINEMAP/bench # this is for bench dir
GWAS_D=GWAS
GENO_D=GENOTYPE

for pop in EUR ; do

    if [ -s "${GWAS_D}/cojo_select_bench/${pop}/${phe}/association/chrAuto.SV.cojo-slct.vid" ]; then

        ## for fine-mapp
        mkdir -p ${FINEMAP_D}/FM/${pop}/${phe}
        mkdir -p ${FINEMAP_D}/FM/${pop}/${phe}/z_file
        mkdir -p ${FINEMAP_D}/FM/${pop}/${phe}/geno_file
        mkdir -p ${FINEMAP_D}/FM/${pop}/${phe}/ld_file
        mkdir -p ${FINEMAP_D}/FM/${pop}/${phe}/output_file
        mkdir -p ${FINEMAP_D}/FM/${pop}/${phe}/stat_pip

        # get z.all
        echo -e "rsid chromosome position allele1 allele2 maf beta se" > ${FINEMAP_D}/FM/${pop}/${phe}/${phe}.all.z

        if [ ${cob} == "B" ]; then
            cat ${GWAS_D}/fastGWA_ori_bench/${pop}/${phe}/${phe}.Merged.${pop}.maf0.01.hwe-6.chrAuto.fastGWA | sed '1d' |\
            awk '{if($7<=0.5) print $2,$1,$3,$4,$5,$7,$11,$12; else print $2,$1,$3,$4,$5,(1-$7),$11*-1,$12}' >> ${FINEMAP_D}/FM/${pop}/${phe}/${phe}.all.z
        else
            cat ${GWAS_D}/fastGWA_ori_bench/${pop}/${phe}/${phe}.Merged.${pop}.maf0.01.hwe-6.chrAuto.fastGWA | sed '1d' |\
            awk '{if($7<=0.5) print $2,$1,$3,$4,$5,$7,$8,$9; else print $2,$1,$3,$4,$5,(1-$7),$8*-1,$9}' >> ${FINEMAP_D}/FM/${pop}/${phe}/${phe}.all.z
        fi

        # gather cojo-slct var cordinates flanking 500K, in total 1M 
        for slt in $(cat ${GWAS_D}/cojo_select_bench/${pop}/${phe}/association/chrAuto.SV.cojo-slct.vid)
        do
            chm=$(cat ${FINEMAP_D}/FM/${pop}/${phe}/${phe}.all.z | grep -w "${slt}" |awk '{print $2}')
            pos=$(cat ${FINEMAP_D}/FM/${pop}/${phe}/${phe}.all.z | grep -w "${slt}" |awk '{print $3}')
            let ups=(pos-500000)
            let downs=(pos+500000)
            if (( $ups <= 0 ))
            then
                ups=0
            fi
            echo ${chm}-${ups}-${downs}-${slt} >> ${FINEMAP_D}/FM/${pop}/${phe}/cordi_regions.txt
        done

        # extract z chunk based on final cojo-slct var cordinates, and select exactly overlapped SNP and SV in genotype data (LD ref).

        for fst in $(cat ${FINEMAP_D}/FM/${pop}/${phe}/cordi_regions.txt)
        do
            chm=$(echo ${fst} | awk -F'-' '{print $1}')
            srt=$(echo ${fst} | awk -F'-' '{print $2}')
            end=$(echo ${fst} | awk -F'-' '{print $3}')
            vid=$(echo ${fst} | awk -F'-' '{print $4}')
            echo -e "rsid chromosome position allele1 allele2 maf beta se" > ${FINEMAP_D}/FM/${pop}/${phe}/z_file/${phe}.${vid}.${chm}-${srt}-${end}.z
            cat ${FINEMAP_D}/FM/${pop}/${phe}/${phe}.all.z | sed '1d' | awk '{if($2=='"$chm"') print $0}' |\
            awk '{if($3>'"$srt"' && $3<'"$end"') print $0}' >> ${FINEMAP_D}/FM/${pop}/${phe}/z_file/${phe}.${vid}.${chm}-${srt}-${end}.z
            cat ${FINEMAP_D}/FM/${pop}/${phe}/z_file/${phe}.${vid}.${chm}-${srt}-${end}.z | sed '1d' | awk '{print $1}' > ${FINEMAP_D}/FM/${pop}/${phe}/z_file/${phe}.${vid}.${chm}-${srt}-${end}.rsid

            # prepare geno and LD files
            plink2 \
            --bfile ${GENO_D}/HGSVC2_HPRCY1_HRC_v3_hg38_merged/${pop}/Merged.${pop}.chr${chm}.maf0.01.hwe-6 \
            --keep ${GENO_D}/infos/zUnRel_ID/UnRel.0.05.Merged.${pop}.maf0.01.hwe-6.chrAuto.grm.id \
            --extract ${FINEMAP_D}/FM/${pop}/${phe}/z_file/${phe}.${vid}.${chm}-${srt}-${end}.rsid \
            --make-bed \
            --threads 2 \
            --memory 30000 \
            --out ${FINEMAP_D}/FM/${pop}/${phe}/geno_file/${phe}.${vid}.${chm}-${srt}-${end}.unRel.geno

            plink \
            --bfile ${FINEMAP_D}/FM/${pop}/${phe}/geno_file/${phe}.${vid}.${chm}-${srt}-${end}.unRel.geno \
            --r square spaces \
            --threads 2 \
            --memory 30000 \
            --out ${FINEMAP_D}/FM/${pop}/${phe}/ld_file/${phe}.${vid}.${chm}-${srt}-${end}.unRel.plkLD
        done


        # prepare master file
        echo 'z;ld;snp;config;cred;log;n_samples' > ${FINEMAP_D}/FM/${pop}/${phe}/${phe}.cojo.master
        for fst in $(cat ${FINEMAP_D}/FM/${pop}/${phe}/cordi_regions.txt)
        do
            chm=$(echo ${fst} | awk -F'-' '{print $1}')
            srt=$(echo ${fst} | awk -F'-' '{print $2}')
            end=$(echo ${fst} | awk -F'-' '{print $3}')
            vid=$(echo ${fst} | awk -F'-' '{print $4}')
            zFile=${FINEMAP_D}/FM/${pop}/${phe}/z_file/${phe}.${vid}.${chm}-${srt}-${end}.z
            ldFile=${FINEMAP_D}/FM/${pop}/${phe}/ld_file/${phe}.${vid}.${chm}-${srt}-${end}.unRel.plkLD.ld
            snpFile=${FINEMAP_D}/FM/${pop}/${phe}/output_file/${phe}.${vid}.${chm}-${srt}-${end}.snp
            configFile=${FINEMAP_D}/FM/${pop}/${phe}/output_file/${phe}.${vid}.${chm}-${srt}-${end}.config
            credFile=${FINEMAP_D}/FM/${pop}/${phe}/output_file/${phe}.${vid}.${chm}-${srt}-${end}.cred
            logFile=${FINEMAP_D}/FM/${pop}/${phe}/output_file/${phe}.${vid}.${chm}-${srt}-${end}.log
            n_samples=$(cat ${GWAS_D}/fastGWA_ori_bench/${pop}/${phe}/${phe}.Merged.${pop}.maf0.01.hwe-6.chrAuto.fastGWA | sed -n '2p' | awk '{print $6}')
            echo "${zFile};${ldFile};${snpFile};${configFile};${credFile};${logFile};${n_samples}" >> ${FINEMAP_D}/FM/${pop}/${phe}/${phe}.cojo.master
        done

        # Run finemap
        finemap_v1.4.1_x86_64 \
        --sss \
        --in-files ${FINEMAP_D}/FM/${pop}/${phe}/${phe}.cojo.master \
        --n-causal-snps 10 \
        --n-threads ${tr} \
        --log

    fi
done

