#!/bin/bash


sn=${SLURM_ARRAY_TASK_ID}

data_set=HGSVC2_HPRCY1

SID=$(cat ${data_set}_assembly/zSample.id.gender| sed -n "${sn}p"|awk '{print $1}')
SEX=$(cat ${data_set}_assembly/zSample.id.gender| sed -n "${sn}p"|awk '{print $2}')

ref_F=hg38.no_alt.fa
h1_F=${data_set}_assembly/${SID}/h1.fasta
h2_F=${data_set}_assembly/${SID}/h2.fasta
PAR_F=hs38.PAR.bed
vcfOut=CKCG_plus3_58_variant/dipcall/ori/${SID}

if [ $SEX == 'male' ];then
run-dipcall -u -t ${tr} -x ${PAR_F} ${SID} ${ref_F} ${h1_F} ${h2_F} > ${vcfOut}/${SID}.mak
cd ${vcfOut}
make -j2 -f ${vcfOut}/${SID}.mak
else
run-dipcall -u -t ${tr} ${SID} ${ref_F} ${h1_F} ${h2_F} > ${vcfOut}/${SID}.mak
cd ${vcfOut}
make -j2 -f ${vcfOut}/${SID}.mak
fi

vcfOut=CKCG_plus3_58_variant/pav/ori/${SID}

cp config.json ${vcfOut}

cd ${vcfOut}

snakemake -c${tr} -s Snakefile  pav_${SID}.vcf.gz

bamOut=${data_set}_alignment_unimap_GRCh38
vcfOut=CKCG_plus3_58_variant/svim/ori/${SID}

unimap -a -x asm5 --cs -r2k -t ${tr} ${ref_F} ${h1_F} > ${bamOut}/${SID}.h1.sam
unimap -a -x asm5 --cs -r2k -t ${tr} ${ref_F} ${h2_F} > ${bamOut}/${SID}.h2.sam
samtools sort -m8G -@ ${tr} -o ${bamOut}/${SID}.h1.sorted.bam ${bamOut}/${SID}.h1.sam
samtools sort -m8G -@ ${tr} -o ${bamOut}/${SID}.h2.sorted.bam ${bamOut}/${SID}.h2.sam
samtools index -@ ${tr} ${bamOut}/${SID}.h1.sorted.bam
samtools index -@ ${tr} ${bamOut}/${SID}.h2.sorted.bam

svim-asm diploid --min_sv_size 50 ${vcfOut} ${bamOut}/${SID}.h1.sorted.bam ${bamOut}/${SID}.h2.sorted.bam ${ref_F}


# preprocess
bcftools reheader --header ${qcOut}/header.txt --threads ${tr} ${vcfOri_F} |\
bcftools norm -Oz --threads ${tr} --fasta-ref ${ref_F} --check-ref s --multiallelics - |\
bcftools norm -Oz --threads ${tr} --rm-dup exact |\
bcftools view -Oz --threads ${tr} --exclude 'ALT="*,." | REF="*,."' |\
bcftools view -Oz --threads ${tr} --include 'CHROM="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"' |\
bcftools annotate -Oz --set-id "dip_${SID}_%CHROM\_%POS" |\
bcftools sort -Oz --temp-dir ${qcOut}/tmp-${SID} --output ${qcOut}/${SID}.dip.bi.norm.exSND.chrAuto.addID.sort.vcf.gz
tabix -p vcf ${qcOut}/${SID}.dip.bi.norm.exSND.chrAuto.addID.sort.vcf.gz

# filter 
bcftools view -Oz --threads ${tr} --targets-file ^${model_region_F} ${qcOut}/${SID}.dip.bi.norm.exSND.chrAuto.addID.sort.vcf.gz |\
bcftools view -Oz --threads ${tr} --targets-file ${cf_bed_F} |\
bcftools sort -Oz --temp-dir ${qcOut}/tmp-${SID} --output ${qcOut}/${SID}.dip.bi.norm.exSND.chrAuto.addID.exCEN.inConFi.sort.vcf.gz
tabix -p vcf ${qcOut}/${SID}.dip.bi.norm.exSND.chrAuto.addID.exCEN.inConFi.sort.vcf.gz


for ff in addID addID.exCEN.inConFi;do
# split into SNP
bcftools view -Oz --threads ${tr} --types snps ${qcOut}/${SID}.dip.bi.norm.exSND.chrAuto.${ff}.sort.vcf.gz |\
bcftools sort -Oz --temp-dir ${qcOut}/zSNP/tmp-${SID} --output ${qcOut}/zSNP/${SID}.dip.bi.norm.exSND.chrAuto.${ff}.SNP.sort.vcf.gz
tabix -p vcf ${qcOut}/zSNP/${SID}.dip.bi.norm.exSND.chrAuto.${ff}.SNP.sort.vcf.gz

# split into InDel
bcftools view -Oz --threads ${tr} --include 'ILEN >= -'"${n1}"' && ILEN <= '"${n1}"'' ${qcOut}/${SID}.dip.bi.norm.exSND.chrAuto.${ff}.sort.vcf.gz |\
bcftools sort -Oz --temp-dir ${qcOut}/zInDel/tmp-${SID} --output ${qcOut}/zInDel/${SID}.dip.bi.norm.exSND.chrAuto.${ff}.InDel.sort.vcf.gz
tabix -p vcf ${qcOut}/zInDel/${SID}.dip.bi.norm.exSND.chrAuto.${ff}.InDel.sort.vcf.gz

# split into SV
bcftools view -Oz --threads ${tr} --include 'ILEN <= -'"${n2}"' || ILEN >= '"${n2}"'' ${qcOut}/${SID}.dip.bi.norm.exSND.chrAuto.${ff}.sort.vcf.gz |\
bcftools sort -Oz --temp-dir ${qcOut}/zSV/tmp-${SID} --output ${qcOut}/zSV/${SID}.dip.bi.norm.exSND.chrAuto.${ff}.SV.sort.vcf.gz
tabix -p vcf ${qcOut}/zSV/${SID}.dip.bi.norm.exSND.chrAuto.${ff}.SV.sort.vcf.gz
done


deal=bi.norm.exSND.chrAuto.addID.exCEN.inConFi

SID=$(cat ${data_set}_assembly/zSample.id.gender| sed -n "${sn}p"|awk '{print $1}')

for vt in SNP InDel;do
mkdir -p integ/${vt}/${SID}
outD=integ/${vt}/${SID}

dip_f=dip/qc/${SID}/z${vt}/${SID}.dip.${deal}.${vt}.sort.vcf.gz
pav_f=pav/qc/${SID}/z${vt}/${SID}.pav.${deal}.${vt}.sort.vcf.gz


echo "${dip_f}" > ${outD}/${SID}-unimap-${vt}_dip_pav.vcfFile
echo "${pav_f}" >> ${outD}/${SID}-unimap-${vt}_dip_pav.vcfFile


bcftools merge -Oz --file-list ${outD}/${SID}-unimap-${vt}_dip_pav.vcfFile --merge none --force-samples |\
bcftools view -Oz --exclude 'GT="./."' --output ${outD}/${SID}_unimap_match-two.${deal}.${vt}.sort.pre.vcf.gz

zcat ${outD}/${SID}_unimap_match-two.${deal}.${vt}.sort.pre.vcf.gz |grep -v '^#' > ${outD}/${SID}_unimap_match-two.${deal}.${vt}.noheader.vcf
zcat ${outD}/${SID}_unimap_match-two.${deal}.${vt}.sort.pre.vcf.gz |grep '^#' > ${outD}/${SID}_unimap_match-two.${deal}.${vt}.header

sed -i '$d' ${outD}/${SID}_unimap_match-two.${deal}.${vt}.header
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${SID}" >> ${outD}/${SID}_unimap_match-two.${deal}.${vt}.header

python3 Match-SGV.py \
${outD}/${SID}_unimap_match-two.${deal}.${vt}.noheader.vcf \
${outD}/${SID}_unimap_match-two.${deal}.${vt}.noheader.pre1.vcf

cat ${outD}/${SID}_unimap_match-two.${deal}.${vt}.noheader.pre1.vcf |awk -v OFS='\t' '{print $3,$4,$5,$6,$7,$8,$9,$10,$1,$2}' > \
${outD}/${SID}_unimap_match-two.${deal}.${vt}.noheader.pre2.vcf

sed -i 's/\//|/g' ${outD}/${SID}_unimap_match-two.${deal}.${vt}.noheader.pre2.vcf


cat ${outD}/${SID}_unimap_match-two.${deal}.${vt}.header ${outD}/${SID}_unimap_match-two.${deal}.${vt}.noheader.pre2.vcf |\
bcftools sort -Oz --temp-dir ${outD}/tmp --output ${outD}/${SID}_unimap_match-two.${deal}.${vt}.sort.vcf.gz

tabix -p vcf ${outD}/${SID}_unimap_match-two.${deal}.${vt}.sort.vcf.gz

done


ref_F=hg38.no_alt.fa
deal=bi.norm.exSND.chrAuto.addID.exCEN.inConFi
ff=addID.exCEN.inConFi
outD=${data_set}_variant/integ/SV/${SID}

for cls in pav-dip-svim;do
cl_1=$(echo ${cls}|awk -F'-' '{print $1}')
cl_2=$(echo ${cls}|awk -F'-' '{print $2}')
cl_3=$(echo ${cls}|awk -F'-' '{print $3}')

cl_1_qc_F=${data_set}_variant/${cl_1}/qc/${SID}/zSV/${SID}.${cl_1}.${deal}.SV.sort.vcf.gz

## exclude use ID not chr_pos
echo -e "." > ${outD}/zFinal_SV_each_caller/${SID}_${cl_1}_unimap_uniq-SV.${deal}.txt

sort ${outD}/${cl_1}-${cl_2}-unimap-SV-${SID}-${ff}.falsepositives.vcf \
${outD}/${cl_1}-${cl_3}-unimap-SV-${SID}-${ff}.falsepositives.vcf |uniq -d|awk '{print $3}' >> ${outD}/zFinal_SV_each_caller/${SID}_${cl_1}_unimap_uniq-SV.${deal}.txt

bcftools view -Oz --include ID!=@${outD}/zFinal_SV_each_caller/${SID}_${cl_1}_unimap_uniq-SV.${deal}.txt ${cl_1_qc_F} |\
bcftools sort -Oz --temp-dir ${outD}/zFinal_SV_each_caller/tmp-${cl_1} --output ${outD}/zFinal_SV_each_caller/${SID}_${cl_1}_unimap_match-two.${deal}.SV.sort.vcf.gz
tabix -p vcf ${outD}/zFinal_SV_each_caller/${SID}_${cl_1}_unimap_match-two.${deal}.SV.sort.vcf.gz
done


