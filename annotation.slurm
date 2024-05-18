#!/bin/bash


GP=HGSVC2_HPRCY1

ANNOT_D=ANNOTATE
SV_vcf_D=${GP}/dotZero


for c in $(seq 1 22);do
# anno SV type and len
truvari anno svinfo -i ${SV_vcf_D}/${GP}.chr${c}.stage2.SV.s5.cr50.zero.sort.vcf.gz -o ${ANNOT_D}/${GP}/SV_base_sets/${GP}.chr${c}.stage2.SV.s5.cr50.addTYPE.addLEN.vcf

# annotSV
${ANNOTSV}/bin/AnnotSV \
-SVminSize 0 \
-SVinputFile ${ANNOT_D}/${GP}/SV_base_sets/${GP}.chr${c}.stage2.SV.s5.cr50.addTYPE.addLEN.vcf \
-outputFile ${ANNOT_D}/${GP}/SV_base_AnnotSV/${GP}.chr${c}.stage2.SV.s5.cr50.AnnotSV.tsv

# anno INFO and drop genotypes
bcftools +fill-tags ${ANNOT_D}/${GP}/SV_base_sets/${GP}.chr${c}.stage2.SV.s5.cr50.addTYPE.addLEN.vcf |\
bcftools view -Ov --drop-genotypes --output ${ANNOT_D}/${GP}/SV_base_sets/${GP}.chr${c}.stage2.SV.s5.cr50.addTYPE.addLEN.addINFO.dropGT.vcf \
&& rm ${ANNOT_D}/${GP}/SV_base_sets/${GP}.chr${c}.stage2.SV.s5.cr50.addTYPE.addLEN.vcf

done


