#!/bin/bash

for spe in main_jt;do 

aimD=zzSim_spe_${spe}
mkdir -p ${aimD}

tr=1
n=${SLURM_ARRAY_TASK_ID}

rep_N=7500

s=$((${n}*5-5+1))
e=$((${n}*5))


tmpDir=/data/jy_${spe}/${n}; rm -rf ${tmpDir}; mkdir -p ${tmpDir}

mkdir -p ${tmpDir}/cv_pool
mkdir -p ${tmpDir}/output

cp For_cv_select/* ${tmpDir}/cv_pool
cp gcta64_2021Jun02 ${tmpDir}/cv_pool

gcta64=${tmpDir}/cv_pool/gcta64_2021Jun02


for sd in 0 ;do
for ARR in STA ; do
for region in fe_en diff_en whole ;do
for c in 2 10 20 ; do

geno_mac2_F=Chrs_22/OKG.EUR.chr${c}.mac2.507

comm_SV_ID=${tmpDir}/cv_pool/chr${c}.comm.SV_${region}
comm_SGV_ID=${tmpDir}/cv_pool/chr${c}.comm.SGV_${region}
rare_SV_ID=${tmpDir}/cv_pool/chr${c}.rare.SV_${region}
rare_SGV_ID=${tmpDir}/cv_pool/chr${c}.rare.SGV_${region}


for ldp in Split; do # there is no difference for 1-1
for ln_mn in 1-1 ;do

ln=$(echo $ln_mn|awk -F'-' '{print $1}')
mn=$(echo $ln_mn|awk -F'-' '{print $2}')

spl=l${ln}m${mn}
coN=$((ln*mn*2))


mgrm=/storage/yangjianLab/baiweiyang/SV_Imputation_Project_s5_2t/17.SIMULATION/FOR_heritability/Imputed_EUR/STA/zFor_LDMS/z${ldp}_${spl}_rmM${sd}sd/chr${c}.grm.files

# cv_N is mean cv_N of SVs
for q2 in 0.2 0.3 0.4; do 
	for cv_N in 30 60 90 ; do 
		for rc in m2080 mC m1090 m3070 mR ; do 
			if [ $rc == 'm2080' ]; then
				if [ $cv_N == '30' ]; then SV_comm=24; SV_rare=6; fi
				if [ $cv_N == '60' ]; then SV_comm=48; SV_rare=12; fi
				if [ $cv_N == '90' ]; then SV_comm=72; SV_rare=18; fi
			fi
			if [ $rc == 'm1090' ]; then
				if [ $cv_N == '30' ]; then SV_comm=27; SV_rare=3; fi
				if [ $cv_N == '60' ]; then SV_comm=54; SV_rare=6; fi
				if [ $cv_N == '90' ]; then SV_comm=81; SV_rare=9; fi
			fi
			if [ $rc == 'm3070' ]; then
				if [ $cv_N == '30' ]; then SV_comm=21; SV_rare=9; fi
				if [ $cv_N == '60' ]; then SV_comm=42; SV_rare=18; fi
				if [ $cv_N == '90' ]; then SV_comm=63; SV_rare=27; fi
			fi
			if [ $rc == 'mC' ]; then
				if [ $cv_N == '30' ]; then SV_comm=30; SV_rare=0; fi
				if [ $cv_N == '60' ]; then SV_comm=60; SV_rare=0; fi
				if [ $cv_N == '90' ]; then SV_comm=90; SV_rare=0; fi
			fi
			if [ $rc == 'mR' ]; then
				if [ $cv_N == '30' ]; then SV_comm=0; SV_rare=30; fi
				if [ $cv_N == '60' ]; then SV_comm=0; SV_rare=60; fi
				if [ $cv_N == '90' ]; then SV_comm=0; SV_rare=90; fi
			fi

			SGV_comm_1=$((${SV_comm}*1)); SGV_rare_1=$((${SV_rare}*1))
			SGV_comm_3=$((${SV_comm}*3)); SGV_rare_3=$((${SV_rare}*3))
			SGV_comm_6=$((${SV_comm}*6)); SGV_rare_6=$((${SV_rare}*6))
			SGV_comm_9=$((${SV_comm}*9)); SGV_rare_9=$((${SV_rare}*9))

			for cvt in SV SGV mix11 mix13 mix16 mix19;do
				pheno_DD=${tmpDir}/pheno_rmM${sd}sd_${region}_chr${c}_${ldp}_${spl}_qsq${q2}_N${cv_N}_${rc}_${cvt}; rm -rf ${pheno_DD}; mkdir -p ${pheno_DD}
				greml_DD=${tmpDir}/greml_rmM${sd}sd_${region}_chr${c}_${ldp}_${spl}_qsq${q2}_N${cv_N}_${rc}_${cvt}; rm -rf ${greml_DD}; mkdir -p ${greml_DD}

				for rep in $(seq ${s} ${e}); do

					if [ ${cvt} == 'SV' ]; then 
						cat ${comm_SV_ID} | shuf -n ${SV_comm} > ${pheno_DD}/rep_${rep}.cv
						cat ${rare_SV_ID} | shuf -n ${SV_rare} >> ${pheno_DD}/rep_${rep}.cv
					fi

					if [ ${cvt} == 'SGV' ]; then # use the same number as SVs
						cat ${comm_SGV_ID} |shuf -n ${SGV_comm_1} > ${pheno_DD}/rep_${rep}.cv
						cat ${rare_SGV_ID} |shuf -n ${SGV_rare_1} >> ${pheno_DD}/rep_${rep}.cv
					fi

					if [ ${cvt} == 'mix11' ]; then
						cat ${comm_SV_ID} | shuf -n ${SV_comm} > ${pheno_DD}/rep_${rep}.cv
						cat ${rare_SV_ID} | shuf -n ${SV_rare} >> ${pheno_DD}/rep_${rep}.cv
						cat ${comm_SGV_ID} |shuf -n ${SGV_comm_1} >> ${pheno_DD}/rep_${rep}.cv
						cat ${rare_SGV_ID} |shuf -n ${SGV_rare_1} >> ${pheno_DD}/rep_${rep}.cv
					fi

					if [ ${cvt} == 'mix13' ]; then
						cat ${comm_SV_ID} | shuf -n ${SV_comm} > ${pheno_DD}/rep_${rep}.cv
						cat ${rare_SV_ID} | shuf -n ${SV_rare} >> ${pheno_DD}/rep_${rep}.cv
						cat ${comm_SGV_ID} |shuf -n ${SGV_comm_3} >> ${pheno_DD}/rep_${rep}.cv
						cat ${rare_SGV_ID} |shuf -n ${SGV_rare_3} >> ${pheno_DD}/rep_${rep}.cv
					fi

					if [ ${cvt} == 'mix16' ]; then
						cat ${comm_SV_ID} | shuf -n ${SV_comm} > ${pheno_DD}/rep_${rep}.cv
						cat ${rare_SV_ID} | shuf -n ${SV_rare} >> ${pheno_DD}/rep_${rep}.cv
						cat ${comm_SGV_ID} |shuf -n ${SGV_comm_6} >> ${pheno_DD}/rep_${rep}.cv
						cat ${rare_SGV_ID} |shuf -n ${SGV_rare_6} >> ${pheno_DD}/rep_${rep}.cv
					fi

					if [ ${cvt} == 'mix19' ]; then
						cat ${comm_SV_ID} | shuf -n ${SV_comm} > ${pheno_DD}/rep_${rep}.cv
						cat ${rare_SV_ID} | shuf -n ${SV_rare} >> ${pheno_DD}/rep_${rep}.cv
						cat ${comm_SGV_ID} |shuf -n ${SGV_comm_9} >> ${pheno_DD}/rep_${rep}.cv
						cat ${rare_SGV_ID} |shuf -n ${SGV_rare_9} >> ${pheno_DD}/rep_${rep}.cv
					fi

					# generate phenotypes
					${gcta64} --simu-qt --simu-hsq ${q2} --threads ${tr} --bfile ${geno_mac2_F} --simu-causal-loci ${pheno_DD}/rep_${rep}.cv --out ${pheno_DD}/rep_${rep}

					# run greml with no constrain
					pheno=${pheno_DD}/rep_${rep}.phen
					# greml joint model
					${gcta64} --reml --reml-no-lrt --reml-no-constrain --mgrm  ${mgrm} --threads ${tr} --pheno ${pheno} --out ${greml_DD}/rep_${rep}.joint

					# summ result
					# joint
					hsq_F=${greml_DD}/rep_${rep}.joint.hsq
					log_F=${greml_DD}/rep_${rep}.joint.log
					elem="${spe},sd${sd},rep${rep},${region},chr${c},all,Joint,qsq${q2},N${cv_N},${rc},${cvt}"
					out_F=${aimD}/${spe}_M${sd}sd_${region}_chr${c}_all_Joint.txt
					if [ ${coN} == '12' ]; then 
						if [ -f "$hsq_F" ]; then 
							cat ${hsq_F}|sed -n "16,27p"|awk '{print $2}'|tr '\n' '\t'|awk -v OFS=',' '{print "sa","'${elem}'","iterF",$1+$2+$3+$4+$5+$6,$7+$8+$9+$10+$11+$12,$1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12}' >> ${out_F}
							its=$(cat ${log_F}|grep -E "^5\s"|awk '{print $1}')
							if [ -z "$its" ]; then
							cat ${log_F}|grep -E "^2\s"|awk -v OFS=',' '{print "su","'${elem}'","its",$3+$4+$5+$6+$7+$8,$9+$10+$11+$12+$13+$14,$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14}' >> ${out_F}
							else
							cat ${log_F}|grep -E "^5\s"|awk -v OFS=',' '{print "su","'${elem}'","its",$3+$4+$5+$6+$7+$8,$9+$10+$11+$12+$13+$14,$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14}' >> ${out_F}
							fi
						else
							cat ${log_F}|grep -E "^5\s"|awk -v OFS=',' '{print "fa","'${elem}'","its",$3+$4+$5+$6+$7+$8,$9+$10+$11+$12+$13+$14,$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14}' >> ${out_F}
						fi
					fi
					if [ ${coN} == '8' ]; then 
						if [ -f "$hsq_F" ]; then 
							cat ${hsq_F}|sed -n "12,19p"|awk '{print $2}'|tr '\n' '\t'|awk -v OFS=',' '{print "sa","'${elem}'","iterF",$1+$2+$3+$4,$5+$6+$7+$8,$1+$2+$3+$4+$5+$6+$7+$8}' >> ${out_F}
							its=$(cat ${log_F}|grep -E "^5\s"|awk '{print $1}')
							if [ -z "$its" ]; then
							cat ${log_F}|grep -E "^2\s"|awk -v OFS=',' '{print "su","'${elem}'","its",$3+$4+$5+$6,$7+$8+$9+$10,$3+$4+$5+$6+$7+$8+$9+$10}' >> ${out_F}
							else
							cat ${log_F}|grep -E "^5\s"|awk -v OFS=',' '{print "su","'${elem}'","its",$3+$4+$5+$6,$7+$8+$9+$10,$3+$4+$5+$6+$7+$8+$9+$10}' >> ${out_F}
							fi
						else
							cat ${log_F}|grep -E "^5\s"|awk -v OFS=',' '{print "fa","'${elem}'","its",$3+$4+$5+$6,$7+$8+$9+$10,$3+$4+$5+$6+$7+$8+$9+$10}' >> ${out_F}
						fi
					fi
					if [ ${coN} == '6' ]; then 
						if [ -f "$hsq_F" ]; then 
							cat ${hsq_F}|sed -n "10,15p"|awk '{print $2}'|tr '\n' '\t'|awk -v OFS=',' '{print "sa","'${elem}'","iterF",$1+$2+$3,$4+$5+$6,$1+$2+$3+$4+$5+$6}' >> ${out_F}
							its=$(cat ${log_F}|grep -E "^5\s"|awk '{print $1}')
							if [ -z "$its" ]; then
							cat ${log_F}|grep -E "^2\s"|awk -v OFS=',' '{print "su","'${elem}'","its",$3+$4+$5,$6+$7+$8,$3+$4+$5+$6+$7+$8}' >> ${out_F}
							else
							cat ${log_F}|grep -E "^5\s"|awk -v OFS=',' '{print "su","'${elem}'","its",$3+$4+$5,$6+$7+$8,$3+$4+$5+$6+$7+$8}' >> ${out_F}
							fi
						else
							cat ${log_F}|grep -E "^5\s"|awk -v OFS=',' '{print "fa","'${elem}'","its",$3+$4+$5,$6+$7+$8,$3+$4+$5+$6+$7+$8}' >> ${out_F}
						fi
					fi
					if [ ${coN} == '4' ]; then 
						if [ -f "$hsq_F" ]; then 
							cat ${hsq_F}|sed -n "8,11p"|awk '{print $2}'|tr '\n' '\t'|awk -v OFS=',' '{print "sa","'${elem}'","iterF",$1+$2,$3+$4,$1+$2+$3+$4}' >> ${out_F}
							its=$(cat ${log_F}|grep -E "^5\s"|awk '{print $1}')
							if [ -z "$its" ]; then
							cat ${log_F}|grep -E "^2\s"|awk -v OFS=',' '{print "su","'${elem}'","its",$3+$4,$5+$6,$3+$4+$5+$6}' >> ${out_F}
							else
							cat ${log_F}|grep -E "^5\s"|awk -v OFS=',' '{print "su","'${elem}'","its",$3+$4,$5+$6,$3+$4+$5+$6}' >> ${out_F}
							fi
						else
							cat ${log_F}|grep -E "^5\s"|awk -v OFS=',' '{print "fa","'${elem}'","its",$3+$4,$5+$6,$3+$4+$5+$6}' >> ${out_F}
						fi
					fi
					if [ ${coN} == '2' ]; then 
						if [ -f "$hsq_F" ]; then 
							cat ${hsq_F}|sed -n "6,7p"|awk '{print $2}'|tr '\n' '\t'|awk -v OFS=',' '{print "sa","'${elem}'","iterF",$1,$2,$1+$2}' >> ${out_F}
							its=$(cat ${log_F}|grep -E "^5\s"|awk '{print $1}')
							if [ -z "$its" ]; then
							cat ${log_F}|grep -E "^2\s"|awk -v OFS=',' '{print "su","'${elem}'","its",$3,$4,$3+$4}' >> ${out_F}
							else
							cat ${log_F}|grep -E "^5\s"|awk -v OFS=',' '{print "su","'${elem}'","its",$3,$4,$3+$4}' >> ${out_F}
							fi
						else
							cat ${log_F}|grep -E "^5\s"|awk -v OFS=',' '{print "fa","'${elem}'","its",$3,$4,$3+$4}' >> ${out_F}
						fi
					fi
				done
			done
		done
	done
done

done
done




done
done
done

done


rm -rf ${tmpDir}
rm -rf ${tmpDir}


done



for spe in main_mg;do  # total scenarios: 2430 * 2

aimD=zzSim_spe_${spe}
mkdir -p ${aimD}

tr=1
n=${SLURM_ARRAY_TASK_ID}

rep_N=7500

s=$((${n}*5-5+1))
e=$((${n}*5))


tmpDir=/data/bwy_${spe}/${n}; rm -rf ${tmpDir}; mkdir -p ${tmpDir}

mkdir -p ${tmpDir}/cv_pool
mkdir -p ${tmpDir}/output

cp For_cv_select/* ${tmpDir}/cv_pool
cp gcta64_2021Jun02 ${tmpDir}/cv_pool

gcta64=${tmpDir}/cv_pool/gcta64_2021Jun02


for sd in 0 ;do
for ARR in STA ; do
for region in fe_en diff_en whole ;do
for c in 2 10 20 ; do

geno_mac2_F=Chrs_22/OKG.EUR.chr${c}.mac2.507

comm_SV_ID=${tmpDir}/cv_pool/chr${c}.comm.SV_${region}
comm_SGV_ID=${tmpDir}/cv_pool/chr${c}.comm.SGV_${region}
rare_SV_ID=${tmpDir}/cv_pool/chr${c}.rare.SV_${region}
rare_SGV_ID=${tmpDir}/cv_pool/chr${c}.rare.SGV_${region}


for ldp in Split; do # there is no difference for 1-1
for ln_mn in 1-1 ;do # now for marginal, no use mgrm but single grm file

ln=$(echo $ln_mn|awk -F'-' '{print $1}')
mn=$(echo $ln_mn|awk -F'-' '{print $2}')

spl=l${ln}m${mn}
coN=1 # now for marginal, no use mgrm but single grm file


SGV_grm=/storage/yangjianLab/baiweiyang/SV_Imputation_Project_s5_2t/17.SIMULATION/FOR_heritability/Imputed_EUR/STA/zFor_LDMS/z${ldp}_${spl}_rmM${sd}sd/chr${c}.SGV.id.lds_1_maf_1
SV_grm=/storage/yangjianLab/baiweiyang/SV_Imputation_Project_s5_2t/17.SIMULATION/FOR_heritability/Imputed_EUR/STA/zFor_LDMS/z${ldp}_${spl}_rmM${sd}sd/chr${c}.SV.id.lds_1_maf_1

# cv_N is mean cv_N of SVs
for q2 in 0.2 0.3 0.4; do 
	for cv_N in 30 60 90 ; do 
		for rc in m2080 mC m1090 m3070 mR ; do 
			if [ $rc == 'm2080' ]; then
				if [ $cv_N == '30' ]; then SV_comm=24; SV_rare=6; fi
				if [ $cv_N == '60' ]; then SV_comm=48; SV_rare=12; fi
				if [ $cv_N == '90' ]; then SV_comm=72; SV_rare=18; fi
			fi
			if [ $rc == 'm1090' ]; then
				if [ $cv_N == '30' ]; then SV_comm=27; SV_rare=3; fi
				if [ $cv_N == '60' ]; then SV_comm=54; SV_rare=6; fi
				if [ $cv_N == '90' ]; then SV_comm=81; SV_rare=9; fi
			fi
			if [ $rc == 'm3070' ]; then
				if [ $cv_N == '30' ]; then SV_comm=21; SV_rare=9; fi
				if [ $cv_N == '60' ]; then SV_comm=42; SV_rare=18; fi
				if [ $cv_N == '90' ]; then SV_comm=63; SV_rare=27; fi
			fi
			if [ $rc == 'mC' ]; then
				if [ $cv_N == '30' ]; then SV_comm=30; SV_rare=0; fi
				if [ $cv_N == '60' ]; then SV_comm=60; SV_rare=0; fi
				if [ $cv_N == '90' ]; then SV_comm=90; SV_rare=0; fi
			fi
			if [ $rc == 'mR' ]; then
				if [ $cv_N == '30' ]; then SV_comm=0; SV_rare=30; fi
				if [ $cv_N == '60' ]; then SV_comm=0; SV_rare=60; fi
				if [ $cv_N == '90' ]; then SV_comm=0; SV_rare=90; fi
			fi

			SGV_comm_1=$((${SV_comm}*1)); SGV_rare_1=$((${SV_rare}*1))
			SGV_comm_3=$((${SV_comm}*3)); SGV_rare_3=$((${SV_rare}*3))
			SGV_comm_6=$((${SV_comm}*6)); SGV_rare_6=$((${SV_rare}*6))
			SGV_comm_9=$((${SV_comm}*9)); SGV_rare_9=$((${SV_rare}*9))

			for cvt in SV SGV mix11 mix13 mix16 mix19;do
				pheno_DD=${tmpDir}/pheno_rmM${sd}sd_${region}_chr${c}_${ldp}_${spl}_qsq${q2}_N${cv_N}_${rc}_${cvt}; rm -rf ${pheno_DD}; mkdir -p ${pheno_DD}
				greml_DD=${tmpDir}/greml_rmM${sd}sd_${region}_chr${c}_${ldp}_${spl}_qsq${q2}_N${cv_N}_${rc}_${cvt}; rm -rf ${greml_DD}; mkdir -p ${greml_DD}

				for rep in $(seq ${s} ${e}); do

					if [ ${cvt} == 'SV' ]; then 
						cat ${comm_SV_ID} | shuf -n ${SV_comm} > ${pheno_DD}/rep_${rep}.cv
						cat ${rare_SV_ID} | shuf -n ${SV_rare} >> ${pheno_DD}/rep_${rep}.cv
					fi

					if [ ${cvt} == 'SGV' ]; then # use the same number as SVs
						cat ${comm_SGV_ID} |shuf -n ${SGV_comm_1} > ${pheno_DD}/rep_${rep}.cv
						cat ${rare_SGV_ID} |shuf -n ${SGV_rare_1} >> ${pheno_DD}/rep_${rep}.cv
					fi

					if [ ${cvt} == 'mix11' ]; then
						cat ${comm_SV_ID} | shuf -n ${SV_comm} > ${pheno_DD}/rep_${rep}.cv
						cat ${rare_SV_ID} | shuf -n ${SV_rare} >> ${pheno_DD}/rep_${rep}.cv
						cat ${comm_SGV_ID} |shuf -n ${SGV_comm_1} >> ${pheno_DD}/rep_${rep}.cv
						cat ${rare_SGV_ID} |shuf -n ${SGV_rare_1} >> ${pheno_DD}/rep_${rep}.cv
					fi

					if [ ${cvt} == 'mix13' ]; then
						cat ${comm_SV_ID} | shuf -n ${SV_comm} > ${pheno_DD}/rep_${rep}.cv
						cat ${rare_SV_ID} | shuf -n ${SV_rare} >> ${pheno_DD}/rep_${rep}.cv
						cat ${comm_SGV_ID} |shuf -n ${SGV_comm_3} >> ${pheno_DD}/rep_${rep}.cv
						cat ${rare_SGV_ID} |shuf -n ${SGV_rare_3} >> ${pheno_DD}/rep_${rep}.cv
					fi

					if [ ${cvt} == 'mix16' ]; then
						cat ${comm_SV_ID} | shuf -n ${SV_comm} > ${pheno_DD}/rep_${rep}.cv
						cat ${rare_SV_ID} | shuf -n ${SV_rare} >> ${pheno_DD}/rep_${rep}.cv
						cat ${comm_SGV_ID} |shuf -n ${SGV_comm_6} >> ${pheno_DD}/rep_${rep}.cv
						cat ${rare_SGV_ID} |shuf -n ${SGV_rare_6} >> ${pheno_DD}/rep_${rep}.cv
					fi

					if [ ${cvt} == 'mix19' ]; then
						cat ${comm_SV_ID} | shuf -n ${SV_comm} > ${pheno_DD}/rep_${rep}.cv
						cat ${rare_SV_ID} | shuf -n ${SV_rare} >> ${pheno_DD}/rep_${rep}.cv
						cat ${comm_SGV_ID} |shuf -n ${SGV_comm_9} >> ${pheno_DD}/rep_${rep}.cv
						cat ${rare_SGV_ID} |shuf -n ${SGV_rare_9} >> ${pheno_DD}/rep_${rep}.cv
					fi

					# generate phenotypes
					${gcta64} --simu-qt --simu-hsq ${q2} --threads ${tr} --bfile ${geno_mac2_F} --simu-causal-loci ${pheno_DD}/rep_${rep}.cv --out ${pheno_DD}/rep_${rep}

					# run greml with no constrain
					pheno=${pheno_DD}/rep_${rep}.phen
					# greml joint model
					${gcta64} --reml --reml-no-lrt --reml-no-constrain --grm  ${SGV_grm} --threads ${tr} --pheno ${pheno} --out ${greml_DD}/rep_${rep}.SGV.marginal
					${gcta64} --reml --reml-no-lrt --reml-no-constrain --grm  ${SV_grm} --threads ${tr} --pheno ${pheno} --out ${greml_DD}/rep_${rep}.SV.marginal

					# summ result
					# marginal
					hsq_SGV_F=${greml_DD}/rep_${rep}.SGV.marginal.hsq
					log_SGV_F=${greml_DD}/rep_${rep}.SGV.marginal.log

					hsq_SV_F=${greml_DD}/rep_${rep}.SV.marginal.hsq
					log_SV_F=${greml_DD}/rep_${rep}.SV.marginal.log
					
					elem="${spe},sd${sd},rep${rep},${region},chr${c},all,Marginal,qsq${q2},N${cv_N},${rc},${cvt}"

					out_SGV_F=${aimD}/${spe}_M${sd}sd_${region}_chr${c}_all_Marginal_SGV.txt
					out_SV_F=${aimD}/${spe}_M${sd}sd_${region}_chr${c}_all_Marginal_SV.txt
					
					# coN always is 1
					if [ ${coN} == '1' ]; then 
						if [ -f "$hsq_SGV_F" ]; then 
							cat ${hsq_SGV_F}|sed -n "5p"|awk '{print $2}'|awk -v OFS=',' '{print "sa","'${elem}'","iterF",$1}' >> ${out_SGV_F}
							its=$(cat ${log_SGV_F}|grep -E "^5\s"|awk '{print $1}')
							if [ -z "$its" ]; then
							cat ${log_SGV_F}|grep -E "^2\s"|awk -v OFS=',' '{print "su","'${elem}'","its",$3}' >> ${out_SGV_F}
							else
							cat ${log_SGV_F}|grep -E "^5\s"|awk -v OFS=',' '{print "su","'${elem}'","its",$3}' >> ${out_SGV_F}
							fi
						else
							cat ${log_SGV_F}|grep -E "^5\s"|awk -v OFS=',' '{print "fa","'${elem}'","its",$3}' >> ${out_SGV_F}
						fi
						if [ -f "$hsq_SV_F" ]; then 
							cat ${hsq_SV_F}|sed -n "5p"|awk '{print $2}'|awk -v OFS=',' '{print "sa","'${elem}'","iterF",$1}' >> ${out_SV_F}
							its=$(cat ${log_SV_F}|grep -E "^5\s"|awk '{print $1}')
							if [ -z "$its" ]; then
							cat ${log_SV_F}|grep -E "^2\s"|awk -v OFS=',' '{print "su","'${elem}'","its",$3}' >> ${out_SV_F}
							else
							cat ${log_SV_F}|grep -E "^5\s"|awk -v OFS=',' '{print "su","'${elem}'","its",$3}' >> ${out_SV_F}
							fi
						else
							cat ${log_SV_F}|grep -E "^5\s"|awk -v OFS=',' '{print "fa","'${elem}'","its",$3}' >> ${out_SV_F}
						fi
					fi
				done
			done
		done
	done
done

done
done




done
done
done

done


rm -rf ${tmpDir}
rm -rf ${tmpDir}


done



tmpD=/data/bwy/main_jt_ext/${sn}

rm -rf ${tmpD}
mkdir -p ${tmpD}

meta_F=meta.270

region=$(cat ${meta_F}|sed -n "${sn}p"|awk '{print $1}')
c=$(cat ${meta_F}|sed -n "${sn}p"|awk '{print $2}')
rc=$(cat ${meta_F}|sed -n "${sn}p"|awk '{print $3}')
cvt=$(cat ${meta_F}|sed -n "${sn}p"|awk '{print $4}')

for spe in main_jt;do  # total scenarios: 2430
aimD=zzSim_spe_${spe}
for sd in 0 ;do
for ARR in STA ; do # only run STA
for ldp in Split; do # 1 bin (no bin) also used Split since it does not matter
for ln_mn in 1-1 ;do
ln=$(echo $ln_mn|awk -F'-' '{print $1}')
mn=$(echo $ln_mn|awk -F'-' '{print $2}')
spl=l${ln}m${mn}
for q2 in 0.2 0.3 0.4 ; do 
for cv_N in 30 60 90 ; do 

ori_out_F=${aimD}/${spe}_M${sd}sd_${region}_chr${c}_all_Joint.txt

sce_out_F=${tmpD}/Joint_${spe}_sd${sd}_${region}_chr${c}_all_Joint_qsq${q2}_N${cv_N}_${rc}_${cvt}.txt_rank_by_its
sce_out_F_sa=${tmpD}/Joint_${spe}_sd${sd}_${region}_chr${c}_all_Joint_qsq${q2}_N${cv_N}_${rc}_${cvt}.txt_sa
sce_out_F_getPos=${tmpD}/Joint_${spe}_sd${sd}_${region}_chr${c}_all_Joint_qsq${q2}_N${cv_N}_${rc}_${cvt}.txt_rank_by_its_getPos

cat ${ori_out_F}|grep -w ${spe}|grep -w sd${sd}|grep -w ${region}|grep -w chr${c}|grep -w qsq${q2}|grep -w N${cv_N}|grep -w ${rc}|grep -w ${cvt}|grep -w its|\
sort -t ',' -k15 -g|awk -v OFS=',' '{print $0,NR}' > ${sce_out_F}

cat ${ori_out_F}|grep -w ${spe}|grep -w sd${sd}|grep -w ${region}|grep -w chr${c}|grep -w qsq${q2}|grep -w N${cv_N}|grep -w ${rc}|grep -w ${cvt}|grep -w sa > ${sce_out_F_sa}


TL=$(cat ${sce_out_F}|tail -n1|awk -F',' '{print $17}')

for rank in $(seq 1 ${TL}); do

rep=$(cat ${sce_out_F}|sed -n "${rank}p"|awk -F',' '{print $4}')
FoS=$(cat ${sce_out_F}|sed -n "${rank}p"|awk -F',' '{print $1}')

if [ $FoS == 'su' ]; then
cat ${sce_out_F_sa}|grep -w ${rep}|awk -F',' '{print $0}' >> ${sce_out_F_getPos}
fi

if [ $FoS == 'fa' ]; then
echo -e "fa,main_jt,sd${sd},${rep},${region},chr${c},all,Joint,qsq${q2},N${cv_N},${rc},${cvt},iterR,,," >> ${sce_out_F_getPos}
fi

done

done
done
done
done
done
done
done

mv ${tmpD}/*txt_rank_by_its_getPos ${aimD}/zSce


rm -rf ${tmpD}
rm -rf ${tmpD}




tmpD=/data/bwy/main_mg_ext_tibu1/${sn}

rm -rf ${tmpD}
mkdir -p ${tmpD}

meta_F=meta.540

region=$(cat ${meta_F}|sed -n "${sn}p"|awk '{print $1}')
c=$(cat ${meta_F}|sed -n "${sn}p"|awk '{print $2}')
out=$(cat ${meta_F}|sed -n "${sn}p"|awk '{print $3}')
rc=$(cat ${meta_F}|sed -n "${sn}p"|awk '{print $4}')
cvt=$(cat ${meta_F}|sed -n "${sn}p"|awk '{print $5}')

for spe in main_mg;do  # total scenarios: 2430 * 2
aimD=zzSim_spe_${spe}
for sd in 0 ;do
for ARR in STA ; do # only run STA
for ldp in Split; do # 1 bin (no bin) also used Split since it does not matter
for ln_mn in 1-1 ;do
ln=$(echo $ln_mn|awk -F'-' '{print $1}')
mn=$(echo $ln_mn|awk -F'-' '{print $2}')
spl=l${ln}m${mn}
for q2 in 0.2 0.3 0.4 ; do 
for cv_N in 30 60 90 ; do 

ori_out_F=${aimD}/${spe}_M${sd}sd_${region}_chr${c}_all_${out}.txt

sce_out_F=${tmpD}/${out}_${spe}_sd${sd}_${region}_chr${c}_all_Marginal_qsq${q2}_N${cv_N}_${rc}_${cvt}.txt_rank_by_its
sce_out_F_sa=${tmpD}/${out}_${spe}_sd${sd}_${region}_chr${c}_all_Marginal_qsq${q2}_N${cv_N}_${rc}_${cvt}.txt_sa
sce_out_F_getPos=${tmpD}/${out}_${spe}_sd${sd}_${region}_chr${c}_all_Marginal_qsq${q2}_N${cv_N}_${rc}_${cvt}.txt_rank_by_its_getPos

cat ${ori_out_F}|grep -w ${spe}|grep -w sd${sd}|grep -w ${region}|grep -w chr${c}|grep -w qsq${q2}|grep -w N${cv_N}|grep -w ${rc}|grep -w ${cvt}|grep -w its|\
sort -t ',' -k14 -g|awk -v OFS=',' '{print $0,NR}' > ${sce_out_F}

cat ${ori_out_F}|grep -w ${spe}|grep -w sd${sd}|grep -w ${region}|grep -w chr${c}|grep -w qsq${q2}|grep -w N${cv_N}|grep -w ${rc}|grep -w ${cvt}|grep -w sa > ${sce_out_F_sa}


TL=$(cat ${sce_out_F}|tail -n1|awk -F',' '{print $15}')

for rank in $(seq 1 ${TL}); do

rep=$(cat ${sce_out_F}|sed -n "${rank}p"|awk -F',' '{print $4}')
FoS=$(cat ${sce_out_F}|sed -n "${rank}p"|awk -F',' '{print $1}')

if [ $FoS == 'su' ]; then
cat ${sce_out_F_sa}|grep -w ${rep}|awk -F',' '{print $0}' >> ${sce_out_F_getPos}
fi

if [ $FoS == 'fa' ]; then
echo -e "fa,main_mg,sd${sd},${rep},${region},chr${c},all,Marginal,qsq${q2},N${cv_N},${rc},${cvt},iterR," >> ${sce_out_F_getPos}
fi

done

done
done
done
done
done
done
done

mv ${tmpD}/*txt_rank_by_its_getPos ${aimD}/zSce


rm -rf ${tmpD}
rm -rf ${tmpD}
