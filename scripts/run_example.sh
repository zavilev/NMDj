#!/bin/bash
set -e
##################### config #########################################

#path to directory containing STAR
star_path='/opt/bio/bin'

#number of threads for parallel processing
threads=8 

#directory of pyIPSA workflow
pyipsa_dir='/gss/home/l.zavileisky.ext/Programs/pyIPSA'

#name of pyIPSA conda environment (with Snakemake and pyIPSA dependencies)
ipsa_env="ipsa"

#<path to conda>/etc/profile.d/conda.sh
condapath="/gss/home/l.zavileisky.ext/Programs/miniconda3/etc/profile.d/conda.sh"

################ end of config ######################################

#1 prepare reference
gtf='../example/input_data/chr1_filtered.gtf'
stringtie_gtf='../example/input_data/chr1_stringtie.gtf'
genome='../example/input_data/chr1.fa'
outgenome='../example/rsem_ref'
ref="${outgenome}/chr1"

zcat "${gtf}.gz" > "${gtf}"
zcat "${stringtie_gtf}.gz" > "${stringtie_gtf}"
zcat "${genome}.gz" > "${genome}"

mkdir -p ${outgenome}
rsem-prepare-reference -p ${threads} --star --star-path ${star_path} --gtf ${gtf} ${genome} ${ref}

echo "reference prepaired"

#2 simulate reads
model='../example/input_data/SRR1319322.model'
expr='../example/input_data/SRR1319322.chr1.isoforms.results'
theta=0.28
N=1000000
simdir='../example/simulation/'

mkdir -p ${simdir}
rsem-simulate-reads ${ref} \
                    ${model} \
                    ${expr} \
                    ${theta} ${N} "${simdir}SRR1319322"

echo "reads simulated"

#3 align simulated reads
fq1=${simdir}SRR1319322_1.fq
fq2=${simdir}SRR1319322_2.fq
stardir="../example/bam/"

mkdir -p ${stardir}
${star_path}/STAR --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 \
     --alignMatesGapMax 1000000 --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outFilterType BySJout \
     --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --sjdbScore 1 --outSAMstrandField intronMotif \
     --chimSegmentMin 15 --limitBAMsortRAM 38000000000 --runThreadN ${threads} --genomeDir ${outgenome} --readFilesIn ${fq1} ${fq2} \
     --outFileNamePrefix ${stardir} --outStd Log
        
echo "simulated reads aligned to genome"
        
#4 run pyIPSA
wd=$(pwd)
ipsaconfig="${wd}/../example/input_data/pyIPSA_config.yaml"
export ipsainp="${wd}/../example/bam"
export ipsaout="${wd}/../example/ipsa_output"
envsubst < ${ipsaconfig} > ${ipsaconfig}.tmp && mv ${ipsaconfig}.tmp ${ipsaconfig}

source ${condapath}
conda activate ${ipsa_env}
cd ${pyipsa_dir}
snakemake --configfile ${ipsaconfig} --cores ${threads}
cd ${wd}
conda deactivate

echo "pyIPSA finished"
echo "starting NMDj"

#5 run NMDj with simulated counts
printf '%s\n' ${ipsaout}/?6/* > ../example/input_data/ipsa_files.txt
./NMDj.py -g ${gtf} -n -r 'tag:MANE_Select' -p ../example/nmdj_output/ -q ../example/input_data/ipsa_files.txt --threads ${threads}

#6 run NMDj with novel transcripts from StringTie
./NMDj.py -g ${stringtie_gtf} -o -a ${gtf} -G ${genome} -p ../example/nmdj_output_stringtie/ --threads ${threads}

#7 compare NMDj PSI with ground truth simulated PSI
./plot_psi.py "${simdir}SRR1319322.sim.isoforms.results" "../example/nmdj_output" "../example/NMDj_PSI_vs_ground_truth.png"

echo $(date)
