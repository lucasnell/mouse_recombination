SISSRs can be downloaded from [here](http://sissrs.rajajothi.com/)

```bash
export SISSRs=${HOME}/SISSRs/sissrs.pl
export Name=13R_Y_Dq.bed
export NameI=Bi_Y_Dq.bed
export Pval=0.01


cd /lustre1/lan/musChIP/bed
# Including input file and p-value cutoff == 0.01

$SISSRs -i $Name -o ${Name/.bed/_Pv}${Pval/0./}.bsites \
-a -F 1000 -s 2700000000 -E 10 -p ${Pval} -b ${NameI}























# ========================================================================================
# ========================================================================================
# ========================================================================================

#               OLD STUFF (pre-December 2015)

# ========================================================================================
# ========================================================================================
# ========================================================================================


#PBS -S /bin/bash
#PBS -q batch
#PBS -N SISSRs
#PBS -l nodes=1:ppn=6:AMD
#PBS -l mem=12gb
#PBS -l walltime=10:00:00

export BED=/usr/local/apps/bedtools/latest/bin/bedtools
export SISSRs=${HOME}/SISSRs/sissrs.pl
export DIR=/lustre1/lan/Mus_ChIP/bam

export NAMES=(B6_anti-DMC1 cast_H3K4me3 pwd_H3K4me3 wsb_H3K4me3)

for name in ${NAMES[@]}
do
cd ${DIR}/${name}
$BED bamtobed -i ${name}_q20.bam >> ${name}_q20.bed
done


for name in ${NAMES[@]}
do
cd ${DIR}/${name}
mv ${name}_q20.bed ../q20_beds
done







export BED=/usr/local/apps/bedtools/latest/bin/bedtools
export SISSRs=${HOME}/SISSRs/sissrs.pl
export DIR=/lustre1/lan/Mus_ChIP/bam
export FILES=(B6_anti-DMC1 cast_H3K4me3 pwd_H3K4me3 wsb_H3K4me3)


cd ${DIR}/q20_beds


for f in ${FILES[@]}
do
$SISSRs -i ${f}_q20.bed -o ${f}_q20.bsites -a -F 1000 -s 2700000000 -E 10
done

-s 171031299
-s 89532639


# For on the desktop...
cd ~/Google_Drive/notes_scripts/R/data/bsites
for f in *.bsites
do
tail -n+58 $f | awk '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' > ${f/bsites/txt}
done

export NAMES=(cast_H9q.bsites pwd_H9q.bsites wsb_H9q.bsites)
for f in ${NAMES[@]}
do
tail -n+58 $f | awk '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' > ${f/bsites/txt}
done

# # THEN IN R...
# PeakDir = "~/Google_Drive/notes_scripts/R/data/bsites"
# peaks = c('/cast_H9q.txt', '/pwd_H9q.txt', '/wsb_H9q.txt')
# #####c('/wsb_Hq.txt', '/pwd_Hq.txt', '/cast_Hq.txt', '/B6_Hq.txt', '/B6_Dq.txt')
# AllPeaks = data.frame(matrix(nrow=0, ncol=6))
# colnames(AllPeaks) = c('Strain', 'cStart', 'cEnd', 'NumTags', 'Fold', 'Pval')
# for (p in peaks[1:3]){
#     temp_Peaks = read.table(paste0(PeakDir, p),
#         header = F, sep = '\t', comment.char = '')
#     colnames(temp_Peaks) = c('cStart', 'cEnd', 'NumTags', 'Fold', 'Pval')
#     temp_Peaks = temp_Peaks[! is.na(temp_Peaks$Pval), ]
#     temp_Peaks$Strain = unlist(strsplit(gsub('\\/', '', p), '_'))[1]
#     AllPeaks = rbind(AllPeaks, temp_Peaks)
# }
# write.table(AllPeaks, file = paste0(PeakDir, '/X9q.txt'), sep = '\t', row.names = F)
######  For B6 *.bsites
# for (p in peaks[4:5]){
#     temp_Peaks = read.table(paste0(PeakDir, p),
#         header = F, sep = '\t', comment.char = '')
#     # cStart    cEnd    NumTags    Fold    p-value
#     colnames(temp_Peaks) = c('cStart', 'cEnd', 'NumTags', 'Fold', 'Pval')
#     temp_Peaks = temp_Peaks[! is.na(temp_Peaks$Pval), ]
#     temp_Peaks$Strain = unlist(strsplit(gsub('\\/', '', p), '_'))[1]
#     write.table(temp_Peaks, file = paste0(PeakDir, p),
#         sep = '\t', row.names = F)
# }







# INPUT ACCESSIONS:
#	cast
# SRR1561130
#	pwd
# SRR1561134
#	wsb
# SRR1561132
#   B6 (PAIRED)
# SRR404052
# SRR404053
# SRR404054

export SRA=${HOME}/sratoolkit/bin

for s in SRR1561130 SRR1561134 SRR1561132
do $SRA/prefetch $s
done


for s in *.sra
do $SRA/fastq-dump --split-files --gzip $s
done

#--------------------
#     GETTING B6 INPUT DNA FASTQ
#--------------------

#PBS -S /bin/bash
#PBS -q batch
#PBS -N inputB6_fq
#PBS -l nodes=1:ppn=1:AMD
#PBS -l mem=1gb
#PBS -l walltime=10:00:00
export SRA=${HOME}/sratoolkit/bin
for s in SRR404052 SRR404053 SRR404054
do ${SRA}/fastq-dump --split-files --gzip $s
done





#--------------------
#     Aligned B6 BAM to BED
#       with 'XS:i' filtering
#--------------------

#PBS -S /bin/bash
#PBS -q batch
#PBS -N B6x_bed
#PBS -l nodes=1:ppn=1:AMD
#PBS -l mem=2gb
#PBS -l walltime=10:00:00

export BED=/usr/local/apps/bedtools/latest/bin/bedtools

cd /lustre1/lan/Mus_ChIP/bam

for b in B6*x.bam
do
    $BED bamtobed -i $b >> ${b/.bam/.bed}
    mv ${b/.bam/.bed} ./bed
done




#--------------------
#     SISSRs on all bed files
#	(done interactively)
#--------------------


export SISSRs=${HOME}/SISSRs/sissrs.pl

# These lists were created in ChIP_SISSRs.py
export NAMES=(B6_D.bed B6_H.bed cast_H.bed pwd_H.bed wsb_H.bed \
B6_Dq.bed B6_Hq.bed cast_Hq.bed pwd_Hq.bed wsb_Hq.bed B6_Dx.bed \
B6_Hx.bed cast_Hx.bed pwd_Hx.bed wsb_Hx.bed)

export NAMESI=(B6_DI.bed B6_HI.bed cast_HI.bed pwd_HI.bed wsb_HI.bed \
B6_DIq.bed B6_HIq.bed cast_HIq.bed pwd_HIq.bed wsb_HIq.bed B6_DIx.bed \
B6_HIx.bed cast_HIx.bed pwd_HIx.bed wsb_HIx.bed)

# Bash indexes start at zero
for i in {0..14}
do
$SISSRs -i ${NAMES[i]} -o ${NAMES[i]/.bed/.bsites} \
-a -F 1000 -s 2700000000 -E 10 -b ${NAMESI[i]}
done

#-s 171031299
#-s 89532639


# For chrX_mm9 files...
export SISSRs=${HOME}/SISSRs/sissrs.pl
export NAMES=(cast_H9q pwd_H9q wsb_H9q)
export NAMESI=(cast_H9Iq pwd_H9Iq wsb_H9Iq)

for i in {0..2}
do
$SISSRs -i ${NAMES[i]}.bed -o ${NAMES[i]}.bsites \
-a -F 1000 -s 2700000000 -E 10 -b ${NAMESI[i]}.bed
done





# Now less-strictly on B6 aligned to chrY and MAPQ >=20 filtered
# This is to identify "cold" spots
# Default P-value threshold = 0.001


export SISSRs=${HOME}/SISSRs/sissrs.pl
export Name=B_Y_Dq.bed
export NameI=Bi_Y_Dq.bed
# If changing below, make sure it has exactly two decimal places
export Pval=0.01


cd /lustre1/lan/musChIP/bed


$SISSRs -i $Name -o ${Name/.bed/_Pv}${Pval/0./}.bsites \
-a -F 1000 -s 2700000000 -E 10 -p ${Pval} -b ${NameI}

# cmd+k
cat ${Name/.bed/_Pv}${Pval/0./}.bsites | wc -l

cat ${Name/.bed/}.bsites | wc -l


# Just a test (removing -E and -F args)
# $SISSRs -i $Name -o ${Name/.bed/_Pv}${Pval/0./}_NoEF.bsites \
# -a -s 2700000000 -p ${Pval} -b ${NameI}
#
# cat ${Name/.bed/_Pv}${Pval/0./}_NoEF.bsites | wc -l
```
