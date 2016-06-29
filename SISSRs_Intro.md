SISSRs: Site Identification from Short Sequence Reads
======

Used for calling peaks in ChIPseq
-----

Pipeline by Lucas Nell, June 2016
--------

# Software
SISSRs and its manual can be downloaded from [here](http://sissrs.rajajothi.com/).
If you'd rather be fancy and use the Terminal to do all the downloading and extracting, 
here is how you'd do that:

```bash
# Make SISSRs folder and move working directory into it
mkdir ~/SISSRs
cd ~/SISSRs

# Download compressed tar archive from SISSRs website
wget http://dir.nhlbi.nih.gov/papers/lmi/epigenomes/sissrs/sissrs_v1.4.tar.gz

# Extract all files from it (there should only be 2)
tar -xvzf sissrs_v1.4.tar.gz

# Remove the archive
rm sissrs_v1.4.tar.gz
```


# Data

### Input files

SISSRs uses [BED files](https://genome.ucsc.edu/FAQ/FAQformat.html) as input. 
To make these files, I first aligned ChIPseq reads (in FASTQ format) to the mouse 
reference genome, resulting in BAM files. 
Some of these BAM files were then filtered, others not; the file name reveals what
filtering, if any, was done to a particular file (see below). 
These BAM files were then converted to BED files using `bedtools bamtobed`:

```bash
bedtools bamtobed -i FILE_NAME.bam > FILE_NAME.bed
```

BED files can be found in the directory
`/Volumes/MW_18TB/Lucas_Nell/lan/musChIP/bed`. I primarily used BED files that were
filtered for mapping quality (MAPQ) ≥ 20.

*`MAPQ = -10 * log10(p)`*

"... where *p* is an estimate of the probability that the alignment does not 
correspond to the read's true point of origin." 
([`bowtie2` manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml))

So MAPQ ≥ 20 should mean that there is a probability of ≤ 0.01 that that read's 
alignment is not in the correct location.



### Output files

The SISSRs algorithm produced files with the extension `.bsites`. These files are located
in the same directory as the BED files and have the same naming scheme.



### File naming scheme

The naming scheme for all files is as follows:

```
<strain(s)>_<chromosome>[<build>]_<antibody>[<input>][<filter>][_<P-value threshold>]
```
*Note*: Square brackets indicate name-parts not present in all files.

### Name descriptions:

- __`strain(s)`__
    + Mouse strain: B6, CAST, PWD, or WSB
    + If multiple strains are represented by a file, they are all put here, separated 
      by '-' (e.g., `CAST-PWD-WSB`)
- __`chromosome`__
    + Chromosome aligned to: 'X' or 'Y'
    + 'Y' is male-specific Y from Soh et al. (2014) in *Cell* 
      ([link](http://dx.doi.org/10.1016/j.cell.2014.09.052))
- __`build`__
    + Not present in files aligned to build mm10 (the newest mouse genome version)
    + 'b9' --> it was aligned to build mm9
- __`antibody`__
    + 'H' --> anti-H3K4me3
    + 'D' --> anti-DMC1
- __`input`__
    + 'i' --> Only present for input DNA samples
- __`filter`__
    + Not present in unfiltered files
    + 'x' --> 'XS:i'-tag filtering
    + 'q' --> MAPQ >= 20 filtering
- __`P-value threshold`__
    + Only present in *.bsites files (output from SISSRs algorithm)
    + 'Pv' indicates a p-value threshold other than default (default = 0.001, or 0.1%)
    + Number following 'Pv' indicates p-value threshold in percent
        Example: '_Pv10' for 10% threshold

__Example:__

File `B6_Xb9_Dx.bam` would be...

- strain B6
- aligned to chr X
- build mm9
- anti-DMC1 ChIP
- filtered for 'XS:i' tag.



# Running SISSRs

Assuming you've downloaded and extracted the SISSRs perl file (`sissrs.pl`) into a 
folder named  `SISSRs` in your home directory (`/Users/mwlab/SISSRs` on the lab 
desktop), you can directly call `sissrs.pl` by first running the following in your 
script:

```bash
export PATH=${PATH}:~/SISSRs
```

Now you can run the SISSRs algorithm in the directory for BED and bsites files.
The following code uses many of the parameters also used in Brick et al. (2012) 
([link](http://dx.doi.org/10.1038/nature11089)), plus a p-value cutoff of 0.01. 
It also assumes the above naming scheme to auto-generate the input DNA filename
and the output bsites filename. Lastly, I'm adding an object `output_dir` so that you 
don't actually make new files in the original directory.
This is because you are not yet worthy of trust.

Comments below explain what some lines of code are doing.

```bash
export PATH=${PATH}:~/SISSRs

cd /Volumes/MW_18TB/Lucas_Nell/lan/musChIP/bed

chip_bed=B6_Y_Dq.bed
p_val=0.01

# Using ticks to create a variable from the stdout from a command
# Piping to sed 's/string2/string2/' replaces any "string1" in the stdout to "string2"
input_dna_bed=`echo ${chip_bed} | sed 's/_D/_Di/g; s/_H/_Hi/g'`

# The commands ${var/%x/y} and ${var/#x/y} replaces "x" at the beginning and end of the
# variable var with "y"
# If still confused, try the following command:
# var="xmx"; echo ${var/%x/y}; echo ${var/#x/y}
output_name=${chip_bed/%.bed/}_Pv${p_val/#0./}.bsites

output_dir=~/Desktop

# Running SISSRs
sissrs.pl -i ${chip_bed} -o ${output_dir}/${output_name} \
    -s 2700000000 -a -F 1000 -E 10 -p ${p_val} -b ${input_dna_bed}
```


# SISSRs reference

Raja Jothi, Suresh Cuddapah, Artem Barski, Kairong Cui, Keji Zhao.
Genome-wide identification of in vivo protein-DNA binding sites from ChIP-Seq data.
*Nucleic Acids Research*, __36__(16):5221-31, 2008.
