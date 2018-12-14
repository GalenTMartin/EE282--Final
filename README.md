## Summary
The purpose of this project is to locate mCHH islands within several grass genomes. 

## Background
  In plants, cytosines can be methylated in three sequence contexts: CG, CHG, and CHH where H = any base other than G. One main purpose of methylation is to silence transposable elements (TEs), which are mostly deleterious. Plants methylate cytosines within TEs in all three of the sequence contexts, and this appears to keep them in a heterochromatic state where they aren’t expressed and can’t proliferate. Strangely, genes are also often highly methylated with the opposite outcome. Highly methylated genes are generally moderately expressed across a broad range of tissues. However, while methylated TEs have cytosines methylated in all three contexts, genes are generally only highly methylated in the CG context. 
  
  Looking across the genome, one can imagine long sections of TEs tightly wrapped in heterochromatin with CG, CHG, and CHH cytosines methylated, and genic sections in euchromatin with only CG cytosines methylated. An additional epigenetic feature in this genomic landscape was recently discovered in maize, where genes usually have small regions (~100 bp) dense with methylated CHH cytosines before the transcription start site and after the terminator. Presence of these regions near genes, called mCHH islands, has been associated with more highly expressed genes in maize (Li et al. 2015). mCHH islands have also been found in many other species as well, but it was unclear whether the association between island presence and expression was universal (Niederhuth et al. 2016).
  
  A post doc in my lab, Dr. Danelle Seymour, is currently writing a paper examining gene body methylation across 8 grass species. So, I’m using her bisulfite sequencing and RNA-seq data to investigate mCHH islands and try to associate them with expression. At the moment, I'm still trying to figure out the best way to define these islands (both in terms of what proportion of CHH sites need to be methylated and how to associate these regions with genes) but for this project I have used the definition which has been used previously in mCHH island literature. In this rigid (and arbitrary) definition, mCHH islands are 100 bp regions with >25% mCHH (Gent et al., 2013; Li et al. 2015; Niederhuth et al. 2016).
  
  The information that I've included is only for maize. After writing these scripts, I had a talk with Danelle and Brandon, and we decided to try a Hidden Markov Model program which Danelle's ph.d. lab had previously used to identify differentially methylated regions. I decided it wouldn't really be fruitful to expand these scripts to the other seven species, but it could easily be done by changing file names within the following scripts.

### Method
#### Data
The data provided by Danelle consists of two files for every species.

The first shows all cytosines with sufficient coverage: 
```
chr  pos  strand  mC  umC  context  sequence
1    75   +       0   4    CHH      CCT
1    76   +       0   4    CHH      CTT
1    80   +       0   7    CHH      CCA
1    81   +       2   5    CHH      CAT
1    88   +       0   20   CHH      CAC
1    90   +       0   21   CHH      CAT
1    94   +       1   22   CHH      CTT
1    103  +       18  5    CHG      CTG
1    113  +       0   24   CHH      CCA
...
```
and the second shows methylated cytosines:
```
chr  pos  strand  mC  umC  context  sequence  p.value               p.adjust
1    81   +       2   5    CHH      CAT       0.000946744365673008  0.0217129118518571
1    103  +       18  5    CHG      CTG       1.63463016667941e-31  6.6695237535131e-29
1    132  +       5   19   CHG      CCG       5.62267294804554e-06  0.000152457277098441
1    133  +       24  0    CG       CGT       3.72398857354409e-48  1.66547138217842e-44
1    149  +       22  1    CG       CGG       7.59874346557903e-43  1.61123143135357e-39
1    168  +       18  6    CG       CGT       3.39689855065004e-31  1.36911496888313e-28
1    184  +       10  8    CG       CGA       7.00587326203953e-16  4.43728835678191e-14
1    219  -       3   0    CG       CGA       1.17862606852008e-06  3.41036682420103e-05
1    222  +       2   4    CHG      CAG       0.00174650251359567   0.0397254180690472
...
```
### Step 1: separating CHH cytosines into chromosomes
First, I turned full C and mC files into files with Cs in CHH context and mCs in CHH context using grep:
```
cd /bio/galentm/CHHislands/data/processed/
mkdir chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10

cd ../raw/
grep CHH P25_mC_cov2_BY_05.txt > /bio/galentm/CHHislands/data/processed/P25_unmCHH.txt
grep CHH P25_mC_cov2_BY_05.txt > /bio/galentm/CHHislands/data/processed/P25_mCHH.txt
```
Next, I broke the files up into separate chromosome files **containing only a list of positions** using the following code:
```
cd /bio/galentm/CHHislands/data/processed/
awk '{if ($1 == "1") print $2;}' P25_unmCHH.txt > /bio/galentm/CHHislands/data/processed/chr1/P25_chr1_CHH.txt
awk '{if ($1 == "2") print $2;}' P25_unmCHH.txt > /bio/galentm/CHHislands/data/processed/chr2/P25_chr2_CHH.txt
awk '{if ($1 == "3") print $2;}' P25_unmCHH.txt > /bio/galentm/CHHislands/data/processed/chr3/P25_chr3_CHH.txt
awk '{if ($1 == "4") print $2;}' P25_unmCHH.txt > /bio/galentm/CHHislands/data/processed/chr4/P25_chr4_CHH.txt
awk '{if ($1 == "5") print $2;}' P25_unmCHH.txt > /bio/galentm/CHHislands/data/processed/chr5/P25_chr5_CHH.txt
awk '{if ($1 == "6") print $2;}' P25_unmCHH.txt > /bio/galentm/CHHislands/data/processed/chr6/P25_chr6_CHH.txt
awk '{if ($1 == "7") print $2;}' P25_unmCHH.txt > /bio/galentm/CHHislands/data/processed/chr7/P25_chr7_CHH.txt
awk '{if ($1 == "8") print $2;}' P25_unmCHH.txt > /bio/galentm/CHHislands/data/processed/chr8/P25_chr8_CHH.txt
awk '{if ($1 == "9") print $2;}' P25_unmCHH.txt > /bio/galentm/CHHislands/data/processed/chr9/P25_chr9_CHH.txt
awk '{if ($1 == "10") print $2;}' P25_unmCHH.txt > /bio/galentm/CHHislands/data/processed/chr10/P25_chr10_CHH.txt
```
The code shown above separates all cytosines into separate chromosome files within separate chromosome directories. The following chunk separates all methylated cytosines into separate chromosome files within separate chromosome directories:
```
cd /bio/galentm/CHHislands/data/processed/
awk '{if ($1 == "1") print $2;}' P25_mCHH.txt > /bio/galentm/CHHislands/data/processed/P25_chr1_mCHH.txt
awk '{if ($1 == "2") print $2;}' P25_mCHH.txt > /bio/galentm/CHHislands/data/processed/P25_chr2_mCHH.txt
awk '{if ($1 == "3") print $2;}' P25_mCHH.txt > /bio/galentm/CHHislands/data/processed/P25_chr3_mCHH.txt
awk '{if ($1 == "4") print $2;}' P25_mCHH.txt > /bio/galentm/CHHislands/data/processed/P25_chr4_mCHH.txt
awk '{if ($1 == "5") print $2;}' P25_mCHH.txt > /bio/galentm/CHHislands/data/processed/P25_chr5_mCHH.txt
awk '{if ($1 == "6") print $2;}' P25_mCHH.txt > /bio/galentm/CHHislands/data/processed/P25_chr6_mCHH.txt
awk '{if ($1 == "7") print $2;}' P25_mCHH.txt > /bio/galentm/CHHislands/data/processed/P25_chr7_mCHH.txt
awk '{if ($1 == "8") print $2;}' P25_mCHH.txt > /bio/galentm/CHHislands/data/processed/P25_chr8_mCHH.txt
awk '{if ($1 == "9") print $2;}' P25_mCHH.txt > /bio/galentm/CHHislands/data/processed/P25_chr9_mCHH.txt
awk '{if ($1 == "10") print $2;}' P25_mCHH.txt > /bio/galentm/CHHislands/data/processed/P25_chr10_mCHH.txt
```
From here on, I'll focus on chromosome 1 so that I don't have to write things 10 times.

### Step 2: breaking chromosomes into 100 bp tiles, and finding mCHH frequency within those tiles
For the next part, I passed the data into R and used the hist() function to find frequencies of non-menthylated and mCHHs within 100 bp tiles:
```
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

#sink() creates file from output
sink(file = "/dfs1/bio/galentm/CHHislands/data/processed/chr1/P25_chr1_CHH_binned.txt")
setwd("/dfs1/bio/galentm/CHHislands/data/processed/chr1/")

#max print has to be changed for this to work, since to file contains many rows. I just made it arbitrarily huge
options(max.print=9999999999)

#reads data (list of all CHH positions) into a df
CHH <- read.table("P25_chr1_CHH.txt", header = FALSE)
#extracts a vector from said df
CHHvec <- CHH[[1]]

#br is just a vector containing a sequence of numbers spaced 100 apart. 
#its components will serve as breaks within the hist() function below
#307039400 is just the position# of the last CHH in the chromosome (found using bash tail function) rounded to the next 100
br <- seq(0,307039400,by=100)

#The hist() function finds the number of CHHs within each of the 100 bp tiles
freq <- hist(CHHvec, breaks = br, include.lowest = TRUE, plot = FALSE)
imp <- cbind(freq[[1]],freq[[2]])
#colnames(imp) <- c("interval", "count")
imp

#I chose to save output as a df because I found that it makes import back into R easier later on. 
```
I used a similar R script to bin the mCHH data:
```
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
sink(file = "/dfs1/bio/galentm/CHHislands/data/processed/P25_chr1_mCHH_binned.txt")
setwd("/dfs1/bio/galentm/CHHislands/data/processed")
options(max.print=999999999)

CHH <- read.table("P25_chr1_mCHH.txt", header = FALSE)
CHHvec <- CHH[[1]]

br <- seq(0,307039400,by=100)
freq <- hist(CHHvec, breaks = br, include.lowest = TRUE, plot = FALSE)
imp <- cbind(freq[[1]],freq[[2]])
#colnames(imp) <- c("interval", "count")
imp
```
Next, I manipulated the data in the terminal to get the R outputs from above to be simple lists of frequencies:
```
awk ' { print $3 } ' P25_chr1_CHH_binned.txt > P25_chr1_justCHHfreqs.txt
awk ' { print $3 } ' P25_chr1_mCHH_binned.txt > P25_chr1_justmCHHfreqs.txt
```
To get percentages of mCHH in each 100 bp tile, I divided the mCHH frequencies by CHH frequencies in R:
```
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

#sink() creates file from output
sink(file = "/dfs1/bio/galentm/CHHislands/data/processed/chr1/P25_chr1_CHHpercents_R.txt")
setwd("/dfs1/bio/galentm/CHHislands/data/processed/chr1")
options(max.print=999999999)

CHH <- read.table("P25_chr1_justCHHfreqs.txt", header = FALSE)
mCHH <- read.table("P25_chr1_justmCHHfreqs.txt", header = FALSE)

#data as vectors
CHHvec <- CHH[[1]]
mCHHvec <- mCHH[[1]]

mCHHpercents <- mCHHvec/CHHvec
#mCHHpercents <- as.data.frame(mCHHpercents)
#mCHHpercents
br <- seq(0,307039400,by=100)+1

#creates final file with two columns, one with start position of tile/bin and one with mCHH percent
final <- cbind(br,mCHHpercents)
final
```
To find mCHH islands, I used the bash terminal:
```
awk ' { print $2, $3 } ' P25_chr1_CHHpercents_R.txt \
| grep ^[^m] \
| grep -v NaN \
| awk ' $2 > 0.25 { print } ' > mCHHislands.txt
```
My final output is found in the file "mCHHislands.txt," which contains two columns. One lists starting position of the 100 bp tile in question, and the other contains % mCHH

As mentioned previously, this turned out to not be the method by which we'll end up identifying these islands, so I never went through the process of actually associating them with genes. In principle, one could do this by using BEDtools and looking for overlap between the regions identified here and ~1 kb regions before and after genes in the annotation .gff files.

