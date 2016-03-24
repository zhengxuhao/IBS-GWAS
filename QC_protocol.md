#QUALITY CONTROL (QC) PROTOCOL FOR IBS GWAS DATA

##Reference 

1. Anderson, C. A. et al. Data quality control in genetic case-control association studies. Nat. Protoc. 5, 1564–1573 (2010).
2. Coleman, J. R. I. et al. Quality control, imputation and analysis of genome-wide genotyping data from the Illumina HumanCoreExome microarray. Brief. Funct. Genomics elv037 (2015). doi:10.1093/bfgp/elv037. [Script here](https://github.com/JoniColeman/gwas_scripts)

##Contact

**Tenghao Zheng**, tenghao.zheng@ki.se 

##Data Input

Plink binary format (BED, BIM & FAM)

##Tool used in this protocol:

Computer workstation with Unix or Linux operating system
 
  * Computer workstation with Unix or Linux operating system
  
  * [PLINK 1.9](https://www.cog-genomics.org/plink2)
  
  * SMARTPCA function in [EIGENSOFT version 6.0.1](http://www.hsph.harvard.edu/alkes-price/software/) for running PCA 
  
  * [R](http://cran.r-project.org/): for statistical in data analysis and graphing 

##PART A: INDIVIDUAL SAMPLE QC


###1. Identification of individuals with discordant sex information


  In shell, type:

```
plink --bfile raw-GWA-data --check-sex --out raw-GWA-data 
grep PROBLEM raw-GWA-data.sexcheck > raw-GWA-data.sexprobs 
awk '{$1=$1 "\t";$2= $2 "\t"; print}' raw-GWA-data.sexprobs |cut -f1-2 > fail-sexcheck-qc.txt
```
 
 File “fail-sexcheck-qc.txt” contains family IDs and individual IDs of all these individuals to remove.    


###2. Identification of individuals with elevated missing data rates or outlying heterozygosity rate

 We will remove

 * Individuals with elevated missing data rates (__missing>0,03__) 
 * Individuals with outlying heterozygosity rate (__out of mean +/- 3 SD__)

 In shell, type:
```
plink --bfile raw-GWA-data --missing --out raw-GWA-data 
plink --bfile raw-GWA-data --het --out raw-GWA-data 
```
 Then run rscript “QC_imiss_het.R” in R environment

 A graph **raw-GWA-data.imiss-vs-het.pdf** will be generated for checking missing calling rates and heterozygosity rate of all individuals.  
 
 A file **fail-imisshet-qc.txt** will be generated containing family IDs and individual IDs of all these individuals to remove.    

###3. Identification of duplicated or related individuals

**A. Pruning data to reduce the number of markers to reduce computational complexity**

 * Using a window of **50** variants and a shift of **5** variants between windows, with an r2 cut-off of **0.2**:  
 * Strong LD region will be excluded based on “high-LD-regions.txt”.
 
In shell, type:
```
plink \
--file raw-GWA-data \
--exclude high-LD-regions.txt \
--range \
--indep-pairwise 50 5 0.2 \
--out raw-GWA-data
```

**B. Generating IBS matrix**

 In shell, type:
```
plink \
--bfile raw-GWA-data \
--extract raw-GWA-data.prune.in \
--genome \
--out raw-GWA-data
```

 Then we will identify all pairs of individuals with **IBS > 0.185** which outputs the ID of the individual from the pair with lowest call rate (**the one with high call rate in a pair will be kept**).  

 Note:  The expectation is that IBD = 1 for duplicates or monozygotic twins, IBD = 0.5 for first-degree relatives, IBD = 0.25 for second-degree relatives and IBD = 0.125 for third-degree relatives.  It is typical to remove one individual from each pair with an IBD value of > 0.1875, which is halfway between the expected IBD for third- and second-degree relatives.

 In shell, type:
```
perl run-IBD-QC.pl raw-GWA-data

```
 A file **fail-IBD-QC.txt** will be generated to exclude these samples from downstream analyses.

 
###4. Plot individuals on components drawn from the HapMap reference populations to assess likely ancestry groupings.

Required: 
 
* HapMap reference population (hapmap3CEU.CHB.JPT.YRI.b37) has already been transformed to Hg19 build, including CEU CHB JPT YRI populations.  [script here](https://github.com/Wall-Facer/IBS-GWAS/blob/master/Build_HapMapb37_ref).
* hapmap3r2_CEU.CHB.JPT.YRI.no-at-cg-snps.txt (list of all HapMap3 markers without a/t & c/g markers)
* smarpca, convertf functions in EIGENSOFT 6.0.1 package
 
####A. excluding  A/T and C/G SNPs in HapMap3 data

```
plink \
--bfile raw-GWA-data \
--extract hapmap3r2_CEU.CHB.JPT.YRI.no-at-cg-snps.txt \
--make-bed \
--out raw-GWA-data.no-at-cg-snps
```
 
 
####B. Merge data with HapMap3 population

```
plink \
--bfile raw-GWA-data.no-at-cg-snps \
--bmerge hapmap3CEU.CHB.JPT.YRI.b37 \
--extract raw-GWA-data.prune.in \
--make-bed \
--out raw-GWA-data.Hapmap.merged
```

**If failed, run following (remove "_###_"):**

```
plink --bfile hapmap3CEU.CHB.JPT.YRI.b37 --flip raw-GWA-data.Hapmap.merged-merge.missnp --make-bed --out hapmap3CEU.CHB.JPT.YRI.b37_trial

plink --bfile raw-GWA-data.no-at-cg-snps --bmerge hapmap3CEU.CHB.JPT.YRI.b37_trial --make-bed --extract raw-GWA-data.pruned.in --out raw-GWA-data.Hapmap.merged
```



####C. Run Smarpca
 
```
convertf -p <(printf "genotypename: raw-GWA-data.Hapmap.merged.bed
snpname: raw-GWA-data.Hapmap.merged.bim
indivname: raw-GWA-data.Hapmap.merged.fam
outputformat: EIGENSTRAT
genotypeoutname: raw-GWA-data.Hapmap.merged.eigenstratgeno
snpoutname: raw-GWA-data.Hapmap.merged.snp
indivoutname: raw-GWA-data.Hapmap.merged.ind")

```

Followed by:

```
smartpca.perl \
-i raw-GWA-data.Hapmap.merged.eigenstratgeno \
-a raw-GWA-data.Hapmap.merged.snp \
-b raw-GWA-data.Hapmap.merged.ind \
-o raw-GWA-data.Hapmap.merged.pca \
-p raw-GWA-data.Hapmap.merged.plot \
-e raw-GWA-data.Hapmap.merged.eval \
-l raw-GWA-data.Hapmap.merged.log \
-m 0 \
-t 10 \
-k 10

```

Then check **raw-GWA-data.Hapmap.merged.plot.pdf** to decide whether to remove outliers or not, create **"fail-outliers-QC.txt"**.   
 
###5. Population stratification by principal component analysis

[Please find details here.](https://github.com/Wall-Facer/IBS-GWAS/blob/master/PCA_4_pop_stratification.md)

###6. Removal of all individuals failing sample QC

In shell, type:
```
cat fail-* | sort -k1 | uniq > fail-qc-inds.txt
```

The file **fail-qc-inds.txt** should now contain a list of unique individuals failing the previous QC steps. 

To remove them from the data set, type the following command at the shell prompt:

```
plink \
--bfile raw-GWA-data \
--remove fail-qc-inds.txt \
--make-bed \
--out GWA-data-ind-clean
```

##PART B: SNP QC

###1. Identification of all SNPs with an excessive missing data rate

```
plink \
--bfile clean-inds-GWA-data \
--missing \
--out clean-inds-GWA-data
```

Then to plott a histogram of the missing genotype rate 

In R environment,  run Rscript **“lmiss-hist.R”**.

This will create the graph **“clean-inds-GWA-data.lmiss.pdf”**.  

Examine the plot to decide a reasonable threshold at which to exclude SNPs based on elevated missing data.  

In general, we set **0.05** for this value.


###2. Test SNPs for different genotype call rates between cases and controls

In shell, type
```
plink \
--bfile clean-inds-GWA-data \
--test-missing \
--out clean-inds-GWA-data
```

To highlight all SNPs with significant differences in case and control call rates (**p<10-5**), we run the following script:  

```
perl run-diffmiss-qc.pl clean-inds-GWA-data
```

The command creates a file called **“fail-diffmiss-qc.txt”**, which can be used to exclude these SNPs from downstream association analyses.

###3.  MAF <0.01

###4.  HWE<10-7 (only within controls)

###5. Removal of all SNPs failing QC

In shell, type:
```
plink \
--bfile clean-inds-GWA-data \
--exclude fail-diffmiss-qc.txt \
--geno 0.05 \
–-maf 0.01 \
--hwe 0.0000001 \
--make-bed \
--out clean-GWA-data
```

Then we will get the clean data after QC for downstream imputation or association analyses: “clean-GWA-data.bed”, “clean-GWA-data.bim” and “clean-GWA-data.fam”. 
