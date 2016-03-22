# IBS-GWAS

GWAS on multinational IBS case-control cohorts, about creating standardized protocol for quality control (QC), imputation and association analysis.


##QUALITY CONTROL (QC) PROTOCOL FOR IBS GWAS DATA

**Please note**: In this protocol, we are going to perform QC on a ‘per-individual’ basis before conducting QC on a ‘per-marker’ basis to maximize the number of markers remaining in the study (suggested by Nature Protocols. 2010, DOI: 10.1038/nprot.2010.116, PMID: 21085122)


####Data

In Plink format (PED, MAP or BED, BIM & FAM)

####Tool used in this protocol:

  * [PLINK 1.9](https://www.cog-genomics.org/plink2)
  
  * SMARTPCA function in [EIGENSOFT version 6.0.1](http://www.hsph.harvard.edu/alkes-price/software/) for running PCA 
  
  * [R](http://cran.r-project.org/): for statistical in data analysis and graphing 

###PART A: INDIVIDUAL SAMPLE QC


####1. Identification of individuals with discordant sex information


  In shell, type:

```
plink --bfile raw-GWA-data --check-sex --out raw-GWA-data 
grep PROBLEM raw-GWA-data.sexcheck > raw-GWA-data.sexprobs 
awk '{$1=$1 "\t";$2= $2 "\t"; print}' raw-GWA-data.sexprobs |cut -f1-2 > fail-sexcheck-qc.txt
```
 
 File “fail-sexcheck-qc.txt” contains family IDs and individual IDs of all these individuals to remove.    


####2. Identification of individuals with elevated missing data rates or outlying heterozygosity rate

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

####3. Identification of duplicated or related individuals

**A. Pruning data to reduce the number of markers to reduce computational complexity**

 * Using a window of 50 variants and a shift of 5 variants between windows, with an r2 cut-off of 0.2:  
 * Strong LD region will be excluded based on “high-LD-regions.txt”.
 
In shell, type:
```
 plink --file raw-GWA-data --exclude high-LD-regions.txt --range --indep-pairwise 50 5 0.2 --out raw-GWA-data
```

**B. Generating IBS matrix**

 In shell, type:
```
 plink --bfile raw-GWA-data --extract raw-GWA-data.prune.in --genome --out raw-GWA-data
```

 Then we will identify all pairs of individuals with IBS > 0.185 which outputs the ID of the individual from the pair with lowest call rate (the one with high call rate in a pair will be kept).  

 Note:  The expectation is that IBD = 1 for duplicates or monozygotic twins, IBD = 0.5 for first-degree relatives, IBD = 0.25 for second-degree relatives and IBD = 0.125 for third-degree relatives.  It is typical to remove one individual from each pair with an IBD value of > 0.1875, which is halfway between the expected IBD for third- and second-degree relatives.

 In shell, type:
```
 perl run-IBD-QC.pl raw-GWA-data

```
 A file **fail-IBD-QC.txt** will be generated to exclude these samples from downstream analyses.


####4. Population stratification by principal component analysis in EIGENSOFT 6.0.1 package


Note: Run EIGENSOFT using **LD-pruned binary files**


**A. Convert Plink Bfiles to EIGENSOFT format using CONVERTF**

```
convertf -p <(printf "genotypename: raw-GWA-data.bed
snpname: raw-GWA-data.bim
indivname: raw-GWA-data.fam
outputformat: EIGENSTRAT
genotypeoutname: raw-GWA-data_pop_strat.eigenstratgeno
snpoutname: raw-GWA-data_pop_strat.snp
indivoutname: raw-GWA-data_pop_strat.ind")
```

**B. Run SmartPCA to check population stratification by principal component analysis**

```
smartpca.perl \
-i raw-GWA-data_pop_strat.eigenstratgeno \
-a raw-GWA-data_pop_strat.snp \
-b raw-GWA-data_pop_strat.ind \
-o raw-GWA-data_pop_strat.pca \
-p raw-GWA-data_pop_strat.plot \
-e raw-GWA-data_pop_strat.eval \
-l raw-GWA-data_pop_strat.log \
-m 0 \
-t 100 \
-k 100 \
-s 6
```
**Interpertation**:
```
-i example.geno  : genotype file in any format
-a example.snp   : snp file in any format 
-b example.ind   : indiv file in any format 
-k k             : (Default is 10) number of principal components to output
-o example.pca   : output file of principal components.  Individuals removed
                   as outliers will have all values set to 0.0 in this file.
-p example.plot  : prefix of output plot files of top 2 principal components.
                   (labeling individuals according to labels in indiv file)
-e example.eval  : output file of all eigenvalues
-l example.log   : output logfile
-m maxiter       : (Default is 5) maximum number of outlier removal iterations.
                   To turn off outlier removal, set -m 0.
-t topk          : (Default is 10) number of principal components along which 
                   to remove outliers during each outlier removal iteration.
-s sigma         : (Default is 6.0) number of standard deviations which an
                   individual must exceed, along one of topk top principal
		   components, in order to be removed as an outlier.

OPTIONAL FLAGS:
-w poplist       : compute eigenvectors using populations in poplist only,
                   where poplist is an ASCII file with one population per line
-y plotlist      : output plot will include populations in plotlist only, 
                   where plotlist is an ASCII file with one population per line
-z badsnpname    : list of SNPs which should be excluded from the analysis
-q YES/NO        : If set to YES, assume that there is a single population and
                   the population field contains real-valued phenotypes.
                   (Corresponds to qtmode parameter in smartpca program.)
                   The default value for this parameter is NO.

```


**C. (Optional) Plot top 2 PCs by *ploteig* function**

A pca plot (in PDF file) will be generated in above smartpca process, running on the top 2 PCs.  If that fails to work, run plotpig function to  

```
plotpig \
-i raw-GWA-data_pop_strat.pca.evec \
-c 1:2 \
-p Case:Control \
-x \
-o raw-GWA-data_pop_strat.plot
```


**D. Check STATISTICAL SIGNFICANCE of each principal component by twstats**

```
> twstats \
-t twtable \
-i raw-GWA-data_pop_strat.eval \
-o GWA-data_pop_strat_twout
```

 Corrected only PCs with P values <0.05 

**E. smarteigenstrat.perl: run EIGENSTRAT stratification correction.**  

 **Optional** Run evec2pca function to tranfer *.evec* file to *.pca* (if .pca file was not generated in smartpca process)

 ```
 evec2pca.perl 100 raw-GWA-data_pop_strat.pca.evec raw-GWA-data_pop_strat.ind raw-GWA-data_pop_strat.pca
 ```

 Then run smarteigenstrat program:

``` 
smarteigenstrat.perl \
-i raw-GWA-data_pop_strat.eigenstratgeno \
-a raw-GWA-data_pop_strat.snp \
-b raw-GWA-data_pop_strat.ind \
-q NO \
-p raw-GWA-data_pop_strat.pca.evec \
-k 10 \
-o raw-GWA-data_pop_strat_cor.chisq \
-l raw-GWA-data_pop_strat_cor.log
```
**Interpretation**:

```
-i example.geno : genotype file in any format
-a example.snp : snp file in any format 
-b example.ind : individual file in any format
  We note that phenotype information will be contained in example.ind,
  either as Case/Control labels or quantitative phenotypes if -q set to YES.
-q YES/NO : If set to YES, use quantitative phenotypes in example.ind.
  If -q is set to YES, the third column of the input individual file
  in EIGENSTRAT format (or sixth column of input individual file in PED format)
  should be real numbers.  The value -100.0 signifies "missing data".
  If -q is set to NO, these values should be "Case" or "Control".
  The default value for the -q parameter is NO.     
-p example.pca : input file of principal components (output of smartpca.perl)
-k 1 : (Default is 10) number of principal components along which to
  correct for stratification.  Note that l must be less than or equal to
  the number of principal components reported in the file example.pca.
-o example.chisq : chisq association statistics.  File contains log of
  flags to eigenstrat program, followed by one line per SNP:
  The first entry of each line is Armitage chisq statistic (Armitage, 1955)
    defined as NSAMPLES x (correlation between genotype and phenotype)^2.
    If the set of individuals with genotype and phenotype both valid
    is monomorphic for either genotype or phenotype, then NA is reported.
  The second entry of each line is the EIGENSTRAT chisq statistic, defined as
    (NSAMPLES-l-1) x (corr between adjusted_genotype and adjusted_phenotype)^2.
    If the set of individuals with genotype and phenotype both valid
    is monomorphic for either genotype or phenotype, then NA is reported.
  Note: even if l=0, there is a tiny difference between the two statistics
    because Armitage uses NSAMPLES while we use NSAMPLES-1, which we
    consider to be appropriate.
-l example.log : standard output file
```

**F. gc.perl: apply Genomic Control to the association statistics computed by EIGENSTRAT**

```
 gc.perl raw-GWA-data_pop_strat_cor.chisq raw-GWA-data_pop_strat_cor_report
```

 Check the lamda value before and after correction, then change QC number for correction to adjust.

**G. Remove outliers**

 only apply to large datasets, for small dataset we go step 5 to remove possible outliers
 
 
####5. Plot individuals on components drawn from the HapMap reference populations to assess likely ancestry groupings.


####6. Removal of all individuals failing sample QC

In shell, type:
>cat fail-* | sort -k1 | uniq > fail-qc-inds.txt

The file **fail-qc-inds.txt** should now contain a list of unique individuals failing the previous QC steps. 

To remove them from the data set, type the following command at the shell prompt:

>plink --bfile raw-GWA-data --remove fail-qc-inds.txt --make-bed --out clean-inds-GWA-data


##PART B: SNP QC

1. Identification of all SNPs with an excessive missing data rate

plink --bfile clean-inds-GWA-data --missing --out clean-inds-GWA-data


Then to plott a histogram of the missing genotype rate 

In R environment,  run Rscript “lmiss-hist.R”.

This will create the graph “clean-inds-GWA-data.lmiss.pdf”.  Examine the plot to decide a reasonable threshold at which to exclude SNPs based on elevated missing data.  

In general, we set 0.05 for this value.
2. Test SNPs for different genotype call rates between cases and controls

In shell, type

plink --bfile clean-inds-GWA-data --test-missing --out clean-inds-GWA-data

To highlight all SNPs with significant differences in case and control call rates (p<10-5), we run the following script:  

In shell, type:

perl run-diffmiss-qc.pl clean-inds-GWA-data

The command creates a file called “fail-diffmiss-qc.txt”, which can be used to exclude these SNPs from downstream association analyses.

3.  MAF <0.01
4.  HWE<10-7
5. Removal of all SNPs failing QC

In shell, type:

plink --bfile clean-inds-GWA-data --exclude fail-diffmiss-qc.txt --geno 0.05 –maf 0.01 --hwe 0.0000001 --make-bed --out clean-GWA-data


Then we will get the clean data after QC for downstream imputation or association analyses: “clean-GWA-data.bed”, “clean-GWA-data.bim” and “clean-GWA-data.fam”. 
