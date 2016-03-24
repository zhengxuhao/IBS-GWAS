##Population stratification by principal component analysis in EIGENSOFT 6.0.1 package

Note: Run EIGENSOFT using **LD-pruned binary files**


###A. Convert Plink Bfiles to EIGENSOFT format using CONVERTF

```
convertf -p <(printf "genotypename: raw-GWA-data.bed
snpname: raw-GWA-data.bim
indivname: raw-GWA-data.fam
outputformat: EIGENSTRAT
genotypeoutname: raw-GWA-data.eigenstratgeno
snpoutname: raw-GWA-data.snp
indivoutname: raw-GWA-data.ind")
```

###B. Run SmartPCA to check population stratification by principal component analysis

```
smartpca.perl \
-i raw-GWA-data.eigenstratgeno \
-a raw-GWA-data.snp \
-b raw-GWA-data.ind \
-o raw-GWA-data_pop_strat.pca \
-p raw-GWA-data_pop_strat.plot \
-e raw-GWA-data_pop_strat.eval \
-l raw-GWA-data_pop_strat.log \
-m 0 \
-t 10 \
-k 10 \
-s 6
```
  


###C. (**Optional**) Plot top 2 PCs by *ploteig* function

A pca plot (in PDF file) will be generated in above smartpca process, running on the top 2 PCs.  If that fails to work, run plotpig function to  

```
plotpig \
-i raw-GWA-data_pop_strat.pca.evec \
-c 1:2 \
-p Case:Control \
-x \
-o raw-GWA-data_pop_strat.plot
```


###D. Check STATISTICAL SIGNFICANCE of each principal component by twstats

```
twstats \
-t twtable \
-i raw-GWA-data_pop_strat.eval \
-o GWA-data_pop_strat_twout
```

 Check P vaule of each PC in file "GWA-data_pop_strat_twout", corrected only PCs with P values <0.05. 

###E. smarteigenstrat.perl: run EIGENSTRAT stratification correction.

**Optional** Run evec2pca function to tranfer *.evec* file to *.pca* (if .pca file was not generated in smartpca process)

```
evec2pca.perl 100 raw-GWA-data_pop_strat.pca.evec raw-GWA-data.ind raw-GWA-data_pop_strat.pca
```

Then run smarteigenstrat program:

``` 
smarteigenstrat.perl \
-i raw-GWA-data.eigenstratgeno \
-a raw-GWA-data.snp \
-b raw-GWA-data.ind \
-q NO \
-p raw-GWA-data_pop_strat.pca.evec \
-k 10 \
-o raw-GWA-data_pop_strat_cor.chisq \
-l raw-GWA-data_pop_strat_cor.log
```



###F. gc.perl: apply Genomic Control to the association statistics computed by EIGENSTRAT

```
gc.perl raw-GWA-data_pop_strat_cor.chisq raw-GWA-data_pop_strat_cor_report
```

Check the lamda value before and after correction in output report file, then change QC number for correction to adjust if necessary.


###G. Remove outliers

only apply to large datasets, for small dataset we go step 4 to remove possible outliers
 
Re-run smartpca function
