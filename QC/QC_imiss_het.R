imiss=read.table("TEST.major.imiss",h=T)
imiss$logF_MISS = log10(imiss[,6])
het=read.table("TEST.major.het",h=T)
het$meanHet = (het$N.NM. - het$O.HOM.)/het$N.NM.
library("geneplotter")
colors  <- densCols(imiss$logF_MISS,het$meanHet)
pdf("TEST.major.imiss-vs-het.pdf")
plot(imiss$logF_MISS,het$meanHet, pch=20, xlab="Proportion of missing genotypes", ylab="Heterozygosity rate",axes=T)
axis(2,at=c(0,0.05,0.10,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5),tick=T)
axis(1,at=c(-3,-2,-1,0),labels=c(0.001,0.01,0.1,1))
abline(h=mean(het$meanHet)-(3*sd(het$meanHet)),col="RED",lty=2)
abline(h=mean(het$meanHet)+(3*sd(het$meanHet)),col="RED",lty=2)
abline(v=-1.522879, col="RED", lty=2)


subset1 <- subset(het,het$meanHet>mean(het$meanHet)+(3*sd(het$meanHet))
                  |het$meanHet<mean(het$meanHet)-(3*sd(het$meanHet)))
subset2 <- subset(imiss,imiss$logF_MISS>log10(0.03))
subset1 <- subset1[,1:2]
subset2 <- subset2[,1:2]
subset <- rbind(subset1, subset2)

subset <- subset[!duplicated(subset$IID),]
write.table(x = subset, file = "fail-imisshet-qc.txt", sep = "\t", row.names = F,col.names = F)
