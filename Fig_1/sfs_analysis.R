## this script is to compare the empirical SFS to the theoretical neutral expectation
## hebin
## 25 nov 2013

## read in data
x <- read.table(file=pipe("cut -f1-12 PTC_freq_from_sup_table_v6.csv"),head=TRUE)

## down-sample the dataset to 150 lines.
x$total <- x$ref.strain.num + x$alt.strain.num + x$other.strain.num
x1 <- subset(x, total >= 150, c("chr","pos","ref.strain.num","alt.strain.num","other.strain.num","total"))
x1$resample <- apply(x1[,3:5],1,function(z) rhyper(1,z[2],z[1]+z[3],150))
x2 <- subset(x1, resample>0)

## calculate the expected SFS -- the # of i/150 alleles is expected to be theta/i
## which means the proportion is (1/i) / sum(1/i with i from 1 to 149)
a <- sum(1/1:149)
exp.sfs <- nrow(x2) * 1/1:149 / a
obs.sfs <- tabulate(x2$resample,149)
res <- data.frame(Exp=exp.sfs,Obs=obs.sfs)
## you can change the way you want to group the variants below
res$group <- cut(1:149, br=c(0,1,15,150),labels=c("Singleton","<=10%",">10%"))

## output
cat("Observed Site Frequency Spectrum...\n")
tapply(res$Obs, res$group, sum) 
cat("Expected Site Frequency Spectrum...\n")
tapply(res$Exp, res$group, sum)
