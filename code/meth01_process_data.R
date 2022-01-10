#' Load required libraries (ignore warning messages)
library(stringi)
library(magrittr)
library(data.table)
library(svd)
library(ewastools)

# setwd('/cloud/project')

#' ## Importing the data
#' Metadata stored in `data/pheno.csv` is imported using `fread` from the
#' data.table package and stored in the object `pheno`
pheno = fread("data/pheno.csv")

# Take a look at the metadata
pheno

#' Column `gsm` is the sample ID (GEO acession) as well as the first part of the
#' file name. For each sample there are two files, `$$$_Red.idat.gz` and
#' `$$$_Grn.idat.gz`, containing fluorescence intensities in the red and green
#' color channel, respectively. Import the data using `read_idats`. Provide only
#' the `$$$` part of the file name. Store in object `meth`.
meth = read_idats("data/" %s+% pheno$gsm)


#' `meth` is a list. Take a look at the some of the list elements
names(meth)

#' Which platform was used.
meth$platform

#' `manifest` is a data.table with information about all the probes on the chip
meth$manifest[4001:4010]

#' For example, `probe_id` is a unique identifier for each probe. The combination
#' of `chr` and `mapinfo` gives the genomic coordinates of the targeted loci.
#' The column `probe_type` refers to the type of loci that are targeted by these
#' probes. Here, all the probes are CpG sites, i.e., "cg".
#'
#' **QUESTION:** What other type of loci besides CpG sites are there?
#'
#' The column `channel` tells us the color channel and Infinium probe design type.
#' "Grn" and "Red" imply Infinium Type I, "Both" implies Type II.
#'
#' `addressU` and `addressM` are the unique identifier of the beads (a bead is not
#' is not the same as a probe).
#'
#' **QUESTION:** Why is the column `addressM` missing for entries with
#' `channel=="Both"`?
#'
#' Some columns (`index`,`OOBi`,`Ui`,`Mi`) are used for the internal workings of
#' the `ewastools` package and not of interest for the user


#' Fluorescence intensities observed at the bead locations are stores in
#' `M` (methylated) and `U` (unmethylated) (`U`).
dim(meth$M)
#'The matrix has 485,577 rows, one for each probe, and 35 columns, one for each
#' sample

#' Let's look at these fluorescence intensities for three probes in three samples.
meth$U[201:203,1:3]
meth$M[201:203,1:3]

#' **QUESTION:** Comparing corresponding entries across `U` and `M` above, what can
#' you infer about these loci?
#'
#' Because of the random assembly of the chips, the copy number can vary for each
#' bead, and the values in `U` and `M` actually are averages. Copy numbers are
#' stored in `V` (corresponding to `U`) and `N` (corresponding to `M`). For some
#' beads the copy number can be zero.
dim(meth$N)
meth$N[201:203,1:3]

#' **QUESTION:** How many probes are missing in the first sample?
#'
#' **QUESTION:** For type II probes, the entries in `N` and `V` are always the same. Why?
#'
#' The `meth` list includes some for elements, among them matrices and a manifest
#' for control probes. The control probes don't target CpG sites, but are used to
#' monitor the various experimental steps or for preprocessing. We will make use
#' of them below for quality control.








# ------------------------------------------------------------------------------
#' ## Quality control

#' Create a flag for samples that we want to exclude from the final dataset
pheno[,exclude:=FALSE]

#' ### Control metrics
#' The first quality check evaluates 17 control metrics which are describe in the
#' BeadArray Controls Reporter Software Guide from Illumina.

#' `control_metrics()` calculates the control metrics from the control probes
#' `sample_failure()` evaluates control metrics against thresholds

meth %>% control_metrics %>% sample_failure -> pheno$failed
#' Here we are making use of the pipe operator (%>%) from the `magrittr` package.
#' The line above is equivalent to (i.e., could be alternatively written as)
# pheno$failed = sample_failure(control_metrics(meth))

#' **QUESTION:** How many samples failed this check?


# ------------------------
#' ### Detection p-values

#' `detectionP()` calculates detecton p-values and stores them as in a new list
#' element `detP`. Detection p-values should be based on raw data, i.e., before
#' preprocessing such as dye-bias correction. We will use `mask()` later to drop
#' those observations with p-values above a specified cut-off
meth = ewastools::detectionP(meth)

#' Comparison of detection p-values for probes targeting the Y chromosome across
#' males and females
#' Get the indices of the probes
chrY = meth$manifest[chr=='Y',index]
#' There are 416 such probes
length(chrY)
#' Retrieve the corresponding p-values
detP = meth$detP[chrY,]
#' Count those below 0.01 as detected
detP = colSums(detP<0.01,na.rm=TRUE)

boxplot(split(detP,pheno$sex),ylab="# of detected Y chromosome probes")
split(detP,pheno$sex) %>% sapply(median)

#' Almost all of the 416 chromosome probes are called detected in the males whereas
#' among females 60 probes are detected on average.
#'
#' **QUESTION:** Excluding the Y chromosome and missing probes, what is the percentage
#' of undetected (0.01 cut-off) probes in the first sample. (Use the `table`) function.



# ------------------------
#' ### Dye-bias correction
#' Infinium BeadChips use two fluorescent dyes that are linked to the nucleotides
#' used in the the single-base extension step. A and T nucleotides use are linked
#' with a red dye (the red color channel), G and C nucleotides are linked with a
#' green dye (green color channel). Uncorrected data usually feature higher
#' intensities in the red color channel, the so-called dye bias. For probes of
#' Infinium type II design, which use separate color channels to measure the
#' methylated and unmethylated signal, this results in a shifted distribution of
#' betavalues. (Probes of Infinium design type I are not affected, as they measure
#' both signals in the same color channel.)
#'
#' `correct_dye_bias()` adjusts the green color channel using the red color channel
#' as reference.

meth %>% dont_normalize                      -> with_bias
meth %>% correct_dye_bias %>% dont_normalize -> corrected

#' Looking at the same loci as before we can check whether dye-bias correction
#' changed the beta values
#' Each probe is of a different type/color channel
meth$manifest$channel[201:203]

with_bias[201:203,1:3] %>% round(4)
corrected[201:203,1:3] %>% round(4)

#' **QUESTION:** Why are the beta values for cg06091566 unchanged?


#' Plotting beta values for heterozygous SNPs, we can observe the dye bias as a
#' deviation from 0.5. For the corrected data, the middle peak aligns with 0.5
snps = meth$manifest[probe_type=="rs" & channel=="Both"]$index

plot (density(with_bias[snps,14],na.rm=TRUE,bw=0.1),col=1)
lines(density(corrected[snps,14],na.rm=TRUE,bw=0.1),col=2)
abline(v=0.5,lty=3)
legend("topleft",col=1:2,legend=c("with_bias","corrected"),lwd=1)



meth %>% correct_dye_bias -> meth
rm(corrected,with_bias)



# ------------------------
#' ### Sex check
#' `check_sex()` computes the normalized average total fluorescence intensity of
#' probes targeting the X and Y chromosomes, and `predict_sex()` infers the sex of
#' the sample donors based on this information. This check should be applied after
#' dye-bias correction but before undetected probes are masked.
pheno[,c("X","Y"):=check_sex(meth)]
pheno[,predicted_sex := predict_sex(X,Y,which(sex=="m"),which(sex=="f"))]

#' Are there instances where the sex inferred from the methylation data does not
#' match the information in the recorded meta data?
pheno[sex!=predicted_sex,.(gsm,sex,predicted_sex)]

#' This is indeed the case for the sample with the ID GSM2260573. This most
#' likely means that this sample was mislabeled.

#' We flag this sample for exclusion
pheno[sex!=predicted_sex,exclude:=TRUE]

#' Let's also plot these data
plot(Y~X,data=pheno,type="n")
text(Y~X,labels=sex,col=ifelse(sex=="m",2,1),data=pheno)
#' The mislabeled sample can be easily spotted as it falls into the wrong cluster.

#' **QUESTION:** What is the sample ID of the one outlier in above plot.
#'

# ------------------------
#' Before the next QC steps, we should mask undetected probes. Here we use a
#' cut-off of 0.01
meth = mask(meth,0.01)

#' Now we can compute a "clean" matrix of beta values. Undetected probes are
#' set to NA
beta = dont_normalize(meth) 



# ------------------------
#' ### SNP outliers
#' `snp_outliers()` returns the average log odds of belonging to the outlier
#' component across all SNP probes. Here a cut-off of -4 is used
#' greater than -4 for exclusion.
snps = meth$manifest[probe_type=="rs"]$index
genotypes = call_genotypes(beta[snps,],learn=FALSE)
pheno$outlier = snp_outliers(genotypes)

stripchart(pheno$outlier,method="jitter",pch=4)
abline(v=-4,lty="dotted",col=2)

#' **QUESTION:** What is the ID of the sample failing this check? Compare the answer
#' to the answer of the previous question. What can we infer from this?
#'

#' Mark the failing sample for exclusion
pheno[outlier > -4,exclude:=TRUE]



# ------------------------
#' ### Check for duplicates
#' The function `enumerate_sample_donors()` can be used to find samples with the 
#' same genetic fingerprint
pheno$donor_id = enumerate_sample_donors(genotypes)

#' Find donor IDs represented more than once
pheno[,n:=.N,by=donor_id]
pheno[n>1,.(gsm,donor_id)]

pheno[gsm=="GSM2260543",exclude:=TRUE] # drop duplicate

#' **QUESTION:** How many different donors are there?
#'

# ------------------------
#' ### Principal component analysis
#' PCA is a popular feature reduction method: it projects high-dimensional data
#' onto a lower-dimensional representation while trying to retain as much
#' variability as possible.

set.seed(292846330)

#' We will apply PCA to the matrix of beta but without the probes targeting the
#' allosomes. Exclude SNP probes as well.
chrXY =  meth$manifest[ chr %in% c("X","Y") | probe_type == "rs"]$index
pcs = beta[-chrXY,]
pcs = pcs - rowMeans(pcs)
pcs = na.omit(pcs)
pcs = t(pcs)
pcs = trlan.svd(pcs,neig=2) # Just the first two principal components
pcs = pcs$u

pheno$pc1 = pcs[,1]
pheno$pc2 = pcs[,2]

plot(pc2~pc1,col=ifelse(sex=="m",2,1),data=pheno)
text(pc2~pc1,labels=pheno[34]$gsm,data=pheno[34],pos=2,offset=1,col=2)

#' GSM2219539 is actually a lung tissue sample from another GEO dataset. It
#' dominates the first principal component and should be excluded as it could
#' drastically change the results of the downstream analysis.

pheno[gsm=="GSM2219539",exclude:=TRUE]

#' PCA may be applied iteratively. After excluding samples that manifest as
#' outliers, repeating PCA can give very different principal components.
#'


# ------------------------
#' ### Leukocyte composition
#' This quality check will only apply in case of blood samples (blood is, however,
#' one of the most commonly studied tissues). The function `estimateLC()`
#' implements the Houseman method (doi.org/10.1186/1471-2105-13-86) to predict the
#' leukocyte composition. The user has the choice between various sets of model
#' parameters trained on various reference datasets (see `?estimateLC` for a list
#' of options). The function operates on the matrix of beta values.


#' we are using the Reinius reference dataset
LC = estimateLC(beta,ref="Reinius")

#' `LC` is a data.table. We can easily add to `pheno`
pheno = cbind(pheno,LC)

#' Plot the proportions of granulocytes versus (the column named `GR`)
plot(sort(pheno$GR),ylim=c(0,1))

#' We observe two outliers, a sample with a very low proportion of granulocytes,
#' and one with a proportion close to 100%

pheno[which.min(GR),.(gsm,exclude)]
#' This is the lung tissue sample from before

pheno[which.max(GR),.(gsm,exclude)]
#' This is actually a sample of purified granulocytes. Flag for exclusion.
pheno[gsm=="GSM1185585",exclude:=TRUE]



#' Let's confirm that the leukocyte composition is a confounded with smoking status
#' LC stratified by smoking status
LC$smoker = pheno$smoker
LC = melt(LC,value.name="proportion",variable.name="cell_type",id.vars="smoker")

boxplot(proportion ~ smoker+cell_type,LC,col=1:2,
	main="Cell type distribution by smoking status",xaxt="n")
axis(1,at=seq(from=1.5, to=11.5,by=2),adj=1,labels=unique(LC$cell_type))
legend("topright",c("Non-smoker","Smoker"),pch=15,bty='n',col=1:2)


# ------------------------
#' ### Clean up

#' Keep only samples that passed all quality checks
keep = which(!pheno$exclude)
pheno = pheno[keep]
beta  = beta[,keep]

#' Keep only the columns relevant for downstream analysis
pheno = pheno[,.(gsm,smoker,sex,CD4,CD8,NK,MO,GR,B)]

#' Also keep a copy of the manfist
manifest = copy(meth$manifest)

#' Write data to disk. (Do not acutally run the next line. Only when
#' running the script locally, you should un-comment the next line
#' in order to save the results to be re-used in the following scripts).
# save(pheno,manifest,beta,file="data/processed.rda")
