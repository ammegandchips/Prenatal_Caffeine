###################################################
#      MATERNAL CAFFEINE ANALYSIS PLAN CODE       # 
###################################################
####################################################################################################################################################################
# The following R code will allow you to complete all the EWAS requested in the maternal caffeine analysis plan.
# If you have insufficient data to complete one or more of the EWAS, you can just skip those models.
# The code also produces .csv files summarising the variables included in the EWASs.
# You shouldn't have to rewrite or add to the following code, unless otherwise stated.
# There are just two inputs required for this analysis:
# 1) pheno: a dataframe containing all the "phenotype" data needed for this project. 
#    Each row is a sample (individual) and each column is a different variable. 
#    Necessary variable names are: "mat.caff","mat.coff","mat.tea","mat.cola","mat.bmi","sex","mat.ses","mat.age","mat.smoking","parity","gest.age","bcell","cd14","cd4t","cd8t","gran","nk"
#    If these columns are named differently in your dataset, please rename the columns accordingly
#    Details on how to code these variables are provided in the analysis plan.
# 2) meth: a matrix of methylation illumina beta values. Each column is a sample and each row is a probe on the array (450k or EPIC). 
#    Column names must correspond to the sample.id column in pheno.
####################################################################################################################################################################

###################################################
#              PACKAGES AND FUNCTIONS             #
###################################################

# Load required packages (if these are not already installed, you will have to install them as the first step)
library(sva)
library(tableone)
library(matrixStats)
library(limma)

# Setup the necessary functions
## Function to remove outliers using the IQR*3 (Tukey) method
IQR.removal <- function(meth.matrix){
  rowIQR <- rowIQRs(meth.matrix, na.rm = T)
  row2575 <- rowQuantiles(meth.matrix, probs = c(0.25, 0.75), na.rm = T)
  maskL <- meth.matrix < row2575[,1] - 3 * rowIQR 
  maskU <- meth.matrix > row2575[,2] + 3 * rowIQR 
  meth.matrix[maskL] <- NA
  meth.matrix[maskU] <- NA
  meth.matrix
}

## Function to generate surrogate variables and merge them with the phenotype data (used to adjust for batch)
SVA.generate <- function(meth.matrix, pheno.data, variable.of.interest, model.covariates,n.sv){
  intersecting.samples <- intersect(pheno.data$sample.id,colnames(meth.matrix))
  pheno.data <- na.omit(pheno.data[which(pheno.data$sample.id %in% intersecting.samples),unique(c("sample.id",variable.of.interest,model.covariates))])
  meth.matrix <- meth.matrix[,match(pheno.data$sample.id,colnames(meth.matrix))]
  k = which(is.na(meth.matrix), arr.ind=TRUE)
  meth.matrix[k] = rowMedians(meth.matrix, na.rm=TRUE)[k[,1]]
  mod = model.matrix(reformulate(paste0("pheno.data$",colnames(pheno.data[-1]))))
  mod0 = mod[,-grep(paste0(variable.of.interest,collapse="|"),colnames(mod))]
  sva.ret = sva(meth.matrix, mod=mod, mod0=mod0, n.sv=n.sv)
  SVs = as.data.frame(sva.ret$sv)
  colnames(SVs) <-paste0("sv",1:ncol(SVs))
  cbind(pheno.data,SVs)
}
## Function to run EWAS
ewas.function <-  function(meth.matrix, pheno.data, variable.of.interest){   
  meth.matrix <- meth.matrix[,match(pheno.data$sample.id,colnames(meth.matrix))]
  model.covariates <- colnames(pheno.data)[-which(colnames(pheno.data) %in% c(variable.of.interest,"sample.id"))]
  des = model.matrix(reformulate(paste0("pheno.data$",c(variable.of.interest,model.covariates))))
  fit = lmFit(meth.matrix, des)
  fit.ebayes = eBayes(fit)
  n = rowSums(!is.na(meth.matrix))
  se = (sqrt(fit.ebayes$s2.post) * fit.ebayes$stdev.unscaled[,grep(paste0(variable.of.interest,collapse="|"),colnames(fit.ebayes$stdev.unscaled))])
  res = data.frame(n=n,
                   coef=fit.ebayes$coefficient[,grep(paste0(variable.of.interest,collapse="|"),colnames(fit.ebayes$coefficient))],
                   se=se,
                   p=fit.ebayes$p.value[,grep(paste0(variable.of.interest,collapse="|"),colnames(fit.ebayes$p.value))])
  res
}

###################################################
#                  ANALYSIS                       #
###################################################

# Set initial parameters
study <- "ALSPAC" #change to your study identifier
timepoint <- "birth" #change depending on the age of the children with methylation samples. Can be "birth", "early_childhood", "late_childhood", "adolescence" or "adult"
cell.names <- tolower(names(cell.counts.cord))[-1]
traits.and.covariates <- c("mat.caff","mat.coff","mat.tea","mat.cola","mat.bmi","sex","mat.ses","mat.age","mat.smoking","parity","gest.age")
covariates <- c("mat.bmi","mat.ses","mat.age","mat.smoking","parity", cell.names)

# Check phenotype data
pheno <- complete.pheno.cord
for(i in 1:length(c("sample.id",traits.and.covariates,cell.names))) {
  print(ifelse(c("sample.id",traits.and.covariates,cell.names)[i] %in% colnames(pheno)==FALSE,
               paste("CAUTION: the variable called",c("sample.id",traits.and.covariates,cell.names)[i],"is missing from pheno"),
               paste("variable called",c("sample.id",traits.and.covariates,cell.names)[i],"is present in pheno")))
}

# Load methylation data (meth)
meth <- meth.matrix.cord

# IQR*3 method to remove outliers (if this has not already been applied to your data)
log.iqr <- data.frame(cpgs = row.names(meth),NAs.before.IQR3 = rowSums(is.na(meth)))
meth <- IQR.removal(meth)
log.iqr$NAs.after.IQR3 <- rowSums(is.na(meth))
save(log.iqr, file=paste0(study,".matcaff.logIQR.",timepoint,".Rdata"))

# Create caff.any 
pheno$caff.any <- NA
pheno$caff.any[which(pheno$mat.caff==0)]<-0
pheno$caff.any[which(pheno$mat.caff>0)]<-1

# Generate surrogate variables for technical batch and merge with pheno data to create the phenotype dataframes for the mutually adjusted models
pheno.covs.any <- SVA.generate(meth, pheno, variable.of.interest = "caff.any", model.covariates = c(covariates,"sex"),n.sv=20)
pheno.gest.any <- SVA.generate(meth, pheno, variable.of.interest = "caff.any", model.covariates = c(covariates,"gest.age","sex"),n.sv=20)
pheno.covs.caff <- SVA.generate(meth, pheno, variable.of.interest = "mat.caff", model.covariates = c(covariates,"sex"),n.sv=20)
pheno.gest.caff <- SVA.generate(meth, pheno, variable.of.interest = "mat.caff", model.covariates = c(covariates,"gest.age","sex"),n.sv=20)
pheno.covs.coff <- SVA.generate(meth, pheno, variable.of.interest = "mat.coff", model.covariates = c(covariates,"sex"),n.sv=20)
pheno.gest.coff <- SVA.generate(meth, pheno, variable.of.interest = "mat.coff", model.covariates = c(covariates,"gest.age","sex"),n.sv=20)
pheno.covs.tea <- SVA.generate(meth, pheno, variable.of.interest = "mat.tea", model.covariates = c(covariates,"sex"),n.sv=20)
pheno.gest.tea <- SVA.generate(meth, pheno, variable.of.interest = "mat.tea", model.covariates = c(covariates,"gest.age","sex"),n.sv=20)
pheno.covs.cola <- SVA.generate(meth, pheno, variable.of.interest = "mat.cola", model.covariates = c(covariates,"sex"),n.sv=20)
pheno.gest.cola <- SVA.generate(meth, pheno, variable.of.interest = "mat.cola", model.covariates = c(covariates,"gest.age","sex"),n.sv=20)

#Create the phenotype dataframes for the sex stratified EWASs
pheno.covs.any.boys.only <- pheno.covs.any[which(pheno.covs.any$sex == 0),]
pheno.covs.any.girls.only <- pheno.covs.any[which(pheno.covs.any$sex == 1),]
pheno.gest.any.boys.only <- pheno.gest.any[which(pheno.gest.any$sex == 0),]
pheno.gest.any.girls.only <- pheno.gest.any[which(pheno.gest.any$sex == 1),]
pheno.covs.caff.boys.only <- pheno.covs.caff[which(pheno.covs.caff$sex == 0),]
pheno.covs.caff.girls.only <- pheno.covs.caff[which(pheno.covs.caff$sex == 1),]
pheno.gest.caff.boys.only <- pheno.gest.caff[which(pheno.gest.caff$sex == 0),]
pheno.gest.caff.girls.only <- pheno.gest.caff[which(pheno.gest.caff$sex == 1),]

# Test associations between maternal caffeine and cell types
cells <- pheno.covs.caff[,which(colnames(pheno.covs.caff) %in% cell.names)]
cells.res <- t(apply(cells,2,function(x) summary(lm(x ~ pheno.covs.caff$mat.caff))$coef[2,]))
write.csv(cells.res,file=paste0(study,".matcaff.cells.res.summary.",timepoint,".csv"))

cells <- pheno.covs.any[,which(colnames(pheno.covs.any) %in% cell.names)]
cells.res <- t(apply(cells,2,function(x) summary(lm(x ~ pheno.covs.any$caff.any))$coef[2,]))
write.csv(cells.res,file=paste0(study,".anycaff.cells.res.summary.",timepoint,".csv"))

# Run each EWAS
ewas.res.min.any <- ewas.function(meth, pheno.covs.any[,!colnames(pheno.covs.any) %in% c("sex",covariates)], variable.of.interest = "caff.any")
ewas.res.covs.any <- ewas.function(meth, pheno.covs.any[,!colnames(pheno.covs.any) == "sex"], variable.of.interest = "caff.any")
ewas.res.covs.any.boys.only <- ewas.function(meth, pheno.covs.any.boys.only[,!colnames(pheno.covs.any.boys.only) == "sex"], variable.of.interest = "caff.any")
ewas.res.covs.any.girls.only <- ewas.function(meth, pheno.covs.any.girls.only[,!colnames(pheno.covs.any.girls.only) == "sex"], variable.of.interest = "caff.any")
ewas.res.gest.any <- ewas.function(meth, pheno.gest.any[,!colnames(pheno.gest.any) == "sex"], variable.of.interest = "caff.any")

ewas.res.min.caff <- ewas.function(meth, pheno.covs.caff[,!colnames(pheno.covs.caff) %in% c("sex",covariates)], variable.of.interest = "mat.caff")
ewas.res.covs.caff <- ewas.function(meth, pheno.covs.caff[,!colnames(pheno.covs.caff) == "sex"], variable.of.interest = "mat.caff")
ewas.res.covs.caff.boys.only <- ewas.function(meth, pheno.covs.caff.boys.only[,!colnames(pheno.covs.caff.boys.only) == "sex"], variable.of.interest = "mat.caff")
ewas.res.covs.caff.girls.only <- ewas.function(meth, pheno.covs.caff.girls.only[,!colnames(pheno.covs.caff.girls.only) == "sex"], variable.of.interest = "mat.caff")
ewas.res.gest.caff <- ewas.function(meth, pheno.gest.caff[,!colnames(pheno.gest.caff) == "sex"], variable.of.interest = "mat.caff")

ewas.res.covs.coff <- ewas.function(meth, pheno.covs.coff[,!colnames(pheno.covs.coff) == "sex"], variable.of.interest = "mat.coff")
ewas.res.covs.tea <- ewas.function(meth, pheno.covs.tea[,!colnames(pheno.covs.tea) == "sex"], variable.of.interest = "mat.tea")
ewas.res.covs.cola <- ewas.function(meth, pheno.covs.cola[,!colnames(pheno.covs.cola) == "sex"], variable.of.interest = "mat.cola")

#Add columns with info on caffiene user status
pheno.gest.caff$user <-NA
pheno.gest.caff$user[which(pheno.gest.caff$mat.caff>0)]<-"user"
pheno.gest.caff$user[which(pheno.gest.caff$mat.caff==0)]<-"non user"

# Summarise pheno data and save summaries as .csv files
caff.tableone.stratified <- as.data.frame(print(CreateTableOne(data=pheno.gest.caff[,-1],factorVars=c("mat.smoking","mat.ses","parity","sex"),strata=pheno.gest.caff$user)),stringsAsFactors=FALSE)
caff.tableone <- as.data.frame(print(CreateTableOne(data=pheno.gest.caff[,-1],factorVars=c("mat.smoking","mat.ses","parity","sex")),stringsAsFactors=FALSE))

coff.tableone <- as.data.frame(print(CreateTableOne(data=pheno.gest.coff[,-1],factorVars=c("mat.smoking","mat.ses","parity","sex")),stringsAsFactors=FALSE))
tea.tableone <- as.data.frame(print(CreateTableOne(data=pheno.gest.tea[,-1],factorVars=c("mat.smoking","mat.ses","parity","sex")),stringsAsFactors=FALSE))
cola.tableone <- as.data.frame(print(CreateTableOne(data=pheno.gest.cola[,-1],factorVars=c("mat.smoking","mat.ses","parity","sex")),stringsAsFactors=FALSE))

write.csv(caff.tableone.stratified,file=paste0(study,".matcaff.stratified.summary.",timepoint,".csv"))
write.csv(caff.tableone,file=paste0(study,".matcaff.summary.",timepoint,".csv"))
write.csv(coff.tableone,file=paste0(study,".matcoff.summary.",timepoint,".csv"))
write.csv(tea.tableone,file=paste0(study,".mattea.summary.",timepoint,".csv"))
write.csv(cola.tableone,file=paste0(study,".matcola.summary.",timepoint,".csv"))

# Save EWAS results as an Rdata file
save(list=intersect(ls(),
                    c("ewas.res.min.any",
                      "ewas.res.covs.any",
                      "ewas.res.gest.any",
                      "ewas.res.covs.any.boys.only",
                      "ewas.res.covs.any.girls.only",
                      "ewas.res.min.caff",
                      "ewas.res.covs.caff",
                      "ewas.res.gest.caff",
                      "ewas.res.covs.caff.boys.only",
                      "ewas.res.covs.caff.girls.only",
                      "ewas.res.covs.coff",
                      "ewas.res.covs.tea",
                      "ewas.res.covs.cola"
                      )),
     file=paste0(study,".caffeine.ewasresults.",timepoint,".Rdata"))

