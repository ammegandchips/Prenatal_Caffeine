# Script written by Laura Schellhas
# November 2019

######################################################################################################### 
#                                             								                                          #
# 			                      DMRff analysis for maternal caffeine EWAS   				                      #
#					      							                                                                        	#			
#########################################################################################################

# For this analysis you will need:
#1 the maternal caffeine EWAS results from your cohort named 'YOURSTUDY.caffeine.ewasresults.birth.Rdata'
#2 the filtered annotation data file that I have sent you called '450k.annotation.filtered.Rdata'
#3 the methylation data from your cohort
#4 the phenofile that you used for your caffeine EWAS analysis

############################################### 
# 1) Set the initial parameters 
###############################################
#change to your study identifier
study <- "ALSPAC" 

###############################################
# 2) Load required packages & functions 
###############################################
#devtools
#install.packages('devtools')
library(devtools)

#dmrff 
# if you are installing the dmrff package for the first time run:
# if (!require(dmrff)) {
#   library(devtools)
#   install_github("perishky/dmrff")
#   library(dmrff)
# }
library(dmrff)

#matrixStats
#install.packages('matrixStats')
library(matrixStats)

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

##################################################################################################################################################################
#3) load the ewas summary statistics (that is the EWAS results file from your cohort, e.g. the ewas results for alspac: "alspac.caffeine.ewasresults.birth.Rdata") 
##################################################################################################################################################################
load(paste0(study, ".caffeine.ewasresults.birth.Rdata"))

##########################################################
#4) load the annotation Rdata file that I have sent you
##########################################################
load('450k.annotation.filtered.RData')

###############################################
#5) load your cohort's methylation data (meth)
###############################################
meth <- meth.matrix.cord

#####################################################################
#6) load the phenofile that you used for your caffeine EWAS analysis
#####################################################################
pheno <- complete.pheno.cord

# Remove missing data (sample.id = methylation ID)
pheno <- na.omit(pheno[c("sample.id","mat.bmi","sex","mat.ses","mat.age","mat.smoking","parity","gest.age")])

#####################################################################
#7) Match meth to only people included in the phenofile
#####################################################################
meth<-meth[,match(pheno$sample.id,colnames(meth))]

#######################################################################################
#8) IQR*3 method to remove outliers (if this has not already been applied to your data)
#######################################################################################
meth <- IQR.removal(meth)

###############################################
#9) Define EWAS models
###############################################
models<-list(ewas.res.covs.caff, 
             ewas.res.covs.caff.boys.only,
             ewas.res.covs.caff.girls.only,
             ewas.res.covs.coff,
             ewas.res.covs.tea,
             ewas.res.covs.cola,
             ewas.res.covs.any)

names(models)<-c('ewas.res.covs.caff', 
                 'ewas.res.covs.caff.boys.only',
                 'ewas.res.covs.caff.girls.only',
                 'ewas.res.covs.coff',
                 'ewas.res.covs.tea',
                 'ewas.res.covs.cola',
                 'ewas.res.covs.any')

###############################################
#10) Add the annotation data to the models
###############################################
# add CpG name as column to models 
cpg.name <- function(df) {
  df$name <- rownames(df)
  df
}
models<- lapply(models, cpg.name)

# merge models with annotation data by cpg name column
merge.anno<-function(df) {
  merge(df, y=annotation, by='name',all.y=T)
}
models<-lapply(models,merge.anno)

######################################################################################
#11) create idx to sort models after genomic position (idx is the same for all models)
######################################################################################
idx <- order(models$ewas.res.covs.caff$chromosome, models$ewas.res.covs.caff$position)

###############################################
#12) run dmrff
###############################################
dmrff.res<-function(df,methylation.matrix){
  dmrff.pre(df$coef[idx],df$se[idx],methylation.matrix[idx,], df$chromosome[idx], 
            df$position[idx], maxsize=20, verbose=T)
}

# sex stratified models
#please be careful to select on sex based on how males and females have been coded in your EWAS analysis (e.g. sex= 0/1 or "male"/"female") 
YOURSTUDY.dmrff.res.boys <- dmrff.res(df=models$ewas.res.covs.caff.boys.only, methylation.matrix=meth[,pheno$sex==0]) #male
YOURSTUDY.dmrff.res.girls <- dmrff.res(df=models$ewas.res.covs.caff.girls.only, methylation.matrix=meth[,pheno$sex==1]) #female
YOURSTUDY.dmrff.res.sex<-list("YOURSTUDY.dmrff.res.boys"=YOURSTUDY.dmrff.res.boys,"YOURSTUDY.dmrff.res.girls"=YOURSTUDY.dmrff.res.girls)

#other models
YOURSTUDY.dmrff.res.others<-lapply(models[-which(names(models)%in% c("ewas.res.covs.caff.girls.only", "ewas.res.covs.caff.boys.only"))], dmrff.res, methylation.matrix=meth)

#all models
YOURSTUDY.dmrff.res <- c(YOURSTUDY.dmrff.res.others,YOURSTUDY.dmrff.res.sex)

###############################################
#13) save results as Rdata file
###############################################
save(YOURSTUDY.dmrff.res, file=paste0(study,'.dmrff.res.RData'))

##############################################
# 		        THANK YOU !!!		               #
##############################################
