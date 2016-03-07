###########################################################################
#
#   Immuneering Corporation
#
#   SOFTWARE COPYRIGHT NOTICE AGREEMENT
#   This software and its documentation are copyright (2014) by the
#   Immuneering Corporation. All rights are reserved.
#
#   This software is supplied without any warranty or guaranteed support
#   whatsoever. Immuneering Corporation cannot be responsible for its use,
#   misuse, or functionality.
#
#
###########################################################################

library(utils)
library(lme4)
library(limma)
library(sva)


getLabels <- function() {
  download.file(url='ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE61nnn/GSE61901/matrix/GSE61901_series_matrix.txt.gz',destfile = 'GSE61901_series_matrix.txt.gz',mode = 'wb')
  data <- read.delim(file='GSE61901_series_matrix.txt.gz',header=F,skip=32)
  SAMPLE <- as.character(unlist(data[2,-1]));
  CLASS <- as.character(unlist(data[26,-1]));
  SLOT <- as.character(unlist(data[12,-1]));
  SLOT <- gsub(pattern = 'array_address: ',replacement = '',x = SLOT)
  BATCH <- as.character(unlist(data[10,-1]));
  BATCH <- gsub(pattern = 'batch: ','B',BATCH)
  data.frame(SAMPLE=factor(SAMPLE),CLASS=factor(CLASS),SLOT=factor(SLOT),BATCH=factor(BATCH))
}

download.file(url = 'http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE61901&format=file&file=GSE61901%5Fnon%2Dnormalized%2Etxt%2Egz',destfile = 'GSE61901_non-normalized.txt.gz',mode = 'wb')
data <- read.delim(file='GSE61901_non-normalized.txt.gz',header=T,skip=5,row.names=1,as.is=T,check.names=F)
data <- data[,!(colnames(data) %in% c('Detection Pval'))]
data <- normalizeQuantiles(data)

labels <- getLabels()
labels <- labels[match(colnames(data),labels$SLOT),]
labels$chip <- factor(substr(as.character(labels$SLOT),1,nchar(as.character(labels$SLOT))-2))
combat_data <- ComBat(data,batch = labels$BATCH,mod = model.matrix(~labels$CLASS))
combat_data_averaged <- avearrays(combat_data,ID = labels$chip)
labels_averaged <- unique(labels[,c(2,4,5)]);
data_averaged <- avearrays(data,ID = labels$chip)

ranova <- function(X,lab,batch,sampid) {
  model <- aov(X~batch+lab+Error(1/sampid))
  unlist(summary(model))[['Error: Within.Pr(>F)2']]
}

oneanova <- function(X,lab) {
  lab <- factor(as.character(lab))
  model <- glm(X~lab)
  coef(summary(model))[2,'Pr(>|t|)']
}

twowayanova_blocking <- function(X,batch,lab) {
  batch <- factor(as.character(batch));
  lab <- factor(as.character(lab));
  model <- glm(X~batch+lab)
  coef(summary(model))[18,'Pr(>|t|)']
}

lmem <- function(X,lab,batch,sampid) {
  # Need to use ML estimates since REML estimates are not valid
  # when fixed effects change
  # http://stats.stackexchange.com/questions/41123/reml-vs-ml-stepaic
  model1 <- lmer(X~0+batch+(1|sampid),REML = F)
  model2 <- lmer(X~0+batch+lab+(1|sampid),REML=F);
  comp <- anova(model1,model2,test='Chisq');
  comp[['Pr(>Chisq)']][2]
}

limma_blocking <- function(inmat,lab) {
  class <- factor(make.names(as.character(lab$CLASS)))
  batch <- factor(make.names(as.character(lab$BATCH)))
  design <- model.matrix(~0+class+batch)
  fit <- lmFit(inmat,design)
  cm <- makeContrasts(DPvsQ=classGA.DP-classGA.Q,levels=design)
  fit2 <- contrasts.fit(fit, cm)
  limma.fit2 <- eBayes(fit2)
  limma.fit2$p.value
}

limma_noblocking <- function(inmat,lab) {
  class <- factor(make.names(as.character(lab$CLASS)))
  design <- model.matrix(~0+class)
  fit <- lmFit(inmat,design)
  cm <- makeContrasts(DPvsQ=classGA.DP-classGA.Q,levels=design)
  fit2 <- contrasts.fit(fit, cm)
  limma.fit2 <- eBayes(fit2)
  limma.fit2$p.value
}

limma_dupcorr <- function(inmat,lab) {
  batch <- factor(make.names(as.character(lab$BATCH)));
  class <- factor(make.names(as.character(lab$CLASS)))
  samearray <- lab$chip
  design <- model.matrix(~0+class+batch)
  corfit <- duplicateCorrelation(inmat,design,block=samearray)
  fit <- lmFit(inmat,design,block=samearray,correlation=corfit$consensus)
  cm <- makeContrasts(DPvsQ=classGA.DP-classGA.Q,levels=design)
  fit2 <- contrasts.fit(fit, cm)
  limma.fit2 <- eBayes(fit2)
  limma.fit2$p.value
}

write.batch.sample.distribution <- function() {
  write.csv(as.matrix(table(unique_samples$labels.BATCH,unique_samples$labels.CLASS)),file='/workspace/copaxone/stage2/rna/ur_02112016_nygaard_data_balance/sample_table.csv',quote=F)
}

dp_vs_q <- list();
dp_vs_q$mat <- data[,labels$CLASS %in% c('GA DP','GA Q')]
dp_vs_q$labels <- labels[labels$CLASS %in% c('GA DP','GA Q'),]

ranova_out <- apply(dp_vs_q$mat,1,ranova,lab=dp_vs_q$labels$CLASS,batch=dp_vs_q$labels$BATCH,sampid=dp_vs_q$labels$chip)
lme4_out <- apply(dp_vs_q$mat,1,lmem,lab=dp_vs_q$labels$CLASS,batch=dp_vs_q$labels$BATCH,sampid=dp_vs_q$labels$chip)
limma_out <- limma_dupcorr(dp_vs_q$mat,lab = dp_vs_q$labels)

dp_vs_q_averaged <- list();
dp_vs_q_averaged$mat <- combat_data_averaged[,labels_averaged$CLASS %in% c('GA DP','GA Q')]
dp_vs_q_averaged$labels <- labels_averaged[labels_averaged$CLASS %in% c('GA DP','GA Q'),]

dp_vs_q_nocombat_averaged <- list();
dp_vs_q_nocombat_averaged$mat <- data_averaged[,labels_averaged$CLASS %in% c('GA DP','GA Q')]
dp_vs_q_nocombat_averaged$labels <- labels_averaged[labels_averaged$CLASS %in% c('GA DP','GA Q'),]

limma_noblocking_out <- limma_noblocking(inmat = dp_vs_q_averaged$mat,lab = dp_vs_q_averaged$labels)
oneanova_out <- apply(dp_vs_q_averaged$mat,1,oneanova,lab=dp_vs_q_averaged$labels$CLASS)

limma_blocking_out <- limma_blocking(inmat = dp_vs_q_nocombat_averaged$mat,lab = dp_vs_q_nocombat_averaged$lab)
twowayanova_blocking_out <- apply(dp_vs_q_nocombat_averaged$mat,1,twowayanova_blocking,batch=dp_vs_q_nocombat_averaged$lab$BATCH,lab=dp_vs_q_nocombat_averaged$lab$CLASS)


result_summary <- data.frame("Quantile Normalization, no averaging of technical replicates, utilize CHIP as blocking variable"=c(length(which(p.adjust(limma_out,method='BH') < 0.05)),NA,length(which(p.adjust(lme4_out,method='BH') < 0.05)),length(which(p.adjust(ranova_out,method='BH') < 0.05))),
           "Quantile normalization (averaging technical replicates), apply ComBat including treatments as covariates"=c(length(which(p.adjust(limma_noblocking_out,method='BH') < 0.05)),length(which(p.adjust(oneanova_out,method='BH') < 0.05)),NA,NA),
           "Quantile normalization (averaging technical replicates), utilize CHIP as blocking variable"=c("11 (Nygaard et al)",length(which(p.adjust(twowayanova_blocking_out,method='BH') < 0.05)),NA,NA))

print(result_summary)