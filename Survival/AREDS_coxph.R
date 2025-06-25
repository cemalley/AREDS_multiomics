set.seed(20)
library(data.table);library(MASS);library(Hmisc);library(survival)
load('~/Metabolomics/AREDS/Data_prepared_AREDS_release4_AllImputed.RData')

results <- data.table()

for (i in 2:1656){

  id = colnames(data)[i]
  
  mydata <- data.table(PARENT_SAMPLE_NAME=data$PARENT_SAMPLE_NAME, as.data.table(data)[,get(id)])
  mydata <- merge(metadata, mydata, by='PARENT_SAMPLE_NAME')
  colnames(mydata)[ncol(mydata)] <- id
  
  mydata[,start:=min(TIME_POINT), by=CLIENT_SAMPLE_ID]
  mydata[,end:=max(TIME_POINT), by=CLIENT_SAMPLE_ID]
  mydata <- mydata[order(CLIENT_SAMPLE_ID, TIME_POINT)]
  mydata[, yearnum := seq_len(.N), by = CLIENT_SAMPLE_ID]
  removeIDs <- unique(mydata[yearnum==1 & LateAMD_person==1,CLIENT_SAMPLE_ID])
  mydata <- mydata[CLIENT_SAMPLE_ID %nin% removeIDs]
  mydata <- mydata[start != end,]
  mydata$trt_ax <- as.factor(mydata$trt_ax)
  mydata$trt_zc <- as.factor(mydata$trt_zc)
  mydata$trt_axzc <- as.factor(mydata$trt_axzc)
  mydata$trt_pl <- as.factor(mydata$trt_pl)
  mydata$male <- as.character(mydata$male)
  mydata$edu <- as.character(mydata$edu)
  mydata$smoked <- as.character(mydata$smoked)
  
  setDT(mydata)
  
  # Update start and end columns
  mydata[, `:=`(start = c(0, TIME_POINT[-.N]),  # start is lagged TIME_POINT
                end = TIME_POINT), 
         by = CLIENT_SAMPLE_ID]
  
  fit <- coxph(formula(paste0('Surv(time=start, time2=end, event=LateAMD_person) ~ ',id,' + age + male + BMI + edu + smoked + trt_ax + trt_zc + trt_axzc + trt_pl')), data=mydata)
  
  summary(fit)
  
  ORs <- exp(summary(fit)$coefficients[id,1] + qnorm(c(0.025,0.5,0.975)) * summary(fit)$coefficients[id,3])
  
  myresults <- data.table(
    id       = id,
    estimate = signif(coef(summary(fit))[id, "coef"], digits = 4),
    OR = signif(ORs[2], digits = 4),
    CIlow = signif(ORs[1], digits=4),
    CIhigh = signif(ORs[3], digits=4),
    pvalue   = signif(coef(summary(fit))[id, "Pr(>|z|)"], digits = 4),
    age_p = signif(coef(summary(fit))['age', "Pr(>|z|)"], digits = 4),
    sex_p = signif(coef(summary(fit))['male1', "Pr(>|z|)"], digits = 4),
    smk_p = signif(coef(summary(fit))['smoked1', "Pr(>|z|)"], digits = 4),
    bmi_p = signif(coef(summary(fit))['BMI', "Pr(>|z|)"], digits = 4),
    edu2_p = signif(coef(summary(fit))['edu2', "Pr(>|z|)"], digits = 4),
    edu3_p = signif(coef(summary(fit))['edu3', "Pr(>|z|)"], digits = 4),
    trt_ax_p = signif(coef(summary(fit))['trt_ax1', "Pr(>|z|)"], digits = 4),
    trt_zc_p = signif(coef(summary(fit))['trt_zc1', "Pr(>|z|)"], digits = 4),
    trt_axzc_p = signif(coef(summary(fit))['trt_axzc1', "Pr(>|z|)"], digits = 4)
  )
  
  print(paste0(i, '(',id,')'))
  
  results <- rbind(myresults, results)
}

chemical_names$CHEM_ID_NEW <- paste0('CHEM_',chemical_names$CHEM_ID)
results <- merge(results, chemical_names[,c('CHEMICAL_NAME','CHEM_ID','CHEM_ID_NEW', 'SUPER_PATHWAY', 'SUB_PATHWAY','HMDB','KEGG')], by.x='id', by.y='CHEM_ID_NEW', all=T)
results <- results[SUPER_PATHWAY !='Xenobiotics' & SUPER_PATHWAY!='' & SUPER_PATHWAY !='Partially Characterized Molecules',] # 960
results[,padj:=p.adjust(pvalue, method='fdr', n=length(pvalue))]

results[,trt_ax_padj:=p.adjust(trt_ax_p, method='fdr', n=length(trt_ax_p))]
results[,trt_zc_padj:=p.adjust(trt_zc_p, method='fdr', n=length(trt_zc_p))]
results[,trt_axzc_padj:=p.adjust(trt_axzc_p, method='fdr', n=length(trt_axzc_p))]

fwrite(results, '~/Metabolomics/AREDS/Release4/Coxph_lateAAMD_053025/AREDS_CoxPH_LateAMDperson_agemalebismokededutrt_endogenousImputed.csv')

nrow(results) # 960 endogenous
hist(results$pvalue)
hist(results$padj)
nrow(results[pvalue < 0.05,]) # 211
nrow(results[padj < 0.2,]) # 184
nrow(results[padj < 0.05,]) # 77
nrow(results[padj < 0.05 & SUPER_PATHWAY !='' & SUPER_PATHWAY %nin% c("Partially Characterized Molecules",'Xenobiotics')]) # 77


# get missingness per metabolite-----
# Calculate missingness percent for columns 2:1656
missingness <- data[, lapply(.SD, function(x) mean(is.na(x)) * 100), .SDcols = 2:1656]

# Convert to long format
missingness_long <- melt(missingness, measure.vars = names(missingness),
                         variable.name = "Metabolite", value.name = "MissingPercent")

head(missingness_long)
