set.seed(20)
library(data.table);library(MASS);library(Hmisc);library(survival)
setwd('~/Metabolomics/AREDS2/')
load('AREDS2_data_prepared.RData')

results <- data.table()

for (i in 2:1656){
  
  id = colnames(data)[i]
  
  mydata <- data.table(PARENT_SAMPLE_NAME=data$PARENT_SAMPLE_NAME, as.data.table(data)[,get(id)])
  mydata <- merge(metadata, mydata, by='PARENT_SAMPLE_NAME')
  colnames(mydata)[ncol(mydata)] <- id
  
  mydata[,start:=min(TIME_POINT), by=CLIENT_SAMPLE_ID]
  mydata[,end:=max(TIME_POINT), by=CLIENT_SAMPLE_ID]
  mydata <- mydata[order(CLIENT_SAMPLE_ID, TIME_POINT)]
  mydata <- mydata[start!=end & !is.na(edu) & !is.na(age), ]
  mydata$treat <- as.factor(mydata$treat)
  mydata$male <- as.character(mydata$male)
  mydata$edu <- as.character(mydata$edu)
  mydata$smoked <- as.character(mydata$smoked)
  
  
  # Convert to data.table for easy manipulation
  setDT(mydata)
  
  # Update start and end columns
  mydata[, `:=`(start = c(0, TIME_POINT[-.N]),  # start is lagged TIME_POINT
                end = TIME_POINT), 
         by = CLIENT_SAMPLE_ID]
  
  
  fit <- coxph(formula(paste0('Surv(time=start, time2=end, event=LateAMD_person) ~ ',id,' + age + male + edu + smoked + treat')), data=mydata)
  
  summary(fit)
  
  ORs <- exp(summary(fit)$coefficients[id,'coef'] + qnorm(c(0.025,0.5,0.975)) * summary(fit)$coefficients[id,'se(coef)'])
  
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
    edu2_p = signif(coef(summary(fit))['edu2', "Pr(>|z|)"], digits = 4),
    edu3_p = signif(coef(summary(fit))['edu3', "Pr(>|z|)"], digits = 4),
    trt_lutzea_p = signif(coef(summary(fit))['treat2', "Pr(>|z|)"], digits = 4),
    trt_dhaepa_p = signif(coef(summary(fit))['treat3', "Pr(>|z|)"], digits = 4),
    trt_combo_p = signif(coef(summary(fit))['treat4', "Pr(>|z|)"], digits = 4)
  )
  
  print(paste0(i, '(',id,')'))
  
  results <- rbind(myresults, results)
}


results <- merge(results, chemical_names[,c('CHEMICAL_NAME','CHEM_ID','CHEM_ID_NEW', 'SUPER_PATHWAY', 'SUB_PATHWAY','HMDB','KEGG')], by.x='id', by.y='CHEM_ID_NEW', all=T)
#results <- results[SUPER_PATHWAY !='Xenobiotics' & !is.na(SUPER_PATHWAY) & SUPER_PATHWAY !='Partially Characterized Molecules',]
results[,padj:=p.adjust(pvalue, method='fdr', n=length(pvalue))]

fwrite(results, '~/Metabolomics/AREDS2/Coxph_lateAMD/AREDS2_CoxPH_LateAMDperson_agemalesmokededutreat_ALL.csv')
results <- fread('~/Metabolomics/AREDS2/Coxph_lateAMD/AREDS2_CoxPH_LateAMDperson_agemalesmokededutreat_ALL.csv')
hist(results$pvalue)
hist(results$padj)
nrow(results) # 948; 1655
nrow(results[pvalue < 0.05,]) # 251; 392
nrow(results[padj < 0.2,]) # 262; 376
nrow(results[padj < 0.05,]) # 99; 103

# pop stats----
metadata[,start:=min(TIME_POINT) & AnyAMD_person==0, by=CLIENT_SAMPLE_ID]
metadata[,end:=max(TIME_POINT) & AnyAMD_person==1, by=CLIENT_SAMPLE_ID]
metadata <- metadata[order(CLIENT_SAMPLE_ID, TIME_POINT)]
metadata[, yearnum := seq_len(.N), by = CLIENT_SAMPLE_ID]
mean(metadata$yearnum) # 2.202368
range(metadata$yearnum)# up to 5 years
sd(metadata$yearnum) # 1.028184

people <- unique(metadata[,c('CLIENT_SAMPLE_ID', 'TIME_POINT' , 'yearnum','start', 'end', 'AnyAMD_person','IntAMD_person','LateAMD_person', 'age', 'male', 'edu','smoked', 'treat')])

people <- setorder(setDT(people), -yearnum)[, head(.SD, 1), keyby = CLIENT_SAMPLE_ID]

# find the numbers/% of eyes that progressed in each of the 3 categories

# Total number of people
total_people <- nrow(people)

# Calculate numbers, percentages, and standard deviations for each category
categories <- c("AnyAMD_person", "IntAMD_person", "LateAMD_person")
results <- data.frame(Category = categories, Number = NA, Percentage = NA, SD_Percentage = NA)

for (i in 1:length(categories)) {
  category <- categories[i]
  num_progressed <- sum(people[[category]])
  perc_progressed <- (num_progressed / total_people) * 100
  
  # Standard deviation calculation for binary data
  p <- num_progressed / total_people
  sd_percentage <- sqrt(p * (1 - p) / total_people) * 100
  
  # Store results
  results$Number[i] <- num_progressed
  results$Percentage[i] <- perc_progressed
  results$SD_Percentage[i] <- sd_percentage
}

# Display the results
print(results)

# expanded table with the other covariables----
# Load necessary library
library(dplyr)

# Calculate summary statistics
summary_table <- people %>%
  summarise(
    # Male Sex
    Male = paste0(sum(male == 1), " (", round(mean(male) * 100, 2), "%)"),
    Female = paste0(sum(male == 0), " (", round((1 - mean(male)) * 100, 2), "%)"),
    
    # Education Levels
    High_School_or_Less = paste0(sum(edu == 1, na.rm=T), " (", round(sum(edu == 1, na.rm=T) / n() * 100, 2), "%)"),
    Some_college = paste0(sum(edu == 2, na.rm=T), " (", round(sum(edu == 2, na.rm=T) / n() * 100, 2), "%)"),
    Bach_higher = paste0(sum(edu == 3, na.rm=T), " (", round(sum(edu == 3, na.rm=T) / n() * 100, 2), "%)"),
    
    # Smoking
    Smoked_Yes = paste0(sum(smoked == 1), " (", round(mean(smoked) * 100, 2), "%)"),
    Smoked_No = paste0(sum(smoked == 0), " (", round((1 - mean(smoked)) * 100, 2), "%)"),
    
    # Treatment Categories
    Treatment_1 = paste0(sum(treat == 1), " (", round(sum(treat == 1) / n() * 100, 2), "%)"),
    Treatment_2 = paste0(sum(treat == 2), " (", round(sum(treat == 2) / n() * 100, 2), "%)"),
    Treatment_3 = paste0(sum(treat == 3), " (", round(sum(treat == 3) / n() * 100, 2), "%)"),
    Treatment_4 = paste0(sum(treat == 4), " (", round(sum(treat == 4) / n() * 100, 2), "%)")
  )

# Transpose the summary table for better readability
summary_table_t <- as.data.frame(t(summary_table))
colnames(summary_table_t) <- "Count (Percentage)"

# Display the table
summary_table_t
