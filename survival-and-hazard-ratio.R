# Human pluripotent stem cells identify molecular targets of trisomy 12 in chronic lymphocytic leukemia patients
# Jennifer C. Reid et al. Cell Reports, 2021

# Load required packages
install.packages('survival')
install.packages('survminer')
install.packages('dplyr')
library(survival)
library(survminer)
library(dplyr)

# SURVIVAL ANALYSIS
# Load data for IGHV-U (N=22) patients used for analysis 
X <- read.delim("IGHV-U.csv", sep=",", header=TRUE, stringsAsFactors = FALSE, row.names=1)
glimpse(X)
#Rows: 22
#$ DISEASE                	<chr> "cll", "cll", "cll", "cll", "cll", "cll", "cll", "mbl", "mbl", "mbl", "cll",…
#$ GENETICS               	<chr> "trisomy12", "normal", "trisomy12", "trisomy12", "normal",…
#$ IGHV                   		<int> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
#$ Progression_0no_1yes   	<int> 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1…
#$ Months_Til_Progression	<dbl> 42.5753, 29.1945, 11.4740, 31.5288, 70.3562, 35.9671, 72.4…
#$ LE_group               	<chr> "LE_above", "LE_above", "LE_above", "LE_below", "LE_below"…
#surv_objectX <- Surv(time = X$Months_Til_Progression, event = X$Progression_0no_1yes)
fit1X <- survfit(surv_objectX ~ LE_group, data = X)
ggsurvplot(fit1X, data = X, pval = TRUE, risk.table=TRUE, surv.median.line = "none", xlim = c(0,48), break.time.by = 12)
ggsave("./survplot.pdf", width = 4.5, height = 5)

fit2 <- survfit(surv_objectX ~ GENETICS, data = X)
ggsurvplot(fit2, data = X, pval = TRUE, xlim = c(0,48), break.time.by = 12)
ggsave("./survplot_cytogenetics.pdf", width = 4.5, height = 5)

fit3 <- survfit(surv_objectX ~ DISEASE, data = X)
ggsurvplot(fit3, data = X, pval = TRUE, xlim = c(0,48), break.time.by = 12)
ggsave("./survplot_disease.pdf", width = 4.5, height = 5)

# COX PROPORTIONAL HAZARDS MODEL
# an HR > 1 indicates an increased risk of death

fit.coxph <- coxph(surv_objectX ~ DISEASE + GENETICS + LE_group, data = X)
ggforest(fit.coxph, data = X)