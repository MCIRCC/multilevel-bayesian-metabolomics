##########
##Purpose: plotting for MIDAS variable importance paper
##Date: 03/23/20
##########

require(tidyverse)
require(ggpubr)

#load data. three main sources of data for MIDAS project
#Simulated data with no missingness (nm), Simulated data with missingess (m), 'real data' analysis
no_missing_file = 'plots_results_no_missing.csv'
miss_file ='plots_results_with_missing.csv'
df_no_missing = read.csv(no_missing_file)
df_miss = read.csv(miss_file)


### True positive rate (TPR) = power ###
### False discovery rate (FDR) ###
### F1 = harmonic mean of precision and recall ###
### Average magnitude error measures how different actual vs. predicted Beta coefficient ###


#No-missingness plots. 
#Completed for simulation when the fraction of metabolites signficant == 0.4
frac_sig = 0.4

#Power = TPR @ Frac_Sig of 0.4. No missingness
pA = ggplot(df_no_missing %>% filter (FRAC_SIG == frac_sig)) + geom_line(aes(x=N_SAMPLE,y=tpr,color=label))
pA = pA + theme_classic() + 
  labs(x= "Sample Size per group",
       y= "\n Power",
       color = "Model") + scale_color_manual(labels = c("Bayes", "Benjamini–Hochberg", "Bonferonni", "Raw"),
                                               values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7")) + 
  scale_x_continuous(breaks=c(25, 50, 75, 100, 125, 150)) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
  
print(pA)

#False discovery rate = FDR @ Frac_Sig of 0.4. No missingness
pB = ggplot(df_no_missing %>% filter (FRAC_SIG == frac_sig)) + geom_line(aes(x=N_SAMPLE,y=fdr,color=label))
pB = pB + theme_classic() + 
  labs(x= "Sample Size per group",
       y= "False \n Discovery Rate",
       color = "Model") + scale_color_manual(labels = c("Bayes", "Benjamini–Hochberg", "Bonferonni", "Raw"),
                                             values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7")) + 
  geom_hline(yintercept=0.05, linetype="dashed", color = "black") + 
  scale_x_continuous(breaks=c(25, 50, 75, 100, 125, 150))
print(pB)

#F1. Harmonic Mean of Precision and Recall
#Mulitiplying orignial value calculated by Chris. May need to remove in future scripting
# pC = ggplot(df_no_missing %>% filter (FRAC_SIG == frac_sig)) + geom_line(aes(x=N_SAMPLE,y=2*f1,color=label))
# pC = pC + theme_classic() + 
#   labs(x= "Sample Size per group",
#        y= "F1 Statistic",
#        color = "Model") + scale_color_manual(labels = c("Bayes", "Benjamini–Hochberg", "Bonferonni", "Raw"),
#                                              values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7")) 
# print(pC)

#Average Exaggeration Ratio
pC = ggplot(df_no_missing %>% filter (FRAC_SIG == frac_sig)) + geom_line(aes(x=N_SAMPLE,y=avg_mag_error,color=label))
pC = pC + theme_classic() + 
  labs(x= "Sample Size per group",
       y= "Average \n Exaggeration Ratio",
       color = "Model") + scale_color_manual(labels = c("Bayes", "Benjamini–Hochberg", "Bonferonni", "Raw"),
                                             values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7")) + 
  geom_hline(yintercept = 1, linetype='dashed', color = 'black') +
  scale_y_continuous(breaks=c(0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)) + 
  scale_x_continuous(breaks=c(25, 50, 75, 100, 125, 150))
print(pC)

fig_nm = ggarrange(pA, pB, pC, labels = c("A", "B", "C"),
          common.legend = TRUE, legend = "top",
          ncol = 1, nrow = 3)

ggsave("figure_nm.tiff", plot = last_plot(),scale = 0.85, 
       width = 5.5, height = 8.5, units = c("in", "cm", "mm"),
       dpi = 300)

#Plots with missing data introduced into the simulation
#Power (TPR) vs. missing rate
qA = ggplot(df_miss) + geom_line(aes(x=MISSING_RATE,y=tpr,color=label))
qA = qA + theme_classic() + 
  labs(x= "Missing Rate", y= "\n Power", color = "Model") + 
  scale_color_manual(labels = c("Bayes", "Benjamini–Hochberg", "Bonferonni", "Raw"),
                                             values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
print(qA)

#FDR vs. missing rate
qB = ggplot(df_miss) + geom_line(aes(x=MISSING_RATE,y=fdr,color=label))
qB = qB + theme_classic() + 
  labs(x= "Missing Rate", y= "False \n Discovery Rate", color = "Model") + 
  scale_color_manual(labels = c("Bayes", "Benjamini–Hochberg", "Bonferonni", "Raw"),
                     values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7")) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "black") 
print(qB)

#F1 stat vs. missing rate
# qC = ggplot(df_miss) + geom_line(aes(x=MISSING_RATE,y=2*f1,color=label))
# qC = qC + theme_classic() + 
#   labs(x= "Missing Rate", y= "F1 Statistic", color = "Model") + 
#   scale_color_manual(labels = c("Bayes", "Benjamini–Hochberg", "Bonferonni", "Raw"),
#                      values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7")) 
# print(qC)

#Average mag error vs. missing rate
qC = ggplot(df_miss) + geom_line(aes(x=MISSING_RATE,y=avg_mag_error,color=label))
qC = qC + theme_classic() + 
  labs(x= "Missing Rate", y= "Average \n Exaggeration Ratio", color = "Model") + 
  scale_color_manual(labels = c("Bayes", "Benjamini–Hochberg", "Bonferonni", "Raw"),
                     values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7")) +
  geom_hline(yintercept = 1, linetype='dashed', color = 'black') + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),
                     limits = c(0.9,1.65), breaks=c(1.0, 1.2, 1.4, 1.6))
  
print(qC)

fig_q = ggarrange(qA, qB, qC, labels = c("A", "B", "C"),
                   common.legend = TRUE, legend = "top",
                   ncol = 1, nrow = 3)

ggsave("figure_miss.tiff", plot = last_plot(),scale = 0.85, 
       width = 5.5, height = 8.5, units = c("in", "cm", "mm"),
       dpi = 300)


#Supplementary figure? Started as S1 --> S2
#How do power, FDR, F1, and Avg_Mag_Error change with different fraction of metaoblites significant 
#set sample size to max = 150 
#power
n_samp = 150 
sA = ggplot(df_no_missing %>% filter (N_SAMPLE == n_samp)) + geom_line(aes(x=FRAC_SIG,y=tpr,color=label))
sA = sA + theme_classic() + 
  labs(x= "Fraction of Signficant Metabolites",
       y= "\n Power",
       color = "Model") + scale_color_manual(labels = c("Bayes", "Benjamini–Hochberg", "Bonferonni", "Raw"),
                                             values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7")) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
print(sA)

#FDR
sB = ggplot(df_no_missing %>% filter (N_SAMPLE == n_samp)) + geom_line(aes(x=FRAC_SIG,y=fdr,color=label))
sB = sB + theme_classic() + 
  labs(x= "Fraction of Signficant Metabolites",
       y= "False \n Discovery Rate",
       color = "Model") + scale_color_manual(labels = c("Bayes", "Benjamini–Hochberg", "Bonferonni", "Raw"),
                                             values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7")) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "black") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
  
print(sB)

#F1
# sC = ggplot(df_no_missing %>% filter (N_SAMPLE == n_samp)) + geom_line(aes(x=FRAC_SIG,y=2*f1,color=label))
# sC = sC + theme_classic() + 
#   labs(x= "Fraction of \n Signficant Metabolites",
#        y= "F1 Statistic",
#        color = "Model") + scale_color_manual(labels = c("Bayes", "Benjamini–Hochberg", "Bonferonni", "Raw"),
#                                              values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7")) 
# print(sC)

#Avg_Mag_Error
sC = ggplot(df_no_missing %>% filter (N_SAMPLE == n_samp)) + geom_line(aes(x=FRAC_SIG,y=avg_mag_error,color=label))
sC = sC + theme_classic() + 
  labs(x= "Fraction of Signficant Metabolites",
       y= "Average \n Exaggeration Ratio",
       color = "Model") + scale_color_manual(labels = c("Bayes", "Benjamini–Hochberg", "Bonferonni", "Raw"),
                                             values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7")) +
  geom_hline(yintercept = 1, linetype='dashed', color = 'black') +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),
                     limits = c(0.70,1.15), breaks=c(0.8, 0.9, 1.0, 1.1))
print(sC)

#
fig_s = ggarrange(sA, sB, sC, labels = c("A", "B", "C"),
                  common.legend = TRUE, legend = "top",
                  ncol = 1, nrow = 3)

ggsave("figure_s1.tiff", plot = last_plot(),scale = 0.85, 
       width = 5.5, height = 8.5, units = c("in", "cm", "mm"),
       dpi = 300)




#this started as fig s2 but became s1
#frac_sig =0.7; no missingness
frac_sig_new =0.7

tA = ggplot(df_no_missing %>% filter (FRAC_SIG == frac_sig_new)) + geom_line(aes(x=N_SAMPLE,y=tpr,color=label))
tA = tA + theme_classic() + 
  labs(x= "Sample Size per group",
       y= "\n Power",
       color = "Model") + scale_color_manual(labels = c("Bayes", "Benjamini–Hochberg", "Bonferonni", "Raw"),
                                             values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7")) + 
  scale_x_continuous(breaks=c(25, 50, 75, 100, 125, 150)) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))

print(tA)

tB = ggplot(df_no_missing %>% filter (FRAC_SIG == frac_sig_new)) + geom_line(aes(x=N_SAMPLE,y=fdr,color=label))
tB = tB + theme_classic() + 
  labs(x= "Sample Size per group",
       y= "False \n Discovery Rate",
       color = "Model") + scale_color_manual(labels = c("Bayes", "Benjamini–Hochberg", "Bonferonni", "Raw"),
                                             values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7")) + 
  geom_hline(yintercept=0.05, linetype="dashed", color = "black") + 
  scale_x_continuous(breaks=c(25, 50, 75, 100, 125, 150))
print(tB)

tC = ggplot(df_no_missing %>% filter (FRAC_SIG == frac_sig_new)) + geom_line(aes(x=N_SAMPLE,y=avg_mag_error,color=label))
tC = tC + theme_classic() + 
  labs(x= "Sample Size per group",
       y= "Average \n Exaggeration Ratio",
       color = "Model") + scale_color_manual(labels = c("Bayes", "Benjamini–Hochberg", "Bonferonni", "Raw"),
                                             values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7")) + 
  geom_hline(yintercept = 1, linetype='dashed', color = 'black') +
  scale_y_continuous(breaks=c(0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)) + 
  scale_x_continuous(breaks=c(25, 50, 75, 100, 125, 150))
print(tC)

fig_s2 = ggarrange(tA, tB, tC, labels = c("A", "B", "C"),
                   common.legend = TRUE, legend = "top",
                   ncol = 1, nrow = 3)

ggsave("figure_s2.tiff", plot = last_plot(),scale = 0.85, 
       width = 5.5, height = 8.5, units = c("in", "cm", "mm"),
       dpi = 300)
