
library(tidyverse)
#install_github("mdhishamshaikh/ViralProduction_R")

####1.0 Creating a Monte Carlo simulation dataframe####

work_dir_simu<- paste0(getwd(), "/results/002_simulation/")
# 
# 
# # Setting the length of experiments
# simu_length<- 1000
# {
# 
# 
#   set.seed(2023)
#   simu_df<- data.frame(matrix(ncol = 5, nrow = 0))
#   colnames(simu_df)<- c("Timepoint", "c_Viruses", #only c_Viruses
#                         "Replicate", "Station_Number",
#                         "Sample_Type")
#   simu_df<- simu_df %>%  mutate(across(!Sample_Type, as.numeric)) %>%
#     mutate(across(Sample_Type, as.character))
# 
# 
#   for (df_no in 1:simu_length){
# 
#     A<- data.frame(Count = runif(n=12, min = 2, max = 3))
#     A$Timepoint<- rep(1:6,2)
#     A$SD_percent<- runif(n = 12, min = 0, max = 100) #UP TO 100% ERROR
#     A$SD<- (A$Count*A$SD_percent)/100
# 
#     #Incorporate this in simulation dataset generation.
#     # generate_positive_numbers <- function(mean_val, std_dev) {
#     #   +     while(TRUE) {
#     #     +         numbers <- rnorm(3, mean = mean_val, sd = std_dev)
#     #     +         if(all(numbers > 0)) return(numbers)
#     #     +     }
# 
#     B<- apply(A, 1, function(x){x[1]+x[4]*abs(scale(rnorm(3)))} )
#     B<- as.data.frame(B)
#     colnames(B)<- rep(1:6,2)
#     B<- B%>% pivot_longer(cols = everything(), names_to = "Timepoint",
#                           values_to = "c_Viruses") %>%
#       mutate(c_Viruses = as.numeric(c_Viruses))
#     B<- B%>%
#       arrange(Timepoint) %>%
#       mutate(Replicate = rep(1:6, 6))%>%
#       mutate(Sample_Type = rep(c(rep("VP",3), rep("VPC",3)),6))
# 
#     B$Station_Number = df_no
# 
# 
#     B<- B%>% mutate(across(!Sample_Type, as.numeric))
# 
#     try(simu_df<- simu_df %>%full_join(B))
#     rm(A,B, df_no)
#   }
# }
# 
# 
# 
# #these columns are required for viralprod but not for the analyses
# cols_to_add<- c("c_Bacteria", "c_HNA", "c_LNA", "c_V1", "c_V2", "c_V3")
# 
# for(col in cols_to_add) {
#   simu_df[[col]] <- runif(nrow(simu_df)) #random number so the fit isn't perfect during calculations
# }
# 
# simu_df[["Location"]] <- "Simulation"
# simu_df[["Depth"]] <- 1
# 
# endpoint_simu_df<- data.frame(Station_Number = 1:simu_length,
#                               Endpoint = sample(3:5, simu_length, replace = T))
# #3-5 because at least 3 points are required, and T6 BEP would be the same as T6 so it is redundant
# 
# #combine this with the simu_df
# simu_df<- full_join(simu_df, endpoint_simu_df)
# 
# #writing it as a csv
# write.csv(simu_df, file = paste0(work_dir_simu, "simulation_dataframe/simulation_df.csv"),
#           row.names = F)
# 


#2.0 Viral production calculations ####

library(viralprod)
work_dir_simu<- paste0(getwd(), "/results/002_simulation/")

#reading simulation csv
simu_df<- read.csv(paste0(work_dir_simu, "simulation_dataframe/simulation_df.csv"))


#need to provide an abundance file. Creating one where we assume abundance to be viral count
#of treatment VP, replicate 1, at timepoint T1

simu_abundance_df<- simu_df %>% dplyr::filter(Timepoint == 1,
                          Replicate ==1,
                          Sample_Type == "VP") %>%
  mutate(Station_Number = as.integer(Station_Number)) %>%
  rename(Total_Bacteria = c_Bacteria,
         Total_Viruses = c_Viruses) %>%
  select(Station_Number, Total_Bacteria, Total_Viruses)


#Some checks before viralprod

vp_check_populations(simu_df)
vp_class_ori_abu(simu_abundance_df) 
class(simu_abundance_df)
class(vp_class_ori_abu(simu_abundance_df))
vp_class_count_data(simu_df)
#passed

#Running viralprod.The third step fails. There is still an unresolved bug in the visualize step, that is no necessary for further analysis
vp_end_to_end(data = simu_df,
              original_abundances = simu_abundance_df,
              #methods = c(1:3,7,10,12),
             #SR_calc = F,
            # BP_endpoint = F,
              write_output = T,
              output_dir = "./results/002_simulation/simulation_viralprod_results")



#3.0 Comparison between methods ####

simu_vp <- read.csv("./results/002_simulation/simulation_viralprod_results/vp_results_ALL.csv")
unique(simu_vp$VP_Method)


#Filtering total viruses, and excluding 'VPC' as it is both lytic and lysogenic.
simu_vp <- simu_vp %>%
  dplyr::filter(Population == "c_Viruses",
                Sample_Type != "VPC") %>%
  mutate(VP_Method = as.factor(VP_Method))


#3.1 between linear methods ####

#Selcting LM methods, and the entire time range (T1_T6)                
lm_vp<- simu_vp %>%
  dplyr::filter(VP_Method == c("LM_ALLPOINTS", "LM_SR_AVG", "LM_AR"),
                Time_Range == "T1_T6")
str(lm_vp)
summary(lm_vp)
unique(lm_vp$VP_Method)
levels(lm_vp$VP_Method)

#Kruskal Wallis

ggpubr::ggboxplot(lm_vp, x = "VP_Method", y = "VP")

kruskal.test(as.numeric(VP) ~ VP_Method, data = lm_vp) #chi-squared = 0.0072784, df = 2, p-value = 0.9964
#no difference between the methods

rstatix::kruskal_effsize(as.numeric(VP) ~ VP_Method, data = lm_vp) #effsize -0.000199
#the difference is minuscule

ggbetweenstats(lm_vp,
               x = VP_Method,
               y = VP,
               type = "nonparametric",
               pairwise.display = "ns")

#3.2 between LM and VIPCAL 

lm_vpcl_vp<- simu_vp %>%
  dplyr::filter(VP_Method == c("LM_ALLPOINTS", "VPCL_AR_DIFF_SE"),
                Time_Range == "T1_T6")
#Kruskal Wallis

ggpubr::ggboxplot(lm_vpcl_vp, x = "VP_Method", y = "VP")

kruskal.test(as.numeric(VP) ~ VP_Method, data = lm_vpcl_vp) #chi-squared = 348.85, df = 1, p-value < 2.2e-16
#SIGNIFICANT difference between the methods

rstatix::kruskal_effsize(as.numeric(VP) ~ VP_Method, data = lm_vpcl_vp) #effsize 0.174 
#the difference is LARGE

ggbetweenstats(lm_vpcl_vp,
               x = VP_Method,
               y = VP,
               type = "nonparametric",
               pairwise.display = "ns")

#3.2 between VIPCAL-SE methods VPCL_AR_DIFF_SE vs VPCL_AR_DIFF_LMER_SE

vpse_vp<- simu_vp %>%
  dplyr::filter(VP_Method == c("VPCL_AR_DIFF_SE", "VPCL_AR_DIFF_LMER_SE"),
                Time_Range == "T1_T6")
#Kruskal Wallis

ggpubr::ggboxplot(vpse_vp, x = "VP_Method", y = "VP")

kruskal.test(as.numeric(VP) ~ VP_Method, data = vpse_vp) #chi-squared = 0.053269, df = 1, p-value = 0.8175
#no difference between the methods

rstatix::kruskal_effsize(as.numeric(VP) ~ VP_Method, data = vpse_vp) #effsize -0.000474
#the difference is minuscule

ggbetweenstats(vpse_vp,
               x = VP_Method,
               y = VP,
               type = "nonparametric",
               pairwise.display = "ns")


#3.3 between VIPCAL vs VIPCAL-SE methods VPCL_AR_DIFF vs VPCL_AR_DIFF_SE

vpcl_vp<- simu_vp %>%
  dplyr::filter(VP_Method == c("VPCL_AR_DIFF", "VPCL_AR_DIFF_SE"),
                Time_Range == "T1_T6",
                VP != 0)
#Kruskal Wallis

ggpubr::ggboxplot(vpcl_vp, x = "VP_Method", y = "VP")

kruskal.test(as.numeric(VP) ~ VP_Method, data = vpcl_vp) #chi-squared = 329.94, df = 1, p-value < 2.2e-16
#significant difference between the methods

rstatix::kruskal_effsize(as.numeric(VP) ~ VP_Method, data = vpcl_vp) #effsize 0.224 
#the difference is large

ggbetweenstats(vpcl_vp,
               x = VP_Method,
               y = VP,
               type = "nonparametric",
               pairwise.display = "ns")
