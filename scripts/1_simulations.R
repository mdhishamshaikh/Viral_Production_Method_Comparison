#Aim: To generate simulation dataset to comparatively assess viral production analyses methods.



#We'll use the viral production being developed on Github for analyses.
#Input as a simualtion dataframe catered to the simualtion function.



source("./scripts/viral_production_Step2.R")
#I will have to remove some input functions from the second script for now.


#Creating a simulation dataframe

####1.0 Creating a dataframe####
# #Set number of dataframes you'd like to create
simu_length<- 1000
{
  
  
  set.seed(2023)
  simu_df<- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(simu_df)<- c("Timepoint", "c_Viruses", "Replicate", "Station_Number", "Sample_Type")
  simu_df<- simu_df %>%  mutate(across(!Sample_Type, as.numeric)) %>%
    mutate(across(Sample_Type, as.character))
  
  
  for (df_no in 1:simu_length){
    
    A<- data.frame(Count = runif(n=12, min = 2, max = 3))
    A$Timepoint<- rep(1:6,2)
    A$SD_percent<- runif(n = 12, min = 0, max = 100) #UP TO 100% ERROR
    A$SD<- (A$Count*A$SD_percent)/100
    
   
    B<- apply(A, 1, function(x){x[1]+x[4]*abs(scale(rnorm(3)))} )
    B<- as.data.frame(B)
    colnames(B)<- rep(1:6,2)
    B<- B%>% pivot_longer(cols = everything(), names_to = 'Timepoint',
                          values_to = 'c_Viruses')
    B<- B%>%
      arrange(Timepoint) %>%
      mutate(Replicate = rep(1:6, 6))%>%
      mutate(Sample_Type = rep(c(rep('VP',3), rep('VPC',3)),6))
    
    B$Station_Number = df_no
    
    
    B<- B%>% mutate(across(!Sample_Type, as.numeric))
    
    try(simu_df<- simu_df %>%full_join(B))
    rm(A,B, df_no)
  }
}

write.csv(simu_df, file = "data/simulation_df.csv")



#Adding additional column to fit vpR tool

simu_df<- read.csv("./data/simulation_df.csv")

cols_to_add<- c('c_Bacteria', 'c_HNA', 'c_LNA', 'c_V1', 'c_V2', 'c_V3')

for(col in cols_to_add) {
  simu_df[[col]] <- runif(nrow(simu_df))
}
cols_to_add<- c("Location", "Depth")

for(col in cols_to_add) {
  simu_df[[col]] <- 1
}


#Run simulation

calc_VP(data = simu_df,
        output_dir = "results",
        SR_calc = F,
        bp_endpoint = F)

simu_vp<- read.csv("./results/vp_calc_ALL.csv")
unique(simu_vp$VP_Type)



simu_vp <- simu_vp %>% filter(Population == 'c_Viruses') %>%
  mutate(Station_Number = as.numeric(Station_Number))

write.csv(simu_vp, "./results/simu_vp_filtered.csv")

