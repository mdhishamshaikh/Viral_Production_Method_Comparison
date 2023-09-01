library(tidyverse)


cr_df<- read.csv("NJ2020.csv")
cr_df<- cr_df %>% select(c("Location", "Expt_No", "Sample_Type", "Timepoint",
                  "Replicate", "c_Bacteria", "c_Viruses", "VBR"))
cr_df$BV<-  cr_df$c_Bacteria * cr_df$c_Viruses

cr_list<- list()

for (loc in unique(cr_df$Location)){
  for (expt in unique(cr_df$Expt_No)){
    for(type in unique(cr_df$Sample_Type)){
      for (rep in unique(cr_df$Replicate)){
        
        df10<- cr_df %>%
          filter(Location == loc,
                 Expt_No == expt,
                 Sample_Type == type,
                 Replicate == rep)
        
        for(time in unique(cr_df$Timepoint)){
          
          cr<-  (df10 %>% filter(Timepoint == time))$BV/(df10 %>% filter(Timepoint == 0))$BV
          
          cr_value<- c(loc, expt, type, rep, time, cr)
          
          cr_list[[length(cr_list)+1]]<- cr_value
        }
          
      }
    }
  }
}

cr_opt<- data.frame(t(sapply(cr_list, c)))

dplyr::bind_rows(cr_list)

