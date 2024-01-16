#Bacterial production affects lytic viral production, and also has an impact on
#lysogenic viral production calculations.


#### 1. Linear regression between bacterial counts from bacterial and viral FCM processing ####

#1.1 Extracting viral and bacterial counts ####
project_title<- "Ba_vs_Vi"
work_dir<- paste0(getwd(), "/results/004_bacterial_production_endpoint/", project_title, "/")
source("scripts/004_bacterial_production_endpoint/1_importing_fcs_metadata_source.R")

#Creates necessary directories
set_up_vp_count()

#Processes metadata file and ensures all required variables are present
metadata_processing(paste0(getwd(),"/data/metadata/NJ2020_Ba_vs_Vi.xlsx"), extension = ".xlsx", project_title = "Ba_vs_Vi")

#Directory where raw FCS data is present
source_directory <- paste0(getwd(), "/data/raw_fcs_data/")

#Imports FCS files to working data directory
import_matching_files(source_dir = source_directory, dest_dir = work_dir, metadata = metadata)

#Adding gates
{
  metadata = metadata
  b_ssc = c(1.1, 2.0, 2.5, 3.2, 3.7, 3.7, 2.0, 0.6)
  b_fl1 = c(1.9, 1.4, 1.4,1.8, 2.8, 3.7, 3.2, 2.75)
  v_ssc = c(-0.25, 1.2)
  v_fl1 = c(-0.1, 1.7)
  hna_ssc = c(0.6, 3.5)
  hna_fl1 = c(2.15, 3.5)
  lna_ssc = c(1.0, 3.5)
  lna_fl1 = c(1.4,2.15)  
  v1_ssc = c(-0.1, 0.90)
  v1_fl1 = c(-0.1, 0.8)
  v2_ssc = c(-0.1, 0.90)
  v2_fl1 = c(0.8, 1.25)
  v3_ssc = c(-0.1, 1.3)
  v3_fl1 = c(1.25, 1.7)
}


#Populating gates
populate_gate_df()

#Extracting plots
get_bv_plots()
#only examined the total bacterial gate for this purpose. All look good.

#Extracting Stats
get_bv_stats()


#Import stats csv
counts<- read.csv(paste0(work_dir, "results/", project_title, "_counts.csv"))
str(counts)

#Add metadata
combine_metadata_counts(metadata_df = metadata, counts_df = counts)

#Calculate TE (Control)
calc_TE(df = counts_metadata)

#Adjust counts with TE to get per mL values
adjust_TE(counts_metadata_df = counts_metadata_TE)


#1.2 Perform Linear Regression ####

#Import counts_per_mL file
counts_per_mL<- read.csv(paste0(work_dir, "results/", project_title, "_per_mL.csv"))

ba_vs_vi_df<- counts_per_mL%>%
  select(Staining_Protocol, Station_Number, Sample_Type, c_Bacteria) %>%
  #group_by(Staining_Protocol, Station_Number, Sample_Type) %>%
  pivot_wider(names_from = Staining_Protocol,
              values_from = c_Bacteria)

ggplot(data = ba_vs_vi_df,
       aes(x = Bacteria,
           y = Viruses))+
  geom_point()+
  geom_smooth(method = 'lm')

ba_vs_vi_model<- lm(Viruses ~ Bacteria, data = ba_vs_vi_df)
summary(ba_vs_vi_model)

#not that great of a fit.


#2.0 Extracting NJ2020  viral and bacterial counts ####

#2.1 Viral reduction Assay ####
project_title<- "NJ2020"
work_dir<- paste0(getwd(), "/results/004_bacterial_production_endpoint/", project_title, "/")
source("scripts/004_bacterial_production_endpoint/1_importing_fcs_metadata_source.R")

#Creates necessary directories
set_up_vp_count()

#Processes metadata file and ensures all required variables are present
metadata_processing(paste0(getwd(),"/data/metadata/NJ2020_Metadata.xlsx"), extension = ".xlsx", project_title = "Ba_vs_Vi")

#Directory where raw FCS data is present
source_directory <- paste0(getwd(), "/data/raw_fcs_data/")
#Imports FCS files to working data directory
import_matching_files(source_dir = source_directory, dest_dir = work_dir, metadata = metadata)

#Adding gates
{ #default
  metadata = metadata
  b_ssc = c(1.1, 2.0, 2.5, 3.2, 3.8, 3.8, 2.0, 0.6)
  b_fl1 = c(1.9, 1.4, 1.4,1.8, 2.8, 3.7, 3.2, 2.75)
  v_ssc = c(-0.25, 1.2)
  v_fl1 = c(-0.1, 1.7)
  hna_ssc = c(0.6, 3.8)
  hna_fl1 = c(2.15, 3.5)
  lna_ssc = c(1.0, 3.8)
  lna_fl1 = c(1.4,2.15)  
  v1_ssc = c(-0.1, 0.90)
  v1_fl1 = c(-0.1, 0.8)
  v2_ssc = c(-0.1, 0.90)
  v2_fl1 = c(0.8, 1.25)
  v3_ssc = c(-0.1, 1.3)
  v3_fl1 = c(1.25, 1.7)
}


#Populating gates
rm(gate_df) #from previous step
populate_gate_df()

#Extracting plots
#get_bv_plots()

#These parameters didn't fit all, so some gate adjustments were made.

#Samples 94:177, 181:491 need- V3 lower boundary lowered by 0.05, 
#V2 upper and lower boundary lowered by 0.05, 
#V1 upper boundary lowered by 0.05
populate_gate_df(sample_range = c(94:177, 181:491),
                 g2_v3_fl1 = c(1.20, 1.7),
                 g2_v2_fl1 = c(0.75, 1.20),
                 g2_v1_fl1 = c(-0.1, 0.75)#from default
)

#Samples 1:92,178:180 need - V3 lower boundary lowered by 0.05
#V2 upper boundary lowered by 0.05
populate_gate_df(sample_range = c(1:92,178:180),
                 g2_v3_fl1 = c(1.20, 1.7),
                 g2_v2_fl1 = c(0.8, 1.20)#from default
)

#Samples 492:600 need - V3 lower boundary lowered by 0.1
#V2 upper boundary lowered by 0.1, lower boundary lowered by 0.05
#V1 upper boundary lowered by 0.05
populate_gate_df(sample_range = c(492:600),
                 g2_v3_fl1 = c(1.25, 1.6),
                 g2_v2_fl1 = c(0.75, 1.15),
                 g2_v1_fl1 = c(-0.1, 0.75)#from default
)

#Samples 1:200, 218:226 need - HNA lower boundary raised by 0.15
#LNA upper boundary raised by 0.15
populate_gate_df(sample_range = c(1:200, 218:226),
                 g2_hna_fl1 = c(2.30, 3.5),
                 g2_lna_fl1 = c(1.4,2.30)#from default
)

#Samples 253:491 need - HNA lower boundary raised by 0.1
#LNA upper boundary raised by 0.15
populate_gate_df(sample_range = c(1:200, 218:226),
                 g2_hna_fl1 = c(2.25, 3.5),
                 g2_lna_fl1 = c(1.4,2.25)#from default
)

#Samples 492:600 need - HNA lower boundary lowered by 0.05
#LNA upper boundary lowered by 0.05
populate_gate_df(sample_range = c(1:200, 218:226),
                 g2_hna_fl1 = c(2.10, 3.5),
                 g2_lna_fl1 = c(1.4,2.10)#from default
)

#Extracting plots again after the adjustment
get_bv_plots()



#Extracting Stats
get_bv_stats()


#Import stats csv
counts<- read.csv(paste0(work_dir, "results/", project_title, "_counts.csv"))
str(counts)

#Add metadata
combine_metadata_counts(metadata_df = metadata, counts_df = counts)

#Calculate TE (Control)
calc_TE(df = counts_metadata)

#Adjust counts with TE to get per mL values
adjust_TE(counts_metadata_df = counts_metadata_TE)





#2.2 NJ2020 Bacterial and viral abundance ####
project_title<- "NJ2020_abundance"
work_dir<- paste0(getwd(), "/results/004_bacterial_production_endpoint/", project_title, "/")
source("scripts/004_bacterial_production_endpoint/1_importing_fcs_metadata_source.R")

#Creates necessary directories
set_up_vp_count()

#Processes metadata file and ensures all required variables are present
metadata_processing(paste0(getwd(),"/data/metadata/NJ2020_Abundance_Metadata.xlsx"), extension = ".xlsx", project_title = "Ba_vs_Vi")

#Directory where raw FCS data is present
source_directory <- paste0(getwd(), "/data/raw_fcs_data/")


#Imports FCS files to working data directory
import_matching_files(source_dir = source_directory, dest_dir = work_dir, metadata = metadata)

#Adding gates
{ #default
  metadata = metadata
  b_ssc = c(1.1, 2.0, 2.5, 3.2, 3.8, 3.8, 2.0, 0.6)
  b_fl1 = c(1.9, 1.5, 1.5, 2.0, 2.8, 3.4, 3.2, 2.75)
  v_ssc = c(-0.25, 1.2)
  v_fl1 = c(-0.1, 1.7)
  hna_ssc = c(0.6, 3.8)
  hna_fl1 = c(2.15, 3.5)
  lna_ssc = c(1.0, 3.8)
  lna_fl1 = c(1.5,2.15)  
  v1_ssc = c(-0.1, 0.90)
  v1_fl1 = c(-0.1, 0.75)
  v2_ssc = c(-0.1, 0.90)
  v2_fl1 = c(0.75, 1.2)
  v3_ssc = c(-0.1, 1.125)
  v3_fl1 = c(1.2, 1.65)
}


#Populating gates
rm(gate_df) #from previous step
populate_gate_df()

#Extracting plots
#get_bv_plots()

#These parameters didn't fit all, so some gate adjustments were made.

#Samples 1:2 need-   HNA lower boundary lowered by 0.05
#LNA upper boundary lowered by 0.05
populate_gate_df(sample_range = c(1:2),
                 g2_hna_fl1 = c(2.1, 3.5),
                 g2_lna_fl1 = c(1.5,2.1)#from default
)

#Samples 13:16 need -all lower boundary raised by 0.1
#LNA upper boundary raised by 0.15
populate_gate_df(sample_range = c(13:16),
                 adj_y_all_by = - 0.05#from default
)
#Extracting plots again after the adjustment
get_bv_plots()



#Extracting Stats
get_bv_stats()


#Import stats csv
counts<- read.csv(paste0(work_dir, "results/", project_title, "_counts.csv"))
str(counts)

#Add metadata
combine_metadata_counts(metadata_df = metadata, counts_df = counts)

#Calculate TE (Control)
calc_TE(df = counts_metadata)

#Adjust counts with TE to get per mL values
adjust_TE(counts_metadata_df = counts_metadata_TE)



#3.0 Relative contact/collision rates ####

#Importing per mL .csv file
cr_df<- read.csv(paste0(work_dir, "results/", project_title, "_per_mL.csv"))
cr_df<- cr_df %>% 
  select(c("Location", "Station_Number", "Sample_Type", "Timepoint",
                           "Replicate", "c_Bacteria", "c_Viruses", "VBR")) %>%
  dplyr::filter(Sample_Type != '0.22')

#Relative contact rate is a function of bacteria x viruses

cr_df$BV<-  cr_df$c_Bacteria * cr_df$c_Viruses

#creating empty contact rate list
cr_list<- list()

#Populating cr_list with relative contact rate
for (loc in unique(cr_df$Location)){
  for (expt in unique(cr_df$Station_Number)){
    for(type in unique(cr_df$Sample_Type)){
      for (rep in unique(cr_df$Replicate)){
        
        df10<- cr_df %>%
          dplyr::filter(Location == loc,
                               Station_Number == expt,
                               Sample_Type == type,
                               Replicate == rep)
        
        for(time in unique(cr_df$Timepoint)){
          
          cr<-  (df10 %>% dplyr::filter(Timepoint == time))$BV/(df10 %>% dplyr::filter(Timepoint == 0))$BV
          print(cr)
          cr_value<- c(loc, expt, type, rep, time, cr)
          
          cr_list[[length(cr_list)+1]]<- cr_value
        }
        
      }
    }
  }
}

cr_opt<- data.table(data.table::transpose(as.data.table(cr_list)))
colnames(cr_opt)<- c("Location", "Station_Number", "Sample_Type","Replicate", "Timepoint",
                     "rate")

cr_opt$rate <- as.numeric(cr_opt$rate)
cr_opt$Timepoint <- as.numeric(cr_opt$Timepoint)

#Calculating mean relative contact rates
cr_opt <- cr_opt %>% 
  group_by(Location, Station_Number, Sample_Type, Timepoint) %>%
  summarise(rate_mean = mean(rate), rate_se = plotrix::std.error(rate)) %>%
  drop_na()

#Output .csv
write.csv(cr_opt, 
          file = paste0(work_dir, "/results/", project_title, "_relative_contact_rates.csv"),
          row.names = F)

#3.0 Viral production calculations for NJ2020 ####
library(devtools)

#Install viral production calculator from Github
#install_github("mdhishamshaikh/ViralProduction_R")
library(viralprod)



#Import FCM count csv
counts_per_mL<- read.csv("./results/004_bacterial_production_endpoint/NJ2020/results/NJ2020_per_mL.csv")
str(counts_per_mL)
counts_per_mL<- counts_per_mL %>%
  dplyr::filter(counts_per_mL$Sample_Type != '0.22')
#Importing NJ2020 abundance

abundance<- read.csv(paste0("./results/004_bacterial_production_endpoint/NJ2020_abundance/results/NJ2020_abundance_per_mL.csv"))
str(abundance)
abundance <- abundance %>%
  mutate(Station_Number = as.integer(Station_Number)) %>%
  rename(Total_Bacteria = c_Bacteria,
         Total_Viruses = c_Viruses)


#Some checks
vp_check_populations(counts_per_mL)
vp_class_ori_abu(abundance) 


class(abundance)
class(vp_class_ori_abu(abundance))

vp_class_count_data(counts_per_mL) #passed


#Running viralprod

vp_end_to_end(data = counts_per_mL,
              original_abundances = abundance,
              write_output = T,
              output_dir = "./results/004_bacterial_production_endpoint/NJ2020/results/viralprod_results")

#Step 3 visualize did not work


                     