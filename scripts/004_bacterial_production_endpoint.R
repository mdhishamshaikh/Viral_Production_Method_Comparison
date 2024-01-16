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
metadata_processing(paste0(getwd(),"/data/metadata/NJ_2020_Ba_vs_Vi.xlsx"), extension = ".xlsx", project_title = "Ba_vs_Vi")

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


