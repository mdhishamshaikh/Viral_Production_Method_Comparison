library(tidyverse)
library(readxl)
library(ggsci)
library(ggpubr)
library(cowplot)
library(UpSetR)

ls_df <- read_excel("data/VP_Assay_WoS_Literature_Survey.xlsx", sheet = "ls_data")
ls_df<- ls_df %>% 
  select(Year_Published,
         Lytic, Lysogenic,
         LM, VIPCAL) %>%
  arrange(Year_Published)

#### 1.0 Assays published ####

lm_df_assays<- ls_df %>% 
  dplyr::filter(LM == 1) %>%
  group_by(Year_Published) %>%
  summarise(Lytic_LM = sum(Lytic), Lysogenic_LM = sum(Lysogenic))

vpcl_df_assays<- ls_df %>% 
  dplyr::filter(VIPCAL == 1) %>%
  group_by(Year_Published) %>%
  summarise(Lytic_VIPCAL = sum(Lytic), Lysogenic_VIPCAL = sum(Lysogenic))

cumu_df_assays <- full_join(lm_df_assays ,vpcl_df_assays)%>%
  arrange(Year_Published)%>%
  mutate(across(everything(), ~replace_na(.x,0)))%>%
  mutate(across(!Year_Published, cumsum))%>%
  pivot_longer(cols = 2:5,
               names_to = 'Type',
               values_to = 'Assay_No')%>%
  mutate(
    VP_Type = str_extract(Type, "Lytic|Lysogenic"),
    Analyses_Type = str_extract(Type, "LM|VIPCAL")
  )

ggplot(data = cumu_df_assays,
       aes(x = Year_Published,
           y = Assay_No,
           shape = VP_Type,
           col = Analyses_Type))+
  scale_color_manual(values = c('orange', 'black'))+
  geom_point()+
  geom_line()+
  theme_classic()


#### 2.0 Studies published ####

lm_df_lytic_studies<- ls_df %>% 
  dplyr::filter(LM == 1, Lytic != 0) %>%
  group_by(Year_Published) %>%
  summarise(Lytic_LM = n())
lm_df_lysogenic_studies<- ls_df %>% 
  dplyr::filter(LM == 1, Lysogenic != 0) %>%
  group_by(Year_Published) %>%
  summarise(Lysogenic_LM = n())

vpcl_df_lytic_assays<- ls_df %>% 
  dplyr::filter(VIPCAL == 1, Lytic != 0) %>%
  group_by(Year_Published) %>%
  summarise(Lytic_VIPCAL = n())
vpcl_df_lysogenic_assays<- ls_df %>% 
  dplyr::filter(VIPCAL == 1, Lysogenic != 0) %>%
  group_by(Year_Published) %>%
  summarise(Lysogenic_VIPCAL = n())


cumu_df_studies <- full_join(lm_df_lytic_studies ,lm_df_lysogenic_studies)%>%
  full_join(vpcl_df_lytic_assays) %>%
  full_join(vpcl_df_lysogenic_assays) %>%
  arrange(Year_Published)%>%
  mutate(across(everything(), ~replace_na(.x,0)))%>%
  mutate(across(!Year_Published, cumsum))%>%
  pivot_longer(cols = 2:5,
               names_to = 'Type',
               values_to = 'Study_No')%>%
  mutate(
    VP_Type = str_extract(Type, "Lytic|Lysogenic"),
    Analyses_Type = str_extract(Type, "LM|VIPCAL")
  )

studies_plot_base <- ggplot(data = cumu_df_studies,
       aes(x = Year_Published,
           y = Study_No,
           shape = VP_Type,
           col = Analyses_Type))+
  scale_color_manual(values = c('orange', 'black'))+
  geom_line()+
  geom_point(size = 3)+
  theme_classic() +
  labs(x = "Year",
       y = "Number of studies",
       color = "Data analyses\nmethods",
       shape = "Viral production\ntypes")+
  theme(axis.title = element_text(size = 16,
                                  face = 'bold'),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14,
                                    face = 'bold'))+ 
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  )
  

studies_plot_base

get_legend(studies_plot_base)

studies_plot<- ggdraw(studies_plot_base +
         theme(legend.position = "none"))+
  draw_plot(get_legend(studies_plot_base), x = 0.2, y = 0.5, width = 0.3, height = 0.4)

ggsave(filename = "ls_studies_plot.svg",
       plot = studies_plot,
       device = "svg",
       width = 5,
       height = 5)



#Studies Venn diagram 

upset_df<- ls_df%>%
mutate_at(vars(Lytic, Lysogenic, LM, VIPCAL), ~ifelse(. > 1, 1, .)) %>%
  select(-c(Year_Published))


#Venn
upset(
  upset_df, 
  sets = c("LM", "VIPCAL", "Lytic", "Lysogenic"), 
  sets.bar.color = "#56B4E9",
  matrix.color = "#D55E00",
  keep.order = TRUE,
  order.by = "freq"
)

#Authors vs data analyses type

#Importing data

ls_df <- read_excel("data/VP_Assay_WoS_Literature_Survey.xlsx", sheet = "ls_data")
ls_df$doi <- tolower(ls_df$doi) 

wos_df<- read_excel("data/VP_Assay_WoS_Literature_Survey.xlsx", sheet = "wos_metadata")
wos_df$DOI...52 <- tolower(wos_df$DOI...52) #Some here might be capitalized

#I will combine the author records from wos_df to ls_df based on doi

comb_df<- full_join(ls_df, wos_df,by = c('doi' = 'DOI...52'))

#Extracting the necessary columns
comb_df <- comb_df %>% 
  select(Year_Published,doi, 'Author Full Names',
         Lytic, Lysogenic,
         LM, VIPCAL) %>%
  mutate_at(vars(Lytic, Lysogenic, LM, VIPCAL), ~ifelse(. > 1, 1, .)) %>% #as we are examining number of studies and not assays
  arrange(Year_Published) %>%
  separate_rows(`Author Full Names`, sep = "; ") %>%
  mutate(Author_Last_Name = word(`Author Full Names`, 1, sep = ", ")) %>%
  select(Year_Published, doi, Author_Last_Name, Lytic, Lysogenic, LM, VIPCAL)


authors_df_sum<- comb_df %>% 
  group_by(Author_Last_Name) %>%
  summarize(
    Lytic = sum(Lytic),
    Lysogenic = sum(Lysogenic),
    LM = sum(LM),
    VIPCAL = sum(VIPCAL)
  )

unique_authors <- unique(comb_df$Author_Last_Name)
author_colors <- rainbow(length(unique_authors))  # This will generate a spectrum of colors. Adjust as needed.
names(author_colors) <- unique_authors

upset(
  comb_df, 
  sets = c("LM", "VIPCAL", "Lytic", "Lysogenic"), 
  sets.bar.color = "#56B4E9",
  matrix.color = "#D55E00",
  keep.order = TRUE,
  order.by = "freq"
) + scale_fill_manual(values = author_colors)

