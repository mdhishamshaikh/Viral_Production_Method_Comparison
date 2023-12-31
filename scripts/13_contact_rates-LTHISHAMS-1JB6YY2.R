library(tidyverse)
library(data.table)
library(cowplot)
library(colorspace)


cr_df<- cr_df <- read.csv("C:/Users/hisham.shaikh/OneDrive - UGent/Projects/FCM_R/ViralProduction_R/NJ2020.csv")





cr_df<- cr_df %>% select(c("Location", "Expt_No", "Sample_Type", "Timepoint",
                  "Replicate", "c_Bacteria", "c_Viruses", "VBR"))

abundance<- read.csv("C:/Users/hisham.shaikh/OneDrive - UGent/Projects/Microbial_Abundances/Microbial_Abundances_NJ2020_PE477_PE486/Microbial Abundances/Microbial_Abundances/Microbial_Abundances/NJ2020_abundance.csv")

abundance1<- abundance %>%
  select(Location, Expt_No, Total_Bacteria, Total_Viruses, VBR) %>%
  rename(c_Bacteria = Total_Bacteria,
         c_Viruses = Total_Viruses )

loc<- 'NJ2020'
expt_no<- 1:7
type<- c('VP', 'VPC', '0.22')
time<- -3
rep<- 1:3
abundance2<- data.frame(expand.grid(loc, expt_no, type, time, rep))
colnames(abundance2)<- c('Location', 'Expt_No', 'Sample_Type',
                         'Timepoint', 'Replicate')




abundance3<- full_join(abundance2, abundance1)


cr_df<- full_join(abundance3, cr_df)


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
          print(cr)
          cr_value<- c(loc, expt, type, rep, time, cr)
          
          cr_list[[length(cr_list)+1]]<- cr_value
        }
          
      }
    }
  }
}


cr_opt<- data.table(data.table::transpose(as.data.table(cr_list)))
colnames(cr_opt)<- c("Location", "Expt_No", "Sample_Type","Replicate", "Timepoint",
             "rate")
cr_opt$rate <- as.numeric(cr_opt$rate)
#Warnings are because of o.22 samples that have negative amount of bacteria. This could be due to TE corection.
#We either get call the bacterial count as 0 and calculate for 0.22 or exclude it from further analyses.
#For now, we'll exclude it.
cr_opt <- cr_opt %>% filter(Sample_Type != '0.22')


#Let's visualise this data #####

#First we summarise the replicates
cr_opt<- cr_opt %>% filter(Timepoint != -3)

cr_opt <- cr_opt %>% 
  group_by(Location, Expt_No, Sample_Type, Timepoint) %>%
  summarise(rate_mean = mean(rate), rate_se = plotrix::std.error(rate))

cr_opt2<- cr_opt %>% ungroup() %>%
  select(-rate_se) %>%
  pivot_wider(names_from = Sample_Type,
              values_from = rate_mean) %>%
  mutate(cr_diff = VP - VPC)


cr_opt[cr_opt$Expt_No == 1 & cr_opt$Timepoint == 9,]$Timepoint<- '17'
cr_opt[cr_opt$Expt_No == 2 & cr_opt$Timepoint == 9,]$Timepoint<- '17.5'
cr_opt[cr_opt$Expt_No == 1 & cr_opt$Timepoint == 12,]$Timepoint<- '20'
cr_opt[cr_opt$Expt_No == 2 & cr_opt$Timepoint == 12,]$Timepoint<- '20.5'
cr_opt$Timepoint<- as.numeric(cr_opt$Timepoint)

cr_opt2[cr_opt2$Expt_No == 1 & cr_opt2$Timepoint == 9,]$Timepoint<- '17'
cr_opt2[cr_opt2$Expt_No == 2 & cr_opt2$Timepoint == 9,]$Timepoint<- '17.5'
cr_opt2[cr_opt2$Expt_No == 1 & cr_opt2$Timepoint == 12,]$Timepoint<- '20'
cr_opt2[cr_opt2$Expt_No == 2 & cr_opt2$Timepoint == 12,]$Timepoint<- '20.5'
cr_opt2$Timepoint<- as.numeric(cr_opt2$Timepoint)
cols<- c("VP" = '#d32f27ff',
         "VPC" = '#4888a2ff',
         "Diff" = '#edb81dff')
shapes<- c("VP" = 21,
           "VPC" = 24,
           "Diff" = 19)
shapes<- c("VP" = 15,
           "VPC" = 17,
           "Diff" = 19)


# Now let's overlap bacterial production values here
bep_df <- read.csv("C:/Users/hisham.shaikh/OneDrive - UGent/Projects/FCM_R/ViralProduction_R/V5000/vp_calc_bp.csv")

bep_df <- bep_df %>%
  select(Location, Expt_No,  Time_Range) %>%
  unique() %>%
  mutate(endpoint = as.numeric(sub(".*\\T0_T?\\s*(\\d+)", "\\1", Time_Range)) )#function to extract the endpoint

bep_df <- merge(cr_opt2, bep_df)


 
 cr_plot1<- ggplot()+
   geom_vline(data = bep_df,
              aes(xintercept = endpoint),
                  colour = '#7E6148FF',
              linewidth =1.5,
              alpha = 0.5,
              show.legend = F)+
   geom_point(data = cr_opt,
              aes(x = as.numeric(Timepoint),
                  y = rate_mean,
                  fill = Sample_Type,
                  shape = Sample_Type,
                  colour = Sample_Type),
              show.legend =  T)+
   geom_line(data = cr_opt2,
             aes(x = as.numeric(Timepoint),
                 y = cr_diff),
             size = 1,
             alpha = 0.5,
             show.legend = F) +
   theme_bw()+
   facet_wrap(. ~ Expt_No, ncol = 2
              # scales = 'free',
   )+
   scale_fill_manual(values = cols)+
   scale_color_manual(values = cols#, breaks = 'VP', labels = 'Bacterial Endpoint'
                      )+
   scale_shape_manual(values = shapes)+
   theme(strip.background = element_rect(fill = lighten("#323D5E",
                                                        amount = 0.0) ,
                                         color = NA),
         strip.text = element_text(face = 'bold',
                                   color = 'white',
                                   size = 10),
         panel.border = element_rect(linewidth = 2),
         panel.background = element_rect(fill = NA),
         #legend.position = c(0.75, 0.1),
         legend.title = element_text(face = 'bold',
                                     size = 10),
         legend.text = element_text(size = 9),
         axis.title = element_text(face = 'bold',
                                   size = 10),
         axis.text = element_text(size = 10))+
   labs(x = 'Timepoint (in hours)',
        y = 'Mean Relative Collision Rate')  +
   guides(color = guide_legend(title="Treatment"),
          shape = guide_legend(title="Treatment"),
          fill= guide_legend(title="Treatment"))

cr_plot1 
 
legend1<- ggpubr::as_ggplot( ggpubr::get_legend(cr_plot1))

cr_plot<- cr_plot1 +theme(legend.position = 'none',
                          panel.grid = element_blank())

cr_plot2<- ggplot()+
  geom_vline(data = bep_df,
             aes(xintercept = endpoint),
             colour = '#7E6148FF',
             linewidth =1.5,
             alpha = 0.5,
             show.legend = T)+
  geom_point(data = cr_opt,
             aes(x = as.numeric(Timepoint),
                 y = rate_mean,
                 #fill = Sample_Type,
                 #shape = Sample_Type,
                 colour = Sample_Type
                 ),
             show.legend =  F)+
  # geom_line(data = cr_opt2,
  #           aes(x = as.numeric(Timepoint),
  #               y = cr_diff),
  #           size = 1,
  #           alpha = 0.5,
  #           show.legend = F) +
  theme_bw()+
  facet_wrap(. ~ Expt_No, ncol = 2
             # scales = 'free',
  )+
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols, breaks = 'VP', labels = 'Bacterial\nEndpoint'
  )+
  scale_shape_manual(values = shapes)+
  theme(strip.background = element_rect(fill = 'gray',
                                        color = NA),
        strip.text = element_text(face = 'bold'),
        panel.border = element_rect(linewidth = 2),
        panel.background = element_rect(fill = NA),
        #legend.position = c(0.75, 0.1),
        legend.title = element_text(face = 'bold',
                                    size = 10),
        legend.text = element_text(size = 9),
        axis.title = element_text(face = 'bold',
                                  size = 15),
        axis.text = element_text(size = 15))+
  labs(x = 'Timepoint (in hours)',
       y = 'Mean Relative Collision Rate')  +
  guides(color = guide_legend(title = ''),
         #shape = guide_legend(title="Treatment"),
         #fill= guide_legend(title="Treatment")
         )

cr_plot2
 
legend2<- ggpubr::as_ggplot( ggpubr::get_legend(cr_plot2))


cr_plot3<- ggplot()+
  # geom_vline(data = bep_df,
  #            aes(xintercept = endpoint),
  #            colour = '#7E6148FF',
  #            linewidth =1.5,
  #            alpha = 0.5,
  #            show.legend = T)+
  geom_point(data = cr_opt,
             aes(x = as.numeric(Timepoint),
                 y = rate_mean,
                 #fill = Sample_Type,
                 #shape = Sample_Type,
                 colour = Sample_Type
             ),
             show.legend =  F)+
  geom_line(data = cr_opt2,
            aes(x = as.numeric(Timepoint),
                y = cr_diff),
            size = 1,
            alpha = 0.5,
            show.legend = T) +
  theme_bw()+
  facet_wrap(. ~ Expt_No, ncol = 2
             # scales = 'free',
  )+
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols, breaks = 'VP', labels = 'Differences in\ncollision rates\nbetween VP\nand VPC\ntreatments'
  )+
  scale_shape_manual(values = shapes)+
  theme(strip.background = element_rect(fill = lighten("#323D5E",
                                                       amount = 0.1) ,
                                        color = NA),
        strip.text = element_text(face = 'bold',
                                  color = 'white'),
        panel.border = element_rect(linewidth = 2),
        panel.background = element_rect(fill = NA),
        #legend.position = c(0.75, 0.1),
        legend.title = element_text(face = 'bold',
                                    size = 10),
        legend.text = element_text(size = 9),
        axis.title = element_text(face = 'bold',
                                  size = 15),
        axis.text = element_text(size = 15))+
  labs(x = 'Timepoint (in hours)',
       y = 'Mean Relative Collision Rate')  +
  guides(color = guide_legend(title = ''),
         #shape = guide_legend(title="Treatment"),
         #fill= guide_legend(title="Treatment")
  )

cr_plot3

legend3<- ggpubr::as_ggplot( ggpubr::get_legend(cr_plot3))

legends<- cowplot::plot_grid(plotlist = list(legend1, legend2, legend3),
                   ncol= 1,
                   rel_widths = c(1, 1,  1))
legends2<- plot_grid( NULL, legends,
          rel_widths = c(2,3),
          nrow = 2)



cr_plot_final<- cr_plot + NULL+ legends +
  plot_layout(width = c(8, 2,2))
cr_plot_final



cr_plot_final<- plot_grid(cr_plot, legends,
                               nrow = 1,
                               rel_widths = c( 2, 1))


cr_plot_final


ggsave(file = 'collision_rate_ww.png', width = 10,
       height = 18, unit = 'cm')








