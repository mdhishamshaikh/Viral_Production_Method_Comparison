
source("./scripts/viral_production_Step2.R")
source("./scripts/3_visualization_source.R")


vis_df<- read.csv("results/simu_vp_filtered.csv")

#### 1.0 Comparison of the two original methods


vis_df <- vis_df  %>%
  mutate(Station_Number = as.numeric(Station_Number))
str(vis_df)
unique(vis_df$Population)
unique(vis_df$VP_Type)


#1.0 Comparison of original two methods ####

#1.1 Set up graphic ####
#1.2 Mock plots to compare methodological difference between LM_SR_AVG (LM_3) and VPCL_AR_DIFF (VPCL_4) ####
#


library(tidyverse)
library(ggsci)
{
  
  df<- data.frame(Timepoint = rep(c('T1', 'T2', 'T3',
                                    'T4', 'T5', 'T6'),6),
                  Sample_Type = c(rep('VP', 18),
                                  rep('VPC', 18)),
                  Replicate = c(rep(1, 6), rep(2, 6), rep(3, 6),
                                rep(4, 6), rep(5, 6), rep(6, 6)),
                  Count = c(2.0, 0.7, 5.0, 6.9, 8.7, 5.0,
                            2.1, 0.3, 5.3, 6.3, 8.3, 5.3,
                            2.9, 1.9, 5.8, 7.2, 7.9, 5.1,
                            1.2, 2.1, 2.0, 6.2, 8.4, 3.0,
                            1.3, 2.4, 1.9, 6.4, 8.2, 3.4,
                            0.8, 2.7, 1.6, 6.1, 8.8, 3.5)
  )
  
  
  df <- df %>% group_by(Sample_Type, Replicate, Timepoint)
  df$Sample_Type<- factor(df$Sample_Type, levels = c('VP', 'VPC', 'Diff'))
  #Mean and SE
  df2<- df%>% group_by(Sample_Type, Timepoint) %>%
    summarise(mean = mean(Count), se = plotrix::std.error(Count)) %>%
    as.data.frame()%>%
    mutate(Timepoint= as.factor(Timepoint))%>%
    mutate(Sample_Type= as.factor(Sample_Type))
  
  df3_a<- df2 %>%
    select(-se)%>%
    pivot_wider(names_from = Sample_Type,
                values_from = c(mean))%>%
    mutate(Diff = VPC - VP)%>%
    pivot_longer(cols = -Timepoint,
                 names_to = 'Sample_Type',
                 values_to = 'mean')
  
  
  df3_b<- df2 %>%
    select(-mean)%>%
    pivot_wider(names_from = Sample_Type,
                values_from = c(se))%>%
    mutate(Diff = VPC - VP)%>%
    pivot_longer(cols = -Timepoint,
                 names_to = 'Sample_Type',
                 values_to = 'se')
  df3<- merge(df3_a, df3_b, by.x = c('Timepoint', 'Sample_Type'),
              by.y = c('Timepoint', 'Sample_Type'))
  
  df3$Sample_Type<- factor(df3$Sample_Type, levels = c('VP', 'VPC', 'Diff'))
  
  
}


cols<- c("VP" = '#e8384fff',
         "VPC" = '#4178bcff',
         "Diff" = '#ebb81cff')
shape<- c("VP" = 15,
          "VPC" = 17,
          "Diff" = 19)

#linear background plot with VP and VPC
bp<- ggplot(data = df, aes(x = as.factor(Timepoint), y = Count, color = Sample_Type, fill = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = 0.5,
             size =2.0,
             position = position_jitter(width = 0.1))+
  scale_shape_manual(values = c(15, 17))+
  theme_classic()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  scale_color_manual(values = cols)+
  theme(legend.position = 'none',
        axis.text = element_text(size = 20, 
                                 face = 'bold'),
        axis.title = element_blank())+
  ylim(0,10)

bp




#####MAIN PLOTS
theme_mockplots<- theme(legend.position = 'none',
                        panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
                        panel.background = element_rect(fill = NA),
                        axis.text = element_text(size = 10, 
                                                 face = 'bold'),
                        axis.title = element_blank())

#lm - vp

lm_1<- ggplot(data = df %>% filter(Sample_Type == 'VP'), aes(x = as.factor(Timepoint), y = Count, color = Sample_Type, fill = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = .8,
             size =2.0,
             position = position_jitter(width = 0.1))+
  scale_shape_manual(values = c(15, 17))+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  scale_color_manual(values = cols)+
  theme_mockplots+
  ylim(-4,10)+
  geom_smooth(aes(group = Replicate), method = 'lm', se = F, alpha = 0.1, size = 0.5)

#LM - VPC
lm_2<- ggplot(data = df %>% filter(Sample_Type == 'VPC'), aes(x = as.factor(Timepoint), y = Count, color = Sample_Type, fill = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = 0.8,
             size =2.0,
             position = position_jitter(width = 0.1))+
  scale_shape_manual(values = c(15, 17))+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  scale_color_manual(values = cols)+
  theme_mockplots+
  ylim(-4,10)+
  geom_smooth(aes(group = Replicate), method = 'lm', se = F, alpha = 0.1, size = 0.5)

#LM -VP VPC
lm_3<- ggplot(data = df , aes(x = as.factor(Timepoint), y = Count, color = Sample_Type, fill = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = 0.8,
             size =2.0,
             position = position_jitter(width = 0.1))+
  scale_shape_manual(values = c(15, 17))+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  scale_color_manual(values = cols)+
  theme_mockplots+
  ylim(-4,10)+
  geom_smooth(aes(group = Replicate), method = 'lm', se = F, alpha = 0.1, size = 0.5)



#vpcl - vp

vpcl_1<-ggplot(data = df %>% filter(Sample_Type == 'VP'), aes(x = as.factor(Timepoint), y = Count, color = Sample_Type, fill = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = 0.2,
             size =2.0,
             position = position_jitter(width = 0.1))+
  scale_shape_manual(values = shape)+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  scale_color_manual(values = cols)+
  theme_mockplots+
  ylim(-4,10)+
  geom_point(data = df3 %>% filter(Sample_Type == 'VP') ,aes(x = as.factor(Timepoint), y = mean, 
                                                             color = Sample_Type, fill = Sample_Type, 
                                                             shape = Sample_Type), alpha = 1, size =2.0) +
  geom_line(data = df3 %>% filter(Sample_Type == 'VP'),
            aes(y = mean, group = Sample_Type),alpha = 1, size = 0.5)

#vpcl -vpc

vpcl_2<-ggplot(data = df %>% filter(Sample_Type == 'VPC'), aes(x = as.factor(Timepoint), y = Count, color = Sample_Type, fill = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = 0.2,
             size = 2.0,
             position = position_jitter(width = 0.1))+
  scale_shape_manual(values = shape)+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  scale_color_manual(values = cols)+
  theme_mockplots+
  ylim(-4,10)+
  geom_point(data = df3 %>% filter(Sample_Type == 'VPC') ,aes(x = as.factor(Timepoint), y = mean, 
                                                              color = Sample_Type, fill = Sample_Type, 
                                                              shape = Sample_Type), alpha = 1, size = 2.0) +
  geom_line(data = df3 %>% filter(Sample_Type == 'VPC'),
            aes(y = mean, group = Sample_Type),alpha = 1, size = 0.5)

#vpcl - vp vpc diff

vpcl_3<- ggplot(data = df, aes(x = as.factor(Timepoint), y = Count, color = Sample_Type, fill = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = 0.2,
             size = 2.0,
             position = position_jitter(width = 0.1))+
  scale_shape_manual(values = shape)+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  scale_color_manual(values = cols)+
  theme_mockplots+
  ylim(-4,10)+
  geom_point(data = df3 ,aes(x = as.factor(Timepoint), y = mean, 
                             color = Sample_Type, fill = Sample_Type, 
                             shape = Sample_Type), alpha = 1, size = 2.0) +
  geom_line(data = df3 %>% filter(Sample_Type != 'Diff'),
            aes(y = mean, group = Sample_Type),alpha = 0.2, size = 0.5)+
  geom_line(data = df3 %>% filter(Sample_Type == 'Diff'),
            aes(y = mean, group = Sample_Type),alpha = 1, size = 0.5)


plot_grid(lm_1, lm_2, lm_3,
          vpcl_1, vpcl_2, vpcl_3,
          nrow = 2) #PLOT: LM VS VIPCAL MOCKPLOTS




og_plots_list<- list(lm_1, vpcl_1, lm_2, vpcl_2, lm_3, vpcl_3)
og_plots_list_name<- c('LM_1', 'VPCL_1',
                       'LM_2', 'VPCL_2',
                       'LM_3', 'VPCL_3')

for (i in 1:6){
  print(i)
  
  ggsave(filename = paste0(og_plots_list_name[i], ".svg"),
         plot = og_plots_list[[i]],
         device = 'svg',
         width = 4,
         height = 4,
         dpi = 800,
         units = 'in')
  
}










#
 




#Linear Regression
violin_df<- simu_vp%>% filter(VP_Type == 'LM_SR_AVG'|  VP_Type == 'VPCL_AR_Diff_No_SE')
violin_df$VP_Type2<- 'A'
ggplot(df6,
       aes(y = VP,
           x = Sample_Type,
           fill = VP_Type))+
  #geom_point()+
  geom_split_violin(nudge = 0.02,
                    color = 'transparent')+
  theme_bw()+
  scale_fill_manual(values = c("black", "#ff9c00ff"))+
  scale_color_manual(values<- NA)+
  theme(axis.line = element_line(linewidth = 1),
        axis.text = element_text(face = 'bold'),
        axis.title = element_text(face = 'bold'),
        legend.box.background = element_rect(color = 'black'),
        legend.box.margin = margin(t = 1, l = 1),
        legend.title = element_text(face = 'bold'),
        legend.position = c(0.9, 0.9))+
  labs(x = 'Treatment',
       y = 'Viral Production')

#1.3.1 Linear Regression

vis_df_og <- vis_df %>%
  filter(VP_Type %in% c('LM_SR_AVG', 'VPCL_AR_DIFF'))

unique(vis_df_og$VP_Type)       

#Check for assumptions
#Normality test
# Levene's Test 
car::leveneTest(VP ~ as.factor(VP_Type), data = vis_df_og)
#p_value = < 2.2e-16. Variance between the two groups is not equal.

# Bartlett's Test
bartlett.test(VP ~ as.factor(VP_Type), data = vis_df_og)
#p_value = < 2.2e-16. Variance between the two groups is not equal.

# Shapiro-Wilk Test (for normality)
shapiro.test(sample(vis_df_og$VP, 5000, replace = FALSE))
#p_value = < 2.2e-16. Not normal.

# Q-Q- plot
qqnorm(vis_df_og$VP)
qqline(vis_df_og$VP)

ggplot(vis_df_og, aes(sample = VP)) + 
  geom_qq() +
  geom_qq_line() +
  theme_bw()

#Create box plots that show distribution of weight loss for each group
boxplot(VP ~ as.factor(VP_Type), xlab='VP_Type', ylab='VP', data=vis_df_og)

#wilcox test #non-normal test
wilcox.test(VP ~ as.factor(VP_Type), data = vis_df_og)
#W = 34050432, p-value < 2.2e-16 Significantly different.

#Linear Regression

vis_df_og_lm<- vis_df_og %>%
  select(-c(abs_VP, VP_SE,
            VP_R_Squared)) %>%
  group_by(Location, Station_Number,
           Depth, Time_Range,
           Population, Sample_Type) %>%
  pivot_wider(names_from = VP_Type,
              values_from = VP) 

#Stats

#Total length of observations
length(vis_df_og_lm$LM_SR_AVG)
#15000

#Cases when LM and VPCL are zero
length(which(vis_df_og_lm$LM_SR_AVG == 0))
#0 0%
length(which( vis_df_og_lm$VPCL_AR_DIFF == 0))
#2056 13.70667%

#Cases when LM and VPCL are less than zero
length(which(vis_df_og_lm$LM_SR_AVG < 0))
#7416 49.44%
length(which( vis_df_og_lm$VPCL_AR_DIFF < 0))
#0 0%


#Cases where LM is higher than VPCL
length(which(vis_df_og_lm$LM_SR_AVG > vis_df_og_lm$VPCL_AR_DIFF))
#511 3.406667%

#Cases when VPCL is higher than LM
length(which(vis_df_og_lm$LM_SR_AVG < vis_df_og_lm$VPCL_AR_DIFF))
#13243 88.28667%

#Cases when VPCL is equal to LM
length(which(vis_df_og_lm$LM_SR_AVG ==vis_df_og_lm$VPCL_AR_DIFF))
#1246 8.306667%



ggplot(data = vis_df_og_lm,
       aes(x = LM_SR_AVG,
           y = VPCL_AR_DIFF,
           col = Sample_Type))+
  geom_point(alpha = 0.1)+
  scale_color_manual(values = cols)+
  theme_bw()+
  geom_abline(slope = 1)+
  geom_smooth(method = 'lm')
  

