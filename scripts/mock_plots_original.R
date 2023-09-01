library(tidyverse)

#### 1.0 Make mock plots to explain LM and VIPCAL methods.

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




#####MAIN PLOTS####
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
          nrow = 2) #PLOT: LM VS VIPCAL MOCKPLOTS####



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


# Add the Y-axis titles
titles <- ggdraw() + 
  draw_label("Linear Regression", x = 0.1, y = 0.8, hjust = 0, vjust = 0.5, fontface = 'bold') + 
  draw_label("VIPCAL", x = 0.1, y = 0.3, hjust = 0, vjust = 0.5, fontface = 'bold') + 
  theme(plot.margin = margin(0, 0, 0, 7 * 1.2),  # Adjust the margin to make it more centered
        plot.background = element_rect(fill = "white"))

# For top titles (VP, VPC, and Diff)
vp_label <- ggdraw() + 
  draw_label("VP", x = 0.2, y = 0.5, fontface = 'bold') + 
  theme(plot.margin = margin(0, 3.8 * 1.2, 0, 0),  # Adjust the margin to make it centered
        plot.background = element_rect(fill = "white"))

vpc_label <- ggdraw() + 
  draw_label("VPC", x = 0.5, y = 0.5, fontface = 'bold') + 
  theme(plot.background = element_rect(fill = "white"))

diff_label <- ggdraw() + 
  draw_label("Diff", x = 0.8, y = 0.5, fontface = 'bold') + 
  theme(plot.margin = margin(0, 0, 0, 3.8 * 1.2),  # Adjust the margin to make it centered
        plot.background = element_rect(fill = "white"))

# Combine everything
plot_grid(vp_label, lm_1, vpcl_1, 
          ncol = 1,
          axis = "tblr")

full_plot <- plot_grid(NULL, titles, lm_1, lm_2, lm_3, vp_label, vpc_label, diff_label,  vpcl_1, vpcl_2, vpcl_3, 
                       ncol = 3, align = 'v', axis = 'l', rel_widths = c(1, 4, 4, 4))

print(full_plot)

