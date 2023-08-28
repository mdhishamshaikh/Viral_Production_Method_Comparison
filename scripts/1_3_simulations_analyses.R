library(tidyverse)
library(ggstatsplot)
library(ggsci)
library(cowplot)
install.packages(parallel)
unique(output2$VP_Type)

make_a_csv <- read.csv("C:/Users/hisham.shaikh/OneDrive - UGent/Projects/FCM_R/ViralProduction_R/simulations/simu_output_df.csv")


#### 1. Linear Model Comparison ####
lm_vp_type <- c('LM_AP', 'LM_SR_AVG',
                'LM_AR', 'LM_AR_Diff',
                'LM_AR_Diff_LMER')

lm_output<- make_a_csv %>%
  filter(VP_Type  %in% c('LM_AP', 'LM_SR_AVG',
                         'LM_AR', 'LM_AR_Diff',
                         'LM_AR_Diff_LMER'))
unique(lm_output$VP_Type)

vpcl_output<- make_a_csv %>%
  filter(!(VP_Type  %in% c('LM_AP', 'LM_SR_AVG',
                           'LM_AR', 'LM_AR_Diff',
                           'LM_AR_Diff_LMER')))

unique(vpcl_output$VP_Type)
        
ggplot(data = lm_output,
       aes(x = VP_Type,
           y = VP,
           fill = VP_Type))+
  geom_violin()+
  geom_boxplot()

lm_comp<- ggbetweenstats(
    data = lm_output,
    x = VP_Type,
    y = VP,
    type = "nonparametric", # ANOVA or Kruskal-Wallis
    plot.type = "violin",
    pairwise.comparisons = TRUE,
    pairwise.display = "s",
    centrality.plotting = FALSE,
    bf.message = FALSE
  )+labs(x = 'VP Data Analyses Type',
         y = 'Viral Production Rate',
    title = "Kruskal-Wallis Test: LM VP")


 ggbetweenstats(
   data = lm_output,
   x = VP_Type,
   y = VP,
   violin.args = list(linewidth = 1,
                      fill = NA),
   type = 'nonparametric',
   plot.type = 'violin',
   point.args = list(alpha= 0, position = position_jitterdodge(jitter.width = 1)),
   pairwise.display = 'all',
   centrality.type = 'parametric',
   centrality.plotting = F,
   centrality.label.args = list(alpha = 0),
   centrality.point.args = list(color = 'black',
                                size = 3)
   
 )+
   geom_violin(data = lm_output,
               aes(x = VP_Type,
                   y = VP,
                   fill = VP_Type))+
   geom_boxplot(size = 1,
                outlier.alpha = 0,
                width = 0.1)+
   scale_fill_lancet()
 
 

a<-  ggbetweenstats(
   data = lm_output,
   x = VP_Type,
   y = VP,
   violin.args = list(linewidth = 1,
                      fill = NA),
   type = 'nonparametric',
   plot.type = 'violin',
   point.args = list(alpha= 0, position = position_jitterdodge(jitter.width = 1)),
   pairwise.display = 's',
   centrality.type = 'parametric',
   centrality.plotting = F,
   centrality.label.args = list(alpha = 0),
   centrality.point.args = list(color = 'black',
                                size = 3)
   
 )+
   geom_violin(data = lm_output,
               aes(x = VP_Type,
                   y = VP,
                   fill = VP_Type))+
   geom_boxplot(size = 1,
                outlier.alpha = 0,
                width = 0.1)+
   scale_fill_npg()+
   theme_classic()+
   theme(legend.position = 'none')
 
 
 b<- ggbetweenstats(
   data = lm_output,
   x = VP_Type,
   y = VP_SE,
   violin.args = list(linewidth = 1,
                      fill = NA),
   type = 'nonparametric',
   plot.type = 'violin',
   point.args = list(alpha= 0, position = position_jitterdodge(jitter.width = 1)),
   pairwise.display = 's',
   centrality.type = 'parametric',
   centrality.plotting = F,
   centrality.label.args = list(alpha = 0),
   centrality.point.args = list(color = 'black',
                                size = 3)
   
 )+
   geom_violin(data = lm_output,
               aes(x = VP_Type,
                   y = VP_SE,
                   fill = VP_Type))+
   geom_boxplot(size = 1,
                outlier.alpha = 0,
                width = 0.1)+
   scale_fill_npg()+
   theme_classic()+
   theme(legend.position = 'none')
 
 extract_stats(a)
 extract_stats(b)
 extract_subtitle(a)
 
 
 combine_plots(list(a,b))
 
 
 
 c<-  ggbetweenstats(
   data = vpcl_output,
   x = VP_Type,
   y = VP,
   violin.args = list(linewidth = 1,
                      fill = NA),
   type = 'nonparametric',
   plot.type = 'violin',
   point.args = list(alpha= 0, position = position_jitterdodge(jitter.width = 1)),
   pairwise.display = 's',
   centrality.type = 'parametric',
   centrality.plotting = F,
   centrality.label.args = list(alpha = 0),
   centrality.point.args = list(color = 'black',
                                size = 3)
   
 )+
   geom_violin(data = vpcl_output,
               aes(x = VP_Type,
                   y = VP,
                   fill = VP_Type))+
   geom_boxplot(size = 1,
                outlier.alpha = 0,
                width = 0.1)+
   scale_fill_npg()+
   theme_classic()+
   theme(legend.position = 'none')
 
 
 d<- ggbetweenstats(
   data = vpcl_output,
   x = VP_Type,
   y = VP_SE,
   violin.args = list(linewidth = 1,
                      fill = NA),
   type = 'nonparametric',
   plot.type = 'violin',
   point.args = list(alpha= 1, position = position_jitterdodge(jitter.width = 1)),
   pairwise.display = 's',
   centrality.type = 'parametric',
   centrality.plotting = F,
   centrality.label.args = list(alpha = 0),
   centrality.point.args = list(color = 'black',
                                size = 3)
   
 )+
   geom_violin(data = vpcl_output,
               aes(x = VP_Type,
                   y = VP_SE,
                   fill = VP_Type))+
   geom_boxplot(size = 1,
                outlier.alpha = 0,
                width = 0.1)+
   scale_fill_npg()+
   theme_classic()+
   theme(legend.position = 'none')
 
 extract_stats(c)
 extract_stats(d)
 extract_subtitle(c)
 
 
 combine_plots(list(c,d))
 
 
 
 e<-  ggbetweenstats(
   data = vpcl_output %>%
     filter(VP_Type %in% c('VPCL_AR_Diff_No_SE',
                           'VPCL_AR_Diff_LMER_SE')),
   x = VP_Type,
   y = VP,
   violin.args = list(linewidth = 1,
                      fill = NA),
   type = 'nonparametric',
   plot.type = 'violin',
   point.args = list(alpha= 0, position = position_jitterdodge(jitter.width = 1)),
   pairwise.display = 'ns',
   centrality.type = 'parametric',
   centrality.plotting = F,
   centrality.label.args = list(alpha = 0),
   centrality.point.args = list(color = 'black',
                                size = 3)
   
 )+
   geom_violin(data =  vpcl_output %>%
                 filter(VP_Type %in% c('VPCL_AR_Diff_No_SE',
                                       'VPCL_AR_Diff_LMER_SE')),
               aes(x = VP_Type,
                   y = VP,
                   fill = VP_Type))+
   geom_boxplot(size = 1,
                outlier.alpha = 0,
                width = 0.05)+
   scale_fill_npg()+
   theme_classic()+
   theme(legend.position = 'none')
 
 
 f<- ggbetweenstats(
   data =  vpcl_output %>%
     filter(VP_Type %in% c('VPCL_AR_Diff_No_SE',
                           'VPCL_AR_Diff_LMER_SE')),
   x = VP_Type,
   y = VP_SE,
   violin.args = list(linewidth = 1,
                      fill = NA),
   type = 'nonparametric',
   plot.type = 'violin',
   point.args = list(alpha= 1, position = position_jitterdodge(jitter.width = 1)),
   pairwise.display = 'ns',
   centrality.type = 'parametric',
   centrality.plotting = F,
   centrality.label.args = list(alpha = 0),
   centrality.point.args = list(color = 'black',
                                size = 3)
   
 )+
   geom_violin(data =  vpcl_output %>%
                 filter(VP_Type %in% c('VPCL_AR_Diff_No_SE',
                                       'VPCL_AR_Diff_LMER_SE')),
               aes(x = VP_Type,
                   y = VP_SE,
                   fill = VP_Type))+
   geom_boxplot(size = 1,
                outlier.alpha = 0,
                width = 0.05)+
   scale_fill_npg()+
   theme_classic()+
   theme(legend.position = 'none')
 
 extract_stats(e)
 extract_stats(f)
 extract_subtitle(e)
 
 
 combine_plots(list(e,f))
 
 
 
 ggbetweenstats(data =  vpcl_output %>%
                  filter(VP_Type %in% c('VPCL_AR_Diff_No_SE',
                                        'VPCL_AR_Diff_LMER_SE')),
                x= VP_Type,
                y = VP,
                pairwise.display = 'ns')
 
 simu_output_df <- read.csv("C:/Users/hisham.shaikh/OneDrive - UGent/Projects/FCM_R/ViralProduction_R/simulations/simu_output_df.csv")
 
 simu_output_df<- simu_output_df %>%
   filter(VP_Type %in% c('VPCL_AR_Diff_No_SE',
                         'VPCL_AR_Diff_LMER_SE',
                         'LM_AP'))
simu_output_df$VP_Type<-   factor(simu_output_df$VP_Type, levels = c('VPCL_AR_Diff_No_SE',
                     'VPCL_AR_Diff_LMER_SE',
                     'LM_AP'))
 
 ggbetweenstats(data =  simu_output_df ,
                x= VP_Type,
                y = VP,
                pairwise.display = 's')
 
 
 
 ggbetweenstats(data =  vpcl_output %>%
                  filter(VP_Type %in% c('VPCL_AR_Diff_No_SE',
                                        'VPCL_AR_Diff_LMER_SE')),
                x= VP_Type,
                y = VP,
                PA)
 
 
 
 
 geom_boxplot()+
   scale_fill_manual(values = NULL)
 
 ggplot(data = lm_output,
           aes(x = VP_Type,
               y = VP,
               color = Sample_Type,
               shape = Sample_Type))+
   ggbeeswarm::geom_beeswarm(cex = 0.3,
                             alpha = 0.2)
   
 
 geom_point(stat = 'identity',
              position = position_jitterdodge(),
              alpha = 0.1)+
   geom_violin(aes(fill = NULL),
               draw_quantiles = c(0.25, 0.5, 0.75))
 
 
 df6<- simu_output_df %>% 
   select(-c(VP_SE, VP_R_Squared, n))%>%group_by(Expt_No, Sample_Type)%>%
   filter(VP_Type %in% c('VPCL_AR_Diff_No_SE',
                         'VPCL_AR_Diff_LMER_SE'))
 df6$Sample_Type<- factor(df6$Sample_Type,
                             levels = c("VP",
                                        "VPC",
                                        "Diff"))
 df7<- simu_output_df %>% 
   select(-c(VP_SE, VP_R_Squared, n))%>%group_by(Expt_No, Sample_Type)%>%
   filter(VP_Type %in% c('VPCL_AR_Diff_No_SE',
                         'VPCL_AR_Diff_LMER_SE',
                         'LM_AP')) %>% pivot_wider(names_from = VP_Type, values_from = VP)
 
 
 lm_vp_fit<- summary(lm(VPCL_AR_Diff_No_SE ~ VPCL_AR_Diff_LMER_SE, data = df7 %>% filter(Sample_Type == 'VP')))
 lm_vp_fit$coefficients[[2]]
 
 lm_vpc_fit<- summary(lm(VPCL_AR_Diff_No_SE ~ VPCL_AR_Diff_LMER_SE, data = df7 %>% filter(Sample_Type == 'VPC')))
 lm_vpc_fit$coefficients[[2]]
 
 lm_diff_fit<- summary(lm(VPCL_AR_Diff_No_SE ~ VPCL_AR_Diff_LMER_SE, data = df7 %>% filter(Sample_Type == 'Diff')))
 lm_diff_fit$coefficients
 
 lm_fit<- summary(lm(VPCL_AR_Diff_No_SE ~ VPCL_AR_Diff_LMER_SE, data = df7))
 lm_fit$coefficients[[2]]
 
 
 
 
 vp_se_comp_legend<- ggplot(data = df7, aes(x = VPCL_AR_Diff_No_SE,
                                     y = VPCL_AR_Diff_LMER_SE,
                                     color = Sample_Type,
                                     fill = Sample_Type,
                                     shape = Sample_Type))+
   geom_point(alpha = 0.2)+
   theme_classic()+
   geom_abline(slope = 1, linewidth = 1)+
   scale_color_npg()+
   labs(x = 'VIPCAL',
        y = 'VIPCAL-SE')+
  geom_smooth(method = 'lm', linewidth =2, alpha = 0.2)+
   scale_color_manual(values = cols)+
   scale_fill_manual(values = cols)+
   scale_shape_manual(values = shape)+
   theme(axis.line = element_line(linewidth = 1),
         axis.text = element_text(face = 'bold', size = 15),
         axis.title = element_text(face = 'bold', size = 15),
         legend.box.background = element_rect(color = 'black'),
         legend.box.margin = margin(t = 1, l = 1),
         legend.title = element_text(face = 'bold'))+
   # xlim(c(-0.1, 2.2))+
   # ylim(c(-0.1, 2.2))+
   labs(fill = 'Treatment',
        color = 'Treatment',
        shape = 'Treatment')
 
vp_se_comp<-  vp_se_comp_legend +
   theme(legend.position = 'none')
 
 

 
 
 ggdraw()+
   draw_plot(vp_se_comp)+
   draw_plot(ggpubr::as_ggplot(ggpubr::get_legend(vp_se_comp_legend)),
             x = 0.35, y = -0.15)
 
ggpubr::as_ggplot(ggpubr::get_legend(vp_se_comp))

#split violins 
violin_df<- simu_output_df%>% filter(VP_Type == 'LM_SR_AVG'|  VP_Type == 'VPCL_AR_Diff_No_SE')
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

 
  