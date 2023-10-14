
source("./scripts/viral_production_Step2.R")
source("./scripts/3_visualization_source.R")
####Geom Split Violin function####
{
  GeomSplitViolin <- ggplot2::ggproto(
    "GeomSplitViolin",
    ggplot2::GeomViolin,
    draw_group = function(self,
                          data,
                          ...,
                          # add the nudge here
                          nudge = 0,
                          draw_quantiles = NULL) {
      data <- transform(data,
                        xminv = x - violinwidth * (x - xmin),
                        xmaxv = x + violinwidth * (xmax - x))
      grp <- data[1, "group"]
      newdata <- plyr::arrange(transform(data,
                                         x = if (grp %% 2 == 1) xminv else xmaxv),
                               if (grp %% 2 == 1) y else -y)
      newdata <- rbind(newdata[1, ],
                       newdata,
                       newdata[nrow(newdata), ],
                       newdata[1, ])
      newdata[c(1, nrow(newdata)-1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
      
      # now nudge them apart
      newdata$x <- ifelse(newdata$group %% 2 == 1,
                          newdata$x - nudge,
                          newdata$x + nudge)
      
      if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
        
        stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
        
        quantiles <- ggplot2:::create_quantile_segment_frame(data,
                                                             draw_quantiles)
        aesthetics <- data[rep(1, nrow(quantiles)),
                           setdiff(names(data), c("x", "y")),
                           drop = FALSE]
        aesthetics$alpha <- rep(1, nrow(quantiles))
        both <- cbind(quantiles, aesthetics)
        quantile_grob <- ggplot2::GeomPath$draw_panel(both, ...)
        ggplot2:::ggname("geom_split_violin",
                         grid::grobTree(ggplot2::GeomPolygon$draw_panel(newdata, ...),
                                        quantile_grob))
      }
      else {
        ggplot2:::ggname("geom_split_violin",
                         ggplot2::GeomPolygon$draw_panel(newdata, ...))
      }
    }
  )
  geom_split_violin <- function(mapping = NULL,
                                data = NULL,
                                stat = "ydensity",
                                position = "identity",
                                # nudge param here
                                nudge = 0,
                                ...,
                                draw_quantiles = NULL,
                                trim = TRUE,
                                scale = "area",
                                na.rm = FALSE,
                                show.legend = NA,
                                inherit.aes = TRUE) {
    
    ggplot2::layer(data = data,
                   mapping = mapping,
                   stat = stat,
                   geom = GeomSplitViolin,
                   position = position,
                   show.legend = show.legend,
                   inherit.aes = inherit.aes,
                   params = list(trim = trim,
                                 scale = scale,
                                 # don't forget the nudge
                                 nudge = nudge,
                                 draw_quantiles = draw_quantiles,
                                 na.rm = na.rm,
                                 ...))
  }
  
}


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
shapes<- c("VP" = 15,
          "VPC" = 17,
          "Diff" = 19)
# shapes<- c("VP" = 22,
#            "VPC" = 24,
#            "Diff" = 23)

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
                        axis.title = element_blank(),
                        panel.grid = element_blank())

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
  scale_shape_manual(values = shapes)+
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
  scale_shape_manual(values = shapes)+
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
  scale_shape_manual(values = shapes)+
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


mockplots_og<- plot_grid(lm_1, lm_2, lm_3,
          vpcl_1, vpcl_2, vpcl_3,
          nrow = 2) #PLOT: LM VS VIPCAL MOCKPLOTS


mockplots_og

ggsave(filename = "mockplots_og.svg",
       plot = mockplots_og,
       dpi = 1000,
       width = 150 ,
       height = 100,
       unit = 'mm' )

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
wilcox.test(VP ~ as.factor(VP_Type), data = vis_df_og %>% filter(Time_Range != 'T1_T2'))
#W = 34050432, p-value < 2.2e-16 Significantly different.

#Linear Regression

vis_df_og_lm<- vis_df_og  %>%
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


#Cases where VPCL is lesser than LM
length(which(vis_df_og_lm$LM_SR_AVG > vis_df_og_lm$VPCL_AR_DIFF))
#511 3.406667%

#Cases when VPCL is higher than LM
length(which(vis_df_og_lm$LM_SR_AVG < vis_df_og_lm$VPCL_AR_DIFF))
#13243 88.28667%

#Cases when VPCL is equal to LM
length(which(vis_df_og_lm$LM_SR_AVG ==vis_df_og_lm$VPCL_AR_DIFF))
#1246 8.306667%

vis_df_og_lm$Sample_Type <- factor(vis_df_og_lm$Sample_Type,
                                   levels = c("VP", "VPC", "Diff"))



ggplot(data = vis_df_og_lm %>% arrange((Sample_Type)),
       aes(x = LM_SR_AVG,
           y = VPCL_AR_DIFF,
           col = Sample_Type,
           shape = Sample_Type))+
  geom_point(alpha = 0.1)+
  #geom_smooth(method = 'lm')+
  scale_color_manual(values = cols)+
  scale_shape_manual(values = shapes)+
  theme_bw()+
  #geom_abline(slope = 1)+
  xlab("Linear Regression (LM-3)")+
  ylab("VIPCAL (VPCL-4)")+
  theme(strip.background = element_rect(fill = lighten("#323d5e",
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
        axis.text = element_text(size = 10),
        panel.grid = element_blank() )+
  #facet_grid(Time_Range ~ Sample_Type)
  facet_wrap(vars(Time_Range),
             nrow = 1)


#Dataframe for linear vs VIP[CAL comparison

df_lm_vpcl<- vis_df_og_lm %>%
  filter(Time_Range != 'T1_T2')

annotations <-  df_lm_vpcl %>%
  group_by(!!sym("Sample_Type")) %>%
  do(model = lm(as.formula(paste("VPCL_AR_DIFF", "~", "LM_SR_AVG")), data = .)) %>%
  rowwise() %>%
  mutate(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    r2 = summary(model)$r.squared,
    eq1 = sprintf("y == %.2f*x + %.2f", slope, intercept), 
    eq2 = sprintf("R^2 == %.2f", r2),
    eq_x = min(vis_df_og_lm$LM_SR_AVG), 
    eq_y1 = max(vis_df_og_lm$VPCL_AR_DIFF), # y position for the first equation
    eq_y2 = max(vis_df_og_lm$VPCL_AR_DIFF) - 0.5  # slightly lower y position for the second equation
  ) %>%
  select(!!sym("Sample_Type"), eq1, eq2, eq_x, eq_y1, eq_y2)



lm_vs_vpcl<- ggplot(data = df_lm_vpcl,
       aes(x = LM_SR_AVG,
           y = VPCL_AR_DIFF,
           col = Sample_Type,
           shape = Sample_Type))+
  geom_point(alpha = 0.05)+
  geom_smooth(method = 'lm', se = T)+
  #annotate("text", x = eq_x, y = eq_y, label = eq, hjust = 0, vjust = 1) +
  scale_fill_manual(values = colorspace::lighten(cols, 0.2))+
  scale_colour_manual(values = cols,
                      breaks = c('VP', 'VPC', 'Diff'),
                      name= 'Treatment') +
  scale_shape_manual(values = shapes,
                     breaks = c('VP', 'VPC', 'Diff'),
                     name= 'Treatment')+
  theme_bw()+
  #geom_abline(slope = 1)+
  xlab("Linear Regression (LM-3)")+
  ylab("VIPCAL (VPCL-4)")+
  theme(strip.background = element_rect(fill = lighten("#323d5e",
                                                       amount = 0.0) ,
                                        color = NA),
        strip.text = element_text(face = 'bold',
                                  color = 'white',
                                  size = 10),
        panel.border = element_rect(linewidth = 2),
        panel.background = element_blank(),
        #legend.position = c(0.75, 0.1),
        legend.title = element_text(face = 'bold',
                                    size = 10),
        legend.text = element_text(size = 9),
        axis.title = element_text(face = 'bold',
                                  size = 10),
        axis.text = element_text(size = 10),
        panel.grid = element_blank())+
  coord_fixed(ratio = 1)+
  facet_wrap(vars(Sample_Type),
             nrow = 1) +
  geom_text(data = annotations, aes(x = eq_x, y = eq_y1, label = eq1),
            hjust = 0, vjust = 1, parse = TRUE, inherit.aes = FALSE) +
  geom_text(data = annotations, aes(x = eq_x, y = eq_y2, label = eq2),
            hjust = 0, vjust = 1, parse = TRUE, inherit.aes = FALSE)
lm_vs_vpcl

ggsave(filename = "1_LM_VS_VPCL_Regression.svg",
       plot = lm_vs_vpcl,
       dpi = 1000,
       width = 150 ,
       unit = 'mm' )


#### 12 Methods Violin plot ####


vis_df<- read.csv("results/simu_vp_filtered.csv")
vis_df <- vis_df  %>%
  mutate(Station_Number = as.numeric(Station_Number),
         VP_Type = factor(VP_Type,
                          levels = c("LM_AP", "LM_SR_AVG",
                                     "LM_AR", "LM_AR_DIFF",
                                     "LM_AR_DIFF_LMER", "VPCL_SR_AVG",
                                     "VPCL_AR", "VPCL_AR_SE",
                                     "VPCL_AR_DIFF", "VPCL_AR_DIFF_SE",
                                     "VPCL_AR_DIFF_LMER", "VPCL_AR_DIFF_LMER_SE"),
                          labels = c("LM-1", "LM-2", "LM-3", "LM-4", "LM-5",
                                     "VPCL-1", "VPCL-2", "VPCL-3", "VPCL-4", "VPCL-5", "VPCL-6", "VPCL-7")))
str(vis_df)
unique(vis_df$Population)
unique(vis_df$VP_Type)


method12<- ggplot(data = vis_df %>%
                    filter(Time_Range != "T1_T2"),
       aes(x = VP_Type,
           y = VP,
           fill = VP_Type)) +
  geom_violin()+
  scale_fill_manual(values = c("#e19c33", "#ff9c00", "#ba4a2e", "#4787a1",
                               "#423a66", "#514759", "#3f5428", "#971e28",
                               "#85b0a3", "#de5168", "#e87d7f", "#514759"))+
  geom_boxplot(width = 0.05,
               outlier.shape = NA, 
               outlier.alpha = 0.00,
               fill = "white")+
  xlab("Viral Production Analyses Methods")+
  ylab("Viral Production Rate")+
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust = 1),
        panel.border = element_rect(linewidth = 2),
        #legend.position = c(0.75, 0.1),
        legend.title = element_text(face = 'bold',
                                    size = 10),
        legend.text = element_text(size = 9),
        axis.title = element_text(face = 'bold',
                                  size = 10),
        axis.text = element_text(size = 10,
                                 face = 'bold'),
        legend.position = 'none')

method12

ggsave(filename = "12_methods.svg",
       plot = method12,
       dpi = 1000,
       width = 250 ,
       height = 125,
       unit = 'mm' )



#Perfeom Kruskal Wallis test between VIPCAL and VIPCAL-SE for VP, VPC and Diff 
v_df <- vis_df %>%
  filter(Time_Range != "T1_T2")

kruskal.test(VP ~ VP_Type, data = v_df)

library(FSA)
dunn_test<- dunnTest(VP ~ VP_Type, data = v_df,
                     method = 'holm')
dunn_output<- as.data.frame(as.list(dunn_test)[[2]])

dunn2<- dunn_output

colnames(dunn2)



dunn2<- dunn2 %>% separate(Comparison, c("Comp1", "Comp2"), " - ")



dunn2<- dunn2 %>%
  arrange(Comp1, Comp2)

corr_p_plot<- ggplot(data = dunn2 , aes(x = Comp1, y = Comp2, fill = P.adj, color = P.adj))+
  geom_tile()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust = 1),
        axis.text = element_text(size = 10,
                                 face = 'bold'),
        axis.title = element_blank(),
        legend.position = c(0.8, 0.3),
        legend.title = element_text(size = 8,
                                    face = 'bold'),
        legend.text = element_text(size = 8),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
        panel.grid = element_blank())+
  labs(fill = 'P(adjusted)',
       color = 'P(adjusted)')+
  scale_fill_gradient(high = '#ea7e82ff', low = '#52495aff')+
  scale_color_gradient(high = '#ea7e82ff', low = '#52495aff')+
  coord_fixed(ratio = 1)

corr_p_plot

ggsave(filename = "corr_p_plot.svg",
       plot = corr_p_plot,
       dpi = 1000,
       width = 125 ,
       height = 125,
       unit = 'mm' )


####VIPCAL vs VIPCAL-SE####

vpcl_df<- vis_df %>% 
  filter(VP_Type %in% c("VPCL-4",
                        "VPCL-7")) %>%
  mutate(Sample_Type = factor(Sample_Type, levels = c("VP", "VPC", "Diff"))) %>%
  filter(Time_Range != "T1_T2")

#Split Violin
vpcl_vs <- ggplot(data = vpcl_df,
       aes(x = Sample_Type,
           y = VP,
           fill = VP_Type)) +
  geom_split_violin(nudge = 0.025,
                    colour = NA)+
  scale_fill_manual(values = c("#de5168", "#514759"),
                    name = "Viral Production\nData Analyses Methods")+
  xlab("Treatments")+
  ylab("Viral Production Rates (per L per h)")+
  theme_bw()+
  theme(panel.border = element_rect(linewidth = 2),
        panel.background = element_blank(),
        legend.position = "top",
        legend.title = element_text(face = 'bold',
                                    size = 8),
        legend.text = element_text(size = 7),
        axis.title = element_text(face = 'bold',
                                  size = 10),
        axis.text = element_text(size = 10),
        panel.grid = element_blank())

vpcl_vs

ggsave(filename = "vpcl_vs.svg",
       plot = vpcl_vs,
       dpi = 1000,
       width = 100,
       height = 120,
       unit = 'mm' )

#Linear Regression

vpcl_df_lm<- vpcl_df %>%
  select(-c( "abs_VP", "VP_SE", "VP_R_Squared")) %>%
  pivot_wider(names_from = VP_Type,
              values_from = VP) %>%
  mutate(Sample_Type = factor(Sample_Type,
                              levels = c("VP", "VPC", "Diff")))

vpcl_lm_plot<- ggplot(data = vpcl_df_lm,
       aes(x = `VPCL-4`,
           y = `VPCL-7`,
           color = Sample_Type,
           shape = Sample_Type))+
  geom_point(alpha = 0.1)+
  scale_color_manual(values = cols,
                     name = "Treatment")+
  scale_shape_manual(values = shapes,
                     name = "Treatment")+
  geom_smooth(method = 'lm',
              alpha =  0.2)+
  geom_abline(slope = 1, linewidth = 1) + 
  xlab("VIPCAL (VPCL-4)")+
  ylab("VIPCAL-SE (VPCL-7)")+
  theme_bw()+
  theme(panel.border = element_rect(linewidth = 2),
        panel.background = element_blank(),
        legend.position = c(0.9,0.2),
        legend.title = element_text(face = 'bold',
                                    size = 8),
        legend.text = element_text(size = 7),
        axis.title = element_text(face = 'bold',
                                  size = 10),
        axis.text = element_text(size = 10),
        panel.grid = element_blank())+
  coord_fixed(ratio = 0.9)

vpcl_lm_plot


ggsave(filename = "vpcl_lm_plot.svg",
       plot = vpcl_lm_plot,
       dpi = 1000,
       width = 100,
       height = 100,
       unit = 'mm' )

#Perfeom Kruskal Wallis test between VIPCAL and VIPCAL-SE for VP, VPC and Diff 
v_df <- make_a_csv
v_df$VP_Type<- factor(v_df$VP_Type, levels = c('LM_AP', 'LM_SR_AVG',
                                             'LM_AR', 'LM_AR_Diff',
                                             'LM_AR_Diff_LMER', "VPCL_SR_AVG" ,
                                             "VPCL_AR_No_SE" , "VPCL_AR_SE" ,
                                             "VPCL_AR_Diff_No_SE" ,
                                             "VPCL_AR_Diff_SE"  ,
                                             "VPCL_AR_Diff_LMER_No_SE", "VPCL_AR_Diff_LMER_SE"
),
labels = c('LM-1', 'LM-2', 'LM-3', 'LM-4', 'LM-5',
           'VPCL-1', 'VPCL-2', 'VPCL-3',
           'VPCL-4', 'VPCL-5', 'VPCL-6',
           'VPCL-7'))


kruskal.test(VP ~ VP_Type, data = v_df)

library(FSA)
dunn_test<- dunnTest(VP ~ VP_Type, data = v_df,
         method = 'holm')
dunn_output<- as.data.frame(as.list(dunn_test)[[2]])

dunn2<- dunn_output

colnames(dunn2)



dunn2<- dunn2 %>% separate(Comparison, c("Comp1", "Comp2"), " - ")



dunn2<- dunn2 %>%
  arrange(Comp1, Comp2)

corr_p_plot<- ggplot(data = dunn2 , aes(x = Comp1, y = Comp2, fill = P.adj, color = P.adj))+
  geom_tile()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 10,
                                 face = 'bold'),
        axis.title = element_blank(),
        legend.position = c(0.8, 0.3),
        legend.title = element_text(size = 8,
                                    face = 'bold'),
        legend.text = element_text(size = 8),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
        panel.grid = element_blank())+
  labs(fill = 'P(adjusted)',
       color = 'P(adjusted)')+
  scale_fill_gradient(high = '#ea7e82ff', low = '#52495aff')+
  scale_color_gradient(high = '#ea7e82ff', low = '#52495aff')

corr_p_plot



dunn2<- dunn2 %>%
  select(-c(Z, P.unadj)) %>%
  pivot_wider(names_from = Comp2,
              values_from = P.adj) %>%
  column_to_rownames(var = 'Comp1')%>%
  as.matrix()

ggcorrplot::ggcorrplot(dunn2)



