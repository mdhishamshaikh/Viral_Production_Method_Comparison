library(tidyverse)
library(ggsci)

df<- data.frame(Timepoint = rep(1:6,6),
                Sample_Type = c(rep('VP', 18),
                rep('VPC', 18)),
                Replicate = c(rep(1, 6), rep(2, 6), rep(3, 6),
                              rep(4, 6), rep(5, 6), rep(6, 6)),
                Count = c(1, 2, 5, 6, 7, 4,
                          1.1, 2.3, 5.3, 6.3, 7.5, 4.3,
                          0.9, 1.9, 5.8, 5.2, 6.9, 4.1,
                          1.2, 3.1, 4.0, 6.2, 5.4, 3.0,
                          1.3, 3.4, 3.9, 6.4, 5.9, 3.4,
                          0.8, 3.7, 3.6, 6.1, 5.8, 3.5)
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
  


bp<- ggplot(data = df, aes(x = as.factor(Timepoint), y = Count, color = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = 0.5)+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  scale_color_lancet()+
  theme(legend.position = 'bottom')
  
bp1<- ggplot(data = df, aes(x = as.factor(Timepoint), y = Count, color = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = 0.1)+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments',  shape = 'Treatments')+
  scale_color_lancet()+
  theme(legend.position = 'bottom')

bp2<- ggplot(data = df2, aes(x = as.factor(Timepoint), y = mean, color = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = 0.5)+
  geom_errorbar(data = df2, aes(ymin=mean-se, ymax=mean+se), width=.2)+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  scale_color_lancet()+
  theme(legend.position = 'bottom')

bp2<- ggplot(data = df3, aes(x = as.factor(Timepoint), y = mean, color = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = 0.5)+
  geom_errorbar(data = df3, aes(ymin=mean-se, ymax=mean+se), width=.2)+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  scale_color_lancet()+
  theme(legend.position = 'bottom')

bp3<- ggplot(data = df3, aes(x = as.factor(Timepoint), y = mean, color = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = 0.5)+
  #geom_errorbar(data = df3, aes(ymin=mean-se, ymax=mean+se), width=.2)+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  scale_color_lancet()+
  theme(legend.position = 'bottom')

#0_Data
bp

#1_ LM_AP
bp +
  geom_smooth(aes(group = Sample_Type) , method = 'lm', se = F)

#2_LM_SR_AVG
bp+
  geom_smooth(aes(group = Replicate), method = 'lm', se = F, alpha = 0.1)

#3_LM_AR
bp1 +
  geom_point(data = df2, aes(x = as.factor(Timepoint), y = mean, shape = Sample_Type))+
  geom_smooth(data = df2, aes(x = as.factor(Timepoint), y = mean, group = Sample_Type), method = 'lm', se = F, alpha = 1)

#4_LM_AR_Diff_5_LM_aR_Diff_LMER

bp1+
  geom_point(data = df3, aes(x = as.factor(Timepoint), y = mean, shape = Sample_Type))+
  geom_smooth(data = df3, aes(x = as.factor(Timepoint), y = mean, group = Sample_Type), method = 'lm', se = F, alpha = 1)+
  scale_colour_manual(breaks = c('VP','VPC','Diff'),
                      values = c('#00468BFF', '#ED0000FF', '#42B540FF'))+
  scale_shape_discrete(breaks = c('VP','VPC','Diff'))


#6_VPCL_SR_AVG

bp +
  geom_line(aes(group = Replicate), lwd = 0.7)

#7_VPCL_AR_No_SE

bp1 +
  geom_point(data = df2, aes(x = as.factor(Timepoint), y = mean, shape = Sample_Type))+
  geom_line(data = df2, aes(x = as.factor(Timepoint), y = mean, group = Sample_Type), lwd = 0.7)

#8_VPCL_AR_SE

bp1 +
  geom_point(data = df2, aes(x = as.factor(Timepoint), y = mean, shape = Sample_Type))+
  geom_errorbar(data = df2, aes(x = as.factor(Timepoint), y = mean, shape = Sample_Type,ymin=mean-se, ymax=mean+se), width=.2)+
  geom_line(data = df2, aes(x = as.factor(Timepoint), y = mean, group = Sample_Type), lwd = 0.7)
  

#9_VPCL_AR_Diff_No_SE

bp3 +
  geom_point(data = df, aes(x = as.factor(Timepoint), y = Count, shape = Sample_Type), alpha = 0.1)+ 
  geom_line(aes(x =  as.factor(Timepoint), y = mean, group = Sample_Type), lwd = 0.7)

#10_VPCL_AR_Diff_SE

bp2 +
  geom_point(data = df, aes(x = as.factor(Timepoint), y = Count, shape = Sample_Type), alpha = 0.1)+
  geom_line(aes(x =  as.factor(Timepoint), y = mean, group = Sample_Type), lwd = 0.7)


ggplot(data = df , aes(x = as.factor(Timepoint), y = Count, color = Sample_Type))+
  geom_point(alpha = 0.5)+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  
  theme(legend.position = 'bottom')+
    geom_smooth(method = 'lm', aes(group = Sample_Type), se = F)+
  scale_color_lancet()





bp2

bp3

bp  +
  geom_line(aes(group = Replicate))
