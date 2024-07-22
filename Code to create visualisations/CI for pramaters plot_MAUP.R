ci_data_MAUP_intercept<- data.frame(level= c("Mesh",
                                   "SA1","SA2",
                                   "SA3", "SA4"),
                          Median = c(-0.45,-0.0037,-0.0052,-0.0005,-0.0018),
                          lower= c(-0.47,-0.0178,-0.02,-0.017,-0.0222),
                          upper= c(-0.43,0.0103,0.0093,0.0164,0.0198),
                          Median_1 = c(0.08,0.002,0.007,0.0125,0.0108),
                          lower_1= c(0.07,0.0013,0.0035,0.0042,0.0029),
                          upper_1= c(0.09,0.0032,0.0132,0.0277,0.0435),
                          Median_2 = c(0.0013,0.0009,0.0031,0.0048,0.019),
                          lower_2= c(0.0007,0.0006,0.0014,0.0018,0.0028),
                          upper_2= c(0.0023,0.0015,0.0058,0.0097,0.0247))

library(tidyverse)
p1<-ggplot(ci_data_MAUP_intercept, aes(level, Median)) + geom_point() + 
  geom_errorbar(aes(ymin = lower, ymax = upper))+labs(x="Aggregate level",
                                                      y= expression(alpha))+
  theme_bw()
p2<-ggplot(ci_data_MAUP_intercept, aes(level, Median_1)) + geom_point() + 
  geom_errorbar(aes(ymin = lower_1, ymax = upper_1))+labs(x="Aggregate level",
                                                          y= expression(sigma[u]^2))+
  theme_bw()
p3<-ggplot(ci_data_MAUP_intercept, aes(level, Median_2)) + geom_point() + 
  geom_errorbar(aes(ymin = lower_2, ymax = upper_2))+labs(x="Aggregate level",
                                                          y=expression(sigma[v]^2))+
  theme_bw()
library(ggpubr)
ggarrange(p1,p2,p3,nrow = 1,ncol=3)

ci_data_MAUP_beta<- data.frame(level= c("Mesh",
                                             "SA1","SA2",
                                             "SA3", "SA4"),
                                    Median_1 = c(-0.33,0.12,0.19,0.10,0.13),
                                    lower_1= c(-0.55,0.09,0.17,0.06,0.05),
                                    upper_1= c(-0.10,0.15,0.23,0.15,0.22),
                               Median_2 = c(-1.71,-0.07,-0.13,-0.05,-0.08),
                               lower_2= c(-1.95,-0.11,-0.17,-0.11,-0.21),
                               upper_2= c(-1.47,-0.03,-0.09,0.017,0.04),
                               Median_3 = c(-0.61,-0.12,-0.23,-0.06,-0.14),
                               lower_3= c(-0.83,-0.17,-0.27,-0.13,-0.29),
                               upper_3= c(-0.39,-0.08,-0.19,0.005,0.008),
                               Median_4 = c(-0.22,-0.18,-0.29,-0.15,-0.19),
                               lower_4= c(-0.44,-0.24,-0.33,-0.23,-0.33),
                               upper_4= c(-0.01,-0.15,-0.24,-0.08,-0.06),
                               Median_5 = c(0.02,-0.28,-0.46,-0.28,-0.31),
                               lower_5= c(-0.20,-0.32,-0.51,-0.34,-0.46),
                               upper_5= c(0.24,-0.22,-0.41,-0.199,-0.15),
                               Median_6 = c(0.26,0.002,0.0015,0.004,0.008),
                               lower_6= c(0.23,0.001,0.008,0.002,0.002),
                               upper_6= c(0.28,0.003,0.003,0.009,0.024),
                               Median_7 = c(0.001,0.0009,0.0012,0.003,0.005),
                               lower_7= c(0.0007,0.0006,0.0007,0.001,0.002),
                               upper_7= c(0.0021,0.0015,0.0019,0.005,0.012))

p1<-ggplot(ci_data_MAUP_beta, aes(level, Median_1)) + geom_point() + 
  geom_errorbar(aes(ymin = lower_1, ymax = upper_1))+labs(x="Aggregate level",
                                                      y= expression(alpha))+
  theme_bw()
p2<-ggplot(ci_data_MAUP_beta, aes(level, Median_2)) + geom_point() + 
  geom_errorbar(aes(ymin = lower_2, ymax = upper_2))+labs(x="Aggregate level",
                                                          y= expression(beta[1]))+
  theme_bw()
p3<-ggplot(ci_data_MAUP_beta, aes(level, Median_3)) + geom_point() + 
  geom_errorbar(aes(ymin = lower_3, ymax = upper_3))+labs(x="Aggregate level",
                                                          y= expression(beta[2]))+
  theme_bw()
p4<-ggplot(ci_data_MAUP_beta, aes(level, Median_4)) + geom_point() + 
  geom_errorbar(aes(ymin = lower_4, ymax = upper_4))+labs(x="Aggregate level",
                                                          y= expression(beta[3]))+
  theme_bw()
p5<-ggplot(ci_data_MAUP_beta, aes(level, Median_5)) + geom_point() + 
  geom_errorbar(aes(ymin = lower_5, ymax = upper_5))+labs(x="Aggregate level",
                                                          y= expression(beta[4]))+
  theme_bw()
p6<-ggplot(ci_data_MAUP_beta, aes(level, Median_6)) + geom_point() + 
  geom_errorbar(aes(ymin = lower_6, ymax = upper_6))+labs(x="Aggregate level",
                                                          y= expression(sigma[u]^2))+
  theme_bw()
p7<-ggplot(ci_data_MAUP_beta, aes(level, Median_7)) + geom_point() + 
  geom_errorbar(aes(ymin = lower_7, ymax = upper_7))+labs(x="Aggregate level",
                                                          y= expression(sigma[v]^2))+
  theme_bw()
library(ggpubr)
ggarrange(p1,p2,p3,p4,p5,p6,p7,nrow = 2,ncol=4)
