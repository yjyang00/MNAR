RMSE.MB05 = readRDS('finalmethod05/RMSE_MB05_NEW.rds')
RMSE.MH05 = readRDS('finalmethod05/RMSE_MH05.rds')
RMSE.MMH05 = readRDS('finalmethod05/RMSE_MMH05.rds')
RMSE.IH05 = readRDS('finalmethod05/RMSE_IH05.rds')

RMSE.MB = readRDS('finalmethod2/RMSE_MB_NEW.rds')
RMSE.MH = readRDS('finalmethod2/RMSE_MH.rds')
RMSE.MMH = readRDS('finalmethod2/RMSE_MMH.rds')
RMSE.IH = readRDS('finalmethod2/RMSE_IH.rds')

hybrid05 = readRDS('finalmethod05/hybrid_RES05_NEW.rds')
MH05 = readRDS('finalmethod05/MH_RES05.rds')
MMH05 = readRDS('finalmethod05/MMH_RES05.rds')
IH05 = readRDS('finalmethod05/IH_RES05.rds')

hybrid = readRDS('finalmethod2/hybrid_RES_NEW.rds')
MH = readRDS('finalmethod2/MH_RES.rds')
MMH = readRDS('finalmethod2/MMH_RES.rds')
IH = readRDS('finalmethod2/IH_RES.rds')

CI.hybrid05 = c(hybrid05$CI.MDAM$CI.T, hybrid05$CI.MDAM$CI.joint, hybrid05$CI.MDAM$CI.cond)
CI.MH05 = c(MH05$CI.MDAM$CI.T, MH05$CI.MDAM$CI.joint, MH05$CI.MDAM$CI.cond)
CI.MMH05 = c(MMH05$CI.MDAM$CI.T, MMH05$CI.MDAM$CI.joint, MMH05$CI.MDAM$CI.cond)
CI.IH05 = c(IH05$CI.MDAM$CI.T, IH05$CI.MDAM$CI.joint, IH05$CI.MDAM$CI.cond)

CI.hybrid = c(hybrid$CI.MDAM$CI.T, hybrid$CI.MDAM$CI.joint, hybrid$CI.MDAM$CI.cond)
CI.MH = c(MH$CI.MDAM$CI.T, MH$CI.MDAM$CI.joint, MH$CI.MDAM$CI.cond)
CI.MMH = c(MMH$CI.MDAM$CI.T, MMH$CI.MDAM$CI.joint, MMH$CI.MDAM$CI.cond)
CI.IH = c(IH$CI.MDAM$CI.T, IH$CI.MDAM$CI.joint, IH$CI.MDAM$CI.cond)

#### plot
df_RMSE_2 = data.frame(CI = c(CI.hybrid, CI.MH, CI.MMH, CI.IH), 
                     RMSE = c(RMSE.MB$rRMSE_T.MB, RMSE.MB$rRMSE_joint.MB, RMSE.MB$rRMSE_cond.MB,
                              RMSE.MH$rRMSE_T.MH, RMSE.MH$rRMSE_joint.MH, RMSE.MH$rRMSE_cond.MH,
                              RMSE.MMH$rRMSE_T.MMH, RMSE.MMH$rRMSE_joint.MMH, RMSE.MMH$rRMSE_cond.MMH,
                              RMSE.IH$rRMSE_T.IH, RMSE.IH$rRMSE_joint.IH, RMSE.IH$rRMSE_cond.IH))
df_RMSE_2$Method = c(rep("MB", 18), rep("MH", 18), rep("MMH", 18), rep("IH", 18))
df_RMSE_2$Type = c(rep("total", 6), rep("joint", 4), rep("conditional", 8))

df_RMSE_05 = data.frame(CI = c(CI.hybrid05, CI.MH05, CI.MMH05, CI.IH05), 
                       RMSE = c(RMSE.MB05$rRMSE_T.MB, RMSE.MB05$rRMSE_joint.MB, RMSE.MB05$rRMSE_cond.MB,
                                RMSE.MH05$rRMSE_T.MH, RMSE.MH05$rRMSE_joint.MH, RMSE.MH05$rRMSE_cond.MH,
                                RMSE.MMH05$rRMSE_T.MMH, RMSE.MMH05$rRMSE_joint.MMH, RMSE.MMH05$rRMSE_cond.MMH,
                                RMSE.IH05$rRMSE_T.IH, RMSE.IH05$rRMSE_joint.IH, RMSE.IH05$rRMSE_cond.IH))
df_RMSE_05$Method = c(rep("MB", 18), rep("MH", 18), rep("MMH", 18), rep("IH", 18))
df_RMSE_05$Type = c(rep("total", 6), rep("joint", 4), rep("conditional", 8))

library(latex2exp)

# MMH Vs IH
plot1 = ggplot(df_RMSE_2 %>% filter(Method %in% c('MMH','IH')), aes(x = RMSE, y = CI, color=Method, shape=Type))+
  geom_point()+
  labs(x=TeX("Relative RMSE: $theta_1$=-2.0"), y="Coverage")+
  theme_classic()

plot2 = ggplot(df_RMSE_05%>% filter(Method %in% c('MMH','IH')), aes(x = RMSE, y = CI, color=Method, shape=Type))+
  geom_point()+
  labs(x=TeX("Relative RMSE: $theta_1$=-0.5"), y="Coverage")+
  theme_classic() + theme(legend.position = "bottom")

plot_legend = lemon::g_legend(plot2 + guides(colour = guide_legend(nrow = 1)) +  guides(shape = guide_legend(nrow = 1)))

gridExtra::grid.arrange(gridExtra::arrangeGrob(plot1 + theme(legend.position="none"), 
                                               plot2 + theme(legend.position="none"),
                                               nrow = 2),
                        plot_legend, nrow = 2, heights = c(5.5, 1))

# MMH, MH, MB
plot3 = ggplot(df_RMSE_2 %>% filter(Method !='IH'), aes(x = RMSE, y = CI, color=Method, shape=Type))+
  geom_point()+
  labs(x=TeX("Relative RMSE: $theta_1$=-2.0"), y="Coverage")+
  theme_classic()

plot4 = ggplot(df_RMSE_05%>% filter(Method != 'IH'), aes(x = RMSE, y = CI, color=Method, shape=Type))+
  geom_point()+
  labs(x=TeX("Relative RMSE: $theta_1$=-0.5"), y="Coverage")+
  theme_classic() + theme(legend.position = "bottom")

plot_legend = lemon::g_legend(plot4 + guides(colour = guide_legend(nrow = 1)) +  guides(shape = guide_legend(nrow = 1)))

gridExtra::grid.arrange(gridExtra::arrangeGrob(plot3 + theme(legend.position="none"), 
                                               plot4 + theme(legend.position="none"),
                                               nrow = 2),
                        plot_legend, nrow = 2, heights = c(5.5, 1))

# together

library(ggplot2)
library(gridExtra)
library(lemon)

# Create the plots
# custom_colors = c('MMH' = rgb(0, 114, 178, maxColorValue = 255),
#                    'IH' = rgb(213, 94, 0, maxColorValue = 255),
#                    'MH' = rgb(0, 158, 115, maxColorValue = 255),
#                    'MB' = rgb(204, 121, 167, maxColorValue = 255))
custom_colors = c('MMH' = '#e02b35',
                  'MH' = '#e9c716',
                  'IH' = '#1a80bb',
                  'MB' = '#082a54')


plot1 <- ggplot(df_RMSE_2 %>% filter(Method %in% c('MMH','IH')), aes(x = RMSE, y = CI, color=Method, shape=Type))+
  geom_point()+
  scale_color_manual(values = custom_colors) +
  labs(x=TeX("Relative RMSE: $theta_1$=-2.0"), y="Coverage")+
  theme_classic()

plot2 <- ggplot(df_RMSE_05 %>% filter(Method %in% c('MMH','IH')), aes(x = RMSE, y = CI, color=Method, shape=Type))+
  geom_point()+
  scale_color_manual(values = custom_colors) +
  labs(x=TeX("Relative RMSE: $theta_1$=-0.5"), y="Coverage")+
  theme_classic() + theme(legend.position = "bottom")

plot_legend1 <- lemon::g_legend(plot2 + guides(colour = guide_legend(nrow = 1)) +  guides(shape = guide_legend(nrow = 1)))

plot3 <- ggplot(df_RMSE_2 %>% filter(Method !='IH'), aes(x = RMSE, y = CI, color=Method, shape=Type))+
  geom_point()+
  scale_color_manual(values = custom_colors) +
  labs(x=TeX("Relative RMSE: $theta_1$=-2.0"), y="Coverage")+
  theme_classic()

plot4 <- ggplot(df_RMSE_05 %>% filter(Method != 'IH'), aes(x = RMSE, y = CI, color=Method, shape=Type))+
  geom_point()+
  scale_color_manual(values = custom_colors) +
  labs(x=TeX("Relative RMSE: $theta_1$=-0.5"), y="Coverage")+
  theme_classic() + theme(legend.position = "bottom")

plot_legend2 <- lemon::g_legend(plot4 + guides(colour = guide_legend(nrow = 1)) +  guides(shape = guide_legend(nrow = 1)))

# Arrange the plots and legends
arrangement1 <- gridExtra::grid.arrange(
  gridExtra::arrangeGrob(plot1 + theme(legend.position="none"), 
                         plot2 + theme(legend.position="none"),
                         nrow = 2),
  plot_legend1, nrow = 2, heights = c(5.5, 1))

arrangement2 <- gridExtra::grid.arrange(
  gridExtra::arrangeGrob(plot3 + theme(legend.position="none"), 
                         plot4 + theme(legend.position="none"),
                         nrow = 2),
  plot_legend2, nrow = 2, heights = c(5.5, 1))

# Combine the two arrangements into one grid with two columns
final_plot <- gridExtra::grid.arrange(arrangement1, arrangement2, ncol=2)

# Display the final plot
print(final_plot)


