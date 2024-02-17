##### Inform the directory 
setwd("D:/Exportações/Chap3/Sim1")
getwd()

##### Some packages
library(ggplot2)
library(ggpubr)
# library(grid)
# library(gridExtra)
library(ggrepel) # Para colocar r?tulos/etiquetas/nomes das esta??es nos mapas

smp_size = 1000
N = 17
K = 100

##### Time series for gauged sites
Y = read.table("Porcentagem_1/Y.txt",
               head = FALSE,
               sep = ",")
colnames(Y) = c("Time", "Site", "Resp.1", "Resp.2")
e_Y_MA = read.table("Porcentagem_1/e_Y.txt",
                      head = FALSE,
                      sep = ",")
colnames(e_Y_MA) = c("Iteration", "Time", "Site", "Resp.1", "Resp.2")
e_Y_MI = read.table("Porcentagem_1/e_Y.txt",
                      head = FALSE,
                      sep = ",")
colnames(e_Y_MI) = c("Iteration", "Time", "Site", "Resp.1", "Resp.2")

n = 10

MA_stats_Y1 = data.frame(1:K, matrix(data = 0, nrow = K, ncol = 5))
MA_stats_Y2 = data.frame(1:K, matrix(data = 0, nrow = K, ncol = 5))
colnames(MA_stats_Y1) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(MA_stats_Y2) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
MI_stats_Y1 = data.frame(1:K, matrix(data = 0, nrow = K, ncol = 5))
MI_stats_Y2 = data.frame(1:K, matrix(data = 0, nrow = K, ncol = 5))
colnames(MI_stats_Y1) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(MI_stats_Y2) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
for(t in 1:K) {
    set_MA = subset(x = e_Y_MA,
                    select = c("Iteration", "Resp.1", "Resp.2"),
                    subset = Site == n & Time == t)
    MA_stats_Y1[t,"LI"] = quantile(x = set_MA[1:smp_size, "Resp.1"],
                                   probs = 0.025)
    MA_stats_Y2[t,"LI"] = quantile(x = set_MA[1:smp_size, "Resp.2"],
                                   probs = 0.025)
    MA_stats_Y1[t,"Mean"] = mean(set_MA[1:smp_size, "Resp.1"])
    MA_stats_Y2[t,"Mean"] = mean(set_MA[1:smp_size, "Resp.2"])
    MA_stats_Y1[t,"True"] = subset(x = Y,
                                   select = c("Resp.1"),
                                   subset = Site == n & Time == t)
    MA_stats_Y2[t,"True"] = subset(x = Y,
                                   select = c("Resp.2"),
                                   subset = Site == n & Time == t)
    MA_stats_Y1[t,"LS"] = quantile(x = set_MA[1:smp_size, "Resp.1"],
                                   probs = 0.975)
    MA_stats_Y2[t,"LS"] = quantile(x = set_MA[1:smp_size, "Resp.2"],
                                   probs = 0.975)
    set_MI = subset(x = e_Y_MI,
                    select = c("Iteration", "Resp.1", "Resp.2"),
                    subset = Site == n & Time == t)
    MI_stats_Y1[t,"LI"] = quantile(x = set_MI[1:smp_size, "Resp.1"],
                                   probs = 0.025)
    MI_stats_Y2[t,"LI"] = quantile(x = set_MI[1:smp_size, "Resp.2"],
                                   probs = 0.025)
    MI_stats_Y1[t,"Mean"] = mean(set_MI[1:smp_size, "Resp.1"])
    MI_stats_Y2[t,"Mean"] = mean(set_MI[1:smp_size, "Resp.2"])
    MI_stats_Y1[t,"True"] = subset(x = Y,
                                   select = c("Resp.1"),
                                   subset = Site == n & Time == t)
    MI_stats_Y2[t,"True"] = subset(x = Y,
                                   select = c("Resp.2"),
                                   subset = Site == n & Time == t)
    MI_stats_Y1[t,"LS"] = quantile(x = set_MI[1:smp_size, "Resp.1"],
                                   probs = 0.975)
    MI_stats_Y2[t,"LS"] = quantile(x = set_MI[1:smp_size, "Resp.2"],
                                   probs = 0.975) 
}

##### Coverage

MA_stats_Y1[,"Coverage"] = MA_stats_Y1[,"True"] == MA_stats_Y1[,"Mean"]
table(MA_stats_Y1[,"Coverage"])

MI_stats_Y1[,"Coverage"] = MI_stats_Y1[,"True"] == MI_stats_Y1[,"Mean"]
table(MI_stats_Y1[,"Coverage"])

MA_stats_Y2[,"Coverage"] = MA_stats_Y2[,"True"] == MA_stats_Y2[,"Mean"]
table(MA_stats_Y2[,"Coverage"])

MI_stats_Y2[,"Coverage"] = MI_stats_Y2[,"True"] == MI_stats_Y2[,"Mean"]
table(MI_stats_Y2[,"Coverage"])


# 1 to 15. 16 to 30, 31 to 45, 46 to 60, 61 to 75, and 76 to 90 
t_ini <- 1
t_fim <- 20

MA_data_1 = subset(x = MA_stats_Y1,
                   subset = Time >= t_ini & Time <= t_fim)
MA_data_2 = subset(x = MA_stats_Y2,
                   subset = Time >= t_ini & Time <= t_fim)
MI_data_1 = subset(x = MI_stats_Y1,
                   subset = Time >= t_ini & Time <= t_fim)
MI_data_2 = subset(x = MI_stats_Y2,
                   subset = Time >= t_ini & Time <= t_fim)
MA_resp1 = ggplot(data = MA_data_1,
                  aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(title = expression(paste(italic(M)[A])),
       y = "Values") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  ylim(min(MA_data_1[,2:5], MI_data_1[,2:5]),
       max(MA_data_1[,2:5], MI_data_1[,2:5])) +
  scale_x_continuous(breaks = seq(t_ini, t_fim, 2))
MI_resp1 = ggplot(data = MI_data_1,
                  aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(title = expression(paste(italic(M)[I])),
       y = "Values") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  ylim(min(MA_data_1[,2:5], MI_data_1[,2:5]),
       max(MA_data_1[,2:5], MI_data_1[,2:5])) +
  scale_x_continuous(breaks = seq(t_ini, t_fim, 2))
MA_resp2 = ggplot(data = MA_data_2,
                  aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(title = expression(paste(italic(M)[A])),
       y = "Values") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  ylim(min(MA_data_2[,2:5], MI_data_2[,2:5]),
       max(MA_data_2[,2:5], MI_data_2[,2:5])) +
  scale_x_continuous(breaks = seq(t_ini, t_fim, 2))
MI_resp2 = ggplot(data = MI_data_2,
                  aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(title = expression(paste(italic(M)[I])),
       y = "Values") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  ylim(min(MA_data_2[,2:5], MI_data_2[,2:5]),
       max(MA_data_2[,2:5], MI_data_2[,2:5])) +
  scale_x_continuous(breaks = seq(t_ini, t_fim, 2))

Fig_resp1 = ggplot(data = MI_data_1,
                   aes(x = Time, y = Mean)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("blue", "transparent")) +
  geom_errorbar(aes(ymin = LI, ymax = LS, color = Coverage),
                 linetype = 2) +
  geom_line(aes(y = True), linetype = 1, color = "black",
            linewidth = 0.75) +
  labs(title = "Resposta 1",
       y = "Values") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position = "none") +
  ylim(min(MA_data_1[,2:5], MI_data_1[,2:5]),
       max(MA_data_1[,2:5], MI_data_1[,2:5])) +
  scale_x_continuous(breaks = seq(t_ini - 1, t_fim, 10))

Fig_resp2 = ggplot(data = MI_data_2,
                   aes(x = Time, y = Mean)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("blue", "transparent")) +
  geom_errorbar(aes(ymin = LI, ymax = LS, color = Coverage),
                linetype = 2) +
  geom_line(aes(y = True), linetype = 1, color = "black",
            linewidth = 0.75) +
  labs(title = "Resposta 2",
       y = "Values") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position = "none") +
  ylim(min(MA_data_2[,2:5], MI_data_2[,2:5]),
       max(MA_data_2[,2:5], MI_data_2[,2:5])) +
  scale_x_continuous(breaks = seq(t_ini - 1, t_fim, 10))

# Landscape 4 x 6 in.
ggarrange(Fig_resp1, Fig_resp2,
          ncol = 1, nrow = 2, legend = "none")
ggarrange(MA_resp2, MI_resp2,
          ncol = 1, nrow = 2, legend = "none")
