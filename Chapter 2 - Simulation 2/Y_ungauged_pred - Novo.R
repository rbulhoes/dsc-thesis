##### Inform the directory 
setwd("D:/Exportações/Chap2/Sim2")
getwd()

##### Some packages
library(ggplot2)
library(ggpubr)
# library(grid)
# library(gridExtra)
library(ggrepel) # Para colocar r?tulos/etiquetas/nomes das esta??es nos mapas

smp_size <- 1000
N <- 16
N_i <- 3
N_tot <- N + N_i
K <- 100
K_p <- 10
K_tot <- K + K_p

##### Time series for gauged sites
Y = read.table("MA/Y.txt",
               head = FALSE,
               sep = ",")
colnames(Y) = c("Time", "Site", "Resp.1", "Resp.2")
Y_i = read.table("MA/Y_i.txt",
                 head = FALSE,
                 sep = ",")
colnames(Y_i) = c("Time", "Site", "Resp.1", "Resp.2")
Y_p_i = read.table("MA/Y_p_i.txt",
                   head = FALSE,
                   sep = ",")
colnames(Y_p_i) = c("Time", "Site", "Resp.1", "Resp.2")
S = read.table("MA/S.txt",
               head = FALSE,
               sep = ",")
colnames(S) = c("Site", "Longitude", "Latitude")
S_i = read.table("MA/S_i.txt",
                 head = FALSE,
                 sep = ",")
colnames(S_i) = c("Site", "Longitude", "Latitude")
Site = 1:N_tot

##### Time series for ungauged sites
e_Y_i_MA = read.table("MA/e_Y_i.txt",
                      head = FALSE,
                      sep = ",")
colnames(e_Y_i_MA) = c("Iteration", "Time", "Site", "Resp.1", "Resp.2")
e_Y_i_MI = read.table("MI/e_Y_i.txt",
                      head = FALSE,
                      sep = ",")
colnames(e_Y_i_MI) = c("Iteration", "Time", "Site", "Resp.1", "Resp.2")
e_Y_p_i_MI = read.table("MI/e_Y_p_i.txt",
                        head = FALSE,
                        sep = ",")
colnames(e_Y_p_i_MI) = c("Iteration", "Time", "Site", "Resp.1", "Resp.2")
e_Y_p_i_MA = read.table("MA/e_Y_p_i.txt",
                        head = FALSE,
                        sep = ",")
colnames(e_Y_p_i_MA) = c("Iteration", "Time", "Site", "Resp.1", "Resp.2")

n_i <- 2

MA_stats_Y1_i = data.frame(1:(K+K_p), matrix(data = 0, nrow = K+K_p, ncol = 5))
MA_stats_Y2_i = data.frame(1:(K+K_p), matrix(data = 0, nrow = K+K_p, ncol = 5))
colnames(MA_stats_Y1_i) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(MA_stats_Y2_i) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
MI_stats_Y1_i = data.frame(1:(K+K_p), matrix(data = 0, nrow = K+K_p, ncol = 5))
MI_stats_Y2_i = data.frame(1:(K+K_p), matrix(data = 0, nrow = K+K_p, ncol = 5))
colnames(MI_stats_Y1_i) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(MI_stats_Y2_i) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
for(t in 1:(K+K_p)) {
  if(t <= K) {
    set_MA = subset(x = e_Y_i_MA,
                    select = c("Iteration", "Resp.1", "Resp.2"),
                    subset = Site == n_i & Time == t)
    MA_stats_Y1_i[t,"LI"] = quantile(x = set_MA[1:smp_size, "Resp.1"],
                                     probs = 0.025)
    MA_stats_Y2_i[t,"LI"] = quantile(x = set_MA[1:smp_size, "Resp.2"],
                                     probs = 0.025)
    MA_stats_Y1_i[t,"Mean"] = mean(set_MA[1:smp_size, "Resp.1"])
    MA_stats_Y2_i[t,"Mean"] = mean(set_MA[1:smp_size, "Resp.2"])
    MA_stats_Y1_i[t,"True"] = subset(x = Y_i,
                                     select = c("Resp.1"),
                                     subset = Site == n_i & Time == t)
    MA_stats_Y2_i[t,"True"] = subset(x = Y_i,
                                     select = c("Resp.2"),
                                     subset = Site == n_i & Time == t)
    MA_stats_Y1_i[t,"LS"] = quantile(x = set_MA[1:smp_size, "Resp.1"],
                                     probs = 0.975)
    MA_stats_Y2_i[t,"LS"] = quantile(x = set_MA[1:smp_size, "Resp.2"],
                                     probs = 0.975)
    set_MI = subset(x = e_Y_i_MI,
                    select = c("Iteration", "Resp.1", "Resp.2"),
                    subset = Site == n_i & Time == t)
    MI_stats_Y1_i[t,"LI"] = quantile(x = set_MI[1:smp_size, "Resp.1"],
                                     probs = 0.025)
    MI_stats_Y2_i[t,"LI"] = quantile(x = set_MI[1:smp_size, "Resp.2"],
                                     probs = 0.025)
    MI_stats_Y1_i[t,"Mean"] = mean(set_MI[1:smp_size, "Resp.1"])
    MI_stats_Y2_i[t,"Mean"] = mean(set_MI[1:smp_size, "Resp.2"])
    MI_stats_Y1_i[t,"True"] = subset(x = Y_i,
                                     select = c("Resp.1"),
                                     subset = Site == n_i & Time == t)
    MI_stats_Y2_i[t,"True"] = subset(x = Y_i,
                                     select = c("Resp.2"),
                                     subset = Site == n_i & Time == t)
    MI_stats_Y1_i[t,"LS"] = quantile(x = set_MI[1:smp_size, "Resp.1"],
                                     probs = 0.975)
    MI_stats_Y2_i[t,"LS"] = quantile(x = set_MI[1:smp_size, "Resp.2"],
                                     probs = 0.975) 
  } else {
    set_MA = subset(x = e_Y_p_i_MA,
                    select = c("Iteration", "Resp.1", "Resp.2"),
                    subset = Site == n_i & Time == t - K)
    MA_stats_Y1_i[t,"LI"] = quantile(x = set_MA[1:smp_size, "Resp.1"],
                                     probs = 0.025)
    MA_stats_Y2_i[t,"LI"] = quantile(x = set_MA[1:smp_size, "Resp.2"],
                                     probs = 0.025)
    MA_stats_Y1_i[t,"Mean"] = mean(set_MA[1:smp_size, "Resp.1"])
    MA_stats_Y2_i[t,"Mean"] = mean(set_MA[1:smp_size, "Resp.2"])
    MA_stats_Y1_i[t,"True"] = subset(x = Y_p_i,
                                     select = c("Resp.1"),
                                     subset = Site == n_i & Time == t - K)
    MA_stats_Y2_i[t,"True"] = subset(x = Y_p_i,
                                     select = c("Resp.2"),
                                     subset = Site == n_i & Time == t - K)
    MA_stats_Y1_i[t,"LS"] = quantile(x = set_MA[1:smp_size, "Resp.1"],
                                     probs = 0.975)
    MA_stats_Y2_i[t,"LS"] = quantile(x = set_MA[1:smp_size, "Resp.2"],
                                     probs = 0.975)
    set_MI = subset(x = e_Y_p_i_MI,
                    select = c("Iteration", "Resp.1", "Resp.2"),
                    subset = Site == n_i & Time == t - K)
    MI_stats_Y1_i[t,"LI"] = quantile(x = set_MI[1:smp_size, "Resp.1"],
                                     probs = 0.025)
    MI_stats_Y2_i[t,"LI"] = quantile(x = set_MI[1:smp_size, "Resp.2"],
                                     probs = 0.025)
    MI_stats_Y1_i[t,"Mean"] = mean(set_MI[1:smp_size, "Resp.1"])
    MI_stats_Y2_i[t,"Mean"] = mean(set_MI[1:smp_size, "Resp.2"])
    MI_stats_Y1_i[t,"True"] = subset(x = Y_p_i,
                                     select = c("Resp.1"),
                                     subset = Site == n_i & Time == t - K)
    MI_stats_Y2_i[t,"True"] = subset(x = Y_p_i,
                                     select = c("Resp.2"),
                                     subset = Site == n_i & Time == t - K)
    MI_stats_Y1_i[t,"LS"] = quantile(x = set_MI[1:smp_size, "Resp.1"],
                                     probs = 0.975)
    MI_stats_Y2_i[t,"LS"] = quantile(x = set_MI[1:smp_size, "Resp.2"],
                                     probs = 0.975)     
  }
}

##### Coverage

MA_stats_Y1_i[,"Coverage"] <- (MA_stats_Y1_i[,"True"] < MA_stats_Y1_i[,"LI"]) |
  (MA_stats_Y1_i[,"True"] > MA_stats_Y1_i[,"LS"])
table(MA_stats_Y1_i[,"Coverage"])

MI_stats_Y1_i[,"Coverage"] <- (MI_stats_Y1_i[,"True"] < MI_stats_Y1_i[,"LI"]) |
  (MI_stats_Y1_i[,"True"] > MI_stats_Y1_i[,"LS"])
table(MI_stats_Y1_i[,"Coverage"])

MA_stats_Y2_i[,"Coverage"] <- (MA_stats_Y2_i[,"True"] < MA_stats_Y2_i[,"LI"]) |
  (MA_stats_Y2_i[,"True"] > MA_stats_Y2_i[,"LS"])
table(MA_stats_Y2_i[,"Coverage"])

MI_stats_Y2_i[,"Coverage"] <- (MI_stats_Y2_i[,"True"] < MI_stats_Y2_i[,"LI"]) |
  (MI_stats_Y2_i[,"True"] > MI_stats_Y2_i[,"LS"])
table(MI_stats_Y2_i[,"Coverage"])

# 1 to 55, and 56 to 110 
t_ini <- 56
t_fim <- 110

MA_data_1 <- subset(x = MA_stats_Y1_i,
                    subset = Time >= t_ini & Time <= t_fim)
MA_data_2 <- subset(x = MA_stats_Y2_i,
                    subset = Time >= t_ini & Time <= t_fim)
MI_data_1 <- subset(x = MI_stats_Y1_i,
                    subset = Time >= t_ini & Time <= t_fim)
MI_data_2 <- subset(x = MI_stats_Y2_i,
                    subset = Time >= t_ini & Time <= t_fim)
MA_resp1 <- ggplot(data = MA_data_1, aes(x = Time, y = True)) +
  geom_ribbon(aes(ymin = LI, ymax = LS), alpha = 0.30, fill = "lightgreen") +
  geom_line(aes(y = Mean), color = "gold", size = 0.45) +
  geom_point(size = 0.65, aes(color = Coverage)) +
  scale_color_manual(breaks = c(TRUE, FALSE), values = c("red", "blue")) +
  labs(title = expression(paste(italic(M)[A])),
       y = "Values") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  ylim(min(MA_data_1[,2:5], MI_data_1[,2:5]),
       max(MA_data_1[,2:5], MI_data_1[,2:5])) +
  scale_x_continuous(breaks = seq(t_ini, t_fim, 10))
MI_resp1 <- ggplot(data = MI_data_1, aes(x = Time, y = True)) +
  geom_ribbon(aes(ymin = LI, ymax = LS), alpha = 0.30, fill = "lightgreen") +
  geom_line(aes(y = Mean), color = "gold", size = 0.45) +
  geom_point(size = 0.65, aes(color = Coverage)) +
  scale_color_manual(breaks = c(TRUE, FALSE), values = c("red", "blue")) +
  labs(title = expression(paste(italic(M)[I])),
       y = "Values") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  ylim(min(MA_data_1[,2:5], MI_data_1[,2:5]),
       max(MA_data_1[,2:5], MI_data_1[,2:5])) +
  scale_x_continuous(breaks = seq(t_ini, t_fim, 10))
MA_resp2 <- ggplot(data = MA_data_2, aes(x = Time, y = True)) +
  geom_ribbon(aes(ymin = LI, ymax = LS), alpha = 0.30, fill = "lightgreen") +
  geom_line(aes(y = Mean), color = "gold", size = 0.45) +
  geom_point(size = 0.65, aes(color = Coverage)) +
  scale_color_manual(breaks = c(TRUE, FALSE), values = c("red", "blue")) +
  labs(title = expression(paste(italic(M)[A])),
       y = "Values") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  ylim(min(MA_data_2[,2:5], MI_data_2[,2:5]),
       max(MA_data_2[,2:5], MI_data_2[,2:5])) +
  scale_x_continuous(breaks = seq(t_ini, t_fim, 10))
MI_resp2 <- ggplot(data = MI_data_2, aes(x = Time, y = True)) +
  geom_ribbon(aes(ymin = LI, ymax = LS), alpha = 0.30, fill = "lightgreen") +
  geom_line(aes(y = Mean), color = "gold", size = 0.45) +
  geom_point(size = 0.65, aes(color = Coverage)) +
  scale_color_manual(breaks = c(TRUE, FALSE), values = c("red", "blue")) +
  labs(title = expression(paste(italic(M)[I])),
       y = "Values") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  ylim(min(MA_data_2[,2:5], MI_data_2[,2:5]),
       max(MA_data_2[,2:5], MI_data_2[,2:5])) +
  scale_x_continuous(breaks = seq(t_ini, t_fim, 10))

# Landscape US Legal
ggarrange(MA_resp1, MI_resp1,
          ncol = 1, nrow = 2, legend = "none")
ggarrange(MA_resp2, MI_resp2,
          ncol = 1, nrow = 2, legend = "none")
