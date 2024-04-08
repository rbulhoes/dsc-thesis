##### Inform the directory 
setwd("D:/Exportações/Chap3/Sim2")
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
Y = read.table("MI/Y.txt",
               head = FALSE,
               sep = ",")
colnames(Y) = c("Time", "Site", "Resp.1", "Resp.2")
e_Y_MA = read.table("MA/e_Y.txt",
                    head = FALSE,
                    sep = ",")
colnames(e_Y_MA) = c("Iteration", "Time", "Site", "Resp.1", "Resp.2")
e_Y_MI = read.table("MI/e_Y.txt",
                    head = FALSE,
                    sep = ",")
colnames(e_Y_MI) = c("Iteration", "Time", "Site", "Resp.1", "Resp.2")

##### Predicted time series for gauged sites
Y_p = read.table("MA/Y_p.txt",
                 head = FALSE,
                 sep = ",")
colnames(Y_p) = c("Time", "Site", "Resp.1", "Resp.2")
e_Y_p_MA = read.table("MA/e_Y_p.txt",
                      head = FALSE,
                      sep = ",")
colnames(e_Y_p_MA) = c("Iteration", "Time", "Site", "Resp.1", "Resp.2")
e_Y_p_MI = read.table("MI/e_Y_p.txt",
                      head = FALSE,
                      sep = ",")
colnames(e_Y_p_MI) = c("Iteration", "Time", "Site", "Resp.1", "Resp.2")

n <- 14 # 14 (explored in my thesis), ? (histogram)

MA_stats_Y1 = data.frame(1:(K+K_p), matrix(data = 0, nrow = K+K_p, ncol = 5))
MA_stats_Y2 = data.frame(1:(K+K_p), matrix(data = 0, nrow = K+K_p, ncol = 5))
colnames(MA_stats_Y1) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(MA_stats_Y2) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
MI_stats_Y1 = data.frame(1:(K+K_p), matrix(data = 0, nrow = K+K_p, ncol = 5))
MI_stats_Y2 = data.frame(1:(K+K_p), matrix(data = 0, nrow = K+K_p, ncol = 5))
colnames(MI_stats_Y1) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(MI_stats_Y2) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
for(t in 1:(K+K_p)) {
  if(t <= K) {
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
  } else {
    set_MA = subset(x = e_Y_p_MA,
                    select = c("Iteration", "Resp.1", "Resp.2"),
                    subset = Site == n & Time == t - K)
    MA_stats_Y1[t,"LI"] = quantile(x = set_MA[1:smp_size, "Resp.1"],
                                   probs = 0.025)
    MA_stats_Y2[t,"LI"] = quantile(x = set_MA[1:smp_size, "Resp.2"],
                                   probs = 0.025)
    MA_stats_Y1[t,"Mean"] = mean(set_MA[1:smp_size, "Resp.1"])
    MA_stats_Y2[t,"Mean"] = mean(set_MA[1:smp_size, "Resp.2"])
    MA_stats_Y1[t,"True"] = subset(x = Y_p,
                                   select = c("Resp.1"),
                                   subset = Site == n & Time == t - K)
    MA_stats_Y2[t,"True"] = subset(x = Y_p,
                                   select = c("Resp.2"),
                                   subset = Site == n & Time == t - K)
    MA_stats_Y1[t,"LS"] = quantile(x = set_MA[1:smp_size, "Resp.1"],
                                   probs = 0.975)
    MA_stats_Y2[t,"LS"] = quantile(x = set_MA[1:smp_size, "Resp.2"],
                                   probs = 0.975)
    set_MI = subset(x = e_Y_p_MI,
                    select = c("Iteration", "Resp.1", "Resp.2"),
                    subset = Site == n & Time == t - K)
    MI_stats_Y1[t,"LI"] = quantile(x = set_MI[1:smp_size, "Resp.1"],
                                   probs = 0.025)
    MI_stats_Y2[t,"LI"] = quantile(x = set_MI[1:smp_size, "Resp.2"],
                                   probs = 0.025)
    MI_stats_Y1[t,"Mean"] = mean(set_MI[1:smp_size, "Resp.1"])
    MI_stats_Y2[t,"Mean"] = mean(set_MI[1:smp_size, "Resp.2"])
    MI_stats_Y1[t,"True"] = subset(x = Y_p,
                                   select = c("Resp.1"),
                                   subset = Site == n & Time == t - K)
    MI_stats_Y2[t,"True"] = subset(x = Y_p,
                                   select = c("Resp.2"),
                                   subset = Site == n & Time == t - K)
    MI_stats_Y1[t,"LS"] = quantile(x = set_MI[1:smp_size, "Resp.1"],
                                   probs = 0.975)
    MI_stats_Y2[t,"LS"] = quantile(x = set_MI[1:smp_size, "Resp.2"],
                                   probs = 0.975)
  }
}


##### Coverage

MA_stats_Y1[,"Coverage"] <- (MA_stats_Y1[,"True"] < MA_stats_Y1[,"LI"]) |
  (MA_stats_Y1[,"True"] > MA_stats_Y1[,"LS"])
table(MA_stats_Y1[,"Coverage"])

MI_stats_Y1[,"Coverage"] <- (MI_stats_Y1[,"True"] < MI_stats_Y1[,"LI"]) |
  (MI_stats_Y1[,"True"] > MI_stats_Y1[,"LS"])
table(MI_stats_Y1[,"Coverage"])

MA_stats_Y2[,"Coverage"] <- (MA_stats_Y2[,"True"] < MA_stats_Y2[,"LI"]) |
  (MA_stats_Y2[,"True"] > MA_stats_Y2[,"LS"])
table(MA_stats_Y2[,"Coverage"])

MI_stats_Y2[,"Coverage"] <- (MI_stats_Y2[,"True"] < MI_stats_Y2[,"LI"]) |
  (MI_stats_Y2[,"True"] > MI_stats_Y2[,"LS"])
table(MI_stats_Y2[,"Coverage"])

# 0 to 36, 37 to 73, and 74 to 110 
t_ini <- 74
t_fim <- 110

MA_data_1 <- subset(x = MA_stats_Y1,
                    subset = Time >= t_ini & Time <= t_fim)
MA_data_2 <- subset(x = MA_stats_Y2,
                    subset = Time >= t_ini & Time <= t_fim)
MI_data_1 <- subset(x = MI_stats_Y1,
                    subset = Time >= t_ini & Time <= t_fim)
MI_data_2 <- subset(x = MI_stats_Y2,
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
  scale_x_continuous(breaks = seq(t_ini - 1, t_fim, 10))
MI_resp1 <- ggplot(data = MI_data_1, aes(x = Time, y = True)) +
  geom_ribbon(aes(ymin = LI, ymax = LS), alpha = 0.30, fill = "lightgreen") +
  geom_line(aes(y = Mean), color = "gold", size = 0.45) +
  geom_point(size = 0.65, aes(color = Coverage)) +
  scale_color_manual(breaks = c(TRUE, FALSE), values = c("red", "blue")) +
  labs(title = expression(paste(italic(M)["I"])),
       y = "Values") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  ylim(min(MA_data_1[,2:5], MI_data_1[,2:5]),
       max(MA_data_1[,2:5], MI_data_1[,2:5])) +
  scale_x_continuous(breaks = seq(t_ini - 1, t_fim, 10))
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
  scale_x_continuous(breaks = seq(t_ini - 1, t_fim, 10))
MI_resp2 <- ggplot(data = MI_data_2, aes(x = Time, y = True)) +
  geom_ribbon(aes(ymin = LI, ymax = LS), alpha = 0.30, fill = "lightgreen") +
  geom_line(aes(y = Mean), color = "gold", size = 0.45) +
  geom_point(size = 0.65, aes(color = Coverage)) +
  scale_color_manual(breaks = c(TRUE, FALSE), values = c("red", "blue")) +
  labs(title = expression(paste(italic(M)["I"])),
       y = "Values") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  ylim(min(MA_data_2[,2:5], MI_data_2[,2:5]),
       max(MA_data_2[,2:5], MI_data_2[,2:5])) +
  scale_x_continuous(breaks = seq(t_ini - 1, t_fim, 10))

# Landscape US Legal
ggarrange(MA_resp1, MI_resp1,
          ncol = 1, nrow = 2, legend = "none")
ggarrange(MA_resp2, MI_resp2,
          ncol = 1, nrow = 2, legend = "none")

## Comparing ranges

for (t in 1:K) {
  # range_Y1_MA <- MA_stats_Y1[t, "LS"] - MA_stats_Y1[t, "LI"]
  # range_Y1_MI <- MA_stats_Y1[t, "LS"] - MA_stats_Y1[t, "LI"]
  # test_Y1 <- range_Y1_MA > range_Y1_MI
  diff_Y1_MA <- abs(MA_stats_Y1[t, "Mean"] - MA_stats_Y1[t, "True"])
  diff_Y1_MI <- abs(MI_stats_Y1[t, "Mean"] - MI_stats_Y1[t, "True"])
  test_Y1 <- diff_Y1_MA < diff_Y1_MI
  # range_Y2_MA <- MA_stats_Y2[t, "LS"] - MA_stats_Y2[t, "LI"]
  # range_Y2_MI <- MA_stats_Y2[t, "LS"] - MA_stats_Y2[t, "LI"]
  # test_Y2 <- range_Y2_MA > range_Y2_MI
  diff_Y2_MA <- abs(MA_stats_Y2[t, "Mean"] - MA_stats_Y2[t, "True"])
  diff_Y2_MI <- abs(MI_stats_Y2[t, "Mean"] - MI_stats_Y2[t, "True"])
  test_Y2 <- diff_Y2_MA < diff_Y2_MI
  if (test_Y1 & test_Y2) {
    print(paste(t, "yes"))
  } else {
    print(paste(t, "no"))
  }
}


## Histograms for the observed time t <- t0, the gauged site n <- 14,
## and both response variables, by model (anisotropic and isotropic models)
t_Y1 <- 14
t_Y2 <- 2
# ggarrange(MA_resp1, MA_resp2, ncol = 1, nrow = 2, legend = "none")

set_hist_MA_Y1 <- subset(x = e_Y_MA,
                         subset = Time == t_Y1 & Site == n,
                         select = Resp.1)
set_hist_MI_Y1 <- subset(x = e_Y_MI,
                         subset = Time == t_Y1 & Site == n,
                         select = Resp.1)

minimo1 <- min(min(set_hist_MA_Y1), min(set_hist_MI_Y1))
maximo1 <- max(max(set_hist_MA_Y1), max(set_hist_MI_Y1))

hist_MA_Y1 <- ggplot(set_hist_MA_Y1, aes(x = Resp.1)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = Resp.1,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 15) +
  geom_vline(xintercept = MA_stats_Y1[t_Y1, "True"], 
             colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = MA_stats_Y1[t_Y1, "Mean"], 
             colour = "gold") +
  geom_vline(xintercept = MA_stats_Y1[t_Y1, "LI"],
             linetype="dashed") +
  geom_vline(xintercept = MA_stats_Y1[t_Y1, "LS"],
             linetype="dashed") +
  labs(x = "Response 1", y = "%") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  xlim(minimo1, maximo1) + ylim(0, 25)

hist_MI_Y1 <- ggplot(set_hist_MI_Y1, aes(x = Resp.1)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = Resp.1,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 15) +
  geom_vline(xintercept = MI_stats_Y1[t_Y1, "True"], 
             colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = MI_stats_Y1[t_Y1, "Mean"], 
             colour = "gold") +
  geom_vline(xintercept = MI_stats_Y1[t_Y1, "LI"],
             linetype="dashed") +
  geom_vline(xintercept = MI_stats_Y1[t_Y1, "LS"],
             linetype="dashed") +
  labs(x = "Response 1", y = "%") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  xlim(minimo1, maximo1) + ylim(0, 25)

set_hist_MA_Y2 <- subset(x = e_Y_MA,
                         subset = Time == t_Y2 & Site == n,
                         select = Resp.2)
set_hist_MI_Y2 <- subset(x = e_Y_MI,
                         subset = Time == t_Y2 & Site == n,
                         select = Resp.2)

minimo2 <- min(min(set_hist_MA_Y2), min(set_hist_MI_Y2))
maximo2 <- max(max(set_hist_MA_Y2), max(set_hist_MI_Y2))

hist_MA_Y2 <- ggplot(set_hist_MA_Y2, aes(x = Resp.2)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = Resp.2,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 15) +
  geom_vline(xintercept = MA_stats_Y2[t_Y2, "True"], 
             colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = MA_stats_Y2[t_Y2, "Mean"], 
             colour = "gold") +
  geom_vline(xintercept = MA_stats_Y2[t_Y2, "LI"],
             linetype="dashed") +
  geom_vline(xintercept = MA_stats_Y2[t_Y2, "LS"],
             linetype="dashed") +
  labs(x = "Response 2", y = "%") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  xlim(minimo2, maximo2) + ylim(0, 25)

hist_MI_Y2 <- ggplot(set_hist_MI_Y2, aes(x = Resp.2)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = Resp.2,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 15) +
  geom_vline(xintercept = MI_stats_Y2[t_Y2, "True"], 
             colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = MI_stats_Y2[t_Y2, "Mean"], 
             colour = "gold") +
  geom_vline(xintercept = MI_stats_Y2[t_Y2, "LI"],
             linetype="dashed") +
  geom_vline(xintercept = MI_stats_Y2[t_Y2, "LS"],
             linetype="dashed") +
  labs(x = "Response 2", y = "%") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  xlim(minimo2, maximo2) + ylim(0, 25)

# Salvar 4 x 6 in. - Landscape
ggarrange(hist_MA_Y1, hist_MI_Y1, ncol = 2, nrow = 1)
ggarrange(hist_MA_Y2, hist_MI_Y2, ncol = 2, nrow = 1)
