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
K <- 100 # Number of times
K_p <- 10
K_tot <- K + K_p

beta_AM <- read.table("MA/beta.txt",
                      head = FALSE,
                      sep = ",")
colnames(beta_AM) <- c("Time", "beta00", "beta01", "beta10", "beta11")
beta_IM <- read.table("MI/beta.txt",
                      head = FALSE,
                      sep = ",")
colnames(beta_IM) <- c("Time", "beta00", "beta01", "beta10", "beta11")
e_beta_AM <- read.table("MA/e_beta.txt",
                        head = FALSE,
                        sep = ",")
colnames(e_beta_AM) <- c("Iteration", "Time",
                         "beta00", "beta01", "beta10", "beta11")
e_beta_IM <- read.table("MI/e_beta.txt",
                        head = FALSE,
                        sep = ",")
colnames(e_beta_IM) <- c("Iteration", "Time",
                         "beta00", "beta01", "beta10", "beta11")
beta_p_AM <- read.table("MA/beta_p.txt",
                        head = FALSE,
                        sep = ",")
colnames(beta_p_AM) <- c("Time", "beta00", "beta01", "beta10", "beta11")
beta_p_IM <- read.table("MI/beta_p.txt",
                        head = FALSE,
                        sep = ",")
colnames(beta_p_IM) <- c("Time", "beta00", "beta01", "beta10", "beta11")
e_beta_p_AM <- read.table("MA/e_beta_p.txt",
                          head = FALSE,
                          sep = ",")
colnames(e_beta_p_AM) <- c("Iteration", "Time",
                           "beta00", "beta01", "beta10", "beta11")
e_beta_p_IM <- read.table("MI/e_beta_p.txt",
                          head = FALSE,
                          sep = ",")
colnames(e_beta_p_IM) <- c("Iteration", "Time",
                           "beta00", "beta01", "beta10", "beta11")

MA_stats_b00 <- data.frame(0:K_tot, matrix(data = 0, nrow = K_tot+1, ncol = 5))
MA_stats_b01 <- data.frame(0:K_tot, matrix(data = 0, nrow = K_tot+1, ncol = 5))
MA_stats_b10 <- data.frame(0:K_tot, matrix(data = 0, nrow = K_tot+1, ncol = 5))
MA_stats_b11 <- data.frame(0:K_tot, matrix(data = 0, nrow = K_tot+1, ncol = 5))
colnames(MA_stats_b00) <- c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(MA_stats_b01) <- c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(MA_stats_b10) <- c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(MA_stats_b11) <- c("Time", "LI", "Mean", "True", "LS", "Coverage")
MI_stats_b00 <- data.frame(0:K_tot, matrix(data = 0, nrow = K_tot+1, ncol = 5))
MI_stats_b01 <- data.frame(0:K_tot, matrix(data = 0, nrow = K_tot+1, ncol = 5))
MI_stats_b10 <- data.frame(0:K_tot, matrix(data = 0, nrow = K_tot+1, ncol = 5))
MI_stats_b11 <- data.frame(0:K_tot, matrix(data = 0, nrow = K_tot+1, ncol = 5))
colnames(MI_stats_b00) <- c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(MI_stats_b01) <- c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(MI_stats_b10) <- c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(MI_stats_b11) <- c("Time", "LI", "Mean", "True", "LS", "Coverage")
for(t in 0:K_tot) {
  if(t <= K) {
    set_MA = subset(x = e_beta_AM,
                    subset = Time == t,
                    select = c("Iteration", "beta00", "beta01",
                               "beta10", "beta11"))
    MA_stats_b00[t+1,"LI"] = quantile(x = set_MA[1:smp_size, "beta00"],
                                      probs = 0.025)
    MA_stats_b01[t+1,"LI"] = quantile(x = set_MA[1:smp_size, "beta01"],
                                      probs = 0.025)
    MA_stats_b10[t+1,"LI"] = quantile(x = set_MA[1:smp_size, "beta10"],
                                      probs = 0.025)
    MA_stats_b11[t+1,"LI"] = quantile(x = set_MA[1:smp_size, "beta11"],
                                      probs = 0.025)
    MA_stats_b00[t+1,"Mean"] = mean(set_MA[1:smp_size, "beta00"])
    MA_stats_b01[t+1,"Mean"] = mean(set_MA[1:smp_size, "beta01"])
    MA_stats_b10[t+1,"Mean"] = mean(set_MA[1:smp_size, "beta10"])
    MA_stats_b11[t+1,"Mean"] = mean(set_MA[1:smp_size, "beta11"])
    MA_stats_b00[t+1,"True"] = subset(x = beta_AM,
                                      subset = Time == t,
                                      select = c("beta00"))
    MA_stats_b01[t+1,"True"] = subset(x = beta_AM,
                                      subset = Time == t,
                                      select = c("beta01"))
    MA_stats_b10[t+1,"True"] = subset(x = beta_AM,
                                      subset = Time == t,
                                      select = c("beta10"))
    MA_stats_b11[t+1,"True"] = subset(x = beta_AM,
                                      subset = Time == t,
                                      select = c("beta11"))
    MA_stats_b00[t+1,"LS"] = quantile(x = set_MA[1:smp_size, "beta00"],
                                      probs = 0.975)
    MA_stats_b01[t+1,"LS"] = quantile(x = set_MA[1:smp_size, "beta01"],
                                      probs = 0.975)
    MA_stats_b10[t+1,"LS"] = quantile(x = set_MA[1:smp_size, "beta10"],
                                      probs = 0.975)
    MA_stats_b11[t+1,"LS"] = quantile(x = set_MA[1:smp_size, "beta11"],
                                      probs = 0.975)
  } else {
    set_MA = subset(x = e_beta_p_AM,
                    subset = Time == t - K,
                    select = c("Iteration", "beta00", "beta01",
                               "beta10", "beta11"))
    MA_stats_b00[t+1,"LI"] = quantile(x = set_MA[1:smp_size, "beta00"],
                                      probs = 0.025)
    MA_stats_b01[t+1,"LI"] = quantile(x = set_MA[1:smp_size, "beta01"],
                                      probs = 0.025)
    MA_stats_b10[t+1,"LI"] = quantile(x = set_MA[1:smp_size, "beta10"],
                                      probs = 0.025)
    MA_stats_b11[t+1,"LI"] = quantile(x = set_MA[1:smp_size, "beta11"],
                                      probs = 0.025)
    MA_stats_b00[t+1,"Mean"] = mean(set_MA[1:smp_size, "beta00"])
    MA_stats_b01[t+1,"Mean"] = mean(set_MA[1:smp_size, "beta01"])
    MA_stats_b10[t+1,"Mean"] = mean(set_MA[1:smp_size, "beta10"])
    MA_stats_b11[t+1,"Mean"] = mean(set_MA[1:smp_size, "beta11"])
    MA_stats_b00[t+1,"True"] = subset(x = beta_p_AM,
                                      subset = Time == t - K,
                                      select = c("beta00"))
    MA_stats_b01[t+1,"True"] = subset(x = beta_p_AM,
                                      subset = Time == t - K,
                                      select = c("beta01"))
    MA_stats_b10[t+1,"True"] = subset(x = beta_p_AM,
                                      subset = Time == t - K,
                                      select = c("beta10"))
    MA_stats_b11[t+1,"True"] = subset(x = beta_p_AM,
                                      subset = Time == t - K,
                                      select = c("beta11"))
    MA_stats_b00[t+1,"LS"] = quantile(x = set_MA[1:smp_size, "beta00"],
                                      probs = 0.975)
    MA_stats_b01[t+1,"LS"] = quantile(x = set_MA[1:smp_size, "beta01"],
                                      probs = 0.975)
    MA_stats_b10[t+1,"LS"] = quantile(x = set_MA[1:smp_size, "beta10"],
                                      probs = 0.975)
    MA_stats_b11[t+1,"LS"] = quantile(x = set_MA[1:smp_size, "beta11"],
                                      probs = 0.975)
  }
}
for(t in 0:K_tot) {
  if(t <= K) {
    set_MI = subset(x = e_beta_AM,
                    subset = Time == t,
                    select = c("Iteration", "beta00", "beta01",
                               "beta10", "beta11"))
    MI_stats_b00[t+1,"LI"] = quantile(x = set_MI[1:smp_size, "beta00"],
                                      probs = 0.025)
    MI_stats_b01[t+1,"LI"] = quantile(x = set_MI[1:smp_size, "beta01"],
                                      probs = 0.025)
    MI_stats_b10[t+1,"LI"] = quantile(x = set_MI[1:smp_size, "beta10"],
                                      probs = 0.025)
    MI_stats_b11[t+1,"LI"] = quantile(x = set_MI[1:smp_size, "beta11"],
                                      probs = 0.025)
    MI_stats_b00[t+1,"Mean"] = mean(set_MI[1:smp_size, "beta00"])
    MI_stats_b01[t+1,"Mean"] = mean(set_MI[1:smp_size, "beta01"])
    MI_stats_b10[t+1,"Mean"] = mean(set_MI[1:smp_size, "beta10"])
    MI_stats_b11[t+1,"Mean"] = mean(set_MI[1:smp_size, "beta11"])
    MI_stats_b00[t+1,"True"] = subset(x = beta_AM,
                                      subset = Time == t,
                                      select = c("beta00"))
    MI_stats_b01[t+1,"True"] = subset(x = beta_AM,
                                      subset = Time == t,
                                      select = c("beta01"))
    MI_stats_b10[t+1,"True"] = subset(x = beta_AM,
                                      subset = Time == t,
                                      select = c("beta10"))
    MI_stats_b11[t+1,"True"] = subset(x = beta_AM,
                                      subset = Time == t,
                                      select = c("beta11"))
    MI_stats_b00[t+1,"LS"] = quantile(x = set_MI[1:smp_size, "beta00"],
                                      probs = 0.975)
    MI_stats_b01[t+1,"LS"] = quantile(x = set_MI[1:smp_size, "beta01"],
                                      probs = 0.975)
    MI_stats_b10[t+1,"LS"] = quantile(x = set_MI[1:smp_size, "beta10"],
                                      probs = 0.975)
    MI_stats_b11[t+1,"LS"] = quantile(x = set_MI[1:smp_size, "beta11"],
                                      probs = 0.975)
  } else {
    set_MI = subset(x = e_beta_p_AM,
                    subset = Time == t - K,
                    select = c("Iteration", "beta00", "beta01",
                               "beta10", "beta11"))
    MI_stats_b00[t+1,"LI"] = quantile(x = set_MI[1:smp_size, "beta00"],
                                      probs = 0.025)
    MI_stats_b01[t+1,"LI"] = quantile(x = set_MI[1:smp_size, "beta01"],
                                      probs = 0.025)
    MI_stats_b10[t+1,"LI"] = quantile(x = set_MI[1:smp_size, "beta10"],
                                      probs = 0.025)
    MI_stats_b11[t+1,"LI"] = quantile(x = set_MI[1:smp_size, "beta11"],
                                      probs = 0.025)
    MI_stats_b00[t+1,"Mean"] = mean(set_MI[1:smp_size, "beta00"])
    MI_stats_b01[t+1,"Mean"] = mean(set_MI[1:smp_size, "beta01"])
    MI_stats_b10[t+1,"Mean"] = mean(set_MI[1:smp_size, "beta10"])
    MI_stats_b11[t+1,"Mean"] = mean(set_MI[1:smp_size, "beta11"])
    MI_stats_b00[t+1,"True"] = subset(x = beta_p_AM,
                                      subset = Time == t - K,
                                      select = c("beta00"))
    MI_stats_b01[t+1,"True"] = subset(x = beta_p_AM,
                                      subset = Time == t - K,
                                      select = c("beta01"))
    MI_stats_b10[t+1,"True"] = subset(x = beta_p_AM,
                                      subset = Time == t - K,
                                      select = c("beta10"))
    MI_stats_b11[t+1,"True"] = subset(x = beta_p_AM,
                                      subset = Time == t - K,
                                      select = c("beta11"))
    MI_stats_b00[t+1,"LS"] = quantile(x = set_MI[1:smp_size, "beta00"],
                                      probs = 0.975)
    MI_stats_b01[t+1,"LS"] = quantile(x = set_MI[1:smp_size, "beta01"],
                                      probs = 0.975)
    MI_stats_b10[t+1,"LS"] = quantile(x = set_MI[1:smp_size, "beta10"],
                                      probs = 0.975)
    MI_stats_b11[t+1,"LS"] = quantile(x = set_MI[1:smp_size, "beta11"],
                                      probs = 0.975)
  }
  
}

MA_stats_b00[,"Coverage"]=MA_stats_b00[,"True"] >= MA_stats_b00[,"LI"] &
                            MA_stats_b00[,"True"] <= MA_stats_b00[,"LS"]
table(MA_stats_b00[,"Coverage"])

MI_stats_b00[,"Coverage"]=MI_stats_b00[,"True"] >= MI_stats_b00[,"LI"] &
                            MI_stats_b00[,"True"] <= MI_stats_b00[,"LS"]
table(MI_stats_b00[,"Coverage"])

MA_stats_b01[,"Coverage"]=MA_stats_b01[,"True"] >= MA_stats_b01[,"LI"] &
  MA_stats_b01[,"True"] <= MA_stats_b01[,"LS"]
table(MA_stats_b01[,"Coverage"])

MI_stats_b01[,"Coverage"]=MI_stats_b01[,"True"] >= MI_stats_b01[,"LI"] &
  MI_stats_b01[,"True"] <= MI_stats_b01[,"LS"]
table(MI_stats_b01[,"Coverage"])

MA_stats_b10[,"Coverage"]=MA_stats_b10[,"True"] >= MA_stats_b10[,"LI"] &
  MA_stats_b10[,"True"] <= MA_stats_b10[,"LS"]
table(MA_stats_b10[,"Coverage"])

MI_stats_b10[,"Coverage"]=MI_stats_b10[,"True"] >= MI_stats_b10[,"LI"] &
  MI_stats_b10[,"True"] <= MI_stats_b10[,"LS"]
table(MI_stats_b10[,"Coverage"])

MA_stats_b11[,"Coverage"]=MA_stats_b11[,"True"] >= MA_stats_b11[,"LI"] &
  MA_stats_b11[,"True"] <= MA_stats_b11[,"LS"]
table(MA_stats_b11[,"Coverage"])

MI_stats_b11[,"Coverage"]=MI_stats_b11[,"True"] >= MI_stats_b11[,"LI"] &
  MI_stats_b11[,"True"] <= MI_stats_b11[,"LS"]
table(MI_stats_b11[,"Coverage"])

plot_MA_beta00 = ggplot(data = MA_stats_b00,
                          aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["0,1,"][italic(t)])) +
  labs(title = expression(italic(M)[A])) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_MI_beta00 = ggplot(data = MI_stats_b00,
                          aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["0,1,"][italic(t)])) +
  labs(title = expression(italic(M)[I])) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_MA_beta01 = ggplot(data = MA_stats_b01,
                          aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["0,2,"][italic(t)])) +
  labs(title = expression(italic(M)[A])) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_MI_beta01 = ggplot(data = MI_stats_b01,
                          aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["0,2,"][italic(t)])) +
  labs(title = expression(italic(M)[I])) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_MA_beta10 = ggplot(data = MA_stats_b10,
                          aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["1,1,"][italic(t)])) +
  labs(title = expression(italic(M)[A])) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_MI_beta10 = ggplot(data = MI_stats_b10,
                          aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["1,1,"][italic(t)])) +
  labs(title = expression(italic(M)[I])) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_MA_beta11 = ggplot(data = MA_stats_b11,
                          aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["1,2,"][italic(t)])) +
  labs(title = expression(italic(M)[A])) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_MI_beta11 = ggplot(data = MI_stats_b11,
                          aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["1,2,"][italic(t)])) +
  labs(title = expression(italic(M)[I])) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

# PDF in US Legal - Landscape
ggarrange(plot_MA_beta00, plot_MI_beta00,
          ncol = 1, nrow = 2, legend = "none")

ggarrange(plot_MA_beta01, plot_MI_beta01,
          ncol = 1, nrow = 2, legend = "none")

ggarrange(plot_MA_beta10, plot_MI_beta10,
          ncol = 1, nrow = 2, legend = "none")

ggarrange(plot_MA_beta11, plot_MI_beta11,
          ncol = 1, nrow = 2, legend = "none")

