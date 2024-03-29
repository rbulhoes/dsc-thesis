##### Inform the directory 
setwd("D:/Exportações/Chap2/Sim2")
getwd()

##### Some packages
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(ggrepel) # Para colocar r?tulos/etiquetas/nomes das esta??es nos mapas
library(scoringutils)

smp_size <- 1000
N <- 16
N_i <- 3
N_tot <- N + N_i
K <- 100
K_p <- 10
K_tot <- K + K_p

Y <- read.table("MI/Y.txt",
                head = FALSE,
                sep = ",")
colnames(Y) <- c("Time", "Site", "Resp.1", "Resp.2")
Y_i <- read.table("MI/Y_i.txt",
                  head = FALSE,
                  sep = ",")
colnames(Y_i) <- c("Time", "Site", "Resp.1", "Resp.2")
S <- read.table("MI/S.txt",
                head = FALSE,
                sep = ",")
colnames(S) <- c("Site", "Longitude", "Latitude")
S_i <- read.table("MI/S_i.txt",
                  head = FALSE,
                  sep = ",")
colnames(S_i) <- c("Site", "Longitude", "Latitude")
Site <- 1:N_tot

##### Maps for each t
t <- 20
data_t <- data.frame(Site,
                     rbind(subset(x = S,
                                  select = c("Longitude", "Latitude")),
                           subset(x = S_i,
                                  select = c("Longitude", "Latitude"))),
                     rbind(subset(x = Y,
                                  select = c("Resp.1", "Resp.2"),
                                  subset = Time == t),
                           subset(x = Y_i,
                                  select = c("Resp.1", "Resp.2"),
                                  subset = Time == t)))

R1_t <- ggplot() +
  geom_point(data = data_t, size = 7, aes(x = Longitude,
                                          y = Latitude,
                                          color = Resp.1)) +
  scale_colour_gradient(low = "gold", high = "red", na.value = NA) +
  geom_text_repel(data = data_t,
                  aes(x = Longitude, y = Latitude, label = Site),
                  size = 3) +
  labs(title = paste("Response 1, Time =", t)) +
  theme(plot.title = element_text(hjust = 0.5))
R2_t <- ggplot() +
  geom_point(data = data_t, size = 7, aes(x = Longitude,
                                          y = Latitude,
                                          color = Resp.2)) +
  scale_colour_gradient(low = "gold", high = "red", na.value = NA) +
  geom_text_repel(data = data_t,
                  aes(x = Longitude, y = Latitude, label = Site),
                  size = 3) +
  labs(title = paste("Response 2, Time =", t)) + 
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(R1_t, R2_t, ncol = 2, nrow = 1)

##### Geographic region
library(ggrepel) # Para colocar r?tulos/etiquetas/nomes das esta??es nos mapas
Sites <- 1:N_tot
Category <- factor(x = c(rep(1, 2),
                         rep(2, N - 2),
                         rep(3, N_i)),
                   levels = 1:3,
                   labels = c("Anchor", "Non-anchor", "Interpolation"))
GR <- data.frame(Sites, Category, rbind(S[1:N, 2:3], S_i[1:N_i, 2:3]))
ggplot() + 
  geom_point(data = GR, size = 2, aes(x = Longitude,
                                      y = Latitude,
                                      color = Category)) +
  geom_text_repel(data = GR,
                  aes(x = Longitude, y = Latitude, label = Sites),
                  size = 3.5) +
  theme(legend.position = 'right',
        text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = 1)
# Landscape 5 x 7 in.



##### Time series for ungauged sites
e_Y_i_MA <- read.table("MA/e_Y_i.txt",
                       head = FALSE,
                       sep = ",")
colnames(e_Y_i_MA) <- c("Iteration", "Time", "Site", "Resp.1", "Resp.2")
e_Y_i_MI <- read.table("MI/e_Y_i.txt",
                       head = FALSE,
                       sep = ",")
colnames(e_Y_i_MI) <- c("Iteration", "Time", "Site", "Resp.1", "Resp.2")

PMSE <- function(K_obs, N_i_obs, sample, pop) {
  aux_R1 <- NULL
  for (t in 1:K_obs) {
    for (n in 1:N_i_obs) {
      est <- subset(x = sample,
                    select = "Resp.1",
                    subset = Time == t & Site == n)
      true <- subset(x = pop,
                     select = "Resp.1",
                     subset = Time == t & Site == n)
      aux_R1 <- c(aux_R1, (mean(est[,1]) - as.numeric(true))**2)
    }
  }
  aux_R2 <- NULL
  for (t in 1:K_obs) {
    for (n in 1:N_i_obs) {
      est <- subset(x = sample,
                    select = "Resp.2",
                    subset = Time == t & Site == n)
      true <- subset(x = pop,
                     select = "Resp.2",
                     subset = Time == t & Site == n)
      aux_R2 <- c(aux_R2, (mean(est[,1]) - as.numeric(true))**2)
    }
  }
  aux <- c(aux_R1, aux_R2)
  PMSE <- mean(aux)
  return(PMSE)
}

PMSE_MA <- PMSE(K_obs = K, N_i_obs = N_i, sample = e_Y_i_MA, pop = Y_i)
PMSE_MI <- PMSE(K_obs = K, N_i_obs = N_i, sample = e_Y_i_MI, pop = Y_i)

# 1 to 3
n_i <- 2

MA_stats_Y1_i <- data.frame(1:K, matrix(data = 0, nrow = K, ncol = 5))
MA_stats_Y2_i <- data.frame(1:K, matrix(data = 0, nrow = K, ncol = 5))
colnames(MA_stats_Y1_i) <- c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(MA_stats_Y2_i) <- c("Time", "LI", "Mean", "True", "LS", "Coverage")
MI_stats_Y1_i <- data.frame(1:K, matrix(data = 0, nrow = K, ncol = 5))
MI_stats_Y2_i <- data.frame(1:K, matrix(data = 0, nrow = K, ncol = 5))
colnames(MI_stats_Y1_i) <- c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(MI_stats_Y2_i) <- c("Time", "LI", "Mean", "True", "LS", "Coverage")
for(t in 1:K) {
  set_MA <- subset(x = e_Y_i_MA,
                   select = c("Iteration", "Resp.1", "Resp.2"),
                   subset = Site == n_i & Time == t)
  MA_stats_Y1_i[t,"LI"] <- quantile(x = set_MA[1:smp_size, "Resp.1"],
                                    probs = 0.025)
  MA_stats_Y2_i[t,"LI"] <- quantile(x = set_MA[1:smp_size, "Resp.2"],
                                    probs = 0.025)
  MA_stats_Y1_i[t,"Mean"] <- mean(set_MA[1:smp_size, "Resp.1"])
  MA_stats_Y2_i[t,"Mean"] <- mean(set_MA[1:smp_size, "Resp.2"])
  MA_stats_Y1_i[t,"True"] <- subset(x = Y_i,
                                    select = c("Resp.1"),
                                    subset = Site == n_i & Time == t)
  MA_stats_Y2_i[t,"True"] <- subset(x = Y_i,
                                    select = c("Resp.2"),
                                    subset = Site == n_i & Time == t)
  MA_stats_Y1_i[t,"LS"] <- quantile(x = set_MA[1:smp_size, "Resp.1"],
                                    probs = 0.975)
  MA_stats_Y2_i[t,"LS"] <- quantile(x = set_MA[1:smp_size, "Resp.2"],
                                    probs = 0.975)
  set_MI <- subset(x = e_Y_i_MI,
                   select = c("Iteration", "Resp.1", "Resp.2"),
                   subset = Site == n_i & Time == t)
  MI_stats_Y1_i[t,"LI"] <- quantile(x = set_MI[1:smp_size, "Resp.1"],
                                    probs = 0.025)
  MI_stats_Y2_i[t,"LI"] <- quantile(x = set_MI[1:smp_size, "Resp.2"],
                                    probs = 0.025)
  MI_stats_Y1_i[t,"Mean"] <- mean(set_MI[1:smp_size, "Resp.1"])
  MI_stats_Y2_i[t,"Mean"] <- mean(set_MI[1:smp_size, "Resp.2"])
  MI_stats_Y1_i[t,"True"] <- subset(x = Y_i,
                                    select = c("Resp.1"),
                                    subset = Site == n_i & Time == t)
  MI_stats_Y2_i[t,"True"] <- subset(x = Y_i,
                                    select = c("Resp.2"),
                                    subset = Site == n_i & Time == t)
  MI_stats_Y1_i[t,"LS"] <- quantile(x = set_MI[1:smp_size, "Resp.1"],
                                    probs = 0.975)
  MI_stats_Y2_i[t,"LS"] <- quantile(x = set_MI[1:smp_size, "Resp.2"],
                                    probs = 0.975) 
}

##### Diagnostics

# ECP
MA_stats_Y1_i[,"Coverage"] = MA_stats_Y1_i[,"True"] >= MA_stats_Y1_i[,"LI"] &
  MA_stats_Y1_i[,"True"] <= MA_stats_Y1_i[,"LS"]
table(MA_stats_Y1_i[,"Coverage"])
round(table(MA_stats_Y1_i[,"Coverage"])/sum(table(MA_stats_Y1_i[,"Coverage"])),
      digits = 3)*100

MI_stats_Y1_i[,"Coverage"] = MI_stats_Y1_i[,"True"] >= MI_stats_Y1_i[,"LI"] &
  MI_stats_Y1_i[,"True"] <= MI_stats_Y1_i[,"LS"]
table(MI_stats_Y1_i[,"Coverage"])
round(table(MI_stats_Y1_i[,"Coverage"])/sum(table(MI_stats_Y1_i[,"Coverage"])),
      digits = 3)*100

MA_stats_Y2_i[,"Coverage"] = MA_stats_Y2_i[,"True"] >= MA_stats_Y2_i[,"LI"] &
  MA_stats_Y2_i[,"True"] <= MA_stats_Y2_i[,"LS"]
table(MA_stats_Y2_i[,"Coverage"])
round(table(MA_stats_Y2_i[,"Coverage"])/sum(table(MA_stats_Y2_i[,"Coverage"])),
      digits = 3)*100

MI_stats_Y2_i[,"Coverage"] = MI_stats_Y2_i[,"True"] >= MI_stats_Y2_i[,"LI"] &
  MI_stats_Y2_i[,"True"] <= MI_stats_Y2_i[,"LS"]
table(MI_stats_Y2_i[,"Coverage"])
round(table(MI_stats_Y2_i[,"Coverage"])/sum(table(MI_stats_Y2_i[,"Coverage"])),
      digits = 3)*100

# IS
mean(interval_score(true_values = MA_stats_Y1_i[,"True"],
                    lower = MA_stats_Y1_i[,"LI"],
                    upper = MA_stats_Y1_i[,"LS"],
                    interval_range = 95))
mean(interval_score(true_values = MI_stats_Y1_i[,"True"],
                    lower = MI_stats_Y1_i[,"LI"],
                    upper = MI_stats_Y1_i[,"LS"],
                    interval_range = 95))
mean(interval_score(true_values = MA_stats_Y2_i[,"True"],
                    lower = MA_stats_Y2_i[,"LI"],
                    upper = MA_stats_Y2_i[,"LS"],
                    interval_range = 95))
mean(interval_score(true_values = MI_stats_Y2_i[,"True"],
                    lower = MI_stats_Y2_i[,"LI"],
                    upper = MI_stats_Y2_i[,"LS"],
                    interval_range = 95))
