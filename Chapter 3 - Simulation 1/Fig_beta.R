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
MV1 = 200; MV2 = 200; MV4 = 200

beta_MV1 = read.table("1MV/beta.txt",
                      head = FALSE,
                      sep = ",")
colnames(beta_MV1) = c("Time", "beta00", "beta01", "beta10", "beta11")
beta_MV2 = read.table("2MV/beta.txt",
                      head = FALSE,
                      sep = ",")
colnames(beta_MV2) = c("Time", "beta00", "beta01", "beta10", "beta11")
beta_MV4 = read.table("4MV/beta.txt",
                      head = FALSE,
                      sep = ",")
colnames(beta_MV4) = c("Time", "beta00", "beta01", "beta10", "beta11")
e_beta_TMV1 = read.table("1MV/e_beta.txt",
                         head = FALSE,
                         sep = ",")
colnames(e_beta_TMV1) = c("Iteration", "Time",
                          "beta00", "beta01", "beta10", "beta11")
e_beta_TMV2 = read.table("2MV/e_beta.txt",
                         head = FALSE,
                         sep = ",")
colnames(e_beta_TMV2) = c("Iteration", "Time",
                          "beta00", "beta01", "beta10", "beta11")
e_beta_TMV4 = read.table("4MV/e_beta.txt",
                         head = FALSE,
                         sep = ",")
colnames(e_beta_TMV4) = c("Iteration", "Time",
                          "beta00", "beta01", "beta10", "beta11")

TMV1_stats_b00 = data.frame(0:MV1, matrix(data = 0, nrow = MV1+1, ncol = 5))
TMV1_stats_b01 = data.frame(0:MV1, matrix(data = 0, nrow = MV1+1, ncol = 5))
TMV1_stats_b10 = data.frame(0:MV1, matrix(data = 0, nrow = MV1+1, ncol = 5))
TMV1_stats_b11 = data.frame(0:MV1, matrix(data = 0, nrow = MV1+1, ncol = 5))
colnames(TMV1_stats_b00) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(TMV1_stats_b01) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(TMV1_stats_b10) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(TMV1_stats_b11) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
TMV2_stats_b00 = data.frame(0:MV2, matrix(data = 0, nrow = MV2+1, ncol = 5))
TMV2_stats_b01 = data.frame(0:MV2, matrix(data = 0, nrow = MV2+1, ncol = 5))
TMV2_stats_b10 = data.frame(0:MV2, matrix(data = 0, nrow = MV2+1, ncol = 5))
TMV2_stats_b11 = data.frame(0:MV2, matrix(data = 0, nrow = MV2+1, ncol = 5))
colnames(TMV2_stats_b00) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(TMV2_stats_b01) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(TMV2_stats_b10) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(TMV2_stats_b11) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
TMV4_stats_b00 = data.frame(0:MV4, matrix(data = 0, nrow = MV4+1, ncol = 5))
TMV4_stats_b01 = data.frame(0:MV4, matrix(data = 0, nrow = MV4+1, ncol = 5))
TMV4_stats_b10 = data.frame(0:MV4, matrix(data = 0, nrow = MV4+1, ncol = 5))
TMV4_stats_b11 = data.frame(0:MV4, matrix(data = 0, nrow = MV4+1, ncol = 5))
colnames(TMV4_stats_b00) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(TMV4_stats_b01) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(TMV4_stats_b10) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(TMV4_stats_b11) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
for(t in 0:MV1) {
  set_TMV1 = subset(x = e_beta_TMV1,
                    subset = Time == t,
                    select = c("Iteration", "beta00", "beta01",
                               "beta10", "beta11"))
  TMV1_stats_b00[t+1,"LI"] = quantile(x = set_TMV1[1:smp_size, "beta00"],
                                      probs = 0.025)
  TMV1_stats_b01[t+1,"LI"] = quantile(x = set_TMV1[1:smp_size, "beta01"],
                                      probs = 0.025)
  TMV1_stats_b10[t+1,"LI"] = quantile(x = set_TMV1[1:smp_size, "beta10"],
                                      probs = 0.025)
  TMV1_stats_b11[t+1,"LI"] = quantile(x = set_TMV1[1:smp_size, "beta11"],
                                      probs = 0.025)
  TMV1_stats_b00[t+1,"Mean"] = mean(set_TMV1[1:smp_size, "beta00"])
  TMV1_stats_b01[t+1,"Mean"] = mean(set_TMV1[1:smp_size, "beta01"])
  TMV1_stats_b10[t+1,"Mean"] = mean(set_TMV1[1:smp_size, "beta10"])
  TMV1_stats_b11[t+1,"Mean"] = mean(set_TMV1[1:smp_size, "beta11"])
  TMV1_stats_b00[t+1,"True"] = subset(x = beta_MV1,
                                      subset = Time == t,
                                      select = c("beta00"))
  TMV1_stats_b01[t+1,"True"] = subset(x = beta_MV1,
                                      subset = Time == t,
                                      select = c("beta01"))
  TMV1_stats_b10[t+1,"True"] = subset(x = beta_MV1,
                                      subset = Time == t,
                                      select = c("beta10"))
  TMV1_stats_b11[t+1,"True"] = subset(x = beta_MV1,
                                      subset = Time == t,
                                      select = c("beta11"))
  TMV1_stats_b00[t+1,"LS"] = quantile(x = set_TMV1[1:smp_size, "beta00"],
                                      probs = 0.975)
  TMV1_stats_b01[t+1,"LS"] = quantile(x = set_TMV1[1:smp_size, "beta01"],
                                      probs = 0.975)
  TMV1_stats_b10[t+1,"LS"] = quantile(x = set_TMV1[1:smp_size, "beta10"],
                                      probs = 0.975)
  TMV1_stats_b11[t+1,"LS"] = quantile(x = set_TMV1[1:smp_size, "beta11"],
                                      probs = 0.975)
}
for(t in 0:MV2) {
  set_TMV2 = subset(x = e_beta_TMV2,
                    subset = Time == t,
                    select = c("Iteration", "beta00", "beta01",
                               "beta10", "beta11"))
  TMV2_stats_b00[t+1,"LI"] = quantile(x = set_TMV2[1:smp_size, "beta00"],
                                      probs = 0.025)
  TMV2_stats_b01[t+1,"LI"] = quantile(x = set_TMV2[1:smp_size, "beta01"],
                                      probs = 0.025)
  TMV2_stats_b10[t+1,"LI"] = quantile(x = set_TMV2[1:smp_size, "beta10"],
                                      probs = 0.025)
  TMV2_stats_b11[t+1,"LI"] = quantile(x = set_TMV2[1:smp_size, "beta11"],
                                      probs = 0.025)
  TMV2_stats_b00[t+1,"Mean"] = mean(set_TMV2[1:smp_size, "beta00"])
  TMV2_stats_b01[t+1,"Mean"] = mean(set_TMV2[1:smp_size, "beta01"])
  TMV2_stats_b10[t+1,"Mean"] = mean(set_TMV2[1:smp_size, "beta10"])
  TMV2_stats_b11[t+1,"Mean"] = mean(set_TMV2[1:smp_size, "beta11"])
  TMV2_stats_b00[t+1,"True"] = subset(x = beta_MV2,
                                      subset = Time == t,
                                      select = c("beta00"))
  TMV2_stats_b01[t+1,"True"] = subset(x = beta_MV2,
                                      subset = Time == t,
                                      select = c("beta01"))
  TMV2_stats_b10[t+1,"True"] = subset(x = beta_MV2,
                                      subset = Time == t,
                                      select = c("beta10"))
  TMV2_stats_b11[t+1,"True"] = subset(x = beta_MV2,
                                      subset = Time == t,
                                      select = c("beta11"))
  TMV2_stats_b00[t+1,"LS"] = quantile(x = set_TMV2[1:smp_size, "beta00"],
                                      probs = 0.975)
  TMV2_stats_b01[t+1,"LS"] = quantile(x = set_TMV2[1:smp_size, "beta01"],
                                      probs = 0.975)
  TMV2_stats_b10[t+1,"LS"] = quantile(x = set_TMV2[1:smp_size, "beta10"],
                                      probs = 0.975)
  TMV2_stats_b11[t+1,"LS"] = quantile(x = set_TMV2[1:smp_size, "beta11"],
                                      probs = 0.975)
}
for(t in 0:MV4) {
  set_TMV4 = subset(x = e_beta_TMV4,
                    subset = Time == t,
                    select = c("Iteration", "beta00", "beta01",
                               "beta10", "beta11"))
  TMV4_stats_b00[t+1,"LI"] = quantile(x = set_TMV4[1:smp_size, "beta00"],
                                      probs = 0.025)
  TMV4_stats_b01[t+1,"LI"] = quantile(x = set_TMV4[1:smp_size, "beta01"],
                                      probs = 0.025)
  TMV4_stats_b10[t+1,"LI"] = quantile(x = set_TMV4[1:smp_size, "beta10"],
                                      probs = 0.025)
  TMV4_stats_b11[t+1,"LI"] = quantile(x = set_TMV4[1:smp_size, "beta11"],
                                      probs = 0.025)
  TMV4_stats_b00[t+1,"Mean"] = mean(set_TMV4[1:smp_size, "beta00"])
  TMV4_stats_b01[t+1,"Mean"] = mean(set_TMV4[1:smp_size, "beta01"])
  TMV4_stats_b10[t+1,"Mean"] = mean(set_TMV4[1:smp_size, "beta10"])
  TMV4_stats_b11[t+1,"Mean"] = mean(set_TMV4[1:smp_size, "beta11"])
  TMV4_stats_b00[t+1,"True"] = subset(x = beta_MV4,
                                      subset = Time == t,
                                      select = c("beta00"))
  TMV4_stats_b01[t+1,"True"] = subset(x = beta_MV4,
                                      subset = Time == t,
                                      select = c("beta01"))
  TMV4_stats_b10[t+1,"True"] = subset(x = beta_MV4,
                                      subset = Time == t,
                                      select = c("beta10"))
  TMV4_stats_b11[t+1,"True"] = subset(x = beta_MV4,
                                      subset = Time == t,
                                      select = c("beta11"))
  TMV4_stats_b00[t+1,"LS"] = quantile(x = set_TMV4[1:smp_size, "beta00"],
                                      probs = 0.975)
  TMV4_stats_b01[t+1,"LS"] = quantile(x = set_TMV4[1:smp_size, "beta01"],
                                      probs = 0.975)
  TMV4_stats_b10[t+1,"LS"] = quantile(x = set_TMV4[1:smp_size, "beta10"],
                                      probs = 0.975)
  TMV4_stats_b11[t+1,"LS"] = quantile(x = set_TMV4[1:smp_size, "beta11"],
                                      probs = 0.975)
}

TMV1_stats_b00[,"Coverage"]=TMV1_stats_b00[,"True"] >= TMV1_stats_b00[,"LI"] &
                            TMV1_stats_b00[,"True"] <= TMV1_stats_b00[,"LS"]
table(TMV1_stats_b00[,"Coverage"])

TMV2_stats_b00[,"Coverage"]=TMV2_stats_b00[,"True"] >= TMV2_stats_b00[,"LI"] &
                            TMV2_stats_b00[,"True"] <= TMV2_stats_b00[,"LS"]
table(TMV2_stats_b00[,"Coverage"])

TMV4_stats_b00[,"Coverage"]=TMV4_stats_b00[,"True"] >= TMV4_stats_b00[,"LI"] &
  TMV4_stats_b00[,"True"] <= TMV4_stats_b00[,"LS"]
table(TMV4_stats_b00[,"Coverage"])

TMV1_stats_b01[,"Coverage"]=TMV1_stats_b01[,"True"] >= TMV1_stats_b01[,"LI"] &
  TMV1_stats_b01[,"True"] <= TMV1_stats_b01[,"LS"]
table(TMV1_stats_b01[,"Coverage"])

TMV2_stats_b01[,"Coverage"]=TMV2_stats_b01[,"True"] >= TMV2_stats_b01[,"LI"] &
  TMV2_stats_b01[,"True"] <= TMV2_stats_b01[,"LS"]
table(TMV2_stats_b01[,"Coverage"])

TMV4_stats_b01[,"Coverage"]=TMV4_stats_b01[,"True"] >= TMV4_stats_b01[,"LI"] &
  TMV4_stats_b01[,"True"] <= TMV4_stats_b01[,"LS"]
table(TMV4_stats_b01[,"Coverage"])

TMV1_stats_b10[,"Coverage"]=TMV1_stats_b10[,"True"] >= TMV1_stats_b10[,"LI"] &
  TMV1_stats_b10[,"True"] <= TMV1_stats_b10[,"LS"]
table(TMV1_stats_b10[,"Coverage"])

TMV2_stats_b10[,"Coverage"]=TMV2_stats_b10[,"True"] >= TMV2_stats_b10[,"LI"] &
  TMV2_stats_b10[,"True"] <= TMV2_stats_b10[,"LS"]
table(TMV2_stats_b10[,"Coverage"])

TMV4_stats_b10[,"Coverage"]=TMV4_stats_b10[,"True"] >= TMV4_stats_b10[,"LI"] &
  TMV4_stats_b10[,"True"] <= TMV4_stats_b10[,"LS"]
table(TMV4_stats_b10[,"Coverage"])

TMV1_stats_b11[,"Coverage"]=TMV1_stats_b11[,"True"] >= TMV1_stats_b11[,"LI"] &
  TMV1_stats_b11[,"True"] <= TMV1_stats_b11[,"LS"]
table(TMV1_stats_b11[,"Coverage"])

TMV2_stats_b11[,"Coverage"]=TMV2_stats_b11[,"True"] >= TMV2_stats_b11[,"LI"] &
  TMV2_stats_b11[,"True"] <= TMV2_stats_b11[,"LS"]
table(TMV2_stats_b11[,"Coverage"])

TMV4_stats_b11[,"Coverage"]=TMV4_stats_b11[,"True"] >= TMV4_stats_b11[,"LI"] &
  TMV4_stats_b11[,"True"] <= TMV4_stats_b11[,"LS"]
table(TMV4_stats_b11[,"Coverage"])

plot_TMV1_beta00 <- ggplot(data = TMV1_stats_b00, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["0,1,"][italic(t)])) +
  labs(title = "Case 1") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_TMV2_beta00 <- ggplot(data = TMV2_stats_b00, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["0,1,"][italic(t)])) +
  labs(title = "Case 2") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_TMV4_beta00 <- ggplot(data = TMV4_stats_b00, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["0,1,"][italic(t)])) +
  labs(title = "Case 3") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_TMV1_beta01 <- ggplot(data = TMV1_stats_b01, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["0,2,"][italic(t)])) +
  labs(title = "Case 1") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_TMV2_beta01 <- ggplot(data = TMV2_stats_b01, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["0,2,"][italic(t)])) +
  labs(title = "Case 2") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_TMV4_beta01 <- ggplot(data = TMV4_stats_b01, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["0,2,"][italic(t)])) +
  labs(title = "Case 3") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_TMV1_beta10 <- ggplot(data = TMV1_stats_b10, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["1,1,"][italic(t)])) +
  labs(title = "Case 1") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_TMV2_beta10 <- ggplot(data = TMV2_stats_b10, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["1,1,"][italic(t)])) +
  labs(title = "Case 2") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_TMV4_beta10 <- ggplot(data = TMV4_stats_b10, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["1,1,"][italic(t)])) +
  labs(title = "Case 3") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_TMV1_beta11 <- ggplot(data = TMV1_stats_b11, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["1,2,"][italic(t)])) +
  labs(title = "Case 1") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_TMV2_beta11 <- ggplot(data = TMV2_stats_b11, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["1,2,"][italic(t)])) +
  labs(title = "Case 2") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_TMV4_beta11 <- ggplot(data = TMV4_stats_b11, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["1,2,"][italic(t)])) +
  labs(title = "Case 3") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

# PDF in US Legal - Landscape
ggarrange(plot_TMV1_beta00, plot_TMV2_beta00, plot_TMV4_beta00,
          ncol = 1, nrow = 3, legend = "none")

ggarrange(plot_TMV1_beta01, plot_TMV2_beta01, plot_TMV4_beta01,
          ncol = 1, nrow = 3, legend = "none")

ggarrange(plot_TMV1_beta10, plot_TMV2_beta10, plot_TMV4_beta10,
          ncol = 1, nrow = 3, legend = "none")

ggarrange(plot_TMV1_beta11, plot_TMV2_beta11, plot_TMV4_beta11,
          ncol = 1, nrow = 3, legend = "none")
