##### Inform the directory 
setwd("D:/Exportações/Chap2/Sim1")
getwd()

##### Some packages
library(ggplot2)
library(ggpubr)
# library(grid)
# library(gridExtra)
library(ggrepel) # Para colocar r?tulos/etiquetas/nomes das esta??es nos mapas

smp_size = 1000
Dez = 10; Cem = 100; Mil = 1000

beta_Dez = read.table("TamanhoDez/beta.txt",
                      head = FALSE,
                      sep = ",")
colnames(beta_Dez) = c("Time", "beta00", "beta01", "beta10", "beta11")
beta_Cem = read.table("TamanhoCem_AF/beta.txt",
                      head = FALSE,
                      sep = ",")
colnames(beta_Cem) = c("Time", "beta00", "beta01", "beta10", "beta11")
beta_Mil = read.table("TamanhoMil/beta.txt",
                      head = FALSE,
                      sep = ",")
colnames(beta_Mil) = c("Time", "beta00", "beta01", "beta10", "beta11")
e_beta_TDez = read.table("TamanhoDez/e_beta.txt",
                         head = FALSE,
                         sep = ",")
colnames(e_beta_TDez) = c("Iteration", "Time",
                          "beta00", "beta01", "beta10", "beta11")
e_beta_TCem = read.table("TamanhoCem_AF/e_beta.txt",
                         head = FALSE,
                         sep = ",")
colnames(e_beta_TCem) = c("Iteration", "Time",
                          "beta00", "beta01", "beta10", "beta11")
e_beta_TMil = read.table("TamanhoMil/e_beta.txt",
                         head = FALSE,
                         sep = ",")
colnames(e_beta_TMil) = c("Iteration", "Time",
                          "beta00", "beta01", "beta10", "beta11")

TDez_stats_b00 = data.frame(0:Dez, matrix(data = 0, nrow = Dez+1, ncol = 5))
TDez_stats_b01 = data.frame(0:Dez, matrix(data = 0, nrow = Dez+1, ncol = 5))
TDez_stats_b10 = data.frame(0:Dez, matrix(data = 0, nrow = Dez+1, ncol = 5))
TDez_stats_b11 = data.frame(0:Dez, matrix(data = 0, nrow = Dez+1, ncol = 5))
colnames(TDez_stats_b00) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(TDez_stats_b01) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(TDez_stats_b10) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(TDez_stats_b11) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
TCem_stats_b00 = data.frame(0:Cem, matrix(data = 0, nrow = Cem+1, ncol = 5))
TCem_stats_b01 = data.frame(0:Cem, matrix(data = 0, nrow = Cem+1, ncol = 5))
TCem_stats_b10 = data.frame(0:Cem, matrix(data = 0, nrow = Cem+1, ncol = 5))
TCem_stats_b11 = data.frame(0:Cem, matrix(data = 0, nrow = Cem+1, ncol = 5))
colnames(TCem_stats_b00) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(TCem_stats_b01) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(TCem_stats_b10) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(TCem_stats_b11) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
TMil_stats_b00 = data.frame(0:Mil, matrix(data = 0, nrow = Mil+1, ncol = 5))
TMil_stats_b01 = data.frame(0:Mil, matrix(data = 0, nrow = Mil+1, ncol = 5))
TMil_stats_b10 = data.frame(0:Mil, matrix(data = 0, nrow = Mil+1, ncol = 5))
TMil_stats_b11 = data.frame(0:Mil, matrix(data = 0, nrow = Mil+1, ncol = 5))
colnames(TMil_stats_b00) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(TMil_stats_b01) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(TMil_stats_b10) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
colnames(TMil_stats_b11) = c("Time", "LI", "Mean", "True", "LS", "Coverage")
for(t in 0:Dez) {
  set_TDez = subset(x = e_beta_TDez,
                    subset = Time == t,
                    select = c("Iteration", "beta00", "beta01",
                               "beta10", "beta11"))
  TDez_stats_b00[t+1,"LI"] = quantile(x = set_TDez[1:smp_size, "beta00"],
                                      probs = 0.025)
  TDez_stats_b01[t+1,"LI"] = quantile(x = set_TDez[1:smp_size, "beta01"],
                                      probs = 0.025)
  TDez_stats_b10[t+1,"LI"] = quantile(x = set_TDez[1:smp_size, "beta10"],
                                      probs = 0.025)
  TDez_stats_b11[t+1,"LI"] = quantile(x = set_TDez[1:smp_size, "beta11"],
                                      probs = 0.025)
  TDez_stats_b00[t+1,"Mean"] = mean(set_TDez[1:smp_size, "beta00"])
  TDez_stats_b01[t+1,"Mean"] = mean(set_TDez[1:smp_size, "beta01"])
  TDez_stats_b10[t+1,"Mean"] = mean(set_TDez[1:smp_size, "beta10"])
  TDez_stats_b11[t+1,"Mean"] = mean(set_TDez[1:smp_size, "beta11"])
  TDez_stats_b00[t+1,"True"] = subset(x = beta_Dez,
                                      subset = Time == t,
                                      select = c("beta00"))
  TDez_stats_b01[t+1,"True"] = subset(x = beta_Dez,
                                      subset = Time == t,
                                      select = c("beta01"))
  TDez_stats_b10[t+1,"True"] = subset(x = beta_Dez,
                                      subset = Time == t,
                                      select = c("beta10"))
  TDez_stats_b11[t+1,"True"] = subset(x = beta_Dez,
                                      subset = Time == t,
                                      select = c("beta11"))
  TDez_stats_b00[t+1,"LS"] = quantile(x = set_TDez[1:smp_size, "beta00"],
                                      probs = 0.975)
  TDez_stats_b01[t+1,"LS"] = quantile(x = set_TDez[1:smp_size, "beta01"],
                                      probs = 0.975)
  TDez_stats_b10[t+1,"LS"] = quantile(x = set_TDez[1:smp_size, "beta10"],
                                      probs = 0.975)
  TDez_stats_b11[t+1,"LS"] = quantile(x = set_TDez[1:smp_size, "beta11"],
                                      probs = 0.975)
}
for(t in 0:Cem) {
  set_TCem = subset(x = e_beta_TCem,
                    subset = Time == t,
                    select = c("Iteration", "beta00", "beta01",
                               "beta10", "beta11"))
  TCem_stats_b00[t+1,"LI"] = quantile(x = set_TCem[1:smp_size, "beta00"],
                                      probs = 0.025)
  TCem_stats_b01[t+1,"LI"] = quantile(x = set_TCem[1:smp_size, "beta01"],
                                      probs = 0.025)
  TCem_stats_b10[t+1,"LI"] = quantile(x = set_TCem[1:smp_size, "beta10"],
                                      probs = 0.025)
  TCem_stats_b11[t+1,"LI"] = quantile(x = set_TCem[1:smp_size, "beta11"],
                                      probs = 0.025)
  TCem_stats_b00[t+1,"Mean"] = mean(set_TCem[1:smp_size, "beta00"])
  TCem_stats_b01[t+1,"Mean"] = mean(set_TCem[1:smp_size, "beta01"])
  TCem_stats_b10[t+1,"Mean"] = mean(set_TCem[1:smp_size, "beta10"])
  TCem_stats_b11[t+1,"Mean"] = mean(set_TCem[1:smp_size, "beta11"])
  TCem_stats_b00[t+1,"True"] = subset(x = beta_Cem,
                                      subset = Time == t,
                                      select = c("beta00"))
  TCem_stats_b01[t+1,"True"] = subset(x = beta_Cem,
                                      subset = Time == t,
                                      select = c("beta01"))
  TCem_stats_b10[t+1,"True"] = subset(x = beta_Cem,
                                      subset = Time == t,
                                      select = c("beta10"))
  TCem_stats_b11[t+1,"True"] = subset(x = beta_Cem,
                                      subset = Time == t,
                                      select = c("beta11"))
  TCem_stats_b00[t+1,"LS"] = quantile(x = set_TCem[1:smp_size, "beta00"],
                                      probs = 0.975)
  TCem_stats_b01[t+1,"LS"] = quantile(x = set_TCem[1:smp_size, "beta01"],
                                      probs = 0.975)
  TCem_stats_b10[t+1,"LS"] = quantile(x = set_TCem[1:smp_size, "beta10"],
                                      probs = 0.975)
  TCem_stats_b11[t+1,"LS"] = quantile(x = set_TCem[1:smp_size, "beta11"],
                                      probs = 0.975)
}
for(t in 0:Mil) {
  set_TMil = subset(x = e_beta_TMil,
                    subset = Time == t,
                    select = c("Iteration", "beta00", "beta01",
                               "beta10", "beta11"))
  TMil_stats_b00[t+1,"LI"] = quantile(x = set_TMil[1:smp_size, "beta00"],
                                      probs = 0.025)
  TMil_stats_b01[t+1,"LI"] = quantile(x = set_TMil[1:smp_size, "beta01"],
                                      probs = 0.025)
  TMil_stats_b10[t+1,"LI"] = quantile(x = set_TMil[1:smp_size, "beta10"],
                                      probs = 0.025)
  TMil_stats_b11[t+1,"LI"] = quantile(x = set_TMil[1:smp_size, "beta11"],
                                      probs = 0.025)
  TMil_stats_b00[t+1,"Mean"] = mean(set_TMil[1:smp_size, "beta00"])
  TMil_stats_b01[t+1,"Mean"] = mean(set_TMil[1:smp_size, "beta01"])
  TMil_stats_b10[t+1,"Mean"] = mean(set_TMil[1:smp_size, "beta10"])
  TMil_stats_b11[t+1,"Mean"] = mean(set_TMil[1:smp_size, "beta11"])
  TMil_stats_b00[t+1,"True"] = subset(x = beta_Mil,
                                      subset = Time == t,
                                      select = c("beta00"))
  TMil_stats_b01[t+1,"True"] = subset(x = beta_Mil,
                                      subset = Time == t,
                                      select = c("beta01"))
  TMil_stats_b10[t+1,"True"] = subset(x = beta_Mil,
                                      subset = Time == t,
                                      select = c("beta10"))
  TMil_stats_b11[t+1,"True"] = subset(x = beta_Mil,
                                      subset = Time == t,
                                      select = c("beta11"))
  TMil_stats_b00[t+1,"LS"] = quantile(x = set_TMil[1:smp_size, "beta00"],
                                      probs = 0.975)
  TMil_stats_b01[t+1,"LS"] = quantile(x = set_TMil[1:smp_size, "beta01"],
                                      probs = 0.975)
  TMil_stats_b10[t+1,"LS"] = quantile(x = set_TMil[1:smp_size, "beta10"],
                                      probs = 0.975)
  TMil_stats_b11[t+1,"LS"] = quantile(x = set_TMil[1:smp_size, "beta11"],
                                      probs = 0.975)
}

TDez_stats_b00[,"Coverage"]=TDez_stats_b00[,"True"] >= TDez_stats_b00[,"LI"] &
                            TDez_stats_b00[,"True"] <= TDez_stats_b00[,"LS"]
table(TDez_stats_b00[,"Coverage"])

TCem_stats_b00[,"Coverage"]=TCem_stats_b00[,"True"] >= TCem_stats_b00[,"LI"] &
                            TCem_stats_b00[,"True"] <= TCem_stats_b00[,"LS"]
table(TCem_stats_b00[,"Coverage"])

TMil_stats_b00[,"Coverage"]=TMil_stats_b00[,"True"] >= TMil_stats_b00[,"LI"] &
  TMil_stats_b00[,"True"] <= TMil_stats_b00[,"LS"]
table(TMil_stats_b00[,"Coverage"])

TDez_stats_b01[,"Coverage"]=TDez_stats_b01[,"True"] >= TDez_stats_b01[,"LI"] &
  TDez_stats_b01[,"True"] <= TDez_stats_b01[,"LS"]
table(TDez_stats_b01[,"Coverage"])

TCem_stats_b01[,"Coverage"]=TCem_stats_b01[,"True"] >= TCem_stats_b01[,"LI"] &
  TCem_stats_b01[,"True"] <= TCem_stats_b01[,"LS"]
table(TCem_stats_b01[,"Coverage"])

TMil_stats_b01[,"Coverage"]=TMil_stats_b01[,"True"] >= TMil_stats_b01[,"LI"] &
  TMil_stats_b01[,"True"] <= TMil_stats_b01[,"LS"]
table(TMil_stats_b01[,"Coverage"])

TDez_stats_b10[,"Coverage"]=TDez_stats_b10[,"True"] >= TDez_stats_b10[,"LI"] &
  TDez_stats_b10[,"True"] <= TDez_stats_b10[,"LS"]
table(TDez_stats_b10[,"Coverage"])

TCem_stats_b10[,"Coverage"]=TCem_stats_b10[,"True"] >= TCem_stats_b10[,"LI"] &
  TCem_stats_b10[,"True"] <= TCem_stats_b10[,"LS"]
table(TCem_stats_b10[,"Coverage"])

TMil_stats_b10[,"Coverage"]=TMil_stats_b10[,"True"] >= TMil_stats_b10[,"LI"] &
  TMil_stats_b10[,"True"] <= TMil_stats_b10[,"LS"]
table(TMil_stats_b10[,"Coverage"])

TDez_stats_b11[,"Coverage"]=TDez_stats_b11[,"True"] >= TDez_stats_b11[,"LI"] &
  TDez_stats_b11[,"True"] <= TDez_stats_b11[,"LS"]
table(TDez_stats_b11[,"Coverage"])

TCem_stats_b11[,"Coverage"]=TCem_stats_b11[,"True"] >= TCem_stats_b11[,"LI"] &
  TCem_stats_b11[,"True"] <= TCem_stats_b11[,"LS"]
table(TCem_stats_b11[,"Coverage"])

TMil_stats_b11[,"Coverage"]=TMil_stats_b11[,"True"] >= TMil_stats_b11[,"LI"] &
  TMil_stats_b11[,"True"] <= TMil_stats_b11[,"LS"]
table(TMil_stats_b11[,"Coverage"])

plot_TDez_beta00 <- ggplot(data = TDez_stats_b00, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["0,1,"][italic(t)])) +
  labs(title = expression(italic(T) == 10)) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10))

plot_TCem_beta00 <- ggplot(data = TCem_stats_b00, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["0,1,"][italic(t)])) +
  labs(title = expression(italic(T) == 100)) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_TMil_beta00 <- ggplot(data = TMil_stats_b00, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["0,1,"][italic(t)])) +
  labs(title = expression(italic(T) == 1000)) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_TDez_beta01 <- ggplot(data = TDez_stats_b01, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["0,2,"][italic(t)])) +
  labs(title = expression(italic(T) == 10)) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10))

plot_TCem_beta01 <- ggplot(data = TCem_stats_b01, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["0,2,"][italic(t)])) +
  labs(title = expression(italic(T) == 100)) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_TMil_beta01 <- ggplot(data = TMil_stats_b01, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["0,2,"][italic(t)])) +
  labs(title = expression(italic(T) == 1000)) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_TDez_beta10 <- ggplot(data = TDez_stats_b10, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["1,1,"][italic(t)])) +
  labs(title = expression(italic(T) == 10)) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10))

plot_TCem_beta10 <- ggplot(data = TCem_stats_b10, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["1,1,"][italic(t)])) +
  labs(title = expression(italic(T) == 100)) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_TMil_beta10 <- ggplot(data = TMil_stats_b10, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["1,1,"][italic(t)])) +
  labs(title = expression(italic(T) == 1000)) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_TDez_beta11 <- ggplot(data = TDez_stats_b11, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["1,2,"][italic(t)])) +
  labs(title = expression(italic(T) == 10)) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10))

plot_TCem_beta11 <- ggplot(data = TCem_stats_b11, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["1,2,"][italic(t)])) +
  labs(title = expression(italic(T) == 100)) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

plot_TMil_beta11 <- ggplot(data = TMil_stats_b11, aes(x = Time, y = True)) +
  geom_point(size = 1.5, aes(color = Coverage)) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values=c("red", "blue")) +
  geom_line(aes(y = LI), linetype = 2) +
  geom_line(aes(y = LS), linetype = 2) +
  geom_line(aes(y = Mean), color = "gold") +
  labs(x = expression(italic(t)),
       y = expression(beta["1,2,"][italic(t)])) +
  labs(title = expression(italic(T) == 1000)) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

# PDF in US Legal - Landscape
ggarrange(plot_TDez_beta00, plot_TCem_beta00, plot_TMil_beta00,
          ncol = 1, nrow = 3, legend = "none")

ggarrange(plot_TDez_beta01, plot_TCem_beta01, plot_TMil_beta01,
          ncol = 1, nrow = 3, legend = "none")

ggarrange(plot_TDez_beta10, plot_TCem_beta10, plot_TMil_beta10,
          ncol = 1, nrow = 3, legend = "none")

ggarrange(plot_TDez_beta11, plot_TCem_beta11, plot_TMil_beta11,
          ncol = 1, nrow = 3, legend = "none")
