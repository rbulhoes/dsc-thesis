##### Inform the directory 
setwd("D:/Exportações/Chap2/Sim1")
getwd()

##### Some packages
library(ggplot2)
library(ggpubr)
# library(grid)
# library(gridExtra)
library(ggrepel) # Para colocar r?tulos/etiquetas/nomes das esta??es nos mapas

smp_size <- 1000

V <- 0.4
phi <- 0.1
Sigma <- matrix(data = c(1, 0.7, 0.7, 1), ncol = 2)

MA_10_VSigma <- read.table("TamanhoDez/e_VSigma.txt", head = FALSE, sep = ",")
colnames(MA_10_VSigma) <- c("Iteration", "VSig11", "VSig12", "VSig22")
MA_10_phi <- read.table(file = "TamanhoDez/e_phi.txt", head = FALSE, sep = ",")
colnames(MA_10_phi) <- c("Iteration", "phi")
MA_10_phiVSig <- data.frame(MA_10_phi$Iteration,
                            MA_10_phi$phi * MA_10_VSigma$VSig11,
                            MA_10_phi$phi * MA_10_VSigma$VSig12,
                            MA_10_phi$phi * MA_10_VSigma$VSig22)
colnames(MA_10_phiVSig) <- c("Iteration", "phiVSig11", "phiVSig12", "phiVSig22")

tr_MA_10_phiVSig11 <- ggplot(data = MA_10_phiVSig[1:smp_size,c(1,2)], 
                             aes(x = Iteration, y = phiVSig11)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[1,1], colour = "purple") +
  geom_hline(yintercept = quantile(MA_10_phiVSig$phiVSig11, 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MA_10_phiVSig$phiVSig11, 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

tr_MA_10_phiVSig12 <- ggplot(data = MA_10_phiVSig[1:smp_size,c(1,3)], 
                             aes(x = Iteration, y = phiVSig12)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[1,2], colour = "purple") +
  geom_hline(yintercept = quantile(MA_10_phiVSig$phiVSig12, 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MA_10_phiVSig$phiVSig12, 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

tr_MA_10_phiVSig22 <- ggplot(data = MA_10_phiVSig[1:smp_size,c(1,4)], 
                             aes(x = Iteration, y = phiVSig22)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[2,2], colour = "purple") +
  geom_hline(yintercept = quantile(MA_10_phiVSig$phiVSig22, 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MA_10_phiVSig$phiVSig22, 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

MA_100_VSigma <- read.table("TamanhoCem_AF/e_VSigma.txt", head = FALSE, sep = ",")
colnames(MA_100_VSigma) <- c("Iteration", "VSig11", "VSig12", "VSig22")
MA_100_phi <- read.table(file = "TamanhoCem_AF/e_phi.txt", head = FALSE, sep = ",")
colnames(MA_100_phi) <- c("Iteration", "phi")
MA_100_phiVSig <- data.frame(MA_100_phi$Iteration,
                             MA_100_phi$phi * MA_100_VSigma$VSig11,
                             MA_100_phi$phi * MA_100_VSigma$VSig12,
                             MA_100_phi$phi * MA_100_VSigma$VSig22)
colnames(MA_100_phiVSig) <- c("Iteration", "phiVSig11", "phiVSig12", "phiVSig22")

tr_MA_100_phiVSig11 <- ggplot(data = MA_100_phiVSig[1:smp_size,c(1,2)], 
                              aes(x = Iteration, y = phiVSig11)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[1,1], colour = "purple") +
  geom_hline(yintercept = quantile(MA_100_phiVSig$phiVSig11, 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MA_100_phiVSig$phiVSig11, 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

tr_MA_100_phiVSig12 <- ggplot(data = MA_100_phiVSig[1:smp_size,c(1,3)], 
                              aes(x = Iteration, y = phiVSig12)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[1,2], colour = "purple") +
  geom_hline(yintercept = quantile(MA_100_phiVSig$phiVSig12, 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MA_100_phiVSig$phiVSig12, 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

tr_MA_100_phiVSig22 <- ggplot(data = MA_100_phiVSig[1:smp_size,c(1,4)], 
                              aes(x = Iteration, y = phiVSig22)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[2,2], colour = "purple") +
  geom_hline(yintercept = quantile(MA_100_phiVSig$phiVSig22, 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MA_100_phiVSig$phiVSig22, 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

MA_1000_VSigma <- read.table("TamanhoMil/e_VSigma.txt", head = FALSE, sep = ",")
colnames(MA_1000_VSigma) <- c("Iteration", "VSig11", "VSig12", "VSig22")
MA_1000_phi <- read.table(file = "TamanhoMil/e_phi.txt", head = FALSE, sep = ",")
colnames(MA_1000_phi) <- c("Iteration", "phi")
MA_1000_phiVSig <- data.frame(MA_1000_phi$Iteration,
                              MA_1000_phi$phi * MA_1000_VSigma$VSig11,
                              MA_1000_phi$phi * MA_1000_VSigma$VSig12,
                              MA_1000_phi$phi * MA_1000_VSigma$VSig22)
colnames(MA_1000_phiVSig) <- c("Iteration", "phiVSig11", "phiVSig12", "phiVSig22")

tr_MA_1000_phiVSig11 <- ggplot(data = MA_1000_phiVSig[1:smp_size,c(1,2)], 
                               aes(x = Iteration, y = phiVSig11)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[1,1], colour = "purple") +
  geom_hline(yintercept = quantile(MA_1000_phiVSig$phiVSig11, 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MA_1000_phiVSig$phiVSig11, 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

tr_MA_1000_phiVSig12 <- ggplot(data = MA_1000_phiVSig[1:smp_size,c(1,3)], 
                               aes(x = Iteration, y = phiVSig12)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[1,2], colour = "purple") +
  geom_hline(yintercept = quantile(MA_1000_phiVSig$phiVSig12, 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MA_1000_phiVSig$phiVSig12, 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

tr_MA_1000_phiVSig22 <- ggplot(data = MA_1000_phiVSig[1:smp_size,c(1,4)], 
                               aes(x = Iteration, y = phiVSig22)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[2,2], colour = "purple") +
  geom_hline(yintercept = quantile(MA_1000_phiVSig$phiVSig22, 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MA_1000_phiVSig$phiVSig22, 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

# Salvar 4 x 6 in. - Landscape
ggarrange(tr_MA_10_phiVSig11, tr_MA_10_phiVSig12, tr_MA_10_phiVSig22,
          ncol = 3, nrow = 1)
ggarrange(tr_MA_100_phiVSig11, tr_MA_100_phiVSig12, tr_MA_100_phiVSig22,
          ncol = 3, nrow = 1)
ggarrange(tr_MA_1000_phiVSig11, tr_MA_1000_phiVSig12, tr_MA_1000_phiVSig22,
          ncol = 3, nrow = 1)

min11 <- min(min(MA_10_phiVSig[,2]), min(MA_100_phiVSig[,2]),
             min(MA_1000_phiVSig[,2]), 0.01)
min12 <- min(min(MA_10_phiVSig[,3]), min(MA_100_phiVSig[,3]),
             min(MA_1000_phiVSig[,3]), 0.01)
min22 <- min(min(MA_10_phiVSig[,4]), min(MA_100_phiVSig[,4]),
             min(MA_1000_phiVSig[,4]), 0.01)
max11 <- max(max(MA_10_phiVSig[,2]), max(MA_100_phiVSig[,2]),
             max(MA_1000_phiVSig[,2]), 0.10)
max12 <- max(max(MA_10_phiVSig[,3]), max(MA_100_phiVSig[,3]),
             max(MA_1000_phiVSig[,3]), 0.08)
max22 <- max(max(MA_10_phiVSig[,4]), max(MA_100_phiVSig[,4]),
             max(MA_1000_phiVSig[,4]), 0.11)

hist_10_phiVSig11 <- ggplot(MA_10_phiVSig, aes(x = phiVSig11)) + 
  geom_histogram(color = "gray50", fill = "snow1", bins = 10, 
                 mapping=aes(x = phiVSig11,
                             y = after_stat(count)/sum(after_stat(count))*100)) +
  geom_vline(xintercept = phi*V*Sigma[1,1], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MA_10_phiVSig[,"phiVSig11"]), colour = "gold") +
  geom_vline(xintercept = quantile(MA_10_phiVSig[,"phiVSig11"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MA_10_phiVSig[,"phiVSig11"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min11, max11), ylim = c(0, 60)) +
  labs(x = expression(paste(phi, italic(V), Sigma["1,1"])), y = "%") +
  theme(text = element_text(size = 20))

hist_100_phiVSig11 <- ggplot(MA_100_phiVSig, aes(x = phiVSig11)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = phiVSig11,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 10) +
  geom_vline(xintercept = phi*V*Sigma[1,1], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MA_100_phiVSig[,"phiVSig11"]), colour = "gold") +
  geom_vline(xintercept = quantile(MA_100_phiVSig[,"phiVSig11"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MA_100_phiVSig[,"phiVSig11"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min11, max11), ylim = c(0, 60)) +
  labs(x = expression(paste(phi, italic(V), Sigma["1,1"])), y = "%") +
  theme(text = element_text(size = 20))

hist_1000_phiVSig11 <- ggplot(MA_1000_phiVSig, aes(x = phiVSig11)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = phiVSig11,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 10) +
  geom_vline(xintercept = phi*V*Sigma[1,1], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MA_1000_phiVSig[,"phiVSig11"]), colour = "gold") +
  geom_vline(xintercept = quantile(MA_1000_phiVSig[,"phiVSig11"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MA_1000_phiVSig[,"phiVSig11"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min11, max11), ylim = c(0, 60)) +
  labs(x = expression(paste(phi, italic(V), Sigma["1,1"])), y = "%") +
  theme(text = element_text(size = 20))

hist_10_phiVSig12 <- ggplot(MA_10_phiVSig, aes(x = phiVSig12)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = phiVSig12,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 10) +
  geom_vline(xintercept = phi*V*Sigma[1,2], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MA_10_phiVSig[,"phiVSig12"]), colour = "gold") +
  geom_vline(xintercept = quantile(MA_10_phiVSig[,"phiVSig12"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MA_10_phiVSig[,"phiVSig12"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min12, max12), ylim = c(0, 60)) +
  labs(x = expression(paste(phi, italic(V), Sigma["1,2"])), y = "%") +
  theme(text = element_text(size = 20))

hist_100_phiVSig12 <- ggplot(MA_100_phiVSig, aes(x = phiVSig12)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = phiVSig12,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 10) +
  geom_vline(xintercept = phi*V*Sigma[1,2], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MA_100_phiVSig[,"phiVSig12"]), colour = "gold") +
  geom_vline(xintercept = quantile(MA_100_phiVSig[,"phiVSig12"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MA_100_phiVSig[,"phiVSig12"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min12, max12), ylim = c(0, 60)) +
  labs(x = expression(paste(phi, italic(V), Sigma["1,2"])), y = "%") +
  theme(text = element_text(size = 20))

hist_1000_phiVSig12 <- ggplot(MA_1000_phiVSig, aes(x = phiVSig12)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = phiVSig12,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 10) +
  geom_vline(xintercept = phi*V*Sigma[1,2], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MA_1000_phiVSig[,"phiVSig12"]), colour = "gold") +
  geom_vline(xintercept = quantile(MA_1000_phiVSig[,"phiVSig12"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MA_1000_phiVSig[,"phiVSig12"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min12, max12), ylim = c(0, 60)) +
  labs(x = expression(paste(phi, italic(V), Sigma["1,2"])), y = "%") +
  theme(text = element_text(size = 20))

hist_10_phiVSig22 <- ggplot(MA_10_phiVSig, aes(x = phiVSig22)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = phiVSig22,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 10) +
  geom_vline(xintercept = phi*V*Sigma[2,2], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MA_10_phiVSig[,"phiVSig22"]), colour = "gold") +
  geom_vline(xintercept = quantile(MA_10_phiVSig[,"phiVSig22"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MA_10_phiVSig[,"phiVSig22"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min22, max22), ylim = c(0, 60)) +
  labs(x = expression(paste(phi, italic(V), Sigma["2,2"])), y = "%") +
  theme(text = element_text(size = 20))

hist_100_phiVSig22 <- ggplot(MA_100_phiVSig, aes(x = phiVSig22)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = phiVSig22,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 10) +
  geom_vline(xintercept = phi*V*Sigma[2,2], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MA_100_phiVSig[,"phiVSig22"]), colour = "gold") +
  geom_vline(xintercept = quantile(MA_100_phiVSig[,"phiVSig22"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MA_100_phiVSig[,"phiVSig22"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min22, max22), ylim = c(0, 60)) +
  labs(x = expression(paste(phi, italic(V), Sigma["2,2"])), y = "%") +
  theme(text = element_text(size = 20))

hist_1000_phiVSig22 <- ggplot(MA_1000_phiVSig, aes(x = phiVSig22)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = phiVSig22,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 10) +
  geom_vline(xintercept = phi*V*Sigma[2,2], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MA_1000_phiVSig[,"phiVSig22"]), colour = "gold") +
  geom_vline(xintercept = quantile(MA_1000_phiVSig[,"phiVSig22"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MA_1000_phiVSig[,"phiVSig22"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min22, max22), ylim = c(0, 60)) +
  labs(x = expression(paste(phi, italic(V), Sigma["2,2"])), y = "%") +
  theme(text = element_text(size = 20))

ggarrange(hist_10_phiVSig11, hist_100_phiVSig11, hist_1000_phiVSig11,
          ncol = 3, nrow = 1)
ggarrange(hist_10_phiVSig12, hist_100_phiVSig12, hist_1000_phiVSig12,
          ncol = 3, nrow = 1)
ggarrange(hist_10_phiVSig22, hist_100_phiVSig22, hist_1000_phiVSig22,
          ncol = 3, nrow = 1)
