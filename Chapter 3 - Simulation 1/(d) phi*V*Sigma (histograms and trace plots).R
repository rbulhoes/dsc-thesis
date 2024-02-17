##### Inform the directory 
setwd("D:/Exportações/Chap3/Sim1")
getwd()

##### Some packages
library(ggplot2)
library(ggpubr)
# library(grid)
# library(gridExtra)
library(ggrepel) # Para colocar r?tulos/etiquetas/nomes das esta??es nos mapas

smp_size <- 1000

V <- 0.1
phi <- 0.3
Sigma <- matrix(data = c(1.00, 0.75, 0.75, 1.00), ncol = 2)

MV1_VSigma <- read.table("1MV/e_VSigma.txt", head = FALSE, sep = ",")
colnames(MV1_VSigma) <- c("Iteration", "VSig11", "VSig12", "VSig22")
MV1_phi <- read.table(file = "1MV/e_phi.txt", head = FALSE, sep = ",")
colnames(MV1_phi) <- c("Iteration", "phi")
MV1_phiVSig <- data.frame(MV1_phi$Iteration,
                            MV1_phi$phi * MV1_VSigma$VSig11,
                            MV1_phi$phi * MV1_VSigma$VSig12,
                            MV1_phi$phi * MV1_VSigma$VSig22)
colnames(MV1_phiVSig) <- c("Iteration", "phiVSig11", "phiVSig12", "phiVSig22")

tr_MV1_phiVSig11 <- ggplot(data = MV1_phiVSig[1:smp_size,c(1,2)], 
                             aes(x = Iteration, y = phiVSig11)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[1,1], colour = "purple") +
  geom_hline(yintercept = quantile(MV1_phiVSig[,"phiVSig11"], 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MV1_phiVSig[,"phiVSig11"], 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

tr_MV1_phiVSig12 <- ggplot(data = MV1_phiVSig[1:smp_size,c(1,3)], 
                             aes(x = Iteration, y = phiVSig12)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[1,2], colour = "purple") +
  geom_hline(yintercept = quantile(MV1_phiVSig[,"phiVSig12"], 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MV1_phiVSig[,"phiVSig12"], 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

tr_MV1_phiVSig22 <- ggplot(data = MV1_phiVSig[1:smp_size,c(1,4)], 
                             aes(x = Iteration, y = phiVSig22)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[2,2], colour = "purple") +
  geom_hline(yintercept = quantile(MV1_phiVSig[,"phiVSig22"], 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MV1_phiVSig[,"phiVSig22"], 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

MV2_VSigma <- read.table("2MV/e_VSigma.txt", head = FALSE, sep = ",")
colnames(MV2_VSigma) <- c("Iteration", "VSig11", "VSig12", "VSig22")
MV2_phi <- read.table(file = "2MV/e_phi.txt", head = FALSE, sep = ",")
colnames(MV2_phi) <- c("Iteration", "phi")
MV2_phiVSig <- data.frame(MV2_phi$Iteration,
                             MV2_phi$phi * MV2_VSigma$VSig11,
                             MV2_phi$phi * MV2_VSigma$VSig12,
                             MV2_phi$phi * MV2_VSigma$VSig22)
colnames(MV2_phiVSig) <- c("Iteration", "phiVSig11", "phiVSig12", "phiVSig22")

tr_MV2_phiVSig11 <- ggplot(data = MV2_phiVSig[1:smp_size,c(1,2)], 
                              aes(x = Iteration, y = phiVSig11)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[1,1], colour = "purple") +
  geom_hline(yintercept = quantile(MV2_phiVSig[,"phiVSig11"], 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MV2_phiVSig[,"phiVSig11"], 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

tr_MV2_phiVSig12 <- ggplot(data = MV2_phiVSig[1:smp_size,c(1,3)], 
                              aes(x = Iteration, y = phiVSig12)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[1,2], colour = "purple") +
  geom_hline(yintercept = quantile(MV2_phiVSig[,"phiVSig12"], 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MV2_phiVSig[,"phiVSig12"], 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

tr_MV2_phiVSig22 <- ggplot(data = MV2_phiVSig[1:smp_size,c(1,4)], 
                              aes(x = Iteration, y = phiVSig22)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[2,2], colour = "purple") +
  geom_hline(yintercept = quantile(MV2_phiVSig[,"phiVSig22"], 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MV2_phiVSig[,"phiVSig22"], 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

MV4_VSigma <- read.table("4MV/e_VSigma.txt", head = FALSE, sep = ",")
colnames(MV4_VSigma) <- c("Iteration", "VSig11", "VSig12", "VSig22")
MV4_phi <- read.table(file = "4MV/e_phi.txt", head = FALSE, sep = ",")
colnames(MV4_phi) <- c("Iteration", "phi")
MV4_phiVSig <- data.frame(MV4_phi$Iteration,
                              MV4_phi$phi * MV4_VSigma$VSig11,
                              MV4_phi$phi * MV4_VSigma$VSig12,
                              MV4_phi$phi * MV4_VSigma$VSig22)
colnames(MV4_phiVSig) <- c("Iteration", "phiVSig11", "phiVSig12", "phiVSig22")

tr_MV4_phiVSig11 <- ggplot(data = MV4_phiVSig[1:smp_size,c(1,2)], 
                               aes(x = Iteration, y = phiVSig11)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[1,1], colour = "purple") +
  geom_hline(yintercept = quantile(MV4_phiVSig[,"phiVSig11"], 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MV4_phiVSig[,"phiVSig11"], 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

tr_MV4_phiVSig12 <- ggplot(data = MV4_phiVSig[1:smp_size,c(1,3)], 
                               aes(x = Iteration, y = phiVSig12)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[1,2], colour = "purple") +
  geom_hline(yintercept = quantile(MV4_phiVSig[,"phiVSig12"], 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MV4_phiVSig[,"phiVSig12"], 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

tr_MV4_phiVSig22 <- ggplot(data = MV4_phiVSig[1:smp_size,c(1,4)], 
                               aes(x = Iteration, y = phiVSig22)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[2,2], colour = "purple") +
  geom_hline(yintercept = quantile(MV4_phiVSig[,"phiVSig22"], 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MV4_phiVSig[,"phiVSig22"], 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

# Salvar 4 x 6 in. - Landscape
ggarrange(tr_MV1_phiVSig11, tr_MV1_phiVSig12, tr_MV1_phiVSig22,
          ncol = 3, nrow = 1)
ggarrange(tr_MV2_phiVSig11, tr_MV2_phiVSig12, tr_MV2_phiVSig22,
          ncol = 3, nrow = 1)
ggarrange(tr_MV4_phiVSig11, tr_MV4_phiVSig12, tr_MV4_phiVSig22,
          ncol = 3, nrow = 1)

min11 <- min(min(MV1_phiVSig[,2]), min(MV2_phiVSig[,2]),
             min(MV4_phiVSig[,2]))
min12 <- min(min(MV1_phiVSig[,3]), min(MV2_phiVSig[,3]),
             min(MV4_phiVSig[,3]))
min22 <- min(min(MV1_phiVSig[,4]), min(MV2_phiVSig[,4]),
             min(MV4_phiVSig[,4]))
max11 <- max(max(MV1_phiVSig[,2]), max(MV2_phiVSig[,2]),
             max(MV4_phiVSig[,2]))
max12 <- max(max(MV1_phiVSig[,3]), max(MV2_phiVSig[,3]),
             max(MV4_phiVSig[,3]))
max22 <- max(max(MV1_phiVSig[,4]), max(MV2_phiVSig[,4]),
             max(MV4_phiVSig[,4]))

hist_MV1_phiVSig11 <- ggplot(MV1_phiVSig, aes(x = phiVSig11)) + 
  geom_histogram(color = "gray50", fill = "snow1", bins = 10, 
                 mapping=aes(x = phiVSig11,
                             y = after_stat(count)/sum(after_stat(count))*100)) +
  geom_vline(xintercept = phi*V*Sigma[1,1], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MV1_phiVSig[,"phiVSig11"]), colour = "gold") +
  geom_vline(xintercept = quantile(MV1_phiVSig[,"phiVSig11"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MV1_phiVSig[,"phiVSig11"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min11, max11), ylim = c(0, 50)) +
  labs(x = expression(paste(phi, italic(V), Sigma["1,1"])), y = "%") +
  theme(text = element_text(size = 20))

hist_MV2_phiVSig11 <- ggplot(MV2_phiVSig, aes(x = phiVSig11)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = phiVSig11,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 10) +
  geom_vline(xintercept = phi*V*Sigma[1,1], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MV2_phiVSig[,"phiVSig11"]), colour = "gold") +
  geom_vline(xintercept = quantile(MV2_phiVSig[,"phiVSig11"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MV2_phiVSig[,"phiVSig11"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min11, max11), ylim = c(0, 50)) +
  labs(x = expression(paste(phi, italic(V), Sigma["1,1"])), y = "%") +
  theme(text = element_text(size = 20))

hist_MV4_phiVSig11 <- ggplot(MV4_phiVSig, aes(x = phiVSig11)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = phiVSig11,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 10) +
  geom_vline(xintercept = phi*V*Sigma[1,1], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MV4_phiVSig[,"phiVSig11"]), colour = "gold") +
  geom_vline(xintercept = quantile(MV4_phiVSig[,"phiVSig11"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MV4_phiVSig[,"phiVSig11"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min11, max11), ylim = c(0, 50)) +
  labs(x = expression(paste(phi, italic(V), Sigma["1,1"])), y = "%") +
  theme(text = element_text(size = 20))

hist_MV1_phiVSig12 <- ggplot(MV1_phiVSig, aes(x = phiVSig12)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = phiVSig12,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 10) +
  geom_vline(xintercept = phi*V*Sigma[1,2], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MV1_phiVSig[,"phiVSig12"]), colour = "gold") +
  geom_vline(xintercept = quantile(MV1_phiVSig[,"phiVSig12"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MV1_phiVSig[,"phiVSig12"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min12, max12), ylim = c(0, 50)) +
  labs(x = expression(paste(phi, italic(V), Sigma["1,2"])), y = "%") +
  theme(text = element_text(size = 20))

hist_MV2_phiVSig12 <- ggplot(MV2_phiVSig, aes(x = phiVSig12)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = phiVSig12,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 10) +
  geom_vline(xintercept = phi*V*Sigma[1,2], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MV2_phiVSig[,"phiVSig12"]), colour = "gold") +
  geom_vline(xintercept = quantile(MV2_phiVSig[,"phiVSig12"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MV2_phiVSig[,"phiVSig12"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min12, max12), ylim = c(0, 50)) +
  labs(x = expression(paste(phi, italic(V), Sigma["1,2"])), y = "%") +
  theme(text = element_text(size = 20))

hist_MV4_phiVSig12 <- ggplot(MV4_phiVSig, aes(x = phiVSig12)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = phiVSig12,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 10) +
  geom_vline(xintercept = phi*V*Sigma[1,2], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MV4_phiVSig[,"phiVSig12"]), colour = "gold") +
  geom_vline(xintercept = quantile(MV4_phiVSig[,"phiVSig12"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MV4_phiVSig[,"phiVSig12"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min12, max12), ylim = c(0, 50)) +
  labs(x = expression(paste(phi, italic(V), Sigma["1,2"])), y = "%") +
  theme(text = element_text(size = 20))

hist_MV1_phiVSig22 <- ggplot(MV1_phiVSig, aes(x = phiVSig22)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = phiVSig22,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 10) +
  geom_vline(xintercept = phi*V*Sigma[2,2], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MV1_phiVSig[,"phiVSig22"]), colour = "gold") +
  geom_vline(xintercept = quantile(MV1_phiVSig[,"phiVSig22"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MV1_phiVSig[,"phiVSig22"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min22, max22), ylim = c(0, 60)) +
  labs(x = expression(paste(phi, italic(V), Sigma["2,2"])), y = "%") +
  theme(text = element_text(size = 20))

hist_MV2_phiVSig22 <- ggplot(MV2_phiVSig, aes(x = phiVSig22)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = phiVSig22,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 10) +
  geom_vline(xintercept = phi*V*Sigma[2,2], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MV2_phiVSig[,"phiVSig22"]), colour = "gold") +
  geom_vline(xintercept = quantile(MV2_phiVSig[,"phiVSig22"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MV2_phiVSig[,"phiVSig22"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min22, max22), ylim = c(0, 60)) +
  labs(x = expression(paste(phi, italic(V), Sigma["2,2"])), y = "%") +
  theme(text = element_text(size = 20))

hist_MV4_phiVSig22 <- ggplot(MV4_phiVSig, aes(x = phiVSig22)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = phiVSig22,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 10) +
  geom_vline(xintercept = phi*V*Sigma[2,2], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MV4_phiVSig[,"phiVSig22"]), colour = "gold") +
  geom_vline(xintercept = quantile(MV4_phiVSig[,"phiVSig22"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MV4_phiVSig[,"phiVSig22"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min22, max22), ylim = c(0, 60)) +
  labs(x = expression(paste(phi, italic(V), Sigma["2,2"])), y = "%") +
  theme(text = element_text(size = 20))

ggarrange(hist_MV1_phiVSig11, hist_MV2_phiVSig11, hist_MV4_phiVSig11,
          ncol = 3, nrow = 1)
ggarrange(hist_MV1_phiVSig12, hist_MV2_phiVSig12, hist_MV4_phiVSig12,
          ncol = 3, nrow = 1)
ggarrange(hist_MV1_phiVSig22, hist_MV2_phiVSig22, hist_MV4_phiVSig22,
          ncol = 3, nrow = 1)
