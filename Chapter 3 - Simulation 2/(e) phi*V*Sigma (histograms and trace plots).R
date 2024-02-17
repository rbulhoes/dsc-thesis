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

V <- 0.60
phi <- 0.40
Sigma <- matrix(data = c(1.00, 0.85, 0.85, 1.00), ncol = 2)

MA_VSigma <- read.table("MA/e_VSigma.txt", head = FALSE, sep = ",")
colnames(MA_VSigma) <- c("Iteration", "VSig11", "VSig12", "VSig22")
MA_phi <- read.table(file = "MA/e_phi.txt", head = FALSE, sep = ",")
colnames(MA_phi) <- c("Iteration", "phi")
MA_phiVSig <- data.frame(MA_phi$Iteration,
                          MA_phi$phi * MA_VSigma$VSig11,
                          MA_phi$phi * MA_VSigma$VSig12,
                          MA_phi$phi * MA_VSigma$VSig22)
colnames(MA_phiVSig) <- c("Iteration", "phiVSig11", "phiVSig12", "phiVSig22")

tr_MA_phiVSig11 <- ggplot(data = MA_phiVSig[1:smp_size,c(1,2)], 
                           aes(x = Iteration, y = phiVSig11)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[1,1], colour = "purple") +
  geom_hline(yintercept = quantile(MA_phiVSig[,"phiVSig11"], 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MA_phiVSig[,"phiVSig11"], 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

tr_MA_phiVSig12 <- ggplot(data = MA_phiVSig[1:smp_size,c(1,3)], 
                           aes(x = Iteration, y = phiVSig12)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[1,2], colour = "purple") +
  geom_hline(yintercept = quantile(MA_phiVSig[,"phiVSig12"], 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MA_phiVSig[,"phiVSig12"], 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

tr_MA_phiVSig22 <- ggplot(data = MA_phiVSig[1:smp_size,c(1,4)], 
                           aes(x = Iteration, y = phiVSig22)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[2,2], colour = "purple") +
  geom_hline(yintercept = quantile(MA_phiVSig[,"phiVSig22"], 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MA_phiVSig[,"phiVSig22"], 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

MI_VSigma <- read.table("MI/e_VSigma.txt", head = FALSE, sep = ",")
colnames(MI_VSigma) <- c("Iteration", "VSig11", "VSig12", "VSig22")
MI_phi <- read.table(file = "MI/e_phi.txt", head = FALSE, sep = ",")
colnames(MI_phi) <- c("Iteration", "phi")
MI_phiVSig <- data.frame(MI_phi$Iteration,
                          MI_phi$phi * MI_VSigma$VSig11,
                          MI_phi$phi * MI_VSigma$VSig12,
                          MI_phi$phi * MI_VSigma$VSig22)
colnames(MI_phiVSig) <- c("Iteration", "phiVSig11", "phiVSig12", "phiVSig22")

tr_MI_phiVSig11 <- ggplot(data = MI_phiVSig[1:smp_size,c(1,2)], 
                           aes(x = Iteration, y = phiVSig11)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[1,1], colour = "purple") +
  geom_hline(yintercept = quantile(MI_phiVSig[,"phiVSig11"], 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MI_phiVSig[,"phiVSig11"], 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

tr_MI_phiVSig12 <- ggplot(data = MI_phiVSig[1:smp_size,c(1,3)], 
                           aes(x = Iteration, y = phiVSig12)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[1,2], colour = "purple") +
  geom_hline(yintercept = quantile(MI_phiVSig[,"phiVSig12"], 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MI_phiVSig[,"phiVSig12"], 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

tr_MI_phiVSig22 <- ggplot(data = MI_phiVSig[1:smp_size,c(1,4)], 
                           aes(x = Iteration, y = phiVSig22)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = phi*V*Sigma[2,2], colour = "purple") +
  geom_hline(yintercept = quantile(MI_phiVSig[,"phiVSig22"], 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MI_phiVSig[,"phiVSig22"], 0.975),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value")

# Salvar 4 x 6 in. - Landscape
ggarrange(tr_MA_phiVSig11, tr_MA_phiVSig12, tr_MA_phiVSig22,
          ncol = 3, nrow = 1)
ggarrange(tr_MI_phiVSig11, tr_MI_phiVSig12, tr_MI_phiVSig22,
          ncol = 3, nrow = 1)

min11 <- min(min(MA_phiVSig[,2]), min(MI_phiVSig[,2]))
min12 <- min(min(MA_phiVSig[,3]), min(MI_phiVSig[,3]))
min22 <- min(min(MA_phiVSig[,4]), min(MI_phiVSig[,4]))
max11 <- max(max(MA_phiVSig[,2]), max(MI_phiVSig[,2]))
max12 <- max(max(MA_phiVSig[,3]), max(MI_phiVSig[,3]))
max22 <- max(max(MA_phiVSig[,4]), max(MI_phiVSig[,4]))

hist_MA_phiVSig11 <- ggplot(MA_phiVSig, aes(x = phiVSig11)) + 
  geom_histogram(color = "gray50", fill = "snow1", bins = 10, 
                 mapping=aes(x = phiVSig11,
                             y = after_stat(count)/sum(after_stat(count))*100)) +
  geom_vline(xintercept = phi*V*Sigma[1,1], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MA_phiVSig[,"phiVSig11"]), colour = "gold") +
  geom_vline(xintercept = quantile(MA_phiVSig[,"phiVSig11"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MA_phiVSig[,"phiVSig11"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min11, max11), ylim = c(0, 60)) +
  labs(x = expression(paste(phi, italic(V), Sigma["1,1"])), y = "%") +
  theme(text = element_text(size = 20))

hist_MI_phiVSig11 <- ggplot(MI_phiVSig, aes(x = phiVSig11)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = phiVSig11,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 10) +
  geom_vline(xintercept = phi*V*Sigma[1,1], colour = "red", linetype="dotdash") +
  geom_vline(xintercept = mean(MI_phiVSig[,"phiVSig11"]), colour = "gold") +
  geom_vline(xintercept = quantile(MI_phiVSig[,"phiVSig11"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MI_phiVSig[,"phiVSig11"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min11, max11), ylim = c(0, 60)) +
  labs(x = expression(paste(phi, italic(V), Sigma["1,1"])), y = "%") +
  theme(text = element_text(size = 20))

hist_MA_phiVSig12 <- ggplot(MA_phiVSig, aes(x = phiVSig12)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = phiVSig12,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 10) +
  geom_vline(xintercept = phi*V*Sigma[1,2], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MA_phiVSig[,"phiVSig12"]), colour = "gold") +
  geom_vline(xintercept = quantile(MA_phiVSig[,"phiVSig12"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MA_phiVSig[,"phiVSig12"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min12, max12), ylim = c(0, 60)) +
  labs(x = expression(paste(phi, italic(V), Sigma["1,2"])), y = "%") +
  theme(text = element_text(size = 20))

hist_MI_phiVSig12 <- ggplot(MI_phiVSig, aes(x = phiVSig12)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = phiVSig12,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 10) +
  geom_vline(xintercept = phi*V*Sigma[1,2], colour = "red", linetype="dotdash") +
  geom_vline(xintercept = mean(MI_phiVSig[,"phiVSig12"]), colour = "gold") +
  geom_vline(xintercept = quantile(MI_phiVSig[,"phiVSig12"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MI_phiVSig[,"phiVSig12"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min12, max12), ylim = c(0, 60)) +
  labs(x = expression(paste(phi, italic(V), Sigma["1,2"])), y = "%") +
  theme(text = element_text(size = 20))

hist_MA_phiVSig22 <- ggplot(MA_phiVSig, aes(x = phiVSig22)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = phiVSig22,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 10) +
  geom_vline(xintercept = phi*V*Sigma[2,2], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MA_phiVSig[,"phiVSig22"]), colour = "gold") +
  geom_vline(xintercept = quantile(MA_phiVSig[,"phiVSig22"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MA_phiVSig[,"phiVSig22"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min22, max22), ylim = c(0, 50)) +
  labs(x = expression(paste(phi, italic(V), Sigma["2,2"])), y = "%") +
  theme(text = element_text(size = 20))

hist_MI_phiVSig22 <- ggplot(MI_phiVSig, aes(x = phiVSig22)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = phiVSig22,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 10) +
  geom_vline(xintercept = phi*V*Sigma[2,2], colour = "red", linetype="dotdash") +
  geom_vline(xintercept = mean(MI_phiVSig[,"phiVSig22"]), colour = "gold") +
  geom_vline(xintercept = quantile(MI_phiVSig[,"phiVSig22"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MI_phiVSig[,"phiVSig22"], 0.975),
             linetype="dashed") +
  coord_cartesian(xlim = c(min22, max22), ylim = c(0, 50)) +
  labs(x = expression(paste(phi, italic(V), Sigma["2,2"])), y = "%") +
  theme(text = element_text(size = 20))

ggarrange(hist_MA_phiVSig11, hist_MI_phiVSig11,
          ncol = 2, nrow = 1)
ggarrange(hist_MA_phiVSig12, hist_MI_phiVSig12,
          ncol = 2, nrow = 1)
ggarrange(hist_MA_phiVSig22, hist_MI_phiVSig22,
          ncol = 2, nrow = 1)
