##### Inform the directory 
setwd("D:/Exportações/Chap2/Sim2")
getwd()

##### Some packages
library(ggplot2)
library(ggpubr)
library(ggrepel) # Para colocar r?tulos/etiquetas/nomes das esta??es nos mapas

smp_size <- 1000

V <- 1.00
Sigma <- matrix(data = c(1.00, 0.80, 0.80, 1.00), ncol = 2)

MA_VSigma <- read.table("MA/e_VSigma.txt", head = FALSE, sep = ",")
colnames(MA_VSigma) <- c("Iteration", "VSig11", "VSig12", "VSig22")
MI_VSigma <- read.table("MI/e_VSigma.txt", head = FALSE, sep = ",")
colnames(MI_VSigma) <- c("Iteration", "VSig11", "VSig12", "VSig22")

tr_MA_VSig11 <- ggplot(data = MA_VSigma[1:smp_size,c(1,2)], 
                       aes(x = Iteration, y = VSig11)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = V*Sigma[1,1], colour = "purple") +
  geom_hline(yintercept = quantile(MA_VSigma[,"VSig11"], 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MA_VSigma[,"VSig11"], 0.975),
             linetype = "dashed") +
  xlab(label = "Sample") + ylab(label = "Sampled value")

tr_MI_VSig11 <- ggplot(data = MI_VSigma[1:smp_size,c(1,2)], 
                       aes(x = Iteration, y = VSig11)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = V*Sigma[1,1], colour = "purple") +
  geom_hline(yintercept = quantile(MI_VSigma[,"VSig11"], 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MI_VSigma[,"VSig11"], 0.975),
             linetype = "dashed") +
  xlab(label = "Sample") + ylab(label = "Sampled value")

tr_MA_VSig12 <- ggplot(data = MA_VSigma[1:smp_size,c(1,3)], 
                       aes(x = Iteration, y = VSig12)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = V*Sigma[1,2], colour = "purple") +
  geom_hline(yintercept = quantile(MA_VSigma[,"VSig12"], 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MA_VSigma[,"VSig12"], 0.975),
             linetype = "dashed") +
  xlab(label = "Sample") + ylab(label = "Sampled value")

tr_MI_VSig12 <- ggplot(data = MI_VSigma[1:smp_size,c(1,3)], 
                       aes(x = Iteration, y = VSig12)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = V*Sigma[1,2], colour = "purple") +
  geom_hline(yintercept = quantile(MI_VSigma[,"VSig12"], 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MI_VSigma[,"VSig12"], 0.975),
             linetype = "dashed") +
  xlab(label = "Sample") + ylab(label = "Sampled value")

tr_MA_VSig22 <- ggplot(data = MA_VSigma[1:smp_size,c(1,4)], 
                       aes(x = Iteration, y = VSig22)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = V*Sigma[2,2], colour = "purple") +
  geom_hline(yintercept = quantile(MA_VSigma[,"VSig22"], 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MA_VSigma[,"VSig22"], 0.975),
             linetype = "dashed") +
  xlab(label = "Sample") + ylab(label = "Sampled value")

tr_MI_VSig22 <- ggplot(data = MI_VSigma[1:smp_size,c(1,4)], 
                       aes(x = Iteration, y = VSig22)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = V*Sigma[2,2], colour = "purple") +
  geom_hline(yintercept = quantile(MI_VSigma[,"VSig22"], 0.025),
             linetype = "dashed") +
  geom_hline(yintercept = quantile(MI_VSigma[,"VSig22"], 0.975),
             linetype = "dashed") +
  xlab(label = "Sample") + ylab(label = "Sampled value")

# Salvar 4 x 6 in. - Landscape
ggarrange(tr_MA_VSig11, tr_MA_VSig12, tr_MA_VSig22,
          tr_MI_VSig11, tr_MI_VSig12, tr_MI_VSig22,
          ncol = 3, nrow = 2)

min11 <- min(0.5, min(MA_VSigma[,2]), min(MI_VSigma[,2]))
min12 <- min(0.3, min(MA_VSigma[,3]), min(MI_VSigma[,3]))
min22 <- min(0.3, min(MA_VSigma[,4]), min(MI_VSigma[,4]))
max11 <- max(max(MA_VSigma[,2]), max(MI_VSigma[,2]))
max12 <- max(max(MA_VSigma[,3]), max(MI_VSigma[,3]))
max22 <- max(max(MA_VSigma[,4]), max(MI_VSigma[,4]))

hist_MA_VSig11 <- ggplot(MA_VSigma, aes(x = VSig11)) + 
  geom_histogram(color = "gray50", fill = "snow1", bins = 10, 
                 mapping=aes(x = VSig11,
                             y = after_stat(count)/sum(after_stat(count))*100)) +
  geom_vline(xintercept = V*Sigma[1,1], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MA_VSigma[,"VSig11"]), colour = "gold") +
  geom_vline(xintercept = quantile(MA_VSigma[,"VSig11"], 0.025),
             linetype = "dashed") +
  geom_vline(xintercept = quantile(MA_VSigma[,"VSig11"], 0.975),
             linetype = "dashed") +
  coord_cartesian(xlim = c(min11, max11), ylim = c(0, 41)) +
  labs(x = expression(paste(italic(V), Sigma["1,1"])), y = "%") +
  theme(text = element_text(size = 20))

hist_MI_VSig11 <- ggplot(MI_VSigma, aes(x = VSig11)) + 
  geom_histogram(color = "gray50", fill = "snow1", bins = 10,
                 mapping=aes(x = VSig11,
                             y = after_stat(count)/sum(after_stat(count))*100)) +
  geom_vline(xintercept = V*Sigma[1,1], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MI_VSigma[,"VSig11"]), colour = "gold") +
  geom_vline(xintercept = quantile(MI_VSigma[,"VSig11"], 0.025),
             linetype = "dashed") +
  geom_vline(xintercept = quantile(MI_VSigma[,"VSig11"], 0.975),
             linetype = "dashed") +
  coord_cartesian(xlim = c(min11, max11), ylim = c(0, 41)) +
  labs(x = expression(paste(italic(V), Sigma["1,1"])), y = "%") +
  theme(text = element_text(size = 20))

hist_MA_VSig12 <- ggplot(MA_VSigma, aes(x = VSig12)) + 
  geom_histogram(color = "gray50", fill = "snow1", bins = 10, 
                 mapping=aes(x = VSig12,
                             y = after_stat(count)/sum(after_stat(count))*100)) +
  geom_vline(xintercept = V*Sigma[1,2], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MA_VSigma[,"VSig12"]), colour = "gold") +
  geom_vline(xintercept = quantile(MA_VSigma[,"VSig12"], 0.025),
             linetype = "dashed") +
  geom_vline(xintercept = quantile(MA_VSigma[,"VSig12"], 0.975),
             linetype = "dashed") +
  coord_cartesian(xlim = c(min12, max12), ylim = c(0, 41)) +
  labs(x = expression(paste(italic(V), Sigma["1,2"])), y = "%") +
  theme(text = element_text(size = 20))

hist_MI_VSig12 <- ggplot(MI_VSigma, aes(x = VSig12)) + 
  geom_histogram(color = "gray50", fill = "snow1", bins = 10,
                 mapping=aes(x = VSig12,
                             y = after_stat(count)/sum(after_stat(count))*100)) +
  geom_vline(xintercept = V*Sigma[1,2], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MI_VSigma[,"VSig12"]), colour = "gold") +
  geom_vline(xintercept = quantile(MI_VSigma[,"VSig12"], 0.025),
             linetype = "dashed") +
  geom_vline(xintercept = quantile(MI_VSigma[,"VSig12"], 0.975),
             linetype = "dashed") +
  coord_cartesian(xlim = c(min12, max12), ylim = c(0, 41)) +
  labs(x = expression(paste(italic(V), Sigma["1,2"])), y = "%") +
  theme(text = element_text(size = 20))

hist_MA_VSig22 <- ggplot(MA_VSigma, aes(x = VSig22)) + 
  geom_histogram(color = "gray50", fill = "snow1", bins = 10, 
                 mapping=aes(x = VSig22,
                             y = after_stat(count)/sum(after_stat(count))*100)) +
  geom_vline(xintercept = V*Sigma[2,2], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MA_VSigma[,"VSig22"]), colour = "gold") +
  geom_vline(xintercept = quantile(MA_VSigma[,"VSig22"], 0.025),
             linetype = "dashed") +
  geom_vline(xintercept = quantile(MA_VSigma[,"VSig22"], 0.975),
             linetype = "dashed") +
  coord_cartesian(xlim = c(min22, max22), ylim = c(0, 41)) +
  labs(x = expression(paste(italic(V), Sigma["2,2"])), y = "%") +
  theme(text = element_text(size = 20))

hist_MI_VSig22 <- ggplot(MI_VSigma, aes(x = VSig22)) + 
  geom_histogram(color = "gray50", fill = "snow1", bins = 10,
                 mapping=aes(x = VSig12,
                             y = after_stat(count)/sum(after_stat(count))*100)) +
  geom_vline(xintercept = V*Sigma[2,2], colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = mean(MI_VSigma[,"VSig22"]), colour = "gold") +
  geom_vline(xintercept = quantile(MI_VSigma[,"VSig22"], 0.025),
             linetype = "dashed") +
  geom_vline(xintercept = quantile(MI_VSigma[,"VSig22"], 0.975),
             linetype = "dashed") +
  coord_cartesian(xlim = c(min22, max22), ylim = c(0, 41)) +
  labs(x = expression(paste(italic(V), Sigma["2,2"])), y = "%") +
  theme(text = element_text(size = 20))

# Salvar 4 x 6 in. - Landscape
ggarrange(hist_MA_VSig11, hist_MA_VSig12, hist_MA_VSig22,
          hist_MI_VSig11, hist_MI_VSig12, hist_MI_VSig22,
          ncol = 3, nrow = 2)
