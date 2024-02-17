##### Inform the directory 
setwd("D:/Exportações/Chap2/Sim2")
getwd()

##### Some packages
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(ggrepel) # Para colocar r?tulos/etiquetas/nomes das esta??es nos mapas

smp_size <- 1000
N <- 16

true_phi <- 0.5

MA_phi <- read.table("MA/e_phi.txt", head = FALSE, sep = ",")
colnames(MA_phi) <- c("Iteration", "Value")
MI_phi <- read.table("MI/e_phi.txt", head = FALSE, sep = ",")
colnames(MI_phi) <- c("Iteration", "Value")

tr_MA_phi <- ggplot(data = MA_phi, 
                    aes(x = Iteration, y = Value)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = true_phi, colour = "purple") +
  geom_hline(yintercept = quantile(MA_phi[,"Value"], 0.025),
             linetype="dashed") +
  geom_hline(yintercept = quantile(MA_phi[,"Value"], 0.975),
             linetype="dashed") +
  xlab(label = "Sample") + ylab(label = "Sampled value")

tr_MI_phi <- ggplot(data = MI_phi, 
                    aes(x = Iteration, y = Value)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = true_phi, colour = "purple") +
  geom_hline(yintercept = quantile(MI_phi[,"Value"], 0.025),
             linetype="dashed") +
  geom_hline(yintercept = quantile(MI_phi[,"Value"], 0.975),
             linetype="dashed") +
  xlab(label = "Sample") + ylab(label = "Sampled value")

minimo <- min(min(MA_phi[,2]), min(MI_phi[,2]))
maximo <- max(max(MA_phi[,2]), max(MI_phi[,2]))

hist_MA_phi <- ggplot(MA_phi, aes(x = Value)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = Value,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 15) +
  geom_vline(xintercept = (true_phi+0.005), colour = "blue", linetype="dotdash") +
  geom_vline(xintercept = (mean(MA_phi[,"Value"])-0.005), colour = "gold") +
  geom_vline(xintercept = quantile(MA_phi[,"Value"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MA_phi[,"Value"], 0.975),
             linetype="dashed") +
  labs(x = expression(phi), y = "%") +
#  labs(title = expression(paste("Histogram - ",phi))) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  xlim(minimo, maximo)

hist_MI_phi <- ggplot(MI_phi, aes(x = Value)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = Value,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 15) +
  geom_vline(xintercept = true_phi, colour = "red", linetype = "dotdash") +
#  geom_vline(xintercept = mod_phi, colour = "green3", linetype = "dotdash") +
  geom_vline(xintercept = mean(MI_phi[,"Value"]), colour = "gold") +
  geom_vline(xintercept = quantile(MI_phi[,"Value"], 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(MI_phi[,"Value"], 0.975),
             linetype="dashed") +
  labs(x = expression(phi), y = "%") +
  #  labs(title = expression(paste("Histogram - ",phi))) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20)) +
  xlim(minimo, maximo)

# Salvar 4 x 6 in. - Landscape
ggarrange(tr_MA_phi, tr_MI_phi, ncol = 2, nrow = 1)
ggarrange(hist_MA_phi, hist_MI_phi, ncol = 2, nrow = 1)
