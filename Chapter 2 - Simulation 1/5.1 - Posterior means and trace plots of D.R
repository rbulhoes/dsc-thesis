##### Inform the directory 
setwd("D:/Exportações/Chap2/Sim1")
getwd()

##### Some packages
library(ggplot2)
library(ggpubr)
# library(grid)
# library(gridExtra)

smp_size <- 1000
N <- 17

D <- read.table(file = "TamanhoDez/D.txt", head = FALSE, sep = ",")
colnames(D) <- c("Site", "Longitude", "Latitude")
S <- read.table("TamanhoDez/S.txt", head = FALSE, sep = ",")
colnames(S) <- c("Site", "Longitude", "Latitude")
e_D_TamanhoDez <- read.table("TamanhoDez/e_D.txt", head = FALSE, sep = ",")
colnames(e_D_TamanhoDez) <- c("Iteration", "Site", "Longitude", "Latitude")
e_D_TamanhoCem <- read.table("TamanhoCem_AF/e_D.txt", head = FALSE, sep = ",")
colnames(e_D_TamanhoCem) <- c("Iteration", "Site", "Longitude", "Latitude")
e_D_TamanhoMil <- read.table("TamanhoMil/e_D.txt", head = FALSE, sep = ",")
colnames(e_D_TamanhoMil) <- c("Iteration", "Site", "Longitude", "Latitude")

##### Geographic region
library(ggrepel) # Para colocar r?tulos/etiquetas/nomes das esta??es nos mapas
Sites <- 1:N
Category <- factor(x = c(rep(x = 1, times = 2), rep(x = 2, times = N - 2)),
                   levels = c(1, 2), labels = c("Anchor", "Non-anchor"))
GR <- data.frame(Sites, Category, S[1:N, 2:3])
## Salvar 5 x 7 in - Landscape, e usar no LaTeX com scale = 0.6.
ggplot() + 
  geom_point(data = GR, size = 2, aes(x = Longitude,
                                      y = Latitude,
                                      color = Category)) +
  geom_text_repel(data = GR, hjust = "right",
                  aes(x = Longitude, y = Latitude, label = Sites),
                  size = 3.5) +
  scale_color_manual(breaks = c("Anchor", "Non-anchor"),
                     values=c("red", "green")) +
  labs(title = NULL) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        text = element_text(size = 20)) +
  coord_fixed(ratio = 1)

TamanhoDez_stats_lon <- matrix(data = 0, nrow = N, ncol = 4)
TamanhoDez_stats_lat <- matrix(data = 0, nrow = N, ncol = 4)
colnames(TamanhoDez_stats_lon) <- c("LI", "Mean", "True", "LS")
colnames(TamanhoDez_stats_lat) <- c("LI", "Mean", "True", "LS")
TamanhoCem_stats_lon <- matrix(data = 0, nrow = N, ncol = 4)
TamanhoCem_stats_lat <- matrix(data = 0, nrow = N, ncol = 4)
colnames(TamanhoCem_stats_lon) <- c("LI", "Mean", "True", "LS")
colnames(TamanhoCem_stats_lat) <- c("LI", "Mean", "True", "LS")
TamanhoMil_stats_lon <- matrix(data = 0, nrow = N, ncol = 4)
TamanhoMil_stats_lat <- matrix(data = 0, nrow = N, ncol = 4)
colnames(TamanhoMil_stats_lon) <- c("LI", "Mean", "True", "LS")
colnames(TamanhoMil_stats_lat) <- c("LI", "Mean", "True", "LS")
for(n in 1:N) {
    TamanhoDez_n = subset(x = e_D_TamanhoDez,
                          select = c("Longitude", "Latitude"),
                          subset = Site == n)
    TamanhoDez_stats_lon[n, "LI"] = quantile(x = TamanhoDez_n[,"Longitude"],
                                       probs = 0.025)
    TamanhoDez_stats_lat[n, "LI"] = quantile(x = TamanhoDez_n[,"Latitude"],
                                       probs = 0.025)
    TamanhoDez_stats_lon[n, "Mean"] = mean(x = TamanhoDez_n[,"Longitude"])
    TamanhoDez_stats_lat[n, "Mean"] = mean(x = TamanhoDez_n[,"Latitude"])
    TamanhoDez_stats_lon[n, "True"] = D[n, 2]
    TamanhoDez_stats_lat[n, "True"] = D[n, 3]
    TamanhoDez_stats_lon[n, "LS"] = quantile(x = TamanhoDez_n[,"Longitude"],
                                       probs = 0.975)
    TamanhoDez_stats_lat[n, "LS"] = quantile(x = TamanhoDez_n[,"Latitude"],
                                       probs = 0.975)
    TamanhoCem_n = subset(x = e_D_TamanhoCem,
                          select = c("Longitude", "Latitude"),
                          subset = Site == n)
    TamanhoCem_stats_lon[n, "LI"] = quantile(x = TamanhoCem_n[,"Longitude"],
                                      probs = 0.025)
    TamanhoCem_stats_lat[n, "LI"] = quantile(x = TamanhoCem_n[,"Latitude"],
                                      probs = 0.025)
    TamanhoCem_stats_lon[n, "Mean"] = mean(x = TamanhoCem_n[,"Longitude"])
    TamanhoCem_stats_lat[n, "Mean"] = mean(x = TamanhoCem_n[,"Latitude"])
    TamanhoCem_stats_lon[n, "True"] = D[n, 2]
    TamanhoCem_stats_lat[n, "True"] = D[n, 3]
    TamanhoCem_stats_lon[n, "LS"] = quantile(x = TamanhoCem_n[,"Longitude"],
                                      probs = 0.975)
    TamanhoCem_stats_lat[n, "LS"] = quantile(x = TamanhoCem_n[,"Latitude"],
                                      probs = 0.975)
    TamanhoMil_n = subset(x = e_D_TamanhoMil,
                          select = c("Longitude", "Latitude"),
                          subset = Site == n)
    TamanhoMil_stats_lon[n, "LI"] = quantile(x = TamanhoMil_n[,"Longitude"],
                                             probs = 0.025)
    TamanhoMil_stats_lat[n, "LI"] = quantile(x = TamanhoMil_n[,"Latitude"],
                                             probs = 0.025)
    TamanhoMil_stats_lon[n, "Mean"] = mean(x = TamanhoMil_n[,"Longitude"])
    TamanhoMil_stats_lat[n, "Mean"] = mean(x = TamanhoMil_n[,"Latitude"])
    TamanhoMil_stats_lon[n, "True"] = D[n, 2]
    TamanhoMil_stats_lat[n, "True"] = D[n, 3]
    TamanhoMil_stats_lon[n, "LS"] = quantile(x = TamanhoMil_n[,"Longitude"],
                                             probs = 0.975)
    TamanhoMil_stats_lat[n, "LS"] = quantile(x = TamanhoMil_n[,"Latitude"],
                                             probs = 0.975)
}

TamanhoDez_stats = data.frame(TamanhoDez_stats_lon, TamanhoDez_stats_lat)
colnames(TamanhoDez_stats) = c("Lon-LI", "Lon-Mean", "Lon-True", "Lon-LS",
                         "Lat-LI", "Lat-Mean", "Lat-True", "Lat-LS")
TamanhoCem_stats = data.frame(TamanhoCem_stats_lon, TamanhoCem_stats_lat)
colnames(TamanhoCem_stats) = c("Lon-LI", "Lon-Mean", "Lon-True", "Lon-LS",
                        "Lat-LI", "Lat-Mean", "Lat-True", "Lat-LS")
TamanhoMil_stats = data.frame(TamanhoMil_stats_lon, TamanhoMil_stats_lat)
colnames(TamanhoMil_stats) = c("Lon-LI", "Lon-Mean", "Lon-True", "Lon-LS",
                               "Lat-LI", "Lat-Mean", "Lat-True", "Lat-LS")

n <- 17 # Use 03, 04, ..., 09, 10, ..., 17
central <- "C:/Users/rodri/OneDrive/Tese/Figures/Chap2/Sim1_tr_D"
string <- ifelse(test = n < 10,
                 yes = sprintf("%02d", n), 
                 no = as.character(n))
formato <- ".pdf"

TamanhoDez_n = subset(x = e_D_TamanhoDez,
                      select = c("Iteration", "Longitude", "Latitude"),
                      subset = Site == n)
TamanhoCem_n = subset(x = e_D_TamanhoCem,
                      select = c("Iteration", "Longitude", "Latitude"),
                      subset = Site == n)
TamanhoMil_n = subset(x = e_D_TamanhoMil,
                      select = c("Iteration", "Longitude", "Latitude"),
                      subset = Site == n)

### Trace plots
tr_TamanhoDez_lon_n <- ggplot(data = TamanhoDez_n, 
                              aes(x = Iteration, y = Longitude)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = D[n,2], colour = "purple") +
  geom_hline(yintercept = TamanhoDez_stats_lon[n,1], linetype = "dashed") +
  geom_hline(yintercept = TamanhoDez_stats_lon[n,4], linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value") 

tr_TamanhoCem_lon_n <- ggplot(data = TamanhoCem_n, 
                              aes(x = Iteration, y = Longitude)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = D[n,2], colour = "purple") +
  geom_hline(yintercept = TamanhoCem_stats_lon[n,1], linetype = "dashed") +
  geom_hline(yintercept = TamanhoCem_stats_lon[n,4], linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value") 

tr_TamanhoMil_lon_n <- ggplot(data = TamanhoMil_n, 
                              aes(x = Iteration, y = Longitude)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = D[n,2], colour = "purple") +
  geom_hline(yintercept = TamanhoMil_stats_lon[n,1], linetype = "dashed") +
  geom_hline(yintercept = TamanhoMil_stats_lon[n,4], linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value") 

tr_TamanhoDez_lat_n <- ggplot(data = TamanhoDez_n, 
                              aes(x = Iteration, y = Latitude)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = D[n,3], colour = "purple") +
  geom_hline(yintercept = TamanhoDez_stats_lat[n,1], linetype = "dashed") +
  geom_hline(yintercept = TamanhoDez_stats_lat[n,4], linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value") 

tr_TamanhoCem_lat_n <- ggplot(data = TamanhoCem_n, 
                              aes(x = Iteration, y = Latitude)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = D[n,3], colour = "purple") +
  geom_hline(yintercept = TamanhoCem_stats_lat[n,1], linetype = "dashed") +
  geom_hline(yintercept = TamanhoCem_stats_lat[n,4], linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value") 

tr_TamanhoMil_lat_n <- ggplot(data = TamanhoMil_n, 
                              aes(x = Iteration, y = Latitude)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = D[n,3], colour = "purple") +
  geom_hline(yintercept = TamanhoMil_stats_lat[n,1], linetype = "dashed") +
  geom_hline(yintercept = TamanhoMil_stats_lat[n,4], linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value") 

ggarrange(tr_TamanhoDez_lon_n, tr_TamanhoDez_lat_n,
          tr_TamanhoCem_lon_n, tr_TamanhoCem_lat_n,
          tr_TamanhoMil_lon_n, tr_TamanhoMil_lat_n,
          ncol = 2, nrow = 3)
  
pdf(file = paste(central, string, "_1_T10", formato, sep = ""), 
    # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4) # The height of the plot in inches
tr_TamanhoDez_lon_n
dev.off()

pdf(file = paste(central, string, "_1_T100", formato, sep = ""), 
    # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4) # The height of the plot in inches
tr_TamanhoCem_lon_n
dev.off()

pdf(file = paste(central, string, "_1_T1000", formato, sep = ""), 
    # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4) # The height of the plot in inches
tr_TamanhoMil_lon_n
dev.off()

pdf(file = paste(central, string, "_2_T10", formato, sep = ""), 
    # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4) # The height of the plot in inches
tr_TamanhoDez_lat_n
dev.off()

pdf(file = paste(central, string, "_2_T100", formato, sep = ""), 
    # The directory you want to save the file in
    width = 6, # The width of the plot in inches.
    height = 4) # The height of the plot in inches
tr_TamanhoCem_lat_n
dev.off()

pdf(file = paste(central, string, "_2_T1000", formato, sep = ""), 
    # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4) # The height of the plot in inches
tr_TamanhoMil_lat_n
dev.off()

