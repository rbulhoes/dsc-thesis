##### Inform the directory 
setwd("D:/Exportações/Chap3/Sim1")
getwd()

##### Some packages
library(ggplot2)
library(ggpubr)
library(coda)
# library(grid)
# library(gridExtra)

smp_size <- 1000
N <- 16

D <- read.table("1MV/D.txt", head = FALSE, sep = ",")
colnames(D) <- c("Site", "Longitude", "Latitude")
S = read.table("1MV/S.txt",
               head = FALSE,
               sep = ",")
colnames(S) = c("Site", "Longitude", "Latitude")
e_D_1MV = read.table("1MV/e_D.txt",
                      head = FALSE,
                      sep = ",")
colnames(e_D_1MV) = c("Iteration", "Site", "Longitude", "Latitude")
e_D_2MV = read.table("2MV/e_D.txt",
                     head = FALSE,
                     sep = ",")
colnames(e_D_2MV) = c("Iteration", "Site", "Longitude", "Latitude")
e_D_4MV = read.table("4MV/e_D.txt",
                            head = FALSE,
                            sep = ",")
colnames(e_D_4MV) = c("Iteration", "Site", "Longitude", "Latitude")

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

MV1_stats_lon = matrix(data = 0, nrow = N, ncol = 4)
MV1_stats_lat = matrix(data = 0, nrow = N, ncol = 4)
colnames(MV1_stats_lon) = c("LI", "Mean", "True", "LS")
colnames(MV1_stats_lat) = c("LI", "Mean", "True", "LS")
MV2_stats_lon = matrix(data = 0, nrow = N, ncol = 4)
MV2_stats_lat = matrix(data = 0, nrow = N, ncol = 4)
colnames(MV2_stats_lon) = c("LI", "Mean", "True", "LS")
colnames(MV2_stats_lat) = c("LI", "Mean", "True", "LS")
MV4_stats_lon = matrix(data = 0, nrow = N, ncol = 4)
MV4_stats_lat = matrix(data = 0, nrow = N, ncol = 4)
colnames(MV4_stats_lon) = c("LI", "Mean", "True", "LS")
colnames(MV4_stats_lat) = c("LI", "Mean", "True", "LS")
for(n in 1:N) {
    MV1_n = subset(x = e_D_1MV,
                          select = c("Longitude", "Latitude"),
                          subset = Site == n)
    MV1_stats_lon[n, "LI"] = quantile(x = MV1_n[,"Longitude"],
                                       probs = 0.025)
    MV1_stats_lat[n, "LI"] = quantile(x = MV1_n[,"Latitude"],
                                       probs = 0.025)
    MV1_stats_lon[n, "Mean"] = mean(x = MV1_n[,"Longitude"])
    MV1_stats_lat[n, "Mean"] = mean(x = MV1_n[,"Latitude"])
    MV1_stats_lon[n, "True"] = D[n, 2]
    MV1_stats_lat[n, "True"] = D[n, 3]
    MV1_stats_lon[n, "LS"] = quantile(x = MV1_n[,"Longitude"],
                                       probs = 0.975)
    MV1_stats_lat[n, "LS"] = quantile(x = MV1_n[,"Latitude"],
                                       probs = 0.975)
    MV2_n = subset(x = e_D_2MV,
                          select = c("Longitude", "Latitude"),
                          subset = Site == n)
    MV2_stats_lon[n, "LI"] = quantile(x = MV2_n[,"Longitude"],
                                      probs = 0.025)
    MV2_stats_lat[n, "LI"] = quantile(x = MV2_n[,"Latitude"],
                                      probs = 0.025)
    MV2_stats_lon[n, "Mean"] = mean(x = MV2_n[,"Longitude"])
    MV2_stats_lat[n, "Mean"] = mean(x = MV2_n[,"Latitude"])
    MV2_stats_lon[n, "True"] = D[n, 2]
    MV2_stats_lat[n, "True"] = D[n, 3]
    MV2_stats_lon[n, "LS"] = quantile(x = MV2_n[,"Longitude"],
                                      probs = 0.975)
    MV2_stats_lat[n, "LS"] = quantile(x = MV2_n[,"Latitude"],
                                      probs = 0.975)
    MV4_n = subset(x = e_D_4MV,
                          select = c("Longitude", "Latitude"),
                          subset = Site == n)
    MV4_stats_lon[n, "LI"] = quantile(x = MV4_n[,"Longitude"],
                                             probs = 0.025)
    MV4_stats_lat[n, "LI"] = quantile(x = MV4_n[,"Latitude"],
                                             probs = 0.025)
    MV4_stats_lon[n, "Mean"] = mean(x = MV4_n[,"Longitude"])
    MV4_stats_lat[n, "Mean"] = mean(x = MV4_n[,"Latitude"])
    MV4_stats_lon[n, "True"] = D[n, 2]
    MV4_stats_lat[n, "True"] = D[n, 3]
    MV4_stats_lon[n, "LS"] = quantile(x = MV4_n[,"Longitude"],
                                             probs = 0.975)
    MV4_stats_lat[n, "LS"] = quantile(x = MV4_n[,"Latitude"],
                                             probs = 0.975)
}

MV1_stats = data.frame(MV1_stats_lon, MV1_stats_lat)
colnames(MV1_stats) = c("Lon-LI", "Lon-Mean", "Lon-True", "Lon-LS",
                         "Lat-LI", "Lat-Mean", "Lat-True", "Lat-LS")
MV2_stats = data.frame(MV2_stats_lon, MV2_stats_lat)
colnames(MV2_stats) = c("Lon-LI", "Lon-Mean", "Lon-True", "Lon-LS",
                        "Lat-LI", "Lat-Mean", "Lat-True", "Lat-LS")
MV4_stats = data.frame(MV4_stats_lon, MV4_stats_lat)
colnames(MV4_stats) = c("Lon-LI", "Lon-Mean", "Lon-True", "Lon-LS",
                               "Lat-LI", "Lat-Mean", "Lat-True", "Lat-LS")

n <- 16 # Use 03, 04, ..., 09, 10, ..., 16
central <- "C:/Users/rodri/OneDrive/Tese/Figures/Chap3/Sim1_tr_D"
string <- ifelse(test = n < 10,
                 yes = sprintf("%02d", n), 
                 no = as.character(n))
formato <- ".pdf"

MV1_n <- subset(x = e_D_1MV,
                select = c("Iteration", "Longitude", "Latitude"),
                subset = Site == n)
MV2_n <- subset(x = e_D_2MV,
                select = c("Iteration", "Longitude", "Latitude"),
                subset = Site == n)
MV4_n <- subset(x = e_D_4MV,
                select = c("Iteration", "Longitude", "Latitude"),
                subset = Site == n)

### Trace plots
tr_MV1_lon_n <- ggplot(data = MV1_n, 
                              aes(x = Iteration, y = Longitude)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = D[n,2], colour = "purple") +
  geom_hline(yintercept = MV1_stats_lon[n,1], linetype = "dashed") +
  geom_hline(yintercept = MV1_stats_lon[n,4], linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value") 

tr_MV2_lon_n <- ggplot(data = MV2_n, 
                              aes(x = Iteration, y = Longitude)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = D[n,2], colour = "purple") +
  geom_hline(yintercept = MV2_stats_lon[n,1], linetype = "dashed") +
  geom_hline(yintercept = MV2_stats_lon[n,4], linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value") 

tr_MV4_lon_n <- ggplot(data = MV4_n, 
                              aes(x = Iteration, y = Longitude)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = D[n,2], colour = "purple") +
  geom_hline(yintercept = MV4_stats_lon[n,1], linetype = "dashed") +
  geom_hline(yintercept = MV4_stats_lon[n,4], linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value") 

tr_MV1_lat_n <- ggplot(data = MV1_n, 
                              aes(x = Iteration, y = Latitude)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = D[n,3], colour = "purple") +
  geom_hline(yintercept = MV1_stats_lat[n,1], linetype = "dashed") +
  geom_hline(yintercept = MV1_stats_lat[n,4], linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value") 

tr_MV2_lat_n <- ggplot(data = MV2_n, 
                              aes(x = Iteration, y = Latitude)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = D[n,3], colour = "purple") +
  geom_hline(yintercept = MV2_stats_lat[n,1], linetype = "dashed") +
  geom_hline(yintercept = MV2_stats_lat[n,4], linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value") 

tr_MV4_lat_n <- ggplot(data = MV4_n, 
                              aes(x = Iteration, y = Latitude)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = D[n,3], colour = "purple") +
  geom_hline(yintercept = MV4_stats_lat[n,1], linetype = "dashed") +
  geom_hline(yintercept = MV4_stats_lat[n,4], linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 15)) +
  xlab(label = "Sample") + ylab(label = "Sampled value") 

ggarrange(tr_MV1_lon_n, tr_MV1_lat_n,
          tr_MV2_lon_n, tr_MV2_lat_n,
          tr_MV4_lon_n, tr_MV4_lat_n,
          ncol = 2, nrow = 3)

pdf(file = paste(central, string, "_1_MV1", formato, sep = ""), 
    # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4) # The height of the plot in inches
tr_MV1_lon_n
dev.off()

pdf(file = paste(central, string, "_1_MV2", formato, sep = ""), 
    # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4) # The height of the plot in inches
tr_MV2_lon_n
dev.off()

pdf(file = paste(central, string, "_1_MV4", formato, sep = ""), 
    # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4) # The height of the plot in inches
tr_MV4_lon_n
dev.off()

pdf(file = paste(central, string, "_2_MV1", formato, sep = ""), 
    # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4) # The height of the plot in inches
tr_MV1_lat_n
dev.off()

pdf(file = paste(central, string, "_2_MV2", formato, sep = ""), 
    # The directory you want to save the file in
    width = 6, # The width of the plot in inches.
    height = 4) # The height of the plot in inches
tr_MV2_lat_n
dev.off()

pdf(file = paste(central, string, "_2_MV4", formato, sep = ""), 
    # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4) # The height of the plot in inches
tr_MV4_lat_n
dev.off()

