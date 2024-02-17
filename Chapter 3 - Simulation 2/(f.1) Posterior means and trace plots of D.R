##### Inform the directory 
setwd("D:/Exportações/Chap3/Sim2")
getwd()

##### Some packages
library(ggplot2)
library(ggpubr)
# library(grid)
# library(gridExtra)

smp_size <- 1000
N <- 16
N_i <- 3
N_tot <- N + N_i

D = read.table("MA/D.txt",
               head = FALSE,
               sep = ",")
colnames(D) = c("Site", "Longitude", "Latitude")
D_i = read.table("MA/D_i.txt",
                 head = FALSE,
                 sep = ",")
colnames(D_i) = c("Site", "Longitude", "Latitude")
S = read.table("MA/S.txt",
               head = FALSE,
               sep = ",")
colnames(S) = c("Site", "Longitude", "Latitude")
S_i = read.table("MA/S_i.txt",
                 head = FALSE,
                 sep = ",")
colnames(S_i) = c("Site", "Longitude", "Latitude")
e_D_MA = read.table("MA/e_D.txt",
                    head = FALSE,
                    sep = ",")
colnames(e_D_MA) = c("Iteration", "Site", "Longitude", "Latitude")

MA_stats_lon = matrix(data = 0, nrow = N, ncol = 3)
MA_stats_lat = matrix(data = 0, nrow = N, ncol = 3)
colnames(MA_stats_lon) = c("LI", "Mean", "LS")
colnames(MA_stats_lat) = c("LI", "Mean", "LS")
for(n in 1:N) {
  MA_n = subset(x = e_D_MA,
                select = c("Longitude", "Latitude"),
                subset = Site == n)
  MA_stats_lon[n, "LI"] = quantile(x = MA_n[,"Longitude"],
                                   probs = 0.025)
  MA_stats_lat[n, "LI"] = quantile(x = MA_n[,"Latitude"],
                                   probs = 0.025)
  MA_stats_lon[n, "Mean"] = mean(x = MA_n[,"Longitude"])
  MA_stats_lat[n, "Mean"] = mean(x = MA_n[,"Latitude"])
  MA_stats_lon[n, "LS"] = quantile(x = MA_n[,"Longitude"],
                                   probs = 0.975)
  MA_stats_lat[n, "LS"] = quantile(x = MA_n[,"Latitude"],
                                   probs = 0.975)
}

MA_stats = data.frame(MA_stats_lon, MA_stats_lat)
colnames(MA_stats) = c("Lon-LI", "Lon-Mean", "Lon-LS",
                       "Lat-LI", "Lat-Mean", "Lat-LS")

##### Geographic region
library(ggrepel) # Para colocar r?tulos/etiquetas/nomes das esta??es nos mapas
Sites = 1:N_tot
Category = factor(x = c(rep(1, 2),
                        rep(2, N - 2),
                        rep(3, N_i)),
                  levels = 1:3,
                  labels = c("Anchor", "Non-anchor", "Interpolation"))
GR = data.frame(Sites, Category, rbind(S[1:N, 2:3], S_i[1:N_i, 2:3]))
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

n <- 16 # Use 03, 04, ..., 09, 10, ..., 16
central <- "C:/Users/rodri/OneDrive/Tese/Figures/Chap3/Sim2_tr_D"
string <- ifelse(test = n < 10,
                 yes = sprintf("%02d", n), 
                 no = as.character(n))
formato <- ".pdf"
MA_n <- subset(x = e_D_MA,
               select = c("Iteration", "Longitude", "Latitude"),
               subset = Site == n)


### Trace plots
tr_MA_lon_n <- ggplot(data = MA_n, 
                      aes(x = Iteration, y = Longitude)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = D[n,2], colour = "purple") +
  geom_hline(yintercept = MA_stats_lon[n,1], linetype = "dashed") +
  geom_hline(yintercept = MA_stats_lon[n,3], linetype = "dashed") +  
  xlab(label = "Sample") + ylab(label = "Sampled value") 

tr_MA_lat_n <- ggplot(data = MA_n, 
                      aes(x = Iteration, y = Latitude)) + 
  geom_line(colour = "gray60") +
  geom_hline(yintercept = D[n,3], colour = "purple") +
  geom_hline(yintercept = MA_stats_lat[n,1], linetype = "dashed") +
  geom_hline(yintercept = MA_stats_lat[n,3], linetype = "dashed") +  
  xlab(label = "Sample") + ylab(label = "Sampled value") 

ggarrange(tr_MA_lon_n, tr_MA_lat_n,
          ncol = 2, nrow = 1)

pdf(file = paste(central, string, "_1", formato, sep = ""),
    # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4) # The height of the plot in inches
tr_MA_lon_n
dev.off()

pdf(file = paste(central, string, "_2", formato, sep = ""),
    # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4) # The height of the plot in inches
tr_MA_lat_n
dev.off()
