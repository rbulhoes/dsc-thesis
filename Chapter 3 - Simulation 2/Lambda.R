Lambda <- matrix(data = c(2.23606798,  -1.78885438, 0, 1.34164079),
                 byrow = TRUE, nrow = 2)

vec11 <- NULL; vec12 <- NULL; vec21 <- NULL; vec22 <- NULL
n <- 3; n_ <- 5
for (k in 1:smp_size) {
  d1nk <- subset(x = e_D_MA, subset = Iteration == k)[n,3]
  d2nk <- subset(x = e_D_MA, subset = Iteration == k)[n,4]
  d1n_k <- subset(x = e_D_MA, subset = Iteration == k)[n_,3]
  d2n_k <- subset(x = e_D_MA, subset = Iteration == k)[n_,4]
  aux <- S[n,2]*S[n_,3] - S[n_,2]*S[n,3]
  now11 <- (d1nk*S[n_,3] - d1n_k*S[n,3]) / aux
  now12 <- (S[n,2]*d1n_k - S[n_,2]*d1nk) / aux
  now21 <- (d2nk*S[n_,3] - d2n_k*S[n,3]) / aux
  now22 <- (S[n,2]*d2n_k - S[n_,2]*d2nk) / aux
  vec11 <- c(vec11, now11)
  vec12 <- c(vec12, now12)
  vec21 <- c(vec21, now21)
  vec22 <- c(vec22, now22)
}

c(quantile(vec11, c(0.025, 0.5, 0.975)), Lambda[1,1])
c(quantile(vec12, c(0.025, 0.5, 0.975)), Lambda[1,2])
c(quantile(vec21, c(0.025, 0.5, 0.975)), Lambda[2,1])
c(quantile(vec22, c(0.025, 0.5, 0.975)), Lambda[2,2])

vec11 <- NULL; vec12 <- NULL; vec21 <- NULL; vec22 <- NULL
for (n in 3:(N-1)) {
  for (n_ in (n+1):N) {
    aux <- S[n,2]*S[n_,3] - S[n_,2]*S[n,3]
    if (aux != 0) {
      for (k in 1:smp_size) {
        #        sub <- subset(x = e_D_MA, subset = Iteration == k & Site %in% c(n, n_),
        #                      select = c("Longitude", "Latitude"))
        sub <- e_D_MA[which(e_D_MA$Iteration == k & e_D_MA$Site %in% c(n, n_)),3:4]
        d1nk <- sub[1, "Longitude"]
        d2nk <- sub[1, "Latitude"]
        d1n_k <- sub[2, "Longitude"]
        d2n_k <- sub[2, "Latitude"]
        now11 <- (d1nk*S[n_,3] - d1n_k*S[n,3]) / aux
        now12 <- (S[n,2]*d1n_k - S[n_,2]*d1nk) / aux
        now21 <- (d2nk*S[n_,3] - d2n_k*S[n,3]) / aux
        now22 <- (S[n,2]*d2n_k - S[n_,2]*d2nk) / aux
        vec11 <- c(vec11, now11)
        vec12 <- c(vec12, now12)
        vec21 <- c(vec21, now21)
        vec22 <- c(vec22, now22)
        print(c(n, n_, k))
      }
    } else {
      invisible()
    }
  }
}

c(quantile(vec11, c(0.025, 0.5, 0.975)), Lambda[1,1])
c(quantile(vec12, c(0.025, 0.5, 0.975)), Lambda[1,2])
c(quantile(vec21, c(0.025, 0.5, 0.975)), Lambda[2,1])
c(quantile(vec22, c(0.025, 0.5, 0.975)), Lambda[2,2])

data11 <- data.frame(1:length(vec11), vec11)
colnames(data11) <- c("Sample", "Value")
data12 <- data.frame(1:length(vec12), vec12)
colnames(data12) <- c("Sample", "Value")
data21 <- data.frame(1:length(vec21), vec21)
colnames(data21) <- c("Sample", "Value")
data22 <- data.frame(1:length(vec22), vec22)
colnames(data22) <- c("Sample", "Value")

hist_Lamb11 <- ggplot(data11, aes(x = Value)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = Value,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 15) +
  geom_vline(xintercept = Lambda[1,1], colour = "blue", linetype = "dotdash") +
  geom_vline(xintercept = mean(vec11), colour = "gold") +
  geom_vline(xintercept = quantile(vec11, 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(vec11, 0.975),
             linetype="dashed") +
  labs(x = expression(Lambda["1,1"]), y = "%") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

hist_Lamb12 <- ggplot(data12, aes(x = Value)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = Value,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 15) +
  geom_vline(xintercept = Lambda[1,2], colour = "blue", linetype = "dotdash") +
  geom_vline(xintercept = mean(vec12), colour = "gold") +
  geom_vline(xintercept = quantile(vec12, 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(vec12, 0.975),
             linetype="dashed") +
  labs(x = expression(Lambda["1,2"]), y = "%") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

hist_Lamb21 <- ggplot(data21, aes(x = Value)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = Value,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 15) +
  geom_vline(xintercept = Lambda[2,1], colour = "blue", linetype = "dotdash") +
  geom_vline(xintercept = mean(vec21), colour = "gold") +
  geom_vline(xintercept = quantile(vec21, 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(vec21, 0.975),
             linetype="dashed") +
  labs(x = expression(Lambda["2,1"]), y = "%") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

hist_Lamb22 <- ggplot(data22, aes(x = Value)) + 
  geom_histogram(color = "gray50",
                 fill = "snow1",
                 mapping=aes(x = Value,
                             y = after_stat(count)/sum(after_stat(count))*100),
                 bins = 15) +
  geom_vline(xintercept = Lambda[2,2], colour = "blue", linetype = "dotdash") +
  geom_vline(xintercept = mean(vec22), colour = "gold") +
  geom_vline(xintercept = quantile(vec22, 0.025),
             linetype="dashed") +
  geom_vline(xintercept = quantile(vec22, 0.975),
             linetype="dashed") +
  labs(x = expression(Lambda["2,2"]), y = "%") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))

# 4 x 6 Landscape
ggarrange(hist_Lamb11, hist_Lamb12, hist_Lamb21, hist_Lamb22,
          ncol = 2, nrow = 2)