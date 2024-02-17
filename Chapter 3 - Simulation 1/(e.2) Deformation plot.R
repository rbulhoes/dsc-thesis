library(Morpho)

deformation <- function(matrix, tarmatrix, ngrid = 0, lwd = 1, bg_S, bg_D,
                        show = c(1:2), lines = TRUE, lcol = 1, pch_S = 19,
                        col1 = 2, col2 = 3, pcaxis = FALSE, pch_D = 19,
                        add = FALSE, wireframe = NULL, margin = 0.2,
                        gridcol = "grey", cex1 = 1, cex2 = 1, linecol, ...) 
{
  k <- dim(matrix)[1]
  x0 <- NULL
  if (ngrid > 1) {
    x2 <- x1 <- c(0:(ngrid - 1)) / ngrid
    x0 <- as.matrix(expand.grid(x1,
                                x2))
    xrange <- diff(range(matrix[, 1]))
    yrange <- diff(range(matrix[, 2]))
    xrange1 <- diff(range(tarmatrix[, 1]))
    yrange1 <- diff(range(tarmatrix[, 2]))
    mean.mat <- c((min(matrix[, 1]) + max(matrix[, 1])) / 2, 
                  (min(matrix[, 2]) + max(matrix[, 2])) / 2)
    cent.mat <- scale(matrix,
                      scale = FALSE,
                      center = mean.mat)
    maxi <- max(c(xrange,
                  yrange,
                  xrange1,
                  yrange1))
    maxi <- (1 + margin) * maxi
    x0 <- maxi * x0
    x0 <- scale(x0,
                scale = FALSE)
    x0[, 2] <- x0[, 2]
    if (pcaxis) 
      space <- eigen(crossprod(cent.mat))$vectors
    else space <- diag(2)
    x0 <- (x0 %*% space)
    x00 <- x0 <- scale(x0,
                       center = -mean.mat,
                       scale = F)
    x0 <- tps3d(x0,
                matrix,
                tarmatrix,
                threads = 1)
  }
  lims <- apply(rbind(matrix,
                      tarmatrix,
                      x0),
                2,
                range)
  if (1 %in% show) {
    if (add) 
      points(matrix,
             col = col1,
             cex = cex1,
             pch = pch_S,
             bg = bg_S)
    else plot(matrix,
              col = col1,
              xlim = lims[, 1],
              ylim = lims[, 2],
              asp = 0,
              pch = pch_S,
              bg = bg_S,
              xlab = "First coordinate",
              ylab = "Second coordinate",
              axes = TRUE, 
              cex = cex1)
    if (!is.null(wireframe)) 
      lineplot(matrix,
               wireframe,
               col = col1,
               pch = pch_S,
               bg = bg_S,
               lwd = lwd)
  }
  if (2 %in% show) {
    if (1 %in% show || add) 
      points(tarmatrix,
             col = col2,
             cex = cex2,
             pch = pch_D,
             bg = bg_D)
    else plot(tarmatrix,
              col = col2,
              xlim = lims[, 1],
              ylim = lims[, 2],
              asp = 0,
              xlab = "First coordinate",
              ylab = "Second coordinate",
              axes = T, 
              pch = pch_D,
              bg = bg_D,
              cex = cex2)
    if (!is.null(wireframe)) 
      lineplot(tarmatrix,
               wireframe,
               col = col2,
               pch = pch_D,
               bg = bg_D,
               lwd = lwd)
  }
  if (lines) {
    linemesh <- list()
    linemesh$vb <- rbind(matrix,
                         tarmatrix)
    linemesh$it <- cbind(1:k, (1:k) + k)
    for (i in 1:nrow(linemesh$it))
      lines(linemesh$vb[linemesh$it[i, ], ],
            lwd = lwd,
            col = linecol,
            lty = 2)
  }
  if (ngrid > 1) {
    myrange <- 0:(ngrid - 1)
    for (i in 0:(ngrid - 1)) {
      lines(x0[(1:ngrid) + (i * ngrid), ],
            col = gridcol)
      lines(x0[(myrange * ngrid) + i + 1, ],
            col = gridcol)
    }
  }
}

# Salvar 4 x 6 in. - Landscape (Sim1_Def_True e Sim1_Def_T10)
par(mar=c(4,4,0.1,0.1), cex = 1.2)
deformation(matrix = as.matrix(S[,2:3]), # regi?o geogr?fica
            tarmatrix = as.matrix(D[,2:3]),
            ngrid = 20,
            pch_S = 25,
            pch_D = 24,
            bg_S = "goldenrod3",
            bg_D = "darkolivegreen1",
            col1 = "black", # col1 - de matrix, col2 - de tarmatrix
            col2 = "black",
            pcaxis = FALSE,
            margin = 0.225,
            lwd = 0, 
            show = c(1:2),
            lines = TRUE,
            linecol = "black",
            lcol = 1, 
            add = FALSE,
            wireframe = NULL,
            gridcol = "lightgrey",
            cex1 = 1.0,
            cex2 = 0.9)
legend("topright", inset = .07, pch = c(25, 24), horiz = TRUE,
       cex = c(0.75, 0.75), pt.bg = c("goldenrod3", "darkolivegreen1"),
       legend = c(as.expression(bquote(italic(S))),
                  as.expression(bquote(italic(D)))))

deformation(matrix = as.matrix(S[,2:3]), # regi?o geogr?fica
            tarmatrix = as.matrix(MV1_stats[1:N, c("Lon-Mean", "Lat-Mean")]),
            ngrid = 20,
            pch_S = 25,
            pch_D = 24,
            bg_S = "goldenrod3",
            bg_D = "darkolivegreen1",
            col1 = "black", # col1 - de matrix, col2 - de tarmatrix
            col2 = "black",
            pcaxis = FALSE,
            margin = 0.225,
            lwd = 0, 
            show = c(1:2),
            lines = TRUE,
            linecol = "black",
            lcol = 1, 
            add = FALSE,
            wireframe = NULL,
            gridcol = "lightgrey",
            cex1 = 1.0,
            cex2 = 0.9)
legend("topright", inset = .07, pch = c(25, 24), horiz = TRUE,
       cex = c(0.75, 0.75), pt.bg = c("goldenrod3", "darkolivegreen1"),
       legend = c(as.expression(bquote(italic(S))),
                  as.expression(bquote(italic(D)))))

deformation(matrix = as.matrix(S[,2:3]), # regi?o geogr?fica
            tarmatrix = as.matrix(MV2_stats[1:N, c("Lon-Mean", "Lat-Mean")]),
            ngrid = 20,
            pch_S = 25,
            pch_D = 24,
            bg_S = "goldenrod3",
            bg_D = "darkolivegreen1",
            col1 = "black", # col1 - de matrix, col2 - de tarmatrix
            col2 = "black",
            pcaxis = FALSE,
            margin = 0.225,
            lwd = 0, 
            show = c(1:2),
            lines = TRUE,
            linecol = "black",
            lcol = 1, 
            add = FALSE,
            wireframe = NULL,
            gridcol = "lightgrey",
            cex1 = 1.0,
            cex2 = 0.9)
legend("topright", inset = .07, pch = c(25, 24), horiz = TRUE,
       cex = c(0.75, 0.75), pt.bg = c("goldenrod3", "darkolivegreen1"),
       legend = c(as.expression(bquote(italic(S))),
                  as.expression(bquote(italic(D)))))

deformation(matrix = as.matrix(S[,2:3]), # regi?o geogr?fica
            tarmatrix = as.matrix(MV4_stats[1:N, c("Lon-Mean", "Lat-Mean")]),
            ngrid = 20,
            pch_S = 25,
            pch_D = 24,
            bg_S = "goldenrod3",
            bg_D = "darkolivegreen1",
            col1 = "black", # col1 - de matrix, col2 - de tarmatrix
            col2 = "black",
            pcaxis = FALSE,
            margin = 0.225,
            lwd = 0, 
            show = c(1:2),
            lines = TRUE,
            linecol = "black",
            lcol = 1, 
            add = FALSE,
            wireframe = NULL,
            gridcol = "lightgrey",
            cex1 = 1.0,
            cex2 = 0.9)
legend("topright", inset = .07, pch = c(25, 24), horiz = TRUE,
       cex = c(0.75, 0.75), pt.bg = c("goldenrod3", "darkolivegreen1"),
       legend = c(as.expression(bquote(italic(S))),
                  as.expression(bquote(italic(D)))))

dev.off()
