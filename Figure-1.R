library(Spectra)
library(MsBackendMgf)


sps <- Spectra("data/caffeine.mgf", source = MsBackendMgf())
Q <- sps[1]
T <- sps[2]

cl <- "#FFFF33"
expand <- 2

png("Figure-1.png", width = 12, height = 4, res = 600, pointsize = 6,
    units = "cm")
par(mfrow = c(1, 3), mar = c(4, 4.3, 1, 0.5), cex.lab = 2, cex.axis = 1.6)

plot_mirror <- function() {
    plotSpectraMirror(Q, T, ppm = 0, lwd = 1)
    text(80, 100, labels = "query", pos = 1, cex = 2)
    text(80, -100, labels = "target", pos = 3, cex = 2)
}

## Left
plot_mirror()
rect(mz(Q)[[1]] - expand,
     rep(0, length(Q)),
     mz(Q)[[1]] + expand,
     intensity(Q)[[1]] + 2 * expand,
     col = paste0(cl, 20), border = cl)
idx <- c(3, 4, 5)
rect(mz(T)[[1]][idx] - expand,
     rep(0, length(idx)),
     mz(T)[[1]][idx] + expand,
     - intensity(T)[[1]][idx] - 2 * expand,
     col = paste0(cl, 20), border = cl)

## Right
plot_mirror()
idx <- c(3, 4, 5)
rect(mz(Q)[[1]][idx] - expand,
     rep(0, length(idx)),
     mz(Q)[[1]][idx] + expand,
     intensity(Q)[[1]][idx] + 2 * expand,
     col = paste0(cl, 20), border = cl)
idx <- 1:5
rect(mz(T)[[1]][idx] - expand,
     rep(0, length(idx)),
     mz(T)[[1]][idx] + expand,
     - intensity(T)[[1]][idx] - 2 * expand,
     col = paste0(cl, 20), border = cl)

## Inner
plot_mirror()
idx <- c(3, 4, 5)
rect(mz(Q)[[1]][idx] - expand,
     rep(0, length(idx)),
     mz(Q)[[1]][idx] + expand,
     intensity(Q)[[1]][idx] + 2 * expand,
     col = paste0(cl, 20), border = cl)
rect(mz(T)[[1]][idx] - expand,
     rep(0, length(idx)),
     mz(T)[[1]][idx] + expand,
     - intensity(T)[[1]][idx] - 2 * expand,
     col = paste0(cl, 20), border = cl)
dev.off()
