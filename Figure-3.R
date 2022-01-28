## A mirror plot.

library(xcms)

std_data <- readMSData("data/mzML/HighIS_Mix07_CE20_POS.mzML",
                       mode = "onDisk")

cwp <- CentWaveParam(snthresh = 10, prefilter = c(3, 4000), ppm = 40,
                     peakwidth = c(2, 8), integrate = 2)
std_data <- findChromPeaks(std_data, param = cwp)
#' Peak refinement
std_data <- refineChromPeaks(
    std_data, MergeNeighboringPeaksParam(expandRt = 3))
pks <- data.frame(peak_id = rownames(chromPeaks(std_data)),
                  chromPeaks(std_data))

std_ions <- read.table("data/std_ions.txt", sep = "\t", header = TRUE)

library(CompoundDb)

#' Load a CompDb database with compound annotation from HMDB
cdb <- CompDb("data/CompDb.Hsapiens.HMDB.5.0.sqlite")
idb <- IonDb(tempfile(), cdb)

#' Insert measured m/z and retention times for ions
idb <- insertIon(idb, std_ions)

library(MetaboAnnotation)
param <- MzRtParam(ppm = 10, toleranceRt = 7)
pks_match <- matchMz(
    pks, ions(idb, c("compound_id", "ion_adduct", "ion_mz", "ion_rt", "name")),
    param = param, mzColname = c("mz", "ion_mz"),
    rtColname = c("rt", "ion_rt"))
pks_match

pks_match <- pks_match[whichQuery(pks_match)]



std_spectra <- chromPeakSpectra(std_data, return.type = "Spectra",
                                peaks = pks_match$peak_id)

#' Define a function to remove low intensity peaks
low_int <- function(x, ...) {
    x > max(x, na.rm = TRUE) * 0.05
}
#' Remove peaks with an intensity below 5% of BPI
std_spectra <- filterIntensity(std_spectra, intensity = low_int)

#' Remove peaks with less than 3 peaks
std_spectra <- std_spectra[lengths(std_spectra) > 2]
std_spectra

#' Define a function to *scale* the intensities
scale_int <- function(x, ...) {
    maxint <- max(x[, "intensity"], na.rm = TRUE)
    x[, "intensity"] <- 100 * x[, "intensity"] / maxint
    x
}
#' *Apply* the function to the data
std_spectra <- addProcessing(std_spectra, scale_int)

spectra_match <- matchSpectra(
    std_spectra, Spectra(cdb),
    param = CompareSpectraParam(ppm = 50, requirePrecursor = FALSE))
spectra_match

caf <- pruneTarget(spectra_match[1])
caf@matches <- caf@matches[1, ]

png("Figure-3.png", width = 8, height = 6, units = "cm",
    res = 600, pointsize = 6)
par(mar = c(4.5, 4.5, 1, 1.5), cex.lab = 2, cex.axis = 1.3, lwd = 1.5)
plotSpectraMirror(caf)
text(x = 100, y = -95, caf$target_compound_id, cex = 2)
dev.off()