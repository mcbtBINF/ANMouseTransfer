# Calculating effect size
rm(list = ls())

library("stringr")
#library("effsize")

# From https://stackoverflow.com/questions/15436702/estimate-cohens-d-for-effect-size
cohens_d <- function(x, y) {
  lx <- length(x)- 1
  ly <- length(y)- 1
  md  <- abs(mean(x) - mean(y))        ## mean difference (numerator)
  csd <- lx * var(x) + ly * var(y)
  csd <- csd/(lx + ly)
  csd <- sqrt(csd)                     ## common sd computation
  
  cd  <- md/csd                        ## cohen's d
}

if (.Platform$OS.type == "windows") {
  baseDir <- c("g:/My Drive/Datasets/CarrollData/Carroll_MouseTransfer/")
  #baseDir <- c("/Google Drive/Datasets/CarrollData/Carroll_MouseTransfer/")
} else {
  #baseDir <- c("/Users/mbrown67/Google Drive/Datasets/CarrollData/Carroll_MouseTransfer/")
  baseDir <- c("My Drive/Datasets/CarrollData/Carroll_MouseTransfer/")  
}

processedDir <- paste0(baseDir, "processed/")
analysisDir <- paste0(baseDir, "analysis/")
metadataDir <- paste0(baseDir, "metadata/")

taxonomicLevels <- c("phylum", "class", "order", "family", "genus")

setwd(processedDir)
fileName <- paste0("genus", "_LogNormwithMetadata.txt")
myT <- read.table(fileName, sep = "\t", header=TRUE, stringsAsFactors = FALSE, na.strings = "n.a.")
setwd(analysisDir)
abundanceIndices <- 51:ncol(myT)
# Only on the forward reads
myT <- myT[grep("_1", myT$comboID),]
# Drop samples with reads < 1000
myT <- myT[myT$depthAtLevel >= 1000, ]

myT.Humans <- myT[grep("donor", myT$Sample.ID), ]
myT.Humans$Mouse.group <- sapply(sapply(strsplit(myT.Humans$Sample.ID, "\\."),"[[",1), function(x) str_sub(x, -2, -1))
humanGroupMeans <- aggregate(myT.Humans[,abundanceIndices], list(myT.Humans$Mouse.group), mean)
HCvT1effectSizes <- vector()
index <- 1
for (aIter in abundanceIndices) {
  HCvT1effectSizes[index] <- cohens_d(myT.Humans[myT.Humans$Mouse.group == "HC", aIter], myT.Humans[myT.Humans$Mouse.group == "T1", aIter])
  index = index + 1
}
HCvT1effectSizes.mean <- mean(HCvT1effectSizes, na.rm = TRUE)

HCvT2effectSizes <- vector()
index <- 1
for (aIter in abundanceIndices) {
  HCvT2effectSizes[index] <- cohens_d(myT.Humans[myT.Humans$Mouse.group == "HC", aIter], myT.Humans[myT.Humans$Mouse.group == "T2", aIter])
  index = index + 1
}
HCvT2effectSizes.mean <- mean(HCvT2effectSizes, na.rm = TRUE)

T1vT2effectSizes <- vector()
index <- 1
for (aIter in abundanceIndices) {
  T1vT2effectSizes[index] <- cohens_d(myT.Humans[myT.Humans$Mouse.group == "T1", aIter], myT.Humans[myT.Humans$Mouse.group == "T2", aIter])
  index = index + 1
}
T1vT2effectSizes.mean <- mean(T1vT2effectSizes, na.rm = TRUE)

myT.Mice <- myT[myT$Mouse.group %in% c("HC", "T1", "T2"),]
myT.Mice <- myT.Mice[myT.Mice$Week == 4,]
mouseGroupMeans <- aggregate(myT.Mice[,abundanceIndices], list(myT.Mice$Mouse.group), mean)

HCvT1effectSizesMouse <- vector()
index <- 1
for (aIter in abundanceIndices) {
  HCvT1effectSizesMouse[index] <- cohens_d(myT.Mice[myT.Mice$Mouse.group == "HC", aIter], myT.Mice[myT.Mice$Mouse.group == "T1", aIter])
  index = index + 1
}
HCvT1effectSizesMouse.mean <- mean(HCvT1effectSizesMouse, na.rm = TRUE)

HCvT2effectSizesMouse <- vector()
index <- 1
for (aIter in abundanceIndices) {
  HCvT2effectSizesMouse[index] <- cohens_d(myT.Mice[myT.Mice$Mouse.group == "HC", aIter], myT.Mice[myT.Mice$Mouse.group == "T2", aIter])
  index = index + 1
}
HCvT2effectSizesMouse.mean <- mean(HCvT2effectSizesMouse, na.rm = TRUE)

T1vT2effectSizesMouse <- vector()
index <- 1
for (aIter in abundanceIndices) {
  T1vT2effectSizesMouse[index] <- cohens_d(myT.Mice[myT.Mice$Mouse.group == "T1", aIter], myT.Mice[myT.Mice$Mouse.group == "T2", aIter])
  index = index + 1
}
T1vT2effectSizesMouse.mean <- mean(T1vT2effectSizesMouse, na.rm = TRUE)
