# Power Analysis with "pwr" package

rm(list = ls())

library("pwr")
# Assuming simple t-test
# Paired T1 vs T2 with human donors
effectSizeVector <- c(0.10, 0.20, 0.23, 0.42, 0.50, 0.80)
sigLevelVector <- c(0.05, 0.10)

T1vsT2List <- list()
index <- 1

for (esv in effectSizeVector){
  for (slv in sigLevelVector){
    T1vsT2List[[index]] <- pwr.t.test(n=40, d=esv, sig.level = slv, type="paired")
    index = index + 1
  }
}

powermatT1vsT2 <- matrix(sapply(T1vsT2List, FUN = function(x) x$power), nrow = length(sigLevelVector), ncol=length(effectSizeVector))
colnames(powermatT1vsT2) <- effectSizeVector
rownames(powermatT1vsT2) <- sigLevelVector

T1vsT2.power80.sig.05 <- pwr.t.test(n=40, sig=0.05, power = 0.80, type="paired")$d
T1vsT2.power80.sig.10 <- pwr.t.test(n=40, sig=0.10, power = 0.80, type="paired")$d

ANvsHCList <- list()
index <- 1

for (esv in effectSizeVector){
  for (slv in sigLevelVector){
    ANvsHCList[[index]] <- pwr.t2n.test(n1=40, n2=15, d=esv, sig.level = slv)
    index = index + 1
  }
}

powermatANvsHC <- matrix(sapply(ANvsHCList, FUN = function(x) x$power), nrow = length(sigLevelVector), ncol=length(effectSizeVector))
colnames(powermatANvsHC) <- effectSizeVector
rownames(powermatANvsHC) <- sigLevelVector

ANvsHC.power80.sig.05 <- pwr.t2n.test(n1=40, n2=15, sig=0.05, power = 0.80)$d
ANvsHC.power80.sig.10 <- pwr.t2n.test(n1=40, n2=15, sig=0.10, power = 0.80)$d
# T1 vs HC; T2 vs HC for human cohort

# With replication in mouse cohort
T1vsT2MouseList <- list()
index <- 1

for (esv in effectSizeVector){
  for (slv in sigLevelVector){
    T1vsT2MouseList[[index]] <- pwr.t.test(n=120, d=esv, sig.level = slv, type="paired")
    index = index + 1
  }
}

powermatT1vsT2Mouse <- matrix(sapply(T1vsT2MouseList, FUN = function(x) x$power), nrow = length(sigLevelVector), ncol=length(effectSizeVector))
colnames(powermatT1vsT2Mouse) <- effectSizeVector
rownames(powermatT1vsT2Mouse) <- sigLevelVector

T1vsT2Mouse.power80.sig.05 <- pwr.t.test(n=120, sig=0.05, power = 0.80, type="paired")$d
T1vsT2Mouse.power80.sig.10 <- pwr.t.test(n=120, sig=0.10, power = 0.80, type="paired")$d

ANvsHCMouseList <- list()
index <- 1

for (esv in effectSizeVector){
  for (slv in sigLevelVector){
    ANvsHCMouseList[[index]] <- pwr.t2n.test(n1=120, n2=45, d=esv, sig.level = slv)
    index = index + 1
  }
}

powermatANvsHCMouse <- matrix(sapply(ANvsHCMouseList, FUN = function(x) x$power), nrow = length(sigLevelVector), ncol=length(effectSizeVector))
colnames(powermatANvsHCMouse) <- effectSizeVector
rownames(powermatANvsHCMouse) <- sigLevelVector

ANvsHCMouse.power80.sig.05 <- pwr.t2n.test(n1=120, n2=45, sig=0.05, power = 0.80)$d
ANvsHCMouse.power80.sig.10 <- pwr.t2n.test(n1=120, n2=45, sig=0.10, power = 0.80)$d
