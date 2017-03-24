#!/usr/bin/env Rscript

library(moments)
library(labstatR)
library(gstat)
library(sp)
library(rgdal)
library(tools)
library(lubridate)
library(parallel)
library(foreach)
library(doParallel)

cl <- makeCluster(detectCores())
registerDoParallel(cl, cores=detectCores())

dataPath <- "/home/sidereus/documents/projects/kriging/data/"
codePath <- "/home/sidereus/vcs/git/github_personal/kriging/src/"
coordPath <- "/home/sidereus/documents/projects/kriging/data/coords/R_semivariogram.csv"
plotPath <- "/home/sidereus/documents/projects/kriging/plot/"

file.sources <- list.files(paste(codePath,"functions",sep=""),pattern="*.R",full.names=T)
sapply(file.sources,source,.GlobalEnv)

files <- list.files(paste(dataPath), pattern="*.csv", full.names=T)
fileNames <- list.files(paste(dataPath), pattern="*.csv")

nstations <- 10 # for local kriging

for (nFile in 1:length(files)) {

    dataFile <- files[nFile]
    type <- file_path_sans_ext(fileNames[nFile])
    print(type)
    data <- inputDataProcessing(dataFile)
    coordinate <- inputCoordinatesProcessing(coordPath)

    list[totaldf,totalcoord] <- dataCoordProcessing(data, coordinate)
    
    foreach (i=1:length(data[,1]), .inorder=FALSE, .packages=c("gstat","lubridate","tools","rgdal","gstat","sp")) %dopar% {

        date <- data[i,1]
        monthNum <- month(as.POSIXlt(date))
        print(date)
        df_media <- totaldf[[i]]
        coord <- totalcoord[[i]]

        bins <- 8
        computed_cutoff <- computeCutoff()
        computed_width <- computed_cutoff/bins

        omnidirectSperimentalVariogram(date, plotPath, type)
        singleVariogram <- c("Lin","Sph","Exp","Gau")
        krigingType <- c("Ordinary", "KED")
        variogramFitting(singleVariogram, krigingType, monthNum, date, type, dataPath, plotPath, nstations)

    }

}
