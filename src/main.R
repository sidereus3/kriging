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

workdir <- "/home/sidereus/documents/projects/kriging/"
codePath <- "/home/sidereus/vcs/git/github_personal/kriging/src/"
dataPath <- paste(workdir,"data/",sep="")
coordPath <- paste(workdir,"data/coords/",sep="")
plotPath <- paste(workdir,"plot/",sep="")
outDataPath <- paste(workdir,"outData/",sep="")

file.sources <- list.files(paste(codePath,"functions",sep=""),pattern="*.R",full.names=T)
sapply(file.sources,source,.GlobalEnv)

files <- list.files(paste(dataPath), pattern="*.csv", full.names=T)
fileNames <- list.files(paste(dataPath), pattern="*.csv")

nstations <- 10 # for local kriging

demToGrid <- read.table(paste(coordPath,"DEM.xyz",sep=""),sep=" ")
colnames(demToGrid) <- c("x","y","quota")
gridded(demToGrid) <- ~x+y
proj4string(demToGrid) <- CRS("+init=epsg:32632 +proj=utm +zone=32N ellps=WGS84")

for (nFile in 1:length(files)) {

    dataFile <- files[nFile]
    type <- file_path_sans_ext(fileNames[nFile])
    print(type)
    data <- inputDataProcessing(dataFile)
    coordinate <- inputCoordinatesProcessing(coordPath)

    list[totaldf,totalcoord] <- dataCoordProcessing(data, coordinate)

    oldMonthNum <- 0

    foreach (i=1:length(data[,1]), .inorder=FALSE, .packages=c("gstat","lubridate","tools","rgdal")) %dopar% {
    #for (i in 1:length(data[,1])) {
        date <- data[i,1]
        monthNum <- month(as.POSIXlt(date))
        print(date)

        ## if (monthNum > oldMonthNum) {
        ##     head <- colnames(data)
        ##     write(head, file=paste(outDataPath,"monthNum",sep=""),sep=",")
        ##     oldMonthNum <- oldMonthNum + 1
        ## } else {
        ##     write(val, file=paste(outDataPath,monthNum,".csv"))
        ## }

        df_media <- totaldf[[i]]
        coord <- totalcoord[[i]]

        bins <- 8
        computed_cutoff <- computeCutoff()
        computed_width <- computed_cutoff/bins

        omnidirectSperimentalVariogram(date, plotPath, type)
        singleVariogram <- c("Lin","Sph","Exp","Gau")
        krigingType <- c("Ordinary", "KED")
        variogramFitting(singleVariogram, krigingType,
                         monthNum, date, type, dataPath,
                         plotPath, nstations, coord)

    }

}
