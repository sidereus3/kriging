#!/usr/bin/env Rscript

library(moments)
library(labstatR)
library(gstat)
library(sp)
library(rgdal)
library(tools)
library(raster)
library(lubridate)
library(parallel)
library(foreach)
library(doParallel)

cl <- makeCluster(detectCores())
registerDoParallel(cl, cores=detectCores())

workdir <- "/home/sidereus/documents/projects/kriging/codeR/"
codePath <- "/home/sidereus/vcs/git/github_personal/kriging/codeR/src/"
dataPath <- paste(workdir,"data/",sep="")
coordPath <- paste(workdir,"data/coords/",sep="")
plotPath <- paste(workdir,"plot/",sep="")
outDataPath <- paste(workdir,"outData/",sep="")

# parses and applies functions.R and list.R
file.sources <- list.files(paste(codePath,"functions",sep=""),pattern="*.R",full.names=T)
sapply(file.sources,source,.GlobalEnv)

# reads csv files and file names
files <- list.files(paste(dataPath), pattern="*.csv", full.names=T)
fileNames <- list.files(paste(dataPath), pattern="*.csv")

nstations <- 10 # for local kriging

# reads DEM raster from xyz file table
dem <- raster(paste(coordPath,"DEM.xyz",sep=""),sep=" ")
# reads the DEM as data frame
demToGrid <- read.table(paste(coordPath,"DEM.xyz",sep=""),sep=" ")

# creates a spatial pixels data frame with proper projection
colnames(demToGrid) <- c("x","y","quota")
gridded(demToGrid) <- ~x+y
proj4string(demToGrid) <- CRS("+init=epsg:32632 +proj=utm +zone=32N ellps=WGS84")

# this function returns a dataframe containing index, elevation and name of each station
stationsIndex <- function() {

    coordinate$name <- row.names(coordinate)
    coordinates(coordinate) <- ~long+lat
    tmpCoord <- latlong2utm32n(coordinate)
    # create the dataframe to host data
    tmpCoord.df <- data.frame(index=coordinates(tmpCoord)[,1],
                              quota=tmpCoord$quota,
                              name=tmpCoord$name)

    for (i in 1:length(tmpCoord.df[,1])) {
        point = matrix(nrow=1,ncol=2)
        point[1,1] = coordinates(tmpCoord)[i,1]
        point[1,2] = coordinates(tmpCoord)[i,2]

        tmpCoord.df[i,1] = cellFromXY(dem, point)
    }

    return(tmpCoord.df)

}

# loops over files, one may have min, max, avg
for (nFile in 1:length(files)) {

    dataFile <- files[nFile]
    # returns file name without extension
    # for plotting name purposes only
    type <- file_path_sans_ext(fileNames[nFile])
    print(type)

    # data preprocessing. This is case specific
    data <- inputDataProcessing(dataFile)

    # sets up the data frame for every semivariogram
    outdata_lin <- matrix(nrow=length(data[,1]),ncol=length(data[1,]))
    outdata_lin <- data.frame(outdata_lin)
    colnames(outdata_lin) <- colnames(data)

    outdata_exp <- matrix(nrow=length(data[,1]),ncol=length(data[1,]))
    outdata_exp <- data.frame(outdata_exp)
    colnames(outdata_exp) <- colnames(data)

    outdata_sph <- matrix(nrow=length(data[,1]),ncol=length(data[1,]))
    outdata_sph <- data.frame(outdata_sph)
    colnames(outdata_sph) <- colnames(data)

    outdata_gau <- matrix(nrow=length(data[,1]),ncol=length(data[1,]))
    outdata_gau <- data.frame(outdata_gau)
    colnames(outdata_gau) <- colnames(data)

    outdata_bes <- matrix(nrow=length(data[,1]),ncol=length(data[1,]))
    outdata_bes <- data.frame(outdata_bes)
    colnames(outdata_bes) <- colnames(data)

    # then makes a copy for ked
    outdata_ked_lin <- outdata_lin
    outdata_ked_exp <- outdata_exp
    outdata_ked_sph <- outdata_sph
    outdata_ked_gau <- outdata_gau
    outdata_ked_bes <- outdata_bes

    # meteo stations coordinates and elevation
    coordinate <- inputCoordinatesProcessing(coordPath)

    # data per date plus just valid values
    list[totaldf,totalcoord] <- dataCoordProcessing(data, coordinate)

    # gets a dataframe with grid position of each meteo station
    saveVal <- stationsIndex()

    #foreach (i=1:length(data[,1]), .inorder=FALSE, .packages=c("gstat","lubridate","tools","rgdal","raster")) %dopar% {
    for (i in 1:length(data[,1])) {
        timestep <- i
        date <- data[i,1]
        monthNum <- month(as.POSIXlt(date))
        print(date)

        df_media <- totaldf[[i]]
        coord <- totalcoord[[i]]

        bins <- 8
        computed_cutoff <- computeCutoff()
        computed_width <- computed_cutoff/bins

        # computation of sperimental variogram + plot
        omnidirectSperimentalVariogram(date, plotPath, type)
        singleVariogram <- c("Lin","Sph","Exp","Gau")
        krigingType <- c("Ordinary", "KED")
        outdata <- variogramFitting(singleVariogram, krigingType,
                                    monthNum, date, type, dataPath,
                                    plotPath, nstations, coord,
                                    outdata_lin, outdata_exp, outdata_sph, outdata_gau,
                                    outdata_ked_lin, outdata_ked_exp, outdata_ked_sph, outdata_ked_gau, outdata_bes, outdata_ked_bes)
        outdata_lin <- outdata[[1]]
        outdata_exp <- outdata[[2]]
        outdata_sph <- outdata[[3]]
        outdata_gau <- outdata[[4]]
        outdata_ked_lin <- outdata[[5]]
        outdata_ked_exp <- outdata[[6]]
        outdata_ked_sph <- outdata[[7]]
        outdata_ked_gau <- outdata[[8]]
        outdata_bes <- outdata[[9]]
        outdata_ked_bes <- outdata[[10]]

    }

    write.csv(outdata_lin, paste(outDataPath,"test",type,"lin.csv",sep=""), sep=",",row.names=FALSE)
    write.csv(outdata_gau, paste(outDataPath,"test",type,"gau.csv",sep=""), sep=",",row.names=FALSE)
    write.csv(outdata_sph, paste(outDataPath,"test",type,"sph.csv",sep=""), sep=",",row.names=FALSE)
    write.csv(outdata_exp, paste(outDataPath,"test",type,"exp.csv",sep=""), sep=",",row.names=FALSE)
    write.csv(outdata_bes, paste(outDataPath,"test",type,"bes.csv",sep=""), sep=",",row.names=FALSE)

    write.csv(outdata_ked_lin, paste(outDataPath,"test",type,"ked_lin.csv",sep=""), sep=",",row.names=FALSE)
    write.csv(outdata_ked_gau, paste(outDataPath,"test",type,"ked_gau.csv",sep=""), sep=",",row.names=FALSE)
    write.csv(outdata_ked_sph, paste(outDataPath,"test",type,"ked_sph.csv",sep=""), sep=",",row.names=FALSE)
    write.csv(outdata_ked_exp, paste(outDataPath,"test",type,"ked_exp.csv",sep=""), sep=",",row.names=FALSE)
    write.csv(outdata_ked_bes, paste(outDataPath,"test",type,"bes_exp.csv",sep=""), sep=",",row.names=FALSE)

}
