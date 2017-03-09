library(moments)
library(labstatR)
library(gstat)
library(sp)
library(rgdal)
library(tools)

dataPath <- "/home/sidereus/documents/projects/kriging/data/"
codePath <- "/home/sidereus/vcs/git/github_personal/kriging/src/"
coordPath <- "/home/sidereus/documents/projects/kriging/data/coords/R_semivariogram.csv"
plotPath <- "/home/sidereus/documents/projects/kriging/plot/"

file.sources <- list.files(paste(codePath,"functions/",sep=""),pattern="*.R",full.names=T)
sapply(file.sources,source,.GlobalEnv)

files <- list.files(paste(dataPath), pattern="*.csv", full.names=T)
fileNames <- list.files(paste(dataPath), pattern="*.csv")

for (nFile in 1:length(files)) {

    dataFile <- files[nFile]
    type <- file_path_sans_ext(fileNames[nFile])
    print(type)
    data <- inputDataProcessing(dataFile)
    coordinate <- inputCoordinatesProcessing(coordPath)

    list[totaldf,totalcoord] <- dataCoordProcessing(data, coordinate)
    
    for (i in 1:12) {

        print(indexToMonth(i))
        df_media <- totaldf[[i]]
        coord <- totalcoord[[i]]

        bins <- 8
        computed_cutoff <- computeCutoff()
        computed_width <- computed_cutoff/bins

        omnidirectSperimentalVariogram(i, type)
        singleVariogram <- c("Lin","Sph","Exp","Gau")
        variogramFitting(singleVariogram, i, type, dataPath, plotPath)

    }

}
