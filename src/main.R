library(moments);library(labstatR);library(gstat);library(sp);library(rgdal)

dataPath <- "data/medie.csv"
coordPath <- "data/R_semivariogram.csv"

source(functions.R)

data <- inputDataProcessing(dataPath)
coordinates <- inputCoordinatesProcessing(coordPath)

totaldf <- list()
totalcoord <- list()

for(i in 1:length(data[1,])){
    tmp <- data.frame(array(0, dim=c(length(data[i,])-1,1)))

    for(index in 1:length(data[1,])-1){
        tmp[index,1] <- as.numeric(data[i,index+1])
    }

    row.names(tmp) <- colnames(data[,-1])
    colnames(tmp) <- c("media")

    df_winter <- merge(coordinate,tmp,0)
    coordinates(df_winter)<-~lat+long

    proj4string(df_winter) <- CRS("+proj=longlat +datum=WGS84")
    res <- spTransform(df_winter, CRS("+proj=utm +zone=32N ellps=WGS84"))

    res_chop <- res
    chop <- 0
    
    for (row in 1:length(res[,1])) {
        if (is.na(res[row,3]$media)) {
            tmpRow <- row - chop
            res_chop <- res_chop[-c(tmpRow),]
            chop <- chop+1
            }
        }
    
    totaldf[[i]] <- res_chop
    totalcoord[[i]] <- data.frame(x=coordinates(res_chop)[,1], y=coordinates(res_chop)[,2])
}

for (i in 1:length(totaldf)) {

    df_media <- totaldf[[i]]
    coord <- totalcoord[[i]]
    type <- "medie"

    bins <- 8
    computed_cutoff <- computeCutoff()
    computed_width <- computed_cutoff/bins

    omnidirectSperimentalVariogram(i, type)
    singleVariogram <- c("Sph","Exp","Gau","Lin")
    variogramFitting(singleVariogram, i, type)

}
