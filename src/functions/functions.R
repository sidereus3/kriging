inputDataProcessing <- function(path) {
    inputData<-read.csv(paste(path,sep=""), header=TRUE,sep=",",stringsAsFactors=FALSE)
    colnames(inputData)<-gsub("X","",colnames(inputData))
    inputData$ID <- NULL
    inputData$posix <- as.POSIXct(strptime(inputData[,1],"%Y-%m-%d %H:%M"), origin="1970-01-01", tz='GMT');
    inputData$dataora <- inputData$posix
    inputData$posix <- NULL

    row.names(inputData) <- inputData$dataora

    return(inputData)
}

inputCoordinatesProcessing <- function(path) {
    inputCoordinates<-read.table(paste(path,sep=""),header=TRUE,sep=",",row.names=3,
                           na.string=NA,stringsAsFactors=FALSE)
    inputCoordinates <- inputCoordinates[order(row.names(inputCoordinates)),]
    inputCoordinates$value <- NULL
    return(inputCoordinates)
}

dataCoordProcessing <- function(data, coordinate) {

    df <- rep(list(),12)
    coord <- rep(list(),12)

    for(i in 1:length(data[1,])){
        tmp <- data.frame(array(0, dim=c(length(data[i,])-1,1)))

        for(index in 1:length(data[1,])-1){
            tmp[index,1] <- as.numeric(data[i,index+1])
        }

        row.names(tmp) <- colnames(data[,-1])
        colnames(tmp) <- c("media")

        df_winter <- merge(coordinate,tmp,0)
        coordinates(df_winter)<-~lat+long

        res <- latlong2utm32n(df_winter)

        res_chop <- res
        chop <- 0

        for (row in 1:length(res[,1])) {
            if (is.na(res[row,3]$media)) {
                tmpRow <- row - chop
                res_chop <- res_chop[-c(tmpRow),]
                chop <- chop+1
            }
        }

        df[[i]] <- res_chop
        coord[[i]] <- data.frame(x=coordinates(res_chop)[,1],
                                 y=coordinates(res_chop)[,2])
    }

    return(list(df = df, coord = coord))
}

computeCutoff <- function() {

    x_min <- coord$x[1]
    y_min <- coord$y[1]
    x_max <- coord$x[1]
    y_max <- coord$y[1]

    for (i in 2:length(coord$x)) {

        x_val <- coord$x[i]
        y_val <- coord$y[i]
        
        x_min <- min(x_min, x_val)
        y_min <- min(y_min, y_val)

        x_max <- max(x_max, x_val)
        y_max <- max(y_max, y_val)
        
    }

    delta_x <- x_max - x_min
    delta_y <- y_max - y_min

    diagonal <- sqrt((delta_x)^2 + (delta_y)^2)

    return(diagonal/3)
    
}
 
omnidirectSperimentalVariogram <- function(date, type) {

    prec.vgm = variogram(media ~ 1, df_media, cutoff=computed_cutoff, width=computed_width)
    pdf(paste("/home/sidereus/documents/projects/kriging/plot/prec_VarAnalisys_",type,"_",date,".pdf",sep=""))
    print(plot(prec.vgm, main=paste("Semivariogramma sperimentale - ",date," - ", type), xlab="distanza (m)", ylab="semivariogramma"))
    dev.off()

}

latlong2utm32n <- function(inputDataFrame) {

    proj4string(inputDataFrame) <- CRS("+proj=longlat +datum=WGS84")
    res <- spTransform(inputDataFrame, CRS("+proj=utm +zone=32N ellps=WGS84"))

    return(res)

}

indexToMonth <- function(index) {

    month <- switch(index,
                    '1'="January",
                    '2'="February",
                    '3'="March",
                    '4'="April",
                    '5'="May",
                    '6'="June",
                    '7'="July",
                    '8'="August",
                    '9'="September",
                    '10'="October",
                    '11'="November",
                    '12'="December")

    return(month)
}

variogramFitting <- function(singleVariogram, lineIndex, type, dataPath, plotPath) {

    date <- indexToMonth(lineIndex)

    for (index in 1:length(singleVariogram)){

        variogramType <- as.character(singleVariogram[index])
        variogramData <- read.csv(paste(dataPath,"variogramData/",type,"_",variogramType,".csv",sep=""),sep=";",header=T,row.names=1, stringsAsFactor=F)

        nugget <- variogramData[which(rownames(variogramData) == "TV_nugget:"), lineIndex]
        sill <- variogramData[which(rownames(variogramData) == "TV_sill:"), lineIndex]
        range <- variogramData[which(rownames(variogramData) == "TV_range:"), lineIndex]

        prec.vgm = variogram(media ~ 1, df_media, cutoff=computed_cutoff, width=computed_width)
        prec.fit = fit.variogram(prec.vgm, model = vgm(sill,variogramType,range,nugget),fit.kappa=T)
        pdf(paste(plotPath,"precFit_OMNI_",variogramType,"_",lineIndex,"_",type,".pdf",sep=""))
        print(plot(prec.vgm,prec.fit,xlab="Distance",ylab="Semivariance",main=paste(variogramType,"Variogram - ",date," - ",type), sub=paste("sill ", prec.fit$psill[2]," - range ", prec.fit$range[2])))
        dev.off()

    }

}
