inputDataProcessing <- function(path) {
    inputData<-read.csv(paste(inputDataPath,sep=""), header=TRUE,sep=",",stringsAsFactors=FALSE)
    colnames(inputData)<-gsub("X","",colnames(inputData))
    inputData$ID <- NULL
    inputData$posix <- as.POSIXct(strptime(inputData[,1],"%Y-%m-%d %H:%M"), origin="1970-01-01", tz='GMT');
    inputData$dataora <- inputData$posix
    inputData$posix <- NULL

    row.names(inputData) <- inputData$dataora

    return(inputData)
}

inputCoordinatesProcessing <- function(path) {
    inputCoordinates<-read.table(paste(coordPath,sep=""),header=TRUE,sep=",",row.names=3,
                           na.string=NA,stringsAsFactors=FALSE)
    inputCoordinates <- inputCoordinates[order(row.names(inputCoordinates)),]
    inputCoordinates$value <- NULL
    return(inputCoordinates)
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

    pdf(paste("plot/prec_VarAnalisys_",type,"_",date,".pdf",sep=""))
    print(plot(prec.vgm, main=paste("Semivariogramma sperimentale - ",date," - ", type), xlab="distanza (m)", ylab="semivariogramma"))
    dev.off()

    return()

}

variogramFitting <- function(singleVariogram, date, type) {



    for (index in 1:length(singleVariogram)){

        variogramType <- as.character(singleVariogram[index])
        prec.vgm = variogram(media ~ 1, df_media, cutoff=computed_cutoff, width=computed_width)
        prec.fit = fit.variogram(prec.vgm, model = vgm(3000,variogramType,10000,10))
        pdf(paste("plot/precFit_OMNI_",variogramType,"_",date,"_",type,".pdf",sep=""))
        print(plot(prec.vgm,prec.fit,xlab="Distance",ylab="Semivariance",main=paste(variogramType,"Variogram - ",date," - ",type)))
        dev.off()

    }

}
