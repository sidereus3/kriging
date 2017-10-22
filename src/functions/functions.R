inputDataProcessing <- function(path) {
    inputData<-read.csv(paste(path,sep=""), header=TRUE,sep=",",stringsAsFactors=FALSE)
    colnames(inputData)<-gsub("X","",colnames(inputData))
    inputData$ID <- NULL
    inputData$posix <- as.POSIXct(strptime(inputData[,1],"%Y-%m-%d %H:%M"), origin="1970-01-01", tz='GMT');
    inputData$dataora <- inputData$posix
    inputData$posix <- NULL

    for (i in 1:length(colnames(inputData))) {
        if (colnames(inputData)[i] == 90290)
            index <- i
    }

    row.names(inputData) <- inputData$dataora
    inputData <- inputData[,-index]

    return(inputData)
}

inputCoordinatesProcessing <- function(path) {
    inputCoordinates<-read.table(paste(path,"R_semivariogram.csv",sep=""),header=TRUE,sep=",",row.names=3,
                           na.string=NA,stringsAsFactors=FALSE)
    inputCoordinates <- inputCoordinates[order(row.names(inputCoordinates)),]
    inputCoordinates$value <- NULL
    return(inputCoordinates)
}

dataCoordProcessing <- function(data, coordinate) {

    df <- rep(list(),length(data[1,]))
    coord <- rep(list(),length(data[1,]))

    for(i in 1:length(data[1,])){
        tmp <- data.frame(array(0, dim=c(length(data[i,])-1,1)))

        for(index in 1:length(data[1,])-1){
            tmp[index,1] <- as.numeric(data[i,index+1])
        }

        row.names(tmp) <- colnames(data[,-1])
        colnames(tmp) <- c("media")

        df_winter <- merge(coordinate,tmp,0)
        coordinates(df_winter)<-~long+lat

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
                                 y=coordinates(res_chop)[,2],
                                 z=res_chop$quota)
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
 
omnidirectSperimentalVariogram <- function(date, plotPath, type) {

    prec.vgm = variogram(media ~ 1, df_media, cutoff=computed_cutoff, width=computed_width)
    pdf(paste(plotPath,"prec_VarAnalisys_",type,"_",date,".pdf",sep=""))
    print(plot(prec.vgm, main=paste("Semivariogramma sperimentale - ",date," - ", type), xlab="distanza (m)", ylab="semivariogramma"))
    dev.off()

}

latlong2utm32n <- function(inputDataFrame) {

    proj4string(inputDataFrame) <- CRS("+init=epsg:4326")
    res <- spTransform(inputDataFrame, CRS("+init=epsg:32632 +proj=utm +zone=32N ellps=WGS84"))

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

computeCrossValidation <- function(krigingType, variogramType, type, prec.fit, date, plotPath, nstations, coord,
                                   data_lin, data_exp, data_sph, data_gau,
                                   data_ked_lin, data_ked_exp, data_ked_sph, data_ked_gau, data_bes, data_ked_bes) {

    for(i in 1:length(krigingType)) {
        kriging <- krigingType[i]
        if (kriging=="Ordinary") {
            prec.cv <- krige.cv(media ~ 1, loc=df_media, model=prec.fit)
            prec.kriged <- krige(media ~ 1, loc=df_media, demToGrid,model=prec.fit)

            tmp <- raster(prec.kriged)
            if (variogramType=="Lin") {
                data_lin[timestep,1] <- as.character(date)
                for (j in 1:length(saveVal$index)) {
                    index <- saveVal$index[j]
                    data_lin[timestep,j+1] <- as.numeric(tmp[index])
                }
            } else if (variogramType=="Gau") {
                data_gau[timestep,1] <- as.character(date)
                for (j in 1:length(saveVal$index)) {
                    index <- saveVal$index[j]
                    data_gau[timestep,j+1] <- as.numeric(tmp[index])
                }
            } else if (variogramType=="Sph") {
                data_sph[timestep,1] <- as.character(date)
                for (j in 1:length(saveVal$index)) {
                    index <- saveVal$index[j]
                    data_sph[timestep,j+1] <- as.numeric(tmp[index])
                }
            } else if (variogramType=="Exp") {
                data_exp[timestep,1] <- as.character(date)
                for (j in 1:length(saveVal$index)) {
                    index <- saveVal$index[j]
                    data_exp[timestep,j+1] <- as.numeric(tmp[index])
                }
            } else if (variogramType=="Bes") {
                data_bes[timestep,1] <- as.character(date)
                for (j in 1:length(saveVal$index)) {
                    index <- saveVal$index[j]
                    data_bes[timestep,j+1] <- as.numeric(tmp[index])
                }
            }

            if (nstations > 0) {
                ## prec.kriged_loc <- krige(media ~ 1, loc=df_media,__missing_xyz__,model=prec.fit, nmax=10)
            }
        } else if (kriging=="KED") {
            prec.cv <- krige.cv(media ~ quota, loc=df_media, model=prec.fit)
            prec.kriged <- krige(media ~ quota, loc=df_media, demToGrid,model=prec.fit)

            tmp <- raster(prec.kriged)
            if (variogramType=="Lin") {
                data_ked_lin[timestep,1] <- as.character(date)
                for (j in 1:length(saveVal$index)) {
                    index <- saveVal$index[j]
                    data_ked_lin[timestep,j+1] <- as.numeric(tmp[index])
                }
            } else if (variogramType=="Gau") {
                data_ked_gau[timestep,1] <- as.character(date)
                for (j in 1:length(saveVal$index)) {
                    index <- saveVal$index[j]
                    data_ked_gau[timestep,j+1] <- as.numeric(tmp[index])
                }
            } else if (variogramType=="Sph") {
                data_ked_sph[timestep,1] <- as.character(date)
                for (j in 1:length(saveVal$index)) {
                    index <- saveVal$index[j]
                    data_ked_sph[timestep,j+1] <- as.numeric(tmp[index])
                }
            } else if (variogramType=="Exp") {
                data_ked_exp[timestep,1] <- as.character(date)
                for (j in 1:length(saveVal$index)) {
                    index <- saveVal$index[j]
                    data_ked_exp[timestep,j+1] <- as.numeric(tmp[index])
                }
            } else if (variogramType=="Bes") {
                data_ked_bes[timestep,1] <- as.character(date)
                for (j in 1:length(saveVal$index)) {
                    index <- saveVal$index[j]
                    data_ked_bes[timestep,j+1] <- as.numeric(tmp[index])
                }
            }

            if (nstations > 0) {
                ## prec.kriged_loc <- krige(media ~ elev, loc=df_media,__missing_xyz__,model=prec.fit, nmax=10)
            }
        }
        
        pdf(paste(plotPath,"observed_predict_",type,"_",date,"_",variogramType,"_krig_", kriging,".pdf",sep=""))
        plot(prec.cv$observed,prec.cv$var1.pred,main=paste("observed_predict_",type,"_",date,"_",variogramType,"_krig_",kriging,sep=""))
        abline(0,1)
        dev.off()

        pdf(paste(plotPath,"krig_",kriging,"_variogram_",variogramType,"_date_",date,".pdf",sep=""))
        spplot(prec.kriged,"var1.pred",contour=T,col.regions=rainbow(100,start=.5,end=.75),main=paste("Kringing temperatura -",date))
        dev.off()
        if (nstations > 0) {
            ## pdf(paste(plotPath,"krig_",kriging,"_variogram_",variogramType,"_localnstat",nstations,".pdf",sep=""))
            ## ssplot(prec.kriged,"var1.pred",contour=T,col.regions=rainbow(100,start=.5,end=.75),main="Kringing precipitazione - inverno 2004")
            ## dev.off()
        }
    }

    dataout <- rep(list(),10)

    dataout[[1]] <- data_lin
    dataout[[2]] <- data_exp
    dataout[[3]] <- data_sph
    dataout[[4]] <- data_gau
    dataout[[5]] <- data_ked_lin
    dataout[[6]] <- data_ked_exp
    dataout[[7]] <- data_ked_sph
    dataout[[8]] <- data_ked_gau
    dataout[[9]] <- data_bes
    dataout[[10]] <- data_ked_bes

    return(dataout)

}

variogramFitting <- function(singleVariogram, krigingType, monthNum, date, type, dataPath, plotPath, nstations, coord,
                             data_lin, data_exp, data_sph, data_gau,
                             data_ked_lin, data_ked_exp, data_ked_sph, data_ked_gau, data_bes, data_ked_bes) {

    monthName <- indexToMonth(monthNum)

    for (index in 1:length(singleVariogram)){

        variogramType <- as.character(singleVariogram[index])
        variogramData <- read.csv(paste(dataPath,"variogramData/",
                                        type,"_",
                                        variogramType,".csv",sep=""),
                                  sep=";",header=T,row.names=1, stringsAsFactor=F)

        nugget <- variogramData[which(rownames(variogramData) == "TV_nugget:"), monthNum]
        sill <- variogramData[which(rownames(variogramData) == "TV_sill:"), monthNum]
        range <- variogramData[which(rownames(variogramData) == "TV_range:"), monthNum]

        prec.vgm = variogram(media ~ 1, df_media,
                             cutoff=computed_cutoff,
                             width=computed_width)
        prec.fit = fit.variogram(prec.vgm,
                                 model = vgm(sill,variogramType,range,nugget),
                                 fit.kappa=T)
        pdf(paste(plotPath,"precFit_OMNI_",
                  variogramType,"_",
                  date,"_",type,".pdf",sep=""))
        print(plot(prec.vgm,prec.fit,xlab="Distance",ylab="Semivariance",
                   main=paste(variogramType,"Variogram - ",date," - ",type),
                   sub=paste("sill ",prec.fit$psill[2]," - range ",prec.fit$range[2])))
        dev.off()

        outdata <- computeCrossValidation(krigingType, variogramType, type, prec.fit, date, plotPath, nstations, coord,
                                          data_lin, data_exp, data_sph, data_gau,
                                          data_ked_lin, data_ked_exp, data_ked_sph, data_ked_gau, data_bes, data_ked_bes)
        data_lin <- outdata[[1]]
        data_exp <- outdata[[2]]
        data_sph <- outdata[[3]]
        data_gau <- outdata[[4]]
        data_ked_lin <- outdata[[5]]
        data_ked_exp <- outdata[[6]]
        data_ked_sph <- outdata[[7]]
        data_ked_gau <- outdata[[8]]
        data_bes <- outdata[[9]]
        data_ked_bes <- outdata[[10]]
    }

    dataout <- rep(list(),10)

    dataout[[1]] <- data_lin
    dataout[[2]] <- data_exp
    dataout[[3]] <- data_sph
    dataout[[4]] <- data_gau
    dataout[[5]] <- data_ked_lin
    dataout[[6]] <- data_ked_exp
    dataout[[7]] <- data_ked_sph
    dataout[[8]] <- data_ked_gau
    dataout[[9]] <- data_bes
    dataout[[10]] <- data_ked_bes

    return(dataout)

}
