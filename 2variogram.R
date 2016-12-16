# --------------------------------------------------------------------------- #
#                      ---- EXERCISE 2 - NEW SECTION ----
# --------------------------------------------------------------------------- #


testToleranceAngle <- function(toleranceAngle, alphaRange) {

    for(index in 1:length(toleranceAngle)){
        angle <- as.numeric(toleranceAngle[index])
                                        # COMPUTATION AND PLOT OF EXPERIMENTAL VARIOGRAM
                                        # EXAMPLE: AVERAGE ANNUAL PRECIPITATION FOR 4 DIRECTIONS
                                        # .................
                                        # variogram command
                                        # no regression ~1
                                        # alpha for a direction variogram
                                        # I'm calculating the directional variogram. If I show that the field is isotropic there is no reason to
                                        # compute in each direction
        prec.dir = variogram(media ~ 1, df_media, width=2000,
                             alpha = alphaRange, tol.hor=angle) #1 vuol dire che non si usano regressori
                                        # df_media: data frame from which $media information are retrieved (both coordinates and values)
                                        # width: the width of subsequent distance intervals into which data point pairs are grouped for semivariance estimates
                                        # alpha: direction in plane (x,y), in positive degrees clockwise from positive y (North):
                                        # alpha=0 for North direction (increasing y), alpha=90 for East direction East (increasing x);
                                        # alpha can be also a vector of directions (see example)
                                        # tol.hor horizontal tolerance angle in degrees -22.5 and +22.5
                             prec.dir
                             pdf(paste("exam/prec_DIRw2000_tol",angle,".pdf",sep=""))
                             print(plot(prec.dir,plot.numbers = TRUE, main=paste("Punti analizzati - tol ",angle,sep="")))
                             print(plot(prec.dir,panel=vgm.panel.xyplot, multipanel = FALSE,main=paste("Semivariogrammi direzionali - tol ",angle,sep=""),xlab="distanza (m)",ylab="semivariogramma",auto.key=TRUE))
                             dev.off()

    }

}

testWidthOmnidirectional <- function(testWidth) {

    for(index in 1:length(testWidth)){
        tmpWidth <- testWidth[index]
        prec.vgm = variogram(media ~ 1, df_media, width=tmpWidth, alpha = 0, tol.hor=90)

        pdf(paste("exam/prec_VarAnalisys_w",tmpWidth,".pdf",sep=""))
            print(plot(prec.vgm, main=paste("Semivariogramma sperimentale - W ",tmpWidth), xlab="distanza (m)", ylab="semivariogramma"))
            dev.off()
    }

}

bestWidth <- function(widthInterval, singleVariogram) {



    for (index in 1:length(singleVariogram)){

        variogramType <- as.character(singleVariogram[index])
        tmpWidth <- widthInterval[1]
        finalWidth <- 0
        minError <- 100000
        while(tmpWidth < widthInterval[2]){

            prec.vgm = variogram(media ~ 1, df_media, width=tmpWidth, alpha = 0, tol.hor=90)
            prec.fit = fit.variogram(prec.vgm, model = vgm(3000,variogramType,10000,10))

            tmpMin <- attributes(prec.fit)$SSErr

            show(tmpWidth)
    
            if (tmpMin < minError) {
                finalWidth <- tmpWidth
                minError <- tmpMin
            }

            tmpWidth = tmpWidth + 100
        }

        pdf(paste("exam/precFit_OMNI_",variogramType,"_w",finalWidth,".pdf",sep=""))
        print(plot(prec.vgm,prec.fit,xlab="Distance",ylab="Semivariance",main=paste(variogramType,"Variogram - w",finalWidth)))
        dev.off()

    }


}

toleranceAngle <- c(22.5,67.5,112.5,157.5)
alphaRange <- c(22.5)

testToleranceAngle(toleranceAngle, alphaRange)

testWidth <- c(2000,4000,5000,6000,7000,8000)

testWidthOmnidirectional(testWidth)

widthInterval <- c(4000,6000)
singleVariogram <- c("Sph","Exp","Gau")

bestWidth(widthInterval, singleVariogram)
