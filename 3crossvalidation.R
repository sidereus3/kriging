# --------------------------------------------------------------------------- #
#                      ---- EXERCISE 3 - NEW SECTION ----
# --------------------------------------------------------------------------- #

# Create Adige DTM (1km) spatial dataframe:
# Loading specified data sets
adige<-read.table("Adige_1000m_wgs84-utm32n.xyz",sep=",")
# Created from 10m resolution DEM and then aggregated at 1 km
# It has been exported to a tabular format: lat, long, elev
colnames(adige)<-c("long","lat","elev")
# We can use gridded command in two ways:
# 1) to return logical (TRUE or FALSE) telling whether the object is gridded or not
# 2) to change a non-gridded (non structured) structure to a gridded one (good plot)
gridded(adige) = ~long+lat # long and lat are the coordinates

# Plot location: ----
#pdf("AdigeDEM_with-stations.pdf")
#plot(adige) # DTM
#points(df_media$long,df_media$lat,col="red",cex=0.5,pch=16) # Measurement stations
# points : is used to add to an existing plot scatter points
#dev.off()


# CROSS VALIDATION PROCEDURE (LEAVE ONE OUT)
# .................
# Leave-one-out cross validation (LOOCV) predicts the value at an observational
# location by leaving out the observed value, using the other points to make the prediction.
# The procedure is then repeated for all the measurement points.
# Semivariogram and interpolation (SK, OK, KED) models should be fitted or predefined.

## do this part with other models (Exp, Gaussian)


computeCrossValidation <- function(variogramType, nmaxStation, krigingType) {

    for(index in 1:length(variogramType)){
        variogramFunction <- variogramType[index]
        prec.vgm = variogram(media~ 1, df_media, width=5400, alpha = 0, tol.hor=90)
        precs.fit = fit.variogram(prec.vgm, model = vgm(3000,variogramFunction,30000,10))

        for(index2 in 1:length(nmaxStation)){
            nmaxVal <- nmaxStation[index2]

            for(index3 in 1:length(krigingType)){
                kriging <- krigingType[index3]
                if(kriging=="Ordinary")
                    prec.cv = krige.cv(media ~ 1, loc = df_media, model = precs.fit, nmax=nmaxVal)
                else if (kriging=="KED")
                    prec.cv = krige.cv(media ~ elev, loc = df_media, model = precs.fit, nmax=nmaxVal)



                res <- prec.cv$residual

                pdf(paste("exam/observed_predict_",variogramFunction,"_nmax",nmaxVal,"_krig",kriging,".pdf",sep=""))
                plot(prec.cv$observed,prec.cv$var1.pred,main=paste("observed_predict_",variogramFunction,"_nmax",nmaxVal,"_krig",kriging,".pdf",sep=""))
                abline(0,1)
                text(220,50,paste("Corr residual",round(cor(prec.cv$var1.pred,prec.cv$residual),digits = 2)))
                text(220,55,paste("Min",round(min(res),digits=2),"- Max",round(max(res),digits=2)))
                text(220,60,paste("Mean",round(mean(res),digits=2)))
                text(220,65,paste("Mean abs",round(mean(abs(res)),digits=2)))
                text(220,70,paste("Var",round(var(res),digits=2)))
                text(220,75,paste("Corr",round(cor(prec.cv$observed,prec.cv$var1.pred),digits=2)))
                dev.off()
                
                pdf(paste("exam/bubble_model",variogramFunction,"_nmax",nmaxVal,"_krig",kriging,".pdf",sep=""))
                print(bubble(prec.cv,"residual",main=paste("Residual - model",variogramFunction,"- nmax",nmaxVal,"- krig",kriging)))
                dev.off()

            }

        }
        
    }

}

variogramType <- c("Sph","Exp","Gau")
nmaxStation <- c(4,8,16,24,32)
krigingType <- c("Ordinary", "KED")

computeCrossValidation(variogramType,nmaxStation,krigingType)

prec.vgm <- variogram(media ~ 1,df_media, width=5900, alpha = 0, tol.hor=90)
precs.fit <- fit.variogram(prec.vgm, model=vgm(10000,"Gau",30000,10))
prec.kriged <- krige(media ~ elev, loc = df_media, adige, model=precs.fit, nmax=32)

pdf("exam/kriging.pdf")
spplot(prec.kriged,"var1.pred",contour=T,col.regions=rainbow(100,start=.5,end=.75),main="Kringing precipitazione - inverno 2004")
dev.off()







df_media <- totaldf[[1]]
prec.vgm <- variogram(media ~ 1,df_media, width=4800, alpha = 0, tol.hor=90)
precs.fit <- fit.variogram(prec.vgm, model=vgm(1000,"Sph",30000,10))
prec.kriged2000 <- krige(media ~ elev, loc = df_media, adige, model=precs.fit, nmax=32)

media[1] <- mean(prec.kriged2000$var1.pred)

df_media <- totaldf[[2]]
prec.vgm <- variogram(media ~ 1,df_media, width=5900, alpha = 0, tol.hor=90)
precs.fit <- fit.variogram(prec.vgm, model=vgm(1000,"Sph",30000,10))
prec.kriged2001 <- krige(media ~ elev, loc = df_media, adige, model=precs.fit, nmax=32)

media[2] <- mean(prec.kriged2001$var1.pred)

df_media <- totaldf[[3]]
prec.vgm <- variogram(media ~ 1,df_media, width=5000, alpha = 0, tol.hor=90)
precs.fit <- fit.variogram(prec.vgm, model=vgm(1000,"Exp",30000,10))
prec.kriged2002 <- krige(media ~ elev, loc = df_media, adige, model=precs.fit, nmax=32)

media[3] <- mean(prec.kriged2002$var1.pred)

df_media <- totaldf[[4]]
prec.vgm <- variogram(media ~ 1,df_media, width=4900, alpha = 0, tol.hor=90)
precs.fit <- fit.variogram(prec.vgm, model=vgm(1000,"Exp",30000,10))
prec.kriged2003 <- krige(media ~ elev, loc = df_media, adige, model=precs.fit, nmax=24)

media[4] <- mean(prec.kriged2003$var1.pred)

df_media <- totaldf[[5]]
prec.vgm <- variogram(media ~ 1,df_media, width=5900, alpha = 0, tol.hor=90)
precs.fit <- fit.variogram(prec.vgm, model=vgm(1000,"Gau",30000,10))
prec.kriged2004 <- krige(media ~ elev, loc = df_media, adige, model=precs.fit, nmax=32)

media[5] <- mean(prec.kriged2004$var1.pred)

df_media <- totaldf[[6]]
prec.vgm <- variogram(media ~ 1,df_media, width=5400, alpha = 0, tol.hor=90)
precs.fit <- fit.variogram(prec.vgm, model=vgm(1000,"Exp",30000,10))
prec.kriged2005 <- krige(media ~ elev, loc = df_media, adige, model=precs.fit, nmax=32)

media[6] <- mean(prec.kriged2005$var1.pred)

df_media <- totaldf[[7]]
prec.vgm <- variogram(media ~ 1,df_media, width=4600, alpha = 0, tol.hor=90)
precs.fit <- fit.variogram(prec.vgm, model=vgm(3000,"Exp",30000,10))
prec.kriged2006 <- krige(media ~ elev, loc = df_media, adige, model=precs.fit, nmax=24)

media[7] <- mean(prec.kriged2006$var1.pred)

df_media <- totaldf[[8]]
prec.vgm <- variogram(media ~ 1,df_media, width=5100, alpha = 0, tol.hor=90)
precs.fit <- fit.variogram(prec.vgm, model=vgm(1000,"Exp",30000,10))
prec.kriged2007 <- krige(media ~ elev, loc = df_media, adige, model=precs.fit, nmax=8)

media[8] <- mean(prec.kriged2007$var1.pred)

df_media <- totaldf[[9]]
prec.vgm <- variogram(media ~ 1,df_media, width=6000, alpha = 0, tol.hor=90)
precs.fit <- fit.variogram(prec.vgm, model=vgm(1000,"Exp",30000,10))
prec.kriged2008 <- krige(media ~ elev, loc = df_media, adige, model=precs.fit, nmax=8)

media[9] <- mean(prec.kriged2008$var1.pred)

df_media <- totaldf[[10]]
prec.vgm <- variogram(media ~ 1,df_media, width=4300, alpha = 0, tol.hor=90)
precs.fit <- fit.variogram(prec.vgm, model=vgm(1000,"Sph",30000,10))
prec.kriged2009 <- krige(media ~ elev, loc = df_media, adige, model=precs.fit, nmax=8)

media[10] <- mean(prec.kriged2009$var1.pred)

year <- c(2000:2009)

pdf("exam/trendPrecipitazione.pdf")
plot(year,media,type="l",main="Trend della precipitazione media invernale dal 2000 al 2009")
points(year,media,col="red")
dev.off()
