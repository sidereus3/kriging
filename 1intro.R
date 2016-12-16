is.leapyear <- function(year){
    #http://en.wikipedia.org/wiki/Leap_year
    return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))
}

library(moments);library(labstatR);library(gstat);library(sp)

# --------------------------------------------------------------------------- #
#                            ---- INTRODUCTION ----
# --------------------------------------------------------------------------- #
#### DATAFRAME preparation ----
### Data and coordinates:
datipioggia<-read.table("PRC-N-TagSum_ANN-MetG_TN-BZ_1970-01-01_2010-01-01.csv",
                        header=TRUE,sep=",",na.string=NA,stringsAsFactors=FALSE)
# stringAsFactors=FALSE : character vectors are not converted to factors
colnames(datipioggia)<-gsub("X","",colnames(datipioggia))
# gsub: R adds an "X" before numeric names, this command substitutes "X" with ""

coordinate<-read.table("Adige_CoordinatePRC.csv",header=TRUE,sep=",",row.names=1,
                       na.string=NA,stringsAsFactors=FALSE)
# header=TRUE : the first row of the table is readed as the header (column names)
# row.names=1 : the first column of the table contain the row names

# Preparing the table: we are splitting (strsplit command) the column 1 (dates)
# in three columns (YYYY, MM, DD) using "-" as separator.
# We are column-binding (cbind) them to the original dataframe.
dp_ts<-cbind(sapply(strsplit(datipioggia[,1],"-"),"[",1),
             sapply(strsplit(datipioggia[,1],"-"),"[",2),
             sapply(strsplit(datipioggia[,1],"-"),"[",3),
             datipioggia)
colnames(dp_ts)[(1:3)]<-c("YYYY","MM","DD")
#View(dp_ts)

#### Selecting data ----
# Selecting a PERIOD between two dates (subset)
dp_ts<-subset(dp_ts,c(dp_ts$Date>="2000-01-01"&dp_ts$Date<="2010-01-01"))
# Seasonal: running the selection over the column "MM", choosing those months
# which are identified by 01, 02, 03 (January, February, March - winter)


dp_ts<-dp_ts[dp_ts$MM %in% c("01","02","03"),]

# Cleaning NAs:
# Selecting those stations in which the number of NAs is less than 10%
NApres<-t(as.data.frame(colSums(is.na(dp_ts))/dim(dp_ts)[1]*100))
# colSums running on is.na calculates the number of NAs for each column of dp_ts
dp_ts<-rbind(NApres,dp_ts)
# Add the vector NApres to the dataframe dp_ts
# Warning! If column number differs, the command doesn't work!
dp_ts<-dp_ts[-1,which(dp_ts[1,]<10)]
# which selects only columns in which NAs are less than 10%
# "-1" exclude the first row (containing NAs %)

# Substituting remaining NAs with "0":
dp_ts[is.na(dp_ts)==TRUE]<-0

# Are NAs still there? (Number of NAs has to be 0!)
length(dp_ts[which(is.na(dp_ts)==TRUE,arr.ind=TRUE)])

# Aggregate by Year ---
annuale<-aggregate(dp_ts[,-(1:4)], list(dp_ts$YYYY),sum,na.rm=FALSE);
# dp_ts[,-(1:4)] : excludes columns which contain date information (YYYY, MM, DD, Date)
# aggregate: is used on dp_ts$YYYY to aggregate (sum) together rows with the same year


mensile <- aggregate(dp_ts[,-(1:4)], list(dp_ts$YYYY,dp_ts$MM),sum,na.rm=FALSE)

colnames(mensile)[(1:2)] <- c("YYYY","MM")

seasonal <- aggregate(mensile[,-(1:2)],list(mensile$YYYY), sum,na.rm=F)

rownames(seasonal)<-seasonal[,1]
# Changing row names with labels present in annuale dataframe

tmpSeasonal = seasonal

## for(index in 1:length(seasonal[,1])){

##     show(seasonal[index,1])
##     if(is.leapyear(as.numeric(seasonal[index,1]))) tmpSeasonal[index,2:length(seasonal[1,])] <- tmpSeasonal[index,2:length(seasonal[1,])]/91
##     else tmpSeasonal[index,2:length(seasonal[1,])] <- tmpSeasonal[index,2:length(seasonal[1,])]/90
    
## }

totaldf <- list()

for(i in 1:length(seasonal[,1])){

    tmp <- data.frame(array(0, dim=c(length(seasonal[i,])-1,1)))

    for(index in 1:length(seasonal[1,])-1){
        tmp[index,1] <- as.numeric(seasonal[i,index+1])
    }

    row.names(tmp) <- colnames(seasonal[,-1])
    colnames(tmp) <- c("media")

    df_winter <- merge(coordinate,tmp,0)
    coordinates(df_winter)<-~long+lat

    totaldf[[i]] <- df_winter
    
}

df_media <- totaldf[[1]]


pdf("exam/bubbles.pdf")
bubble(df_media, "media", col=c("blue"),pch=20,
       key.space="right", main="Precipitazione cumulata - inverno", maxsize=2,key.entries=30*(1:10))
dev.off()

pdf("exam/hist.pdf")
hist(df_media$media, main="Precipitazione cumulata - inverno", breaks=seq(0,300,10),col="gray", freq=F, xlab="Classes")
dev.off()

pdf("exam/ecdf.pdf")
plot(ecdf(df_media$media), do.points = F, verticals = T,ylab="ECDF", xlab="Precipitazione media [mm]")
x1 <- seq(0,300,1)
lines(x1,pnorm(x1,mean=mean(df_media$media),sd=sqrt(var(df_media$media))), lty=3)
dev.off()


## # Annual mean ---
## medialp<-as.data.frame(colMeans(annuale[,-1])) # We do not consider first column
## #mediaseason <- as.data.frame(colMeans(seasonal[,-1]))
## colnames(medialp)<-c("media");summary(medialp);

## #### Create S4 object type (SpatialPointsDataFrame) ----
## # - row names of coordinate table and dataframe should coincide!
## # - dataframe class should be numeric!
## ## S4 for every year - if you want to mantain yearly subdivision the following
## # structure can be used: - optional
## #df_ann<-t(annuale);class(df_ann)<-c("numeric");df_ann<-merge(coordinate,df_ann,0)
## #coordinates(df_ann)<-~long+lat

## ## S4 for average annual precipitation:
## df_media<-merge(coordinate,medialp,0)
## # merge : merge coordinate and medialp dataframes using common row names (parameter "0")
## # which contains the station code (e.g. T0001)
## # it is necessary because the dataframes differ: coordinate contains every stations,
## # "medialp" contains only stations with NAs < 10%

## coordinates(df_media)<-~long+lat
## # transforms lat long information to create a SpatialPointsDataFrame
## #### END! - DATAFRAME preparation ----
