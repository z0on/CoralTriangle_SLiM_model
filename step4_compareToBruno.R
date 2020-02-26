# repace directory location below with the name of the cloned github repo directory
setwd("~/Dropbox/Documents/coral_triangle_SLiM/trianglePaper_finalSubmission/CoralTriangle_SLiM_model/")
library(tidyverse)
library(viridis)
library(rworldmap)
map.polygon <- getMap(resolution = "low")

# computing relative reef responses (cover change relative to pre-warming value), for the main parameter setting

load('~/Dropbox/Documents/coral_triangle_SLiM/coral_triangle_model/version3-4/REv32_sin_8_8_base.RData')
c20=subset(covall,dec==20)$cover
c10=subset(covall,dec==10)$cover
c5=subset(covall,dec==5)$cover
c0=subset(covall,dec==-1)$cover
rr=data.frame(cbind(r50=(c5-c0)/c0, r100=(c10-c0)/c0,r200=(c20-c0)/c0))
load("sst_modified.RData")
sst=cbind(sst,rrv3=rr)
names(sst)

# reading bruno data, subsetting them for Indo-Pacific
bruno=read.csv("Global_Cover_1_5_10_bruno_may28_2019.csv")
bruno=subset(bruno, LONGITUDE>min(sst$lon)-1 & LONGITUDE<max(sst$lon)+1 & LATITUDE>min(sst$lat)-1 & LATITUDE<max(sst$lat)+1)

# assembling table of coordinats for two datasets, SLiM model ("tri", from "coral triangle") and "bruno"
lat=c(sst$lat,bruno$LATITUDE)
lon=c(sst$lon,bruno$LONGITUDE)
treefs=paste("t",sst$reef_id,sep="")
breefs=paste("b",bruno$REEF_ID_GR,sep="")
coords=data.frame(cbind(lat,lon,"reef"=c(treefs,breefs)))
c2=list();i=1
for (r in unique(coords$reef)) {
	c2[[i]]=subset(coords,reef==r)[1,]
	i=i+1
}
c2=do.call(rbind,c2)
row.names(c2)=c2$reef
br=grep("b",c2$reef)
tr=grep("t",c2$reef)
c2$set="bruno"
c2$set[tr]="tri"
head(c2)

#------------ assembling table with predictors and coral cover changes for clusters including BOTH modeled and Bruno 
# NB: lph in the new table = logprophots in sst dataframe = pr05 in the paper

# clustering reefs witin 1 degree
hc=hclust(dist(c2[,c("lat","lon")]),method="complete")
c2$grp=cutree(hc,h=1)

c2$lat=as.numeric(as.character(c2$lat))
c2$lon=as.numeric(as.character(c2$lon))
cc=1
dec=list(); cc=1
for (g in unique(c2$grp)) {
	message(g)
	ss=subset(c2,grp==g)
	if("bruno" %in% ss$set &  "tri" %in% ss$set) {
			lon=mean(ss$lon)
			lat=mean(ss$lat)
			sb=sub("b","",as.character(subset(ss,set=="bruno")$reef))
			brd=subset(bruno,REEF_ID_GR %in% sb)
			if (range(brd$YEAR)[2]-range(brd$YEAR)[1]>10 & range(brd$YEAR)[2]>2000) {
				ye=summary(lm(HARD_COR_P~YEAR,brd))$coefficients[2,1] # beta for year, bruno
				pe=summary(lm(HARD_COR_P~YEAR,brd))$coefficients[2,4] # p-value for year, bruno
				st=sub("t","",as.character(subset(ss,set=="tri")$reef))
				d850=mean(subset(sst,reef_id %in% st)$rrv3.r50)
				d8100=mean(subset(sst,reef_id %in% st)$rrv3.r100)
				d8200=mean(subset(sst,reef_id %in% st)$rrv3.r200)
				lph=mean(subset(sst,reef_id %in% st)$lprophots)
				mt=mean(subset(sst,reef_id %in% st)$meanT)
				ph=mean(subset(sst,reef_id %in% st)$prophots)
				dec[[cc]]=data.frame(cbind(lon,lat,ye,pe,d850,d8100,d8200,lph,mt))
				cc=cc+1
			}
	}
}

dec=do.call(rbind,dec)
str(dec)
dec$pe[dec$pe<1e-3]=1e-3

# ---- map colored by pr05, with Bruno locations as crosses (Fig. 3a)

ggplot() + 
geom_polygon(data = map.polygon, aes(x = long, y = lat, group = group),lwd=0.1,fill="grey80",col="black") +
 xlim(min(sst$lon)-1, max(sst$lon)+1) + 
  ylim(min(sst$lat)-5, max(sst$lat)+3) + 
  coord_equal() +   theme_void() +
  geom_point(data=sst, aes(x=lon, y=lat,color=lprophots,cex=km2))+
  geom_point(data=dec, aes(x=lon, y=lat,size=1),color="cyan",pch=3,cex=2)+
  scale_color_viridis(option="B")

# ---- our results vs bruno results 

# observed cover change vs pr05
summary(lm(ye~lph,dec))
ggplot(dec,aes(10^(lph),ye),color="grey80")+geom_point()+theme_bw()+xlab("pr05")+ylab("real cover\nchange, % / y")+scale_x_log10()+geom_smooth(method="lm")

# observed cover change vs meanTemperature
summary(lm(ye~mt,dec))
ggplot(dec,aes(mt,ye),color="grey80")+geom_point()+theme_bw()+xlab("meanT")+ylab("real cover\nchange, % / y")+geom_smooth(method="lm")

# observed cover change vs modeled cover change, after 50 years of warming
summary(lm(ye~d850,dec))
ggplot(dec,aes(d850,ye),color="grey80")+geom_point()+theme_bw()+xlab("50yrs of RCP 8.5")+ylab("real cover\nchange, % / y")+geom_smooth(method="lm")

# observed cover change vs modeled cover change, after 100 years of warming
summary(lm(ye~d8100,dec))
ggplot(dec,aes(d8100,ye),color="grey80")+geom_point()+theme_bw()+xlab("100yrs of RCP 8.5")+ylab("real cover\nchange, % / y")+geom_smooth(method="lm")

