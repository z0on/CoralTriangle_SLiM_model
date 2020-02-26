

# repace directory location below with the name of the cloned github repo directory
setwd("~/Dropbox/Documents/coral_triangle_SLiM/trianglePaper_finalSubmission/CoralTriangle_SLiM_model/")

library(tidyverse)
library(Rmisc)
library(viridis)
library(gridExtra)
library(rworldmap)

load("Amill_CRA_OutMmat_sst_modified.RData")

# ------ plotting individual coral cover maps (Fig. 2 C-F)

files=c("REv32_8_8_base.RData","REv32_8NM_8_NM.RData")

ress=list();i=1
for (f in files) {
	load(f)
	message(f)
	ress[[length(ress)+1]]=merge(covall,sst,by="reef_id",all.x=T)
	}


decs=c(10,20)
topcover=0.6
#head(sgg)

map.polygon <- getMap(resolution = "low")

plots=list()
for (f in 1:length(ress)) {
	for (d in decs) {
		sgg=subset(ress[[f]],dec==d)
#		head(sgg)
		if (d==0) { showLegend=T} else {showLegend=F}
		showLegend=F
		plots[[length(plots)+1]]=ggplot() +
			geom_polygon(data = map.polygon, aes(x = long, y = lat, group = group),size=0.1,fill="grey80",col="black") +
			xlim(min(sst$lon)-1, max(sst$lon)+1) + 
			ylim(min(sst$lat)-5, max(sst$lat)+3) + 
			coord_equal() + 
			geom_point(data=sgg, aes(lon,lat,color=cover,cex=km2), alpha=0.5,show.legend=showLegend) +
			scale_colour_viridis(option="B",limits=c(0,topcover))+
			theme_void()
			#+ggtitle(paste(files[f],d))
	}
}
pdf("fig2_maps.pdf",width=9,height=9)
do.call("grid.arrange", c(plots, nrow=2))
dev.off()

# ------ plotting movie frames (supplemental movies)

# set input file names. See RawResults/model_run_abbreviations.txt for guidance what each filename identifier (32tri(...)a.clean) means

# currently set for making Supplemental movie 2 (parameter variations)
files=c("32tri8a.clean","32tri8ljma.clean","32tri8lma.clean","32tri8lpa.clean")

readslim=function(filename) {
	g85=read.table(paste("RawResults/",filename,sep=""),sep="\t",skip=19)
	names(g85)=c("g","reef_id","fit","phen","t","nmut","h","age","nad","mort")
	g85$reef_id=g85$reef_id+1
	g85$reef_id=as.factor(g85$reef_id)
	g85=subset(g85,g>5300 & g<=5700)
	g85=merge(g85,sst[,c(1:5,35:37)],by="reef_id",all.x=T)
	gg=g85 %>% mutate(cover=nad/km2)
	gg=gg %>% arrange(g)
	return(gg)
}

ress=list()
for (f in files) {
	message(f)
	ress[[f]]=readslim(f) 
	}

# labels of models runs to appear on movie frames
model=c("main","low juv.mort.","low mut.rate","narrow fit.curve")
library(rworldmap)
library(gridExtra)
map.polygon <- getMap(resolution = "low")

 rescaler = function(x, to = c(0, 1), from = NULL) {
    			 ifelse(x<1.2, 
    			 scales::rescale(x, to = to, from = c(min(x, na.rm = TRUE), 1.2)), 
    			 1)
 }


mxx=c();mnn=c();minn=c();i=1
 for (i in 1:length(files)) { mxx=append(mxx,max(ress[[i]]$cover[1:20])) }
 for (i in 1:length(files)) { mnn=append(mnn,mean(ress[[i]]$cover[1:20])) }
 for (i in 1:length(files)) { minn=append(minn,min(ress[[i]]$cover[1:60])) }
 mnn
 mxx
 topcover= max(mxx)
 for (i in 1:length(files)) {
 	ress[[i]]$cover[ress[[i]]$cover>topcover]=topcover
	 ress[[i]]$cover=rescaler(ress[[i]]$cover)
 }

i=5550
library(viridis)
library(gridExtra)
system("mkdir jpg_jj")
system("mkdir jpg_s")
single.run=1
for (i in c(5480:5700)){
	plo=list()
	txtcol="cyan3"
	if(i>5500) { txtcol="coral"}
	for (f in 1:length(files)) {
		sgg=subset(ress[[files[f]]],g==i)
		plo[[f]]=ggplot() +
			geom_polygon(data = map.polygon, aes(x = long, y = lat, group = group),size=0.1,fill="grey80",col="black") +
			xlim(min(sst$lon)-1, max(sst$lon)+1) + 
			ylim(min(sst$lat)-5, max(sst$lat)+3) + 
			coord_equal() + 
			geom_point(data=sgg, aes(lon,lat,color=cover,size=km2),alpha=0.5) + scale_size_area(max_size=8)+
			scale_colour_viridis(option="B",limits = c(0, topcover), oob = scales::squish)+
			theme_void()+
			annotate("text",
				label=model[f],
				col="grey20",cex=8,x=max(sgg$lon)-13,y=max(sgg$lat))+
			theme(
			  legend.position = "none",
				plot.margin=unit(c(-0.05,-0.05,-0.1,-0.1), "cm"))
#			theme(panel.background = element_rect(fill = "grey40",colour = "grey40",size = 0.5, linetype = "solid"))

		
		if (f==1) { 	plo[[f]]=plo[[f]]+annotate("text",
				label=i-5500,
				col=txtcol,cex=7,x=min(sgg$lon)+1,y=max(sgg$lat))
		}
		if (f==single.run) {
	
			jpeg(filename=paste("jpg_s/",i,".jpg",sep=""),width=320,height=320,quality=90)
#			pdf(file=paste("jpg_s/",i,".pdf",sep=""))
			plot(plo[[f]]+annotate("text",
				label=i-5500,
				col=txtcol,cex=7,x=min(sgg$lon)+5,y=max(sgg$lat)-1) )
			dev.off()
		}
	}
	jpeg(filename=paste("jpg_jj/",i,".jpg",sep=""),width=960,height=960)
#	pdf(file=paste("jpg_jj/",i,".pdf",sep=""))
	grid.arrange(plo[[1]],plo[[2]],plo[[3]],plo[[4]],ncol=2,heights=c(1,1))
	dev.off()
}

#--------- map of pr05 (fig. 3 A without Bruno data sites)

names(sgg)
ggplot() +
			geom_polygon(data = map.polygon, aes(x = long, y = lat, group = group),size=0.1,fill="grey80",col="black") +
			xlim(min(sst$lon)-1, max(sst$lon)+1) + 
			ylim(min(sst$lat)-5, max(sst$lat)+3) + 
			coord_equal() + 
			geom_point(data=sst, aes(lon,lat,color=lprophots,cex=km2), alpha=0.5) +
			scale_colour_viridis(option="B")+
			theme_void()


# ----- variation in coral cover (Fig. S2)

names(sst)
readslim=function(filename) {
	g85=read.table(filename,sep="\t",skip=19)
	names(g85)=c("g","reef_id","fit","phen","t","nmut","h","age","nad","mort")
	g85$reef_id=g85$reef_id+1
	g85$reef_id=as.factor(g85$reef_id)
	g85=subset(g85,g>5300 & g<=5700)
	g85=merge(g85,sst[,c(1:5,35:37)],by="reef_id",all.x=T)
	gg=g85 %>% mutate(cover=nad/km2)
	gg=gg %>% arrange(g)
	return(gg)
}

res1=readslim("RawResults/32tri8a.clean")
res1$cover=apply(cbind(res1$cover,res2$cover,res3$cover),1,mean)
res1$g=res1$g-5500
p=1
sp=subset(res1,reef_id==p)
sp$cover[sp$cover<=0.01]=0.01
lo=loess(cover~g,sp,span=0.25)
sp$Wlo=predict(lo,newdata=sp)
sp$Wlo[sp$Wlo<=0.01]=0.01
plot(cover~g,sp,type="l")
lines(Wlo~g,sp,col="red")

plot((sp$cover/sp$Wlo-1)~sp$g,type="l",col=rgb(0,0,0,alpha=0.03),ylim=c(-1,1),xlim=c(-80,80), xlab="year of warming",ylab="coral cover deviation",mgp=c(2.3,1,0))
for (p in 2:680){
	sp=subset(res1,reef_id==p)
	sp$cover[sp$cover<=0.01]=0.01
	lo=loess(cover~g,sp,span=0.25)
	sp$Wlo=predict(lo,newdata=sp)
	sp$Wlo[sp$Wlo<=0.01]=0.01
	lines((sp$cover/sp$Wlo-1)~sp$g,type="l",col=rgb(0,0,0,alpha=0.03))
}
abline(v=0,lty=3,col="grey70")




