# repace directory location below with the name of the cloned github repo directory
setwd("~/Dropbox/Documents/coral_triangle_SLiM/trianglePaper_finalSubmission/CoralTriangle_SLiM_model/")
load("sst_modified.RData")
#sst$logarea2=log(sst$kms2_reef,10)
mig=as.matrix(read.csv("Amill_CRA_Mmat.csv"))
# head(sst)

#--------------- generating out-migration matrix (for nonWF SLiM)
mign=c()
for (i in 1:nrow(mig)) { mign=cbind(mign,mig[,i]*sst$kms2_reef[i]) }
sumsout=apply(mign,1,sum)
mig2=c()
for (i in 1:nrow(mign)) { mig2=rbind(mig2,mign[i,]/sumsout[i]) }
#save(mig2,sst,file="Amill_CRA_OutMmat_sst_modified.RData")

mig=mig2
minMig=1e-5 # minimal migration rate; values below that will be set to 0

#------ GBR pops selector: if you want to run a small model with just 6 populations somewhere in the GBR. Skip this if redoing the whole thing.

	npops=6
	gbr=which(sst$x_cent<2e+6 & sst$x_cent>700000 & sst$y_cent<(-1.75e+6) & sst$y_cent>(-3e+6))
	sst.gbr=sst[gbr,]
	mig.gbr=mig[gbr,gbr]
    select=sample(1:nrow(sst.gbr),npops,prob=sst.gbr$logarea)
    sst1=sst.gbr[select,]
    mig1=mig.gbr[select,select]
pheatmap(log(mig1+minMig,10))
	sst=sst1
	mig=mig1
	
#------------- setting environment and pop sizes

# temperature
sst$hi=sst$meanT+sst$TST
sst$hi1=sst$meanT+1.5*sst$stdT

# relative reef size
sst$relsize=sst$kms2_reef/max(sst$kms2_reef)

minpop=100 # carrying capacity of smallest reef
psize=sst$kms2_reef
psize=minpop*sst$kms2_reef/min(sst$kms2_reef)

nsites= nrow(sst) # change this to smaller number if you want to run a smaller model
choose=sample(nrow(sst),nsites,prob=sst$relsize)
choose=sort(choose)
migc=as.matrix(mig[choose,choose])
migc[migc<minMig]=0

# ------- making a connectivity figure - map with arcs

migcc=migc
migcc[migc<0.999]=1
migcc[migc<0.1]=0.5
migcc[migc<0.01]=0.25
migcc[migc<1e-3]=0

table(migcc==0)
hist(migcc)

# migration curves
curves=list();cc=1
for (i in 1:(nrow(sst)-1)) {
	for (j in 1:nrow(sst)) {
		if (i==j) { next}	
		if (migcc[i,j]>0) {
			curves[[cc]]=data.frame(x1=sst$lon[i],y1=sst$lat[i],x2=sst$lon[j],y2=sst$lat[j],rate=migcc[i,j])
			cc=cc+1
		}
	}
}
curv=do.call(rbind,curves)

library(rworldmap)
map.polygon <- getMap(resolution = "low")

require(ggplot2)
require(pheatmap)

# pop sizes, temperature, and connectivity map (Fig. 2a)	
conmap=ggplot() + 
geom_polygon(data = map.polygon, aes(x = long, y = lat, group = group),lwd=0.1,fill="grey80",col="black") +
 xlim(min(sst$lon)-1, max(sst$lon)+1) + 
  ylim(min(sst$lat)-5, max(sst$lat)+3) + 
  coord_equal() +   
  geom_curve(data=curv,curvature=-0.2, aes(x=x1,y=y1,xend=x2,yend=y2),size=curv$rate/4,color="cyan3")+
  geom_point(data=sst[choose,], aes(x=lon, y=lat,color=meanT,cex=km2),alpha=0.5)+
  scale_color_gradient( low="blue",high="coral")+
  theme_void()
 
conmap
  
save(conmap,file="connectivity_map_ggplot.RData")

# -------
# delta-temperature maps for RCP 4.5 and 8.5 (Fig. S1)

# RCP 4.5
ggplot() + 
geom_polygon(data = map.polygon, aes(x = long, y = lat, group = group),lwd=0.1,fill="grey80",col="black") +
 xlim(min(sst$lon)-1, max(sst$lon)+1) + 
  ylim(min(sst$lat)-5, max(sst$lat)+3) + 
  coord_equal() +   
  geom_point(data=sst[choose,], aes(x=lon, y=lat,color=DT_RCP45,cex=km2),alpha=0.5)+
scale_color_gradient( low="blue",high="coral",limits=c(min(sst$DT_RCP45),max(sst$DT_RCP85)))+
  theme_void()

# RCP 8.5
 ggplot() + 
geom_polygon(data = map.polygon, aes(x = long, y = lat, group = group),lwd=0.1,fill="grey80",col="black") +
 xlim(min(sst$lon)-1, max(sst$lon)+1) + 
  ylim(min(sst$lat)-5, max(sst$lat)+3) + 
  coord_equal() +   
  geom_point(data=sst[choose,], aes(x=lon, y=lat,color=DT_RCP85,cex=km2),alpha=0.5)+
scale_color_gradient( low="blue",high="coral",limits=c(min(sst$DT_RCP45),max(sst$DT_RCP85)))+
  theme_void()

#--------- writing down migration matrix and population sizes

psize=sst$rel

write.table(migc,file=paste("tri3_migration_",length(choose),".txt",sep=""),quote=F,row.names=F,sep="\t",col.names=F)
write.table(round(sst[choose,"kms2_reef"],0),file=paste("tri3_popsize_",length(choose),".txt",sep=""),quote=F,row.names=F,sep="\t",col.names=F)

sstc=sst[choose,]

pops=seq(1,nsites) # population names
e.local=sstc$meanT-mean(sstc$meanT) # local environmental settting
sin.p=5 # period of sinusoidal climate fluctuations, in reproduction cycles
# ---- rand /  sin settings:
# e.sd= 0.2 # SD of random environmentsl fluctuations (different in each pop)
sin.amp=0.5 # amplitude of sinusoidal climate fluctuations, common to all pops
# sin.rand=0.2 # amplitude of random climate fluctuations, common to all pops
e.sd= 0.0 # SD of random environmentsl fluctuations (different in each pop)
sin.rand=0.0 # amplitude of random climate fluctuations, common to all pops
burnin=5500 # length of burnin period
gmax=6000 # maximum number of generations

e.increment=sstc$DT_RCP85/90 # increment in environmental setting per reproduction cycle (year) past burn-in period

# for no-fluctuation and sinusoidal model, there is only one environmental profile ("a") for all SLiM replicates; for random model there are four different profiles.
for (index in c("a","b","c","d")) {
	message(index)
	dd=lapply(seq_len(gmax),function(i) {
		newgen=e.local+sin(i*2*pi/sin.p)*0.5*(sin.amp+rnorm(1,0,sin.rand))
		newgen=newgen+rnorm(length(newgen),0,e.sd)
		if (i>burnin){
			newgen=newgen+(i-burnin)*e.increment
		}
		return(round(newgen,3))
		}
	)
	envs=data.frame(do.call(cbind,dd))
	write.table(envs,row.names=F,quote=F,sep="\t",file=paste("tri85sin_",index,"_environment.txt",sep=""))
}

# ------ sanity check: let's plot environmental profile over the onset of warming for a random sample of 5 populations

envplot=function(x, popset=sample(length(choose),5), genset=c(5400:5700)) {
	require(ggplot2)
	E=data.frame(t(x))
	te=stack(E)
	names(te)=c("env","pop")
	te$pop=sub("X","",te$pop)
	te$gen=seq(ncol(x))
	ggplot(subset(te,pop %in% popset & gen %in% genset),aes(gen,env,colour=pop))+geom_line()
}
envplot(envs)

