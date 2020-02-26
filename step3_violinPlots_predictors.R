# repace directory location below with the name of the cloned github repo directory
setwd("~/Dropbox/Documents/coral_triangle_SLiM/trianglePaper_finalSubmission/CoralTriangle_SLiM_model/")

library(vegan)
library(tidyverse)

load('sst_modified.RData')

files=as.character(read.table("resfiles")[,1])
meta=do.call(rbind,strsplit(files,"[_.]"))[,3:4]
row.names(meta)=paste(meta[,1],meta[,2],sep=".")
ids=row.names(meta)
gvar=list();pw=list();ages=list();cover=list();resp=list()
for (a in 1:length(files)) { 
	ll=load(files[a])
# averaging coral cover (raws), genetic variaiton (hs), and adult age (ages) across replicates 
	cover[[a]]=covall
	maxdec=max(covall$dec)
	if (a>1 & maxdec<20){
# adding empty cover values if the simulation was shorter than 200 years in future
    adds=subset(cover[[1]],dec>maxdec)
    adds$cover=0
	  cover[[a]]=rbind(cover[[a]],adds)
	  cover[[a]]=cover[[a]][order(cover[[a]]$dec,cover[[a]]$reef_id),]
	}
	# correcting cover measure for "ps" setting, with twice larger carrying capacity
	if (meta[a,2] == "ps") { cover[[a]][,"cover"]=cover[[a]][,"cover"]/2 }
	ages[[a]]=ageall
	covs=data.frame(t(cover[[a]] %>% dplyr::select(reef_id,dec,cover) %>% spread(reef_id,cover) %>% dplyr::select(-dec)))
# pw: pre-warming
	pw[[a]]=apply(covs[,1:20],1,mean)
# changes in cover at 100 and 200 years of warming
	resp[[a]]=data.frame(cbind(r100=(pw[[a]]-covs[,30])/pw[[a]], r200=(pw[[a]]-covs[,40])/pw[[a]]))
}
names(resp)=ids

# ------ pre-warming cover vs pr05 (Fig.S4b)

dat=merge(cover[[which(ids=="8.base")]],sst,by="reef_id",all.x=T)
dat0=subset(dat,dec==0)
str(dat0)

ggplot(dat0,aes(prophots,cover))+geom_point(alpha=0.3)+geom_smooth()+scale_x_log10()+theme_bw()+xlab("pr05")+ylab("pre-warming\ncoral cover")+ylim(0,0.8)

# ---------- violin plots of coral cover (Fig 1)

toplot=c(1:4,8:15)
#toplot=c(1:nrow(meta))
plots=list()
i=1;a=4
for (a in toplot) {
	dat=merge(cover[[a]],sst,by="reef_id",all.x=T)
	head(dat)
	datt=subset(dat,dec %in% c(0,5,10,15,20))
# removing top and bottom 10 reefs in each decade
	datt2=c();d=10
	for (d in c(0,5,10,15,20)) {
		dd=subset(datt,dec==d)
		datt2=data.frame(rbind(datt2,dd[order(dd$cover)[11:670],]))
	}
	datt2$dec=factor(datt2$dec)
	plots[[i]]=ggplot(datt2,aes(dec,cover))+geom_violin(scale="width")+ggtitle(ids[a])+ylim(0,0.7)+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())
	i=i+1
}
library(gridExtra)
pdf(file="violins2.pdf",height=1.5,width=13)
do.call("grid.arrange", c(plots, nrow=1))
dev.off()

# ---------- violin plots of age (not used)

toplot=c(1:4,8:13)
#toplot=c(1:nrow(meta))
plots=list()
i=1;a=4
for (a in toplot) {
  dat=merge(ages[[a]],sst,by="reef_id",all.x=T)
  head(dat)
  datt=subset(dat,dec %in% c(0,5,10,15,20))
  # removing top and bottom 10 reefs in each decade
  datt2=c();d=10
  for (d in c(0,5,10,15,20)) {
    dd=subset(datt,dec==d)
    datt2=data.frame(rbind(datt2,dd[order(dd$age)[11:670],]))
  }
  datt2$dec=factor(datt2$dec)
  head(datt2)
  plots[[i]]=ggplot(datt2,aes(dec,age))+geom_violin(scale="width")+ggtitle(ids[a])+ylim(4,12.5)+theme(axis.title.x=element_blank(),axis.title.y=element_blank())
  i=i+1
}
library(gridExtra)
pdf(file="violins_age.pdf",height=1.5,width=12)
do.call("grid.arrange", c(plots, nrow=1))
dev.off()


# ----- predictors

calc=ids[-c(1,14,15)]
ves1=c();ves2=c();ves1.0=c();ves2.0=c();nves1=c();nves1.0=c();nves2.0=c();nves2=c();mm=c();plo=list();mm0=c()
i=1;m=2
par(mfrow=c(2,3))
for (m in which(ids %in% calc)) {
	print(ids[m])
	
	# calculating reef response relative to pre-warming state
	c20=subset(cover[[m]],dec==20)$cover
	c10=subset(cover[[m]],dec==10)$cover
	c0=subset(cover[[m]],dec==-1)$cover
	rr=data.frame(cbind(r100=(c10-c0)/c0, r200=(c20-c0)/c0))
	dat=cbind(sst,rr)
	dat=dat[dat$r100!=Inf & !is.na(dat$r100) & dat$r200!=Inf & !is.na(dat$r200),]
	dat$r100[dat$r100>1]=1
	dat$r200[dat$r200>1]=1
	dat$dt=dat$DT_RCP85	
	
	# binary response variable: does the reef decline below 50% of original cover?
	dat$p100=as.numeric(dat$r100>-0.5) # (after 100 years)
	dat$p200=as.numeric(dat$r200>-0.5) # (after 200 years)

	# making plots of correspondence between the rate of reef decline and predictors
	nbins=20 # number of bins to put reefs in, depending on predictor variable 
	props1=props2=c();lpwin=seq(-3,0,length.out=nbins)  # for lprophots, which is the same as pr05
#	props1=props2=c();lpwin=seq(min(dat$meanT),max(dat$meanT),length.out=nbins) # for present-day temperature
	for (k in 1:(length(lpwin)-1)) {
		sp=subset(dat,lprophots>=lpwin[k] & lprophots<lpwin[k+1]) # for pr05
#		sp=subset(dat,meanT>=lpwin[k] & meanT<lpwin[k+1]) # for temperature
		props1=append(props1,length(which(sp$r100>-0.5))/nrow(sp))
		props2=append(props2,length(which(sp$r200>-0.5))/nrow(sp))
	}
	prop.declined=data.frame(cbind(lp=lpwin[1:length(lpwin)-1],props1,props2))
	# declines over 100 years are cyan, declines over 200 years are red
	plo[[length(plo)+1]]=ggplot(prop.declined) + geom_point(aes(x=lp, y=props1),colour="cyan3",alpha=0.5) + geom_smooth(aes(x=lp, y=props1),method = "glm", method.args = list(family = "binomial"), se = FALSE,colour="cyan3") +geom_point(aes(x=lp, y=props2),colour="coral",alpha=0.5) + geom_smooth(aes(x=lp, y=props2),method = "glm", method.args = list(family = "binomial"), se = FALSE,colour="coral") +theme_bw()+ggtitle(meta[m,2])+ylim(0,1)

# all three precictors, after 100 years
head(dat)
	aa=anova(glm(p100~lprophots+meanT+dt,family="binomial",dat))
	print(summary(glm(p100~lprophots+meanT+dt,family="binomial",dat)))
	preds1=row.names(aa)
	devtot=aa[1,4]
	varexp1=aa$Deviance/devtot

# all three precictors, after 200 years
	aa=anova(glm(p200~lprophots+meanT+dt,family="binomial",dat))
	print(summary(glm(p200~lprophots+meanT+dt,family="binomial",dat)))
	preds2=row.names(aa)
	devtot=aa[1,4]
	varexp2=aa$Deviance/devtot

# precictors without pr05, after 100 years
	aa=anova(glm(p100~meanT+dt,family="binomial",dat))
	print(summary(glm(p100~meanT+dt,family="binomial",dat)))
	preds1.0=row.names(aa)
	devtot=aa[1,4]
	varexp1.0=aa$Deviance/devtot

# precictors without pr05, after 100 years
	aa=anova(glm(p200~meanT+dt,family="binomial",dat))
	print(summary(glm(p200~meanT+dt,family="binomial",dat)))
	preds2.0=row.names(aa)
	devtot=aa[1,4]
	varexp2.0=aa$Deviance/devtot

# combining "deviance explained" results
	ves1=append(ves1,varexp1[-1])
	ves1.0=append(ves1.0,varexp1.0[-1])
	nves1=append(nves1,preds1[-1])
	nves1.0=append(nves1.0,preds1.0[-1])
	ves2=append(ves2,varexp2[-1])
	ves2.0=append(ves2.0,varexp2.0[-1])
	nves2=append(nves2,preds2[-1])
	nves2.0=append(nves2,preds2.0[-1])
	mm=append(mm,rep(ids[m],length(varexp1[-1])))
	mm0=append(mm0,rep(meta[m,2],length(varexp1.0[-1])))
}
library(gridExtra)
do.call("grid.arrange", c(plo, nrow=3))

# putting together "variance explained" table

vexp1=data.frame(cbind(var=ves1,f=nves1,m=mm))
vexp2=data.frame(cbind(var=ves2,f=nves2,m=mm))
vexp1.0=data.frame(cbind(var=ves1.0,f=nves1.0,m=mm0))
vexp2.0=data.frame(cbind(var=ves2.0,f=nves1.0,m=mm0))
vexp1$var=as.numeric(as.character(vexp1$var))
vexp1.0$var=as.numeric(as.character(vexp1.0$var))
vexp2$var=as.numeric(as.character(vexp2$var))
vexp2.0$var=as.numeric(as.character(vexp2.0$var))

# setting colors - same color for same factor

library(RColorBrewer)
myColors <- brewer.pal(4,"BrBG")
names(myColors) <- levels(vexp1$f)
colScale <- scale_fill_manual(name = "f",values = myColors)

# plot deviance explained in all models (not used)

p1=ggplot(subset(vexp1,f!="Residuals"),aes(m,var,fill=f))+geom_bar(stat="identity")+theme_bw()+ylim(0,1)+colScale+theme(axis.text.x = element_text(angle = 45,hjust=1))
p2=ggplot(subset(vexp2,f!="Residuals"),aes(m,var,fill=f))+geom_bar(stat="identity")+theme_bw()+ylim(0,1)+colScale+theme(axis.text.x = element_text(angle = 45,hjust=1))
p1.0=ggplot(subset(vexp1.0,f!="Residuals"),aes(m,var,fill=f))+geom_bar(stat="identity")+theme_bw()+ylim(0,1)+colScale+theme(axis.text.x = element_text(angle = 45,hjust=1))
p2.0=ggplot(subset(vexp2.0,f!="Residuals"),aes(m,var,fill=f))+geom_bar(stat="identity")+theme_bw()+ylim(0,1)+colScale+theme(axis.text.x = element_text(angle = 45,hjust=1))
grid.arrange(p1,p1.0,p2,p2.0,nrow=1)

# plot only the "main setting" model - Fig 3b

baseExpl=data.frame(rbind(subset(vexp1,m=="8.base"),subset(vexp2,m=="8.base"),subset(vexp1.0,m=="base"),subset(vexp2.0,m=="base")))
baseExpl$yr=c(rep("100",length(preds1)-1),rep("200",length(preds1)-1),rep("100",length(preds1.0)-1),rep("200",length(preds1.0)-1))
baseExpl$set=c(rep("with pr05",2*(length(preds1)-1)),rep("without pr05",2*(length(preds1.0)-1)))
ggplot(baseExpl,aes(yr,var,fill=f))+geom_bar(stat="identity")+theme_bw()+colScale+theme(axis.text.x = element_text(angle = 45,hjust=1))+xlab("year of warming")+ylab("deviance explained")+ylim(0,0.7)+facet_wrap(~set)

# ----- pr05 vs meanT plot (Fig S3a)

ggplot(sst,aes(meanT,prophots))+geom_point(alpha=0.3)+theme_bw()+scale_y_log10()+ylab("pr05")+xlab("pre-warm T")

# ----------- how much fitness drops under different fitness curves
 sd=0.5
 dnorm(1,0,sd)/dnorm(0,0,sd)

#----- plotting individual coral cover tracks (not used in the last version of manuscript)

load('sst_modified.RData')
plots=list()
row.names(meta)=c(1:nrow(meta))
toplot=c(2,8,7,5,12,)
toplot=c(1:nrow(meta))

ids
i=1
for (a in toplot) {
	dat=merge(cover[[a]],sst,by="reef_id",all.x=T)
# plots[[a]]=ggplot(dat,aes(dec,cover,group=reef_id))+geom_line(alpha=0.1)+theme_bw()+theme(legend.position="none")+ylim(0,0.75)+xlab("decade")+ylab("coral cover")+ggtitle(ids[a])
 plots[[i]]=ggplot(dat,aes(dec,cover,group=reef_id))+geom_line(alpha=0.03)+theme_bw()+theme(legend.position="none")+ggtitle(ids[a])
 #+ylim(0,0.75)+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())
  	i=i+1
 }
 library(gridExtra)
do.call("grid.arrange", c(plots, nrow=3))


