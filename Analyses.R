##########################
# Body condition analyses
##########################

# load required libraries
library(tidyverse);library(nlme);
library(MuMIn);library(stringi)
library(spdep);library(INLA)
library(marcoUtils);library(mgcv)
library(stringi);library(gtools)
library(XLConnect);library(ncf)
library(RColorBrewer)

# source required functions and info for labeling graphs
marcofunctions<-list.files("/mnt/data1tb/Dropbox/Andrea/ndvi/scripts/functions",full.names=TRUE)
for (f in 1:length(marcofunctions)) {source(marcofunctions[f])}

# read in files
lfv <- read.table(file = "/mnt/data1tb/Dropbox/Andrea/ndvi/vulturemodels/LFVdata.txt", header = T, dec = ",")
wbv<- read.table(file = "/mnt/data1tb/Dropbox/Andrea/ndvi/vulturemodels/WBVdata.txt", header = T, dec = ",")
# insert names
lfv$english<-"Lappet-faced vulture"
wbv$english<-"White-backed vulture"

# exclude outliers
lfv<-subset(lfv,NDVI_1m < 0.27)
wbv<-subset(wbv,NDVI_1m < 0.35)

# rescale body condition index
bind_rows(lfv,wbv)    %>%
  dplyr::select(-c(RingCode,Ring,Wing,Mass,Species))  %>%
  group_by(english)  %>%
 mutate(ScaledMassIndex=scale(ScaledMassIndex)) -> dat1

# duplicated coordinates
unique(dat1[,c("X","Y")]) %>%
  mutate(site.no=1:length(X)) %>%
  inner_join(dat1) %>%
  group_by(site.no) %>%
  summarise(dup=n()) %>%
  filter(dup==1) %>%
  inner_join(unique(dat1[,c("X","Y")]) %>%
               mutate(site.no=1:length(X)) %>%
               inner_join(dat1)) %>%
  dplyr::select(-c(dup)) ->dat2

# rescale predictors
dat2 %>% 
  group_by(english) %>%
  tidyr::gather(variable,value,-c(site.no,english,SiteID,X,Y,Date,ScaledMassIndex))  %>%
  group_by(site.no,variable,SiteID,english,X,Y,Date,ScaledMassIndex) %>%
  summarise(value=mean(value)) %>%
  group_by(english,variable) %>%
  mutate(value=scale(value))  %>%
  spread(variable,value) ->dat3

# prepare wbv data
subset(dat3,english=="White-backed vulture") ->wbv1
coordinates(wbv1)<-~X+Y
proj4string(wbv1)<-latlon
wbv1<-spTransform(wbv1,CRS(ml))
as.data.frame(coordinates(wbv1)) %>%
  rename(e=X,n=Y) ->coorsml
data.frame(coorsml,as.data.frame(wbv1)) %>%
  mutate(NDVI_36m2=NDVI_36m^2,NDVI_24m2=NDVI_24m^2,NDVI_12m2=NDVI_12m^2,
         NDVI_3m2=NDVI_3m^2,NDVI_1m2=NDVI_1m^2,NDVI_36m3=NDVI_36m^3,
         NDVI_24m3=NDVI_24m^3,NDVI_12m3=NDVI_12m^3,NDVI_3m3=NDVI_3m^3,
         NDVI_1m3=NDVI_1m^3,PA_cover2=PA_cover^2,PA_cover3=PA_cover^3,
         Year2=Year^2,Year3=Year^3) ->wbv1    

wbv1$e<-wbv1$e/1000
wbv1$n<-wbv1$n/1000

# prepare lfv data
subset(dat3,english=="Lappet-faced vulture") ->lfv1

coordinates(lfv1)<-~X+Y
proj4string(lfv1)<-latlon
lfv1<-spTransform(lfv1,CRS(ml))
as.data.frame(coordinates(lfv1)) %>%
  rename(e=X,n=Y) ->coorsml
data.frame(coorsml,as.data.frame(lfv1)) %>%
  mutate(NDVI_36m2=NDVI_36m^2,NDVI_24m2=NDVI_24m^2,NDVI_12m2=NDVI_12m^2,
         NDVI_3m2=NDVI_3m^2,NDVI_1m2=NDVI_1m^2,NDVI_36m3=NDVI_36m^3,
         NDVI_24m3=NDVI_24m^3,NDVI_12m3=NDVI_12m^3,NDVI_3m3=NDVI_3m^3,
         NDVI_1m3=NDVI_1m^3,PA_cover2=PA_cover^2,PA_cover3=PA_cover^3,
         Year2=Year^2,Year3=Year^3) ->lfv1    

lfv1$e<-lfv1$e/1000
lfv1$n<-lfv1$n/1000

# Lapped-faced vulture
# construct mesh
mesh1<-inla.mesh.2d(as.matrix(lfv1[,c("e","n")]),max.edge =20,cutoff=40)
# define weight factors
A5<-inla.spde.make.A(mesh1,loc=as.matrix(lfv1[,c("e","n")]))
# define the spde
spde<-inla.spde2.matern(mesh1,alpha=2)
# define spatial field
w.index<-inla.spde.make.index(name="w",n.spde = spde$n.spde,
                              n.group=1,n.repl=1)
# define the stack
lfvstack<-inla.stack(tag="fit",data=list(y=lfv1$ScaledMassIndex),
                     A=list(1,1,A5),effects=list(Intercept=rep(1,dim(lfv1)[1]),
                                         X=lfv1[,c(names(lfv1))],w=w.index))

# Fit additive models for Lappet-faced vulture

# linear terms + siteID + gaussian random field
lfvres2<-fitaddmod(indat=lfvstack,ranef='+f(w,model=spde)+f(SiteID,model="iid")',
                   quad=FALSE)

save(lfvres2,file="/mnt/data1tb/Dropbox/Andrea/ndvi/resultstochoose/Rfiles/lfvres2")

# Interaction models for Lappet-faced vulture (PA_cover *NDVI) (=input raw dataframes!)
lfvres3<-fitintmod(indat=lfv1,ranef='+f(w,model=spde)+f(SiteID,model="iid")',
          inplot=lfv)

save(lfvres3,file="/mnt/data1tb/Dropbox/Andrea/ndvi/resultstochoose/Rfiles/lfvres3")

# White-backed vulture
# construct mesh
mesh1<-inla.mesh.2d(as.matrix(wbv1[,c("e","n")]),max.edge =20,cutoff=40)
# define weight factors
A5<-inla.spde.make.A(mesh1,loc=as.matrix(wbv1[,c("e","n")]))
# define the spde
spde<-inla.spde2.matern(mesh1,alpha=2)
# define spatial field
w.index<-inla.spde.make.index(name="w",n.spde = spde$n.spde,
                              n.group=1,n.repl=1)
# define the stack
wbvstack<-inla.stack(tag="fit",data=list(y=wbv1$ScaledMassIndex),
                     A=list(1,1,A5),effects=list(Intercept=rep(1,dim(wbv1)[1]),
                                                 X=wbv1[,c(names(wbv1))],w=w.index))

# linear terms + siteID + gaussian random field
wbvres2<-fitaddmod(indat=wbvstack,ranef='+f(w,model=spde)+f(SiteID,model="iid")',
                   quad=FALSE)
save(wbvres2,file="/mnt/data1tb/Dropbox/Andrea/ndvi/resultstochoose/Rfiles/wbvres2")

# ---Coefficient tables
wbvres2[[1]] %>%
    group_by(modelset) %>%
    filter(variable!="Intercept")  %>%
    dplyr::select(c(mean,variable,modelset,X0.025quant,X0.975quant)) %>%
    rename(lowerCI=X0.025quant,higherCI=X0.975quant,model.name=modelset) %>%
    mutate(Species="White-backed vulture") %>%
    dplyr::select(Species,model.name,variable,mean,lowerCI,higherCI) ->lin.coef.wbv

writeWorksheetToFile(data=lin.coef.wbv,file="/mnt/data1tb/Dropbox/Andrea/ndvi/resultstochoose/coefs/coefs1.xlsx",
                     sheet = "Sheet1", header = TRUE,startCol=1,
                     startRow=1,styleAction =XLC$"STYLE_ACTION.NONE")

# lappet-faced v.: additive models 
lfvres2[[1]] %>%
    group_by(modelset) %>%
    filter(variable!="Intercept")  %>%
    dplyr::select(c(mean,variable,modelset,X0.025quant,X0.975quant)) %>%
    rename(lowerCI=X0.025quant,higherCI=X0.975quant,model.name=modelset) %>%
    mutate(Species="Lappet-faced vulture") %>%
    dplyr::select(Species,model.name,variable,mean,lowerCI,higherCI) ->lin.coef.lfv

writeWorksheetToFile(data=lin.coef.lfv,file="/mnt/data1tb/Dropbox/Andrea/ndvi/resultstochoose/coefs/coefs1.xlsx",
                     sheet = "Sheet1", header = TRUE,startCol=1,
                     startRow=16,styleAction =XLC$"STYLE_ACTION.NONE")

# lappet-faced v.: interaction models
lfvres3[[1]] %>%
    group_by(modelset) %>%
    filter(variable!="Intercept")  %>%
    dplyr::select(c(mean,variable,modelset,X0.025quant,X0.975quant)) %>%
    rename(lowerCI=X0.025quant,higherCI=X0.975quant,model.name=modelset) %>%
    mutate(Species="Lappet-faced vulture") %>%
    dplyr::select(Species,model.name,variable,mean,lowerCI,higherCI) ->int.coef.lfv

writeWorksheetToFile(data=int.coef.lfv,file="/mnt/data1tb/Dropbox/Andrea/ndvi/resultstochoose/coefs/coefs1.xlsx",
                     sheet = "Sheet1", header = TRUE,startCol=8,
                     startRow=1,styleAction =XLC$"STYLE_ACTION.NONE")

# --- Spatial autocorrelation analysis
# additive models 
spautowrap(indat=lfvres2[[3]],group="model",coor=lfv1[,c("e","n")])  %>%
    mutate(species="Lappet-faced vulture") %>%
    bind_rows(spautowrap(indat=wbvres2[[3]],group="model",coor=wbv1[,c("e","n")]) %>%
                  mutate(species="White-backed vulture")) %>%
    ggplot(data=.,aes(x=distance,y=correlation))+facet_grid(model~species)+
    geom_line(size=0.9)+theme_bw()+ylim(-1,1)+xlab("Distance (km)")+
    ylab("Correlation")+theme(text=element_text(size=12,colour="black"),axis.text=element_text(colour="black"))+
    theme(strip.text.x=element_text(size=12,face="bold"),
          strip.text.y=element_text(size=6.5,face="bold")) ->additivesac

ggsave(additivesac,file="/mnt/data1tb/Dropbox/Andrea/ndvi/resultstochoose/figuresNOGRF/additivesac.png",
       width=6,height=8,dpi=400)

# Interaction models
spautowrap(indat=lfvres3[[3]],group="model",coor=lfv1[,c("e","n")])  %>%
    inner_join(modlookup2) %>%
    mutate(species="Lappet-faced vulture") %>%
    ggplot(data=.,aes(x=distance,y=correlation))+facet_grid(model.lab~species,scale="fixed")+
    geom_line(size=0.9)+theme_bw()+ylim(-1,1)+xlab("Distance (km)")+
    ylab("Correlation")+theme(text=element_text(size=12,colour="black"),axis.text=element_text(colour="black"))+
    theme(strip.text.x=element_text(size=12,face="bold"),
          strip.text.y=element_text(size=8,face="bold")) ->intsac

# combine plots together and save results
combined<-plot_grid(additivesac,intsac,ncol = 2)
ggsave(combined,file="/mnt/data1tb/Dropbox/Andrea/ndvi/resultstochoose/figures/SACall.png",
       width=9,height=8,dpi=400)

# --- Posterior distributions
# Lappet-faced vulture
lfvres2[[2]] %>%
    filter(var.name!="Intercept") %>%
    inner_join(modlookup) %>%
    inner_join(varlookup)  %>%
    rename(variable.lab=label) %>%
bind_rows(lfvres3[[2]] %>%
    inner_join(varlookup1) %>%
    rename(variable.lab=model.lab) %>%
    inner_join(modlookupint))  %>%
    mutate(model.lab=factor(model.lab,levels=c(unique(model.lab)))) %>%
ggplot(data=.,aes(x=x,y=y))+facet_wrap(model.lab~variable.lab,scale="free",ncol=4)+
    geom_line(size=0.9)+theme_bw()+geom_vline(xintercept =0, linetype="dotted")+
    ylab("Density")+xlab("Predictor")+theme(strip.text=element_text(face="bold",size=9,
    colour="black"),axis.text =element_text(face="bold",size=10,colour="black"),
    axis.title=element_text(face="bold",size=12,colour="black")) ->post.lfv
# save plot
ggsave(post.lfv,file="/mnt/data1tb/Dropbox/Andrea/ndvi/resultstochoose/figures/posteriorsLFV.png",
       width=8,height=8.5,dpi=400)

# White-backed vulture
wbvres2[[2]] %>%
    filter(var.name!="Intercept") %>%
    inner_join(modlookup) %>%
    inner_join(varlookup) %>%
    ggplot(data=.,aes(x=x,y=y))+facet_wrap(model.lab~label,scale="free")+
    geom_line(size=0.9)+theme_bw()+geom_vline(xintercept =0, linetype="dotted")+
    ylab("Density")+xlab("Predictor")+theme(strip.text=element_text(face="bold",size=10,
    colour="black"),axis.text =element_text(face="bold",size=12,colour="black"),
    axis.title=element_text(face="bold",size=12,colour="black")) ->post.wbv.add
# save plot
ggsave(post.wbv.add,file="/mnt/data1tb/Dropbox/Andrea/ndvi/resultstochoose/figures/posteriorsWBVadditivemodels.png",
       width=8,height=8.5,dpi=400)

# ---scatterplots 

# create dataframe with values
lfv.sct<-scattervalues(indat=lfv,ranef='+f(w,model=spde)+f(SiteID,model="iid")')
wbv.sct<-scattervalues(indat=wbv,ranef='+f(w,model=spde)+f(SiteID,model="iid")')

# PA cover
# fitted values 
lfv.sct[[2]] %>%
    filter(pred.name=="PA_cover") %>%
    mutate(species="Lappet-faced vulture") %>%
    bind_rows(wbv.sct[[2]] %>%
                  filter(pred.name=="PA_cover") %>%
                  mutate(species="White-backed vulture")) -> pa.fitted
# coefficients
lfv.sct[[1]] %>%
    filter(pred.name=="PA_cover") %>%
    mutate(species="Lappet-faced vulture") %>%
    bind_rows(wbv.sct[[1]] %>%
                  filter(pred.name=="PA_cover") %>%
                  mutate(species="White-backed vulture")) -> pa.coefs
# raw data
lfv %>%
    mutate(species="Lappet-faced vulture") %>%
    bind_rows(wbv %>% 
                  mutate(species="White-backed vulture")) %>%
    dplyr::select(species,ScaledMassIndex,SiteID,PA_cover) ->pa.raw

# scatterplot
ggplot(data=pa.fitted,aes(x=predictor,y=pred))+
            geom_point(data=pa.raw,aes(y=ScaledMassIndex,x=PA_cover),
            size=4,col="black",alpha=0.3)+theme_bw()+
            geom_abline(data=pa.coefs,aes(intercept=Int,slope=Slope),size=1)+
            geom_ribbon(aes(ymin=min,ymax=max),linetype=2,alpha= 0.3)+
            facet_grid(~species)+
theme(strip.text.x = element_text(size=15, face="bold"),
      strip.text.y = element_text(size=12, face="bold"))+
    scale_x_continuous(breaks = seq(from=0,to=100, by = 25))+
    theme(text = element_text(size=15),axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black")) +xlab("PA cover")+
    ylab("Scaled Body Mass Index")->pa.scat
# save plot
ggsave(filename="/mnt/data1tb/Dropbox/Andrea/ndvi/resultstochoose/figures/MassPAJune18.png",
       plot =pa.scat,width=10,height=8,dpi=400)

# NDVI
# fitted values
lfv.sct[[2]] %>%
    filter(pred.name!="PA_cover") %>%
    mutate(species="Lappet-faced vulture") %>%
    bind_rows(wbv.sct[[2]] %>%
            filter(pred.name!="PA_cover") %>%
            mutate(species="White-backed vulture")) %>%
        left_join(data.frame(pred.name=c("NDVI_1m","NDVI_3m","NDVI_12m","NDVI_24m","NDVI_36m"),
            ndvilabel=factor(c("1 month","3 months","1 year"," 2 years"," 3 years"),
            levels=c("1 month","3 months","1 year"," 2 years"," 3 years"))))-> ndvi.fitted
# coefficients
lfv.sct[[1]] %>%
    filter(pred.name!="PA_cover") %>%
    mutate(species="Lappet-faced vulture") %>%
    bind_rows(wbv.sct[[1]] %>%
    filter(pred.name!="PA_cover") %>%
    mutate(species="White-backed vulture")) %>%
    left_join(data.frame(pred.name=c("NDVI_1m","NDVI_3m","NDVI_12m","NDVI_24m","NDVI_36m"),
        ndvilabel=factor(c("1 month","3 months","1 year"," 2 years"," 3 years"),
    levels=c("1 month","3 months","1 year"," 2 years"," 3 years"))))-> ndvi.coefs
# raw data
lfv %>%
    mutate(Species="Lappet-faced vulture") %>%
    bind_rows(wbv %>% 
                  mutate(Species="White-backed vulture")) %>%
    dplyr::select(english,ScaledMassIndex,SiteID,NDVI_1m,NDVI_3m,NDVI_12m,NDVI_24m,NDVI_36m,PA_cover) %>%
    tidyr::gather(ndviname,ndvivalue,-c(english,ScaledMassIndex,SiteID))  %>%
    inner_join(data.frame(ndviname=c("NDVI_1m","NDVI_3m","NDVI_12m","NDVI_24m","NDVI_36m"),
                          ndvilabel=factor(c("1 month","3 months","1 year"," 2 years"," 3 years"),
                                           levels=c("1 month","3 months","1 year"," 2 years"," 3 years")))) %>%
    dplyr::rename(species=english)->combined.m2

# scatterplot
ggplot(data=ndvi.fitted,aes(x=predictor,y=pred))+
            geom_point(data=combined.m2,aes(y=ScaledMassIndex,x=ndvivalue),
            size=2,col="black",alpha=0.1)+theme_bw()+
    facet_grid(ndvilabel~species,scales="free",space="free")+
            geom_abline(data=ndvi.coefs,aes(intercept=Int,slope=Slope),size=0.8)+
            geom_ribbon(aes(ymin=min,ymax=max),linetype=2,alpha= 0.3)+
    theme(text = element_text(size=15),axis.text = element_text(colour="black"))+
theme(strip.text= element_text(size=12, face="bold"))+ylab("Scaled Mass Index")+xlab("NDVI")->ndvi.scat
# save plot 
ggsave(filename="/mnt/data1tb/Dropbox/Andrea/ndvi/resultstochoose/figures/MassNDVIJune18.png",
       plot =ndvi.scat,width=7,height=8,dpi=400)

# ---Fancy interaction plot (NDVI*PA_cover)
cols <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))
lfvres3[[4]] %>%
    inner_join(modlookup1)  %>%
    ggplot(data=.,aes(x=NDVI,y=PA_cover,z=pred))+
    geom_raster(aes(fill=pred))+facet_wrap(~model.lab,scale="free",ncol=2)+
    scale_fill_gradientn(colours = cols(30))+theme_bw()+
    labs(fill = "Scaled Body Mass Index")+
    ylab("Protected Area Cover")+xlab("NDVI")->intplot

ggsave(filename="/mnt/data1tb/Dropbox/Andrea/ndvi/resultstochoose/figures/Interaction1.png",
       plot =intplot,width=7.5,height=6.5,dpi=400)