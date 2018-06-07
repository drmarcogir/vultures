#' Fit Bayesian spatial regression to body condition data
#' This is only a wrapper function that takes only one input.
#' @indat= input stack data for INLA model
#' @quad= linear or quadratic terms?
#' @ranef= random effects character vector

scattervalues<-function(indat,ranef){
    # store results
    slopeintres<-NULL
    fitted<-NULL
    var.l<-c("PA_cover","NDVI_1m","NDVI_3m","NDVI_12m","NDVI_24m","NDVI_36m")
    #----prepare data
    coordinates(indat)<-~X+Y
    proj4string(indat)<-latlon
    indat1<-spTransform(indat,CRS(ml))
    as.data.frame(coordinates(indat1)) %>%
        rename(e=X,n=Y) ->coorsml
    # select right columns 
    data.frame(coorsml,as.data.frame(indat1)) %>%
        mutate(Year_sc=as.numeric(scale(Year))) %>%
        dplyr::select(ScaledMassIndex,PA_cover,NDVI_1m,NDVI_3m,NDVI_12m,NDVI_24m,NDVI_36m,Year_sc,e,n,SiteID) %>%
        mutate(e=e/1000,n=n/1000) ->indat1
    # rescale response
    indat1 %>%    
        mutate(ScaledMassIndex=ScaledMassIndex/1000)  -> indat1
    for (i in 1:length(var.l)){
        indat1$original<-indat1[,var.l[i]]
        # rescaling info for variable i
        var.i.att1<-attr(scale(indat1[,var.l[i]]),'scaled:scale')
        var.i.att2<-attr(scale(indat1[,var.l[i]]),'scaled:center')
        # rescale variable i
        indat1[,var.l[i]]<-as.numeric(scale(indat1[,var.l[i]]))
        # construct mesh
        mesh1<-inla.mesh.2d(as.matrix(indat1[,c("e","n")]),max.edge =20,cutoff=40)
        # define weight factors
        A5<-inla.spde.make.A(mesh1,loc=as.matrix(indat1[,c("e","n")]))
        # define the spde
        spde<-inla.spde2.matern(mesh1,alpha=2)
        # define spatial field
        w.index<-inla.spde.make.index(name="w",n.spde = spde$n.spde,
                                      n.group=1,n.repl=1)
        # define the model fitting stack
        stackfit<-inla.stack(tag="fit",data=list(y=indat1$ScaledMassIndex),
                         A=list(1,1,A5),effects=list(Intercept=rep(1,dim(indat1)[1]),
                            X=indat1[,c("Year_sc","PA_cover","SiteID",var.l)],w=w.index))
        
        # create prediction dataframe
        preddf1<-data.frame(varint=seq(from=min(indat1[,var.l[i]]),to=max(indat1[,var.l[i]]),
        length.out = 100))
        colnames(preddf1)[1]<-var.l[i]
        # define the model prediction stack
        stackpred<-inla.stack(tag="predict",data=list(y=NA),
                             A=list(1,1),effects=list(Intercept=rep(1,dim(preddf1)[1]),
                            X=preddf1))
        # combine stacks together
        all.stacks<-inla.stack(stackfit,stackpred)
        
        # formula
        formfit<-as.formula(paste('y~-1+Intercept+',var.l[i],ranef,sep=""))
        # fit model
        mod.inla<-try(inla(formfit,family="gaussian",data=inla.stack.data(all.stacks),
                           control.compute = list(dic=TRUE),
                           control.predictor=list(A=inla.stack.A(all.stacks),compute=TRUE)))
        # extract coefficients
        slopeint1<-data.frame(Int=(mod.inla$summary.fixed[,1][1]-mod.inla$summary.fixed[,1][2]*mean(indat1$original)/sd(indat1$original))*1000,
                              Slope=((mod.inla$summary.fixed[,1][2]/sd(indat1$original))*1000),
                              pred.name=var.l[i])
        slopeintres<-rbind(slopeint1,slopeintres)
        # fitted values
         # fitted
        I0 = inla.stack.index(all.stacks, "predict")
        E0 = mod.inla$summary.fitted.values[I0$data,]
        preddf<-data.frame(pred=E0$mean*1000,min=E0[,3]*1000,max=E0[,5]*1000,
                           predictor=preddf1[,1]*var.i.att1+var.i.att2,pred.name=var.l[i])
        fitted<-rbind(preddf, fitted)
    }
    # put final results together 
       finalres<-vector("list",2)
       names(finalres)<-c("coefs","fitted")
       finalres[[1]]<-slopeintres
       finalres[[2]]<-fitted
       return(finalres)
}