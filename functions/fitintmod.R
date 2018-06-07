#' Fit Bayesian spatial regression to body condition data
#' This is only a wrapper function that takes only one input.
#' @indat= input stack data for INLA model
#' @quad= linear or quadratic terms?
#' @ranef= random effects character vector



fitintmod<-function(indat,ranef,inplot){
    # variables
    var.l<-c("NDVI_1m","NDVI_3m","NDVI_12m","NDVI_24m","NDVI_36m")
    
#----numerical results
    # construct mesh
    mesh1<-inla.mesh.2d(as.matrix(indat[,c("e","n")]),max.edge =20,cutoff=40)
    # define weight factors
    A5<-inla.spde.make.A(mesh1,loc=as.matrix(indat[,c("e","n")]))
    # define the spde
    spde<-inla.spde2.matern(mesh1,alpha=2)
    # define spatial field
    w.index<-inla.spde.make.index(name="w",n.spde = spde$n.spde,
                                  n.group=1,n.repl=1)
    # define the stack
    stackfit<-inla.stack(tag="fit",data=list(y=indat$ScaledMassIndex),
                         A=list(1,1,A5),effects=list(Intercept=rep(1,dim(indat)[1]),
                                               X=indat[,c(names(indat))],w=w.index))
    # store results
    results<-NULL
    results.post<-NULL
    results.res<-NULL
    predictions<-NULL
    # fit models
    for (i in 1:length(var.l)){
        print(paste("variable of interest-->",var.l[i],sep=""))
        # create list of variables
        # formula
        formfit<-as.formula(paste('y~-1+Intercept+',paste(var.l[i],"*PA_cover+Year",sep=""),
                                  ranef,sep=""))
        # fit model
        mod.inla<-try(inla(formfit,family="gaussian",data=inla.stack.data(stackfit),
                           control.compute = list(dic=TRUE),
                           control.predictor=list(A=inla.stack.A(stackfit),compute=TRUE)))
        if(class(mod.inla)=="try-error"){
            next
        } else {
            # summary of fixed effects
            fixeddf<-data.frame(mod.inla$summary.fixed,variable=row.names(mod.inla$summary.fixed),
                                DIC=mod.inla$dic$dic,model=paste(var.l[i],"*PA_cover+Year",sep=""),
                                modelset=var.l[i])
            results<-rbind(fixeddf,results)
            # posterior distributions
            for (h in 1:length(mod.inla$marginals.fixed)){
                tmpmarg<-data.frame(mod.inla$marginals.fixed[[h]],
                                    var.name=names(mod.inla$marginals.fixed)[h],
                                    model=paste(var.l[i],"*PA_cover+Year",sep=""),modelset=var.l[i])
                results.post<-rbind(tmpmarg,results.post)
                
            }
            # residuals 
            I0 = inla.stack.index(stackfit, "fit")
            E0 = mod.inla$summary.fitted.values[I0$data,]
            resdf<-data.frame(res=stackfit$data$data$y-E0$mean,
                              model=paste(var.l[i],"*PA_cover+Year",sep=""),
                              modelset=var.l[i])
            results.res<-rbind(resdf,results.res)
        }
}
    #----interaction plot values
    coordinates(inplot)<-~X+Y
    proj4string(inplot)<-latlon
    inplot1<-spTransform(inplot,CRS(ml))
    as.data.frame(coordinates(inplot1)) %>%
        rename(e=X,n=Y) ->coorsml
    # select right columns 
    data.frame(coorsml,as.data.frame(inplot1)) %>%
        mutate(Year_sc=as.numeric(scale(Year))) %>%
        dplyr::select(ScaledMassIndex,PA_cover,NDVI_1m,NDVI_3m,NDVI_12m,NDVI_24m,NDVI_36m,Year_sc,e,n,SiteID) %>%
        mutate(e=e/1000,n=n/1000) ->inplot1
# store rescaling information for variables that remain the same
    resp.att1<-attr(scale(inplot1$ScaledMassIndex),'scaled:scale')
    resp.att2<-attr(scale(inplot1$ScaledMassIndex),'scaled:center')
    pa.att1<-attr(scale(inplot1$PA_cover),'scaled:scale')
    pa.att2<-attr(scale(inplot1$PA_cover),'scaled:center')
    # rescale variables that remain the same
    inplot1 %>%    
        mutate(PA_cover=as.numeric(scale(PA_cover))) %>%
        mutate(ScaledMassIndex=as.numeric(scale(ScaledMassIndex))) -> inplot1
    
for (i in 1:length(var.l)){
    # rescaling info for variable i
    var.i.att1<-attr(scale(inplot1[,var.l[i]]),'scaled:scale')
    var.i.att2<-attr(scale(inplot1[,var.l[i]]),'scaled:center')
    # rescale variable i
    inplot1[,var.l[i]]<-as.numeric(scale(inplot1[,var.l[i]]))
    # construct mesh
    mesh1<-inla.mesh.2d(as.matrix(inplot1[,c("e","n")]),max.edge =20,cutoff=40)
    # define weight factors
    A5<-inla.spde.make.A(mesh1,loc=as.matrix(inplot1[,c("e","n")]))
    # define the spde
    spde<-inla.spde2.matern(mesh1,alpha=2)
    # define spatial field
    w.index<-inla.spde.make.index(name="w",n.spde = spde$n.spde,
                                  n.group=1,n.repl=1)
    # define the model fitting stack
    stackfit<-inla.stack(tag="fit",data=list(y=inplot1$ScaledMassIndex),
                         A=list(1,1,A5),effects=list(Intercept=rep(1,dim(inplot1)[1]),
                            X=inplot1[,c("Year_sc","PA_cover","SiteID",var.l)],w=w.index))
        # create prediction dataframe
        preddf<-data.frame(PA_cover=seq(from=min(inplot1$PA_cover),to=max(inplot1$PA_cover),
        length.out = 30),Year_sc=seq(from=min(inplot1$Year_sc),to=max(inplot1$Year_sc),
        length.out = 30),varint=seq(from=min(inplot1[,var.l[i]]),to=max(inplot1[,var.l[i]]),
        length.out = 30))
        colnames(preddf)[3]<-var.l[i]
        preddf1<-expand.grid(preddf)
        # define the model prediction stack
        stackpred<-inla.stack(tag="predict",data=list(y=NA),
                             A=list(1,1),effects=list(Intercept=rep(1,dim(preddf1)[1]),
                            X=preddf1))
        # combine stacks together
        all.stacks<-inla.stack(stackfit,stackpred)
        # formula
        formfit<-as.formula(paste('y~-1+Intercept+',paste("Year_sc+",var.l[i],"+",
                var.l[i],"*PA_cover",sep=""),ranef,sep=""))
        # fit model
        mod.inla<-try(inla(formfit,family="gaussian",data=inla.stack.data(all.stacks),
                           control.compute = list(dic=TRUE),
                           control.predictor=list(A=inla.stack.A(all.stacks),compute=TRUE)))

        # fitted
        I0 = inla.stack.index(all.stacks, "predict")
        E0 = mod.inla$summary.fitted.values[I0$data,]
        respred<-data.frame(pred=E0$mean,preddf1,modelset=var.l[i])
        # backtransform everything on original scale
        respred$pred<-respred$pred*resp.att1+resp.att2
        respred$PA_cover<-respred$PA_cover*pa.att1+pa.att2
        respred[,var.l[i]]<-respred[,var.l[i]]*var.i.att1+var.i.att2
        colnames(respred)[4]<-c("NDVI")
        # bind results!
        predictions<-rbind(respred,predictions)
    } 
    globresults<-vector("list",4)
    names(globresults)<-c("coefs","posteriors","residuals","predictions")
    globresults[[1]]<-results
    globresults[[2]]<-results.post
    globresults[[3]]<-results.res
    globresults[[4]]<-predictions
    return(globresults)
} # end of function    
