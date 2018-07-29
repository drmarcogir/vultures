#' Fit Bayesian spatial regression to body condition data
#' This is only a wrapper function that takes only one input.
#' @indat= input stack data for INLA model
#' @quad= linear or quadratic terms?
#' @ranef= random effects character vector


fitaddmod<-function(indat,quad=FALSE,ranef){
    if(quad==TRUE){
        # store results
        results<-NULL
        results.post<-NULL
        results.res<-NULL
        # list of predictors
        var.l<-c("PA_cover","NDVI_1m","NDVI_3m","NDVI_12m","NDVI_24m","NDVI_36m")
        for (i in 1:length(var.l)){
            print(paste("variable of interest-->",var.l[i],sep=""))
            # create list of variables
            var.l1<-c("Year",paste("Year",2,sep=""),var.l[i],paste(var.l[i],2,sep=""))
            # formula
            formfit<-as.formula(paste('y~-1+Intercept+',paste(var.l1,collapse="+"),
                                      ranef,sep=""))
            # fit model
            mod.inla<-try(inla(formfit,family="gaussian",data=inla.stack.data(indat),
                               control.compute = list(dic=TRUE),
                               control.predictor=list(A=inla.stack.A(indat),compute=TRUE)))
            if(class(mod.inla)=="try-error"){
                next
            } else {
                # summary of fixed effects
                fixeddf<-data.frame(mod.inla$summary.fixed,variable=row.names(mod.inla$summary.fixed),
                                    DIC=mod.inla$dic$dic,model=paste(var.l1,collapse="+"),
                                    modelset=var.l[i])
                results<-rbind(fixeddf,results)
                # posterior distributions
                for (h in 1:length(mod.inla$marginals.fixed)){
                    tmpmarg<-data.frame(mod.inla$marginals.fixed[[h]],
                                        var.name=names(mod.inla$marginals.fixed)[h],
                                        model=paste(var.l1,collapse="+"),modelset=var.l[i])
                    results.post<-rbind(tmpmarg,results.post)
                    
                }
                # residuals 
                I0 = inla.stack.index(indat, "fit")
                E0 = mod.inla$summary.fitted.values[I0$data,]
                resdf<-data.frame(res=indat$data$data$y-E0$mean,
                                  model=paste(var.l1,collapse="+"),
                                  modelset=var.l[i])
                results.res<-rbind(resdf,results.res)
            }
        }
        # Year model 
        # formula
        formfit<-as.formula(paste('y~-1+Intercept+Year+Year2',
                                  ranef,sep=""))
        # fit model
        mod.inla<-try(inla(formfit,family="gaussian",data=inla.stack.data(indat),
                           control.compute = list(dic=TRUE),
                           control.predictor=list(A=inla.stack.A(indat),compute=TRUE)))
        # summary of fixed effects
        fixeddf<-data.frame(mod.inla$summary.fixed,variable=row.names(mod.inla$summary.fixed),
                            DIC=mod.inla$dic$dic,model="Year2",
                            modelset="Year+Year2")
        results<-rbind(fixeddf,results)
        # posterior distributions
        for (h in 1:length(mod.inla$marginals.fixed)){
            tmpmarg<-data.frame(mod.inla$marginals.fixed[[h]],
                                var.name=names(mod.inla$marginals.fixed)[h],
                                model="Year",modelset="Year")
            results.post<-rbind(tmpmarg,results.post)
            
        }
        # residuals 
        I0 = inla.stack.index(indat, "fit")
        E0 = mod.inla$summary.fitted.values[I0$data,]
        resdf<-data.frame(res=indat$data$data$y-E0$mean,
                          model=paste(var.l1,collapse="+"),
                          modelset=var.l[i])
        results.res<-rbind(resdf,results.res)
    }
    ####### end of quadratic models #######    
    if(quad==FALSE){
        # store results
        results<-NULL
        results.post<-NULL
        results.res<-NULL
        # list of predictors
        var.l<-c("PA_cover","NDVI_1m","NDVI_3m","NDVI_12m","NDVI_24m","NDVI_36m")
        for (i in 1:length(var.l)){
            print(paste("variable of interest-->",var.l[i],sep=""))
            # create list of variables
            var.l1<-paste("Year","+",var.l[i],sep="")
            # formula
            formfit<-as.formula(paste('y~-1+Intercept+',paste(var.l1,collapse="+"),
                                      ranef,sep=""))
            # fit model
            mod.inla<-try(inla(formfit,family="gaussian",data=inla.stack.data(indat),
                               control.compute = list(dic=TRUE),
                               control.predictor=list(A=inla.stack.A(indat),compute=TRUE)))
            if(class(mod.inla)=="try-error"){
                next
            } else {
                # summary of fixed effects
                fixeddf<-data.frame(mod.inla$summary.fixed,variable=row.names(mod.inla$summary.fixed),
                                    DIC=mod.inla$dic$dic,model=var.l1,
                                    modelset=var.l[i])
                results<-rbind(fixeddf,results)
                # posterior distributions
                for (h in 1:length(mod.inla$marginals.fixed)){
                    tmpmarg<-data.frame(mod.inla$marginals.fixed[[h]],
                                        var.name=names(mod.inla$marginals.fixed)[h],
                                        model=var.l1,modelset=var.l[i])
                    results.post<-rbind(tmpmarg,results.post)
                    
                }
                # residuals 
                I0 = inla.stack.index(indat, "fit")
                E0 = mod.inla$summary.fitted.values[I0$data,]
                resdf<-data.frame(res=indat$data$data$y-E0$mean,
                                  model=var.l1,
                                  modelset=var.l[i])
                results.res<-rbind(resdf,results.res)
            }
        }
        
        # Year model 
        # formula
        formfit<-as.formula(paste('y~-1+Intercept+Year',
                                  ranef,sep=""))
        # fit model
        mod.inla<-try(inla(formfit,family="gaussian",data=inla.stack.data(indat),
                           control.compute = list(dic=TRUE),
                           control.predictor=list(A=inla.stack.A(indat),compute=TRUE)))
        # summary of fixed effects
        fixeddf<-data.frame(mod.inla$summary.fixed,variable=row.names(mod.inla$summary.fixed),
                            DIC=mod.inla$dic$dic,model="Year",
                            modelset="Year")
        results<-rbind(fixeddf,results)
        # posterior distributions
        for (h in 1:length(mod.inla$marginals.fixed)){
            tmpmarg<-data.frame(mod.inla$marginals.fixed[[h]],
                                var.name=names(mod.inla$marginals.fixed)[h],
                                model="Year",modelset="Year")
            results.post<-rbind(tmpmarg,results.post)
            
        }
    }
    ####### end of linear models #######    
    res.l<-vector("list",3)
    names(res.l)<-c("coefs","posteriors","residuals")
    row.names(results)<-1:dim(results)[1]
    res.l[[1]]<-results
    res.l[[2]]<- results.post
    res.l[[3]]<-results.res
    return(res.l)
} # end of function

    