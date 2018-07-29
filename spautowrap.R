#' Wrapper function for calculating Moran's correlograms
#' from dataframe containing residuals
#' @indat= input dataframe
#' @group= grouping variable (i.e. model)
#' @coor= two-column dataframe containing coordinates (easting and northing)

spautowrap<-function(indat,group,coor){
    mod.l<-as.character(unique(indat[,group]))
    results<-NULL
    for (mod.i in mod.l){
        # subset residuals of interest
        tmp<-indat[indat[,group] %in% mod.i,]
        # compute spatial correlogram
        mod.c<-correlog(x=coor[,1], y=coor[,2],
                        z=tmp$res,resamp=0,increment=10)
        # wrap results into dataframe
        resdf<-data.frame(distance=mod.c$mean.of.class,correlation=mod.c$correlation,
                          model=unique(tmp$model),modelset=unique(tmp$modelset))
        results<-rbind(resdf,results)
        
    }
    return(results)
}    
        