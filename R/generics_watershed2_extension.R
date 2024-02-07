# Author: Emanuele Cordano
# Date : October 2023
# Version 1.0
# License GPL v3

setMethod("watershed2", signature(p="SpatRaster",pp_offset="integer"), 
          function(p,pp_offset,filename="", ...) { 
            ##  v <- match.arg(unique(v), c("aspect", "flowdir", "roughness", "slope", "TPI", "TRI"))
            ##  unit <- match.arg(unit, c("degrees", "radians"))
            ##  opt <- spatOptions(filename, ...)
            ##  seed <- ifelse("flowdirection" %in% v, .seed(), 0)
            print("watershed")
            opt <- spatOptions(filename, ...)
           ## p@ptr <- p@ptr$watershed2(as.integer(pp_offset-1),opt)
            p@ptr <- p@ptr$watershed2(as.integer(pp_offset-1),opt)
            messages(p, "watershed2") ## EC 20210318
            return(p)
            ## p@ptr <- uu
            ##messages(p, "watershed2")
          }
          
)

#### END EC 20210702

#### EC 20220809
setMethod("pitfinder2", signature(p="SpatRaster"), 
          function(p,filename="", ...) { 
            
            opt <- spatOptions(filename, ...)
           # p@ptr <- p@ptr$pitfinder2(opt)# EC 20231026
            p@ptr <- p@ptr$pitfinder2(opt)
            messages(p, "pitfinder2") ## EC 20210318
            return(p)
            ## p@ptr <- uu
            ##messages(p, "watershed2")
          }
          
)
#### END EC 20220809


#### EC 20231031
setMethod("NIDP2", signature(p="SpatRaster"), 
          function(p,filename="", ...) { 
            
            opt <- spatOptions(filename, ...)
         
            p@ptr <- p@ptr$NIDP2(opt)
            messages(p, "NIDP2") ## EC 20231031
            return(p)
         
          }
          
)
#### END EC 20231031

#### EC 20231104
setMethod("flowAccu2", signature(p="SpatRaster"), 
          function(p,filename="", ...) { 
            
            opt <- spatOptions(filename, ...)
            
            p@ptr <- p@ptr$flowAccu2(opt)
            messages(p, "flowAccu2") ## EC 20231104
            return(p)
            
          }
          
)
#### END EC 20231104


#### EC 2023114
setMethod("flowAccu2_weight", signature(p="SpatRaster",weight="SpatRaster"), 
          function(p,weight,filename="", ...) { 
            
            opt <- terra:::spatOptions(filename,...)
            print("ba")
            p@ptr <- p@ptr$flowAccu2_weight(weight@ptr,opt)
            messages(p, "flowAccu2_weight") ## EC 20231104
            return(p)
            
          }
          
)
#### END EC 20231104






