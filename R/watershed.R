# Author: Emanuele Cordano
# Date : October 2023
# Version 1.0
# License GPL v3

setMethod("watershed", signature(x="SpatRaster"), 
    function(x, pourpoint, filename="", ...) { 
        opt <- spatOptions(filename, ...)		
		cell <- cellFromXY(x, pourpoint)
		if (is.na(cell)) error("watershed", "pourpoint not on raster")
        x@pntr <- x@pntr$watershed2(as.integer(cell-1), opt)
        messages(x, "watershed") ## EC 20210318
    }
)

setMethod("pitfinder", signature(x="SpatRaster"), 
    function(x,pits_on_boundary=TRUE,filename="", ...) { 
        opt <- spatOptions(filename, ...)
        x@pntr <- x@pntr$pitfinder2(as.integer(pits_on_boundary),opt)

        messages(x, "pitfinder") ## EC 20210318
    }
)

setMethod("NIDP", signature(x="SpatRaster"), 
    function(x, filename="", ...) { 
        opt <- spatOptions(filename, ...)
        x@pntr <- x@pntr$NIDP2(opt)
        messages(x, "NIDP") ## EC 20231031
    }
)

setMethod("flowAccumulation", signature(x="SpatRaster"), 
    function(x, weight=NULL, filename="", ...) { 
		opt <- spatOptions(filename, ...)
		if (is.null(weight)) {      
			x@pntr <- x@pntr$flowAccu2(opt)
	    } else {
			x@pntr <- x@pntr$flowAccu2_weight(weight@pntr, opt)
		} 
		messages(x, "flowAccumulation") 
    }      
)

setMethod("flowDir", signature(x="SpatRaster"), 
	function(x, lambda=0.5, deviation_type=c("ltd","lad"), max_iters=10^6, filename="", ...) { 
		## http://www.idrologia.unimore.it/orlandini/web-archive/seminars/nyc-2008-2.pdf
		## ltd least transverse deviation
		## lad least angular deviation
<<<<<<< HEAD
	  if (length(deviation_type)<1) deviation_type <- "ltd"
	  if (grepl(" ", deviation_type[1])) {
	    deviation_type <- deviation_type[1] |>
	      trimws() |>              # remove leading/trailing whitespace (base R)
	      strsplit("\\s+") |>      # split into words using regex (base R)
	      unlist() |>              # flatten list into a character vector
	      substr(1, 1) |>          # extract first character of each word
	      paste(collapse = "")     # join initials into a single string
	  } else {
	    deviation_type <- deviation_type[1]
	  }
		deviation_type <- match.arg(tolower(deviation_type), c("ltd", "lad"))
		
		use_lad <- deviation_type=="lad"
=======
		deviation_type <- match.arg(tolower(deviation_type), c("ltd", "lad"))
		use_lad <- deviation_type == "lad"
>>>>>>> 4fb7fade7d313e14147261b00d649a0219e4397d
		opt <- spatOptions(filename, ...)
		x@pntr <- x@pntr$d8ltd(lambda, use_lad, max_iters, opt)
		names(x) <- sprintf("flowdir_%s_l=%s", deviation_type[1], as.character(lambda))
		messages(x, "flowDir") ## EC 20210318
	}
)


setMethod("pitfiller", signature(x="SpatRaster"), 
	function(x,pit=NULL,flowdir=NULL,niter=10,lambda=0,deviation_type="lad",max_iters=10^6,U=1,D=300,beta=0.9,theta_exp=0.5,filename="",...) { 
	 
	  
	  deviation_type <- match.arg(tolower(deviation_type), c("ltd", "lad"))
	  
		if (is.null(flowdir)) flowdir <- terrain(x,"flowdir") 
		if (is.null(pitfinder)) pit <- pitfinder(flowdir) 
		use_lad=1
		if (deviation_type=="ltd") use_lad=0
		flowdir[pit>0] <- 0 
		opt <- spatOptions(filename, ...)
		##  uselad=0
		print(pit)
		x@pntr <- x@pntr$pitfillerm(pit@pntr,flowdir@pntr,niter,lambda,use_lad,max_iters,U,D,beta,theta_exp,opt)
		messages(x, "pitfiller") ## EC 20210318
	}
)


# SpatRaster  SpatRaster::pitfillerm(SpatRaster pits,SpatRaster flowdirs,int niter, double lambda,int use_lad,
#				    double U,double D,double beta,double theta_exp, // see // see reference doi:10.1016/j.advwatres.2006.11.016)    
# SpatOptions &opt) {
