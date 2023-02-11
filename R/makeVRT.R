# Author: Robert J. Hijmans
# Date :  October 2010, January 2022
# Version 1.0
# Licence GPL v3


gdalDType <- function(dtype) {
	dps <- c("INT2S", "INT4S", "INT1U", "INT2U", "INT4U", "FLT4S", "FLT8S")
	if (!(dtype %in% dps)) {
		stop(paste(dtype, "is not a valid data type. Should be one of:", paste(dps, collapse=", ")))
	}
	bytesize <- as.integer(substr(dtype, 4, 4))
	size <- bytesize * 8
	type <- substr(dtype, 1, 3)
	if (type == "INT") {
		type <- "Int"
		if (size == 64) {
			size <- 32
			warning("8 byte integer values not supported by GDAL, changed to 4 byte integer values")
		}
		if (substr(dtype, 5, 5) == "U") {
			if (size == 8) {
				return(c("Byte", 1))
			} else {
				type <- paste("U", type, sep="")
			}
		}
	} else {
		type <- "Float"
	}
	return(c(paste0(type, size), bytesize))
}


makeVRT <- function(filename, nrow, ncol, nlyr=1, extent, xmin, ymin, xres, yres=xres, xycenter=TRUE, crs="+proj=longlat", lyrnms="", datatype, NAflag=NA, bandorder="BIL", byteorder="LSB", toptobottom=TRUE, offset=0, scale=1) {

	stopifnot(length(filename)==1)
	stopifnot(file.exists(filename))
	if (tolower(tools::file_ext(filename)) == "vrt") {
		stop("cannot (over)write a vrt header for a vrt file")
	}
	lyrnms <- rep(lyrnms, length.out=nlyr)
	fvrt <- paste0(filename, ".vrt")

	if (missing(datatype)) {
		bytes <- file.info(filename)$size / (3601 * 3601)
		if (bytes == 1) {
			datatype <- "INT1U"
		} else if (bytes == 2) {
			datatype <- "INT2U"
		} else if (bytes == 4) {
			datatype <- "FLT4S"
		} else if (bytes == 8) {
			datatype <- "FLT8S"
		}
	}

	gd <- gdalDType(datatype[1])
	datatype <- gd[1]
	pixsize <- as.integer(gd[2])

	if (bandorder[1] == "BIL") {
		pixoff <- pixsize
		lineoff <- pixsize * ncol * nlyr
		imgoff <- ((1:nlyr)-1) * ncol * pixsize
	} else if (bandorder[1] == "BSQ") {
		pixoff <- pixsize
		lineoff <- pixsize * ncol
		imgoff <- ((1:nlyr)-1) *  nrow*ncol * pixsize
	} else if (bandorder[1] == "BIP") {
		pixoff <- pixsize * nlyr
		lineoff <- pixsize * ncol * nlyr
		imgoff <- (1:nlyr)-1
	} else {
		stop("unknown bandorder")
	}

	stopifnot(byteorder[1] %in% c("LSB", "MSB"))
	if (toptobottom[1]) { rotation <- 0 } else { rotation <- 180 }
	res <- abs(c(xres, yres))
	if (missing(extent)) {
		if (xycenter) {
			xmin <- xmin - res[1]/2
			ymin <- ymin - res[2]/2
		}
		ymax <- ymin + nrow * res[2]
	} else {
		xmin <- xmin(extent)
		ymax <- ymax(extent)
	}

	f <- file(fvrt, "w")
	cat('<VRTDataset rasterXSize="', ncol, '" rasterYSize="', nrow, '">\n' , sep = "", file = f)
	cat('<GeoTransform>', xmin, ', ', res[1], ', ', rotation, ', ', ymax, ', ', 0.0, ', ', -1*res[2], '</GeoTransform>\n', sep = "", file = f)

	if (! is.na(crs) ) {
		cat('<SRS>', crs ,'</SRS>\n', sep = "", file = f)
	}

	for (i in nlyr) {
		cat('\t<VRTRasterBand dataType="', datatype, '" band="', i, '" subClass="VRTRawRasterBand">\n', sep = "" , file = f)
		cat('\t\t<Description>', lyrnms[i], '</Description>\n', sep = "", file = f)
		cat('\t\t<SourceFilename relativetoVRT="1">', basename(filename), '</SourceFilename>\n', sep = "", file = f)
		cat('\t\t<ImageOffset>', imgoff[i], '</ImageOffset>\n', sep = "", file = f)
		cat('\t\t<PixelOffset>', pixoff, '</PixelOffset>\n', sep = "", file = f)
		cat('\t\t<LineOffset>', lineoff, '</LineOffset>\n', sep = "", file = f)
		cat('\t\t<ByteOrder>', byteorder, '</ByteOrder>\n', sep = "", file = f)
		if (!is.na(NAflag)) {
			cat('\t\t<NoDataValue>', NAflag, '</NoDataValue>\n', sep = "", file = f)
		}
		if (isTRUE(offset != 0) || isTRUE(scale != 1)) {
			cat('\t\t<Offset>', offset, '</Offset>\n', sep = "", file = f)
			cat('\t\t<Scale>', scale, '</Scale>\n', sep = "", file = f)
		}
		cat('\t</VRTRasterBand>\n', sep = "", file = f)
	}
	cat('</VRTDataset>\n', sep = "", file = f)
	close(f)
	return(fvrt)
}

# a = makeVRT(ff[1], 3601, 3601, 1, xmin=37, ymin=37, xres=1/3600, lyrnms="aspect", datatype="INT2U", byteorder="MSB")
