# Author: Robert J. Hijmans
# Date : October 2011
# October 2011
# version 1
# License GPL v3



setMethod("%in%", signature(x='SpatRaster', table='ANY'),
	function(x, table) {
		app(x, function(i) i %in% table)
	}
)


setMethod("match", signature(x='Raster', table='ANY', nomatch='ANY', incomparables='ANY'),
	function(x, table, nomatch, incomparables) {
		app(x, function(i) match(i, table, nomatch, incomparables))
	}
)

