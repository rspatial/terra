.geocolors <- function(name, ...) {
	if (name == "aspect") {
		#0% black, 50% white, 100% black
		pal <- grDevices::colorRampPalette(c("black", "white", "black"), ...)
		cols <- cbind(1:360, pal(360))
	} else if (name == "aspectclr") {
		pal <- grDevices::colorRampPalette(c("yellow", "green", "cyan", "red", "yellow"), ...)
		cols <- cbind(1:360, pal(360))
		cols[1,2] <- "#FFFFFF"
	}
}
