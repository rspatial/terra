
extentFromCells <- function (object, cells)
{
    cells <- stats::na.omit(unique(round(cells)))
    cells <- cells[cells > 0 & cells <= ncell(object)]
    if (length(cells) < 1) {
        stop("no valid cells")
    }
    r <- res(object)
    dx <- r[1] * c(-0.5, 0.5)
    dy <- r[2] * c(-0.5, 0.5)
    ext(range(xFromCell(object, cells)) + dx, range(yFromCell(object,
        cells)) + dy)
}

cellsFromExtent <- function (object, extent, expand = FALSE)
{
    object <- rast(object)
    extent <- align(ext(extent), object)
    innerBox <- intersect(ext(object), extent)
    if (is.null(innerBox)) {
        return(NULL)
    }
    srow <- rowFromY(object, ymax(innerBox) - 0.5 * yres(object))
    erow <- rowFromY(object, ymin(innerBox) + 0.5 * yres(object))
    scol <- colFromX(object, xmin(innerBox) + 0.5 * xres(object))
    ecol <- colFromX(object, xmax(innerBox) - 0.5 * xres(object))
    if (expand) {
        srow <- srow - round((ymax(extent) - ymax(innerBox))/yres(object))
        erow <- erow + round((ymin(innerBox) - ymin(extent))/yres(object))
        scol <- scol - round((xmin(innerBox) - xmin(extent))/xres(object))
        ecol <- ecol + round((xmax(extent) - xmax(innerBox))/xres(object))
    }
    return(cellFromRowColCombine(object, srow:erow, scol:ecol))
}



