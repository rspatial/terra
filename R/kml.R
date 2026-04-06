# KML ExtendedData for OGR "KML" driver (e.g. Windows) vs "LIBKML" (many Linux builds).
# GDAL KML often only fills Name/Description; remaining fields sit in XML. See #1954.

.full_kml_filename <- function(x) {
	x <- enc2utf8(trimws(x[1]))
	if (grepl("\\.(kml|kmz)$", x, ignore.case = TRUE)) {
		x
	} else {
		NULL
	}
}


.vect_kml_source_for_xml <- function(x) {
	x <- enc2utf8(trimws(x[1]))
	if (grepl("^/vsicurl/", x)) {
		sub("^/vsicurl/", "", x)
	} else {
		x
	}
}


# local-name() XPaths work for all common KML xmlns layouts; the XML package then
# warns about default namespaces — suppressed inside .read_kml_extended_table.
.kml_has_geometry <- function(pm) {
	length(XML::getNodeSet(
		pm,
		paste(
			".//*[local-name()='Polygon' or local-name()='Point' or local-name()='LineString'",
			"or local-name()='MultiGeometry' or local-name()='Track'",
			"or local-name()='MultiTrack' or local-name()='Model']"
		)
	)) > 0
}


.parse_kml_placemark_attrs <- function(pm) {
	out <- list()
	getv <- function(nodes) {
		if (length(nodes)) XML::xmlValue(nodes[[1]]) else ""
	}
	nodes <- XML::getNodeSet(pm, "./*[local-name()='name']")
	if (length(nodes)) {
		out[["Name"]] <- getv(nodes)
	}
	nodes <- XML::getNodeSet(pm, "./*[local-name()='description']")
	if (length(nodes)) {
		out[["Description"]] <- getv(nodes)
	}
	for (sd in XML::getNodeSet(pm, ".//*[local-name()='SimpleData']")) {
		nm <- XML::xmlGetAttr(sd, "name")
		if (!is.na(nm) && nzchar(nm)) {
			out[[nm]] <- XML::xmlValue(sd)
		}
	}
	for (d in XML::getNodeSet(pm, ".//*[local-name()='Data']")) {
		nm <- XML::xmlGetAttr(d, "name")
		if (is.na(nm) || !nzchar(nm)) {
			next
		}
		vn <- XML::getNodeSet(d, "./*[local-name()='value']")
		out[[nm]] <- if (length(vn)) XML::xmlValue(vn[[1]]) else ""
	}
	out
}


.read_kml_extended_table <- function(filename) {
	suppressWarnings({
		doc <- XML::xmlParse(filename)
		pms <- XML::getNodeSet(doc, "//*[local-name()='Placemark']")
		if (length(pms) < 1) {
			return(NULL)
		}
		rows <- list()
		for (i in seq_along(pms)) {
			if (.kml_has_geometry(pms[[i]])) {
				rows[[length(rows) + 1L]] <- .parse_kml_placemark_attrs(pms[[i]])
			}
		}
		if (length(rows) < 1) {
			return(NULL)
		}
		alln <- unique(unlist(lapply(rows, names)))
		mat <- matrix(NA_character_, nrow = length(rows), ncol = length(alln))
		colnames(mat) <- alln
		for (i in seq_along(rows)) {
			r <- rows[[i]]
			for (nm in names(r)) {
				j <- match(nm, alln)
				if (!is.na(j)) {
					mat[i, j] <- r[[nm]]
				}
			}
		}
		out <- as.data.frame(mat, stringsAsFactors = FALSE)
	})
	out
}


# Unzip KMZ to a temp dir; return list(path=inner_kml, dir=tempdir) or NULL.
# dir must be unlinked after parsing. Prefer doc.kml (Google Earth convention).
.kml_extract_from_kmz <- function(kmz_path) {
	if (!file.exists(kmz_path)) {
		return(NULL)
	}
	zl <- try(utils::unzip(kmz_path, list = TRUE), silent = TRUE)
	if (inherits(zl, "try-error") || !is.data.frame(zl) || !length(zl$Name)) {
		return(NULL)
	}
	inzip <- zl$Name[grepl("\\.kml$", zl$Name, ignore.case = TRUE)]
	if (!length(inzip)) {
		return(NULL)
	}
	td <- tempfile("terra_kmz_")
	dir.create(td)
	st <- try(utils::unzip(kmz_path, exdir = td), silent = TRUE)
	if (inherits(st, "try-error")) {
		unlink(td, recursive = TRUE)
		return(NULL)
	}
	fls <- list.files(td, pattern = "\\.kml$", full.names = TRUE,
		recursive = TRUE, ignore.case = TRUE)
	if (!length(fls)) {
		unlink(td, recursive = TRUE)
		return(NULL)
	}
	docf <- fls[tolower(basename(fls)) == "doc.kml"]
	kmlf <- if (length(docf)) docf[[1]] else sort(fls)[1]
	list(path = kmlf, dir = td)
}


.try_kml_extended_attributes <- function(v, x, kml.extended = NULL) {
	if (!inherits(v, "SpatVector")) {
		return(v)
	}
	if (is.null(.full_kml_filename(x))) {
		return(v)
	}
	if (isFALSE(kml.extended)) {
		return(v)
	}
	# Always attempt merge when kml.extended is NULL or TRUE. LIBKML (common on
	# Linux) exposes many KML structure fields as columns (often NA); the plain
	# KML driver on Windows does not — merging from XML keeps results portable.
	if (!requireNamespace("XML", quietly = TRUE)) {
		if (isTRUE(kml.extended) || ncol(v) <= 2) {
			warn("vect", "install package 'XML' to use kml.extended; or set kml.extended=FALSE")
		}
		return(v)
	}

	fn <- .vect_kml_source_for_xml(x)
	tf <- NULL
	parse_path <- fn
	if (grepl("^https?://", fn, ignore.case = TRUE)) {
		url_base <- sub("[?#].*$", "", fn)
		ext <- if (grepl("\\.kmz$", url_base, ignore.case = TRUE)) {
			".kmz"
		} else {
			".kml"
		}
		tf <- tempfile(fileext = ext)
		on.exit(unlink(tf), add = TRUE)
		st <- try(utils::download.file(fn, tf, mode = "wb", quiet = TRUE), silent = TRUE)
		if (inherits(st, "try-error") || !isTRUE(st == 0L) || !file.exists(tf)) {
			warn("vect", "kml.extended: failed to download KML/KMZ for XML parsing")
			return(v)
		}
		parse_path <- tf
	} else {
		nx <- suppressWarnings(normalizePath(fn, winslash = "/", mustWork = FALSE))
		parse_path <- if (!is.na(nx) && nzchar(nx) && file.exists(nx)) nx else fn
		if (!file.exists(parse_path)) {
			return(v)
		}
	}

	kmz_work <- NULL
	if (grepl("\\.kmz$", parse_path, ignore.case = TRUE)) {
		kmz_work <- .kml_extract_from_kmz(parse_path)
		if (is.null(kmz_work)) {
			warn("vect", "kml.extended: could not extract KML from KMZ")
			return(v)
		}
		parse_path <- kmz_work$path
		on.exit(unlink(kmz_work$dir, recursive = TRUE), add = TRUE)
	}

	tab <- try(.read_kml_extended_table(parse_path), silent = TRUE)
	if (inherits(tab, "try-error") || is.null(tab)) {
		return(v)
	}
	if (nrow(tab) != nrow(v)) {
		warn("vect", paste0(
			"kml.extended: ", nrow(tab), " KML placemark(s) vs ", nrow(v),
			" GDAL feature(s); attributes not replaced"
		))
		return(v)
	}
	names(tab) <- make.names(names(tab), unique = TRUE)
	tab[] <- lapply(tab, utils::type.convert, as.is = FALSE)
	values(v) <- tab
	v
}


.vect_kml_merge_extended <- function(v, x, kml.extended, proxy, what) {
	if (isTRUE(proxy) || nzchar(what)) {
		return(v)
	}
	.try_kml_extended_attributes(v, x, kml.extended)
}
