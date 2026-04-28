# KML ExtendedData parsing (#1954) when OGR KML driver omits fields

if (!requireNamespace("XML", quietly = TRUE)) {
	exit_file("No XML package")
}

{ # one block for on.exit 

	kml <- c(
		"<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
		"<kml xmlns=\"http://www.opengis.net/kml/2.2\"><Document>",
		"<Placemark><name>PN</name>",
		"<ExtendedData><SchemaData>",
		"<SimpleData name=\"id\">7</SimpleData>",
		"<SimpleData name=\"camp_id\">c1</SimpleData>",
		"</SchemaData></ExtendedData>",
		"<Polygon><outerBoundaryIs><LinearRing>",
		"<coordinates>0,0,0 1,0,0 1,1,0 0,1,0 0,0,0</coordinates>",
		"</LinearRing></outerBoundaryIs></Polygon>",
		"</Placemark></Document></kml>"
	)
	f <- tempfile(fileext = ".kml")
	on.exit(unlink(f), add = TRUE)
	writeLines(kml, f, useBytes = TRUE)

	tab <- terra:::.read_kml_extended_table(f)
	expect_true(is.data.frame(tab))
	expect_equal(nrow(tab), 1L)
	expect_true("Name" %in% names(tab))
	expect_true(any(grepl("id", names(tab), ignore.case = TRUE)))

	# KMZ wrap (needs working zip in utils::zip)
	zf <- tempfile(fileext = ".kmz")
	on.exit(unlink(zf), add = TRUE)
	writeLines(kml, f, useBytes = TRUE)
	tr <- try(utils::zip(zf, f), silent = TRUE)
	if (!inherits(tr, "try-error") && file.exists(zf)) {
		ex <- terra:::.kml_extract_from_kmz(zf)
		expect_true(!is.null(ex))
		expect_true(file.exists(ex$path))
		tab2 <- terra:::.read_kml_extended_table(ex$path)
		expect_equal(nrow(tab2), 1L)
		unlink(ex$dir, recursive = TRUE)
	}
}
