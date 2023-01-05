// Copyright (c) 2018-2023  Robert J. Hijmans
//
// This file is part of the "spat" library.
//
// spat is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// spat is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with spat. If not, see <http://www.gnu.org/licenses/>.

#include "spatVector.h"
#include "string_utils.h"
#include <stdexcept>
#include "NA.h"


#ifdef useGDAL

#include "file_utils.h"
#include "ogrsf_frmts.h"



GDALDataset* SpatVector::write_ogr(std::string filename, std::string lyrname, std::string driver, bool append, bool overwrite, std::vector<std::string> options) {

    GDALDataset *poDS = NULL;
	if (filename != "") {
		if (file_exists(filename)) {
			if ((!overwrite) && (!append)) {
				setError("file exists. Use 'overwrite=TRUE' to overwrite it");
				return(poDS);
			} else {
				options.push_back("OVERWRITE=YES");
			}
		} else {
			append = false;
		}
		if (nrow() == 0) {
			setError("no geometries to write");
			return(poDS);
		}
	}

	if (append) {

		#if GDAL_VERSION_MAJOR < 3
			setError("GDAL >= 3 required for inserting layers into an existing file");
			return(poDS);
		#endif

		poDS = static_cast<GDALDataset*>(GDALOpenEx(filename.c_str(), GDAL_OF_VECTOR | GDAL_OF_UPDATE,
				NULL, NULL, NULL ));

		std::vector<std::string> lyrnms;

		size_t n = poDS->GetLayerCount();
		for (size_t i=0; i<n; i++) {
			OGRLayer *poLayer = poDS->GetLayer(i);
			if (poLayer != NULL) {
				lyrnms.push_back((std::string)poLayer->GetName());
			}
		}
		if (is_in_vector(lyrname, lyrnms)) {
			if (!overwrite) {
				setError("layer exists. Use 'overwrite=TRUE' to overwrite it");
				return(poDS);
			} else {
				options.push_back("OVERWRITE=YES");
			}
		}
	} else {
		GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName( driver.c_str() );
		if( poDriver == NULL )  {
			setError( driver + " driver not available");
			return poDS;
		}
		char **papszMetadata;
		papszMetadata = poDriver->GetMetadata();
		if (!CSLFetchBoolean( papszMetadata, GDAL_DCAP_VECTOR, FALSE)) {
			setError(driver + " is not a vector format");
			return poDS;
		}
		if (!CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE)) {
			setError("cannot create a "+ driver + " dataset");
			return poDS;
		}
		poDS = poDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL );
	}

    if( poDS == NULL ) {
        setError("Creation of output dataset failed" );
        return poDS;
    }

	OGRwkbGeometryType wkb;
	SpatGeomType geomtype = geoms[0].gtype;
	if (geomtype == points) {
		wkb = wkbPoint;
	} else if (geomtype == lines) {
		wkb = wkbMultiLineString;
	} else if (geomtype == polygons) {
		wkb = wkbMultiPolygon;
	} else {
        setError("this geometry type is not supported: " + type());
        return poDS;
	}

	std::string s = srs.wkt;


	OGRSpatialReference *SRS = NULL;
	if (s != "") {
		SRS = new OGRSpatialReference;
		OGRErr err = SRS->SetFromUserInput(s.c_str());
#if GDAL_VERSION_NUM >= 2050000
		SRS->SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
#endif
		if (err != OGRERR_NONE) {
			setError("crs error");
			delete SRS;
			return poDS;
		}
	}

	size_t nGroupTransactions = 0;

    OGRLayer *poLayer;
	char** papszOptions = NULL;
	if (options.size() > 0) {
		for (size_t i=0; i<options.size(); i++) {
			std::vector<std::string> gopt = strsplit(options[i], "=");
			if (gopt.size() == 2) {
				if (gopt[0] == "nGroupTransactions") {
					try  {
						nGroupTransactions = std::stoi(gopt[1]);
					} catch (std::invalid_argument &e)  {
						nGroupTransactions = 0;
					}
				} else {
					papszOptions = CSLSetNameValue(papszOptions, gopt[0].c_str(), gopt[1].c_str() );
				}
			}
		}
		// papszOptions = CSLSetNameValue( papszOptions, "ENCODING", "UTF-8" );
    }
	poLayer = poDS->CreateLayer(lyrname.c_str(), SRS, wkb, papszOptions);
	CSLDestroy(papszOptions);
    if( poLayer == NULL ) {
        setError( "Layer creation failed" );
        return poDS;
    }
//	if (SRS != NULL) SRS->Release();
	if (SRS != NULL) OSRDestroySpatialReference(SRS);


	std::vector<std::string> nms = get_names();
	std::vector<std::string> tps = df.get_datatypes();
	OGRFieldType otype;
	int nfields = nms.size();
	size_t ngeoms = size();


	for (int i=0; i<nfields; i++) {

		OGRFieldSubType eSubType = OFSTNone;
		if (tps[i] == "double") {
			otype = OFTReal;
		} else if (tps[i] == "long") {
			otype = OFTInteger64;
		} else if (tps[i] == "bool") {
			otype = OFTInteger;
			eSubType = OFSTBoolean;
		} else {
			otype = OFTString;
		}

		OGRFieldDefn oField(nms[i].c_str(), otype);
		oField.SetSubType(eSubType);
		if (otype == OFTString) {
			size_t w = 10;
			w = std::max(w, df.strwidth(i));
			oField.SetWidth(w); 
		}
		if( poLayer->CreateField( &oField ) != OGRERR_NONE ) {
			setError( "Field creation failed for: " + nms[i]);
			return poDS;
		}
	}

	// use a single transaction as in sf
	// makes a big difference for gpkg by avoiding many INSERTs
	bool can_do_transaction = poDS->TestCapability(ODsCTransactions); // == TRUE);
	bool transaction = false;
	if (can_do_transaction) {
		transaction = (poDS->StartTransaction() == OGRERR_NONE);
		if (! transaction) {
			setError("transaction failed");
			return poDS;
		}
	}
	// chunks

	if (nGroupTransactions == 0) {
		nGroupTransactions = 50000;
	}
	size_t gcntr = 0;
	long longNA = NA<long>::value;

	for (size_t i=0; i<ngeoms; i++) {

		OGRFeature *poFeature;
        poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
		for (int j=0; j<nfields; j++) {
			if (tps[j] == "double") {
				double dval = df.getDvalue(i, j);
				if (!std::isnan(dval)) {
					poFeature->SetField(j, df.getDvalue(i, j));
				}
			} else if (tps[j] == "long") {
				long ival = df.getIvalue(i, j);
				if (ival != longNA) {
					poFeature->SetField(j, (GIntBig)ival);
				}
			} else if (tps[j] == "bool") {
				poFeature->SetField(j, df.getBvalue(i, j));
			} else if (tps[j] == "time") {
				SpatTime_t tval = df.getTvalue(i, j);
				if (tval != longNA) {
					poFeature->SetField(j, (GIntBig)tval);
				}
			} else if (tps[j] == "factor") {
				SpatFactor f = df.getFvalue(i, j);
				if (f.v[0] != 0) {
					std::string s = f.getLabel(0);
					poFeature->SetField(j, f.getLabel(0).c_str());
				}
			} else {
				std::string s = df.getSvalue(i, j);
				if (s != df.NAS) {
					poFeature->SetField(j, df.getSvalue(i, j).c_str());
				}
			}
		}
		//r++;

// points -- also need to do multi-points
		OGRPoint pt;
		if (wkb == wkbPoint) {
			if (!std::isnan(geoms[i].parts[0].x[0])) {
				pt.setX( geoms[i].parts[0].x[0] );
				pt.setY( geoms[i].parts[0].y[0] );
			}
			poFeature->SetGeometry( &pt );

// lines
		} else if (wkb == wkbMultiLineString) {
			OGRMultiLineString poGeom;
			for (size_t j=0; j<geoms[i].size(); j++) {
				OGRLineString poLine = OGRLineString();
				for (size_t k=0; k<geoms[i].parts[j].size(); k++) {
					if (!std::isnan(geoms[i].parts[j].x[k])) {
						pt.setX(geoms[i].parts[j].x[k]);
						pt.setY(geoms[i].parts[j].y[k]);
						poLine.setPoint(k, &pt);
					}
				}
				if (poGeom.addGeometry(&poLine) != OGRERR_NONE ) {
					setError("cannot add line");
					return poDS;
				}
			}
			if (poFeature->SetGeometry( &poGeom ) != OGRERR_NONE) {
				setError("cannot set geometry");
				return poDS;
			}

// polygons
		} else if (wkb == wkbMultiPolygon) {
			SpatGeom g = getGeom(i);
			OGRMultiPolygon poGeom;
			for (size_t j=0; j<g.size(); j++) {
				OGRLinearRing poRing;
				SpatPart p = g.getPart(j);
				for (size_t k=0; k<p.size(); k++) {
					if (!std::isnan(p.x[k])) {
						pt.setX(p.x[k]);
						pt.setY(p.y[k]);
						poRing.setPoint(k, &pt);
					}
				}
				OGRPolygon polyGeom;
				if (polyGeom.addRing(&poRing) != OGRERR_NONE ) {
					setError("cannot add ring");
					return poDS;
				}

				if (p.hasHoles()) {
					for (size_t h=0; h < p.nHoles(); h++) {
						SpatHole hole = p.getHole(h);
						OGRLinearRing poHole;
						for (size_t k=0; k<hole.size(); k++) {
							pt.setX(hole.x[k]);
							pt.setY(hole.y[k]);
							poHole.setPoint(k, &pt);
						}
						if (polyGeom.addRing(&poHole) != OGRERR_NONE ) {
							setError("cannot add hole");
							return poDS;
						}
					}
				}
				poGeom.addGeometry( &polyGeom);
				//closeRings
			}

			//OGRMultiPolygon* mGeom = poGeom.toMultiPolygon();
			if (poFeature->SetGeometry( &poGeom ) != OGRERR_NONE) {
				setError("cannot set geometry");
				return poDS;
			}
		} else {
			setError("Only points, lines and polygons are currently supported");
			return poDS;
		}

		if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE ) {
			setError("Failed to create feature");
			return poDS;
        }
        OGRFeature::DestroyFeature( poFeature );
		gcntr++;
		if (transaction && (gcntr == nGroupTransactions)) {
			if (poDS->CommitTransaction() != OGRERR_NONE) {
				poDS->RollbackTransaction();
				setError("transaction commit failed");
			}
			gcntr = 0;
			transaction = (poDS->StartTransaction() == OGRERR_NONE);
			if (! transaction) {
				setError("transaction failed");
				return poDS;
			}
		}
    }
	if (transaction && (gcntr>0) && (poDS->CommitTransaction() != OGRERR_NONE)) {
		poDS->RollbackTransaction();
		setError("transaction commit failed");
	}
	return poDS;
}



bool SpatVector::write(std::string filename, std::string lyrname, std::string driver, bool append, bool overwrite, std::vector<std::string> options) {

	if (nrow() == 0) {
		addWarning("nothing to write");
		return false;
	}

	GDALDataset *poDS = write_ogr(filename, lyrname, driver, append, overwrite, options);
    if (poDS != NULL) GDALClose( poDS );
	if (hasError()) {
		return false;
	}
	return true;

}

GDALDataset* SpatVector::GDAL_ds() {
	return write_ogr("", "layer", "Memory", false, true, std::vector<std::string>());
}


#include <fstream>

bool SpatDataFrame::write_dbf(std::string filename, bool overwrite, SpatOptions &opt) {
// filename is here "raster.tif"
// to write "raster.tif.vat.dbf"

	if (filename != "") {
		if (file_exists(filename) & (!overwrite)) {
			setError("file exists. Use 'overwrite=TRUE' to overwrite it");
			return(false);
		}
		if (nrow() == 0) {
			setError("nothing to write");
			return(false);
		}
	}

	std::string fbase = tempFile(opt.get_tempdir(), opt.pid, "");
	std::string f = fbase + ".shp";
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName( "ESRI Shapefile" );
    GDALDataset *poDS = NULL;
    poDS = poDriver->Create(f.c_str(), 0, 0, 0, GDT_Unknown, NULL );
    if( poDS == NULL ) {
        setError("Creation of output dataset failed" );
        return false;
    }

	OGRwkbGeometryType wkb = wkbPoint;
	OGRSpatialReference *SRS = NULL;
    OGRLayer *poLayer;
    poLayer = poDS->CreateLayer("dbf", SRS, wkb, NULL );
    if( poLayer == NULL ) {
        setError( "Layer creation failed" );
        return false;
    }
	std::vector<std::string> nms = get_names();
	std::vector<std::string> tps = get_datatypes();
	OGRFieldType otype;
	int nfields = nms.size();

	for (int i=0; i<nfields; i++) {
		if (tps[i] == "double") {
			otype = OFTReal;
		} else if (tps[i] == "long") {
			otype = OFTInteger64;
		} else {
			otype = OFTString;
		}

		OGRFieldDefn oField(nms[i].c_str(), otype);
		if (otype == OFTString) {
			oField.SetWidth(50); // needs to be computed
		}
		if( poLayer->CreateField( &oField ) != OGRERR_NONE ) {
			setError( "Field creation failed for: " + nms[i]);
			return false;
		}
	}

	for (size_t i=0; i<nrow(); i++) {

		OGRFeature *poFeature;
        poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
		for (int j=0; j<nfields; j++) {
			if (tps[j] == "double") {
				poFeature->SetField(j, getDvalue(i, j));
			} else if (tps[j] == "long") {
				poFeature->SetField(j, (GIntBig) getIvalue(i, j));
			} else {
				poFeature->SetField(j, getSvalue(i, j).c_str());
			}
		}

		OGRPoint pt;
		pt.setX( 0.0 );
		pt.setY( 0.0 );
		poFeature->SetGeometry( &pt );

		if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE ) {
			setError("Failed to create feature");
			return false;
        }

        OGRFeature::DestroyFeature( poFeature );
    }
    GDALClose( poDS );
	f = fbase + ".dbf";
	filename += ".vat.dbf";
	// c++17 has file_copy
    std::ifstream  src(f.c_str(), std::ios::binary);
    std::ofstream  dst(filename.c_str(),  std::ios::binary);
    dst << src.rdbuf();

	filename.erase(filename.length()-3);
	filename += "cpg";
	std::ofstream cpg;
	cpg.open (filename.c_str());
	cpg << "UTF-8";
	cpg.close();

	return true;
}


bool SpatVector::delete_layers(std::string filename, std::vector<std::string> layers, bool return_error) {

	if (filename == "") {
		setError("empty filename");
		return false;
	}
	if (!file_exists(filename)) {
		setError("file does not exist");
		return false;
	}
	if (layers.size() == 0) return(true);

    GDALDataset *poDS = static_cast<GDALDataset*>(GDALOpenEx(filename.c_str(), GDAL_OF_VECTOR | GDAL_OF_UPDATE,
				NULL, NULL, NULL ));

    if( poDS == NULL ) {
        setError("Cannot open or update this dataset" );
        return false;
    }

	std::string fails;

	size_t n = poDS->GetLayerCount();
	for (int i =(n-1); i > 0; i--) {
		size_t m = layers.size();
		if (m == 0) break;

		OGRLayer *poLayer = poDS->GetLayer(i);
		if (poLayer == NULL) continue;
		std::string lname = poLayer->GetName();
		for (size_t j=0; j<m; j++) {
			if (lname.compare(layers[j]) == 0) {
				OGRErr err = poDS->DeleteLayer(i);
				if (err == OGRERR_UNSUPPORTED_OPERATION) {
					setError("Deleting layer not supported for this file (format / driver)");
					GDALClose(poDS);
					return(false);
				}
				if (err != OGRERR_NONE) {
					if (fails.size() > 0) {
						fails += ", " + layers[j];
					} else {
						fails = layers[j];
					}
				}
				layers.erase(layers.begin() + j);
				break;
			}
		}
	}
	GDALClose(poDS);
	if (layers.size() > 0) {
		fails += concatenate(layers, ", ");
	}
	if (fails.size() > 0) {
		if (return_error) {
			setError("deleting failed for: " + fails);
		} else {
			addWarning("deleting failed for: " + fails);
		}
	}
	return true;
}

#endif

