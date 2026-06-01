// Copyright (c) 2018-2026  Robert J. Hijmans
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

// Round-tripping a SpatNetwork through GDAL's Geographic Network Model
// (GNM). See https://gdal.org/en/stable/user/gnm_data_model.html.
//
// On disk a network is a directory (driver "GNMFile") or a database
// (driver "GNMDatabase") containing two kinds of layers:
//   - "class" layers (here: a `nodes` point layer and an `edges` line
//     layer) that hold the geometry,
//   - "system" layers (`graph`, `meta`, `features`) that hold the
//     incidence list, network metadata and the GFID lookup. We let GNM
//     manage the system layers; we only fill the class layers and call
//     ConnectFeatures() to register edges.

#include "spatNetwork.h"

// Detect whether the GDAL build that we are linking against ships the
// GNM headers AND is recent enough for the API we rely on. GNM has
// been part of GDAL since 2.0, but a few of the helpers stabilized later.

#include "gdal_version.h"
#if defined(__has_include)
#  if __has_include(<gnm.h>) && defined(GDAL_VERSION_NUM) \
      && GDAL_VERSION_NUM >= GDAL_COMPUTE_VERSION(3, 0, 0)
#    define TERRA_HAS_GNM 1
#  endif
#endif

#if TERRA_HAS_GNM
#include "ogrsf_frmts.h"
#include "gdal_priv.h"
#include "ogr_spatialref.h"
#include "cpl_conv.h"
#include "cpl_string.h"
#include "cpl_vsi.h"
#include <gnm.h>
#include <gnmgraph.h>
#endif

#include <map>
#include <utility>


// Split `path` into (parent_dir, basename). `path` may use forward or
// backslashes; the result always uses the input's slash style. If `path`
// has no separator (`./mynet` or just `mynet`), parent becomes ".".
#if TERRA_HAS_GNM
static void split_parent_basename(const std::string &path, std::string &parent, std::string &base) {
	std::string::size_type p = path.find_last_of("/\\");
	if (p == std::string::npos) {
		parent = ".";
		base = path;
	} else if (p == 0) {
		// "/foo" -> parent "/" base "foo"
		parent = path.substr(0, 1);
		base   = path.substr(1);
	} else {
		parent = path.substr(0, p);
		base   = path.substr(p + 1);
	}
	// Strip a trailing slash on `parent` for consistency, except when
	// `parent` is just "/".
	if (parent.size() > 1) {
		char c = parent[parent.size() - 1];
		if (c == '/' || c == '\\') parent.resize(parent.size() - 1);
	}
}
#endif


bool SpatNetwork::write_gnm(std::string filename,
                            std::string driver_name,
                            std::vector<std::string> options) {
#if !TERRA_HAS_GNM
	(void)filename; (void)driver_name; (void)options;
	setError("terra needs GDAL 3 for GNM support");
	return false;
#else
	GDALAllRegister();
	if (driver_name.empty()) driver_name = "GNMFile";
	GDALDriver *poDriver =
		GetGDALDriverManager()->GetDriverByName(driver_name.c_str());
	if (!poDriver) {
		setError("GDAL driver '" + driver_name + "' is not available");
		return false;
	}

	// GNM's FormPath joins (pszFilename, GNM_MD_NAME) into the actual
	// network directory. So we pass the *parent* of the user-supplied
	// path as pszFilename, and use the basename as the network name.
	// The parent must exist; create it if it doesn't.
	std::string parent_dir, net_name;
	split_parent_basename(filename, parent_dir, net_name);
	if (net_name.empty()) {
		setError("could not derive a network name from '" + filename + "'");
		return false;
	}
	if (parent_dir != "." && parent_dir != "/") {
		VSIStatBufL sStat;
		if (VSIStatL(parent_dir.c_str(), &sStat) != 0) {
			if (VSIMkdirRecursive(parent_dir.c_str(), 0755) != 0) {
				setError("could not create parent directory '" + parent_dir + "'");
				return false;
			}
		}
	}

	char **papszOptions = nullptr;
	papszOptions = CSLAddNameValue(papszOptions, GNM_MD_NAME, net_name.c_str());
	papszOptions = CSLAddNameValue(papszOptions, "OVERWRITE", "YES");
	// GNM requires an SRS at create time. If we don't have one, hand it
	// a placeholder that GDAL will accept as "unknown".
	std::string srs_str = srs.wkt;
	if (srs_str.empty()) srs_str = "LOCAL_CS[\"unknown\"]";
	papszOptions = CSLAddNameValue(papszOptions, GNM_MD_SRS, srs_str.c_str());
	for (size_t i = 0; i < options.size(); i++) {
		const std::string &kv = options[i];
		std::string::size_type p = kv.find('=');
		if (p != std::string::npos) {
			std::string k = kv.substr(0, p);
			std::string v = kv.substr(p + 1);
			papszOptions = CSLAddNameValue(papszOptions, k.c_str(), v.c_str());
		}
	}

	GDALDataset *poDS = poDriver->Create(parent_dir.c_str(), 0, 0, 0,
	                                     GDT_Unknown, papszOptions);
	CSLDestroy(papszOptions);
	if (!poDS) {
		setError("could not create GNM network at '" + filename + "'");
		return false;
	}
	GNMGenericNetwork *poNet = dynamic_cast<GNMGenericNetwork *>(poDS);
	if (!poNet) {
		GDALClose(poDS);
		setError("created dataset is not a GNMGenericNetwork");
		return false;
	}

	// GNMFileNetwork::Create does NOT reliably write the `_gnm_srs.prj`
	// sidecar even though its Open path needs it (Open calls
	// LoadNetworkSrs, which reads that file). Write it ourselves so
	// that the dataset can be read back -- by terra, by ogrinfo, or by
	// any other GDAL consumer. We just persist the WKT we already
	// have: `srs.wkt` is the terra-side WKT2 string we used in
	// GNM_MD_SRS above, so going through OGRSpatialReference here
	// would only round-trip the same text through GDAL.
	{
		std::string srs_path = filename + "/_gnm_srs.prj";
		VSILFILE *fp = VSIFOpenL(srs_path.c_str(), "w");
		if (fp) {
			if (!srs_str.empty()) {
				VSIFWriteL(srs_str.data(), srs_str.size(), 1, fp);
			}
			VSIFCloseL(fp);
		}
	}

	OGRSpatialReference oSRS;
	OGRSpatialReference *pSRS = nullptr;
	if (!srs.wkt.empty()) {
		if (oSRS.SetFromUserInput(srs.wkt.c_str()) == OGRERR_NONE) {
			pSRS = &oSRS;
		}
	}

	// --- nodes layer ---------------------------------------------------
	OGRLayer *poNodes = poNet->CreateLayer("nodes", pSRS, wkbPoint, nullptr);
	if (!poNodes) {
		GDALClose(poDS);
		setError("could not create 'nodes' layer in GNM dataset");
		return false;
	}
	std::vector<GNMGFID> node_gfid(node_x.size(), -1);
	int gnm_fid_idx_n = poNodes->GetLayerDefn()->GetFieldIndex("gnm_fid");
	for (size_t i = 0; i < node_x.size(); i++) {
		OGRFeature *f = OGRFeature::CreateFeature(poNodes->GetLayerDefn());
		OGRPoint pt(node_x[i], node_y[i]);
		f->SetGeometry(&pt);
		if (poNodes->CreateFeature(f) == OGRERR_NONE) {
			if (gnm_fid_idx_n >= 0) node_gfid[i] = f->GetFieldAsInteger64(gnm_fid_idx_n);
			else node_gfid[i] = f->GetFID();
		}
		OGRFeature::DestroyFeature(f);
	}

	// --- edges layer ---------------------------------------------------
	OGRLayer *poEdges = poNet->CreateLayer("edges", pSRS, wkbLineString, nullptr);
	if (!poEdges) {
		GDALClose(poDS);
		setError("could not create 'edges' layer in GNM dataset");
		return false;
	}
	{
		OGRFieldDefn fd("length", OFTReal);
		OGRErr e = poEdges->CreateField(&fd);
		(void) e;
	}
	if (weighted) {
		OGRFieldDefn fd("weight", OFTReal);
		OGRErr e = poEdges->CreateField(&fd);
		(void) e;
	}
	int gnm_fid_idx_e = poEdges->GetLayerDefn()->GetFieldIndex("gnm_fid");
	int len_idx_w = poEdges->GetLayerDefn()->GetFieldIndex("length");
	int w_idx_w   = poEdges->GetLayerDefn()->GetFieldIndex("weight");

	std::vector<GNMGFID> edge_gfid(edge_from.size(), -1);
	for (size_t i = 0; i < edge_from.size(); i++) {
		OGRFeature *f = OGRFeature::CreateFeature(poEdges->GetLayerDefn());
		OGRLineString ls;
		for (size_t j = 0; j < edge_x[i].size(); j++) {
			ls.addPoint(edge_x[i][j], edge_y[i][j]);
		}
		f->SetGeometry(&ls);
		if (len_idx_w >= 0 && i < edge_length.size()) {
			f->SetField(len_idx_w, edge_length[i]);
		}
		if (weighted && w_idx_w >= 0 && i < edge_weight.size()) {
			f->SetField(w_idx_w, edge_weight[i]);
		}
		if (poEdges->CreateFeature(f) == OGRERR_NONE) {
			if (gnm_fid_idx_e >= 0) edge_gfid[i] = f->GetFieldAsInteger64(gnm_fid_idx_e);
			else edge_gfid[i] = f->GetFID();
		}
		OGRFeature::DestroyFeature(f);
	}

	// --- connections (incidence list) ----------------------------------
	GNMDirection dir_const = directed ? GNM_EDGE_DIR_SRCTOTGT
	                                  : GNM_EDGE_DIR_BOTH;
	for (size_t i = 0; i < edge_from.size(); i++) {
		if (edge_gfid[i] < 0) continue;
		size_t fi = edge_from[i], ti = edge_to[i];
		if (fi >= node_gfid.size() || ti >= node_gfid.size()) continue;
		GNMGFID f = node_gfid[fi];
		GNMGFID t = node_gfid[ti];
		if (f < 0 || t < 0) continue;
		double cost = (weighted && i < edge_weight.size()) ? edge_weight[i]
		            : (i < edge_length.size() ? edge_length[i] : 1.0);
		double invc = cost;     // GNM uses the same unit either way
		poNet->ConnectFeatures(f, t, edge_gfid[i], cost, invc, dir_const);
	}

	GDALClose(poDS);
	return true;
#endif
}


bool SpatNetwork::read_gnm(std::string filename) {
#if !TERRA_HAS_GNM
	(void)filename;
	setError("terra needs GDAL 3 for GNM support");
	return false;
#else
	GDALAllRegister();
	GDALDataset *poDS = (GDALDataset *) GDALOpenEx(filename.c_str(),
		GDAL_OF_GNM | GDAL_OF_READONLY,
		nullptr, nullptr, nullptr);
	if (!poDS) {
		// Some GDAL builds open GNM only when the GNM open option is the
		// only mode bit; others require GDAL_OF_VECTOR alongside. Try
		// the alternate combination before giving up.
		poDS = (GDALDataset *) GDALOpenEx(filename.c_str(),
			GDAL_OF_GNM | GDAL_OF_VECTOR | GDAL_OF_READONLY,
			nullptr, nullptr, nullptr);
	}
	if (!poDS) {
		setError("could not open GNM dataset at " + filename);
		return false;
	}
	GNMGenericNetwork *poNet = dynamic_cast<GNMGenericNetwork *>(poDS);
	if (!poNet) {
		// Some GDAL builds with mismatched RTTI between the link unit
		// and the loaded driver fail dynamic_cast even on valid GNM
		// datasets. Fall back to checking the driver name and casting
		// statically.
		GDALDriver *drv = poDS->GetDriver();
		const char *drv_name = drv ? drv->GetDescription() : "";
		if (drv_name && (EQUAL(drv_name, "GNMFile") || EQUAL(drv_name, "GNMDatabase"))) {
			poNet = static_cast<GNMGenericNetwork *>(poDS);
		}
	}
	if (!poNet) {
		GDALClose(poDS);
		setError("dataset at " + filename + " is not a GNMGenericNetwork");
		return false;
	}

	// Reset *this.
	node_x.clear();   node_y.clear();
	edge_from.clear(); edge_to.clear();
	edge_x.clear();   edge_y.clear();
	edge_source.clear();
	edge_length.clear();
	edge_weight.clear();
	weighted = false;
	directed = false;
	edge_df = SpatDataFrame();
	extent = SpatExtent();

	const OGRSpatialReference *poSRS = poNet->GetSpatialRef();
	if (poSRS) {
		char *wkt = nullptr;
		poSRS->exportToWkt(&wkt);
		if (wkt && *wkt) {
			std::string m;
			srs.set(std::string(wkt), m);
			CPLFree(wkt);
		} else if (wkt) {
			CPLFree(wkt);
		}
	}
	// Fallback: if the GNM-managed SRS came back empty, read the
	// `_gnm_srs.prj` sidecar directly. Stock GDAL builds occasionally
	// fail to populate m_oSRS during Open even though the file is
	// present.
	if (srs.wkt.empty()) {
		std::string srs_path = filename + "/_gnm_srs.prj";
		VSILFILE *fp = VSIFOpenL(srs_path.c_str(), "r");
		if (fp) {
			std::string buf;
			char chunk[4096];
			size_t n;
			while ((n = VSIFReadL(chunk, 1, sizeof(chunk), fp)) > 0) {
				buf.append(chunk, n);
			}
			VSIFCloseL(fp);
			// Strip trailing whitespace.
			while (!buf.empty() &&
			       (buf.back() == '\n' || buf.back() == '\r' ||
			        buf.back() == ' '  || buf.back() == '\t')) {
				buf.pop_back();
			}
			if (!buf.empty()) {
				std::string m;
				srs.set(buf, m);
			}
		}
	}

	// --- locate the class layers ---------------------------------------
	OGRLayer *poNodes = poNet->GetLayerByName("nodes");
	OGRLayer *poEdges = poNet->GetLayerByName("edges");
	if (!poNodes || !poEdges) {
		// Fallback: take the first non-system point layer and the first
		// non-system line layer.
		for (int li = 0; li < poNet->GetLayerCount(); li++) {
			OGRLayer *L = poNet->GetLayer(li);
			if (!L) continue;
			const char *nm = L->GetName();
			if (!nm || nm[0] == '_') continue;     // skip system layers
			OGRwkbGeometryType gt = wkbFlatten(L->GetGeomType());
			if (!poNodes && gt == wkbPoint)      poNodes = L;
			else if (!poEdges && gt == wkbLineString) poEdges = L;
		}
	}
	if (!poNodes || !poEdges) {
		GDALClose(poDS);
		setError("could not find a node (point) and an edge (line) layer "
		         "in the GNM dataset");
		return false;
	}

	// --- read nodes ----------------------------------------------------
	std::map<GNMGFID, size_t> gfid_to_node;
	{
		int idx = poNodes->GetLayerDefn()->GetFieldIndex("gnm_fid");
		poNodes->ResetReading();
		OGRFeature *f;
		while ((f = poNodes->GetNextFeature()) != nullptr) {
			GNMGFID gfid = (idx >= 0) ? f->GetFieldAsInteger64(idx)
			                          : f->GetFID();
			OGRGeometry *g = f->GetGeometryRef();
			if (g && wkbFlatten(g->getGeometryType()) == wkbPoint) {
				OGRPoint *pt = (OGRPoint *) g;
				gfid_to_node[gfid] = node_x.size();
				node_x.push_back(pt->getX());
				node_y.push_back(pt->getY());
			}
			OGRFeature::DestroyFeature(f);
		}
	}

	// --- read edges (geometry + cached length / weight) ---------------
	struct EdgeRec {
		std::vector<double> xs, ys;
		double length = 0.0;
		double weight = 0.0;
		bool has_weight = false;
		bool consumed = false;
	};
	std::map<GNMGFID, EdgeRec> gfid_to_edge;
	{
		int idx_g = poEdges->GetLayerDefn()->GetFieldIndex("gnm_fid");
		int idx_l = poEdges->GetLayerDefn()->GetFieldIndex("length");
		int idx_w = poEdges->GetLayerDefn()->GetFieldIndex("weight");
		poEdges->ResetReading();
		OGRFeature *f;
		while ((f = poEdges->GetNextFeature()) != nullptr) {
			GNMGFID gfid = (idx_g >= 0) ? f->GetFieldAsInteger64(idx_g)
			                            : f->GetFID();
			OGRGeometry *g = f->GetGeometryRef();
			if (g && wkbFlatten(g->getGeometryType()) == wkbLineString) {
				OGRLineString *ls = (OGRLineString *) g;
				EdgeRec er;
				er.xs.reserve(ls->getNumPoints());
				er.ys.reserve(ls->getNumPoints());
				for (int j = 0; j < ls->getNumPoints(); j++) {
					er.xs.push_back(ls->getX(j));
					er.ys.push_back(ls->getY(j));
				}
				if (idx_l >= 0) er.length = f->GetFieldAsDouble(idx_l);
				if (idx_w >= 0) {
					er.has_weight = true;
					er.weight = f->GetFieldAsDouble(idx_w);
				}
				gfid_to_edge[gfid] = std::move(er);
			}
			OGRFeature::DestroyFeature(f);
		}
	}

	// --- read the system "graph" layer (incidence list) ---------------
	// GNM does NOT expose the graph layer through GetLayerByName --
	// it loads the connections internally into a GNMGraph and the layer
	// itself is hidden from the public dataset API. So we open the
	// underlying file (`_gnm_graph.<ext>`) ourselves; the file naming
	// is part of GNM's public on-disk format.
	OGRLayer *poGraph = nullptr;
	GDALDataset *poGraphDS = nullptr;
	{
		// Discover the storage extension by listing the network dir.
		std::string graph_path;
		char **papszFiles = VSIReadDir(filename.c_str());
		if (papszFiles) {
			for (int i = 0; papszFiles[i] != nullptr; i++) {
				// CPLGetBasename strips the directory and extension from
				// a path. CPLGetBasenameSafe (GDAL >= 3.10) returns a
				// std::string-like wrapper, but CPLGetBasename works
				// across our entire supported GDAL range and is fine
				// for a one-shot equality check.
				const char *bn = CPLGetBasename(papszFiles[i]);
				if (bn && EQUAL(bn, "_gnm_graph")) {
					graph_path = filename + "/" + papszFiles[i];
					break;
				}
			}
			CSLDestroy(papszFiles);
		}
		if (!graph_path.empty()) {
			poGraphDS = (GDALDataset *) GDALOpenEx(graph_path.c_str(),
				GDAL_OF_VECTOR | GDAL_OF_READONLY,
				nullptr, nullptr, nullptr);
			if (poGraphDS && poGraphDS->GetLayerCount() > 0) {
				poGraph = poGraphDS->GetLayer(0);
			}
		}
	}
	if (!poGraph) {
		if (poGraphDS) GDALClose(poGraphDS);
		GDALClose(poDS);
		setError("the GNM dataset does not expose its graph (connections) layer");
		return false;
	}
	OGRFeatureDefn *gd = poGraph->GetLayerDefn();
	int s_idx = gd->GetFieldIndex("source");
	int t_idx = gd->GetFieldIndex("target");
	int c_idx = gd->GetFieldIndex("connector");
	int dir_idx = gd->GetFieldIndex("direction");
	if (s_idx < 0 || t_idx < 0 || c_idx < 0) {
		GDALClose(poDS);
		setError("graph layer is missing source/target/connector fields");
		return false;
	}

	bool any_directed = false, any_undirected = false;
	poGraph->ResetReading();
	OGRFeature *f;
	while ((f = poGraph->GetNextFeature()) != nullptr) {
		GNMGFID s = f->GetFieldAsInteger64(s_idx);
		GNMGFID t = f->GetFieldAsInteger64(t_idx);
		GNMGFID c = f->GetFieldAsInteger64(c_idx);
		int dir_v = (dir_idx >= 0) ? f->GetFieldAsInteger(dir_idx) : 0;
		OGRFeature::DestroyFeature(f);

		std::map<GNMGFID, size_t>::iterator its = gfid_to_node.find(s);
		std::map<GNMGFID, size_t>::iterator itt = gfid_to_node.find(t);
		std::map<GNMGFID, EdgeRec>::iterator ite = gfid_to_edge.find(c);
		if (its == gfid_to_node.end() || itt == gfid_to_node.end()
		    || ite == gfid_to_edge.end()) continue;

		// GNM may emit one row per direction for a bidirectional edge.
		// Each edge GFID is unique -- dedupe by setting `consumed`.
		if (ite->second.consumed) continue;
		ite->second.consumed = true;

		edge_from.push_back(its->second);
		edge_to.push_back(itt->second);
		edge_x.push_back(std::move(ite->second.xs));
		edge_y.push_back(std::move(ite->second.ys));
		edge_source.push_back(-1);
		edge_length.push_back(ite->second.length);
		if (ite->second.has_weight) {
			edge_weight.push_back(ite->second.weight);
			weighted = true;
		}

		if (dir_v == GNM_EDGE_DIR_BOTH) any_undirected = true;
		else                            any_directed   = true;
	}
	if (poGraphDS) GDALClose(poGraphDS);
	GDALClose(poDS);

	directed = any_directed && !any_undirected;

	if (weighted && edge_weight.size() != edge_from.size()) {
		// Inconsistent: not every edge had a `weight` value.
		edge_weight.clear();
		weighted = false;
	}

	// If lengths weren't stored, recompute them from the geometry so
	// callers always get sensible values.
	bool any_length = false;
	for (size_t i = 0; i < edge_length.size(); i++) {
		if (edge_length[i] > 0.0) { any_length = true; break; }
	}
	if (!any_length) compute_edge_lengths();

	computeExtent();
	return true;
#endif
}
