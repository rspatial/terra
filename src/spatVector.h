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

#ifndef SPATVECTOR_GUARD
#define SPATVECTOR_GUARD

#include "spatDataframe.h"

#ifdef useGDAL
#include "gdal_priv.h"
#endif


enum SpatGeomType { points, lines, polygons, null};


class SpatHole {
	public:
		virtual ~SpatHole(){}

		std::vector<double> x, y;
		SpatExtent extent;
		//constructors
		SpatHole();
		SpatHole(std::vector<double> X, std::vector<double> Y);
		//methods
		size_t size() { return x.size(); }
		bool empty() { return x.empty(); }
};

class SpatPart {
	public:
		virtual ~SpatPart(){}
		std::vector<double> x, y; //, z;
		std::vector< SpatHole > holes; // polygons only
		SpatExtent extent;

		//constructors
		SpatPart();
		SpatPart(std::vector<double> X, std::vector<double> Y);
		SpatPart(double X, double Y);

		//methods
		size_t size() { return x.size(); }
		bool empty() { return x.empty(); }
		//holes, polygons only
		bool addHole(std::vector<double> X, std::vector<double> Y);
		bool addHole(SpatHole h);
		SpatHole getHole(unsigned i) { return( holes[i] ) ; }
		bool hasHoles() { return !holes.empty();}
		unsigned nHoles() { return holes.size();}
		size_t ncoords();
		bool is_CCW();
};


class SpatGeom {
	public:
		//constructors
		SpatGeom();
		SpatGeom(SpatGeomType g);
		SpatGeom(SpatPart p, SpatGeomType type);
		virtual ~SpatGeom(){}

		SpatGeomType gtype = null;
		std::vector<SpatPart> parts;
		SpatExtent extent;

		//methods
		bool unite(SpatGeom g);
		bool addPart(SpatPart p);
		bool addHole(SpatHole h);
		bool setPart(SpatPart p, unsigned i);
		bool reSetPart(SpatPart p);
		SpatPart getPart(unsigned i);
		//double area_plane();
		//double area_lonlat(double a, double f);
		//double length_plane();
		//double length_lonlat(double a, double f);
		unsigned size() { return parts.size(); };
		bool empty() { return parts.empty(); };
		void remove_duplicate_nodes(int digits);
		size_t ncoords();
		std::vector<std::vector<double>> coordinates();
		void computeExtent();
		void reserve(size_t n) {
			parts.reserve(n);
		}
};


class SpatVectorCollection;

class SpatVector {

	public:
		std::vector<SpatGeom> geoms;
		SpatExtent extent;
		SpatDataFrame df;
		//std::vector<std::string> crs;
		SpatSRS srs;
		bool is_proxy = false;
		std::string read_query = "";
		std::vector<double> read_extent;
		std::string source = "";
		std::string source_layer = "";
		size_t geom_count = 0;
		
		SpatVector();
		//SpatVector(const SpatVector &x);
		SpatVector(SpatGeom g);
		SpatVector(SpatExtent e, std::string crs);
		SpatVector(std::vector<double> x, std::vector<double> y, SpatGeomType g, std::string crs);
		SpatVector(std::vector<std::string> wkt);
		virtual ~SpatVector(){}

		SpatGeom window; // for point patterns, must be polygon

		std::vector<std::string> get_names();
		void set_names(std::vector<std::string> s);
		unsigned nrow();
		unsigned ncol();
		unsigned nxy();

		SpatVector deepCopy() {return *this;}

		SpatExtent getExtent();
//		bool is_geographic();
		bool is_lonlat();
		bool could_be_lonlat();
		std::string type();
		SpatGeomType getGType(std::string &type);

		//std::vector<std::string> getCRS();
		//void setCRS(std::vector<std::string> _crs);

		bool setSRS(std::string _srs) {
			std::string msg;
			if (!srs.set(_srs, msg)){
				addWarning("Cannot set SRS to vector: "+ msg);
				return false;
			}
			return true;	
		}

		std::string getSRS(std::string x) {
			return srs.get(x);
		}

		SpatGeom getGeom(unsigned i);
		bool addGeom(SpatGeom p);
		bool setGeom(SpatGeom p);
		bool replaceGeom(SpatGeom p, unsigned i);
		std::vector<std::vector<double>> getGeometry();
		SpatDataFrame getGeometryDF();
		std::vector<std::string> getGeometryWKT();
		void computeExtent();

		size_t nparts(bool holes);

		size_t ncoords();
		std::vector<std::vector<double>> coordinates();

		SpatVector project(std::string crs, bool partial);
		std::vector<double> project_xy(std::vector<double> x, std::vector<double> y, std::string fromCRS, std::string toCRS);

		SpatVector subset_cols(int i);
		SpatVector subset_cols(std::vector<int> range);
		SpatVector subset_rows(int i);
		SpatVector subset_rows(std::vector<int> range);
		SpatVector subset_rows(std::vector<unsigned> range);
		SpatVector remove_rows(std::vector<unsigned> range);

		void setGeometry(std::string type, std::vector<unsigned> gid, std::vector<unsigned> part, std::vector<double> x, std::vector<double> y, std::vector<unsigned> hole);
		void setPointsGeometry(std::vector<double> &x, std::vector<double> &y);
		void setPointsDF(SpatDataFrame &x, std::vector<unsigned> geo, std::string crs, bool keepgeom);

		std::vector<double> area(std::string unit, bool transform, std::vector<double> mask);

		void reserve(size_t n);
		std::vector<double> length();
		std::vector<double> distance(SpatVector x, bool pairwise, std::string unit);
		std::vector<double> pointdistance(const std::vector<double>& px, const std::vector<double>& py, const std::vector<double>& sx, const std::vector<double>& sy, bool pairwise, double m, bool lonlat);

//		std::vector<double> pointdistance_seq(const std::vector<double>& px, const std::vector<double>& py, double m, bool lonlat);


		std::vector<double> distance(bool sequential, std::string unit);
		std::vector<double> linedistLonLat(SpatVector x);

		std::vector<std::vector<size_t>> knearest(size_t k);

		size_t size();
		bool empty();
		SpatVector as_lines();
		SpatVector as_points(bool multi, bool skiplast=false);
		SpatVector remove_holes();
		SpatVector get_holes();
		SpatVector set_holes(SpatVector x, size_t i);
		SpatVector remove_duplicate_nodes(int digits);
		
		bool read(std::string fname, std::string layer, std::string query, std::vector<double> extent, SpatVector filter, bool as_proxy, std::string what);
		
		bool write(std::string filename, std::string lyrname, std::string driver, bool append, bool overwrite, std::vector<std::string>);

		void make_CCW();
		
#ifdef useGDAL
		GDALDataset* write_ogr(std::string filename, std::string lyrname, std::string driver, bool append, bool overwrite, std::vector<std::string> options);
		GDALDataset* GDAL_ds();
		bool read_ogr(GDALDataset *&poDS, std::string layer, std::string query, std::vector<double> extent, SpatVector filter, bool as_proxy, std::string what);
		SpatVector fromDS(GDALDataset *poDS);
		bool ogr_geoms(std::vector<OGRGeometryH> &ogrgeoms, std::string &message);		
		bool delete_layers(std::string filename, std::vector<std::string> layers, bool return_error);		
		std::vector<std::string> layer_names(std::string filename);	
#endif

// attributes
		std::vector<double> getDv(unsigned i);
		std::vector<long> getIv(unsigned i);
		std::vector<std::string> getSv(unsigned i);
		std::vector<unsigned> getItype();
		std::vector<unsigned> getIplace();

		void add_column(unsigned dtype, std::string name) {
			df.add_column(dtype, name);
		};
		template <typename T>
		bool add_column(std::vector<T> x, std::string name) {
			return df.add_column(x, name);
		}
		bool add_column_bool(std::vector<int> x, std::string name) {
			return df.add_column_bool(x, name);
		}
		bool add_column_time(std::vector<SpatTime_t> x, std::string name, std::string step, std::string zone) {
			return df.add_column_time(x, name, step, zone);
		}
		bool add_column_factor(SpatFactor x, std::string name) {
			return df.add_column(x, name);
		}

		void remove_df() {
			SpatDataFrame empty;
			df = empty;
		};
		
		bool set_df(SpatDataFrame x) {
			if (x.nrow() != nrow()) {
				setError("nrow dataframe does not match nrow geometry");
				return false;
			}
			df = x;
			return true;
		};


		bool remove_column(std::string field) {
			return df.remove_column(field);
		};
		bool remove_column(int i) {
			return df.remove_column(i);
		};
		std::vector<std::string> get_datatypes() {
			return df.get_datatypes();
		}

		SpatMessages msg;
		void setError(std::string s) { msg.setError(s); }
		void addWarning(std::string s) { msg.addWarning(s); }
		bool hasError() { return msg.has_error; }
		bool hasWarning() { return msg.has_warning; }
		std::vector<std::string> getWarnings() { return msg.getWarnings();}
		std::string getError() { return msg.getError();}

		SpatVector append(SpatVector x, bool ignorecrs);
		SpatVector disaggregate(bool segments);
		SpatVector shift(double x, double y);
		SpatVector rescale(double fx, double fy, double x0, double y0);
		SpatVector transpose();
		SpatVector flip(bool vertical);	
		SpatVector rotate(double angle, std::vector<double> x0, std::vector<double> y0);
		SpatVector normalize_longitude();
		SpatVector rotate_longitude(double longitude, bool left);

		std::vector<std::vector<double>> linesNA();
		std::vector<std::vector<std::vector<double>>> linesList();
		std::vector<std::vector<std::vector<std::vector<double>>>> polygonsList();

//ogr 
		std::vector<bool> is_valid();
		SpatVector make_valid();

//geos
		SpatVector make_valid2();

		std::vector<bool> geos_isvalid();
		std::vector<std::string> geos_isvalid_msg();
		std::vector<std::string> wkt();
		std::vector<std::string> wkb();
		std::vector<std::string> hex();
		SpatVector from_hex(std::vector<std::string> x, std::string srs);
		SpatVector make_nodes();
		SpatVector polygonize();
		SpatVector normalize();
		SpatVector boundary();
		SpatVector line_merge();
		SpatVector simplify(double tolerance, bool preserveTopology);
		SpatVector shared_paths(bool index);
		SpatVector shared_paths(SpatVector x, bool index);
		SpatVector snap(double tolerance);
		SpatVector snapto(SpatVector y, double tolerance);
		SpatVector thin(double threshold);

		SpatVector allerretour();
		SpatVectorCollection bienvenue();
		SpatVector aggregate(bool dissolve);
		SpatVector aggregate(std::string field, bool dissolve);

		SpatVector buffer2(std::vector<double> d, unsigned quadsegs);
		SpatVector buffer3(std::vector<double> d, unsigned quadsegs);
//		SpatVector buffer4(double d);
//		SpatVector bufferclip(std::vector<double> d, std::string jointype, double miter_limit, int precision, double arc_tolerance);
		
		SpatVector buffer(std::vector<double> d, unsigned quadsegs, std::string capstyle, std::string joinstyle, double mitrelimit, bool singlesided);
		SpatVector buffer_lonlat(std::string vt, std::vector<double> d, unsigned quadsegs);


		SpatVector point_buffer(std::vector<double>	d, unsigned quadsegs, bool no_multipolygons, bool wrap);

		SpatVector centroid(bool check_lonlat);
		SpatVector point_on_surface(bool check_lonlat);

		SpatVector crop(SpatExtent e, bool wrap);
		SpatVector crop(SpatVector e);
		SpatVector voronoi(SpatVector bnd, double tolerance, int onlyEdges);		
		SpatVector delaunay(double tolerance, int onlyEdges);		
		SpatVector hull(std::string htype, std::string by="");
		SpatVector intersect(SpatVector v, bool values);
		SpatVector unite(SpatVector v);
		SpatVector unite();
		SpatVector erase_agg(SpatVector v);
		SpatVector erase(SpatVector v);
		SpatVector erase(bool sequential);
		SpatVector elongate(double length, bool flat);
		SpatVector mask(SpatVector x, bool inverse);
		SpatVector gaps();		
		SpatVector cover(SpatVector v, bool identity, bool expand);
		SpatVectorCollection split(std::string field);
		SpatVector symdif(SpatVector v);
		SpatVector set_precision(double gridSize);
		std::vector<std::vector<unsigned>> index_2d(SpatVector v);
		std::vector<std::vector<unsigned>> index_sparse(SpatVector v);

		std::vector<std::vector<double>> which_relate(SpatVector v, std::string relation, bool narm);
		std::vector<std::vector<double>> which_relate(std::string relation, bool narm);
		std::vector<bool> is_related(SpatVector v, std::string relation);
//		std::vector<int> relate(SpatVector v, std::string relation);
		std::vector<int> relate(SpatVector v, std::string relation, bool prepared, bool index);
		std::vector<int> relate(std::string relation, bool symmetrical);
		std::vector<int> relateFirst(SpatVector v, std::string relation);
		std::vector<unsigned> equals_exact(SpatVector v, double tol);
		std::vector<unsigned> equals_exact(bool symmetrical, double tol);

		std::vector<double> geos_distance(SpatVector v, bool parallel, std::string fun);
		std::vector<double> geos_distance(bool sequential, std::string fun);

		SpatVector nearest_point(SpatVector v, bool parallel);
		SpatVector nearest_point();
		std::vector<int> nearest_geometry(SpatVector v);
		
		SpatVector sample(unsigned n, std::string method, unsigned seed);
		SpatVector sample_geom(std::vector<unsigned> n, std::string method, unsigned seed);

		SpatVector clearance();
		SpatVector width();

		SpatVector unaryunion();

		SpatVector cbind(SpatDataFrame d);
		void fix_lonlat_overflow();
		SpatVector cross_dateline(bool &fixed);
		SpatVector densify(double interval, bool adjust, bool ignorelonlat);
		SpatVector round(int digits);
		std::vector<unsigned> nullGeoms();
		std::vector<bool> naGeoms();
				
};



class SpatVectorCollection {


	public:
		virtual ~SpatVectorCollection(){}
		SpatVectorCollection();
		SpatVectorCollection(std::string filename, std::string layer, std::string query, std::vector<double> extent, SpatVector filter);
		
		
		SpatVectorCollection deepCopy() { return *this; }
		bool read(std::string fname, std::string layer, std::string query, std::vector<double> extent, SpatVector filter);
		
		bool read_ogr(GDALDataset *&poDS, std::string layer, std::string query, std::vector<double> extent, SpatVector filter);

//		SpatVectorCollection create(std::string filename);

		std::vector<SpatVector> v;
		std::vector<std::string> names;
		std::vector<std::string> getNames() { return names;}
		bool setNames(std::vector<std::string> nms, bool make_valid=false);

		SpatMessages msg;
		void setError(std::string s) { msg.setError(s); }
		void addWarning(std::string s) { msg.addWarning(s); }
		bool hasError() { return msg.has_error; }
		bool hasWarning() { return msg.has_warning; }
		std::vector<std::string> getWarnings() { return msg.getWarnings();}
		std::string getError() { return msg.getError();}

		size_t size() { return v.size(); }
		bool empty() { return v.empty(); }
		void reserve(size_t n) { v.reserve(n); names.reserve(n); }
		void resize(size_t n) { v.resize(n); names.resize(n); }
		void push_back(SpatVector x) {
			v.push_back(x); 
			names.push_back("");
		};
		bool replace(SpatVector x, size_t i) { 
			if (i < size()) {
				v[i] = x; 
				return true;
			} else {
				return false;
			}
		}
		SpatVectorCollection subset(std::vector<size_t> i) { 
			SpatVectorCollection out;
			for (size_t j=0; j<size(); j++) {
				if (i[j] < size()) {
					out.push_back(v[i[j]]); 
				} 
			}
			return out;
		}

		SpatVector get(size_t i) { 
			SpatVector out;
			out.msg = msg;
			if (size() == 0) {
				out.addWarning("empty SpatVector");
				return out;
			}
			if (i < size()) {
				out = v[i];
			} else {
				out.setError("invalid index");
			}
			return out;
		}
		
		SpatVector append();
		SpatVectorCollection from_hex_col(std::vector<std::string> x, std::string srs);
		
};



class SpatVectorProxy {
	public:
		SpatVector v;
		SpatVectorProxy(){}
		virtual ~SpatVectorProxy(){}
		SpatVectorProxy deepCopy() {return *this;}
		SpatVector query_filter(std::string query, std::vector<double> extent, SpatVector filter);
};


#endif // SPATVECTOR_GUARD
