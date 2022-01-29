// Copyright (c) 2018-2020  Robert J. Hijmans
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

//#include "spatBase.h"
#include "spatDataframe.h"
//#include "spatMessages.h"

#ifdef useGDAL
#include "gdal_priv.h"
#endif


enum SpatGeomType { points, lines, polygons, unknown };


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
		//holes, polygons only
		bool addHole(std::vector<double> X, std::vector<double> Y);
		bool addHole(SpatHole h);
		SpatHole getHole(unsigned i) { return( holes[i] ) ; }
		bool hasHoles() { return holes.size() > 0;}
		unsigned nHoles() { return holes.size();}
};


class SpatGeom {
	public:
		//constructors
		SpatGeom();
		SpatGeom(SpatGeomType g);
		SpatGeom(SpatPart p);
		virtual ~SpatGeom(){}

		SpatGeomType gtype = unknown;
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
		void remove_duplicate_nodes(int digits);
		size_t ncoords();
		std::vector<std::vector<double>> coordinates();

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
		void reserve(size_t n) {
			geoms.reserve(n);
		}	

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

		size_t ncoords();
		std::vector<std::vector<double>> coordinates();

		SpatVector project(std::string crs);

		SpatVector subset_cols(int i);
		SpatVector subset_cols(std::vector<int> range);
		SpatVector subset_rows(int i);
		SpatVector subset_rows(std::vector<int> range);
		SpatVector subset_rows(std::vector<unsigned> range);
		SpatVector remove_rows(std::vector<unsigned> range);

		void setGeometry(std::string type, std::vector<unsigned> gid, std::vector<unsigned> part, std::vector<double> x, std::vector<double> y, std::vector<unsigned> hole);
		void setPointsGeometry(std::vector<double> &x, std::vector<double> &y);
		void setPointsDF(SpatDataFrame &x, std::vector<unsigned> geo, std::string crs);

		std::vector<double> area(std::string unit, bool transform, std::vector<double> mask);

		std::vector<double> length();
		std::vector<double> distance(SpatVector x, bool pairwise);
		std::vector<double> distance(bool sequential);
		std::vector<double> linedistLonLat(SpatVector pts);

		std::vector<std::vector<size_t>> knearest(size_t k);

		size_t size();
		SpatVector as_lines();
		SpatVector as_points(bool multi, bool skiplast=false);
		SpatVector remove_holes();
		SpatVector get_holes();
		SpatVector set_holes(SpatVector x, size_t i);
		SpatVector remove_duplicate_nodes(int digits);

		bool read(std::string fname, std::string layer, std::string query, std::vector<double> extent, SpatVector filter, bool as_proxy);
		
		bool write(std::string filename, std::string lyrname, std::string driver, bool append, bool overwrite, std::vector<std::string>);
		
#ifdef useGDAL
		GDALDataset* write_ogr(std::string filename, std::string lyrname, std::string driver, bool append, bool overwrite, std::vector<std::string> options);
		GDALDataset* GDAL_ds();
		bool read_ogr(GDALDataset *poDS, std::string layer, std::string query, std::vector<double> extent, SpatVector filter, bool as_proxy);
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

		void remove_df() {
			SpatDataFrame empty;
			df = empty;
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
		std::string getWarnings() { return msg.getWarnings();}
		std::string getError() { return msg.getError();}

		SpatVector append(SpatVector x, bool ignorecrs);
		SpatVector disaggregate();
		SpatVector shift(double x, double y);
		SpatVector rescale(double fx, double fy, double x0, double y0);
		SpatVector transpose();
		SpatVector flip(bool vertical);	
		SpatVector rotate(double angle, double x0, double y0);

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
		SpatVector shared_paths();
		SpatVector snap(double tolerance);
		SpatVector snapto(SpatVector y, double tolerance);

		SpatVector allerretour();
		SpatVectorCollection bienvenue();
		SpatVector aggregate(bool dissolve);
		SpatVector aggregate(std::string field, bool dissolve);

        SpatVector buffer(std::vector<double> d, unsigned quadsegs);
		SpatVector point_buffer(std::vector<double>	 d, unsigned quadsegs, bool no_multipolygons);

		SpatVector centroid(bool check_lonlat);
		SpatVector crop(SpatExtent e);
		SpatVector crop(SpatVector e);
		SpatVector voronoi(SpatVector e, double tolerance, int onlyEdges);		
		SpatVector delauny(double tolerance, int onlyEdges);		
		SpatVector hull(std::string htype, std::string by="");
		SpatVector intersect(SpatVector v);
		SpatVector unite(SpatVector v);
		SpatVector unite();
		SpatVector erase_agg(SpatVector v);
		SpatVector erase(SpatVector v);
		SpatVector erase();
		SpatVector mask(SpatVector x, bool inverse);
		SpatVector gaps();		
		SpatVector cover(SpatVector v, bool identity);
		SpatVectorCollection split(std::string field);
		SpatVector symdif(SpatVector v);
		std::vector<bool> is_related(SpatVector v, std::string relation);
		std::vector<int> relate(SpatVector v, std::string relation);
		std::vector<int> relate(std::string relation, bool symmetrical);
		std::vector<int> relateFirst(SpatVector v, std::string relation);
		std::vector<double> geos_distance(SpatVector v, bool parallel);
		std::vector<double> geos_distance(bool sequential);


		SpatVector nearest_point(SpatVector v, bool parallel);
		SpatVector nearest_point();
		SpatVector sample(unsigned n, std::string method, unsigned seed);
		SpatVector sample_geom(std::vector<unsigned> n, std::string method, unsigned seed);

		SpatVector clearance();
		SpatVector width();

		SpatVector unaryunion();

		SpatVector cbind(SpatDataFrame d);
		void fix_lonlat_overflow();
		SpatVector cross_dateline(bool &fixed);
		SpatVector densify(double interval, bool adjust);
		SpatVector round(int digits);
};



class SpatVectorCollection {

	private:
		std::vector<SpatVector> v;

	public:
		virtual ~SpatVectorCollection(){}
		SpatVectorCollection deepCopy() { return *this; }

		SpatMessages msg;
		void setError(std::string s) { msg.setError(s); }
		void addWarning(std::string s) { msg.addWarning(s); }
		bool hasError() { return msg.has_error; }
		bool hasWarning() { return msg.has_warning; }
		std::string getWarnings() { return msg.getWarnings();}
		std::string getError() { return msg.getError();}

		size_t size() { return v.size(); }
		void push_back(SpatVector x) { v.push_back(x); };
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
			if (i < size()) {
				out = v[i];
			} else {
				out.setError("invalid index");
			}
			return out;
		}
		
		SpatVector append();
		
};



class SpatVectorProxy {
	public:
		SpatVector v;
		SpatVectorProxy(){}
		virtual ~SpatVectorProxy(){}
		SpatVectorProxy deepCopy() {return *this;}
		SpatVector query_filter(std::string query, std::vector<double> extent, SpatVector filter);
};

