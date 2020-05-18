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

#include "spatBase.h"
#include "spatDataframe.h"
#include "spatMessages.h"

#ifdef useGDAL
#include "gdal_priv.h"
#endif


enum SpatGeomType { points, lines, polygons, unknown };

class SpatHole {
	public:
		std::vector<double> x, y;
		SpatExtent extent;
		SpatHole();
		SpatHole(std::vector<double> X, std::vector<double> Y);
		size_t size() { return x.size(); }	
};

class SpatPart {
	public:
		std::vector<double> x, y; //, z;
		SpatExtent extent;
		SpatPart();
		SpatPart(std::vector<double> X, std::vector<double> Y);
		SpatPart(double X, double Y);
		size_t size() { return x.size(); }
		// for POLYGONS only
		std::vector< SpatHole > holes;
		bool addHole(std::vector<double> X, std::vector<double> Y);
		bool addHole(SpatHole h);
		SpatHole getHole(unsigned i) { return( holes[i] ) ; }
		bool hasHoles() { return holes.size() > 0;}
		unsigned nHoles() { return holes.size();}
};


class SpatGeom {
	public:
		SpatGeomType gtype;
		std::vector<SpatPart> parts;
		SpatExtent extent;
		SpatGeom();
		SpatGeom(SpatGeomType g);
		SpatGeom(SpatPart p);
		bool unite(SpatGeom g);
		bool addPart(SpatPart p);
		bool addHole(SpatHole h);
		bool setPart(SpatPart p, unsigned i);
		SpatPart getPart(unsigned i);
		double area_plane();
		double area_lonlat(double a, double f);
		double length_plane();
		double length_lonlat(double a, double f);
		unsigned size() { return parts.size(); };
};


class SpatLayer {
	public:
		std::vector<SpatGeom> geoms;
		SpatExtent extent;
		SpatDataFrame df;
		//std::vector<std::string> crs;
		SpatSRS srs;
};

class SpatVector {

	public:

		SpatVector();
		//SpatVector(const SpatVector &x);
		SpatVector(SpatGeom g);
		SpatVector(SpatExtent e, std::string crs);

		SpatLayer lyr;
		std::vector<SpatLayer> lyrs;

		std::vector<std::string> get_names();
		void set_names(std::vector<std::string> s);
		unsigned nrow();
		unsigned ncol();
		unsigned nxy();

		SpatExtent getExtent();
		bool is_lonlat();
		bool could_be_lonlat();
		std::string type();
		SpatGeomType getGType(std::string &type);

		//std::vector<std::string> getCRS();
		//void setCRS(std::vector<std::string> _crs);


		bool setSRS(std::string _srs) {
			std::string msg;
			if (!lyr.srs.set(_srs, msg)){
				addWarning("Cannot set SRS to vector: "+ msg);
				return false;
			}
			return true;	
		}

/*
#ifdef useGDAL	
		bool setSRS(OGRSpatialReference *poSRS, std::string &msg) {
			if (!lyr.srs.set(poSRS, msg)){
				addWarning("Cannot set SRS to vector: "+ msg);
				return false;
			}
			return true;				
		}
#endif		
*/

		std::string getSRS(std::string x) {
			return lyr.srs.get(x);
		}
		//std::string getPRJ();
		//void setPRJ(std::string PRJ);

		SpatGeom getGeom(unsigned i);
		bool addGeom(SpatGeom p);
		bool setGeom(SpatGeom p);
		SpatDataFrame getGeometryDF();
		std::vector<std::vector<double>> coordinates();

		SpatVector project(std::string crs);
		//std::vector<std::vector<double>> test(std::vector<double> x, std::vector<double> y, std::string fromcrs, std::string tocrs);

		SpatVector subset_cols(int i);
		SpatVector subset_cols(std::vector<int> range);
		SpatVector subset_rows(int i);
		SpatVector subset_rows(std::vector<int> range);

		void setGeometry(std::string type, std::vector<unsigned> gid, std::vector<unsigned> part, std::vector<double> x, std::vector<double> y, std::vector<unsigned> hole);

		std::vector<double> area();
		std::vector<double> length();
		SpatDataFrame distance(SpatVector x, bool pairwise);
		SpatDataFrame distance();

		size_t size();
		SpatVector as_lines();
		SpatVector as_points();

		bool read(std::string fname);
		
		bool write(std::string filename, std::string lyrname, std::string driver, bool overwrite);
#ifdef useGDAL
		GDALDataset* write_ogr(std::string filename, std::string lyrname, std::string driver, bool overwrite);
		GDALDataset* GDAL_ds();
		bool read_ogr(GDALDataset *poDS);
		SpatVector fromDS(GDALDataset *poDS);
#endif

// attributes
		std::vector<double> getDv(unsigned i);
		std::vector<long> getIv(unsigned i);
		std::vector<std::string> getSv(unsigned i);
		std::vector<unsigned> getItype();
		std::vector<unsigned> getIplace();

		void add_column(unsigned dtype, std::string name) {
			lyr.df.add_column(dtype, name);
		};
		template <typename T>
		bool add_column(std::vector<T> x, std::string name) {
			return lyr.df.add_column(x, name);
		}

		bool remove_column(std::string field) {
			return lyr.df.remove_column(field);
		};
		bool remove_column(int i) {
			return lyr.df.remove_column(i);
		};


		SpatMessages msg;
		void setError(std::string s) { msg.setError(s); }
		void addWarning(std::string s) { msg.addWarning(s); }
		bool hasError() { return msg.has_error; }
		bool hasWarning() { return msg.has_warning; }

		SpatVector point_buffer(double d, unsigned quadsegs);

		std::vector<bool> is_valid();
		SpatVector make_valid();
        SpatVector buffer(double d, unsigned segments, unsigned capstyle);

		SpatVector aggregate(std::string field, bool dissolve);
		SpatVector disaggregate();
//geos
        SpatVector buffer2(double d, unsigned segments, unsigned capstyle);
		SpatVector intersect(SpatVector v);
};

/*

class SpatVector {
	public:
		std::vector<SpatLayer> lyrs;
		unsigned layertypes();
		SpatMessages msg;
		void setError(std::string s) { msg.setError(s); }
		void addWarning(std::string s) { msg.addWarning(s); }
};

*/

