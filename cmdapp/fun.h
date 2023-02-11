#include <string>

SpatRaster make_raster() {
    SpatExtent e(-180,180,-90,90);
    SpatRaster r(180, 365, 1, e, "");
    std::vector<double> v;
    for (size_t i = 0; i < r.ncell(); i++) {
       v.push_back(i+1);
    }
    SpatOptions opt;
    r.setValues(v, opt);
    return r;
}


SpatRaster aggregate(std::vector<std::string> args) {
    SpatOptions opts;
    SpatRaster out;
    if (args.size() < 5) {
	out.setError("usage: aggregate input output fact fun=mean narm=true");
	return(out);
    }

    SpatRaster input(args[2], {-1}, {""}, {""}, {""});
    opts.set_filenames({args[3]});
    // need to split the fact arguments
    // for now assume a single argument
    std::vector<unsigned> fact(1);
    fact[0] = stoi(args[4]);
    std::string fun;
    if (args.size() < 6) {
	fun = "mean";
    } else {
        fun = args[5];
    }
    bool narm;
    if (args.size() < 7) {
       narm = true;
    } else {
       narm = stoi(args[6]);
    }

    out = input.aggregate(fact, fun, narm, opts);
    //std::cout << "out = input.aggregate(" << fact[0] << ", " << fun << "," << narm << ", " << opts.filename << std::endl;
    return(out);
}



    //std::vector<double> x = {-180, -179.9, 190, 0.1};
    //std::vector<double> y = {90, 89.9, 90, 0.1};
    //std::vector<std::vector<double>> out = r.bilinearValues(x, y);
    //std::vector<std::vector<double>> out = r.sampleRandom(10,18);
    //showValues(out);

/*
    std::string f = "c:/temp/agfile.grd";
    opts.filename = f;
    opts.datatype = "FLT4S";
    opts.overwrite = true;

    SpatExtent e(-180,-180,-90,90);
    SpatRaster x(360, 180, 1, e, "");
    std::vector<double> v;
    for (size_t i = 0; i < x.ncell(); i++) {
       v.push_back(i+1);
    }
    x.setValues(v);
    std::vector<unsigned> fact(1);
    fact[0] = 10;
//    fact[1] = 10;
    SpatRaster out = x.aggregate(fact, "mean", true, opts);
    show(out);

    SpatOptions opts;
    SpatExtent e(-180,-150,60,90);
    SpatRaster x(5, 5, 1, e, "");
    std::vector<double> v;
    for (size_t i = 0; i < x.ncell(); i++) {
       v.push_back(i+1);
    }

    std::string f = "c:/temp/file.grd";
    opts.filename = f;
    opts.datatype = "FLT4S";
    opts.overwrite = true;

    x.writeStart(opts);
    x.writeValues(v,0,5,0,5);
    x.writeStop();

    std::vector<double> vv = x.getValues();

    showValues(v);
    std::cout << "\n";

    showValues(x);
    std::cout << "\n";

    show(x);

    SpatRaster out = SpatRaster(f);
    show(out);
    showValues(out);


    SpatExtent e(0,10,0,5);
    SpatRaster r(10, 10, 1, e, "+proj=longlat +datum=WGS84");
    std::vector<double> v(100, NAN);
    v[55] = 1;
    r.setValues(v);
    SpatRaster d = r.distance(opts);
    show(d);


    SpatRaster z(10, 10, 1, e, "");
    v.resize(0);
    for (size_t j = 0; j < 4; j++) {
        for (size_t i=0; i<25; i++) {
            v.push_back(j);
        }
    }
    z.setValues(v);

    std::vector<std::vector<double>> zu = z.unique(false);
    showValues(zu);

    SpatDataFrame d = r.zonal(z, "max", true, opts);
    std::cout << "\n";
    showValues(d.getD(1));


    SpatExtent e(-180,-150,60,90);
    SpatRaster x(30, 30, 1, e, "");
    SpatExtent ee(-100,-50,30,50);
    SpatRaster y(20, 50, 1, ee, "");

    std::vector<double> v;
    for (size_t i = 0; i < x.ncell(); i++) {
       v.push_back(i);
    }
    x.setValues(v);
    v.resize(0);
    for (size_t i = 0; i < y.ncell(); i++) {
       v.push_back(i);
    }
    y.setValues(v);

    show(x);
    show(y);

    SpatRasterCollection rc;
    rc.push_back(x);
    rc.push_back(y);
    SpatRaster mr = rc.merge(opts);

    show(mr);



    SpatOptions opts;
    SpatExtent e(0,10,0,5);
    SpatRaster r(2, 5, 3, e, "");
    std::vector<double> v;
    for (size_t j = 0; j < 6; j++) {
        for (size_t i=0; i<5; i++) {
            v.push_back(i+j);
        }
    }
    r.setValues(v);
    show(r);
    opts.filename = "c:/temp/file.grd";
    opts.datatype = "FLT4S";
    opts.overwrite = true;
    r.writeRaster(opts);
    show(r);
//    showValues(r.source[0].range_min);

  //  SpatRaster out = r.arith(5, ">", opts);
  //  showValues(out);

 //  SpatRaster r("d:/test.grd");
//    show(r);
//     showValues(r);



     SpatVector sv;

    std::vector<unsigned> object = { 1,1,1,1};
    std::vector<unsigned> part = { 1,1,1,1};
    std::vector<double> x = { -180, -140, 10, -140};
    std::vector<double> y = { -20,55,0,-60};
    std::vector<unsigned> hole = { 0, 0, 0, 0, 0};

    sv.setGeometry("polygons", object, part, x, y, hole);
    show(sv);

    SpatDataFrame v = sv.getGeometryDF();
        std::string names = join(v.names, ", ");
        std::cout << "names   : " << names << std::endl;


   	unsigned n = v.ncol();
	std::vector<unsigned> itype = v.itype;
	for (size_t i=0; i < n; i++) {
		if (itype[i] == 0) {
			std::vector<double> dv = v.getD(i);
            std::cout << "double  : " << dv[0] << std::endl;

		} else if (itype[i] == 1) {
			v.getI(i);
		} else if (itype[i] == 2) {
			v.getS(i);
		} else if (itype[i] == 3) {
			std::vector<unsigned> uv = v.getU(i);
           std::cout << "unsigned  : " << uv[0] << std::endl;
		} else {
			v.getB(i);
		}
	}


    //std::vector<std::vector<double>> vv = r.extractVector(sv);




 std::vector<std::vector<double>> out = r.unique(0);
 for (size_t i=0; i<out.size(); i++) {
        for (size_t j=0; j<r.nlyr(); j++) {
            std::cout << out[i][j] << ", ";
        }
        std::cout << "\n";
 }
   //  std::vector<unsigned> fact = {2,2,1};
   //  SpatRaster z = r.disaggregate(fact, opts);
   //  showValues(r);
//     showValues(r);
//     showValues(z);

    SpatExtent e(0,10,0,10);
    SpatRaster r(10, 10, 1, e, "");
    std::vector<double> v;
    for (size_t i=0; i<r.size(); i++) { v.push_back(i); }
    for (size_t i=0; i<r.size(); i=i+5) { v[i] = NAN; }

    r.setValues(v);
    showValues(r);
    show(r);
    std::vector<double> cells = {0,11,31,14};
    std::vector<std::vector<double>> d = r.extractCell(cells);
    showValues(d);
    std::vector<std::vector<double>> xy = r.xyFromCell(cells);
    d = r.extractXY(xy[0], xy[1], "simple");
    showValues(d);



    SpatVector svp = r.as_points(false, true);
    show(svp);
    svp = r.as_points(false, false);
    show(svp);
    svp = r.as_points(true, true);
    show(svp);
    svp = r.as_points(true, false);
    show(svp);
    std::cout << svp.lyr.df.dv[0][0] << ", " << svp.lyr.df.dv[0][99] << "\n";


    std::vector<double> cells = {0, 15};
    std::vector<std::vector<double>> adj = r.adjacent(cells, "bishop", true);
    showValues(adj);
    std::string fn = "c:/temp/test222.grd";
     opts.set_overwrite(true);
     opts.set_filename(fn);
     r.writeRaster(opts);
     SpatRaster x = SpatRaster(fn);
     show(x);
     showValues(x);


    //showValues(b, 2,2);

    SpatExtent e2(0,10,0,10);
    SpatRaster z(6, 6, 1, e2, "");

    //showValues(d, 6, 6);
//    d = r.extractXY(xy[0], xy[1], "bilinear");
//    showValues(d, 6, 6);
//    showValues(r);

    SpatRaster out = r.warp(z, "bilinear", opt);
//    show(out);
    showValues(out);


    std::vector<unsigned> fact = {2,3};
    SpatRaster d = r.disaggregate(fact, opt);
    showValues(d);

    v = out.getValues();
    for (size_t i=0; i<v.size(); i++) { cout << v[i] << " "; }

    SpatExtent e(-180,180,-90,90);
    SpatRaster r(10, 10, 1, e, "");
    show(r);

    SpatExtent x(-180,-144,-90,90);
    SpatOptions opt;
    std::string s = "near";
    SpatRaster a = r.crop(x, s, opt);
    show(a);
    SpatVector v = a.makePolygons(false, false);
    std::vector<double> z = v.area();
    for (size_t i =0; i<z.size(); i++) {
        cout << z[i] << " ";
    }
    cout << "\n ";


    SpatExtent e(0,50,0,50);



    SpatRaster r(10, 10, 1, e, "");
    std::vector<double> v;
    for (size_t i=0; i<r.size(); i++) { v.push_back(NAN); }
    v[23] = 1;
    r.setValues(v);

    SpatRaster x = r.gridDistance("", "", "", false);
    show(x);
    v = x.getValues();
    for (size_t i=0; i <10; i++) {
        for (size_t j=0; j <10; j++) {
             size_t k = i * 10 + j;
             cout << v[k] << ", "   ;
        }
        cout  << "\n" ;
    }


    std::vector<std::vector<double>> rcl(3);
    rcl[0] = {-1,10,50};
    rcl[1] = {10,50,2000};
    rcl[2] = {1,2,3};

    SpatRaster x = r.reclassify(rcl, 0, false, "", "", "", false);
    show(x);
 //   SpatRaster x = r.arith(10, "*", "", "", "", false);
 //   show(x);
 //   std::vector<double> v = x.getValues();
 //   for (size_t i=0; i<v.size(); i++) {
 //           if (v[i] > 10)  std::cout << v[i] << " ";
 //   }


//    SpatRaster g = r.sampleRegular(10);
//    show(g);
	SpatLayer vec;
    std::vector<unsigned> id = {0,1,2,3};
    std::vector<unsigned> part = {0,0,0,0};
    std::vector<double> x = {1,2,3,1};
    std::vector<double> y = {3,1,2,3};
    std::vector<bool> hole = {0,0,0,0};
	vec.setGeometry("point", id, part, x, y, hole);
    unsigned n = vec.size();

    SpatDataFrame g = vec.getGeometryDF();
    cout << g.nrow() << "\n";

	SpatDataFrame out;

	out.add_column(1, "geom");
	out.add_column(1, "part");
	out.add_column(0, "x");
	out.add_column(0, "y");
	out.add_column(1, "hole");

	out.resize(n);

    cout << out.nrow() << "\n";

//    std::vector<double> x = df.dv[0]
//	 for (size_t i=0; i<d.size(); i++) {
//    cout << r.source[0].values[i] << "\n";
//    cout << d[i] << "\n";
// }




    SpatExtent e;
    SpatRaster x, y, z;
    SpatRaster r(3, 3, 1, e, "");
    std::vector<double> v = {1,2,3,4,5,6,7,8,9};
    r.setValues(v);
    show(r);
    std::vector<double> i = {1, 2};
    std::vector<double>  d = r.extractCell(i);


//    SpatRaster y = r.summary("max", true, "", true);
//   show(y);

    cout << NA<short>::value << "\n";
    cout << NA<long>::value << "\n";
    cout << NA<long long>::value << "\n";
    cout << NA<float>::value << "\n";
    cout << NA<double>::value << "\n";
    cout << NA<long double>::value << "\n";
    cout << NA<bool>::value << "\n";

    cout << (int)NA<long double>::value << "\n";


    string f = "C:/temp";
    SpatRaster z(f);
    show(z);
    SpatRaster x = z.cum("sum", false, "", false);


std::vector<unsigned> rcl = {45, 90 ,1};
std::vector<double> e = {-180, 180, -90 ,90};
SpatRaster r(rcl, e, "+proj=longlat +datum=WGS84");
std::vector<double> v(r.ncell());
iota(v.begin(), v.end(), 1);

r.setValues(v);
show(r);
SpatExtent ee(-160, 10, 30, 60);
SpatRaster rc = r.crop(ee);
show(rc);



  void make_unique(std::vector<string> &s);
    void make_valid(std::vector<string> &s);

    std::vector<string> a = {" 1", "1", " 2 ", "b ", "c", "2", "2 ", "e", "e", "b", "b", "c", "d"};
    make_valid(a);
    make_unique(a);
    for (size_t i=0; i<a.size(); i++) {
        cout << a[i] << ", ";
    }
    cout << "\n";


    string f = "C:/soft/R/R-3.5.1/library/terra/external/rlogo.grd";
    SpatRaster z(f);
    show(z);
//std::vector<unsigned> lyrs = {0};
//    z = z.subset( lyrs, "", FALSE );
//    show(z);

    std::vector<unsigned> lyrs = {-1, 0,2,4};
    z = z.subset( lyrs, "", FALSE );
    show(z);



    string filename1 = "c:/temp/test.grd";
    SpatRaster r(filename1);
    show(r);
    r = r.arith(r, "*");
    show(r);
    string f = "c:/temp/write.grd";
    r.writeRaster(f, true);
    SpatRaster x(f);
    show(x);


    string filename3 = "c:/temp/test3.grd";
    SpatRaster rr(filename3);
    SpatRaster x = rr.addSources(rr);
    x = x.addSources(rr);
    x.setNames(std::vector<string>{"a", "b", "c", "d", "e", "f", "g", "h", "i"});
    show(x);
    std::cout << "--" << endl;
    std::vector<unsigned> lyrs = {1,3,5};
    SpatRaster out = x.subset(lyrs, "", FALSE);
    show(out);


    show(r);
  //  x = x.addSources(x);
  //  x.setNames(std::vector<string>{"a", "b", "c", "d"});
    show(x);

    std::cout << "--" << endl;
    std::vector<unsigned> lyrs = {0,5,6};
    std::vector<unsigned> srcs = x.sourcesFromLyrs(lyrs);
    for (size_t i=0; i<lyrs.size(); i++) {
        std::cout << "source: " << srcs[i] << endl;
    }

    SpatRaster out = x.subset({1});
    show(out);

    for (size_t i=0; i<x.nlyr(); i++) {
        std::cout << "source: " << x.sourceFromLyr(i) << endl;
    }
    std::cout << endl;
    std::vector<unsigned> nl = x.nlyrBySource();
    for (size_t i=0; i<nl.size(); i++) {
        std::cout << "layers: " << nl[i] << endl;
    }

    std::cout << endl;
    nl = x.lyrsBySource();
    for (size_t i=0; i<nl.size(); i++) {
        std::cout << "source: " << nl[i] << endl;
    }


    std::vector<unsigned> sub = {1,2};
  //  SpatRaster y = x.subset(sub);
  //  show(y);
    SpatGeomRing g;
    std::vector<double> x = {-180,-140,10,-140};
    std::vector<double> y = {-20,55,0,-60};
    g.set(x, y);
    SpatGeomRings gg;
    gg.addGeom(g);
    SpatPolygons p;
    p.addGeometry(gg);
    SpatExtent e;
    SpatRaster r(90, 180, 1, e, "");
   	SpatRaster out = r.geometry();
	out.writeStart("", true);
  	show(out);
    std::cout << "bs n        : " << out.bs.n << endl;
    std::cout << "bs row[0]   : " << out.bs.row[0] << endl;
    std::cout << "bs nrows[0] : " << out.bs.nrows[0] << endl;
    std::cout << endl;
    SpatRaster rr = r.rasterizePolygons(p, 0, "", true);
    show(rr);

    string filename = "c:/temp/test.grd";
    SpatRaster r(filename);
    r.setNames(std::vector<string>{"a"});
    show(r);
//    SpatRaster x(filename);
    SpatRaster z = r.addSource(r);
    z.setNames(std::vector<string>{"a", "b"});
    show(z);
    z.source.insert(z.source.end(), z.source.begin(), z.source.end());
    show(z);



   if (r.source[0].memory) {std::cout << "mem : true" << endl;} else {std::cout << "mem : false" << endl;}
    if (x.source[0].memory) {std::cout << "mem : true" << endl;} else {std::cout << "mem : false" << endl;}
    show(r);
    show(x);
    SpatRaster m = r.getCopy();
    show(m);
*/
//    system("pause");

