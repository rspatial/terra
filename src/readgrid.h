/*
  C++ Class to read gridfiles
  Robert Hijmans
  January 2008
  r.hijmans@gmail.com
*/

using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>

#include "SimpleIni.h"
#include <string>


class CGrid {
    string grifilename, grdfilename; 
    double xmin, xmax, ymin, ymax, xres, yres;

    void getshortrow(int, int[]);
    void getlongrow(int, int[]);
    void getfloatrow(int, double[]);
    void getdoublerow(int, double[]);
    short getshortcell(int, int);
    long getlongcell(int, int);
    double getdoublecell(int, int);
    float getfloatcell(int, int);

  public:
    string datatype;
    int nrows, ncols;

    CGrid();
    bool setup(string fname);
    bool setfilename(string);

    int getintcell(int, int);
    double getrealcell(int, int);
    double getcell(int, int);
    void getintrow(int, int[]);
    void getrealrow(int, double[]);
    void getrow(int, double[]);
};

CGrid::CGrid()
{
// empty grid
  datatype = "INT2BYTES";
  ncols = 0;
  nrows = 0;
  xmin = 0;
  xmax = 0;
  ymin = 0;
  ymax = 0;
  xres = 0;
  yres = 0;
}


bool CGrid::setfilename(string fname)
{
   grdfilename = fname+".grd";
   grifilename = fname+".gri";
}


bool CGrid::setup(string fname)
{
   setfilename(fname);
   CSimpleIniA ini(TRUE, FALSE, FALSE);

   char ss[grdfilename.length()];
   strcpy(ss, grdfilename.c_str());
   SI_Error rc = ini.LoadFile(ss);
   if (rc < 0) return false;
   else {
     datatype = ini.GetValue("Data", "DataType");
     nrows = atoi(ini.GetValue("GeoReference", "Rows"));
     ncols = atoi(ini.GetValue("GeoReference", "Columns"));
     xmin = atof(ini.GetValue("GeoReference", "MinX"));
     xmax = atof(ini.GetValue("GeoReference", "MaxX"));
     ymin = atof(ini.GetValue("GeoReference", "MinY"));
     ymax = atof(ini.GetValue("GeoReference", "MaxY"));
     xres = (xmax - xmin) / ncols;
     yres = (ymax - ymin) / nrows;

     return true;
   }
}


short CGrid::getshortcell(int rownr, int colnr)
{
  ifstream thefile;
  short value[1];
  const int dsize = 2;
  thefile.open (grifilename.c_str(), ios::in | ios::binary);
  thefile.seekg ( (rownr * ncols + colnr) * dsize, ios::beg);
  thefile.read ((char*)value, sizeof(value)); 
  thefile.close();
  return (int)value[0];
}

long CGrid::getlongcell(int rownr, int colnr)
{
  ifstream thefile;
  long value[1];
  const int dsize = 4;
  thefile.open (grifilename.c_str(), ios::in | ios::binary);
  thefile.seekg ( (rownr * ncols + colnr) * dsize, ios::beg);
  thefile.read ((char*)value, sizeof(value)); 
  thefile.close();
  return (int)value[0];
}

int CGrid::getintcell(int rownr, int colnr)
{
  if (string(datatype) == "INT2BYTES") {return getshortcell(rownr, colnr);}
  else if (string(datatype) == "INT4BYTES") {return getlongcell(rownr, colnr);}
}


float CGrid::getfloatcell(int rownr, int colnr)
{
  ifstream thefile;
  float value[1];
  const int dsize = 4;
  thefile.open (grifilename.c_str(), ios::in | ios::binary);
  thefile.seekg ( (rownr * ncols + colnr) * dsize, ios::beg);
  thefile.read ((char*)value, sizeof(value)); 
  thefile.close();
  return (double)value[0];
}


double CGrid::getdoublecell(int rownr, int colnr)
{
  ifstream thefile;
  double value[1];
  const int dsize = 8;
  thefile.open (grifilename.c_str(), ios::in | ios::binary);
  thefile.seekg ( (rownr * ncols + colnr) * dsize, ios::beg);
  thefile.read ((char*)value, sizeof(value)); 
  thefile.close();
  return (double)value[0];
}


double CGrid::getrealcell(int rownr, int colnr)
{
  if (string(datatype) == "FLT4BYTES") {return getfloatcell(rownr, colnr);}
  else if (string(datatype) == "FLT8BYTES") {return getdoublecell(rownr, colnr);}
}


double CGrid::getcell(int rownr, int colnr)
{
  if (string(datatype, 0, 3) == "INT") { return getintcell(rownr, colnr); }
  else {return getrealcell(rownr, colnr);}
}


void CGrid::getshortrow(int rownr, int vals[])
{
  ifstream thefile;
  int c;
  const int dsize = 2;
  short numbs[ncols]; 

  thefile.open (grifilename.c_str(), ios::in | ios::binary);
  thefile.seekg (rownr * ncols * dsize, ios::beg);
  thefile.read ((char*)numbs, sizeof(numbs));
  for(c=0; c<= ncols-1; c++) { vals[c] = (int)numbs[c]; }
  thefile.close();
}


void CGrid::getlongrow(int rownr, int vals[])
{
  ifstream thefile;
  int c;
  const int dsize = 4;
  long numbs[ncols]; 

  thefile.open (grifilename.c_str(), ios::in | ios::binary);
  thefile.seekg (rownr * ncols * dsize, ios::beg);
  thefile.read ((char*)numbs, sizeof(numbs));
  for(c=0; c<= ncols-1; c++) { vals[c] = (int)numbs[c]; }
  thefile.close();
}


void CGrid::getfloatrow(int rownr, double vals[])
{
  ifstream thefile;
  int c;
  const int dsize = 4;
  float numbs[ncols]; 

  thefile.open (grifilename.c_str(), ios::in | ios::binary);
  thefile.seekg (rownr * ncols * dsize, ios::beg);
  thefile.read ((char*)numbs, sizeof(numbs));
  for(c=0; c<= ncols-1; c++) { vals[c] = (double)numbs[c]; }
  thefile.close();
}


void CGrid::getdoublerow(int rownr, double vals[])
{
  ifstream thefile;
  int c;
  const int dsize = 8;
  double numbs[ncols]; 

  thefile.open (grifilename.c_str(), ios::in | ios::binary);
  thefile.seekg (rownr * ncols * dsize, ios::beg);
  thefile.read ((char*)numbs, sizeof(numbs));
  for(c=0; c<= ncols-1; c++) { vals[c] = (double)numbs[c]; }
  thefile.close();
}


void CGrid::getrealrow(int rownr, double vals[])
{
  if (string(datatype) == "FLT4BYTES") {getfloatrow(rownr, vals);}
  else if (string(datatype) == "FLT8BYTES") {getdoublerow(rownr, vals);}
}


void CGrid::getintrow(int rownr, int vals[])
{
  if (string(datatype) == "INT2BYTES") {getshortrow(rownr, vals);}
  else if (string(datatype) == "INT4BYTES") {getlongrow(rownr, vals);}
}


void CGrid::getrow(int rownr, double vals[])
{
  int c;
  if (string(datatype, 0, 3) == "INT") {
      int ivals[ncols];
      getintrow(rownr, ivals);
      for(c=0; c<= ncols-1; c++) { vals[c] = ivals[c]; }
  }
  else {getrealrow(rownr, vals);}
}

