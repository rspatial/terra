// Copyright (c) 2018-2019  Robert J. Hijmans
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

#ifndef SPATDATE_GUARD
#define SPATDATE_GUARD

#include <vector>
//#include <string>

class SpatDate {
  public:
    long v;
    SpatDate() {v = 0;};
    SpatDate(const int& year, const int& month, const int& day);
    SpatDate(int doy, int year);
    SpatDate(const long& d) { v = d; };
//    SpatDate(const std::string& s);

    std::vector<int> ymd();
    int doy();

    bool is_leap_year();

    SpatDate operator ++(); // prefix
    SpatDate operator ++(int); // postfix
    SpatDate operator --(); // prefix
    SpatDate operator --(int); // postfix
    void operator = (const SpatDate&);
};

SpatDate operator + (SpatDate&, int);
SpatDate operator + (int, SpatDate&);
SpatDate operator - (SpatDate&, int);
int  operator - (SpatDate&, SpatDate&);
bool operator == (const SpatDate&, const SpatDate&); 
bool operator != (const SpatDate&, const SpatDate&);
bool operator < (const SpatDate&, const SpatDate&);
bool operator > (const SpatDate&, const SpatDate&);
bool operator <= (const SpatDate&, const SpatDate&);
bool operator >= (const SpatDate&, const SpatDate&);

bool isleapyear(const int& year) {
  return (year % 4 == 0) && ((year % 400 == 0) || (year % 100 != 0 ));
}

std::vector<int> month_days(int year) {
  if (isleapyear(year)) {
    std::vector<int> md {31,29,31,30,31,30,31,31,30,31,30,31};
    return md;
  } else { 
    std::vector<int> md {31,28,31,30,31,30,31,31,30,31,30,31};
    return md;
  }
}

std::vector<int> SpatDate::ymd() {
  int year, month, day;
  long x = v;
  
  if (x >= 0) {
    long days = 0;
    year = 1969;
    while (days <= x) {
      year++;
      days += (365 + isleapyear(year));
    }
    x = x - days + (365 + isleapyear(year));
    std::vector<int> mdays = month_days(year);
    days = -1;
    month = -1;
    while (days < x) {
      month++;
      days += mdays[month];
    }
    day = x - days + mdays[month];
    month++;
  } else { // (x < 0)
    long days = 0;
    year = 1970;
    while (days > x) {
      year--;
      days -= (365 + isleapyear(year));
    }
    x = x - (days + (365 + isleapyear(year)));
    x = abs(x);
    std::vector<int> mdays = month_days(year);
    days = 0;
    month = 12;
    while (days < x) {
      month--;
      days += mdays[month];
    }
    day = days - x +1;
    month++;
  }
  std::vector<int> result {year, month, day};
  return result;
}


int date_from_ymd(std::vector<int> ymd) {

  int year = 1970;
  long day = -1;
  
  if (ymd[0] > 1970) {
    
    for (int i=0; i<(ymd[0]-1970); i++) {
      day += 365 + isleapyear(year);
      year++;
    }
    std::vector<int> mdays = month_days(year);
    for (int i=0; i<(ymd[1]-1); i++) {
      day += mdays[i];
    }
    day += ymd[2];
    
  } else {

    day = -1;
    for (int i=0; i<(1970-ymd[0]); i++) {
      day -= 365 + isleapyear(year);
      year--;
    }
    std::vector<int> mdays = month_days(year);
    for (int i=0; i<(ymd[1]-1); i++) {
      day += mdays[i];
    }
    day += ymd[2];

  }
  return(day);
}
  
 
bool SpatDate::is_leap_year() {
	std::vector<int> d = ymd();
	return isleapyear( d[0] );
};


SpatDate::SpatDate(const int& year, const int& month, const int& day) {
	std::vector<int> ymd {year, month, day};
	v = date_from_ymd(ymd);
}

SpatDate::SpatDate(int doy, int year) {
	if (doy == 0) doy = 1;
	
	if (doy < 0) {  // this is a weird case, but let's try to handle it
		while ( doy < 0 ) {
			year--;
			doy = doy + 365 + isleapyear(year);
		}
		doy = 365 + isleapyear(year) - doy;

	} else {
		while ( doy > 366 ) {
			doy = doy - 365 - isleapyear(year);
			year++;
		}
	}
	std::vector<int> mdays = month_days(year);
	int month= 0 ;
	while ( doy > mdays[month] ) {
		doy = doy - mdays[month];
		month++;
	}
	SpatDate(year, month+1, doy);
}


int SpatDate::doy() {
	std::vector<int> d = ymd();
  std::vector<int> mdays = month_days(d[0]);

	int doy = 0;
	for (int i=0; i < d[1]-1; i++) {
		doy = doy + mdays[i];
	}
	return( doy + d[2] );
};


bool operator == (const SpatDate& d1,const SpatDate& d2){
	return (d1.v == d2.v);
}

bool operator !=(const SpatDate& d1, const SpatDate& d2){
	return (d1.v != d2.v);
}

inline SpatDate next_date(const SpatDate& d) {
	SpatDate x(d.v+1);
	return x;
};

inline SpatDate previous_date(const SpatDate& d){
	SpatDate x(d.v-1);
	return x;
};

SpatDate SpatDate::operator ++(int){
	SpatDate d = *this;
	*this = next_date(d);
	return d;
}

SpatDate SpatDate::operator ++() {
	*this = next_date(*this);
	return *this;
}

SpatDate SpatDate::operator --(int){ 
	SpatDate d = *this;
	*this = previous_date(d);
	return d;
}

SpatDate SpatDate::operator --(){
	*this = previous_date(*this);
	return *this;
}

bool operator < (const SpatDate& d1, const SpatDate& d2){
	return (d1.v < d2.v);
};

bool operator <=(const SpatDate& d1, const SpatDate& d2){
	return (d1.v <= d2.v);
}

bool operator >=(const SpatDate& d1, const SpatDate& d2) {
	return (d1.v >= d2.v);
}

bool operator > (const SpatDate& d1, const SpatDate& d2) {
	return (d1.v > d2.v);
}

int operator - (SpatDate& d1,  SpatDate& d2) {
	return (d1.v - d2.v);
}

SpatDate operator +(SpatDate &d1, int x){
	return SpatDate(d1.v + x);
}

SpatDate operator +(int x, SpatDate &d1) {
	return SpatDate(d1.v + x);
}

SpatDate operator -(SpatDate &d1, int x){
	return SpatDate(d1.v - x);
}

void SpatDate::operator =(const SpatDate &d){
	this->v = d.v;
}

#endif
