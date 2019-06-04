/*

#include <geos/geom/Polygon.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/CoordinateSequenceFactory.h>
#include <geos/geom/GeometryFactory.h>
#include <iostream>

geos::geom::Polygon* MakeBox(double xmin, double ymin, double xmax, double ymax){
  geos::geom::GeometryFactory factory;
  geos::geom::CoordinateSequence* temp = factory.getCoordinateSequenceFactory()->create((std::size_t) 0, 0);

  temp->add(geos::geom::Coordinate(xmin, ymin));
  temp->add(geos::geom::Coordinate(xmin, ymax));
  temp->add(geos::geom::Coordinate(xmax, ymax));
  temp->add(geos::geom::Coordinate(xmax, ymin));
  //Must close the linear ring or we will get an error:
  //"Points of LinearRing do not form a closed linestring"
  temp->add(geos::geom::Coordinate(xmin, ymin));

  geos::geom::LinearRing *shell=factory.createLinearRing(temp);

  //NULL in this case could instead be a collection of one or more holes
  //in the interior of the polygon
  return factory.createPolygon(shell,NULL);
}

*/
/* 
int main(){
  geos::geom::Polygon* box = MakeBox(0,0,10,10);
  std::cout<<box->getArea()<<std::endl;
}
*/
