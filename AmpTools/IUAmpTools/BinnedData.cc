#include "IUAmpTools/BinnedData.h"
#include <stdio.h>
#include <iostream>

void BinnedData::addPoint(double x, double y, double exl, double exh, double eyl, double eyh){
  datapoint newpoint;
  newpoint.x = x;
  newpoint.y = y;
  newpoint.exl = exl;
  newpoint.exh = exh;
  newpoint.eyl = eyl;
  newpoint.eyh = eyh;
  m_thedata.push_back(newpoint);
}

void BinnedData::addPoint(double x, double y, double ey){
  datapoint newpoint;
  newpoint.x = x;
  newpoint.y = y;
  newpoint.exl = 0;
  newpoint.exh = 0;
  newpoint.eyl = ey;
  newpoint.eyh = ey;  
  m_thedata.push_back(newpoint);
}

void BinnedData::addPoint(double x, double y, double ex, double ey){
  datapoint newpoint;
  newpoint.x = x;
  newpoint.y = y;
  newpoint.exl = ex;
  newpoint.exh = ex;
  newpoint.eyl = ey;
  newpoint.eyh = ey;  
  m_thedata.push_back(newpoint);
}

void BinnedData::addPoint(double x, double y, double ex, double eyl, double eyh){
  datapoint newpoint;
  newpoint.x = x;
  newpoint.y = y;
  newpoint.exl = ex;
  newpoint.exh = ex;
  newpoint.eyl = eyl;
  newpoint.eyh = eyh;  
  m_thedata.push_back(newpoint);
}

datapoint BinnedData::getPoint(unsigned int i){
  datapoint currentPoint;
  if(i<m_thedata.size())
    currentPoint = m_thedata.at(i);
  else{
    cerr << "****** ERROR: requested point in binned data is out of range!! ******" << endl;
    assert(0);
  }
  return currentPoint;
}
