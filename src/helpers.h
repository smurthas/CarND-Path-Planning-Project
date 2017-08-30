#ifndef __HELPERS__
#define __HELPERS__


#include <fstream>
#include <math.h>
//#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>

#include "map.h"
#include <string>
#include <sstream>

using namespace std;

double distance(double x1, double y1, double x2, double y2) {
  return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

// minimum distance from a line segment represented by [(vx, vy), (wx, wy)] to a
// point (px, py)
double min_distance(double vx, double vy, double wx, double wy, double px, double py) {
  const double l2 = pow(distance(vx, vy, wx, wy), 2);
  if (l2 < 0.0001) {
    return distance(vx, vy, px, py);
  }

  double t = ((px - vx) * (wx - vx) + (py - vy) * (wy - vy)) / l2;
  t = max(0.0, min(1.0, t));
  double qx = vx + t * (wx - vx);
  double qy = vy + t * (wy - vy);
  return distance(px, py, qx, qy);
}

Map load_map_from_file(string map_file_) {
  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
    istringstream iss(line);
    double x;
    double y;
    float s;
    float d_x;
    float d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    map_waypoints_x.push_back(x);
    map_waypoints_y.push_back(y);
    map_waypoints_s.push_back(s);
    map_waypoints_dx.push_back(d_x);
    map_waypoints_dy.push_back(d_y);
  }
  Map map(map_waypoints_x, map_waypoints_y, map_waypoints_s);
  return map;
}

std::ostream& operator << (std::ostream& os, const vector<double>& v) {
  os << "[ ";
  for (int i = 0; i < v.size(); i++) {
    os << v[i];
    if (i < v.size() - 1) {
      os  << ", ";
    }
  }
  os << " ]";
  return os;
}


#endif
