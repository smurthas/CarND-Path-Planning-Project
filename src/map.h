#ifndef __MAP__
#define __MAP__

#include <math.h>
#include <vector>

#include "spline.h"

#include "pose.h"

using namespace std;

class Map {

  tk::spline spline_x;
  tk::spline spline_y;

 public:
  Map(const Map& other) {
    this->spline_x = other.spline_x;
    this->spline_y = other.spline_y;
  }

  Map(vector<double> waypoints_x, vector<double> waypoints_y,
      vector<double> waypoints_s) {
    spline_x.set_points(waypoints_s, waypoints_x);
    spline_y.set_points(waypoints_s, waypoints_y);
  }

  Pose get_cartesian(double s, double d) {
    Pose p;
    double x0 = spline_x(s);
    double y0 = spline_y(s);
    double ds = 0.001;
    double x1 = spline_x(s+ds);
    double y1 = spline_y(s+ds);
    double d_heading = atan2(y1-y0, x1-x0) - M_PI/2.0;
    p.x = x0 + cos(d_heading)*d;
    p.y = y0 + sin(d_heading)*d;
    return p;
  }
};

#endif
