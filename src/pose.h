#ifndef __POSE__
#define __POSE__
#include <math.h>

struct Pose {
  double x;
  double y;
  double s;
  double d;

  double y_dot;
  double x_dot;

  double y_ddot;
  double x_ddot;

  double theta;
};

std::ostream& operator << (std::ostream& os, const Pose& p) {
  os << "{ x: " << p.x << ", y: " << p.y << ", s: " << p.s << ", d: " << p.d
     << ", x_dot: " << p.x_dot << ", y_dot: " << p.y_dot
     << ", x_ddot: " << p.x_ddot << ", y_ddot: " << p.y_ddot
     << ", theta: " << p.theta << " }";
  return os;
}

double distance(const Pose& p1, const Pose& p2) {
  double dx = p2.x - p1.x;
  double dy = p2.y - p1.y;
  return sqrt(dx*dx + dy*dy);
}


#endif

