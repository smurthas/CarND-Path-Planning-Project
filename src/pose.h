#ifndef __POSE__
#define __POSE__

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


#endif

