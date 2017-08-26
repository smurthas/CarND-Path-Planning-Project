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

#endif

