#ifndef __MAP__
#define __MAP__

#include <math.h>
#include <vector>

#include "spline.h"
#include "helpers.h"

#include "pose.h"

using namespace std;

class Map {
  double max_s = 6945.554;
  vector<double> waypoints_x;
  vector<double> waypoints_y;
  vector<double> waypoints_s;

  tk::spline spline_x;
  tk::spline spline_y;

 public:
  Map(const Map& other) {
    this->spline_x = other.spline_x;
    this->spline_y = other.spline_y;
    this->waypoints_x = other.waypoints_x;
    this->waypoints_y = other.waypoints_y;
    this->waypoints_s = other.waypoints_s;
  }

  Map(vector<double> waypoints_x, vector<double> waypoints_y,
      vector<double> waypoints_s) {
    // wrap around
    waypoints_s.insert(waypoints_s.begin(), waypoints_s[waypoints_s.size()-1]-max_s);
    waypoints_s.insert(waypoints_s.begin(), waypoints_s[waypoints_s.size()-2]-max_s);
    waypoints_s.push_back(max_s);
    waypoints_s.push_back(max_s + waypoints_s[3]);

    waypoints_x.insert(waypoints_x.begin(), waypoints_x[waypoints_x.size()-1]);
    waypoints_x.insert(waypoints_x.begin(), waypoints_x[waypoints_x.size()-2]);
    waypoints_x.push_back(waypoints_x[2]);
    waypoints_x.push_back(waypoints_x[3]);

    waypoints_y.insert(waypoints_y.begin(), waypoints_y[waypoints_y.size()-1]);
    waypoints_y.insert(waypoints_y.begin(), waypoints_y[waypoints_y.size()-2]);
    waypoints_y.push_back(waypoints_y[2]);
    waypoints_y.push_back(waypoints_y[3]);

    this->waypoints_x = waypoints_x;
    this->waypoints_y = waypoints_y;
    this->waypoints_s = waypoints_s;

    spline_x.set_points(waypoints_s, waypoints_x);
    spline_y.set_points(waypoints_s, waypoints_y);
  }

  double get_d_from_x_y_s(double x, double y, double s) {
    double dx = spline_x(s) - x;
    double dy = spline_y(s) - y;
    return sqrt(dx*dx + dy*dy);
  }

  Pose get_cartesian(double s, double d) {
    while (s > max_s) {
      s -= max_s;
    }
    while (s < 0) {
      s += max_s;
    }
    Pose p;
    double x0 = spline_x(s);
    double y0 = spline_y(s);
    double ds = 0.001;
    double x1 = spline_x(s+ds);
    double y1 = spline_y(s+ds);
    double d_heading = atan2(y1-y0, x1-x0) - M_PI/2.0;
    p.x = x0 + cos(d_heading)*d;
    p.y = y0 + sin(d_heading)*d;
    p.s = s;
    p.d = d;
    return p;
  }

  /*Pose get_from_cartesian(double x, double y) {
    Pose p;
    p.x = x;
    p.y = y;
    int prev = 0;
    int min_dist = 1000000;
    for (int i = 0; i < waypoints_x.size(); i++) {
      int next = (i+1) % waypoints_x.size();
      double d = min_distance(waypoints_x[i], waypoints_y[i],
                              waypoints_x[next], waypoints_y[next],
                              x, y);
      if (d < min_dist) {
        prev = i;
        min_dist = d;
      }
    }

    double sl = waypoints_s[prev];
    double sr;
    if (prev == waypoints_x.size() - 1) {
      sr = max_s;
    } else {
      sr = waypoints_s[prev+1];
    }

    Pose l, r, m;
    double dl, dr, dm;
    double sm;
    do {
      l = get_cartesian(prev_s, 0);
      r = get_cartesian(next_x, 0);
      dl = distance(p, l);
      dr = distance(p, r);
      sm = (sl + sr) / 2.0;
      m = get_cartesian(sm, 0);
      dm = distance(p, m);

      if (dm < dr && dl < dr) {
        // discard dr
        dr = dm;
        r = m;
        sr = sm;
      } else if (dm < dl && dr < dl) {
        // discard dl
        dl = dm;
        l = m;
        sl = sm;
      }
    } while (abs(sl - sr) > 0.01);

    p.d = distance(p, m);
    p.s = sm;

    return p;
  }*/
};

#endif

