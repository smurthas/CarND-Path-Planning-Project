#include <vector>
#include <math.h>
#include <algorithm>

#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/LU"

#include "prediction.h"
#include "map.h"
#include "pose.h"

using namespace std;
//using namespace Eigen;

struct Plan {
  vector<double> x;
  vector<double> y;
  double cost;
};


double eval_poly(vector<double> coeffs, double x) {
  double y = 0;
  for (int i = 0; i < coeffs.size(); i++) {
    y += coeffs[i] * pow(x, i);
  }

  return y;
}

vector<double> JMT(vector< double> start, vector <double> end, double T) {
  double s_i = start[0];
  double s_dot_i = start[1];
  double s_ddot_i = start[2];
  double s_f = end[0];
  double s_dot_f = end[1];
  double s_ddot_f = end[2];
  Eigen::VectorXd esses(3);
  esses << (s_f - (s_i + s_dot_i*T + 0.5*s_ddot_i*T*T)),
           (s_dot_f - (s_dot_i + s_ddot_i*T)),
           (s_ddot_f - s_ddot_i);

  Eigen::MatrixXd Ts(3, 3);
  double T_2 = T*T;
  double T_3 = T_2*T;
  double T_4 = T_3*T;
  double T_5 = T_4*T;
  Ts << T_3,    T_4,    T_5,
      3*T_2,  4*T_3,  5*T_4,
        6*T, 12*T_2, 20*T_3;
  Eigen::VectorXd As = Ts.inverse() * esses;

  return {s_i, s_dot_i, 0.5*s_ddot_i, As[0], As[1], As[2]};
}


class Planner {
  Map map;
  double DT;

  double distance(double x1, double y1, double x2, double y2) {
    return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
  }

  double speed(double x1, double y1, double x2, double y2) {
    return distance(x1, y1, x2, y2) / DT;
  }

  double accel(double x1, double y1, double x2, double y2, double x3, double y3) {
    const double speed1 = speed(x1, y1, x2, y2);
    const double speed2 = speed(x2, y2, x3, y3);
    const double lon_accel = (speed2 - speed1) / DT;
    const double heading1 = atan2(y2-y1, x2-x1);
    const double heading2 = atan2(y3-y2, x3-x2);
    const double yaw_dot = (heading2 - heading1) / DT;
    const double lat_accel = (speed2 + speed1) / 2.0 * yaw_dot;
    return sqrt(lon_accel*lon_accel + lat_accel*lat_accel);
  }


  Plan get_plan(const Pose& pose, double s_end, double d_end, double end_speed, double T) {
    Pose end_pose = map.get_cartesian(s_end, d_end);
    Pose end_pose1 = map.get_cartesian(s_end+end_speed*DT, d_end);
    double end_x_dot = (end_pose1.x - end_pose.x) / DT;
    double end_y_dot = (end_pose1.y - end_pose.y) / DT;

    vector<double> start_x_v = {pose.x, pose.x_dot, pose.x_ddot};
    vector<double> end_x_v = {end_pose.x, end_x_dot, 0.0};

    vector<double> start_y_v = {pose.y, pose.y_dot, pose.y_ddot};
    vector<double> end_y_v = {end_pose.y, end_y_dot, 0.0};


    vector<double> a_x = JMT(start_x_v, end_x_v, T);
    vector<double> a_y = JMT(start_y_v, end_y_v, T);

    Plan p;
    for (double t = DT; t < T; t += DT) {
      p.x.push_back(eval_poly(a_x, t));
      p.y.push_back(eval_poly(a_y, t));
    }

    return p;
  }

  double plan_cost(const Plan& plan) {
    const double count = plan.x.size();
    const double SPEED_LIMIT = 21;
    double speed_cost = 0;
    for (int i = 1; i < count; i++) {
      const double spd = speed(plan.x[i], plan.y[i], plan.x[i-1], plan.y[i-1]);
      //cout << "spd: " << spd << endl;
      if (spd > SPEED_LIMIT) {
        speed_cost += pow(spd - SPEED_LIMIT, 4);
      } else {
        speed_cost += pow(SPEED_LIMIT - spd, 2);
      }
    }

    const double ACCEL_MAX = 10;
    double accel_cost = 0;
    vector<double> accels;
    for (int i = 2; i < count; i++) {
      const double a = accel(plan.x[i-2], plan.y[i-2], plan.x[i-1], plan.y[i-1], plan.x[i], plan.y[i]);
      accels.push_back(a);
      //cout << "accel: " << a << endl;
      const double a_sq = a*a;
      if (a_sq > ACCEL_MAX*ACCEL_MAX) {
        accel_cost += pow(a, 4);
      } else {
        accel_cost += pow(a, 4) / 100.0;
      }
    }

    const double JERK_MAX = 10;
    double jerk_cost = 0;
    for (int i = 1; i < accels.size(); i++) {
      const double jerk = (accels[i] - accels[i-1])/ DT;
      //cout << "jerk: " << jerk << endl;
      //const double j_sq = jerk*jerk;
      /*if (j_sq > JERK_MAX*JERK_MAX ) {
        jerk_cost += j_sq * j_sq;
      } else {*/
        jerk_cost += abs(jerk);
      //}
    }
    speed_cost /= count;
    accel_cost /= count - 1;
    jerk_cost /= count - 2;
    double total_cost = speed_cost + accel_cost + jerk_cost;

    cout << "cost(" << total_cost << "): " << speed_cost << ", " << accel_cost << ", " << jerk_cost << endl;
    return total_cost;
  }


 public:
  Planner(const Map& _map, const double DT) : map(_map) {
    this->DT = DT;
  }

  // need to know:
  //  x, x', x''
  //  y, y', y''
  Plan get_plan(const Pose& pose, const vector<PredictedObject>& predictions) {
    double car_speed = sqrt(pose.x_dot*pose.x_dot + pose.y_dot*pose.y_dot);
    cout << "planning, speed: " << car_speed <<
            ", pose: (" << pose.x << ", " << pose.y << "), (" << pose.s << ")" <<
            ", velo: (" << pose.x_dot << ", " << pose.y_dot << ")" <<
            ", accl: (" << pose.x_ddot << ", " << pose.y_ddot << ")" << endl;

    double min_cost = 100000000000;
    Plan best_plan;
    const double T = 1;
    for (double end_speed = max(0.0, car_speed - 5.0); end_speed <= min(20.0, car_speed+10); end_speed++) {
      cout << "end_speed: " << end_speed << ", ";
      double s_end = pose.s + (car_speed + end_speed) / 2.0 * T;
      double d_end = 6; // TODO: actual lane!

      Plan p = get_plan(pose, s_end, d_end, end_speed, T);
      p.cost = plan_cost(p);
      //Plan p = get_plan(x, y, s_last, s_end, d, end_speed, T, map);
      //double cost = plan_cost(p.x, p.y, car_speed);
      if (p.cost < min_cost) {
        min_cost = p.cost;
        best_plan = p;
      }
    }

    return best_plan;
  }
};

