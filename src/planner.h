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

enum class Intent {
  HOLD,
  PREP_LC_L,
  PREP_LC_R,
  LC_L,
  LC_R,
};
std::ostream& operator << (std::ostream& os, const vector<double>& v) {
  for (int i = 0; i < v.size(); i++) {
    os << v[i];
    if (i < v.size() - 1) {
      os  << ", ";
    }
  }
  return os;
}

std::ostream& operator << (std::ostream& os, const Intent& i) {
  string output = "unknown";
  switch(i) {
    case Intent::HOLD:
      output = "HOLD";
      break;

    case Intent::PREP_LC_L:
      output = "PREP_LC_L";
      break;

    case Intent::PREP_LC_R:
      output = "PREP_LC_R";
      break;

    case Intent::LC_L:
      output = "LC_L";
      break;

    case Intent::LC_R:
      output = "LC_R";
      break;
  }

  os << output;
  return os;
}

struct Plan {
  vector<double> x;
  vector<double> y;
  double cost;
};

std::ostream& operator << (std::ostream& os, const Plan& plan) {
  os << "cost(" << plan.cost << ")" << endl << "points: " << endl;
  double prev_x;
  double prev_y;
  double prev_speed;
  double prev_accel;
  for (int i = 0; i < plan.x.size(); i++) {
    double x = plan.x[i];
    double y = plan.y[i];
    os << "  (" << x << ", " << y << ")";
    double speed = 0;
    double accel = 0;
    if (i > 0) {
      double dx = x - prev_x;
      double dy = y - prev_y;
      speed = sqrt(dx*dx + dy*dy) / 0.02;
      os << ", " << speed << " m/s";
    }
    if (i > 1) {
      accel = (speed - prev_speed) / 0.02;
      os << ", " << accel << " m/s/s";
    }
    if (i > 2) {
      double jerk = (accel - prev_accel) / 0.02;
      os << ", " << jerk << " m/s/s/s";
    }
    //if (i < plan.x.size() - 1) {
      os << endl;
    //}
    prev_x = x;
    prev_y = y;
    prev_speed = speed;
    prev_accel = accel;
  }
  return os;
}

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

  Intent prev_intent = Intent::HOLD;

  double distance(double x1, double y1, double x2, double y2) {
    return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
  }

  double speed(double x1, double y1, double x2, double y2) {
    return distance(x1, y1, x2, y2) / DT;
  }

  double accel(double x1, double y1, double x2, double y2, double x3, double y3) {
    //cout << "1: (" << x1 << ", " << y1 << ") ";
    //cout << "2: (" << x2 << ", " << y2 << ") ";
    //cout << "3: (" << x3 << ", " << y2 << ") ";
    const double speed1 = speed(x1, y1, x2, y2);
    const double speed2 = speed(x2, y2, x3, y3);
    //cout << "s1: " << speed1 << ", s2: " << speed2 << endl;;
    const double lon_accel = (speed2 - speed1) / DT;
    //cout << "lon_accel: " << lon_accel << endl;
    const double heading1 = atan2(y2-y1, x2-x1);
    const double heading2 = atan2(y3-y2, x3-x2);
    const double yaw_dot = (heading2 - heading1) / DT;
    //cout << "yaw_dot: " << yaw_dot << endl;
    const double lat_accel = (speed2 + speed1) / 2.0 * yaw_dot;
    return sqrt(lon_accel*lon_accel + lat_accel*lat_accel);
  }


  Plan get_hold_lane_plan(const Pose& pose, double s_end, double end_speed, double d, double T) {
    //cout << "hold_lane_plan" << endl;
    double s_dot = sqrt(pose.x_dot*pose.x_dot + pose.y_dot*pose.y_dot);
    double s_ddot = sqrt(pose.x_ddot*pose.x_ddot + pose.y_ddot*pose.y_ddot);
    vector<double> start_s_v = {pose.s, s_dot, s_ddot};
    vector<double> end_s_v = {s_end, end_speed, 0.0};

    vector<double> a_s = JMT(start_s_v, end_s_v, T);

    Plan p;
    for (double t = DT; t < T; t += DT) {
      double s = eval_poly(a_s, t);
      Pose _pose = map.get_cartesian(s, d);
      p.x.push_back(_pose.x);
      p.y.push_back(_pose.y);
    }

    return p;
  }

  Plan get_plan(const Pose& pose, double s_end, double d_end, double end_speed, double T) {
    //cout << "s_end: " << s_end << ", d_end: " << d_end << endl;
    Pose end_pose = map.get_cartesian(s_end, d_end);
    Pose end_pose1 = map.get_cartesian(s_end+end_speed*DT, d_end);
    double end_x_dot = (end_pose1.x - end_pose.x) / DT;
    double end_y_dot = (end_pose1.y - end_pose.y) / DT;

    vector<double> start_x_v = {pose.x, pose.x_dot, pose.x_ddot};
    vector<double> end_x_v = {end_pose.x, end_x_dot, 0.0};

    vector<double> start_y_v = {pose.y, pose.y_dot, pose.y_ddot};
    vector<double> end_y_v = {end_pose.y, end_y_dot, 0.0};

    //cout << "start_x_v: " << start_x_v << endl;
    //cout << "end_x_v: " << end_x_v << endl;

    //cout << "start_y_v: " << start_y_v << endl;
    //cout << "end_y_v: " << end_y_v << endl;

    vector<double> a_x = JMT(start_x_v, end_x_v, T);
    vector<double> a_y = JMT(start_y_v, end_y_v, T);

    Plan p;
    for (double t = DT; t < T; t += DT) {
      p.x.push_back(eval_poly(a_x, t));
      p.y.push_back(eval_poly(a_y, t));
    }

    return p;
  }

  double plan_cost(const Plan& plan, double start_s) {
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
      //cout << "a: " << a << endl;
      accels.push_back(a);
      if (a > ACCEL_MAX) {
        accel_cost += pow(a, 8);
      } else {
        accel_cost += pow(a, 4) / 100.0;
      }
    }

    const double JERK_MAX = 20;
    double jerk_cost = 0;
    for (int i = 0; i < accels.size() - 1; i++) {
      const double jerk = abs(accels[i+1] - accels[i])/ DT;
      //cout << "jerk: " << jerk << endl;
      if (jerk > JERK_MAX ) {
        jerk_cost += pow(jerk, 8);
      } else {
        jerk_cost += jerk;
      }
    }

    // CTE cost
    double CTE_cost = 0;
   /* vector<double> CTEs(count);
    double s = start_s;
    for (int i = 0; i < count; i++) {
      double d = map.get_d_from_x_y_s(plan.x[i], plan.y[i], s);
      int lane = max(0, min(2, (int)d/4));
      double lane_d = lane*4.0;
      double CTE = abs(lane_d - d);
      CTE_cost += CTE/100.0;
      if (i < count-1) {
        double dx = plan.x[i+1] - plan.x[i];
        double dy = plan.y[i+1] - plan.y[i];
        s += sqrt(dx*dx + dy*dy);
      }
    }*/

    speed_cost /= count;
    CTE_cost /= count;
    accel_cost /= count - 1;
    jerk_cost /= count - 2;
    double total_cost = speed_cost + CTE_cost + accel_cost + jerk_cost;
    cout << "cost(" << total_cost << "): " << speed_cost << ", " <<
                                              CTE_cost << ", " <<
                                              accel_cost << ", " <<
                                              jerk_cost << endl;
    return total_cost;
  }

  double distance_ahead_in_lane(const double s, const int lane,
      const vector<PredictedObject>& predictions, int time_step) {
    double distance_ahead = 100000;
    //const int lane = pose.d/4;
    for (const PredictedObject& prediction : predictions) {
      Pose car_pose = prediction.poses[time_step];
      const int car_lane = car_pose.d / 4;
      if (car_lane == lane) {
        double ahead = car_pose.s - s;
        if (ahead > 0 && ahead < distance_ahead) {
          distance_ahead = ahead;
        }
      }
    }

    return distance_ahead;
  }

  Intent get_intent(const Pose& pose, const vector<PredictedObject>& predictions) {
    const double buffer = 25;
    const int lane = pose.d / 4;
    Intent intent = prev_intent;
    // state machine for intent
    if (prev_intent == Intent::HOLD) {
      // calculate free time ahead
      double distance_ahead[3] = {
        distance_ahead_in_lane(pose.s, 0, predictions, 0),
        distance_ahead_in_lane(pose.s, 1, predictions, 0),
        distance_ahead_in_lane(pose.s, 2, predictions, 0),
      };
      if (distance_ahead[lane] < buffer) {
        // PREP_LC
        if (lane == 0 && distance_ahead[1] > buffer) {
          intent = Intent::PREP_LC_R;
        } else if (lane == 1) {
          if (distance_ahead[0] > buffer) {
            intent = Intent::PREP_LC_L;
          } else if (distance_ahead[2] > buffer) {
            intent = Intent::PREP_LC_R;
          }
        } else if (lane == 2 && distance_ahead[1] > buffer) {
          intent = Intent::PREP_LC_L;
        }
      }
    }
    return intent;
  }

 public:
  Planner(const Map& _map, const double DT) : map(_map) {
    this->DT = DT;
  }

  // need to know:
  //  x, x', x''
  //  y, y', y''
  Plan get_plan(const Pose& pose, const vector<PredictedObject>& predictions) {
    const int lane = pose.d / 4;
    Intent intent = get_intent(pose, predictions);

    cout << "intent: " << intent << endl;
    Plan plan;
    const double car_speed = sqrt(pose.x_dot*pose.x_dot + pose.y_dot*pose.y_dot);
    if (true || intent == Intent::HOLD) {
      // short JMT to lane center, T is a function of distance from center

      const double min_speed = max(0.0, car_speed - 1.0);
      const double max_speed = min(20.0, car_speed + 1.0);
      const double d_speed = 0.2;
      double min_cost = -1;

      for (double T = 1.0; T<= 1.05; T += 0.1) {
        for (double end_speed = min_speed; end_speed <= max_speed; end_speed += d_speed) {
          double s_end = pose.s + end_speed * T;
          double d_end = 6; // TODO: actual lane!
          cout << "T: " << T << ", car_speed: " << car_speed <<
                                ", end_speed: " << end_speed <<
                                ", s: " << pose.s <<
                                ", s_end: " << s_end <<
                                ": pose: " << pose <<
                                endl;
          Plan p = get_plan(pose, s_end, d_end, end_speed, T);

          //Plan p = get_hold_lane_plan(pose, s_end, end_speed, d_end, T);
          p.cost = plan_cost(p, pose.s);
          //cout << p << endl << endl;
          if (min_cost < 0 || p.cost < min_cost) {
            min_cost = p.cost;
            plan = p;
          }
        }
      }
    } else if (intent == Intent::PREP_LC_R) {
    } else if (intent == Intent::PREP_LC_L) {

    }

    return plan;
  }
};

