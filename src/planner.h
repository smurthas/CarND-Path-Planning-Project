#include <vector>
#include <math.h>
#include <algorithm>

#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/LU"

#include "helpers.h"
#include "prediction.h"
#include "map.h"
#include "pose.h"

using namespace std;
//using namespace Eigen;

enum class Intent {
  HOLD,
  PREP_LC_L,
  PREP_LC_R,
  LC,
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

    case Intent::LC:
      output = "LC";
      break;
  }

  os << output;
  return os;
}

struct Plan {
  vector<Pose> poses;
  double cost;
};

std::ostream& operator << (std::ostream& os, const Plan& plan) {
  os << "cost(" << plan.cost << ")" << endl << "points: " << endl;
  double prev_x;
  double prev_y;
  double prev_speed;
  double prev_accel;
  for (int i = 0; i < plan.poses.size(); i++) {
    double x = plan.poses[i].x;
    double y = plan.poses[i].y;
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
  int goal_lane = 1;
  Plan prev_plan;


  double speed(const Pose& p1, const Pose& p2) {
    return distance(p1, p2) / DT;
  }

  double accel(const Pose& p1, const Pose& p2, const Pose& p3) {
    const double speed1 = speed(p1, p2);
    const double speed2 = speed(p2, p3);
    const double lon_accel = (speed2 - speed1) / DT;
    const double heading1 = atan2(p2.y-p1.y, p2.x-p1.x);
    const double heading2 = atan2(p3.y-p2.y, p3.x-p2.x);
    const double yaw_dot = (heading2 - heading1) / DT;
    const double lat_accel = (speed2 + speed1) / 2.0 * yaw_dot;
    return sqrt(lon_accel*lon_accel + lat_accel*lat_accel);
  }


  Plan get_hold_lane_plan(const Pose& pose, double goal_speed) {
    double s = pose.s;
    double d = pose.d;
    double ds_step = 0;
    double ds_goal = goal_speed*DT;
    if (prev_plan.poses.size() > 0) {
      s = prev_plan.poses[prev_plan.poses.size()-1].s;
      d = prev_plan.poses[prev_plan.poses.size()-1].d;
    }
    if (prev_plan.poses.size() > 1) {
      ds_step = prev_plan.poses[prev_plan.poses.size()-1].s -
                prev_plan.poses[prev_plan.poses.size()-2].s;
    }
    while (prev_plan.poses.size() < 100) {
      if (ds_step < ds_goal) {
        ds_step += 0.001;
      } else if (ds_step > ds_goal) {
        ds_step -= 0.001;
      }
      s += ds_step;
      Pose p = map.get_cartesian(s, d);
      prev_plan.poses.push_back(p);
    }
    return prev_plan;
  }

  Plan get_lane_change_plan(const Pose& pose, int goal_lane, double goal_speed) {
    double s = pose.s;
    double d = pose.d;
    double ds_step = 0;
    double dd_step = 0;
    double goal_d = goal_lane * 4.0 + 2.0;
    cout << "s: " << s << endl;
    cout << "d: " << d << endl;
    cout << "ds_step: " << ds_step << endl;
    cout << "dd_step: " << dd_step << endl;
    cout << "goal_d: " << goal_d << endl;

    if (prev_plan.poses.size() > 0) {
      s = prev_plan.poses[prev_plan.poses.size()-1].s;
      d = prev_plan.poses[prev_plan.poses.size()-1].d;
    }
    if (prev_plan.poses.size() > 1) {
      dd_step = prev_plan.poses[prev_plan.poses.size()-1].d -
                prev_plan.poses[prev_plan.poses.size()-2].d;
      ds_step = prev_plan.poses[prev_plan.poses.size()-1].s -
                prev_plan.poses[prev_plan.poses.size()-2].s;
    }
    cout << "s: " << s << endl;
    cout << "d: " << d << endl;
    cout << "ds_step: " << ds_step << endl;
    cout << "dd_step: " << dd_step << endl;
    const double T = 4.0;
    vector<double> a_d = JMT({d, 0.0, 0.0}, {goal_d, 0.0, 0.0}, T);
    for (double t = DT; t < T; t += DT) {
      double d_step = eval_poly(a_d, t);
      s += ds_step;
      prev_plan.poses.push_back(map.get_cartesian(s, d_step));
    }

    cout << "prev_plan.poses.size(): " << prev_plan.poses.size() << endl;
    return prev_plan;
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
      Pose pose;
      pose.x = eval_poly(a_x, t);
      pose.y = eval_poly(a_y, t);

      p.poses.push_back(pose);
    }

    return p;
  }

  double plan_cost(const Plan& plan, double start_s, const vector<PredictedObject>& predictions) {
    const double count = plan.poses.size();
    const double SPEED_LIMIT = 21;
    double speed_cost = 0;
    for (int i = 1; i < count; i++) {
      const Pose p1 = plan.poses[i-1];
      const Pose p2 = plan.poses[i];
      const double spd = speed(p2, p1);
      //cout << "spd: " << spd << endl;
      if (spd > SPEED_LIMIT) {
        speed_cost += pow(spd - SPEED_LIMIT, 8);
      } else {
        speed_cost += pow(SPEED_LIMIT - spd, 2);
      }
    }

    //double gap_cost = 0.0;
    //for (int i = 0;  )

    const double ACCEL_MAX = 8;
    double accel_cost = 0;
    vector<double> accels;
    for (int i = 2; i < count; i++) {
      const Pose p1 = plan.poses[i-2];
      const Pose p2 = plan.poses[i-1];
      const Pose p3 = plan.poses[i];
      const double a = accel(p1, p2, p3);
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
    vector<double> CTEs(count);
    double s = start_s;
    for (int i = 0; i < count; i++) {
      double d = map.get_d_from_x_y_s(plan.poses[i].x, plan.poses[i].y, s);
      int lane = max(0, min(2, (int)d/4));
      double lane_d = lane*4.0;
      double CTE = abs(lane_d - d);
      CTE_cost += CTE;
      if (i < count-1) {
        s += distance(plan.poses[i], plan.poses[i+1]);
      }
    }

    speed_cost /= count;
    CTE_cost /= count;
    accel_cost /= count - 1;
    jerk_cost /= count - 2;
    double total_cost = speed_cost + CTE_cost + accel_cost + jerk_cost;
    /*cout << "cost(" << total_cost << "): " << speed_cost << ", " <<
                                              CTE_cost << ", " <<
                                              accel_cost << ", " <<
                                              jerk_cost << endl;*/
    return total_cost;
  }

  double speed_ahead_in_lane(const double s, const int lane,
      const vector<PredictedObject>& predictions, int time_step) {
    double d1 = distance_ahead_in_lane(s, lane, predictions, time_step);
    double d2 = distance_ahead_in_lane(s, lane, predictions, time_step+1);
    return (d2-d1) / DT;
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
    const double buffer = 80;
    const int lane = pose.d / 4;
    Intent intent = prev_intent;
    // state machine for intent
    // TODO: change lane if adjacent gap is bigger than current gap
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
          goal_lane = 1;
        } else if (lane == 1) {
          if (distance_ahead[0] > buffer) {
            intent = Intent::PREP_LC_L;
            goal_lane = 0;
          } else if (distance_ahead[2] > buffer) {
            intent = Intent::PREP_LC_R;
            goal_lane = 2;
          }
        } else if (lane == 2 && distance_ahead[1] > buffer) {
          intent = Intent::PREP_LC_L;
          goal_lane = 1;
        }
      }
    } else if (prev_intent == Intent::PREP_LC_L || prev_intent == Intent::PREP_LC_R) {
      intent = Intent::LC;
    } else if (prev_intent == Intent::LC) {
      if (lane == goal_lane) {
        intent = Intent::HOLD;
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
  Plan get_plan(const Pose& pose, const vector<PredictedObject>& predictions, int points_ahead) {
    const int lane = pose.d / 4;
    Intent intent = get_intent(pose, predictions);
    // slice out old points
    while (prev_plan.poses.size() > points_ahead) {
      prev_plan.poses.erase(prev_plan.poses.begin());
    }

    if (intent != prev_intent) {
      cout << "changed intent: " << prev_intent << "→" << intent << endl;
    }

    Plan plan;
    const double car_speed = sqrt(pose.x_dot*pose.x_dot + pose.y_dot*pose.y_dot);
    if (intent == Intent::HOLD) {
      double goal_speed = 20.0;
      double distance_ahead = distance_ahead_in_lane(pose.s, lane, predictions, 0);
      if (distance_ahead < 40) {
        goal_speed = speed_ahead_in_lane(pose.s, lane, predictions, 0);
      }
      plan = get_hold_lane_plan(pose, goal_speed);
    } else if (intent == Intent::PREP_LC_R || intent == Intent::PREP_LC_L) {
      // TODO: wait for a gap before commiting
      plan = get_lane_change_plan(pose, goal_lane, car_speed);
    } else if (intent == Intent::LC) {
      if (lane == goal_lane) {
        double goal_speed = 20.0;
        plan = get_hold_lane_plan(pose, goal_speed);
      } else {
        plan = prev_plan;
      }
    }

    prev_intent = intent;
    prev_plan = plan;

    return plan;
  }
};

