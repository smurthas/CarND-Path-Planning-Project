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

enum class Intent {
  HOLD,
  PREP_LC,
  LC,
};

std::ostream& operator << (std::ostream& os, const Intent& i) {
  string output = "unknown";
  switch(i) {
    case Intent::HOLD:
      output = "HOLD";
      break;

    case Intent::PREP_LC:
      output = "PREP_LC";
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

  int collides(const Plan& plan, const vector<PredictedObject>& predictions) {
    const double MIN_S_GAP = 10;
    const double MIN_D_GAP = 3;
    for (int i = 0; i < plan.poses.size(); i++) {
      const Pose& ego_pose = plan.poses[i];
      for (const PredictedObject& obj : predictions) {
        const Pose& obj_pose = obj.poses[i];
        double ds = abs(obj_pose.s - ego_pose.s);
        double dd = abs(obj_pose.d - ego_pose.d);
        if (ds < MIN_S_GAP && dd < MIN_D_GAP) {
          return obj.id;
        }
      }
    }

    return -1;
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


  Plan get_hold_lane_plan(const Pose& pose, double goal_speed) {
    double s = pose.s;
    double d = pose.d;
    double ds_step = 0;
    double ds_goal = goal_speed*DT;

    // First, inherit existing plan
    Plan plan;
    for (const Pose& p : prev_plan.poses) {
      plan.poses.push_back(p);
    }

    if (plan.poses.size() > 0) {
      s = plan.poses[plan.poses.size()-1].s;
      d = plan.poses[plan.poses.size()-1].d;
    }
    if (plan.poses.size() > 1) {
      ds_step = plan.poses[plan.poses.size()-1].s -
                plan.poses[plan.poses.size()-2].s;
    }
    while (plan.poses.size() < 100) {
      if (ds_step < ds_goal) {
        ds_step += 0.001;
      } else if (ds_step > ds_goal) {
        ds_step -= 0.001;
      }
      s += ds_step;
      Pose p = map.get_cartesian(s, d);
      plan.poses.push_back(p);
    }
    return plan;
  }

  Plan get_lane_change_plan(const Pose& pose, int goal_lane, double goal_speed) {
    double s = pose.s;
    double d = pose.d;
    double ds_step = 0;
    double dd_step = 0;
    double goal_d = goal_lane * 4.0 + 2.0;
    double goal_ds_step = goal_speed*DT;

    // First, inherit existing plan
    Plan plan;
    for (const Pose& p : prev_plan.poses) {
      plan.poses.push_back(p);
      // only keep 10 prev poses during LC
      if (plan.poses.size() > 10) {
        break;
      }
    }

    if (plan.poses.size() > 0) {
      s = plan.poses[plan.poses.size()-1].s;
      d = plan.poses[plan.poses.size()-1].d;
    }
    if (plan.poses.size() > 1) {
      dd_step = plan.poses[plan.poses.size()-1].d -
                plan.poses[plan.poses.size()-2].d;
      ds_step = plan.poses[plan.poses.size()-1].s -
                plan.poses[plan.poses.size()-2].s;
    }

    const double T = 3.0;
    vector<double> a_d = JMT({d, 0.0, 0.0}, {goal_d, 0.0, 0.0}, T);
    for (double t = DT; t < T; t += DT) {
      double d_step = eval_poly(a_d, t);
      s += ds_step;
      plan.poses.push_back(map.get_cartesian(s, d_step));
      if (ds_step < goal_ds_step) {
        ds_step += 0.001;
      }
    }

    // add five lane center points to help us settle into the lane
    for (int i = 0; i < 5; i++) {
      s += ds_step;
      plan.poses.push_back(map.get_cartesian(s, goal_d));
      if (ds_step < goal_ds_step) {
        ds_step += 0.001;
      }
    }

    return plan;
  }

  Intent get_intent(const Pose& pose, const vector<PredictedObject>& predictions) {
    const double buffer = 80;
    const int lane = pose.d / 4;

    Intent intent = prev_intent;

    // calculate free space ahead
    double distance_ahead[3] = {
      distance_ahead_in_lane(pose.s, 0, predictions, 0),
      distance_ahead_in_lane(pose.s, 1, predictions, 0),
      distance_ahead_in_lane(pose.s, 2, predictions, 0),
    };

    if (prev_intent == Intent::HOLD) {
      // if we are closing in on a car, start looking for other lane options
      if (distance_ahead[lane] < buffer) {
        intent = Intent::PREP_LC;
      }
    } else if (prev_intent == Intent::PREP_LC) {
      // if the gap to the car ahead has increased sufficiently, just return to
      // holding our lane
      if (distance_ahead[lane] > buffer + 5) {
        intent = Intent:: HOLD;
      }
    } else if (prev_intent == Intent::LC) {
      // if we've reached the goal lane, switch back to holding our lane
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

  Plan get_plan(const Pose& pose, const vector<PredictedObject>& predictions, int points_ahead) {
    const int lane = pose.d / 4;

    Intent intent = get_intent(pose, predictions);

    // slice out points that we have already passed
    while (prev_plan.poses.size() > points_ahead) {
      prev_plan.poses.erase(prev_plan.poses.begin());
    }

    if (intent != prev_intent) {
      cout << "changed intent: " << prev_intent << "→" << intent << endl;
    }

    double goal_speed = 20.0; // speed limit

    // calculate the distance ahead in all three lanes
    double distance_ahead[3] = {
      distance_ahead_in_lane(pose.s, 0, predictions, 0),
      distance_ahead_in_lane(pose.s, 1, predictions, 0),
      distance_ahead_in_lane(pose.s, 2, predictions, 0),
    };

    double lane_speed = speed_ahead_in_lane(pose.s, lane, predictions, 0);

    // if we are closing in on another car, slow down to just below its speed
    if (distance_ahead[lane] < 20) {
      goal_speed = lane_speed - 1.0;
    } else if (distance_ahead[lane] < 40) {
      goal_speed = lane_speed - 0.5;
    }


    Plan plan;
    const double car_speed = sqrt(pose.x_dot*pose.x_dot + pose.y_dot*pose.y_dot);

    if (intent == Intent::HOLD || intent == Intent::LC) {
      // if we are mid lane change, then this is just planning our path once we
      // arrive in the lane since we keep the previous plan poses
      plan = get_hold_lane_plan(pose, goal_speed);

    } else if (intent == Intent::PREP_LC) {
      // figure out if we should change lanes, only change lanes when an
      // adjacent lane has a longer open distance ahead than the current lane
      const double switch_distance = distance_ahead[lane] + 5;
      if (lane == 0 && distance_ahead[1] > switch_distance) {
        goal_lane = 1;
      } else if (lane == 1) {
        if (distance_ahead[0] > switch_distance) {
          goal_lane = 0;
        } else if (distance_ahead[2] > switch_distance) {
          goal_lane = 2;
        }
      } else if (lane == 2 && distance_ahead[1] > switch_distance) {
        goal_lane = 1;
      }

      if (goal_lane == lane) {
        // if we determine that the current lane is best, stick with that (but
        // look again next time around)
        plan = get_hold_lane_plan(pose, goal_speed);
      } else {
        // generate a candidate plan to change lanes
        plan = get_lane_change_plan(pose, goal_lane, car_speed);

        // check if that plan (nearly) collides with any other cars
        int collision_id = collides(plan, predictions);
        if (collision_id != -1) {
          // if it does, discard the plan and simply hold our lane. We'll try
          // again next tick
          plan = get_hold_lane_plan(pose, goal_speed);
        } else {
          // if no collision adopt the LC plan as the current plan and
          // transition in the LC state
          cout << "made lane chng plan, PREP_LC→LC w goal lane " << goal_lane << endl;

          // do state transition here
          intent = Intent::LC;
        }
      }
    }

    // save this intent and plan for next time
    prev_intent = intent;
    prev_plan = plan;

    return plan;
  }
};

