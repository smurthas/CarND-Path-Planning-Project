#include <cassert>
#include <fstream>
#include <math.h>
#include <chrono>
#include <iostream>
#include <vector>

#include "helpers.h"

#include "map.h"
#include "planner.h"

using namespace std;

int main() {
  cout << "hello!" << endl;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  double DT = 0.02;
  //Map map = load_map_from_file(map_file_);
  vector<double> x_wp = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  vector<double> y_wp = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  vector<double> s_wp = {0.0, 1.41, 2.82, 4.24, 5.66, 7.07};
  Map map(x_wp, y_wp, s_wp);
  Pose origin = map.get_cartesian(0, 0);
  assert(origin.x == 0.0);
  assert(origin.y == 0.0);
  Planner planner(map, DT);
  Predictor predictor(map);

  Pose planning_pose = map.get_cartesian(0.0, 6.0);
  cout << planning_pose << endl;
  assert(int(planning_pose.x*100) == 424);
  assert(int(planning_pose.y*100) == -424);
  planning_pose.x = 0;
  planning_pose.y = 0;
  planning_pose.x_dot = 0;
  planning_pose.y_dot = 0;
  planning_pose.x_ddot = 0;
  planning_pose.y_ddot = 0;
  planning_pose.s = 0;

  vector<PredictedObject> predictions;
  Plan best = planner.get_plan(planning_pose, predictions);

  cout << "best: " << best << endl;

}

