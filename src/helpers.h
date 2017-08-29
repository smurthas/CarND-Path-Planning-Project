#include <fstream>
#include <math.h>
//#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>

#include "map.h"
#include <string>
#include <sstream>

using namespace std;

Map load_map_from_file(string map_file_) {
  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
    istringstream iss(line);
    double x;
    double y;
    float s;
    float d_x;
    float d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    map_waypoints_x.push_back(x);
    map_waypoints_y.push_back(y);
    map_waypoints_s.push_back(s);
    map_waypoints_dx.push_back(d_x);
    map_waypoints_dy.push_back(d_y);
  }
  Map map(map_waypoints_x, map_waypoints_y, map_waypoints_s);
  return map;
}
