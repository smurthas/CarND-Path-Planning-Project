#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
/*#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/LU"*/
#include "json.hpp"

#include "helpers.h"
#include "map.h"
#include "planner.h"

using namespace std;
//using namespace Eigen;

// for convenience
using json = nlohmann::json;
const double DT = 0.02; // 20 ms

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y) {

  double closestLen = 100000; //large number
  int closestWaypoint = 0;

  for(int i = 0; i < maps_x.size(); i++) {
    double map_x = maps_x[i];
    double map_y = maps_y[i];
    double dist = distance(x,y,map_x,map_y);
    if(dist < closestLen) {
      closestLen = dist;
      closestWaypoint = i;
    }

  }

  return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y) {

  int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

  double map_x = maps_x[closestWaypoint];
  double map_y = maps_y[closestWaypoint];

  double heading = atan2( (map_y-y),(map_x-x) );

  double angle = abs(theta-heading);

  if(angle > pi()/4) {
    closestWaypoint++;
  }

  return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y) {
  int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

  int prev_wp;
  prev_wp = next_wp-1;
  if(next_wp == 0) {
    prev_wp  = maps_x.size()-1;
  }

  double n_x = maps_x[next_wp]-maps_x[prev_wp];
  double n_y = maps_y[next_wp]-maps_y[prev_wp];
  double x_x = x - maps_x[prev_wp];
  double x_y = y - maps_y[prev_wp];

  // find the projection of x onto n
  double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
  double proj_x = proj_norm*n_x;
  double proj_y = proj_norm*n_y;

  double frenet_d = distance(x_x,x_y,proj_x,proj_y);

  //see if d value is positive or negative by comparing it to a center point

  double center_x = 1000-maps_x[prev_wp];
  double center_y = 2000-maps_y[prev_wp];
  double centerToPos = distance(center_x,center_y,x_x,x_y);
  double centerToRef = distance(center_x,center_y,proj_x,proj_y);

  if(centerToPos <= centerToRef) {
    frenet_d *= -1;
  }

  // calculate s value
  double frenet_s = 0;
  for(int i = 0; i < prev_wp; i++) {
    frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
  }

  frenet_s += distance(0,0,proj_x,proj_y);

  return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y) {
  int prev_wp = -1;

  while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) )) {
    prev_wp++;
  }

  int wp2 = (prev_wp+1)%maps_x.size();

  double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s-maps_s[prev_wp]);

  double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
  double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

  double perp_heading = heading-pi()/2;

  double x = seg_x + d*cos(perp_heading);
  double y = seg_y + d*sin(perp_heading);

  return {x,y};

}

Pose get_pose(double x, double y, double speed, double yaw) {
  //cout << "get_pose 0" << endl;
  Pose pose = {};
  pose.x = x;
  pose.y = y;
  pose.x_dot = cos(yaw)*speed;
  pose.y_dot = sin(yaw)*speed;
  pose.x_ddot = 0;
  pose.y_ddot = 0;

  return pose;
}

Pose get_pose(vector<double> x, vector<double> y) {
  //cout << "get_pose 1" << endl;
  Pose pose;

  double prev_x_dot = (double(x[1]) - double(x[0])) / DT;
  pose.x = x[2];
  pose.x_dot = (double(x[2]) - double(x[1])) / DT;
  pose.x_ddot = (pose.x_dot - prev_x_dot) / DT;

  double prev_y_dot = (double(y[1]) - double(y[0])) / DT;
  pose.y = y[2];
  pose.y_dot = (double(y[2]) - double(y[1])) / DT;
  pose.y_ddot = (pose.y_dot - prev_y_dot) / DT;

  return pose;
}

int main() {
  uWS::Hub h;
  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  Map map = load_map_from_file(map_file_);
  Planner planner(map, DT);
  Predictor predictor(map);

  double GOAL_SPEED = 2; // m/s
  int mesg_count = 0;
  Pose prev_pose = {};
  Pose prev_prev_pose = {};

  h.onMessage([&GOAL_SPEED,&prev_pose,&prev_prev_pose,&mesg_count,&predictor,&planner,&map](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

          // Main car's localization Data
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
          double car_d = j[1]["d"];
          double car_yaw = deg2rad(j[1]["yaw"]);
          double car_speed = double(j[1]["speed"]) * 0.44704;

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];


          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

          vector<PredictedObject> predictions = predictor.get_predictions(sensor_fusion);

          json msgJson;
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          Pose planning_pose = get_pose(car_x, car_y, car_speed, car_yaw);
          planning_pose.s = car_s;
          planning_pose.d = car_d;
          Plan best = planner.get_plan(planning_pose, predictions, previous_path_x.size());
          //cout << "best: " << best << endl;
          for (const Pose& p : best.poses) {
            next_x_vals.push_back(p.x);
            next_y_vals.push_back(p.y);
          }

          /*if (mesg_count < 100) {
            cout << "mesg_count: " << mesg_count << endl;
            if (mesg_count == 0) {
              double s = car_s;
              double ds_step = 0;
              for (int i = 0; i < 1000; i++) {
                if (ds_step < 0.42) {
                  ds_step += 0.001;
                }
                s += ds_step;
                Pose p = map.get_cartesian(s, car_d);
                next_x_vals.push_back(p.x);
                next_y_vals.push_back(p.y);
              }
            } else {
              cout << "push" << endl;
              for(double x : previous_path_x) {
                next_x_vals.push_back(x);
              }
              for(double y : previous_path_y) {
                next_y_vals.push_back(y);
              }
              cout << "push_ed " << next_y_vals.size() << endl;
            }
            mesg_count++;
          } else {
            next_x_vals = {previous_path_x[0], previous_path_x[1], previous_path_x[2]};
            next_y_vals = {previous_path_y[0], previous_path_y[1], previous_path_y[2]};
            Pose planning_pose = get_pose(next_x_vals, next_y_vals);
            //cout << "planning pose: " <<planning_pose << endl;
            double ds = distance(next_x_vals[0], next_y_vals[0], next_x_vals[2], next_y_vals[2]);
            planning_pose.s = car_s + ds;

            Plan best = planner.get_plan(planning_pose, predictions, previous_path_x.size());
            cout << "best: " << best << endl;
            for (const Pose& p : best.poses) {
              next_x_vals.push_back(p.x);
              next_y_vals.push_back(p.y);
            }
          }*/

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"control\","+ msgJson.dump()+"]";

          //this_thread::sleep_for(chrono::milliseconds(1000));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}

