#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/LU"
#include "json.hpp"
#include "spline.h"

using namespace std;
using namespace Eigen;

// for convenience
using json = nlohmann::json;

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

double distance(double x1, double y1, double x2, double y2) {
  return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
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

vector<double> get_cartesian(double s, double d, tk::spline spline_x, tk::spline spline_y) {
  double x0 = spline_x(s);
  double y0 = spline_y(s);
  double ds = 0.001;
  double x1 = spline_x(s+ds);
  double y1 = spline_y(s+ds);
  double d_heading = atan2(y1-y0, x1-x0) - M_PI/2.0;
  double x = x0 + cos(d_heading)*d;
  double y = y0 + sin(d_heading)*d;
  return { x, y };
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
  VectorXd esses(3);
  esses << (s_f - (s_i + s_dot_i*T + 0.5*s_ddot_i*T*T)),
           (s_dot_f - (s_dot_i + s_ddot_i*T)),
           (s_ddot_f - s_ddot_i);

  MatrixXd Ts(3, 3);
  double T_2 = T*T;
  double T_3 = T_2*T;
  double T_4 = T_3*T;
  double T_5 = T_4*T;
  Ts << T_3,    T_4,    T_5,
      3*T_2,  4*T_3,  5*T_4,
        6*T, 12*T_2, 20*T_3;
  VectorXd As = Ts.inverse() * esses;

  return {s_i, s_dot_i, 0.5*s_ddot_i, As[0], As[1], As[2]};
  /*
  Calculate the Jerk Minimizing Trajectory that connects the initial state
  to the final state in time T.

  INPUTS

  start - the vehicles start location given as a length three array
      corresponding to initial values of [s, s_dot, s_double_dot]

  end   - the desired end state for vehicle. Like "start" this is a
      length three array.

  T     - The duration, in seconds, over which this maneuver should occur.

  OUTPUT
  an array of length 6, each value corresponding to a coefficent in the polynomial
  s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5

  EXAMPLE

  > JMT( [0, 10, 0], [10, 10, 0], 1)
  [0.0, 10.0, 0.0, 0.0, 0.0, 0.0]
 */
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

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

  tk::spline spline_x;
  spline_x.set_points(map_waypoints_s, map_waypoints_x);

  tk::spline spline_y;
  spline_y.set_points(map_waypoints_s, map_waypoints_y);

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&spline_x,&spline_y](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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
          //cout << "car_speed: " << car_speed << endl;

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

          json msgJson;
          vector<double> next_x_vals = {car_x};
          vector<double> next_y_vals = {car_y};

          const double DT = 0.02; // 20 ms
          const double GOAL_SPEED = 15; // m/s


          // TODO: pick a lane
          double d = 6;

          double last_x = car_x;
          double last_x_dot = car_speed*cos(car_yaw);
          double last_x_ddot = 0;
          if (previous_path_x.size() > 2) {
            last_x = previous_path_x[2];
            last_x_dot = (double(previous_path_x[2]) - double(previous_path_x[1])) / DT;
            double prev_x_dot = (double(previous_path_x[1]) - double(previous_path_x[0])) / DT;
            last_x_ddot = (last_x_dot - prev_x_dot) / DT;
            next_x_vals = {previous_path_x[0], previous_path_x[1], previous_path_x[2]};
          } else {
            cout << "WARNING: previous_path_x.size(): " << previous_path_x.size() << endl;
          }

          double last_y = car_y;
          double last_y_dot = car_speed*sin(car_yaw);
          double last_y_ddot = 0;
          if (previous_path_y.size() > 2) {
            last_y = previous_path_y[2];
            last_y_dot = (double(previous_path_y[2]) - double(previous_path_y[1])) / DT;
            double prev_y_dot = (double(previous_path_y[1]) - double(previous_path_y[0])) / DT;
            last_y_ddot = (last_y_dot - prev_y_dot) / DT;
            next_y_vals = {previous_path_y[0], previous_path_y[1], previous_path_y[2]};
          }

          const double pre_planned_time = (next_x_vals.size() - 1) * DT;
          const double T = 1 - pre_planned_time;
          //double s_last = getFrenet(last_x, last_y, car_yaw, map_waypoints_x, map_waypoints_y)[0];
          double s_last = car_s + pre_planned_time*car_speed;
          double s_end = s_last + GOAL_SPEED * T;
          vector<double> end_point = get_cartesian(s_end, d, spline_x, spline_y);
          double end_x = end_point[0];
          double end_y = end_point[1];
          if (distance(last_x, last_y, end_x, end_y) > 3) {
            cout << "boink!" << endl;
          }
          cout << "T: " << T << ", s_last: " << s_last << ", s_end: " << s_end << endl;
          vector<double> end_point1 = get_cartesian(s_end+GOAL_SPEED*DT, d, spline_x, spline_y);
          double end_x1 = end_point1[0];
          double end_y1 = end_point1[1];
          double end_x_dot = (end_x1 - end_x) / DT;
          double end_y_dot = (end_y1 - end_y) / DT;

          vector<double> start_x_v = {last_x, last_x_dot, last_x_ddot};
          vector<double> end_x_v = {end_x, end_x_dot, 0.0};

          vector<double> start_y_v = {last_y, last_y_dot, last_y_ddot};
          vector<double> end_y_v = {end_y, end_y_dot, 0.0};


          vector<double> a_x = JMT(start_x_v, end_x_v, T);
          vector<double> a_y = JMT(start_y_v, end_y_v, T);

          for (double t = DT; t < T; t += DT) {
            next_x_vals.push_back(eval_poly(a_x, t));
            next_y_vals.push_back(eval_poly(a_y, t));
          }

          /*double last_s = car_s;
          double last_heading = car_yaw;
          if (previous_path_y.size() >= 5) {
            int prev_points = 5;
            for (int i = 0; i < prev_points; i++) {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }

            double last_x = previous_path_x[prev_points-1];
            double last_y = previous_path_y[prev_points-1];
            last_s = getFrenet(last_x, last_y, car_yaw, map_waypoints_x, map_waypoints_y)[0];
          }

          for (double s = last_s+0.4; s < last_s+15.0; s += 0.4) {
            double x0 = spline_x(s);
            double y0 = spline_y(s);
            double ds = 0.01;
            double x1 = spline_x(s+ds);
            double y1 = spline_y(s+ds);
            double d_heading = atan2(y1-y0, x1-x0) - M_PI/2.0;
            double x = x0 + cos(d_heading)*d;
            double y = y0 + sin(d_heading)*d;
            cout << "d_heading(" << s << "): " << d_heading << ", cos: " << cos(d_heading) << ", sin: " << sin(d_heading) << ", (" << x << "," << y << ")" << endl;
            next_x_vals.push_back(x);
            next_y_vals.push_back(y);
          }*/

          cout << next_x_vals.size() << ", " << next_y_vals.size() << endl;

          // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
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

