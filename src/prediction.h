#include <vector>

#include "pose.h"

using namespace std;

struct FuturePose {
  double x;
  double y;
  double theta;
  double t;
};

struct PredictedObject {
  int id;
  vector<Pose> poses;
};

class Predictor {
  Map map;
 public:
  Predictor(const Map& _map) : map(_map) {}

  vector<PredictedObject> get_predictions(vector<vector<double> > sensor_input) {
    vector<PredictedObject> predictions;

    for (auto vehicle : sensor_input) {
      PredictedObject car;
      car.id = vehicle[0];

      double ds_dt = sqrt(vehicle[3]*vehicle[3] + vehicle[4]*vehicle[4]);
      double s = vehicle[5];
      double d = vehicle[6];

      Pose now;
      now.x = vehicle[1];
      now.y = vehicle[2];
      now.s = s;
      now.d = d;
      car.poses = {now};
      // just hold lane for now
      for (double t = 0.02; t <= 6; t += 0.02) {
        car.poses.push_back(map.get_cartesian(s + (ds_dt*t), d));
      }
      predictions.push_back(car);
    }

    return predictions;
  }
};

