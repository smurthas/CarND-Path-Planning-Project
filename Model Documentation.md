# Generating paths

The approach this project uses to generating paths involves three states:

1. *HOLD:* hold the current lane, keep a gap to the car in front
2. *PREP_LC:* when a car ahead is driving slowly, this state will be enabled and
the car will continually check the other lanes for gaps.
3. *LC:* this is the state during a lane change.

While the vehicle is holding its current lane, the generated path simply follows
the lanes centerline using its frenet D coordinate, advancing the S coordinate
to maintain the proper speed and distance from the car ahead.

During a lane change, paths are generated using a JMT calculation in frenet
coordinates along the D axis.

