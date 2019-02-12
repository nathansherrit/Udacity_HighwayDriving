#include <fstream> 
#include <math.h> 
#include <uWS/uWS.h> 
#include <chrono> 
#include <iostream> 
#include <thread> 
#include <vector> 
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
#include "Car.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() {
	return M_PI;
}
double deg2rad(double x) {
	return x * pi() / 180;
}
double rad2deg(double x) {
	return x * 180 / pi();
}

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
	auto found_null = s.find("null");
	auto b1 = s.find_first_of("[");
	auto b2 = s.find_first_of("}");
	if (found_null != string::npos) {
		return "";
	}
	else if (b1 != string::npos && b2 != string::npos) {
		return s.substr(b1, b2 - b1 + 2);
	}
	return "";
}

double distance(double x1, double y1, double x2, double y2) {
	return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}
int ClosestWaypoint(double x, double y, const vector <double> & maps_x, const vector <double> & maps_y) {

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for (int i = 0; i < maps_x.size(); i++) {
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x, y, map_x, map_y);
		if (dist < closestLen) {
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;
}

int NextWaypoint(double x, double y, double theta, const vector <double> & maps_x, const vector <double> &maps_y) {
	int closestWaypoint = ClosestWaypoint(x, y, maps_x, maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y - y), (map_x - x));

	double angle = fabs(theta - heading);
	angle = min(2 * pi() - angle, angle);

	if (angle > pi() / 2) {
		closestWaypoint++;
		if (closestWaypoint == maps_x.size()) {
			closestWaypoint = 0;
		}
	}

	return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector <double> getFrenet(double x, double y, double theta, const vector <double> & maps_x, const vector <double> &maps_y) {
	int next_wp = NextWaypoint(x, y, theta, maps_x, maps_y);

	int prev_wp;
	prev_wp = next_wp - 1;
	if (next_wp == 0) {
		prev_wp = maps_x.size() - 1;
	}

	double n_x = maps_x[next_wp] - maps_x[prev_wp];
	double n_y = maps_y[next_wp] - maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x * n_x + x_y * n_y) / (n_x * n_x + n_y * n_y);
	double proj_x = proj_norm * n_x;
	double proj_y = proj_norm * n_y;

	double frenet_d = distance(x_x, x_y, proj_x, proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000 - maps_x[prev_wp];
	double center_y = 2000 - maps_y[prev_wp];
	double centerToPos = distance(center_x, center_y, x_x, x_y);
	double centerToRef = distance(center_x, center_y, proj_x, proj_y);

	if (centerToPos <= centerToRef) {
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for (int i = 0; i < prev_wp; i++) {
		frenet_s += distance(maps_x[i], maps_y[i], maps_x[i + 1], maps_y[i + 1]);
	}

	frenet_s += distance(0, 0, proj_x, proj_y);

	return{
		frenet_s,
		frenet_d
	};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector <double> getXY(double s, double d, const vector <double> & maps_s, const vector <double> & maps_x, const vector <double> &maps_y) {
	int prev_wp = -1;

	while (s > maps_s[prev_wp + 1] && (prev_wp < (int)(maps_s.size() - 1))) {
		prev_wp++;
	}

	int wp2 = (prev_wp + 1) % maps_x.size();

	double heading = atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s - maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp] + seg_s * cos(heading);
	double seg_y = maps_y[prev_wp] + seg_s * sin(heading);

	double perp_heading = heading - pi() / 2;

	double x = seg_x + d * cos(perp_heading);
	double y = seg_y + d * sin(perp_heading);

	return{
		x,
		y
	};

}

bool checkCollisoion(double x1, double y1, double x2, double y2) {
	double width = 3;
	double length = 7;

	if (x1 - length / 2 < x2 + length / 2 && x1 + length / 2 > x2 - length / 2)
		if (y1 - width / 2 < y2 + width / 2 && y1 + width / 2 > y2 - width / 2)
			return true;
	return false;
}

bool isPathValid(Car c, tk::spline s, vector<vector<double>> sensor_fusion, const vector <double> & maps_s, const vector <double> & maps_x, const vector <double> &maps_y, double accel)
{
	double finalTime = 1.0;
	for (int j = 0; j < sensor_fusion.size(); j++)
	{
		double myX = 0;
		double myY = 0;
		vector<double> myFernet;

		double ref_yaw = deg2rad(c.yaw);
		double carVelocity = c.vel;

		double targetVX = sensor_fusion[j][3];
		double targetVY = sensor_fusion[j][4];
		double targetSpeed = sqrt(pow(targetVX, 2) + pow(targetVY, 2));
		double targetD = sensor_fusion[j][6]; // Assume target stays in lane
		double targetS = sensor_fusion[j][5];
		double targetAngle = atan2(targetVY, targetVX);

		for (double t = 0.02; t <= finalTime; t += 0.02)
		{
			double target_x = 60.0;
			double target_y = s(target_x);
			double target_dist = sqrt(target_x*target_x + target_y*target_y);
			double ref_vel_temp = c.vel;
			double N = (target_dist / (.02*carVelocity / 2.24));

			carVelocity += accel*0.02;
			if (carVelocity > 50)
				carVelocity = 50;
			if (carVelocity < 5)
				carVelocity = 5;

			myX = myX + (target_x) / N;
			myY = s(myX);

			targetS = sensor_fusion[j][5] + targetSpeed * t;


			double adjustedX = (myX *cos(ref_yaw) - myY*sin(ref_yaw));
			double adjustedY = (myX *sin(ref_yaw) + myY*cos(ref_yaw));
			adjustedX += c.x;
			adjustedY += c.y;

			myFernet = getFrenet(adjustedX, adjustedY, ref_yaw, maps_x, maps_y);

			if (checkCollisoion(myFernet[0], myFernet[1], targetS, targetD))
			{
				return false;
			}
		}
	}
	return true;
}
vector<double> findPathStats(Car c, int lane, tk::spline s, vector<vector<double>> sensor_fusion, const vector <double> & maps_s, const vector <double> & maps_x, const vector <double> &maps_y, double accel)
{
	double TTC = 0;
	double finalSpeed = c.vel + accel*2.5;
	double finalTime = 2.5;
	bool firstTarget = true;
	double carVelocity = c.vel;

	for (int j = 0; j < sensor_fusion.size(); j++)
	{
		double myX = 0;
		double myY = 0;
		vector<double> myFernet;

		double ref_yaw = deg2rad(c.yaw);
		double carVelocity = c.vel;

		double targetVX = sensor_fusion[j][3];
		double targetVY = sensor_fusion[j][4];
		double targetSpeed = sqrt(pow(targetVX, 2) + pow(targetVY, 2));
		double targetD = sensor_fusion[j][6]; // Assume target stays in lane
		double targetS = sensor_fusion[j][5];
		double targetAngle = atan2(targetVY, targetVX);

		for (double t = 0.02; t <= finalTime; t += 0.02)
		{
			double target_x = 60.0;
			double target_y = s(target_x);
			double target_dist = sqrt(target_x*target_x + target_y*target_y);
			double ref_vel_temp = c.vel;
			double N = (target_dist / (.02*carVelocity / 2.24));

			carVelocity += accel*0.02;
			if (carVelocity > 50)
				carVelocity = 50;
			if (carVelocity < 5)
				carVelocity = 5;

			myX = myX + (target_x) / N;
			myY = s(myX);

			targetS = sensor_fusion[j][5] + targetSpeed * t;


			double adjustedX = (myX *cos(ref_yaw) - myY*sin(ref_yaw));
			double adjustedY = (myX *sin(ref_yaw) + myY*cos(ref_yaw));
			adjustedX += c.x;
			adjustedY += c.y;

			myFernet = getFrenet(adjustedX, adjustedY, ref_yaw, maps_x, maps_y);

			if (checkCollisoion(myFernet[0], myFernet[1], targetS, targetD))
			{				
				return{ -1, 0, t };
			}
		}

		// Find the lane of the target car at the end of the time
		int targetLane = targetD / 4;
		double adjustedS = targetS - myFernet[0];
		int myLane = myFernet[1] / 4;
		carVelocity = c.vel + accel*finalTime;
		if (carVelocity > 50)
			carVelocity = 50;
		if (carVelocity < 5)
			carVelocity = 5;

		bool temp = false;
		// Check if target is in ego lane
		if (targetLane == lane && adjustedS > 0)
		{
			// Check if the time to collision is within 2 seconds
			double ttc = 0;
			double targetVel = sqrt(targetVX*targetVX + targetVY*targetVY);
			ttc = adjustedS / (carVelocity - targetVel);
			if (ttc < 2.0 && ttc > 0)
			{				
				if (ttc < TTC || firstTarget)
				{
					TTC = ttc;
					firstTarget = false;
					temp = true;
					finalSpeed = targetSpeed;
				}
			}
		}	
	}
	if (finalSpeed > 50)
      finalSpeed = 50;
  	if (finalSpeed < 5)
    	finalSpeed = 5;
	return{ TTC, finalSpeed, finalTime };
}

vector<double> findCost(Car c, vector<tk::spline> s, vector<vector<double>> sensor_fusion, const vector <double> & maps_s, const vector <double> & maps_x, const vector <double> &maps_y)
{
	double TTC;
	double finalSpeed = 50;
	vector<double> finalS;
	vector<double> cost;
	double finalTime[] = { 2.5,4,6 };
	double accel[] = { -5, 0, 10 };
	bool invalidPaths = true;
	// Check cost function for each spline		
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			vector<double> TTC_test = findPathStats(c, i, s[i], sensor_fusion, maps_s, maps_x, maps_y, accel[j]);

			TTC = TTC_test[0];
			finalSpeed = TTC_test[1];
			double lastTimeHolder = TTC_test[2];

			if (TTC_test[0] != -1)
			{
				invalidPaths = false;
				double target_x = 30.0;
				double target_y = s[i](target_x);
				double target_dist = sqrt(target_x*target_x + target_y*target_y);
				double N = (target_dist / (.02*c.vel / 2.24));

				// Find the TTC in the valid lanes
				cout << "lane " << i << ", accel " << accel[j] << " is a good lane with TTC = " << TTC << " and max speed = " << finalSpeed << endl;

				if (TTC <= 0)
					cost.push_back(0.001);
				else
					cost.push_back(1 / TTC);

				if (abs(i - c.lane) >= 2) // Add penalty for changing two lanes at once. 			
					cost[i] *= 10;
				else if (c.lane != i) // Adds a penlety for changing the lanes			
					cost[i] *= 5;

				cost[i * 3 + j] *= 1 / finalSpeed;
				cost[i * 3 + j] *= 1 / (20 + accel[j]);
			}
			else
			{
				cost.push_back(999999999);
				cout << "lane " << i << ", accel " << accel[j] << " is a bad lane due to collision" << endl;
			}
		}
	}
	//invalidPaths = true;
	if (invalidPaths)
	{
		cout << "All paths collide" << endl;
		return{ (double)c.lane, -5 };
	}
	for (int i = 0; i < cost.size(); i++)
	{
		cout << "Lane " << i << " Cost = " << cost[i] << endl;
	}

	double minCost = cost[0];
	int minIndex = 0;

	for (int i = 1; i < cost.size(); i++)
	{
		if (cost[i] < minCost)
		{
			minCost = cost[i];
			minIndex = i;
		}
	}

	return{ (double)(minIndex / 3), accel[minIndex % 3] };
}

int lane = 1;
double ref_vel = 0;
int goalLane = 1;
double accel = 10;
double targetTTC = 2.0;
int main()
{
	uWS::Hub h;

	// Load up map values for waypoint's x,y,s and d normalized normal vectors
	vector <double> map_waypoints_x;
	vector <double> map_waypoints_y;
	vector <double> map_waypoints_s;
	vector <double> map_waypoints_dx;
	vector <double> map_waypoints_dy;

	// Waypoint map to read from
	string map_file_ = "../data/highway_map.csv";
	// The max s value before wrapping around the track back to 0
	double max_s = 6945.554;

	ifstream in_map_(map_file_.c_str(), ifstream::in);

	string line;

	while (getline(in_map_, line))
	{
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

	h.onMessage([&map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy](uWS::WebSocket <uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
		// "42" at the start of the message means there's a websocket message event.
		// The 4 signifies a websocket message
		// The 2 signifies a websocket event
		//auto sdata = string(data).substr(0, length);
		//cout << sdata << endl;

		if (length && length > 2 && data[0] == '4' && data[1] == '2') {

			auto s = hasData(data);

			if (s != "") {
				auto j = json::parse(s);

				string event = j[0].get < string >();

				if (event == "telemetry") {
					// j[1] is the data JSON object

					// Main car's localization Data
					double car_x = j[1]["x"];
					double car_y = j[1]["y"];
					double car_s = j[1]["s"];
					double car_d = j[1]["d"];
					double car_yaw = j[1]["yaw"];
					double car_speed = j[1]["speed"];
					double end_path_s = j[1]["end_path_s"];
					double end_path_d = j[1]["end_path_d"];
					// Sensor Fusion Data, a list of all other cars on the same side of the road.
					auto sensor_fusion = j[1]["sensor_fusion"];
					// Previous path data given to the Planner
					auto previous_path_x = j[1]["previous_path_x"];
					auto previous_path_y = j[1]["previous_path_y"];
					Car c = Car(car_x, car_y, car_s, car_d, car_speed, car_yaw);
					cout << "#My Stats " << car_x << ", " << car_y << ", " << car_s << ", " << car_d << ", " << car_speed << ", " << car_yaw << endl;
					bool carInFront = false;
					// Find cars located in from of target

					double ttc = 0;
					for (int i = 0; i < sensor_fusion.size(); i++)
					{
						double targetID = sensor_fusion[i][0];
						double targetX = sensor_fusion[i][1];
						double targetY = sensor_fusion[i][2];
						double targetVX = sensor_fusion[i][3];
						double targetVY = sensor_fusion[i][4];
						double targetS = sensor_fusion[i][5];
						double targetD = sensor_fusion[i][6];
						//cout << "sensor_fusion.push_back(addTarget( " << targetID << ", " << targetX << ", " << targetY << ", " <<targetVX << ", " <<targetVY << ", " <<targetS << ", " <<targetD << ")); " <<endl;
						// Find the lane of the target car
						int targetLane = targetD / 4;
						double adjustedS = targetS - car_s;
						//cout << "Following distance = " << adjustedS << endl;
						// Check if target is in ego lane
						if (goalLane == targetLane && adjustedS > 0)
						{
							// Check if the time to collision is within 2 seconds
							//double ttc = 0;
							double targetVel = sqrt(targetVX*targetVX + targetVY*targetVY);
							ttc = adjustedS / (ref_vel - targetVel);

							if (ttc < 2.0 && ttc > 0)
							{
								cout << "Target in front" << endl;
								carInFront = true;
							}
						}
					}

					vector <vector<double>> ptsx;
					vector <vector<double>> ptsy;

					// Create 3 different paths, Go left, Go right, Stay straight
					int prev_size = previous_path_x.size();
					double ref_x = car_x;
					double ref_y = car_y;
					double ref_yaw = deg2rad(car_yaw);
					for (int i = 0; i < 3; i++)
					{
						if (prev_size < 2) {
							double prev_car_x = car_x - cos(car_yaw);
							double prev_car_y = car_y - sin(car_yaw);

							ptsx.push_back({ prev_car_x });
							ptsy.push_back({ prev_car_y });

							//ptsx[i].push_back(prev_car_x);
							ptsx[i].push_back(car_x);

							//ptsy[i].push_back(prev_car_y);
							ptsy[i].push_back(car_y);
						}
						else {
							ref_x = previous_path_x[prev_size - 1];
							ref_y = previous_path_y[prev_size - 1];

							double ref_x_prev = previous_path_x[prev_size - 2];
							double ref_y_prev = previous_path_y[prev_size - 2];

							ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

							ptsx.push_back({ ref_x_prev });
							ptsy.push_back({ ref_y_prev });

							ptsx[i].push_back(ref_x);
							ptsy[i].push_back(ref_y);
						}

						vector <double> next_wp0 = getXY(car_s + 60, (2 + 4 * i), map_waypoints_s, map_waypoints_x, map_waypoints_y);
						vector <double> next_wp1 = getXY(car_s + 120, (2 + 4 * i), map_waypoints_s, map_waypoints_x, map_waypoints_y);
						vector <double> next_wp2 = getXY(car_s + 180, (2 + 4 * i), map_waypoints_s, map_waypoints_x, map_waypoints_y);

						ptsx[i].push_back(next_wp0[0]);
						ptsx[i].push_back(next_wp1[0]);
						ptsx[i].push_back(next_wp2[0]);

						ptsy[i].push_back(next_wp0[1]);
						ptsy[i].push_back(next_wp1[1]);
						ptsy[i].push_back(next_wp2[1]);
					}

					for (int j = 0; j < ptsx.size(); j++)
					{
						//cout << "Spline " << j << " points" << endl;
						for (int i = 0; i < ptsx[j].size(); i++)
						{
							//cout << "(" << ptsx[j][i] << ", " <<  ptsy[j][i] << ")" << endl;
							double shift_x = ptsx[j][i] - ref_x;
							double shift_y = ptsy[j][i] - ref_y;

							ptsx[j][i] = (shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw));
							ptsy[j][i] = (shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw));

							//cout << "(" << ptsx[j][i] << ", " <<  ptsy[j][i] << ")" << endl;
						}
					}

					tk::spline s1;
					vector<tk::spline> s;

					for (int i = 0; i < 3; i++)
					{
						s1.set_points(ptsx[i], ptsy[i]);
						s.push_back(s1);
					}
					int bestSpline = 0;
					cout << "Current lane = " << c.lane << " Goal lane = " << goalLane << endl;
					bool checkPathValid = isPathValid(c, s[goalLane], sensor_fusion, map_waypoints_s, map_waypoints_x, map_waypoints_y, accel);
					checkPathValid = true;
					if (!checkPathValid)
					{
						cout << "Current path will have colision" << endl;
					}

					if (goalLane != c.lane)//&& checkPathValid)
					{
						cout << "Changing lanes" << endl;
						bestSpline = goalLane;
                      
                      	if(carInFront)
                        {
                          	accel = -5;
                          	cout << "Car in goal lane, slowing down" << endl;
                        }
					}
					else if (carInFront)
					{
						cout << "Target in front of car. Looking for new lane" << endl;
						for (int i = 0; i < sensor_fusion.size(); i++)
						{
							double targetID = sensor_fusion[i][0];
							double targetX = sensor_fusion[i][1];
							double targetY = sensor_fusion[i][2];
							double targetVX = sensor_fusion[i][3];
							double targetVY = sensor_fusion[i][4];
							double targetS = sensor_fusion[i][5];
							double targetD = sensor_fusion[i][6];
							cout << "sensor_fusion.push_back(addTarget( " << targetID << ", " << targetX << ", " << targetY << ", " << targetVX << ", " << targetVY << ", " << targetS << ", " << targetD << ")); " << endl;
						}

						//bestSpline = findCost(c, s, sensor_fusion, map_waypoints_s, map_waypoints_x, map_waypoints_y);
						vector<double> CostFunction = findCost(c, s, sensor_fusion, map_waypoints_s, map_waypoints_x, map_waypoints_y);
						bestSpline = CostFunction[0];
						accel = CostFunction[1];
						goalLane = bestSpline;

						cout << "Best option is spline " << bestSpline << " with accel of " << accel << endl;
					}
					else
					{						
						cout << "No need to change lanes" << endl;
						
						bestSpline = c.lane;
						goalLane = c.lane;

						accel = 10;
						
					}

					//cout << "Best option is spline " << bestSpline << endl;

					//cout << " Getting next vals" << endl;
					vector <double> next_x_vals;
					vector <double> next_y_vals;

					for (int i = 0; i < previous_path_x.size(); i++)
					{
						next_x_vals.push_back(previous_path_x[i]);
						next_y_vals.push_back(previous_path_y[i]);
					}

					double target_x = 30.0;
					double target_y = s[bestSpline](target_x);
					double target_dist = distance(target_x, target_y, 0, 0);

					double x_add_on = 0;

					//cout << "Generating next points" << endl;
					for (int i = 1; i <= 50 - previous_path_x.size(); i++)
					{
						if (goalLane != c.lane)
							ref_vel += accel * 0.02;
						else if (carInFront)
							ref_vel -= 10 * 0.02;
						else
							ref_vel += 10 * 0.02;
						if (ref_vel < 0)
							ref_vel = 0;
						else if (ref_vel > 49)
							ref_vel = 49;


						//cout << i << endl;
						double N = (target_dist / (.02*ref_vel / 2.24));
						double x_point = x_add_on + (target_x) / N;
						double y_point = s[bestSpline](x_point);

						x_add_on = x_point;
						//cout << target_dist << " " <<x_point << " " <<y_point << " " << N << " " << ref_vel << " " <<endl;
						double x_ref = x_point;
						double y_ref = y_point;

						x_point = (x_ref *cos(ref_yaw) - y_ref*sin(ref_yaw));
						y_point = (x_ref *sin(ref_yaw) + y_ref*cos(ref_yaw));

						x_point += ref_x;
						y_point += ref_y;

						next_x_vals.push_back(x_point);
						next_y_vals.push_back(y_point);
					}


					//cout << "path completed" << endl;                  	
					//for(int i = 0; i < next_x_vals.size(); i++)
					//{
					// cout << "point " << i << " ("<<next_x_vals[i] << "," << next_y_vals[i] << ")" << endl; 
					//}
					json msgJson;
					msgJson["next_x"] = next_x_vals;
					msgJson["next_y"] = next_y_vals;

					auto msg = "42[\"control\"," + msgJson.dump() + "]";

					//this_thread::sleep_for(chrono::milliseconds(1000));
					ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

				}
			}
			else {
				// Manual driving
				std::string msg = "42[\"manual\",{}]";
				ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
			}
		}
	});

	// We don't need this since we're not using HTTP but if it's removed the
	// program
	// doesn't compile :-(
	h.onHttpRequest([](uWS::HttpResponse * res, uWS::HttpRequest req, char *data, size_t, size_t) {
		const std::string s = "<h1>Hello world!</h1>";
		if (req.getUrl().valueLength == 1) {
			res->end(s.data(), s.length());
		}
		else {
			// i guess this should be done more gracefully?
			res->end(nullptr, 0);
		}
	});

	h.onConnection([&h](uWS::WebSocket < uWS::SERVER > ws, uWS::HttpRequest req) {
		std::cout << "Connected!!!" << std::endl;
	});

	h.onDisconnection([&h](uWS::WebSocket < uWS::SERVER > ws, int code, char * message, size_t length) {
		ws.close();
		std::cout << "Disconnected" << std::endl;
	});

	int port = 4567;
	if (h.listen(port)) {
		std::cout << "Listening to port " << port << std::endl;
	}
	else {
		std::cerr << "Failed to listen to port" << std::endl;
		return -1;
	}
	h.run();
}