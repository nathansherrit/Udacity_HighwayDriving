#pragma once
class Car
{
public:
	Car();
	Car(double, double, double, double, double, double);
	~Car();

	double x;
   	double y;
	double s;
	double d;
	double vel;
    double yaw;
	int lane;
};

Car::Car(double x_in,double y_in,double s_in, double d_in, double vel_in, double yaw_in)
{
  	x = x_in;
  	y = y_in;
	s = s_in;
	d = d_in;
	vel = vel_in;
  	yaw = yaw_in;
  
	lane = d / 4;
}

Car::~Car()
{
}

