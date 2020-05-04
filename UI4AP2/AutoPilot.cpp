// AutoPilotSelfDev.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "AutoPilot.h"

//wasd123

const auto PI = 3.141592653589793238462643383279502884L;

void AP::hermite_cubic_to_power_cubic(double x1, double f1, double d1, double x2,
	double f2, double d2, double* c0, double* c1, double* c2, double* c3)

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    HERMITE_CUBIC_TO_POWER_CUBIC converts a Hermite cubic to power form.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    13 February 2011
	//
	//  Author:
	//
	//    John Burkardt
	//
	//  Reference:
	//
	//    Fred Fritsch, Ralph Carlson,
	//    Monotone Piecewise Cubic Interpolation,
	//    SIAM Journal on Numerical Analysis,
	//    Volume 17, Number 2, April 1980, pages 238-246.
	//
	//  Parameters:
	//
	//    Input, double X1, F1, D1, the left endpoint, function value
	//    and derivative.
	//
	//    Input, double X2, F2, D2, the right endpoint, function value
	//    and derivative.
	//
	//    Output, double *C0, *C1, *C2, *C3, the power form of the polynomial.
	//
{
	double df;
	double h;

	h = x2 - x1;
	df = (f2 - f1) / h;
	//
	//  Polynomial in terms of X - X1:
	//
	*c0 = f1;
	*c1 = d1;
	*c2 = -(2.0 * d1 - 3.0 * df + d2) / h;
	*c3 = (d1 - 2.0 * df + d2) / h / h;
	//
	//  Shift polynomial to X.
	//
	*c2 = *c2 - x1 * *c3;
	*c1 = *c1 - x1 * *c2;
	*c0 = *c0 - x1 * *c1;
	*c2 = *c2 - x1 * *c3;
	*c1 = *c1 - x1 * *c2;
	*c2 = *c2 - x1 * *c3;

	return;
}

namespace AP {

	Waypoint::Waypoint(double x, double y, double a) : X(x), Y(y), Angle(a) {};

	Waypoint::Waypoint() : X(0), Y(0), Angle(0) {};

	inline bool operator == (const Waypoint& lhs, const Waypoint& rhs)
	{
		if (lhs.X == rhs.X && lhs.Y == rhs.Y && lhs.Angle == rhs.Angle)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	//Tangent basically 
	double Angle2Deriv(double AngleInDegrees)
	{
		double a = AngleInDegrees;
		if (a == 90 || a == 270)
		{
			return NAN;
		}
		else
		{
			return tan(a * (3.141592653589793238462643383279502884197 / 180));
		}
	}

	//uses arc length formula to find distance
	double ArcLengthDistance(SplineFunction TheSplineFunction)
	{   //deriv = 3ax^2 + 2bx + c
		double Ax = 3 * TheSplineFunction.Ax, Bx = 2 * TheSplineFunction.Bx, C = TheSplineFunction.Cx;
		return (TheSplineFunction.PointTwo.X - TheSplineFunction.PointOne.X) * (sqrt(1 + pow(pow(Ax, 2) + Bx + C, 2)));
	}

	//Time given Spline Function and Jerk
	double TimeGivenSFJ(SplineFunction TheSplineFunction, double Jerk)
	{
		/*

		v(t) = interval [0, t] a(t)dt=1/2Jmax*t^2

		x(t) = interval [0, t] v(t)dt=1/2Jmax*t^3

		thus

		time for end position is:

		t = (6*end position/Jmax)^1/3

		and time for maximum speed

		t = (2*velocitymax/Jmax)^1/2

		*/

		double i = ArcLengthDistance(TheSplineFunction);
		return cbrt(6 * ArcLengthDistance(TheSplineFunction) / Jerk);
	}

	//Finds the cubic equation given two poitns
	SplineFunction HermiteFinder(Waypoint PointOne, Waypoint PointTwo)
	{
		// p(x) = c0 + c1 * x + c2 * x^2 + c3 * x^3
		// p(x) = c3 * x^3 + c2 * x^2 + c1 * x + c0
		double* A0 = new double(0), * A1 = new double(0), * A2 = new double(0), * A3 = new double(0);
		hermite_cubic_to_power_cubic(PointOne.X, PointOne.Y, Angle2Deriv(PointOne.Angle), PointTwo.X, PointTwo.Y, Angle2Deriv(PointTwo.Angle), A0, A1, A2, A3);
		return SplineFunction{ PointOne, PointTwo, *A3, *A2, *A1, *A0 };
	}

	//Generates a spline using a path structure
	Spline GenerateSpline(Path ThePath)
	{
		Spline ReturnSpline;
		int NumberOfFunctions = ThePath.size() - 1;
		for (int i = 0; i < NumberOfFunctions; i++)
		{
			SplineFunction Temp = HermiteFinder(ThePath[i], ThePath[i + 1]);
			ReturnSpline.push_back(Temp);
		}
		return ReturnSpline;
	}

	//Velocity calculator
	double Segment::Velocity(double Seconds)
	{
		//deriv 3ax^2 + 2bx + c
		if (Seconds > Time)
		{
			return NAN;
		}
		else
		{
			double Xa = XFunction.Ax, Xb = XFunction.Bx, Xc = XFunction.Cx;
			double Ya = YFunction.Ax, Yb = YFunction.Bx, Yc = YFunction.Cx;
			return sqrt(pow((3 * Xa * pow(Seconds, 2)) + (2 * Xb * Seconds) + Xc, 2) + pow((3 * Ya * pow(Seconds, 2)) + (2 * Yb * Seconds) + Yc, 2));
		}
	}

	//Acceleration Calculation
	double Segment::Acceleration(double Seconds)
	{
		//deriv 6ax + 2b
		if (Seconds > Time)
		{
			return NAN;
		}
		else
		{
			double Xa = XFunction.Ax, Xb = XFunction.Bx;
			double Ya = YFunction.Ax, Yb = YFunction.Cx;
			//return sqrt(pow((6 * Xa * Seconds + 2 * Xb), 2) + pow((6 * Ya * Seconds + 2 * Yb), 2));
			return sqrt(pow(((6 * Xa * Seconds) + (2 * Xb)), 2) + pow(((6 * Ya * Seconds) + (2 * Yb)), 2));
		}
	}

	//Jerk calculation
	double Segment::Jerk(double Seconds)
	{
		//deriv 6a
		if (Seconds > Time)
		{
			return NAN;
		}
		else
		{
			double Xa = XFunction.Ax;
			double Ya = YFunction.Ax;
			return sqrt(pow(6 * Xa, 2) + pow(6 * Ya, 2));
		}
	}

	//Generates a segment.
	Segment GenerateSegment(SplineFunction Function, double Jerk)
	{
		double B0, B1, B2, Time;
		Time = TimeGivenSFJ(Function, Jerk);
		SplineFunction XFunction = HermiteFinder(Waypoint{ 0, 0, Function.PointOne.Angle }, Waypoint{ Time, Function.PointTwo.X - Function.PointOne.X, Function.PointTwo.Angle });
		SplineFunction YFunction = HermiteFinder(Waypoint{ 0, 0, Function.PointOne.Angle }, Waypoint{ Time, Function.PointTwo.Y - Function.PointOne.Y, Function.PointTwo.Angle });
		B0 = Function.Ax; B1 = Function.Bx; B2 = Function.Cx;
		return Segment{ Function, XFunction, YFunction, B0, B1, B2, Time };
	}

	Trajectory GenerateTrajectory(Spline Spline, double Jerk)
	{
		//Trajectory ReturnTrajectory;
		//for (int i = 0; i < Spline.size(); i++)
		//{
		//	Segment Temp = GenerateSegment(Spline[i], Jerk);
		//	ReturnTrajectory.push_back(Temp);
		//}
		//return ReturnTrajectory;
		std::vector<double> a1, a2, x1, x2, y1, y2;
		std::vector<bool> p1rq, p2rq;
		for (int i = 0; i < Spline.size(); i++)
		{
			a1.push_back(Spline[i].PointOne.Angle);
			a2.push_back(Spline[i].PointTwo.Angle);
			x1.push_back(Spline[i].PointOne.X);
			x2.push_back(Spline[i].PointTwo.X);
			y1.push_back(Spline[i].PointOne.Y);
			y2.push_back(Spline[i].PointTwo.Y);
			if (a1[i] < 90 || a1[i] > 270)
			{
				p1rq.push_back(true);
			}
			else
			{
				p1rq.push_back(false);
			}
			if (a2[i] < 90 || a2[i] > 270)
			{
				p2rq.push_back(true);
			}
			else
			{
				p2rq.push_back(false);
			}
		}

		for (int i = 0; i < Spline.size(); i++)
		{
			if (p1rq[i] != p2rq[i])
			{
				if (p1rq[i])
				{
					AP::Waypoint p1{ x1[i], y1[i], a1[i] };
					AP::Waypoint p2{ x1[i] * 1.1, y2[i] * 1.1, 80 };
					AP::Waypoint p3{ x1[i] * 1.1, y2[i] * 1.1, 100 };
					AP::Waypoint p4{ x2[i], y2[i], a2[i] };
					AP::Path smoothPath{ p1, p2, p3, p4 };
					AP::Spline smoothSpline = AP::GenerateSpline(smoothPath);
					Spline.erase(Spline.begin() + i);
					p1rq.erase(p1rq.begin() + i);
					p2rq.erase(p2rq.begin() + i);
					x1.erase(x1.begin() + i);
					x2.erase(x2.begin() + i);
					y1.erase(y1.begin() + i);
					y2.erase(y2.begin() + i);
					a1.erase(a1.begin() + i);
					a2.erase(a2.begin() + i);
					for (int k = 0; k < smoothSpline.size(); k++)
					{
						Spline.emplace(Spline.begin() + i + k, smoothSpline[k]);
						p1rq.emplace(p1rq.begin() + i + k, true);
						p2rq.emplace(p2rq.begin() + i + k, true);
						x1.emplace(x1.begin() + i + k, 0);
						x2.emplace(x2.begin() + i + k, 0);
						y1.emplace(y1.begin() + i + k, 0);
						y2.emplace(y2.begin() + i + k, 0);
						a1.emplace(a1.begin() + i + k, 0);
						a2.emplace(a2.begin() + i + k, 0);
					}
				}
				else if (p2rq[i])
				{
					AP::Waypoint p1{ x1[i], y1[i], a1[i] };
					AP::Waypoint p2{ x1[i] * .6, y1[i] * 1.1, 100 };
					AP::Waypoint p3{ x1[i] * .6, y1[i] * 1.1, 60 };
					AP::Waypoint p4{ x2[i], y2[i], a2[i] };
					AP::Path smoothPath{ p1, p2, p3, p4 };
					AP::Spline smoothSpline = AP::GenerateSpline(smoothPath);
					Spline.erase(Spline.begin() + i);
					p1rq.erase(p1rq.begin() + i);
					p2rq.erase(p2rq.begin() + i);
					x1.erase(x1.begin() + i);
					x2.erase(x2.begin() + i);
					y1.erase(y1.begin() + i);
					y2.erase(y2.begin() + i);
					a1.erase(a1.begin() + i);
					a2.erase(a2.begin() + i);
					for (size_t k = 0; k < smoothSpline.size(); k++)
					{
						Spline.emplace(Spline.begin() + i + k, smoothSpline[k]);
						p1rq.emplace(p1rq.begin() + i + k, true);
						p2rq.emplace(p2rq.begin() + i + k, true);
						x1.emplace(x1.begin() + i + k, 0);
						x2.emplace(x2.begin() + i + k, 0);
						y1.emplace(y1.begin() + i + k, 0);
						y2.emplace(y2.begin() + i + k, 0);
						a1.emplace(a1.begin() + i + k, 0);
						a2.emplace(a2.begin() + i + k, 0);
					}
				}
			}
		}

		Trajectory ReturnTrajectory;
		for (int i = 0; i < Spline.size(); i++)
		{
			Segment Temp = GenerateSegment(Spline[i], Jerk);
			ReturnTrajectory.push_back(Temp);
		}
		return ReturnTrajectory;
	}

	Trajectory TrajectoryGeneration(AP::Path GroupOfWaypoints, double Jerk)
	{
		Trajectory ReturnTrajectory;
		Path modifiedGroupOfWaypoints(GroupOfWaypoints);
		for (size_t t = 0; t < GroupOfWaypoints.size() - 1; t++)
		{
			Waypoint WP1 = GroupOfWaypoints[t], WP2 = GroupOfWaypoints[t + 1];
			std::vector<Waypoint>::const_iterator iterator;
			if ((WP1.Angle < 90 || WP1.Angle > 270) && (WP2.Angle > 90 && WP2.Angle < 270))
			{
				Waypoint NP1(WP1.X * 1.1, WP1.Y * 1.1, 80);
				Waypoint NP2(WP1.X * 1.1, WP1.Y * 1.1, 100);
				iterator = find(modifiedGroupOfWaypoints.begin(), modifiedGroupOfWaypoints.end(), WP1);
				int i = (int)(iterator - modifiedGroupOfWaypoints.begin()) + 1;
				modifiedGroupOfWaypoints.emplace(modifiedGroupOfWaypoints.begin() + i, NP1);
				modifiedGroupOfWaypoints.emplace(modifiedGroupOfWaypoints.begin() + i + 1, NP2);
			}
			else if ((WP1.Angle > 90 && WP1.Angle < 270) && (WP2.Angle < 90 || WP2.Angle > 270))
			{
				Waypoint NP1(WP1.X * 1.1, WP1.Y * 1.1, 100);
				Waypoint NP2(WP1.X * 1.1, WP1.Y * 1.1, 80);
				iterator = find(modifiedGroupOfWaypoints.begin(), modifiedGroupOfWaypoints.end(), WP1);
				int i = (int)(iterator - modifiedGroupOfWaypoints.begin()) + 1;
				modifiedGroupOfWaypoints.emplace(modifiedGroupOfWaypoints.begin() + i, NP1);
				modifiedGroupOfWaypoints.emplace(modifiedGroupOfWaypoints.begin() + i + 1, NP2);
			}
		}
		Spline mySpline = GenerateSpline(modifiedGroupOfWaypoints);
		for (size_t i = 0; i < mySpline.size(); i++)
		{
			Segment Temp = GenerateSegment(mySpline[i], Jerk);
			ReturnTrajectory.push_back(Temp);
		}
		return ReturnTrajectory;
	}

	template<typename Encoder> void AutoTankDrive(TankConfig& TheTankConfig, Encoder& LeftEncoder, Encoder& RightEncoder, double TimeStep, double WidthBetweenWheels, double Jerk, int Ticks, double CircumferenceOfWheel)
	{
		double kP, kI, kD, MaxVelocity, Time;
		TheTankConfig.kP = !1.0 ? kP = TheTankConfig.kP : kP = 1.0; TheTankConfig.kI = !0.0 ? kP = TheTankConfig.kP : kP = 0.0;
		TheTankConfig.kD = !0.0 ? kP = TheTankConfig.kP : kP = 0.0; TheTankConfig.MaxVelocity = !2.0 ? MaxVelocity = TheTankConfig.kP : MaxVelocity = 2.0;
		Time = TheTankConfig.Time;
		int MaxSize = TheTankConfig.LeftTrajectory.size();
		for (int i = 0; i < MaxSize; i++)
		{
			double TotalTime = +TheTankConfig.LeftTrajectory[i].Time;
			if (Time > TotalTime)
			{
				continue;
			}
			else if (Time <= TotalTime)
			{
				double LeftError = TheTankConfig.LeftTrajectory[i].Velocity(Time) - (LeftEncoder.GetVelocity() / Ticks * CircumferenceOfWheel);
				TheTankConfig.LeftValue = (1 / MaxVelocity * LeftError);
				double RightError = TheTankConfig.RightTrajectory[i].Velocity(Time) - (RightEncoder.GetVelocity() / Ticks * CircumferenceOfWheel);
				TheTankConfig.RightValue = (1 / MaxVelocity * RightError);
			}
		}
		TheTankConfig.Time += TimeStep;
	}

	TankConfig GenerateTankConfig(Trajectory OriginalTrajectory, double WidthBetweenWheels, double Jerk)
	{
		TankConfig ReturnConfig;
		Spline* LTemp = new Spline;
		Spline* RTemp = new Spline;
		double radius = WidthBetweenWheels / 2;
		for (int i = 0; i < OriginalTrajectory.size(); i++)
		{
			double x1, y1, x2, y2, a1, a2;
			a1 = OriginalTrajectory[i].Function.PointOne.Angle + 90;
			a2 = OriginalTrajectory[i].Function.PointTwo.Angle + 90;
			x1 = radius * cos(a1 * (PI / 180)) + OriginalTrajectory[i].Function.PointOne.X;
			y1 = radius * sin(a1 * (PI / 180)) + OriginalTrajectory[i].Function.PointOne.Y;
			x2 = radius * cos(a2 * (PI / 180)) + OriginalTrajectory[i].Function.PointTwo.X;
			y2 = radius * sin(a2 * (PI / 180)) + OriginalTrajectory[i].Function.PointTwo.Y;
			Waypoint ModPoint1{ x1, y1, OriginalTrajectory[i].Function.PointOne.Angle }, ModPoint2{ x2, y2, OriginalTrajectory[i].Function.PointTwo.Angle };
			LTemp->push_back(HermiteFinder(ModPoint1, ModPoint2));
		}
		ReturnConfig.LeftTrajectory = GenerateTrajectory(*LTemp, Jerk);
		for (int i = 0; i < OriginalTrajectory.size(); i++)
		{
			double x1, y1, x2, y2, a1, a2;
			a1 = OriginalTrajectory[i].Function.PointOne.Angle - 90;
			a2 = OriginalTrajectory[i].Function.PointTwo.Angle - 90;
			x1 = radius * cos(a1 * (PI / 180)) + OriginalTrajectory[i].Function.PointOne.X;
			y1 = radius * sin(a1 * (PI / 180)) + OriginalTrajectory[i].Function.PointOne.Y;
			x2 = radius * cos(a2 * (PI / 180)) + OriginalTrajectory[i].Function.PointTwo.X;
			y2 = radius * sin(a2 * (PI / 180)) + OriginalTrajectory[i].Function.PointTwo.Y;
			Waypoint ModPoint1{ x1, y1, OriginalTrajectory[i].Function.PointOne.Angle }, ModPoint2{ x2, y2, OriginalTrajectory[i].Function.PointTwo.Angle };
			RTemp->push_back(HermiteFinder(ModPoint1, ModPoint2));
		}
		ReturnConfig.RightTrajectory = GenerateTrajectory(*RTemp, Jerk);
		ReturnConfig.Time = 0.0;
		return ReturnConfig;
	}
}

class TestEncoder
{
public:
	double GetVelocity();
	TestEncoder();
	~TestEncoder();
};

double TestEncoder::GetVelocity()
{
	int min = 1, max = 8;
	//std::random_device ran_device;
	//std::mt19937 generator(ran_device());
	//std::uniform_real_distribution<float> distr(min, max);
	//return distr(generator);
	return min;
}

TestEncoder::TestEncoder()
{
}

TestEncoder::~TestEncoder()
{
}
