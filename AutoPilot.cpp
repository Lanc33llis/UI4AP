// AutoPilotSelfDev.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <functional>
#include <cmath>
#include <string>
#include "AutoPilot.h"

const auto PI = 3.14159265358979323;

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

void hermite_cubic_value(double x1, double f1, double d1, double x2,
	double f2, double d2, int n, double x[], double f[], double d[],
	double s[], double t[])

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    HERMITE_CUBIC_VALUE evaluates a Hermite cubic polynomial.
	//
	//  Discussion:
	//
	//    The input arguments can be vectors.
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
	//    Input, int N, the number of evaluation points.
	//
	//    Input, double X[N], the points at which the Hermite cubic
	//    is to be evaluated.
	//
	//    Output, double F[N], D[N], S[N], T[N], the value and first
	//    three derivatives of the Hermite cubic at X.
	//
{
	double c2;
	double c3;
	double df;
	double h;
	int i;

	h = x2 - x1;
	df = (f2 - f1) / h;

	c2 = -(2.0 * d1 - 3.0 * df + d2) / h;
	c3 = (d1 - 2.0 * df + d2) / h / h;

	for (i = 0; i < n; i++)
	{
		f[i] = f1 + (x[i] - x1) * (d1
			+ (x[i] - x1) * (c2
				+ (x[i] - x1) * c3));
		d[i] = d1 + (x[i] - x1) * (2.0 * c2
			+ (x[i] - x1) * 3.0 * c3);
		s[i] = 2.0 * c2 + (x[i] - x1) * 6.0 * c3;
		t[i] = 6.0 * c3;
	}
	return;
}

double gTrunc(double in, int accuracy)
{
	auto c1 = trunc(in);

	auto c2 = in - c1;
	c2 *= pow(10, accuracy);
	c2 = trunc(c2);
	c2 *= pow(10, -accuracy);

	return c1 + c2;
}

namespace AP
{

	Waypoint::Waypoint(double x, double y, double a) : X(x), Y(y), Angle(a), typeOfFunction(Types::Hermite) {};

	Waypoint::Waypoint() : X(0), Y(0), Angle(0), typeOfFunction(Types::Hermite) {};

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

	inline bool operator == (const SplineFunction& lhs, const SplineFunction& rhs)
	{
		if ((lhs.Ax == rhs.Ax) && (lhs.Bx == rhs.Bx) && (lhs.Cx == rhs.Cx) && (lhs.Dx == rhs.Dx) && 
			(lhs.flag == rhs.flag) && (lhs.PointOne == rhs.PointOne) && (lhs.PointTwo == rhs.PointTwo))
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
			a *= (PI / 180);
			return round(tan(a));
		}
	}

	//uses arc length formula to find distance
	double ArcLengthDistance(SplineFunction TheSplineFunction)
	{   //deriv = 3ax^2 + 2bx + c
		double Ax = 3 * TheSplineFunction.Ax, Bx = 2 * TheSplineFunction.Bx, C = TheSplineFunction.Cx;
		return (TheSplineFunction.PointTwo.X - TheSplineFunction.PointOne.X) * (sqrt(1 + pow(pow(Ax, 2) + Bx + C, 2)));
	}

	Waypoint at(Spline aSpline, size_t index, double x)
	{
		//deriv = 3ax^2 + 2bx + c
		auto& c0 = aSpline[index].Ax, c1 = aSpline[index].Bx, c2 = aSpline[index].Cx, c3 = aSpline[index].Dx;

		auto y = (c0 * x * x * x) + (c1 * x * x) + (c2 * x) + c3;
		auto d = (3 * c0 * x * x) + (2 * c1 * x) + (c2);
		auto a = tan(d * (180 / PI));

		if (a < 0)
		{
			a = abs(a);
			//a = 360 - a;
		}

		if (aSpline[index].PointTwo.X - aSpline[index].PointOne.X < 0 && (int)gTrunc(a, 3) == 0)
		{
			a = 180;
		}
		

		return Waypoint(x, y, a);
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

		return cbrt(6 * ArcLengthDistance(TheSplineFunction) / Jerk);
	}

	//Finds the cubic equation given two poitns
	SplineFunction HermiteFinder(Waypoint PointOne, Waypoint PointTwo)
	{
		// p(x) = c0 + c1 * x + c2 * x^2 + c3 * x^3
		// p(x) = c3 * x^3 + c2 * x^2 + c1 * x + c0
		double* A0 = new double(0), * A1 = new double(0), * A2 = new double(0), * A3 = new double(0);
		hermite_cubic_to_power_cubic(PointOne.X, PointOne.Y, Angle2Deriv(PointOne.Angle), PointTwo.X, PointTwo.Y, Angle2Deriv(PointTwo.Angle), A0, A1, A2, A3);
		return SplineFunction{ PointOne, PointTwo, *A3, *A2, *A1, *A0, Waypoint::Types::Hermite };
		delete A0;
		delete A1;
		delete A2;
		delete A3;
	}

	SplineFunction QuadraticFinder(Waypoint pointOne, Waypoint pointTwo)
	{
		double A0, A1 = pointOne.X, A2 = pointOne.Y, A3 = 0;
		A0 = (pointTwo.Y - A2) / ((pow(pointTwo.X - (A1), 2)));
		return SplineFunction{ pointOne, pointTwo, A0, A1, A2, A3, Waypoint::Types::Quadratic };
	}

	SplineFunction SquareRootFinder(Waypoint pointOne, Waypoint pointTwo)
	{
		double A0, A1 = pointOne.X, A2 = pointOne.Y, A3 = 0;
		A0 = (pointTwo.Y - A2) / (copysign(1, pointTwo.X - A1) * sqrt(abs(pointTwo.X - A1)));
		return SplineFunction{ pointOne, pointTwo, A0, A1, A2, A3, Waypoint::Types::SquareRoot };
	}

	//Generates a spline using a path structure
	Spline GenerateSpline(Path ThePath)
	{
		Spline ReturnSpline;
		int NumberOfFunctions = ThePath.size() - 1;

		for (int i = 0; i < NumberOfFunctions; i++)
		{
			if (ThePath[i].typeOfFunction == Waypoint::Types::Quadratic)
			{
				SplineFunction Temp = QuadraticFinder(ThePath[i], ThePath[i + 1]);
				ReturnSpline.push_back(Temp);
			}
			else if (ThePath[i].typeOfFunction == Waypoint::Types::SquareRoot)
			{
				SplineFunction Temp = SquareRootFinder(ThePath[i], ThePath[i + 1]);
				ReturnSpline.push_back(Temp);
			}
			else
			{
				SplineFunction Temp = HermiteFinder(ThePath[i], ThePath[i + 1]);
				ReturnSpline.push_back(Temp);
			}
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
		SplineFunction XFunction, YFunction;
		if (Function.PointOne.typeOfFunction == Waypoint::Types::Quadratic)
		{
			XFunction = QuadraticFinder(Waypoint{ 0, 0, Function.PointOne.Angle }, Waypoint{ Time, Function.PointTwo.X - Function.PointOne.X, Function.PointTwo.Angle });
			YFunction = QuadraticFinder(Waypoint{ 0, 0, Function.PointOne.Angle }, Waypoint{ Time, Function.PointTwo.Y - Function.PointOne.Y, Function.PointTwo.Angle });
		}
		if (Function.PointOne.typeOfFunction == Waypoint::Types::SquareRoot)
		{
			XFunction = SquareRootFinder(Waypoint{ 0, 0, Function.PointOne.Angle }, Waypoint{ Time, Function.PointTwo.X - Function.PointOne.X, Function.PointTwo.Angle });
			YFunction = SquareRootFinder(Waypoint{ 0, 0, Function.PointOne.Angle }, Waypoint{ Time, Function.PointTwo.Y - Function.PointOne.Y, Function.PointTwo.Angle });
		}
		else
		{
			XFunction = HermiteFinder(Waypoint{ 0, 0, Function.PointOne.Angle }, Waypoint{ Time, Function.PointTwo.X - Function.PointOne.X, Function.PointTwo.Angle });
			YFunction = HermiteFinder(Waypoint{ 0, 0, Function.PointOne.Angle }, Waypoint{ Time, Function.PointTwo.Y - Function.PointOne.Y, Function.PointTwo.Angle });
		}
		B0 = Function.Ax; B1 = Function.Bx; B2 = Function.Cx;
		return Segment{ Function, XFunction, YFunction, B0, B1, B2, Time };
	}

	void interpolater(Waypoint p1, Waypoint p2, Waypoint& np1, Waypoint& np2, bool R2L)
	{
		auto x1 = p1.X, y1 = p1.Y, x2 = p2.X, y2 = p2.Y;
		if (R2L)
		{
			if (y1 < y2)
			{
				if (x1 < x2)
				{

				}
			}
		}
	}

	Path getPoints(Trajectory theTraj)
	{
		Path myPath;
		for (size_t i = 0; i < theTraj.size(); i++)
		{
			Waypoint temp(theTraj[i].Function.PointOne);
			myPath.push_back(temp);
			if (i == theTraj.size() - 1)
			{
				Waypoint temp2(theTraj[i].Function.PointTwo);
				myPath.push_back(temp2);
			}
		}
		return myPath;
	}

	Path getPoints(Spline theSpline)
	{
		Path myPath;
		for (size_t i = 0; i < theSpline.size(); i++)
		{
			Waypoint temp(theSpline[i].PointOne);
			myPath.push_back(temp);
			if (i == theSpline.size() - 1)
			{
				Waypoint temp2(theSpline[i].PointTwo);
				myPath.push_back(temp2);
			}
		}
		return myPath;
	}

	Trajectory TrajectoryGeneration(AP::Path GroupOfWaypoints, double Jerk, bool interpolationOveride)
	{
		Trajectory ReturnTrajectory;
		Path modifiedGroupOfWaypoints(GroupOfWaypoints);
		if (!interpolationOveride)
		{
			for (size_t t = 0; t < GroupOfWaypoints.size() - 1; t++)
			{
				Waypoint WP1 = GroupOfWaypoints[t], WP2 = GroupOfWaypoints[t + 1];
				std::vector<Waypoint>::const_iterator iterator;
				if ((WP1.Angle < 90 || WP1.Angle > 270) && (WP2.Angle > 90 && WP2.Angle < 270))
				{
					Waypoint NP1(WP1.X * 1.1, WP1.Y * 1.1, 80); NP1.typeOfFunction = Waypoint::Types::Hermite; //quad
					Waypoint NP2(WP1.X * 1.1, WP1.Y * 1.1, 100); NP2.typeOfFunction = Waypoint::Types::Hermite; //sqrt
					iterator = find(modifiedGroupOfWaypoints.begin(), modifiedGroupOfWaypoints.end(), WP1);
					int i = (int)(iterator - modifiedGroupOfWaypoints.begin()) + 1;
					modifiedGroupOfWaypoints.emplace(modifiedGroupOfWaypoints.begin() + i, NP1);
					modifiedGroupOfWaypoints.emplace(modifiedGroupOfWaypoints.begin() + i + 1, NP2);
					NP1.Waypoint::~Waypoint();
					NP2.Waypoint::~Waypoint();
				}
				else if ((WP1.Angle > 90 && WP1.Angle < 270) && (WP2.Angle < 90 || WP2.Angle > 270))
				{
					Waypoint NP1(WP1.X * .9, WP1.Y * 1.1, 100); NP1.typeOfFunction = Waypoint::Types::Hermite;
					Waypoint NP2(WP1.X * .9, WP1.Y * 1.1, 60); NP2.typeOfFunction = Waypoint::Types::Hermite;
					iterator = find(modifiedGroupOfWaypoints.begin(), modifiedGroupOfWaypoints.end(), WP1);
					int i = (int)(iterator - modifiedGroupOfWaypoints.begin()) + 1;
					modifiedGroupOfWaypoints.emplace(modifiedGroupOfWaypoints.begin() + i, NP1);
					modifiedGroupOfWaypoints.emplace(modifiedGroupOfWaypoints.begin() + i + 1, NP2);
					NP1.Waypoint::~Waypoint();
					NP2.Waypoint::~Waypoint();
				}
				WP1.Waypoint::~Waypoint();
				WP2.Waypoint::~Waypoint();
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

	Path GenerateNormal(Path thePath, double r, bool outerNorm)
	{
		//auto mySpline = GenerateSpline(thePath);
		//Path tempPath, returnPath;

		//for (size_t i = 0; i < mySpline.size(); i++)
		//{
		//	tempPath.push_back(mySpline[i].PointOne);
		//	if (!(mySpline[i].PointOne.X == mySpline[i].PointTwo.X && mySpline[i].PointOne.Y == mySpline[i].PointTwo.Y))
		//	{
		//		double c0 = mySpline[i].Ax, c1 = mySpline[i].Bx, c2 = mySpline[i].Cx, c3 = mySpline[i].Dx;
		//		std::function<double(double)> splineFunction = [c0, c1, c2, c3](double x) {return ((c0 * pow(x, 3)) + (c1 * pow(x, 2)) + (c2 * x) + c3); };
		//		
		//		auto &p1 = mySpline[i].PointOne;
		//		auto &p2 = mySpline[i].PointTwo;
		//		double domainFourth =( p2.X - p1.X) / 4;

		//		//deriv = 3ax^2 + 2bx + c
		//		double x[] = { mySpline[i].PointOne.X + domainFourth, mySpline[i].PointOne.X + domainFourth + domainFourth };
		//		double y[2];
		//		double d[2];
		//		double d2[2];
		//		double d3[2];


		//		bool is180 = false;

		//		if (p1.Angle == 180 && p2.Angle == 180)
		//			is180 = true;
		//		else 
		//			false;


		//		hermite_cubic_value(p1.X, p1.Y, Angle2Deriv(p1.Angle), p2.X, p2.Y,
		//			Angle2Deriv(p2.Angle), 2, x, y, d, d2, d3);				
		//		//hermite_cubic_value(mySpline[i].PointTwo.X, mySpline[i].PointTwo.Y, Angle2Deriv(mySpline[i].PointTwo.Angle), mySpline[i].PointOne.X, mySpline[i].PointOne.Y,
		//		//	Angle2Deriv(mySpline[i].PointOne.Angle), 2, x, y, d, d2, d3);

		//		auto test = (atan(d[0])) * 180 / PI;

		//		if (d[0] > 0 && d[1] > 0 && !is180)
		//		{
		//			d[0] = round(atan(d[0]) * 180 / PI);
		//			d[1] = round(atan(d[1]) * 180 / PI);
		//		}

		//		if (d[0] < 0)
		//		{
		//			d[0] = 180 + ((atan(d[0]) * 180 / PI));
		//		}

		//		if (d[1] < 0)
		//		{
		//			d[1] = 180 + ((atan(d[1]) * 180 / PI));
		//		}

		//		if (is180)
		//		{
		//			d[0] = 180;
		//			d[1] = 180;
		//		}


		//		Waypoint newPoint(x[0], y[0], d[0]);

		//		tempPath.push_back(newPoint);

		//		newPoint.X = x[1]; newPoint.Y = y[1]; newPoint.Angle = d[1];

		//		tempPath.push_back(newPoint);
		//	}
		//	if (i == mySpline.size() - 1)
		//		tempPath.push_back(mySpline[i].PointTwo);
		//}

		//for (size_t i = 0; i < tempPath.size(); i++)
		//{
		//	double x, y, a, c, s;

		//	a = tempPath[i].Angle + 90;
		//	c = r * (cos(a * (PI / 180)));
		//	s = r * (sin(a * (PI / 180)));

		//	if (outerNorm)
		//	{
		//		x = c + tempPath[i].X;
		//		y = s + tempPath[i].Y;
		//	}
		//	else
		//	{
		//		x = -c + tempPath[i].X;
		//		y = -s + tempPath[i].Y;
		//	}
		//	Waypoint l(x, y, tempPath[i].Angle);
		//	returnPath.push_back(l);
		//}
		//return returnPath;

		auto &oP = thePath; 
		auto oS = GenerateSpline(oP);

		Path nP;
		Path rP;

		for (Waypoint i : oP)
		{
			//double x, y, a, c, s;
			//a = i.Angle + 90;
			//c = r * (cos(a * (PI / 180)));
			//s = r * (sin(a * (PI / 180)));

			//if (outerNorm)
			//{
			//	x = c + i.X;
			//	y = s + i.Y;
			//}
			//else
			//{
			//	x = -c + i.X;
			//	y = -s + i.Y;
			//}

			//nP.push_back(Waypoint(x, y, a - 90));

			//deriv = 3ax ^ 2 + 2bx + c
		}

		for (size_t i = 0; i < nP.size(); i++)
		{
			rP.push_back(nP[i]);

			if (i + 1 == nP.size())
				break;

			if (nP[i] == nP[i + 1])
				continue;

			Waypoint &p1 = nP[i];
			Waypoint &p2 = nP[i + 1];

			double f = (nP[i + 1].X - nP[i].X) / 4;

			double x[] = { nP[i].X + f, nP[i].X + f + f }, y[2], d[2], d1[2], d2[2];

			hermite_cubic_value(p1.X, p1.Y, Angle2Deriv(p1.Angle), p2.X, p2.Y, Angle2Deriv(p2.Angle), 2, x, y, d, d1, d2);

			if (d[0] > 0)
			{
				d[0] = (atan(d[0]) * 180 / PI);
			}

			else if (d[0] < 0)
			{
				d[0] = 180 + ((atan(d[0]) * 180 / PI));
			}

			if (d[1] > 0)
			{
				d[1] = (atan(d[1]) * 180 / PI);
			}

			else if (d[1] < 0)
			{
				d[1] = 180 + ((atan(d[1]) * 180 / PI));
			}

			if (p1.Angle == 180 && p2.Angle == 180)
			{
				d[0] = 180;
				d[1] = 180;
			}

			rP.push_back(Waypoint(x[0], y[0], d[0]));
			rP.push_back(Waypoint(x[1], y[1], d[1]));
		}

		return rP;
	}

	Path FindNormal(Spline theSpline, double r, size_t accuracy, bool outer)
	{
		auto& s = theSpline;
		auto test = getPoints(s);
		Path c;
		for (size_t t = 0; t < theSpline.size(); t++)
		{
			c.push_back(s[t].PointOne);

			auto &p1 = s[t].PointOne;
			auto &p2 = s[t].PointTwo;

			auto d = p2.X - p1.X;

			auto ratio = d / accuracy;

			if (!(s[t].PointOne.X == s[t].PointTwo.X && s[t].PointOne.Y == s[t].PointTwo.Y))
			{
				for (size_t i = 1; i < accuracy; i++)
				{
					auto nx = ratio * i + p1.X;
					auto np = at(s, t, nx);
					c.push_back(np);
				}
			}

			if (t == s.size() - 1)
			{
				c.push_back(s[t].PointTwo);
			}
		}
		return c;
	}

	TankConfig GenerateTankConfig(Trajectory theTraj, double WidthBetweenWheels, double Jerk)
	{
		TankConfig myConfig;
		auto controlPath = getPoints(theTraj);
		auto temp = GenerateSpline(controlPath);
		Path leftWheelPath, rightWheelPath;

		auto r = WidthBetweenWheels / 2;

		//for (size_t i = 0; i < controlPath.size(); i++)
		//{
			//double x, y, a;
			//a = controlPath[i].Angle + 90;
			//x = r * cosf(a * (PI / 180)) + controlPath[i].X;
			//y = r * sinf(a * (PI / 180)) + controlPath[i].Y;
			//Waypoint l(x, y, a - 90);
			//l.typeOfFunction = controlPath[i].typeOfFunction;

			//a = controlPath[i].Angle - 90;
			//x = r * cosf(a * (PI / 180)) + controlPath[i].X;
			//y = r * sinf(a * (PI / 180)) + controlPath[i].Y;
			//Waypoint r(x, y, a + 90);
			//r.typeOfFunction = controlPath[i].typeOfFunction;

		//	double x, y, a, c, s;
		//	a = controlPath[i].Angle + 90;
		//	c = r * cosf(a * (PI / 180));
		//	s = r * sinf(a * (PI / 180));
		//	x = c + controlPath[i].X;
		//	y = s + controlPath[i].Y;
		//	Waypoint l(x, y, controlPath[i].Angle);
		//	leftWheelPath.push_back(l);
		//	x = -c + controlPath[i].X;
		//	y = -s + controlPath[i].Y;
		//	l.X = x; l.Y = y;
		//	rightWheelPath.push_back(l);

		//}

		leftWheelPath = FindNormal(temp, r, 4, true);
		rightWheelPath = FindNormal(temp, r, 4, true);

		myConfig.LeftTrajectory = TrajectoryGeneration(leftWheelPath, Jerk, true);
		myConfig.RightTrajectory = TrajectoryGeneration(rightWheelPath, Jerk, true);

		return myConfig;
	}

//	TankConfig GenerateTankConfig(Trajectory theTraj, double WidthBetweenWheels, double Jerk)
//	{
//		TankConfig myConfig;
//		Path lPath, rPath, trajPath = getPoints(theTraj);
//		auto radius = WidthBetweenWheels / 2;
//
//		for (size_t i = 0; i < theTraj.size(); i++)
//		{
//			double c0 = theTraj[i].Function.Ax, c1 = theTraj[i].Function.Bx, c2 = theTraj[i].Function.Cx, c3 = theTraj[i].Function.Dx;
//			int flag = theTraj[i].Function.flag;
//
//			//gotta find deriv
//
//			std::function<double(double)> Function = [c0, c1, c2, c3, flag](double x)
//			{
//				switch (flag)
//				{
//				case AP::Waypoint::Types::Quadratic:
//					if (flag == (AP::Waypoint::Types::Quadratic | AP::Waypoint::Types::HorizontallyFlipped))
//					{
//
//					}
//
//				case AP::Waypoint::Types::SquareRoot:
//				
//				case AP::Waypoint::Types::Hermite:
//
//				}
//				if (flag == AP::Waypoint::Types::Quadratic)
//				{
//					return ((c0 * pow(x - c1, 2)) + c2);
//				}
//				else if (flag == AP::Waypoint::Types::SquareRoot)
//				{
//					return ((c0 * sqrt(x - c1)) + c2);
//				}
//				else if (flag == AP::Waypoint::Types::Hermite)
//				{
//					return ((c0 * pow(x, 3)) + (c1 * pow(x, 2)) + (c2 * x) + c3);
//				}
//			};
//
//
//			if ((theTraj[i].Function.PointTwo.X - theTraj[i].Function.PointOne.X) < 0)
//			{
//				for (double g = theTraj[i].Function.PointOne.X; g >= theTraj[i].Function.PointTwo.X; g -= .005)
//				{
//					points.push_back(cv::Point(g * xRatio, Function(g) * yRatio));
//				}
//			}
//
//			else
//			{
//				for (double g = theTraj[i].Function.PointOne.X; g <= theTraj[i].Function.PointTwo.X; g += .005)
//				{
//					points.push_back(cv::Point(g * xRatio, Function(g) * yRatio));
//				}
//			}
//
//
//			double x, y, a;
//			a = trajPath[i].Angle + 90;
//			x = radius * cos(a * (PI / 180)) + trajPath[i].X;
//			y = radius * sin(a * (PI / 180)) + trajPath[i].Y;
//			Waypoint l(x, y, a - 90);
//			l.typeOfFunction = trajPath[i].typeOfFunction;
//
//			a = trajPath[i].Angle - 90;
//			x = radius * cos(a * (PI / 180)) + trajPath[i].X;
//			y = radius * sin(a * (PI / 180)) + trajPath[i].Y;
//			Waypoint r(x, y, a + 90);
//			r.typeOfFunction = trajPath[i].typeOfFunction;
//
//			lPath.push_back(l);
//			rPath.push_back(r);
//		}
//
//		myConfig.LeftTrajectory = TrajectoryGeneration(lPath, Jerk);
//		myConfig.RightTrajectory = TrajectoryGeneration(rPath, Jerk);
//
//		return myConfig;
//	}
//
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
