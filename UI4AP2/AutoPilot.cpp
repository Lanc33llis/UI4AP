// AutoPilotSelfDev.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <functional>
#include "AutoPilot.h"

//wasd

const auto PI = 3.141592653589793238462643383279502884L;

void r8vec_bracket3(int n, double t[], double tval, int* left)

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_BRACKET3 finds the interval containing or nearest a given value.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The routine always returns the index LEFT of the sorted array
//    T with the property that either
//    *  T is contained in the interval [ T[LEFT], T[LEFT+1] ], or
//    *  T < T[LEFT] = T[0], or
//    *  T > T[LEFT+1] = T[N-1].
//
//    The routine is useful for interpolation problems, where
//    the abscissa must be located within an interval of data
//    abscissas for interpolation, or the "nearest" interval
//    to the (extreme) abscissa must be found so that extrapolation
//    can be carried out.
//
//    This version of the function has been revised so that the value of
//    LEFT that is returned uses the 0-based indexing natural to C++.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, length of the input array.
//
//    Input, double T[N], an array that has been sorted into ascending order.
//
//    Input, double TVAL, a value to be bracketed by entries of T.
//
//    Input/output, int *LEFT.
//    On input, if 0 <= LEFT <= N-2, LEFT is taken as a suggestion for the
//    interval [ T[LEFT-1] T[LEFT] ] in which TVAL lies.  This interval
//    is searched first, followed by the appropriate interval to the left
//    or right.  After that, a binary search is used.
//    On output, LEFT is set so that the interval [ T[LEFT], T[LEFT+1] ]
//    is the closest to TVAL; it either contains TVAL, or else TVAL
//    lies outside the interval [ T[0], T[N-1] ].
//
{
	int high;
	int low;
	int mid;
	//
	//  Check the input data.
	//
	if (n < 2)
	{
		std::cerr << "\n";
		std::cerr << "R8VEC_BRACKET3 - Fatal error!\n";
		std::cerr << "  N must be at least 2.\n";
		std::exit(1);
	}
	//
	//  If *LEFT is not between 0 and N-2, set it to the middle value.
	//
	if (*left < 0 || n - 2 < *left)
	{
		*left = (n - 1) / 2;
	}
	//
	//  CASE 1: TVAL < T[*LEFT]:
	//  Search for TVAL in (T[I],T[I+1]), for I = 0 to *LEFT-1.
	//
	if (tval < t[*left])
	{
		if (*left == 0)
		{
			return;
		}
		else if (*left == 1)
		{
			*left = 0;
			return;
		}
		else if (t[*left - 1] <= tval)
		{
			*left = *left - 1;
			return;
		}
		else if (tval <= t[1])
		{
			*left = 0;
			return;
		}
		//
		//  ...Binary search for TVAL in (T[I],T[I+1]), for I = 1 to *LEFT-2.
		//
		low = 1;
		high = *left - 2;

		for (; ; )
		{
			if (low == high)
			{
				*left = low;
				return;
			}

			mid = (low + high + 1) / 2;

			if (t[mid] <= tval)
			{
				low = mid;
			}
			else
			{
				high = mid - 1;
			}
		}
	}
	//
	//  CASE 2: T[*LEFT+1] < TVAL:
	//  Search for TVAL in (T[I],T[I+1]) for intervals I = *LEFT+1 to N-2.
	//
	else if (t[*left + 1] < tval)
	{
		if (*left == n - 2)
		{
			return;
		}
		else if (*left == n - 3)
		{
			*left = *left + 1;
			return;
		}
		else if (tval <= t[*left + 2])
		{
			*left = *left + 1;
			return;
		}
		else if (t[n - 2] <= tval)
		{
			*left = n - 2;
			return;
		}
		//
		//  ...Binary search for TVAL in (T[I],T[I+1]) for intervals I = *LEFT+2 to N-3.
		//
		low = *left + 2;
		high = n - 3;

		for (; ; )
		{

			if (low == high)
			{
				*left = low;
				return;
			}

			mid = (low + high + 1) / 2;

			if (t[mid] <= tval)
			{
				low = mid;
			}
			else
			{
				high = mid - 1;
			}
		}
	}
	//
	//  CASE 3: T[*LEFT] <= TVAL <= T[*LEFT+1]:
	//  T is just where the user said it might be.
	//
	else
	{
	}

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

void hermite_cubic_spline_value(int nn, double xn[], double fn[],
	double dn[], int n, double x[], double f[], double d[], double s[],
	double t[])

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    HERMITE_CUBIC_SPLINE_VALUE evaluates a Hermite cubic spline.
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
	//    Input, int NN, the number of data points.
	//
	//    Input, double XN[NN], the coordinates of the data points.
	//    The entries in XN must be in strictly ascending order.
	//
	//    Input, double FN[NN], the function values.
	//
	//    Input, double DN[NN], the derivative values.
	//
	//    Input, int N, the number of sample points.
	//
	//    Input, double X[N], the coordinates of the sample points.
	//
	//    Output, double F[N], the function value at the sample points.
	//
	//    Output, double D[N], the derivative value at the sample points.
	//
	//    Output, double S[N], the second derivative value at the
	//    sample points.
	//
	//    Output, double T[N], the third derivative value at the
	//    sample points.
	//
{
	int i;
	int left;

	left = n / 2;

	for (i = 0; i < n; i++)
	{
		r8vec_bracket3(nn, xn, x[i], &left);

		hermite_cubic_value(xn[left], fn[left], dn[left], xn[left + 1],
			fn[left + 1], dn[left + 1], 1, x + i, f + i, d + i, s + i, t + i);
	}
	return;
}

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
		return SplineFunction{ PointOne, PointTwo, *A3, *A2, *A1, *A0, Waypoint::Types::Hermite };
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
				auto a1 = (100 + 60) / 2;
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
