#pragma once

#include <vector> 

namespace AP
{

	class AutoPilot
	{
	public:
		double MaxJerk;
		AutoPilot(double Jerk);
		~AutoPilot();
	};

	void hermite_cubic_to_power_cubic(double x1, double f1, double d1, double x2, double f2, double d2, double* c0, double* c1, double* c2, double* c3);


	struct Waypoint
	{
		double X, Y, Angle;
		enum Types
		{
			None = 0,
			Quadratic = 1 << 0,
			SquareRoot = 1 << 1,
			Hermite = 1 << 2,
			Exponential = 1 << 3,
			VerticallyFlipped = 1 << 4,
			HorizontallyFlipped = 1 << 5,
		};
		int typeOfFunction;
		Waypoint(double x, double y, double angleInDegrees);
		Waypoint();
	};

	inline bool operator == (const Waypoint& lhs, const Waypoint& rhs);

	typedef std::vector<Waypoint> Path;

	struct SplineFunction
	{
		Waypoint PointOne, PointTwo;
		double Ax, Bx, Cx, Dx;
		int flag;
	};

	typedef std::vector<SplineFunction> Spline;

	struct Segment
	{
		SplineFunction Function, XFunction, YFunction;
		double Velocity(double Seconds); double Acceleration(double Seconds); double Jerk(double Seconds);
		double Ax, Bx, C, Time;
	};

	typedef std::vector<Segment> Trajectory;

	//finds deriv or slope of an angle. Used to find exit angle slope or deriv
	double Angle2Deriv(double AngleInDegrees);

	//Finds the distance of a spline function using arc length
	double ArcLengthDistance(SplineFunction TheSplineFunction);

	//Finds the time needed to traverse the distance given velocity, acceleration, and jerk
	double TimeGivenSFJ(SplineFunction TheSplineFunction, double Jerk);

	SplineFunction HermiteFinder(Waypoint PointOne, Waypoint PointTwo);

	SplineFunction QuadraticFinder(Waypoint pointOne, Waypoint pointTwo);

	SplineFunction SquareRootFinder(Waypoint pointOne, Waypoint pointTwo);

	Spline GenerateSpline(Path ThePath);

	Segment GenerateSegment(SplineFunction function, double Jerk);

	void interpolater(Waypoint &p1, Waypoint &p2, bool R2L);

	Trajectory TrajectoryGeneration(AP::Path GroupOfWaypoints, double Jerk, bool interpolationOveride = false);

	struct TankConfig
	{
		double Time, LeftValue, RightValue, kP, kI, kD, MaxVelocity;
		Trajectory LeftTrajectory, RightTrajectory;
	};

	TankConfig GenerateTankConfig(Trajectory theTraj, double WidthBetweenWheels, double Jerk);

	Path getPoints(Trajectory theTraj);

	template<typename Encoder> void AutoTankDrive(TankConfig& TheTankConfig, Encoder& LeftEncoder, Encoder& RightEncoder, double TimeStep, double WidthBetweenWheels, double Jerk, int Ticks, double CircumferenceOfWheel);
}