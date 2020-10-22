#ifndef __my_math_h_included__
#define __my_math_h_included__
#include <map>

class Probeset;

struct FLOAT_POINT {
public:
	float x, y;
};

struct INT_POINT {
public:
	int x, y;
};

struct DOUBLE_POINT {
public:
	double x, y;
};

struct INFO {
	float m_Value;
	float m_ValueStd;
	int m_Exclude; // -1 for array outlier, 1 when only the first probe is
		// a single outlier, 2 when only the 2nd probe is a single outlier ...
};



int compare_float_point( const void *arg1, const void *arg2 );
int compare_int_point( const void *arg1, const void *arg2 );

double InnerProduct(gsl_vector* v1, gsl_vector* v2, int len, 
					gsl_vector* exclude = NULL, bool mean_corrected = false);

bool Same(gsl_vector* v1, gsl_vector* v2, int len, gsl_vector* exclude = NULL, double delta = 1e-6);


struct compare_probeset_pair
{
	bool operator()( std::pair<Probeset*, double> p1, std::pair<Probeset*, double> p2 ) const 
	{
		if (p1.second > p2.second) {
			return true;
		} else if (p1.first > p2.first) {
			return true;
		} else {
			return false;
		}
	}
};

class LSE_FINDER {
public:
	/*
class LSE_FINDER: importing a pivoting point and an array of points, 
fits a least square line passing the pivoting point; used at tails 
of the "invariant set" (together with running medians for the middle 
part of the "invariant set") to approximate a smoothing spline fit. 
Can be ignored since Eric/Byron already did smoothing spline in C.
*/

	#define MAX_POINT 10000 // increase from 2000 on 10/27/05, but best make it 
	// dynamically allocated
	// fit a line passing (x0, y0)

	float X[MAX_POINT];
	float Y[MAX_POINT];

	int Num;
	double x0, y0; 
	double Alpha, Beta; // the intercept and the slope
	double MeanX, MeanY;
	double StdBeta;
	bool HaveAlpha; // intercept

	void Reset(double x, double y, bool have_alpha = false);

	bool AddPoint(double x, double y);

	void GetMean();

	void GetSlope();
			
	double Predict(double x);
	
};

#endif
