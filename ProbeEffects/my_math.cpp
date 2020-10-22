#include "stdafx.h"
#include "my_math.h"
#include <valarray>
#include <utility>




int compare_float_point( const void *arg1, const void *arg2 ) // only ordered by x
{
	float a = ((FLOAT_POINT *) arg1)->x;
	float b = ((FLOAT_POINT *) arg2)->x;

	if (a > b) 
		return 1;
	else if (a < b)
		return -1;
	else 
		return 0;
}

int compare_int_point( const void *arg1, const void *arg2 ) // only ordered by x
{
	int a = ((INT_POINT *) arg1)->x;
	int b = ((INT_POINT *) arg2)->x;

	return (a - b);
		
}

double InnerProduct(gsl_vector* v1, gsl_vector* v2, int len, gsl_vector* exclude, bool mean_corrected) {
	
	double v1_mean = 0;
	double v2_mean = 0;
	
	if (mean_corrected) {
		
		for (int i = 0; i < len; i ++) {
			if ( exclude == NULL || !gsl_vector_get(exclude, i) ) {
				v1_mean += gsl_vector_get( v1, i);
				v2_mean += gsl_vector_get(v2, i);
			}
		}
		v1_mean /= len;
		v2_mean /= len;
	} 
	
	double sum = 0;
	for (int i = 0; i < len; i ++) {
		if (exclude == NULL || !gsl_vector_get(exclude, i) )
			sum += (gsl_vector_get(v1, i) - v1_mean) * (gsl_vector_get(v2,i) - v2_mean);
	}
	return sum;
}


bool Same(gsl_vector* v1, gsl_vector* v2, int len, gsl_vector* exclude, double delta) {
	// compare two vectors to see if they are close enough
	double sum = 0;
	for (int i = 0; i < len; i ++) {
		if (exclude == NULL || !gsl_vector_get(exclude, i) )
			sum += fabs(gsl_vector_get(v1, i) - gsl_vector_get(v2, i));
	}
	sum /= len;
	return (sum < 1e-6);
}


//bool compare_probeset_pair( std::pair<Probeset*, double>& p1, std::pair<Probeset*, double>& p2 )
//{
//	return (p1.second > p2.second );
//}

void LSE_FINDER::Reset(double x, double y, bool have_alpha) {
		x0 = x;
		y0 = y;
		Num = 0;
		HaveAlpha = have_alpha;
}

bool LSE_FINDER::AddPoint(double x, double y) {
		if (Num == MAX_POINT) {
			//sprintf(prompt_str, "    Exceeds the limit %d in LSE_FINDER; use this many data points to fit slope",
			//	MAX_POINT);
				
			//WritePrompt(prompt_str, false, RED);
			return false;
		}
		X[Num] = float(x - x0);
		Y[Num] = float(y - y0);
		Num ++;
		return true;
		
	}

void LSE_FINDER::GetMean() {
		MeanX = 0, MeanY = 0;

		if (HaveAlpha) {
			for (int i = 0; i < Num; i ++) {
				MeanX += X[i];
				MeanY += Y[i];
			}
			MeanX /= Num;
			MeanY /= Num;
		}

	}


void LSE_FINDER::GetSlope() {
		GetMean();

		double sum_xy = 0, sum_xx = 0;
		int i;
		for (i = 0; i < Num; i ++) {
			sum_xy += (Y[i] - MeanY) * (X[i] - MeanX);
			sum_xx += std::pow(X[i] - MeanX, 2);
		}
		Beta = sum_xy / sum_xx;
		Alpha = MeanY - Beta * MeanX;

		// find sigma
		double RSS = 0;
		for (i = 0; i < Num; i ++) {
			double fitted = Alpha + X[i] * Beta;
			RSS += std::pow(Y[i] - fitted, 2);
		}
		double sig_sqr = RSS / (Num - 2);
		StdBeta = std::sqrt(sig_sqr / sum_xx); // Jennrich p 31

		
	}

double LSE_FINDER::Predict(double x) {
		return y0 + (x - x0) * Beta + .5;
}

