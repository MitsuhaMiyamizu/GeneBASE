#ifndef __DB_VECTOR_h_included__
#define __DB_VECTOR_h_included__

#include <gsl/gsl_statistics_double.h>

#define MISSING -1e6

template <class T>
T partial_sort(T *vec, int len, double quantile) { // No missing data
	
	int left, right, i, j;
	int pos = int(len * quantile);
	if (pos == len)
		pos = len - 1;

	left = 0; 
	right = len - 1;

	T tmp, pivot;
	while(right > left) {
		pivot = vec[right];
		i = left - 1; 
		j = right;
		while (true) {
			while (vec[++i] < pivot);
			while ((j > 0) && vec[--j] > pivot);
			if (i >= j) break;
			
			// swap
			tmp = vec[i];
			vec[i] = vec[j];
			vec[j] = tmp;
		}

		tmp = vec[i];
		vec[i] = vec[right];
		vec[right] = tmp;

		// now i is the postion where vec[j] <= vec[i] if j < i,
		// and vec[j] >=  vec[i] if j > i
		if (i == pos)
			break;
		else if (i > pos) 
			right = i - 1;
		else // if (i < pos)
			left = i + 1;
	}
	return vec[pos];
}



class DB_VECTOR {
public:
	int MaxLen, Len;
	double *Data; 
	double m_Mean;
	
	DB_VECTOR(int max_len) {
		MaxLen = max_len;
		Data = new double[MaxLen];
		Reset();
	}

	~DB_VECTOR() {
		if (Data) {
			delete [] Data;
			Data = NULL;
		}
	}

	void Reset() {
		Len = 0;
		m_Mean = MISSING;
	}
	//void SetTrimmed2Missing(int pct);

	void Add(double value) {
		if (Len == MaxLen) {
		//WritePrompt("  DB_VECTOR error: length limit reached", false, RED);
		} else
			Data[Len ++] = value;
		}

	//double GetTrimmedMean(int pct) { // 12/10/05
	//SetTrimmed2Missing(pct);
	//return GetMean();
	//}

	double GetMean() {
		if (m_Mean == MISSING) {// set as MISSING at Reset()
			//m_Mean = vec_mean(Data, Len);
			m_Mean = gsl_stats_mean ( Data, 1, Len );

		}
		return m_Mean;
	}

	double GetQuantile(double quantile) {
		// 03/03: if Len is <= 5 and quantile 0.8 (for CheckSingleOutlier), 
		// the largest number is returned.

		if (Len > 0)
			return partial_sort(Data, Len, quantile);
		else
			return MISSING;
	}
		
	double GetMedian() {
		return GetQuantile(.5);
	}

	double GetVarAboutMedian () {
		return gsl_stats_variance_m(Data, 1, Len, this->GetMedian() );
	}

	double GetVar() {
		this->GetMean();
		return gsl_stats_variance_m (Data, 1, Len, m_Mean);
	}

	//void Sort() { //10/15/05
		//qsort((void*)Data, Len, sizeof(double), compare_double);
	//}

};
#endif
