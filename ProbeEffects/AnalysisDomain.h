#ifndef __AnalysisDomain_h_included__
#define __AnalysisDomain_h_included__

//**************
// AnalysisDomain has 3 major functions
// 
// 1.  Output background-corrected probe intensities
// 2.  Output gene-level expression indexes
// 3.  Output cross-hyb info
// 4.  Output presence/absence calls
//**************

#include "stdafx.h"
#include "BackgroundModel.h"
#include "TranscriptCluster.h"
#include "chip_data.h"

class AnalysisDomain {

public:
	static AnalysisDomain* GetInstance();
	~AnalysisDomain(void);
	
	bool ModelBackground();
	bool SummarizeExpression();
	bool OutputCrossHybInfo();
	bool ModelPresenceAbsence();
	
private:
	
	AnalysisDomain();
	static bool m_InstanceFlag;
	static AnalysisDomain* m_pTheOne;
	
	bool GetMeanQuantiles();
	void ComputeQuantileNormalizedIntensities( ChipBase* i_pChipBase );
	void DeleteMeanQuantiles(); // Deletes m_pMeanQuantiles;
	void GetIntensityMatrix( gsl_matrix* o_pInten, TranscriptCluster* i_pTranscriptCluster, 
							 vector<CHIP_DATA*>& i_Chips, const bool& i_ExcludeCrossHyb );
	void GetIntensityVector(  gsl_vector* o_pInten, Probe* i_pProbe, vector<CHIP_DATA*>& i_Chips);
	
	double* m_pMeanQuantiles;
	
	
};




#endif

