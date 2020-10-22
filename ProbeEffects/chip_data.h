// CHIP_DATA stores data relevant to each affy CEL file.  
// The class DataDomain manages the CHIP_DATA objects:
// creation, reading in data, normalization, destruction

#ifndef __chip_data_h_included__
#define __chip_data_h_included__

#include "Probe.h"
#include "stdafx.h"
#include "ChipBase.h"

class CHIP_DATA : public ChipBase {
	friend class DataDomain;
	friend class ChipChipDomain;

public:
	CHIP_DATA( const string& i_Id, const int& i_NumCells, const int& i_CellDim );
	virtual ~CHIP_DATA();
	virtual bool SetGCBgd( );
	double GetBkgSubtractedIntensity( Probe* i_pProbe, const gsl_vector* i_pBeta = NULL );
	double GetGCBgdPValue( Probe* i_pProbe );
	double GetGCBgd( Probe* i_pProbe );
	void SetGCBinBackgroundData();
	void SetProbeNormalizationConstant( const gsl_vector* i_pBeta );
	void SetProbeNormalizationConstant();
	double GetBackgroundCorrectedNormalizedIntensity( Probe* i_pProbe );
	
	void SetBackgroundSubtractedNormalizedIntensities( gsl_vector* i_pBeta, const double& i_ScalingConstant );
	void OutputAllBkgdCorrectNormProbes( const gsl_vector* i_pBeta, const double& i_NormalizationConstant, const double& i_SigmaHat, gsl_matrix* i_pCov );
	void OutputAllBkgdCorrectNormProbes();
	//void OutputSelectedBkgdCorrectNormProbes( gsl_vector* i_pBeta, const double& i_NormalizationConstant, const double& i_SigmaHat, gsl_matrix* i_pCov );
	
 private:
		double m_ProbeNormalizationConstant;
	vector< vector<unsigned short> > m_BackgroundData;
	map<int, double> m_VarGCBgd;	// Stores the variance of background probes with GC content between 0 and 25
	void OutputBkgdCorrectNormProbesHeader( ofstream& i_Outfile );
	void OutputBkgdCorrectNormProbesLine( Probe* i_pProbe, gsl_vector* i_pDesignVec, 
												gsl_matrix* i_pCov, 
												gsl_vector* i_pBeta, const double& i_NormalizationConstant, 
												const double& i_SigmaHat, const string& i_ProbesetId, 
												const string& i_ProbesetAnnotationType, const string& i_TcId, 
												ofstream& i_Outfile);
};
#endif
