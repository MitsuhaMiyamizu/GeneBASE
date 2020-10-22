#ifndef __ChipBase_h_included__
#define __ChipBase_h_included__

#include "CELFileData.h"
#include "BackgroundModel.h"


class ChipBase
{
	friend class ChipChipDomain;
	friend class CrossHybDomain; // for quantile normalization
	friend class AnalysisDomain; // for quantile normalization
	
public:
		ChipBase( const string& i_Id, const int& i_NumCells, const int& i_CellDim );
	virtual ~ChipBase();
	bool ReadBinaryCel(const char *file, const bool& i_ReturnCelData = false, affxcel::CCELFileData** o_pCel = NULL);
	bool ModelBackground();
	const BackgroundModel* GetBackgroundModel() const;
	unsigned short GetUnnormalizedIntensity( ProbeBase* i_pProbe ) const;
	unsigned short GetUnnormalizedIntensity( int x, int y ) const;
	string GetId() const;
	string GetFullname() const;
	int GetNumCells() const;
	virtual bool SetGCBgd( ) = 0;
	unsigned short GetMedianGCBkgIntensity( unsigned short i_Val );
	unsigned short GetMedianGCBkgIntensity( ProbeBase* i_pProbe );
	//void ScaleUnnormalized(); 
	//double GetNormalizationConstant() const;
	void GetMATParameterEstimate( gsl_vector* o_pBeta, gsl_vector* o_pLogPM, double& o_Chisq,
		vector<ProbeBase*>::iterator i_ItProbeBegin, 
		vector<ProbeBase*>::iterator i_ItProbeEnd);

protected:
	unsigned short* m_pInten;
	map<unsigned short, unsigned short > m_MedianGCBgd;
	string m_Id;
	int m_NumCells;
	int m_CellDim;
	BackgroundModel* m_pBackgroundModel;
	
	//double m_NormalizationConstant;
};

#endif
