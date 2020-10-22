#ifndef __BackgroundModel_MAT_h_included__
#define __BackgroundModel_MAT_h_included__

//***********
// BackgroundModel_MAT
// BackgroundModel_MAT fits the MAT model using probes stored by BackgroundDomain
//***********

#include "stdafx.h"
#include "BackgroundModel.h"
#include "ChipBase.h"


class BackgroundModel_MAT : public BackgroundModel {
	
public:
	BackgroundModel_MAT( ChipBase* i_pChipBase );
	~BackgroundModel_MAT();
	
	virtual const gsl_vector* GetBackgroundCoefficients() const;
	virtual void OutputCoefficients( std::ofstream& i_Outfile ) const;
	virtual double GetBackgroundPrediction( ProbeBase* i_pProbe ) const;
	double GetBackgroundPredictionStandardError( ProbeBase* i_pProbe ) const;
	
	double GetRSquaredLogScale() const;
	double GetRSquaredOriginalScale() const;
	double GetSigmaHat() const;
	
private:
		void FitModel( ChipBase* i_pChipBase );
		
		gsl_vector* m_pMATCoefficients;
		double m_RSquaredLogScale;
		double m_RSquaredOriginalScale;
		double m_SigmaHat;
		gsl_matrix* m_pCov;
	
};

#endif

