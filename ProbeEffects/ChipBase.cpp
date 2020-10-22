#include "stdafx.h"
#include "ChipBase.h"
#include "CELFileData.h"
using namespace affxcel;
#include "Probe.h"
#include "DB_VECTOR.h"
#include "TranscriptClusterDomain.h"
#include "BackgroundModel_MAT.h"
#include "BackgroundModel_MedianGC.h"
#include "BackgroundModel_None.h"
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_multifit.h>


ChipBase::ChipBase( const string& i_Id, const int& i_NumCells, const int& i_CellDim )
{
	m_Id = i_Id;
	m_NumCells = i_NumCells;
	m_CellDim = i_CellDim;
	m_pBackgroundModel = NULL;
	m_pInten = new unsigned short[ i_NumCells ];
	for ( int i = 0; i < i_NumCells; ++i ) {
		m_pInten[i] = 0; // important for non-square arrays
	}
}

ChipBase::~ChipBase()
{
	delete [] m_pInten;
	m_pInten = NULL;
	delete m_pBackgroundModel;
	m_pBackgroundModel = NULL;
}

//*****************
// Reads Cel file intensities into a unsigned short[ i_NumCells ].  
// We check that m_NumCells is equal to the Cel file's numcells.  
// If the user specifies i_ReturnCelData = true, then the user is responsible for 
// freeing the memory of pCel
//*****************
bool ChipBase::ReadBinaryCel(const char* i_pFile, const bool& i_ReturnCelData, CCELFileData** o_pCel) {
	Parameters* pParameters = Parameters::GetInstance();
	CCELFileData* pCel = new CCELFileData;
	CELFileEntryType entry;

	if (!pCel) {
		Parameters::LOGFILE << "ERROR:  Allocate memory failed in chip_data::ReadBinaryCel()" << endl;
		return false;
	}
	pCel->SetFileName(i_pFile);

	if (!pCel->Read()) {
		Parameters::LOGFILE << "ERROR:  Failed to read Binary CEL file" << endl;
		Parameters::LOGFILE << "Reason:  " << pCel->GetError() << endl;
		delete pCel;
		return false;
	}

	if ( pCel->GetCols()*pCel->GetRows() != this->m_NumCells ) {
		Parameters::LOGFILE << "ERROR:  Cel file number of cells (" << pCel->GetRows() * pCel->GetCols() << ") " <<
		" does not match the numcells for array type " << pParameters->GetArrayType() << endl;
		exit(1);
	}
	
	for (int iy=0; iy< pCel->GetRows(); iy++) { // follow order in text cel file
		for (int ix=0; ix< pCel->GetCols(); ix++) {
			pCel->GetEntry(ix, iy, entry);
			// DChip code
			int idx = Parameters::IDX(ix, iy, this->m_CellDim );
			unsigned short val = static_cast<unsigned short>( entry.Intensity + .5 ); // Add .5 for conversion to int
			m_pInten[idx] = val; 
		}
	}
	if ( i_ReturnCelData ) { 
		*o_pCel = pCel;	
	} else {
		delete pCel;
	}
	return true;
}

bool ChipBase::ModelBackground( )
{
	Parameters* pParameters = Parameters::GetInstance();
	string method = pParameters->GetModelBackgroundMethod();
	if ( method == "mat" ) {
		if ( m_pBackgroundModel == NULL ) { // Only compute the background model one time
			this->m_pBackgroundModel = new BackgroundModel_MAT( this );
		}
	} else if ( method == "median_gc") {
		if ( m_pBackgroundModel == NULL ) {
			this->m_pBackgroundModel = new BackgroundModel_MedianGC( this );
		}
	} else if ( method == "none" ) {
		if ( m_pBackgroundModel == NULL ) {
			this->m_pBackgroundModel = new BackgroundModel_None( this );
		}
	} else {
		Parameters::LOGFILE << "Error:  unknown background model method in ChipBase" << endl;
		exit(1);
	}
	return true;
}

const BackgroundModel* ChipBase::GetBackgroundModel() const
{
	return m_pBackgroundModel;	
}

unsigned short ChipBase::GetUnnormalizedIntensity( ProbeBase* i_pProbe ) const
{
	pair<int, int> posPair = i_pProbe->GetPosition();
	return this->m_pInten[ Parameters::IDX(posPair.first, posPair.second, this->m_CellDim )  ];
}

unsigned short ChipBase::GetUnnormalizedIntensity( int i_x, int i_y ) const
{
	return this->m_pInten[ Parameters::IDX( i_x, i_y, this->m_CellDim ) ];
}

string ChipBase::GetId() const
{
	return this->m_Id;
}

int ChipBase::GetNumCells() const
{
	return this->m_NumCells;
}



unsigned short ChipBase::GetMedianGCBkgIntensity( ProbeBase* i_pProbe )
{
	map<unsigned short, unsigned short>::const_iterator it = this->m_MedianGCBgd.find( i_pProbe->GetGcCount() );
	unsigned short gcbgd;
	if (it == this->m_MedianGCBgd.end() ) {
		gcbgd = 0; 
	} else {
		gcbgd = it->second;
	}
	return gcbgd;
}

unsigned short ChipBase::GetMedianGCBkgIntensity( unsigned short i_Val )
{
	map<unsigned short, unsigned short>::const_iterator it = this->m_MedianGCBgd.find( i_Val );
	unsigned short gcbgd;
	if (it == this->m_MedianGCBgd.end() ) {
		gcbgd = 0; 
	} else {
		gcbgd = it->second;
	}
	return gcbgd;
}


//double ChipBase::GetNormalizationConstant() const
//{
//	return this->m_NormalizationConstant;
//}

void ChipBase::GetMATParameterEstimate( gsl_vector* o_pBeta, gsl_vector* o_pLogPM, double& o_Chisq,
		vector<ProbeBase*>::iterator i_ItProbeBegin, 
		vector<ProbeBase*>::iterator i_ItProbeEnd)
{
	Parameters* pParameters = Parameters::GetInstance();
	int numProbes = distance( i_ItProbeBegin, i_ItProbeEnd );
	int numParameters = pParameters->GetNumParameters_MAT();
	gsl_matrix*	pDesignMatrix = gsl_matrix_calloc( numProbes, numParameters );	
	int selectedProbeCtr = 0;
	for (; i_ItProbeBegin != i_ItProbeEnd; ++i_ItProbeBegin ) {
		ProbeBase* pProbe = *i_ItProbeBegin;
		gsl_vector_view designRowView = gsl_matrix_row (pDesignMatrix, selectedProbeCtr);
		pProbe->SetMATDesignVector( (&designRowView.vector) );
		double pm = static_cast<double>( this->GetUnnormalizedIntensity( pProbe ) );
		double minVal = 1;
		double logPm = log((max(pm, minVal)) );
		gsl_vector_set( o_pLogPM, selectedProbeCtr, logPm);
		++selectedProbeCtr;
	}
	gsl_multifit_linear_workspace* pWorkspace = gsl_multifit_linear_alloc (numProbes, numParameters);
	gsl_matrix* pCov = gsl_matrix_alloc( numParameters, numParameters );
	gsl_multifit_linear (pDesignMatrix, o_pLogPM, o_pBeta, pCov, (&o_Chisq), pWorkspace); 
	
	// Free memory
	gsl_multifit_linear_free( pWorkspace );
	gsl_matrix_free( pDesignMatrix );
	gsl_matrix_free( pCov );
}








