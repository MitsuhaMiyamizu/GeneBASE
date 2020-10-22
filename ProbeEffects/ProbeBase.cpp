#include "stdafx.h"
#include "ProbeBase.h"
#include <gsl/gsl_math.h>

ProbeBase::ProbeBase( const unsigned long& i_Id, const string& i_Sequence  ) : m_ProbeId( i_Id ), m_Sequence( i_Sequence ), 
	m_x( -1 ), m_y( -1 ), m_ExcludeProbe(false), m_pCrossHybInfo( NULL )
{
}

ProbeBase::~ProbeBase(void)
{
	if ( this->m_pCrossHybInfo != NULL ) {
		delete m_pCrossHybInfo;
	}
}

ProbeBase::ProbeBase(const ProbeBase& p) : m_ProbeId( p.m_ProbeId ), m_Sequence( p.m_Sequence ), 
	m_x( p.m_x ), m_y( p.m_y ), m_ExcludeProbe( p.m_ExcludeProbe )
{
	if ( p.m_pCrossHybInfo == NULL ) {
		m_pCrossHybInfo = NULL;
	} else {
		m_pCrossHybInfo = new map< string, pair<int, int> >( *p.m_pCrossHybInfo );
	}
}
ProbeBase& ProbeBase::operator= (const ProbeBase& p)
{
	if (this == &p) return *this;   // Gracefully handle self assignment
	// Put the normal assignment duties here...
	this->m_ProbeId = p.m_ProbeId;
	this->m_Sequence = p.m_Sequence;
	this->m_x = p.m_x;
	this->m_y = p.m_y;
	this->m_ExcludeProbe = p.m_ExcludeProbe;
	if ( p.m_pCrossHybInfo == NULL ) {
		this->m_pCrossHybInfo = NULL;
	} else {
		this->m_pCrossHybInfo = new map< string, pair<int, int> >( *p.m_pCrossHybInfo );
	}
	
	return *this;
}


void ProbeBase::CopyInfo( ProbeBase* i_pProbeBase )
{	
	pair<int, int> posPair = i_pProbeBase->GetPosition();
	this->SetPosition( posPair.first, posPair.second );
	if ( i_pProbeBase->m_pCrossHybInfo != NULL ) {
		this->m_pCrossHybInfo = new map< string, pair<int, int> >( *(i_pProbeBase->m_pCrossHybInfo) );	
	}
}


bool ProbeBase::SetPosition( const int& i_x, const int& i_y ) {
	this->m_x = i_x;
	this->m_y = i_y;
	return true;
}

string ProbeBase::GetSequence() const
{
	return this->m_Sequence;
}

int ProbeBase::GetLength() const
{
	return this->m_Sequence.length();
}

pair<int, int> ProbeBase::GetPosition() const
{
	return make_pair(this->m_x, this->m_y);
}

int ProbeBase::GetGcCount() const
{
	int gcCount = 0;
	int length = this->m_Sequence.length();
	for ( int i = 0; i < length; ++i ) {
		if ( m_Sequence[i] == 'G' || m_Sequence[i] == 'C' ) {
			++gcCount;
		}
	}
	return gcCount;
}

double ProbeBase::GetNucleotideContent( const char& i_Char ) const
{	
	double content = 0;
	int length = this->m_Sequence.length();
	for ( int i = 0; i < length; ++i ) {
		if ( this->m_Sequence[i] == i_Char ) {
			++content;
		}
	}
	return content;
}

//
// Assume that the position is between 0 and this->m_Length
//
double ProbeBase::GetNucleotidePositionIndicator( const char& i_Char, const int& i_Position ) const
{
	if ( this->m_Sequence[i_Position] == i_Char ) {
		return 1;
	} else {
		return 0;
	}	
}

void ProbeBase::SetMATDesignVector( gsl_vector* o_pDesignVector, const bool& i_Binned, const bool& i_UseCopyNumber )
{
	char nucs[4] = {'A','C','G','T' };
	double numA = this->GetNucleotideContent( 'A' );
	double numC = this->GetNucleotideContent( 'C' );
	double numG = this->GetNucleotideContent( 'G' );
	double numT = this->GetNucleotideContent( 'T' );
	int paramCounter = 0;
	if ( i_Binned == false ) { // Including this parameter, design matrix will be of full rank
		gsl_vector_set (o_pDesignVector, 0, numT);
		++paramCounter; // Increments the counter
	}
	for ( int i = 0; i < 25; ++i ) { // probes have length 25
		for ( int j = 0; j < 3; ++j ) {
			gsl_vector_set( o_pDesignVector, (i*3 + j + paramCounter), this->GetNucleotidePositionIndicator( nucs[j], i ) );
		}
	}
	gsl_vector_set( o_pDesignVector, 75 + paramCounter, gsl_pow_2(numA));
	gsl_vector_set( o_pDesignVector, 76 + paramCounter, gsl_pow_2(numC));
	if ( i_Binned == false ) { // Including these parameters, design matrix will be of full rank
		gsl_vector_set( o_pDesignVector, 77 + paramCounter, gsl_pow_2(numG));
		gsl_vector_set( o_pDesignVector, 78 + paramCounter, gsl_pow_2(numT));
	}
}


int h(int i_L, int i_J ) {
	if ( i_L == 0 ) {
		return i_J;
	} else if ( i_L == 1 ) {
		return static_cast<int>( gsl_pow_2( i_J ) );
	} else if ( i_L == 2 ) {
		return static_cast<int>( gsl_pow_3( i_J ) );
	} else if ( i_L == 3) {
		int retVal = static_cast<int>( gsl_pow_3( i_J - 7) );
		if ( retVal > 0 ) {
			return retVal;
		} else {
			return 0; 
		}
	} else if ( i_L == 4 ) {
		int retVal = static_cast<int>( gsl_pow_3( i_J - 15 ) );
		if ( retVal > 0 ) {
			return retVal;
		} else {
			return 0;
		}
	} else {
		return 0;
	}
}

void ProbeBase::AddCrossHybInfo( const string& i_ChTcId, const int& i_Mismatches, const int& i_Indels )
{
	if ( this->m_pCrossHybInfo == NULL ) {
		m_pCrossHybInfo = new map< string, pair<int, int> >;
		m_pCrossHybInfo->insert( make_pair(i_ChTcId, make_pair(i_Mismatches, i_Indels) ) );
	} else {
		map< string, pair<int, int> >::iterator itFind = m_pCrossHybInfo->find( i_ChTcId );
		if ( itFind != m_pCrossHybInfo->end() ) {
			int mm = itFind->second.first;
			int ind = itFind->second.second;
			if ( i_Mismatches + i_Indels < mm + ind ) {
				itFind->second = make_pair( i_Mismatches, i_Indels ); // store the min edit distance btwn a probe and chTc.  
			}
		} else {
			m_pCrossHybInfo->insert( make_pair(i_ChTcId, make_pair(i_Mismatches, i_Indels) ) );
		}
	}
}

unsigned long ProbeBase::GetProbeId() const
{
	return this->m_ProbeId;
}
bool ProbeBase::IsCrossHyb() const
{
	if ( this->m_pCrossHybInfo == NULL || this->m_pCrossHybInfo->size() == 0) {
		return false;
	} else {
		return true;
	}
}
int ProbeBase::GetNumChTcIds() const
{
	if ( this->m_pCrossHybInfo == NULL ) {
		return 0;
	} else {
		return this->m_pCrossHybInfo->size();
	}
}

const map< string, pair<int, int> >* ProbeBase::GetCrossHybInfo() const
{
	return this->m_pCrossHybInfo;
}


void ProbeBase::SetExcludeProbe( const bool& i_val )
{
	this->m_ExcludeProbe = i_val;
}

