#ifndef __ProbeBase_h_included__
#define __ProbeBase_h_included__

class ChipChipProbe;

class ProbeBase
{
public:
	ProbeBase( const unsigned long& i_Id, const string& i_Sequence );
	virtual ~ProbeBase();
	ProbeBase(const ProbeBase& p);	// Copy constructor, due to pointer m_pCrossHybInfo
	ProbeBase& operator= (const ProbeBase& p); // Assignment= operator

	pair<int, int> GetPosition() const;
	string GetSequence() const;
	bool SetPosition(const int& i_x, const int& i_y);
	int GetGcCount() const;
	int GetLength() const;
	void CopyInfo( ProbeBase* i_pProbeBase );
	void AddCrossHybInfo( const string& i_ChTcId, const int& i_Mismatches, const int& i_Indels );
	int GetNumChTcIds() const;
	const map< string, pair<int, int> >* GetCrossHybInfo() const;
	//
	// The number of nucleotides taking value i_Char
	//
	double GetNucleotideContent( const char& i_Char ) const;
	
	//
	// Assume that the position is between 0 and this->m_Length
	//
	double GetNucleotidePositionIndicator( const char& i_Char, const int& i_Position ) const;
	virtual void SetMATDesignVector( gsl_vector* o_pDesignVector, const bool& i_Binned = false, const bool& i_UseCopyNumber = false);
	unsigned long GetProbeId() const;
	bool IsCrossHyb() const;
	//bool ExcludeProbe() const;
	void SetExcludeProbe( const bool& i_val);
	
protected:
	unsigned long m_ProbeId;
	string m_Sequence;
	int m_x; //The x and y coordinates 
	int m_y; //of the probe's location on the array
	bool m_ExcludeProbe; //Whether to exclude the probe from downstream analysis
	map< string, pair< int, int > >* m_pCrossHybInfo;
	
	
};

#endif
