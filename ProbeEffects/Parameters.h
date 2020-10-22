#ifndef __Parameters_h_included__
#define __Parameters_h_included__

#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>

using namespace std;

enum AnnotationConf {CORE, EXTENDED, FULL, FREE, AMBIGUOUS, NONE};

class Parameters
{
public:
	~Parameters(void);
	static Parameters* GetInstance();
	bool LoadParameters( const std::string& i_ParameterInfile );

	// Log file used by all
	static std::ofstream LOGFILE;
	static std::string SEPARATOR;
	static int PROBELENGTH;

	// *****************************
	// Files and folders
	// *****************************
	int GetNumArrays() const;
	// Exon, gene, glue array
	string GetArrayCelFolder() const;
	pair<set<string>::iterator, set<string>::iterator> GetCelFilesIterators();
	string GetPgfFile() const;
	string GetClfFile() const;
	string GetProbesetAnnotation() const;
	string GetCrossHybAnnotation() const;
	AnnotationConf GetAnnotationGeneEstimates() const;
	// Promoter
	pair<set<string>::iterator, set<string>::iterator> GetPromoterCelFilesIterators();
	string GetPromoterBpmapFile() const;
	string GetPromoterArrayCelFolder() const;
	// ChipChip / Tiling
	pair< multimap<string, string>::iterator, multimap<string, string>::iterator> 
		GetChipChipCelFileIterators( const string& i_CharVal );
	pair< multimap<string, string>::iterator, multimap<string, string>::iterator> 
		GetChipChipCelFileIterators();
	pair< map<string, string>::iterator, map<string, string>::iterator > GetChipChipBpmapIterators();
	bool GetChipChipBpmapFile( const string& i_Type, string& i_File ) const;
	string GetChipChipCelFolder() const;
	string GetChipChipBpmapFolder() const;
	

	// promoter input
	static string PROMOTER_BPMAP_FOLDER;

	// get array parameters
	int GetNumCells() const;
	int GetCellDim() const;
	
	// methods for modeling background
	string GetModelBackgroundMethod() const;
	bool TrainModel() const;
	bool OutputCrossHyb() const;
	int GetNumProbesMAT() const;
	string GetArrayType() const;
	bool IncludeCopyNumberMATParameter() const;
	int GetNumParameters_MAT() const;
	
	// presence/absence
	bool CalculatePresenceAbsence() const;
	string GetPresenceAbsenceOutfile() const;
	
	string GetWorkingDirectory() const;
	string GetMATTrainingProbeType() const;
	bool ChipChipTargetsReverse() const;
	string GetOutputModelFit() const;
	bool OutputAllBkgdCorrectNormProbes() const;
	string GetBkgdCorrectNormProbesFile() const;
	bool OutputCelFormat() const;
	string GetNormalizationMethod() const;
	
	//***
	// Summarize gene-level expression
	bool SummarizeExpression() const;
	string GetSummaryMethod() const;
	double GetCorrelationFilterMinCorrelation() const;
	double GetCorrelationFilterMinSd() const;
	double GetCorrelationFilterOutputMasked() const;
	int GetCorrelationFilterEditDistance() const;
	
	string GetCrossHybOutputFile() const;
	string GetCrossHybDir() const;
	string GetTcEstimExprFile() const;
	string GetBkgdCoefFile() const;
	string GetCrossHybType() const;

	// General methods
	int RandomNumber( const int& i_Size );
	gsl_rng* GetRng();
	static unsigned long IDX(int i, int j, int i_CellDim); 
	static pair<int, int> LocsFromIDX(int i_IDX, int i_CellDim );
	static bool MOSTLY_POSITIVE( gsl_vector* i_pVec, int i_Size );
	static string REVERSE_COMPLEMENT( string i_Sequence );
	static string SWAP_ORDER( const string& i_Sequence );
	bool GetTcEstimExpr( map< unsigned long, vector<double>* >& o_GeneExpression );

private:
	static bool m_InstanceFlag;
	static Parameters* m_pTheOne;
	Parameters();
	int m_NUMCELLS;
	int m_CELLDIM;

	// exon
	int m_NUMCELLS_EXON_ARRAY;
	int m_CELLDIM_EXON_ARRAY;
	string m_PROBESET_ANNOTATION;
	//string m_CrossHybAnnotation;

	// gene array
	int m_NUMCELLS_GENE_ARRAY;
	int m_CELLDIM_GENE_ARRAY;
	
	// glue array
	int m_NUMCELLS_GLUE_ARRAY;
	int m_CELLDIM_GLUE_ARRAY;

	// Exon, gene, glue array
	set<string> m_CelFiles;
	string m_ARRAY_CEL_FOLDER;
	string m_CLFFILE;
	string m_PGFFILE;
	//AnnotationConf m_ANNOT_GENE_ESTIMATES;

	// chip chip
	map< string, string > m_ChipChipBpmap;
	multimap<string, string > m_ChipChipCelFiles;
	int m_NUMCELLS_CHIP_CHIP_ARRAY;
	int m_CELLDIM_CHIP_CHIP_ARRAY;
	string m_CHIP_CHIP_CEL_FOLDER;
	string m_CHIP_CHIP_BPMAP_FOLDER;
	bool m_ChipChipTargetsReverse;

	// promoter
	set<string> m_PromoterCelFiles;
	string m_PromoterBpmap;
	string m_PROMOTER_CEL_FOLDER;
	int m_NUMCELLS_PROMOTER_ARRAY;
	int m_CELLDIM_PROMOTER_ARRAY;
	
	//***
	// Background modeling
	string m_MODEL_BACKGROUND_METHOD;
	bool m_TrainModel; // Whether to train a background model for parameter estimates
	bool m_OutputCrossHyb; // Whether to output the cross hyb info
	string m_MATInfile; // the file specifying parameter estimates for the MAT model
	string m_MATTrainingProbeType; // The type of probes used to train the MAT model
	int m_NUM_PROBES_MAT;
	int m_NumParameters_MAT;
	string m_ArrayType; // The type of array we are analyzing
	string m_NormalizationMethod; // The type of normalization method
	bool m_IncludeCopyNumberMATParameter;
	//***
	// Presence/absence modeling
	string m_PresenceAbsenceOutfile;
	bool m_CalculatePresenceAbsence;
	
	//***
	// Summarization of gene-level estimates
	bool m_SummarizeExpression; // Whether to compute gene-level estimates
	string m_SummaryMethod; // correlation_filter, selection, or liwong
	double m_CorrelationFilterMinCorrelation;
	double m_CorrelationFilterMinSd;
	double m_CorrelationFilterOutputMasked;
	int m_CorrelationFilterEditDistance;

	string m_WorkingDirectory;
	string m_OutputModelFit;
	bool m_OutputAllBkgdCorrectNormProbes;
	string m_BkgdCorrectNormProbesFile;
	bool m_OutputCelFormat;
	//***
	// Cross hyb modeling
	string m_CrossHybOutputFile;
	string m_CrossHybDir;
	string m_TcEstimExprFile;
	string m_BkgdCoefFile;
	string m_CrossHybType; 

	gsl_rng* m_pRng; // Random number generator, set up when TranscriptClusterDomain is created.  
	const gsl_rng_type **m_AllGslRngTypes;
	const gsl_rng_type* m_pGslRngType;
};

#endif

