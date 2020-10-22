#include "stdafx.h"
#include "Parameters.h"
#include "StringTokenizer.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multifit.h>

bool Parameters::m_InstanceFlag = false;
Parameters* Parameters::m_pTheOne = NULL;
std::ofstream Parameters::LOGFILE;

int Parameters::PROBELENGTH = 25;
unsigned long Parameters::IDX(int i, int j, int i_CellDim)
{
	unsigned long result = j;
	result *= i_CellDim;
	result += i;
	return result;
}

pair<int, int> Parameters::LocsFromIDX( int i_IDX, int i_CellDim )
{
	int j = i_IDX / i_CellDim; // integer division
	int i = i_IDX - j* i_CellDim;
	return make_pair(i,j);
}

#ifndef WIN32
std::string Parameters::SEPARATOR = "/";
#else
std::string Parameters::SEPARATOR = "\\"; 
#endif

// Promoter input
string Parameters::PROMOTER_BPMAP_FOLDER;

Parameters::Parameters(void)
{
	// Set up a random number generator
	this->m_AllGslRngTypes = gsl_rng_types_setup ();
	this->m_pGslRngType = m_AllGslRngTypes[12];
	gsl_rng_env_setup();
	this->m_pRng = gsl_rng_alloc (  this->m_pGslRngType  ) ;

	// Initialize modeling parameters with default values
	m_TrainModel = false;
	m_OutputCrossHyb = false;
	m_SummarizeExpression = false;
	m_MODEL_BACKGROUND_METHOD = "mat";
	m_NUM_PROBES_MAT = 400000; 
	m_NumParameters_MAT = 80;
	m_ArrayType = "exon"; 
	m_MATTrainingProbeType = "background"; // changed 07/03/08
	m_NormalizationMethod = "core_probe_scaling";
	m_IncludeCopyNumberMATParameter = false;
	m_WorkingDirectory = "";
	m_CrossHybType = "none";
	m_SummaryMethod = "selection";
	m_CorrelationFilterMinCorrelation = 0.7;
	m_CorrelationFilterMinSd = 100;
	m_CorrelationFilterEditDistance = 2;
	m_CorrelationFilterOutputMasked = false;
	m_CalculatePresenceAbsence = false;



	// Initialize array parameters with default values
	m_NUMCELLS = 0;
	m_CELLDIM = 0;
	m_NUMCELLS_CHIP_CHIP_ARRAY = 835396;
	m_CELLDIM_CHIP_CHIP_ARRAY = 914;
	m_NUMCELLS_EXON_ARRAY = 6553600;
	m_CELLDIM_EXON_ARRAY = 2560;
	m_NUMCELLS_PROMOTER_ARRAY = 4691556; 
	m_CELLDIM_PROMOTER_ARRAY = 2166;
	m_NUMCELLS_GENE_ARRAY = 1102500;
	m_CELLDIM_GENE_ARRAY = 1050;
	m_NUMCELLS_GLUE_ARRAY = 4545424;
	m_CELLDIM_GLUE_ARRAY = 2132;
	
	
	// Initialize file names with default values
	m_MATInfile = "";
	m_OutputModelFit = "";
	m_CrossHybOutputFile = "";
	m_CrossHybDir = "";
	m_TcEstimExprFile = "";
	m_BkgdCoefFile = "";
	m_CrossHybType = "";
	m_PresenceAbsenceOutfile = "";
	

	// Fitted values and residuals
	m_OutputAllBkgdCorrectNormProbes = false;
	m_BkgdCorrectNormProbesFile = "";
	m_OutputCelFormat = false;

	// Promoter array
	m_PROMOTER_CEL_FOLDER = "";
	
	// ChipChip array
	m_ChipChipTargetsReverse = true;
	m_CHIP_CHIP_CEL_FOLDER = "";
	m_CHIP_CHIP_BPMAP_FOLDER = "";
	
	// Exon, Gene, glue array
	m_ARRAY_CEL_FOLDER = "";
}

Parameters::~Parameters(void)
{
	gsl_rng_free (this->m_pRng);
}

Parameters* Parameters::GetInstance()
{
	if ( !m_InstanceFlag )
	{
		m_pTheOne = new Parameters();
		m_InstanceFlag = true;
		return m_pTheOne;
	}
	else {
		return m_pTheOne;
	}
}


//
// This method will load in all parameter values from the parameter file 
// specified on the command line.  
// 
//
bool Parameters::LoadParameters( const std::string& i_ParameterInfile )
{
	ifstream infile;
	infile.open( i_ParameterInfile.c_str() );
	if (infile.is_open() == false) {
		return false;
	}

	while (!infile.eof() ) {
		string inputLine;
		getline( infile, inputLine);
		
		if ( inputLine == "[log]" ) {
			while ( !infile.eof() ) {
				string inputParameter;
				getline( infile, inputParameter );
				if ( inputParameter == "" )	break; // go to the next "[...]" 
				StringTokenizer tokenizedInput( inputParameter, "\t" );
				string flag = tokenizedInput.nextToken();
				if ( flag == "logfile" ) {
					Parameters::LOGFILE.open( tokenizedInput.nextToken().c_str() );
					if ( Parameters::LOGFILE.is_open() == false ) {
						cout << "ERROR:  Cannot open log file" << endl;
						return false;
					}
				}
			}
		} else if ( inputLine == "[exon_annotation]" ) {
			while ( !infile.eof() ) {
				string inputParameter;
				getline( infile, inputParameter );
				if ( inputParameter == "" )	break; // go to the next "[...]" 
				StringTokenizer tokenizedInput( inputParameter, "\t" );
				string flag = tokenizedInput.nextToken();	
				if ( flag == "probeset_annotation") {
					this->m_PROBESET_ANNOTATION = tokenizedInput.nextToken();
				} else if ( flag == "pgf_file") {
					this->m_PGFFILE = tokenizedInput.nextToken();
				} else if ( flag == "clf_file") {
					this->m_CLFFILE = tokenizedInput.nextToken();
				} else if ( flag == "num_cells" ) {
					this->m_NUMCELLS_EXON_ARRAY = tokenizedInput.nextIntToken();
				} else if ( flag == "cell_dim" ) {
					this->m_CELLDIM_EXON_ARRAY = tokenizedInput.nextIntToken();
				}
			}
		} else if ( inputLine == "[gene_annotation]" ) {
			while ( !infile.eof() ) {
				string inputParameter;
				getline( infile, inputParameter );
				if ( inputParameter == "" )	break; // go to the next "[...]" 
				StringTokenizer tokenizedInput( inputParameter, "\t" );
				string flag = tokenizedInput.nextToken();	
				if ( flag == "pgf_file") {
					this->m_PGFFILE = tokenizedInput.nextToken();
				} else if ( flag == "clf_file") {
					this->m_CLFFILE = tokenizedInput.nextToken();
				} else if ( flag == "num_cells" ) {
					this->m_NUMCELLS_GENE_ARRAY = tokenizedInput.nextIntToken();
				} else if ( flag == "cell_dim" ) {
					this->m_CELLDIM_GENE_ARRAY = tokenizedInput.nextIntToken();
				}
			}
		} else if ( inputLine == "[glue_annotation]" ) {
			while ( !infile.eof() ) {
				string inputParameter;
				getline( infile, inputParameter );
				if ( inputParameter == "" )	break; // go to the next "[...]" 
				StringTokenizer tokenizedInput( inputParameter, "\t" );
				string flag = tokenizedInput.nextToken();	
				if ( flag == "probeset_annotation") {
					this->m_PROBESET_ANNOTATION = tokenizedInput.nextToken();
				} else if ( flag == "pgf_file") {
					this->m_PGFFILE = tokenizedInput.nextToken();
				} else if ( flag == "clf_file") {
					this->m_CLFFILE = tokenizedInput.nextToken();
				} else if ( flag == "num_cells" ) {
					this->m_NUMCELLS_GLUE_ARRAY = tokenizedInput.nextIntToken();
				} else if ( flag == "cell_dim" ) {
					this->m_CELLDIM_GLUE_ARRAY = tokenizedInput.nextIntToken();
				}
			}
		} else if ( inputLine == "[exon_data]" || inputLine == "[gene_data]" || inputLine == "[glue_data]" ) {
			while ( !infile.eof() ) {
				string inputParameter;
				getline( infile, inputParameter );
				if ( inputParameter == "" )	break; // go to the next "[...]" 
				StringTokenizer tokenizedInput( inputParameter, "\t" );
				string flag = tokenizedInput.nextToken();
				if ( flag == "exon_cel_files" || flag == "gene_cel_files" || flag == "glue_cel_files" ) {
					StringTokenizer tokenizedFiles( tokenizedInput.nextToken(), "," );
					while ( tokenizedFiles.hasMoreTokens() ) {
						string geneFile = tokenizedFiles.nextToken();
						this->m_CelFiles.insert( geneFile );
					}
				} else if ( flag == "folder" ) {
					this->m_ARRAY_CEL_FOLDER = tokenizedInput.nextToken();
				}
			}	
		} else if ( inputLine == "[chip_chip_annotation]" ) {
			while ( !infile.eof() ) {
				string inputParameter;
				getline( infile, inputParameter );
				if ( inputParameter == "" )	break; // go to the next "[...]" 
				StringTokenizer tokenizedInput( inputParameter, "\t" );
				string flag = tokenizedInput.nextToken();
				if ( flag == "folder" ) {
					m_CHIP_CHIP_BPMAP_FOLDER = tokenizedInput.nextToken();
				} else if ( flag == "num_cells" ) {
					this->m_NUMCELLS_CHIP_CHIP_ARRAY = tokenizedInput.nextIntToken();
				} else if ( flag == "cell_dim" ) {
					this->m_CELLDIM_CHIP_CHIP_ARRAY = tokenizedInput.nextIntToken();				
				} else {
					string file = tokenizedInput.nextToken();
					this->m_ChipChipBpmap.insert( make_pair(flag, file) );
				}
			}
		} else if ( inputLine == "[chip_chip_data]" ) {
			while ( !infile.eof() ) {
				string inputParameter;
				getline( infile, inputParameter );
				if ( inputParameter == "" )	break; // go to the next "[...]" 
				StringTokenizer tokenizedInput( inputParameter, "\t" );
				string flag = tokenizedInput.nextToken();
				if ( flag == "folder" ) {
					m_CHIP_CHIP_CEL_FOLDER = tokenizedInput.nextToken();
				} else if ( flag == "target_strand" ) {
					if ( tokenizedInput.nextToken() == "reverse" ) {
						this->m_ChipChipTargetsReverse = true;
					} else {
						this->m_ChipChipTargetsReverse = false;
					}
				} else {
					string file = tokenizedInput.nextToken();
					StringTokenizer tokenizedFiles( file, "," );
					int numFiles = tokenizedFiles.countTokens();
					for ( int i = 0; i < numFiles; ++i ) {
						this->m_ChipChipCelFiles.insert( make_pair(flag, tokenizedFiles.nextToken() ) );
					}
				}
			}
		} else if ( inputLine == "[promoter_data]" ) {
			while ( !infile.eof() ) {
				string inputParameter;
				getline( infile, inputParameter );
				if ( inputParameter == "" )	break; // go to the next "[...]" 
				StringTokenizer tokenizedInput( inputParameter, "\t" );
				string flag = tokenizedInput.nextToken();
				if ( flag == "folder" ) {
					m_PROMOTER_CEL_FOLDER = tokenizedInput.nextToken();
				} else if ( flag == "promoter_cel_files" ) {
					StringTokenizer tokenizedFiles( tokenizedInput.nextToken(), "," );
					while ( tokenizedFiles.hasMoreTokens() ) {
						string promoterFile = tokenizedFiles.nextToken();
						this->m_PromoterCelFiles.insert( promoterFile );
					}
				} 
			}
		} else if ( inputLine == "[promoter_annotation]" ) {
			while( !infile.eof() ) {
				string inputParameter;
				getline( infile, inputParameter );
				if ( inputParameter == "" )	break; // go to the next "[...]" 
				StringTokenizer tokenizedInput( inputParameter, "\t" );
				string flag = tokenizedInput.nextToken();
				if ( flag == "folder" ) {
					Parameters::PROMOTER_BPMAP_FOLDER = tokenizedInput.nextToken();
				} else if ( flag == "bpmap" ) {
					this->m_PromoterBpmap = tokenizedInput.nextToken();
				} else if ( flag == "num_cells" ) {
					this->m_NUMCELLS_PROMOTER_ARRAY = tokenizedInput.nextIntToken();
				} else if ( flag == "cell_dim" ) {
					this->m_CELLDIM_PROMOTER_ARRAY = tokenizedInput.nextIntToken();
				}
			}
		} else if ( inputLine == "[output]" ) {
			while ( !infile.eof() ) {
				string inputParameter;
				getline( infile, inputParameter );
				if ( inputParameter == "" )	break; // go to the next "[...]" 
				StringTokenizer tokenizedInput( inputParameter, "\t" );
				string flag = tokenizedInput.nextToken();
				if ( flag == "output_model_fit" ) {
					this->m_OutputModelFit = tokenizedInput.nextToken();
				} else if ( flag == "output_all_bkgd_correct_norm_probes" ) {
					string outputFittedProbes = tokenizedInput.nextToken();
					if ( outputFittedProbes == "true" ) {
						this->m_OutputAllBkgdCorrectNormProbes = true;
					} else if ( outputFittedProbes == "false" ) {
						this->m_OutputAllBkgdCorrectNormProbes = false;
					} else {
						cout << "ERROR:  Invalid argument to output_all_fitted_exon_probes: " << outputFittedProbes << endl;
						return false;
					}
				} else if ( flag == "bkgd_correct_norm_probes_file" ) {
					this->m_BkgdCorrectNormProbesFile = tokenizedInput.nextToken();
				} else if ( flag == "output_cel_format" ) {
					string outputCelFormat = tokenizedInput.nextToken();
					if ( outputCelFormat == "true" ) {
						this->m_OutputCelFormat = true;
					} else if ( outputCelFormat == "false" ) {
						this->m_OutputCelFormat = false;
					} else {
						cout << "ERROR:  Invalid argument to output_cel_format:  " << outputCelFormat << endl;
						return false;
					}
				} else if ( flag == "dir" ) {
					this->m_WorkingDirectory = tokenizedInput.nextToken();
				} else if ( flag == "cross_hyb_output_file" ) {
					this->m_CrossHybOutputFile = tokenizedInput.nextToken();
				} else if ( flag == "cross_hyb_dir" ) {
					this->m_CrossHybDir = tokenizedInput.nextToken();
				} else if ( flag == "tc_estim_expr_file" ) {
					this->m_TcEstimExprFile = tokenizedInput.nextToken();
				} else if ( flag == "bkgd_coef_file" ) {
					this->m_BkgdCoefFile = tokenizedInput.nextToken();
				} else if ( flag == "presence_absence_output_file" ) {
					this->m_PresenceAbsenceOutfile = tokenizedInput.nextToken();
				}
			}
		} else if ( inputLine == "[model]" ) {
			while ( !infile.eof() ) {
				string inputParameter;
				getline( infile, inputParameter );
				if ( inputParameter == "" )	break; // go to the next "[...]" 
				StringTokenizer tokenizedInput( inputParameter, "\t" );
				string flag = tokenizedInput.nextToken();
				if ( flag == "method" ) {
					this->m_MODEL_BACKGROUND_METHOD = tokenizedInput.nextToken();
					if ( m_MODEL_BACKGROUND_METHOD != "mat" && 
						m_MODEL_BACKGROUND_METHOD != "median_gc" &&
						m_MODEL_BACKGROUND_METHOD != "none" ) 
					{
						cout << "ERROR:  Invalid argument passed to method: " << m_MODEL_BACKGROUND_METHOD << endl;
						return false;
					}
				} else if ( flag == "normalization_method") {
					this->m_NormalizationMethod = tokenizedInput.nextToken();
					if ( m_NormalizationMethod != "core_probe_scaling" &&
						m_NormalizationMethod != "quantile" &&
						m_NormalizationMethod != "none" ) 
					{
						cout << "ERROR:  Invalid argument passed to normalization_method: " << m_NormalizationMethod << endl;
						return false;
					}
				} else if ( flag == "train_model" ) {
					string trainModel = tokenizedInput.nextToken();
					if ( trainModel == "true" ) {
						this->m_TrainModel = true;
					} else if ( trainModel == "false" ) {
						this->m_TrainModel = false;
					} else {
						cout << "ERROR:  Invalid argument passed to train_model: " << trainModel << endl;
						return false;
					}
				} else if ( flag == "output_cross_hyb" ) {
					if ( tokenizedInput.nextToken() == "true" ) {
						this->m_OutputCrossHyb = true;
					} else {
						this->m_OutputCrossHyb = false;
					}
				} else if ( flag == "cross_hyb_type" ) {
					m_CrossHybType = tokenizedInput.nextToken();
					if ( m_CrossHybType != "none" &&
						 m_CrossHybType != "filter") {
						cout << "ERROR:  Invalid argument passed to cross_hyb_type: " << m_CrossHybType << endl;
						return false;
					}
				} else if ( flag == "summarize_expression" ) {
					if ( tokenizedInput.nextToken() == "true" ) {
						this->m_SummarizeExpression = true;
					} else {
						this->m_SummarizeExpression = false;
					}
				} else if ( flag == "summary_method" ) {
					m_SummaryMethod = tokenizedInput.nextToken();
					if ( m_SummaryMethod != "selection" && m_SummaryMethod != "liwong" && 
						 m_SummaryMethod != "deconvolve" && m_SummaryMethod != "correlation_filter" ) {
						cout << "ERROR:  Invalid argument passed to summary_method: " << m_SummaryMethod << endl;
						return false;
					}
				} else if ( flag == "correlation_filter_min_correlation") {
					m_CorrelationFilterMinCorrelation = tokenizedInput.nextFloatToken();
				} else if ( flag == "correlation_filter_min_sd" ) { 
					m_CorrelationFilterMinSd = tokenizedInput.nextFloatToken();					
				} else if ( flag == "correltion_filter_output_masked" ) {
					if ( tokenizedInput.nextToken() == "true" ) {
						this->m_CorrelationFilterOutputMasked = true;
					} else {
						this->m_CorrelationFilterOutputMasked = false;
					}
				} else if ( flag == "correlation_filter_edit_distance" ) {
					m_CorrelationFilterEditDistance = tokenizedInput.nextIntToken();
				} else if ( flag == "mat_infile" ) {
					this->m_MATInfile = tokenizedInput.nextToken();
				} else if ( flag == "array_type" ) {
					this->m_ArrayType = tokenizedInput.nextToken();
					if ( m_ArrayType != "exon" 
						&& m_ArrayType != "tiling"
						&& m_ArrayType != "promoter" 
						&& m_ArrayType != "gene"
						&& m_ArrayType != "glue"
						&& m_ArrayType != "transcriptome" ) {
							cout << "ERROR:  Invalid argument passed to array_type: " << m_ArrayType << endl;
							return false;
						}
				} else if ( flag == "mat_training_probe_type" ) {
					this->m_MATTrainingProbeType = tokenizedInput.nextToken();
					if ( m_MATTrainingProbeType != "full" 
						&& m_MATTrainingProbeType != "extended"
						&& m_MATTrainingProbeType != "core"
						&& m_MATTrainingProbeType != "background" ) {
							cout << "ERROR:  Invalid argument passed to mat_training_probe_type: " << m_MATTrainingProbeType << endl;
							return false;
					}
				} else if ( flag == "include_copy_number_mat_parameter" ) {
					if ( tokenizedInput.nextToken() == "true" ) {
						this->m_IncludeCopyNumberMATParameter = true;
					} else {
						this->m_IncludeCopyNumberMATParameter = false;
					}
				} else if ( flag == "calculate_presence_absence" ) {
					if ( tokenizedInput.nextToken() == "true" ) {
						this->m_CalculatePresenceAbsence = true;
					} else {
						this->m_CalculatePresenceAbsence = false;
					}
				}
			}
		}
	}


	//
	// Now do some error checking on the input values
	//
	if ( Parameters::LOGFILE.is_open() == false ) {
		std::cout << "Cannot open logfile" << endl;
		return false;
	}
	double targetNumCellsExonArray = gsl_pow_2(m_CELLDIM_EXON_ARRAY);
	if ( targetNumCellsExonArray != m_NUMCELLS_EXON_ARRAY ) {
		std::cout << "ERROR:  Invalid arguments to num_cells and cell_dim!" << endl;
		return false;
	}
	string arrayType = this->m_ArrayType;
	if ( arrayType == "exon" ) {
		m_NUMCELLS = this->m_NUMCELLS_EXON_ARRAY;
		m_CELLDIM = this->m_CELLDIM_EXON_ARRAY;
	} else if ( arrayType == "gene" ) {
		m_NUMCELLS = this->m_NUMCELLS_GENE_ARRAY;
		m_CELLDIM = this->m_CELLDIM_GENE_ARRAY;
	} else if ( arrayType == "glue" ) {
		m_NUMCELLS = this->m_NUMCELLS_GLUE_ARRAY;
		m_CELLDIM = this->m_CELLDIM_GLUE_ARRAY;
	} else if ( arrayType == "tiling" ) {
		m_NUMCELLS = this->m_NUMCELLS_CHIP_CHIP_ARRAY;
		m_CELLDIM = this->m_CELLDIM_CHIP_CHIP_ARRAY;
	} else if ( arrayType == "promoter" ) {
		m_NUMCELLS = this->m_NUMCELLS_PROMOTER_ARRAY;
		m_CELLDIM = this->m_CELLDIM_PROMOTER_ARRAY;
	}
	Parameters::LOGFILE << "Success:  Parameters loaded.  " << endl
		<< "Array type: " << this->m_ArrayType << endl
		<< "Array dimensions:  " << m_CELLDIM << endl
		<< "Array total features:  " << m_NUMCELLS << endl;
	
	if ( this->m_TrainModel == true ) { 
		Parameters::LOGFILE << "Train background model: true" << endl;
		Parameters::LOGFILE << "Background model: " << this->m_MODEL_BACKGROUND_METHOD << endl;
		Parameters::LOGFILE << "Background model training probes: " << this->m_MATTrainingProbeType << endl;
		Parameters::LOGFILE << "Output background-corrected intensities: " << this->OutputAllBkgdCorrectNormProbes() << endl;		
	}
	
	if ( this->m_SummarizeExpression == true ) {
		Parameters::LOGFILE << "Summarize gene expression: true" << endl;
		Parameters::LOGFILE << "Summary method: " << this->m_SummaryMethod << endl;
		if ( this->m_SummaryMethod == "correlation_filter" ) {
			Parameters::LOGFILE << "Minimum probe correlation with off-target: " << this->m_CorrelationFilterMinCorrelation << endl;
			Parameters::LOGFILE << "Maximum edit distance between probe and off-target: " << this->m_CorrelationFilterEditDistance << endl;		
		}
	}
	if ( this->m_OutputCrossHyb == true ) {
		Parameters::LOGFILE << "Output probe cross-hybridization info: true" << endl;
	}

	return true;
}
pair<set<string>::iterator, set<string>::iterator> Parameters::GetCelFilesIterators()
{
	return make_pair( this->m_CelFiles.begin(), this->m_CelFiles.end() );
}


pair<set<string>::iterator, set<string>::iterator> Parameters::GetPromoterCelFilesIterators()
{
	return make_pair( this->m_PromoterCelFiles.begin(), this->m_PromoterCelFiles.end() );
}

pair<multimap<string, string>::iterator, multimap<string, string>::iterator> Parameters::GetChipChipCelFileIterators( const string& i_CharVal ) {
	return make_pair( this->m_ChipChipCelFiles.lower_bound( i_CharVal ), this->m_ChipChipCelFiles.upper_bound( i_CharVal ) );
}

pair< multimap<string, string>::iterator, multimap<string, string>::iterator> Parameters::GetChipChipCelFileIterators()
{
	return make_pair( this->m_ChipChipCelFiles.begin(), this->m_ChipChipCelFiles.end() );
}

pair< map<string, string>::iterator, map<string, string>::iterator > Parameters::GetChipChipBpmapIterators()
{
	return make_pair( this->m_ChipChipBpmap.begin(), this->m_ChipChipBpmap.end() );
}






bool Parameters::MOSTLY_POSITIVE( gsl_vector* i_pVec, int i_Size )
{
	int posCount = 0;
	int negCount = 0;
	for ( int i = 0; i < i_Size; ++i) 
	{
		double val = gsl_vector_get(i_pVec, i);
		if ( val > 0 ) {
			++posCount;
		} else if ( val < 0 ) {
			++negCount;
		}
	}
	
	if (negCount > posCount) {
		return false;
	} else {
		return true;
	}
}


int Parameters::RandomNumber( const int& i_Size )
{
	int rand = static_cast<int>( gsl_rng_uniform_int (this->m_pRng, i_Size) );
	return rand;
}	

gsl_rng* Parameters::GetRng()
{
	return this->m_pRng;
}

string Parameters::GetModelBackgroundMethod() const
{
	return this->m_MODEL_BACKGROUND_METHOD;
}

// NumCells and CellDim
int Parameters::GetNumCells() const
{
	return this->m_NUMCELLS;
}
int Parameters::GetCellDim() const
{
	return this->m_CELLDIM;
}


// Model Parameters
bool Parameters::TrainModel() const
{	
	return this->m_TrainModel;
}

bool Parameters::SummarizeExpression() const
{
	return this->m_SummarizeExpression;
}

int Parameters::GetNumProbesMAT() const
{
	return this->m_NUM_PROBES_MAT;
}

bool Parameters::GetChipChipBpmapFile( const string& i_Type, string& i_File ) const
{
	map< string, string >::const_iterator itFind = this->m_ChipChipBpmap.find( i_Type );
	if ( itFind == this->m_ChipChipBpmap.end() ) {
		i_File = "";
		return false;
	} else {
		i_File = itFind->second;
		return true;
	}
}

string Parameters::GetChipChipCelFolder() const
{
	return this->m_CHIP_CHIP_CEL_FOLDER;
}

bool Parameters::ChipChipTargetsReverse() const
{
	return this->m_ChipChipTargetsReverse;
}

bool Parameters::OutputAllBkgdCorrectNormProbes() const
{
	return this->m_OutputAllBkgdCorrectNormProbes;
}
bool Parameters::OutputCelFormat() const
{
	return this->m_OutputCelFormat;
}

string Parameters::GetBkgdCorrectNormProbesFile() const
{
	return this->m_BkgdCorrectNormProbesFile;
}

string Parameters::GetOutputModelFit() const
{
	return this->m_OutputModelFit;
}

string Parameters::GetMATTrainingProbeType() const
{
	return this->m_MATTrainingProbeType;
}

bool Parameters::CalculatePresenceAbsence() const
{
	return this->m_CalculatePresenceAbsence;
}

string Parameters::GetPresenceAbsenceOutfile() const
{
	return this->m_PresenceAbsenceOutfile;
}



string Parameters::GetArrayType() const 
{
	return this->m_ArrayType;
}
string Parameters::GetPromoterBpmapFile() const
{
	return this->m_PromoterBpmap;
}

string Parameters::GetPgfFile() const
{
	return this->m_PGFFILE;
}
string Parameters::GetClfFile() const
{
	return this->m_CLFFILE;
}

string Parameters::GetProbesetAnnotation() const
{
	return this->m_PROBESET_ANNOTATION;
}

string Parameters::GetArrayCelFolder() const
{
	return this->m_ARRAY_CEL_FOLDER;
}

string Parameters::GetChipChipBpmapFolder() const
{
	return this->m_CHIP_CHIP_BPMAP_FOLDER;
}
string Parameters::GetPromoterArrayCelFolder() const
{
	return this->m_PROMOTER_CEL_FOLDER;
}

string Parameters::GetCrossHybOutputFile() const
{	
	return this->m_CrossHybOutputFile;
}
string Parameters::GetCrossHybDir() const
{
	return this->m_CrossHybDir;
}
string Parameters::GetTcEstimExprFile() const
{
	return this->m_TcEstimExprFile;
}

int Parameters::GetNumArrays() const
{
	if ( m_ArrayType == "exon" || m_ArrayType == "gene" || m_ArrayType == "glue" ) {
		return this->m_CelFiles.size();
	} else if ( m_ArrayType == "tiling" ) {
		return this->m_ChipChipCelFiles.size();
	} else if ( m_ArrayType == "promoter" ) {
		return this->m_PromoterCelFiles.size();
	} else {
		return 0;
	}
}
string Parameters::GetNormalizationMethod() const
{
	return m_NormalizationMethod;
}
string Parameters::GetSummaryMethod() const
{
	return m_SummaryMethod;
}
double Parameters::GetCorrelationFilterMinCorrelation() const
{
	return this->m_CorrelationFilterMinCorrelation;
}
double Parameters::GetCorrelationFilterMinSd() const
{
	return this->m_CorrelationFilterMinSd;
}
double Parameters::GetCorrelationFilterOutputMasked() const
{
	return this->m_CorrelationFilterOutputMasked;
}
int Parameters::GetCorrelationFilterEditDistance() const
{
	return this->m_CorrelationFilterEditDistance;
}

string Parameters::REVERSE_COMPLEMENT( string i_Sequence )
{
	int len = i_Sequence.size();
	string revCompSeq = i_Sequence;
	for ( int i = 0; i < len; ++i ) {
		char currChar = i_Sequence[i];
		if ( currChar == 'A' || currChar == 'a' ) {
			revCompSeq[len - i - 1] = 'T';
		} else if ( currChar == 'C' || currChar == 'c' ) {
			revCompSeq[len - i - 1] = 'G';
		} else if ( currChar == 'G' || currChar == 'g' ) {
			revCompSeq[len - i - 1] = 'C';
		} else if ( currChar == 'T' || currChar == 't' ) {
			revCompSeq[len - i - 1] = 'A';
		} else {
			// error!  
			// do nothing 
		}
	}
	return revCompSeq;
}

string Parameters::SWAP_ORDER( const string& i_Sequence )
{
	string newSequence = i_Sequence;
	int len = i_Sequence.length();
	for ( int i = 0; i < len; ++i ) {
		newSequence[i] = i_Sequence[ len - i - 1 ];
	}
	return newSequence;
}

bool Parameters::IncludeCopyNumberMATParameter() const
{
	return this->m_IncludeCopyNumberMATParameter;
}
int Parameters::GetNumParameters_MAT() const
{
	return this->m_NumParameters_MAT;
}

string Parameters::GetWorkingDirectory() const
{
	return this->m_WorkingDirectory;
}
string Parameters::GetBkgdCoefFile() const
{
	return this->m_BkgdCoefFile;
}
bool Parameters::OutputCrossHyb() const
{
	return this->m_OutputCrossHyb;
}
string Parameters::GetCrossHybType() const
{
	return this->m_CrossHybType;
}

bool Parameters::GetTcEstimExpr( map< unsigned long, vector<double>* >& o_GeneExpression )
{
	Parameters* pParameters = Parameters::GetInstance();
	string geneExpressionsFilename = pParameters->GetTcEstimExprFile();
	ifstream geneExpressionsInfile;
	geneExpressionsInfile.open( geneExpressionsFilename.c_str() );
	if ( geneExpressionsInfile.is_open() == false ) {
		Parameters::LOGFILE << "ERROR:  Cannot open gene expresseions file " << geneExpressionsFilename << endl;
		return false;
	}

	string inputLine;
	while( !geneExpressionsInfile.eof() ) {
		getline( geneExpressionsInfile, inputLine );
		if ( inputLine.find("probe") != string::npos || inputLine.find("Transcript") != string::npos) {
			continue; // we got the header line
		}
		StringTokenizer tokenizedInputLine = StringTokenizer( inputLine, "\t" );
		unsigned long tcId = tokenizedInputLine.nextUnsignedLongIntToken();
		vector<double>* pVec = new vector<double>;
		while( tokenizedInputLine.hasMoreTokens() ) {
			double val = tokenizedInputLine.nextFloatToken();
			pVec->push_back( val );
		}
		o_GeneExpression.insert( make_pair( tcId, pVec ) );
	}
	return true;

}

	

