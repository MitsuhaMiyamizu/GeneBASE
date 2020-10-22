// ProbeEffects.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "math.h"
#include "TranscriptClusterDomain.h"
#include "cluster.h"
#include "AnalysisDomain.h"
int main(int argc, char* argv[])
{		
	//********************************
	// Handle command line arguments
	//********************************
	if ( argc == 1) { // Only the default argument is passed
		std::cout << "Options: "  << "\n" 
			<< "-par file specifying parameter values" << "\n" 
			<< endl;
		return 1;
	}
	std::string parameterInfile = "";
	for ( int i = 1; i < argc; ++i ) {
		std::string inputFlag = argv[i];
		++i;
		if ( inputFlag == "-par" ) {
			parameterInfile = argv[i];
		} else {
			std::cout << "Options: "  << "\n" 
				<< "-par file specifying parameter values" << "\n" 
				<< endl;
			return 1;
		}
	}
	

	//********************************
	// Set up the parameters
	//********************************
	Parameters* pParameters = Parameters::GetInstance();
	bool success = pParameters->LoadParameters( parameterInfile );
	if ( success == false ) {
		std::cout << "ERROR:  Cannot load parameters from specified parameter file.  "  << endl;
		return 1;
	}

	string modelMethod = pParameters->GetModelBackgroundMethod(); // Either "MAT" or "ProbeClassify"	
	string arrayType = pParameters->GetArrayType(); // Either exon, tiling, promoter, gene or glue
	AnalysisDomain* pAnalysisDomain = AnalysisDomain::GetInstance();
	//ChipChipDomain* pChipChipDomain = ChipChipDomain::GetInstance( numArrays ); 

	
	
	bool summarizeExpression = pParameters->SummarizeExpression(); // If true, get gene-level estimates
	if ( summarizeExpression == true ) {
		success = pAnalysisDomain->SummarizeExpression();
		//pChipChipDomain->SummarizeExpression(pParameters->GetSummaryMethod(), arrayType );
		if ( success == false ) {
			Parameters::LOGFILE << "ERROR:  Summarize expression failed.  " << endl;
			return 1;
		} else {
			Parameters::LOGFILE << "SUCCESS:  Summarize expression succeeded.  " << endl;
		}
	}
	
	bool trainModel = pParameters->TrainModel(); // If true, we train a model
	if ( trainModel == true ) {
		//*******************************
		// Train the background model using:
		// -- MAT
		// -- MedianGC
		// -- etc.  
		// Output parameter values.  
		//*******************************
		success = pAnalysisDomain->ModelBackground();
		if ( success == false ) {
			Parameters::LOGFILE << "ERROR:  Model background failed.  " << endl;
			return 1;
		} else {
			Parameters::LOGFILE << "SUCCESS:  Model background succeeded.  " << endl;
		}
	} 
		
	bool outputCrossHyb = pParameters->OutputCrossHyb();
	if ( outputCrossHyb == true ) {
		//CrossHybDomain* pCrossHybDomain = CrossHybDomain::GetInstance();
		bool success = pAnalysisDomain->OutputCrossHybInfo();
		//pCrossHybDomain->OutputCrossHybInfo();
		if ( success == false ) {
			Parameters::LOGFILE << "Error:  Output cross hyb info failed.  " << endl;
		} else {
			Parameters::LOGFILE << "Success:  Output cross hyb info succeeded.  " << endl;
		}
	}
	
	bool calculatePresenceAbsence = pParameters->CalculatePresenceAbsence();
	if ( calculatePresenceAbsence == true ) {
		bool success = pAnalysisDomain->ModelPresenceAbsence();
		if ( success == false ) {
			Parameters::LOGFILE << "Error:  Model presence/absence failed.  " << endl;
		} else {
			Parameters::LOGFILE << "Success:  Model presence/absence succeeded.  " << endl;
		}
	}

	return 0;
}










