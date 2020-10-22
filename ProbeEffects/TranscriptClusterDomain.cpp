#include "stdafx.h"
#include "TranscriptClusterDomain.h"
#include "Parameters.h"
#include "StringTokenizer.h"
#include "TranscriptCluster.h"
#include "chip_data.h"
#include "ClfInfo.h"

bool TranscriptClusterDomain::m_InstanceFlag = false;
TranscriptClusterDomain* TranscriptClusterDomain::m_pTheOne = NULL;

TranscriptClusterDomain* TranscriptClusterDomain::GetInstance()
{
	if ( !m_InstanceFlag )
	{
		m_pTheOne = new TranscriptClusterDomain();
		m_InstanceFlag = true;
		return m_pTheOne;
	}
	else {
		return m_pTheOne;
	}
}

TranscriptClusterDomain::TranscriptClusterDomain(void)
{
	// Load transcript cluster info here
	bool success = this->LoadTranscriptClusters(); //this->LoadProbeInfo();
	if (success == false ) {
		exit(1);
	}
	// Other methods for initialization...  
}

bool TranscriptClusterDomain::LoadProbeInfo()
{	
	Parameters* pParameters = Parameters::GetInstance();
	string probeInfoFile = pParameters->GetWorkingDirectory() + "ProbeInfo.peaf";
	ifstream infile;
	infile.open( probeInfoFile.c_str() );
	bool writeProbeInfo = true;
	if ( infile.is_open() ) {
		// check for correct headers
		bool matchingHeader = false;
		string pgfCreateDateHeader;
		getline( infile, pgfCreateDateHeader );

		ifstream pgffile;
		pgffile.open( (pParameters->GetPgfFile()).c_str() ); 
		if (pgffile.is_open() == false) {
			Parameters::LOGFILE << "ERROR:  Cannot open pgf file.  " << endl;
			return false;
		}
		string pgfLine;
		while ( !pgffile.eof() ) {
			getline( pgffile, pgfLine );
			string::size_type pos = pgfLine.find ("#%create_date=",0);
			if ( pos != string::npos ) {
				break;
			}
		}
		if ( pgfLine == pgfCreateDateHeader ) {
			string probesetAnnotationCreateDateHeader;
			getline( infile, probesetAnnotationCreateDateHeader );

			ifstream probesetAnnotationFile;
			probesetAnnotationFile.open( pParameters->GetProbesetAnnotation().c_str() );
			if (probesetAnnotationFile.is_open() == false)
			{
				Parameters::LOGFILE << "ERROR:  Cannot open probeset annotation file.  "  << endl;
				return false;
			}
			string probesetAnnotationLine;
			while ( !probesetAnnotationFile.eof() ) {
				getline( probesetAnnotationFile, probesetAnnotationLine );
				string::size_type pos = probesetAnnotationLine.find( "#%create_date=",0);
				if ( pos != string::npos ) {
					break;
				}
			}
			if ( probesetAnnotationLine == probesetAnnotationCreateDateHeader ) {
				writeProbeInfo = false; // if headers are correct, set writeProbeInfo to false
			}
		}
	}

	if ( writeProbeInfo  ) {
		// Create the "ProbeInfo.peaf" file
		ofstream probeInfoOutfile;
		probeInfoOutfile.open( probeInfoFile.c_str() );
		if ( probeInfoOutfile.is_open() == false ) {
			Parameters::LOGFILE << "ERROR:  Cannot open " << probeInfoFile << "!  Is working directory set correctly in parameter file?  "  << endl;
			return false;
		}
		const string SEP = "\t";
		
		// Read the pgf file
		ifstream pgffile;
		string pgfCreateDateHeader = "";
		pgffile.open( (pParameters->GetPgfFile()).c_str() ); 
		if (pgffile.is_open()) {
			int probesetCounter = 0;
			int probeCounter = 0;
			unsigned long lastProbesetId = 0;
			string lastProbesetType = "";
			while (!pgffile.eof() ) {
				string inputLine;
				getline( pgffile, inputLine);
				if (inputLine[0] == '#'){
					string::size_type pos = inputLine.find ("#%create_date=",0);
					if ( pos != string::npos ) {
						pgfCreateDateHeader = inputLine;
					}
					continue; // Found a comment, to be ignored.
				}
				string::size_type pos1 = inputLine.find("\t", 0);
				if (pos1 != 0){
					// We have a probeset id line
					string probesetId = inputLine.substr(0, pos1);
					lastProbesetId = strtoul( probesetId.c_str(), NULL, 10 );
					string::size_type pos2 = inputLine.find_first_of( "\t\n\r", pos1 + 1 );
					string::size_type len = pos2 - pos1 - 1;
					if (pos2 == string::npos){
						len = inputLine.length() - pos1 - 1;
					}
					lastProbesetType = inputLine.substr(pos1 + 1, len);
					if ( lastProbesetType == "control->chip" ) {
						lastProbesetId = 0;
						lastProbesetType = "";
						continue;
					}
					probesetCounter++;		
				} else {
					string::size_type pos2 = inputLine.find("\t",pos1 + 1);
					if (pos2 == 1) {
						// We have a probe id line
						if (lastProbesetId == 0)
						{
							// Could be because we have a 'control->chip' probeset
							// Don't want to load the corresponding probes.  
							continue;
						}

						string::size_type pos3 = inputLine.find("\t", pos2 + 1);
						string::size_type len = pos3 - pos2 - 1;
						string probeId = inputLine.substr(pos2 + 1, len);
						unsigned long probeIdVal = strtoul ( probeId.c_str(), NULL, 10 ); 
						pos2 = pos3;
						pos3 = inputLine.find("\t", pos2 + 1);
						len = pos3 - pos2 - 1;
						string probeType = inputLine.substr(pos2 + 1, len);
						pos2 = pos3;
						pos3 = inputLine.find("\t", pos2 + 1);
						len = pos3 - pos2 - 1;
						string gcCount = inputLine.substr(pos2 + 1, len);
						pos2 = pos3;
						pos3 = inputLine.find("\t", pos2 + 1);
						len = pos3 - pos2 - 1;
						string probeLength = inputLine.substr(pos2 + 1, len);
						pos2 = pos3;
						pos3 = inputLine.find("\t", pos2 + 1);
						len = pos3 - pos2 - 1;
						string interrogationPosition = inputLine.substr(pos2 + 1, len);
						pos2 = pos3;
						pos3 = inputLine.find_first_of("\t\n\r", pos2 + 1);
						len = pos3 - pos2 - 1;
						if (pos3 == string::npos){
							len = inputLine.length() - pos2 - 1;
						}
						string sequence = inputLine.substr(pos2 + 1, len);

						string level = "none";
						if (lastProbesetType == "control->bgp->genomic" ) {
							level = "genomic_background";
						} else if (lastProbesetType == "control->bgp->antigenomic" ) {
							level = "antigenomic_background";
						}

						ProbeInfo* pProbeInfo = new ProbeInfo( probeIdVal, 0, level, sequence );
						this->m_ProbeInfo.insert( make_pair(lastProbesetId, pProbeInfo) );
					}	
				}
			}
		} else {
			Parameters::LOGFILE << "ERROR" << "\t" << "Cannot open pgf file " << pParameters->GetPgfFile() << " in TranscriptClusterDomain::ReadProbesets.  " << endl;
			return false; // cannot open the pgf file
		}
		Parameters::LOGFILE << "SUCCESS" << "\t" << "Exon array pgf file read successfully." << endl;
		
		// Read the probeset annotation file
		ifstream probesetAnnotationFile;
		string probesetAnnotationCreateDateHeader = "";
		probesetAnnotationFile.open( pParameters->GetProbesetAnnotation().c_str() );
		if (probesetAnnotationFile.is_open() == false){
			Parameters::LOGFILE << "ERROR:  Cannot open probeset annotation file" << endl;
			return false;
		}
		int numTokens; 
		while ( !probesetAnnotationFile.eof() )
		{
			string inputLine;
			getline( probesetAnnotationFile, inputLine);
			
			//Parse each line
			if ( inputLine[0] == '#'){
				string::size_type pos = inputLine.find ("#%create_date=",0);
				if ( pos != string::npos ) {
					probesetAnnotationCreateDateHeader = inputLine;
				}
				// Found a comment.  To be ignored.  
				continue;
			}
			StringTokenizer tokenizedInputLine( inputLine, "," );
			if ( inputLine.find("probeset_id") != string::npos) {
				numTokens = tokenizedInputLine.countTokens(); 
				// Ignore the header line
				continue;
			}
			if ( tokenizedInputLine.countTokens() != numTokens ) {
				continue;
			}
			
			// Remove the leading " and trailing " from each string.  
			// Note that if the file was read using microsoft excel, then 
			// the "'s will already have been removed!  
			string probesetId = tokenizedInputLine.nextToken();
			probesetId = probesetId.substr( 1, (probesetId.length() - 2) );
			string chromosome = tokenizedInputLine.nextToken();
			chromosome = chromosome.substr( 1, (chromosome.length() - 2) );
			string strand = tokenizedInputLine.nextToken();
			strand = strand.substr( 1, (strand.length() - 2) );
			string start = tokenizedInputLine.nextToken();
			start = start.substr( 1, (start.length() - 2) );
			string stop = tokenizedInputLine.nextToken();
			stop = stop.substr( 1, (stop.length() - 2) );
			string probeCountStr = tokenizedInputLine.nextToken();
			int probeCount = atoi( probeCountStr.substr( 1, (probeCountStr.length() - 2) ).c_str() );
			string transcriptClusterId = tokenizedInputLine.nextToken();
			transcriptClusterId = transcriptClusterId.substr( 1, (transcriptClusterId.length() - 2) );
			unsigned long transcriptClusterIdVal = strtoul( transcriptClusterId.c_str(), NULL, 10);
			string exonId = tokenizedInputLine.nextToken();
			exonId = exonId.substr( 1, (exonId.length() - 2) );
			string psrId = tokenizedInputLine.nextToken();
			psrId = psrId.substr( 1, (psrId.length() - 2) );
			string geneAssignment = tokenizedInputLine.nextToken();
			geneAssignment = geneAssignment.substr( 1, (geneAssignment.length() - 2) );
			string mrnaAssignment = tokenizedInputLine.nextToken();
			mrnaAssignment = mrnaAssignment.substr( 1, (mrnaAssignment.length() - 2) );
			string ProbesetTypeStr = tokenizedInputLine.nextToken();
			int ProbesetType = atoi( ProbesetTypeStr.substr( 1, (ProbesetTypeStr.length() - 2) ).c_str() );
			string numberIndependentProbesStr = tokenizedInputLine.nextToken();
			int numberIndependentProbes = atoi( numberIndependentProbesStr.substr( 1, (numberIndependentProbesStr.length() - 2) ).c_str() );
			string numberCrossHybProbesStr = tokenizedInputLine.nextToken();
			int numberCrossHybProbes = atoi( numberCrossHybProbesStr.substr( 1, (numberCrossHybProbesStr.length() - 2) ).c_str() );
			string numberNonoverlappingProbesStr = tokenizedInputLine.nextToken();
			int numberNonoverlappingProbes = atoi( numberNonoverlappingProbesStr.substr( 1, (numberNonoverlappingProbesStr.length() - 2) ).c_str() );
			string level = tokenizedInputLine.nextToken(); // This we care about
			level = level.substr( 1, (level.length() - 2) );

			string boundedStr = tokenizedInputLine.nextToken();
			int bounded = atoi( boundedStr.substr( 1, (boundedStr.length() - 2) ).c_str() );
			string noBoundedEvidenceStr = tokenizedInputLine.nextToken();
			int noBoundedEvidence = atoi( noBoundedEvidenceStr.substr( 1, (noBoundedEvidenceStr.length() - 2) ).c_str() );
			string hasCdsStr = tokenizedInputLine.nextToken();
			int hasCds = atoi( hasCdsStr.substr( 1, (hasCdsStr.length() - 2) ).c_str() );
			string flStr  = tokenizedInputLine.nextToken();
			int fl = atoi( flStr.substr( 1, (flStr.length() - 2) ).c_str() );
			string mrnaStr = tokenizedInputLine.nextToken();
			int mrna = atoi( mrnaStr.substr( 1, (mrnaStr.length() - 2) ).c_str() );
			string estStr = tokenizedInputLine.nextToken();
			int est = atoi( estStr.substr( 1, (estStr.length() - 2) ).c_str() );
			string vegaGeneStr = tokenizedInputLine.nextToken();
			int vegaGene = atoi( vegaGeneStr.substr( 1, (vegaGeneStr.length() - 2) ).c_str() );
			string vegaPseudoGeneStr = tokenizedInputLine.nextToken();
			int vegaPseudoGene = atoi( vegaPseudoGeneStr.substr( 1, (vegaPseudoGeneStr.length() - 2) ).c_str() );
			string ensGeneStr = tokenizedInputLine.nextToken();
			int ensGene = atoi( ensGeneStr.substr( 1, (ensGeneStr.length() - 2) ).c_str() );
			string sgpGeneStr = tokenizedInputLine.nextToken();
			int sgpGene = atoi( sgpGeneStr.substr( 1, (sgpGeneStr.length() - 2) ).c_str() );
			string exoniphyStr = tokenizedInputLine.nextToken();
			int exoniphy = atoi( exoniphyStr.substr( 1, (exoniphyStr.length() - 2) ).c_str() );
			string twinscanStr = tokenizedInputLine.nextToken();
			int twinscan = atoi( twinscanStr.substr( 1, (twinscanStr.length() - 2) ).c_str() );
			string geneidStr = tokenizedInputLine.nextToken();
			int geneid = atoi( geneidStr.substr( 1, (geneidStr.length() - 2) ).c_str() );
			string genscanStr = tokenizedInputLine.nextToken();
			int genscan = atoi( genscanStr.substr( 1, (genscanStr.length() - 2) ).c_str() );
			string genscanSuboptStr = tokenizedInputLine.nextToken();
			int genscanSubopt = atoi( genscanSuboptStr.substr( 1, (genscanSuboptStr.length() - 2) ).c_str() );
			string mouseFlStr = tokenizedInputLine.nextToken();
			int mouseFl = atoi( mouseFlStr.substr( 1, (mouseFlStr.length() - 2) ).c_str() );
			string mouseMrnaStr = tokenizedInputLine.nextToken();
			int mouseMrna = atoi( mouseMrnaStr.substr( 1, (mouseMrnaStr.length() - 2) ).c_str() );
			string ratFlStr = tokenizedInputLine.nextToken();
			int ratFl = atoi( ratFlStr.substr( 1, (ratFlStr.length() - 2) ).c_str() );
			string ratMrnaStr = tokenizedInputLine.nextToken();
			int ratMrna = atoi( ratMrnaStr.substr( 1, (ratMrnaStr.length() - 2) ).c_str() );
			string microRnaRegistryStr = tokenizedInputLine.nextToken();
			int microRnaRegistry = atoi( microRnaRegistryStr.substr( 1, (microRnaRegistryStr.length() - 2) ).c_str() );
			string rnaGeneStr = tokenizedInputLine.nextToken();
			int rnaGene = atoi( rnaGeneStr.substr( 1, (rnaGeneStr.length() - 2) ).c_str() );
			string mitomapStr = tokenizedInputLine.nextToken();
			int mitomap = atoi( mitomapStr.substr( 1, (mitomapStr.length() - 2) ).c_str() );

			unsigned long probesetIdVal = strtoul( probesetId.c_str(), NULL, 10);
			pair<hash_multimap< unsigned long, ProbeInfo* >::iterator, hash_multimap< unsigned long, ProbeInfo* >::iterator > itFindPair = this->m_ProbeInfo.equal_range(probesetIdVal);
			hash_multimap< unsigned long, ProbeInfo* >::iterator itFindBegin = itFindPair.first;
			hash_multimap< unsigned long, ProbeInfo* >::iterator itFindEnd = itFindPair.second;
			for (; itFindBegin != itFindEnd; ++itFindBegin) {
				ProbeInfo* pProbeInfo = itFindBegin->second;
				pProbeInfo->Annotate( transcriptClusterIdVal, level );
			}
		}

		// Output the info
		probeInfoOutfile << pgfCreateDateHeader << endl;
		probeInfoOutfile << probesetAnnotationCreateDateHeader << endl;
		hash_multimap< unsigned long, ProbeInfo* >::iterator itBegin = this->m_ProbeInfo.begin();
		hash_multimap< unsigned long, ProbeInfo* >::iterator itEnd = this->m_ProbeInfo.end();
		for (; itBegin != itEnd; ++itBegin) {
			ProbeInfo* pProbeInfo = itBegin->second;
			probeInfoOutfile << itBegin->first << SEP 
				<< pProbeInfo->m_ProbeId << SEP
				<< pProbeInfo->m_TranscriptClusterId << SEP
				<< pProbeInfo->m_Level << SEP
				<< pProbeInfo->m_Sequence << endl;
			delete pProbeInfo;
			pProbeInfo = NULL;
		}
		probeInfoOutfile.close();
		m_ProbeInfo.clear();		

		// Finally, open the "ProbeInfo.peaf" file
		infile.open( probeInfoFile.c_str() );
	}

	// Process the "ProbeInfo.peaf" file

	// Read and store clf info.  
	bool success = this->ReadClf();
	if ( success == false ) {
		return false;
	}

	// For each line in the file
	unsigned long prevProbesetId = 0;
	Probeset* pPrevProbeset = NULL;
	while ( !infile.eof() ) {
		string inputLine;
		getline( infile, inputLine );
		if ( inputLine[0] == '#' ) {
			// Found a comment
			continue;
		}

		// Create transcript clusters, probesets, probes
		StringTokenizer tokenizedInputLine( inputLine, "\t" );
		unsigned long probesetId = tokenizedInputLine.nextUnsignedLongIntToken();
		unsigned long probeId = tokenizedInputLine.nextUnsignedLongIntToken();
		unsigned long tcId = tokenizedInputLine.nextUnsignedLongIntToken();
		string level = tokenizedInputLine.nextToken();
		string sequence = tokenizedInputLine.nextToken();
		sequence = Parameters::SWAP_ORDER( sequence ); // Exon array probes have sequences 3' to 5' in the pgf
														// We swap the order to follow the 5' to 3' convention

		Probe* pProbe = new Probe( probeId, sequence );
		if ( probeId >= this->m_ClfInfo.size() ) {
			Parameters::LOGFILE << "ERROR:  Probe id not found from Clf file!" << endl;
			delete pProbe;
			continue;
		}
		map<unsigned long, ClfInfo >::iterator itClfInfo = this->m_ClfInfo.find( probeId );
		if ( itClfInfo == m_ClfInfo.end() ) {
			Parameters::LOGFILE << "ERROR:  Mapping error in ClfInfo.  " 
				<< "Probe " << probeId << " has wrong clf annotation." << endl;
		} else {
			int x = itClfInfo->second.GetX();
			int y = itClfInfo->second.GetY();
			bool isPosSet = pProbe->SetPosition(x,y);
			if (isPosSet == false) {
				Parameters::LOGFILE << "Probe " << probeId << " has negative x or y coordinates: " << x << ", " << y << endl;
			}
			bool isCrossHyb = itClfInfo->second.GetIsCrossHyb();
			if ( isCrossHyb == true ) {
				string crossHybInfo = itClfInfo->second.GetChInfo();
				StringTokenizer tokCrossHybInfo = StringTokenizer( crossHybInfo, ";" );
				while ( tokCrossHybInfo.hasMoreTokens() ) {
					string currCrossHybInfo = tokCrossHybInfo.nextToken();
					StringTokenizer tokCurrCrossHybInfo = StringTokenizer( currCrossHybInfo, "," );
					string tcId = tokCurrCrossHybInfo.nextToken();
					int mm = tokCurrCrossHybInfo.nextIntToken();
					int ind = tokCurrCrossHybInfo.nextIntToken();
					pProbe->AddCrossHybInfo( tcId, mm, ind );
				}
			}
		}

		if ( tcId > 0 ) { // Background probes have tcId == 0
			TranscriptCluster* pTranscriptCluster;
			hash_map<unsigned long, TranscriptCluster*>::iterator itFindTc = m_TranscriptClusters.find( tcId );
			if ( itFindTc == m_TranscriptClusters.end() ) {
				pTranscriptCluster = new TranscriptCluster( tcId );
				m_TranscriptClusters.insert( make_pair(tcId, pTranscriptCluster) );
			} else {
				pTranscriptCluster = itFindTc->second;
			}
			Probeset* pProbeset;
			if ( probesetId != prevProbesetId ) {
				// add pProbeset to pTranscriptCluster
				pProbeset = new Probeset(probesetId, "" );
				AnnotationConf annot;
				if ( level == "core" ) {
					annot = CORE;
				} else if ( level == "extended" ) {
					annot = EXTENDED;
				} else if (level == "full" ) {
					annot = FULL;
				} else if ( level == "ambiguous" ) {
					annot = AMBIGUOUS;
				} else if ( level == "free" ) {
					annot = FREE;
				} else {
					annot = NONE;
				}
				pTranscriptCluster->AddProbeset( pProbeset, annot );
				prevProbesetId = probesetId;
				pPrevProbeset = pProbeset;
			} else {
				pProbeset = pPrevProbeset;
			}
			pProbeset->AddProbe( pProbe );
		} else if ( level == "genomic_background" ) {
			int gcCount = pProbe->GetGcCount();
			this->m_GenomicBgdProbes.insert(make_pair(gcCount, pProbe));
		} else if ( level == "antigenomic_background" ) {
			int gcCount = pProbe->GetGcCount();
			this->m_AntigenomicBgdProbes.insert(make_pair(gcCount, pProbe));
		}
	}

	// Clear the clf file
	this->m_ClfInfo.clear();

	Parameters::LOGFILE << "SUCCESS:  Probes, probesets, and transcript clusters created successfully.  " << endl;	
	return true;
}


TranscriptClusterDomain::~TranscriptClusterDomain(void)
{
	// 
	// TBD:  add code here
	//

	hash_map<unsigned long, Probeset*>::iterator itPBegin = this->m_Probesets.begin();
	hash_map<unsigned long, Probeset*>::iterator itPEnd = this->m_Probesets.end();
	for (; itPBegin != itPEnd; ++itPBegin) {
		Probeset* pProbeset = itPBegin->second;
		delete pProbeset;
		pProbeset = NULL;
	}

	hash_map<unsigned long, TranscriptCluster*>::iterator itTCBegin = this->m_TranscriptClusters.begin();
	hash_map<unsigned long, TranscriptCluster*>::iterator itTCEnd = this->m_TranscriptClusters.end();
	for (; itTCBegin != itTCEnd; ++itTCBegin) {
		TranscriptCluster* pTranscriptCluster = itTCBegin->second;
		delete pTranscriptCluster;
		pTranscriptCluster = NULL;
	}



}

bool TranscriptClusterDomain::LoadTranscriptClusters() {
	// First read and store clf info.  
	bool success = this->ReadClf();
	if ( success == false ) {
		return false;
	}
	
	// Second read the probesets from the pgf file.  
	success = this->ReadProbesets();
	if (success == false ) {
		return false;
	}

	Parameters* pParameters = Parameters::GetInstance();
	string arrayType = pParameters->GetArrayType();
	if ( arrayType == "exon" || arrayType == "glue" ) {
		// Next annotate the probesets.  
		// Here TranscriptClusters are created and
		// low conf probesets are stored separately.  
		success = this->AnnotateProbesetsExonArray();
	} else if ( arrayType == "gene" ) {
		success = this->AnnotateProbesetsGeneArray();
	}
	if ( success == false ) {
		return false;
	}



	return true;
}

//
// Reads the probesets and their corresponding probes from the 
// exon array pgf file.  
//
bool TranscriptClusterDomain::ReadProbesets()
{
	Parameters* pParameters = Parameters::GetInstance();
	bool setExcludeProbe = false;
	if ( pParameters->GetCrossHybType() == "filter" ) {
		setExcludeProbe = true;
	}
	ifstream infile;
	infile.open( (pParameters->GetPgfFile()).c_str() ); 
	if (infile.is_open()) {
		// logfile << "Pgf file opened" << endl;
		int probesetCounter = 0;
		Probeset* pLastProbeset = NULL;
		while (!infile.eof() ) {
			string inputLine;
			getline( infile, inputLine);
			if (inputLine[0] == '#')
			{
				// Found a comment, to be ignored.
				continue;
			}

			string::size_type pos1 = inputLine.find("\t", 0);
			if (pos1 != 0)
			{
				// We have a probeset id line
				string probesetId = inputLine.substr(0, pos1);
				string::size_type pos2 = inputLine.find_first_of( "\t\n\r", pos1 + 1 );
				string::size_type len = pos2 - pos1 - 1;
				if (pos2 == string::npos)
				{
					len = inputLine.length() - pos1 - 1;
				}
				string probesetType = inputLine.substr(pos1 + 1, len);
				if ( (probesetType == "control->chip") || (probesetType == "normgene->intron") 
					|| (probesetType == "normgene->exon") || (probesetType == "control->affx") ) {
					pLastProbeset = NULL; // Ignore these probesets
					continue;
				}

				//i_LogFile << "Probeset found with id " << probesetId << " and type " << probesetType << endl;
				unsigned long probesetIdVal = strtoul( probesetId.c_str(), NULL, 10 );
				pLastProbeset = new Probeset(probesetIdVal, probesetType);
				this->m_Probesets.insert(make_pair(probesetIdVal, pLastProbeset));
				probesetCounter++;
							
			} else {
				string::size_type pos2 = inputLine.find("\t",pos1 + 1);
				if (pos2 == 1)
				{
					// We have a probe id line
					string::size_type pos3 = inputLine.find("\t", pos2 + 1);
					string::size_type len = pos3 - pos2 - 1;
					string probeId = inputLine.substr(pos2 + 1, len);
					unsigned long probeIdVal = strtoul( probeId.c_str(), NULL, 10 );
					pos2 = pos3;
					pos3 = inputLine.find("\t", pos2 + 1);
					len = pos3 - pos2 - 1;
					string probeType = inputLine.substr(pos2 + 1, len);
					pos2 = pos3;
					pos3 = inputLine.find("\t", pos2 + 1);
					len = pos3 - pos2 - 1;
					string gcCount = inputLine.substr(pos2 + 1, len); // Convert to an integer?  
					pos2 = pos3;
					pos3 = inputLine.find("\t", pos2 + 1);
					len = pos3 - pos2 - 1;
					string probeLength = inputLine.substr(pos2 + 1, len); // Convert to an integer?
					pos2 = pos3;
					pos3 = inputLine.find("\t", pos2 + 1);
					len = pos3 - pos2 - 1;
					string interrogationPosition = inputLine.substr(pos2 + 1, len); // Convert to an integer?
					pos2 = pos3;
					pos3 = inputLine.find_first_of("\t\n\r", pos2 + 1);
					len = pos3 - pos2 - 1;
					if (pos3 == string::npos)
					{
						len = inputLine.length() - pos2 - 1;
					}
					string sequence = inputLine.substr(pos2 + 1, len);

					//i_LogFile << "Probe id " << probeId << endl;
					if (pLastProbeset == NULL)
					{
						// Could be because we have a 'control->chip', 'control->affx', 
						// 'normgene->intron', or 'normgene->exon' probeset
						// Don't want to load the corresponding probes.  
						continue;
					}
					
					bool isUniqueProbe = pLastProbeset->AddProbe(probeIdVal, sequence);
					if (isUniqueProbe == false)
					{

						//i_LogFile << "Error!  Probe " << probeId << " already exists!  " << endl;
					}
					Probe* pProbe = pLastProbeset->GetProbe(probeId);
					if ( probeIdVal >= this->m_ClfInfo.size() ) {
						Parameters::LOGFILE << "ERROR:  Probe id not found from Clf file!" << endl;
						continue;
					}

					map<unsigned long, ClfInfo >::iterator itClfInfo = this->m_ClfInfo.find(probeIdVal ); // 04/01/08:  changed to a map to handle non-square cel files
					if ( itClfInfo == m_ClfInfo.end() ) {						
						Parameters::LOGFILE << "ERROR:  Mapping error in ClfInfo.  " 
							<< "Probe " << probeId << " has wrong clf annotation." << endl;
					} else {
						int x = itClfInfo->second.GetX();
						int y = itClfInfo->second.GetY();
						bool isPosSet = pProbe->SetPosition(x,y);
						if (isPosSet == false) {
							Parameters::LOGFILE << "Probe " << probeId << " has negative x or y coordinates: " << x << ", " << y << endl;
						}
						bool isCrossHyb = itClfInfo->second.GetIsCrossHyb();
						if ( isCrossHyb == true ) {
							string crossHybInfo = itClfInfo->second.GetChInfo();
							StringTokenizer tokCrossHybInfo = StringTokenizer( crossHybInfo, ";" );
							while ( tokCrossHybInfo.hasMoreTokens() ) {
								string currCrossHybInfo = tokCrossHybInfo.nextToken();
								StringTokenizer tokCurrCrossHybInfo = StringTokenizer( currCrossHybInfo, "," );
								string tcId = tokCurrCrossHybInfo.nextToken();
								int mm = tokCurrCrossHybInfo.nextIntToken();
								int ind = tokCurrCrossHybInfo.nextIntToken();
								pProbe->AddCrossHybInfo( tcId, mm, ind );
							}
							if ( setExcludeProbe == true ) {
								pProbe->SetExcludeProbe(true); // isCrossHyb == true => excludeProbe = true
							}
						}
					}

					// Add pProbe to m_GenomicBgdProbes or m_AntigenomicBgdProbes
					// Note:  probes with other probe types are ignored
					if ( pLastProbeset->IsGenomicBackgroundProbeset() == true )
					{
						this->m_GenomicBgdProbes.insert(make_pair(atoi(gcCount.c_str()), pProbe));
					} else if ( pLastProbeset->IsAntigenomicBackgroundProbeset() == true )
					{
						this->m_AntigenomicBgdProbes.insert(make_pair(atoi(gcCount.c_str()), pProbe));
					}
				}	
			}
		}
	} else {
		Parameters::LOGFILE << "ERROR" << "\t" << "Cannot open pgf file " << pParameters->GetPgfFile() << " in TranscriptClusterDomain::ReadProbesets.  " << endl;
		return false; // cannot open the pgf file
	}
	Parameters::LOGFILE << "SUCCESS" << "\t" << "Exon array pgf file read successfully." << endl;
	this->m_ClfInfo.clear();
	return true;
}


bool TranscriptClusterDomain::AnnotateProbesetsGeneArray()
{
	hash_map<unsigned long, Probeset*>::iterator itPsBegin = this->m_Probesets.begin();
	hash_map<unsigned long, Probeset*>::iterator itPsEnd = this->m_Probesets.end();
	for (; itPsBegin != itPsEnd; ++itPsBegin) {
		Probeset* pProbeset = itPsBegin->second;
		if ( pProbeset->IsAntigenomicBackgroundProbeset() || pProbeset->IsGenomicBackgroundProbeset() ) {
			continue; // Don't create a transcript cluster from this probeset
		} else {
			unsigned long tcId = pProbeset->GetProbesetId();
			TranscriptCluster* pTranscriptCluster = new TranscriptCluster( tcId );
			this->m_TranscriptClusters.insert( make_pair( tcId, pTranscriptCluster ) );
			AnnotationConf annot = CORE; // all probesets assumed to be of annotation core
			pTranscriptCluster->AddProbeset( pProbeset, annot );
		}
	}
	return true;
}

// 
// This method will assign probesets to their corresponding transcript clusters.  
// Also, probesets with low annotational confidence will be added to a map
// of low conf probesets.  
//
bool TranscriptClusterDomain::AnnotateProbesetsExonArray() {
	
	Parameters* pParameters = Parameters::GetInstance();
	ifstream infile;
	infile.open( pParameters->GetProbesetAnnotation().c_str() );
	int numTokens; // Read in from the header file.  Used to verify that an input line is complete.  
	if (infile.is_open() == true)
	{
		//i_LogFile << "Probeset annotation file opened" << endl;
		int numAnnotations = 0;
		int numCoreAnnotations = 0;
		int numExtendedAnnotations = 0;
		int numFullAnnotations = 0;
		int numAmbiguousAnnotations = 0;
		int numFreeAnnotations = 0;
		int numNoAnnotations = 0;

		while ( !infile.eof() )
		{
			string inputLine;
			getline( infile, inputLine);
			
			//Parse each line
			if ( inputLine[0] == '#')
			{
				// Found a comment.  To be ignored.  
				continue;
			}

			StringTokenizer tokenizedInputLine( inputLine, "," );

			if ( inputLine.find("probeset_id") != string::npos) {
				numTokens = tokenizedInputLine.countTokens(); 
				// Ignore the header line
				continue;
			}
			
			if ( tokenizedInputLine.countTokens() != numTokens ) {
				continue;
			}
			
			// Remove the leading " and trailing " from each string.  
			// Note that if the file was read using microsoft excel, then 
			// the "'s will already have been removed!  
			string probesetId = tokenizedInputLine.nextToken();
			probesetId = probesetId.substr( 1, (probesetId.length() - 2) );
			string chromosome = tokenizedInputLine.nextToken();
			chromosome = chromosome.substr( 1, (chromosome.length() - 2) );
			string strand = tokenizedInputLine.nextToken();
			strand = strand.substr( 1, (strand.length() - 2) );
			string start = tokenizedInputLine.nextToken();
			start = start.substr( 1, (start.length() - 2) );
			string stop = tokenizedInputLine.nextToken();
			stop = stop.substr( 1, (stop.length() - 2) );
			string probeCountStr = tokenizedInputLine.nextToken();
			int probeCount = atoi( probeCountStr.substr( 1, (probeCountStr.length() - 2) ).c_str() );
			string transcriptClusterId = tokenizedInputLine.nextToken();
			transcriptClusterId = transcriptClusterId.substr( 1, (transcriptClusterId.length() - 2) );
			string exonId = tokenizedInputLine.nextToken();
			exonId = exonId.substr( 1, (exonId.length() - 2) );
			string psrId = tokenizedInputLine.nextToken();
			psrId = psrId.substr( 1, (psrId.length() - 2) );
			string geneAssignment = tokenizedInputLine.nextToken();
			geneAssignment = geneAssignment.substr( 1, (geneAssignment.length() - 2) );
			string mrnaAssignment = tokenizedInputLine.nextToken();
			mrnaAssignment = mrnaAssignment.substr( 1, (mrnaAssignment.length() - 2) );
			string ProbesetTypeStr = tokenizedInputLine.nextToken();
			int ProbesetType = atoi( ProbesetTypeStr.substr( 1, (ProbesetTypeStr.length() - 2) ).c_str() );
			string numberIndependentProbesStr = tokenizedInputLine.nextToken();
			int numberIndependentProbes = atoi( numberIndependentProbesStr.substr( 1, (numberIndependentProbesStr.length() - 2) ).c_str() );
			string numberCrossHybProbesStr = tokenizedInputLine.nextToken();
			int numberCrossHybProbes = atoi( numberCrossHybProbesStr.substr( 1, (numberCrossHybProbesStr.length() - 2) ).c_str() );
			string numberNonoverlappingProbesStr = tokenizedInputLine.nextToken();
			int numberNonoverlappingProbes = atoi( numberNonoverlappingProbesStr.substr( 1, (numberNonoverlappingProbesStr.length() - 2) ).c_str() );
			string level = tokenizedInputLine.nextToken(); // This we care about
			level = level.substr( 1, (level.length() - 2) );

			//string boundedStr = tokenizedInputLine.nextToken();
			//int bounded = atoi( boundedStr.substr( 1, (boundedStr.length() - 2) ).c_str() );
			//string noBoundedEvidenceStr = tokenizedInputLine.nextToken();
			//int noBoundedEvidence = atoi( noBoundedEvidenceStr.substr( 1, (noBoundedEvidenceStr.length() - 2) ).c_str() );
			//string hasCdsStr = tokenizedInputLine.nextToken();
			//int hasCds = atoi( hasCdsStr.substr( 1, (hasCdsStr.length() - 2) ).c_str() );
			//string flStr  = tokenizedInputLine.nextToken();
			//int fl = atoi( flStr.substr( 1, (flStr.length() - 2) ).c_str() );
			//string mrnaStr = tokenizedInputLine.nextToken();
			//int mrna = atoi( mrnaStr.substr( 1, (mrnaStr.length() - 2) ).c_str() );
			//string estStr = tokenizedInputLine.nextToken();
			//int est = atoi( estStr.substr( 1, (estStr.length() - 2) ).c_str() );
			//string vegaGeneStr = tokenizedInputLine.nextToken();
			//int vegaGene = atoi( vegaGeneStr.substr( 1, (vegaGeneStr.length() - 2) ).c_str() );
			//string vegaPseudoGeneStr = tokenizedInputLine.nextToken();
			//int vegaPseudoGene = atoi( vegaPseudoGeneStr.substr( 1, (vegaPseudoGeneStr.length() - 2) ).c_str() );
			//string ensGeneStr = tokenizedInputLine.nextToken();
			//int ensGene = atoi( ensGeneStr.substr( 1, (ensGeneStr.length() - 2) ).c_str() );
			//string sgpGeneStr = tokenizedInputLine.nextToken();
			//int sgpGene = atoi( sgpGeneStr.substr( 1, (sgpGeneStr.length() - 2) ).c_str() );
			//string exoniphyStr = tokenizedInputLine.nextToken();
			//int exoniphy = atoi( exoniphyStr.substr( 1, (exoniphyStr.length() - 2) ).c_str() );
			//string twinscanStr = tokenizedInputLine.nextToken();
			//int twinscan = atoi( twinscanStr.substr( 1, (twinscanStr.length() - 2) ).c_str() );
			//string geneidStr = tokenizedInputLine.nextToken();
			//int geneid = atoi( geneidStr.substr( 1, (geneidStr.length() - 2) ).c_str() );
			//string genscanStr = tokenizedInputLine.nextToken();
			//int genscan = atoi( genscanStr.substr( 1, (genscanStr.length() - 2) ).c_str() );
			//string genscanSuboptStr = tokenizedInputLine.nextToken();
			//int genscanSubopt = atoi( genscanSuboptStr.substr( 1, (genscanSuboptStr.length() - 2) ).c_str() );
			//string mouseFlStr = tokenizedInputLine.nextToken();
			//int mouseFl = atoi( mouseFlStr.substr( 1, (mouseFlStr.length() - 2) ).c_str() );
			//string mouseMrnaStr = tokenizedInputLine.nextToken();
			//int mouseMrna = atoi( mouseMrnaStr.substr( 1, (mouseMrnaStr.length() - 2) ).c_str() );
			//string ratFlStr = tokenizedInputLine.nextToken();
			//int ratFl = atoi( ratFlStr.substr( 1, (ratFlStr.length() - 2) ).c_str() );
			//string ratMrnaStr = tokenizedInputLine.nextToken();
			//int ratMrna = atoi( ratMrnaStr.substr( 1, (ratMrnaStr.length() - 2) ).c_str() );
			//string microRnaRegistryStr = tokenizedInputLine.nextToken();
			//int microRnaRegistry = atoi( microRnaRegistryStr.substr( 1, (microRnaRegistryStr.length() - 2) ).c_str() );
			//string rnaGeneStr = tokenizedInputLine.nextToken();
			//int rnaGene = atoi( rnaGeneStr.substr( 1, (rnaGeneStr.length() - 2) ).c_str() );
			//string mitomapStr = tokenizedInputLine.nextToken();
			//int mitomap = atoi( mitomapStr.substr( 1, (mitomapStr.length() - 2) ).c_str() );

			// Ensure that we only include probesets with the appropriate level of annotation confidence
			AnnotationConf annotConf;
			if ( level == "core" ) {
				annotConf = CORE;
				numCoreAnnotations++;
			} else if ( level == "extended" ) {
				annotConf = EXTENDED;
				numExtendedAnnotations++;
			} else if ( level == "full" ) {
				annotConf = FULL;
				numFullAnnotations++;
			} else if ( level == "free" ) {
				annotConf = FREE;
				numFreeAnnotations++;
			} else if ( level == "ambiguous" ) {
				annotConf = AMBIGUOUS;
				numAmbiguousAnnotations++;
			} else {
				annotConf = NONE;
				numNoAnnotations++;
			}	

			unsigned long probesetIdVal = strtoul( probesetId.c_str(), NULL, 10);
			hash_map<unsigned long, Probeset*>::iterator it = this->m_Probesets.find(probesetIdVal);
			if (it == this->m_Probesets.end() ) { 
				//i_LogFile << "Error:  Probeset " << probesetId 
					//<< " listed in annotation file, but not found in pgf file! " << endl;
			} else {
				Probeset* pProbeset = it->second;
				pProbeset->Annotate(geneAssignment, start, stop, level, transcriptClusterId); // Can add more annotation later...
				
				if ( annotConf >= FULL ) {
					std::vector<int> gcVec;
					pProbeset->GetGCVec( gcVec );
					this->m_LowConfProbesets.insert( make_pair(gcVec, pProbeset) ); // For modeling unexpressed probeset expression
				}
							
				TranscriptCluster* pTranscriptCluster;
				unsigned long transcriptClusterIdVal = strtoul( transcriptClusterId.c_str(), NULL, 10);
				hash_map<unsigned long, TranscriptCluster*>::iterator it = this->m_TranscriptClusters.find( transcriptClusterIdVal );
				if ( it != this->m_TranscriptClusters.end() ) {
					pTranscriptCluster = it->second;
				} else {
					pTranscriptCluster = new TranscriptCluster( transcriptClusterIdVal );
					this->m_TranscriptClusters.insert( make_pair( transcriptClusterIdVal, pTranscriptCluster ) );
				}
				pTranscriptCluster->AddProbeset( pProbeset, annotConf );

				// 
				// We can filter this->m_Probesets to exclude those probesets with low annotational confidence if necessary.  
				// Note that no probesets will be deleted (or their corresp probes).  
				// Don't think this will be necessary for our purposes...  
				//
				//if ( annotConf <= i_AnnotationConfidence ) { 
				//	numAnnotations++;
				//} else {  // Remove the probeset from m_ProbesetsAnnotConf if the annotation confidence is insufficient
					// Discard the probeset from this->m_Probesets
				//	this->m_Probesets.erase( probesetId );
				//}
				//
			}


		}
		infile.close();
	} else {
		Parameters::LOGFILE << "ERROR" << "\t" << "Cannot open probeset annotation file in TranscriptClusterDomain::AnnotateProbesets.  " << endl;
		return false;
	}
	Parameters::LOGFILE << "SUCCESS" << "\t" << "Success in TranscriptClusterDomain::AnnotateProbesets!" << endl;
	return true;
}

bool TranscriptClusterDomain::ReadClf() {
	ifstream infile;
	infile.open( Parameters::GetInstance()->GetClfFile().c_str());
	if (infile.is_open() == true)
	{
		while (!infile.eof() )
		{
			string inputLine;
			getline( infile, inputLine);
			if (inputLine[0] == '#')
			{
				// Found a comment, to be ignored.  
				continue;
			}
			// Remove the last \r if it exists
			string::size_type pos = inputLine.find_last_of("\r");
			if ( pos != string::npos ) {
				inputLine = inputLine.substr(0, (inputLine.length() - 1) );
			}	
			StringTokenizer tokInputLine = StringTokenizer( inputLine, "\t" );
			unsigned long probeId = tokInputLine.nextUnsignedLongIntToken();
			int x = tokInputLine.nextIntToken();
			int y = tokInputLine.nextIntToken();
			bool isCrossHyb = false;
			string crossHybInfo = "";
			if ( tokInputLine.hasMoreTokens() ) 
			{
				crossHybInfo = tokInputLine.nextToken();
				isCrossHyb = true;
			}
			ClfInfo clfInfo = ClfInfo( x, y, isCrossHyb, crossHybInfo );
			m_ClfInfo.insert( make_pair(probeId, clfInfo) );
		}
	} else {
		Parameters::LOGFILE << "ERROR:  " << "Cannot open exon array clf file!" << endl;
		return false;
	}
	Parameters::LOGFILE << "SUCCESS:  " << "Exon array clf file read successfully"  << endl;
	return true;
}

bool TranscriptClusterDomain::AnnotateCrossHybExonArray() 
{
	// TBD!!!
	// Maybe we should create a probeset annotation that contains this info
	// Save having to go through all tcs one more time
	// Or add the annotation to the clf file
	// That would be better.  
	
	return true;
}

pair<hash_map<unsigned long, TranscriptCluster*>::iterator, hash_map<unsigned long, TranscriptCluster*>::iterator> TranscriptClusterDomain::GetTranscriptClusterIterators()
{
	return make_pair( this->m_TranscriptClusters.begin(), this->m_TranscriptClusters.end() );
}

TranscriptCluster* TranscriptClusterDomain::GetTranscriptCluster( const string& i_TranscriptClusterId )
{
	unsigned long transcriptClusterIdVal = strtoul( i_TranscriptClusterId.c_str(), NULL, 10 );
	TranscriptCluster* pTc = this->GetTranscriptCluster(transcriptClusterIdVal);
	return pTc;
}

TranscriptCluster* TranscriptClusterDomain::GetTranscriptCluster( const unsigned int i_TranscriptClusterId )
{
	hash_map< unsigned long, TranscriptCluster*>::iterator it = this->m_TranscriptClusters.find( i_TranscriptClusterId );
	if ( it != this->m_TranscriptClusters.end() ) {
		return it->second;
	} else {
		return NULL;
	}
}


Probeset* TranscriptClusterDomain::GetProbeset( const string& i_ProbesetId )
{
	unsigned long probesetIdVal = strtoul( i_ProbesetId.c_str(), NULL, 10 );
	hash_map< unsigned long, Probeset*>::iterator it = this->m_Probesets.find( probesetIdVal );
	if ( it != this->m_Probesets.end() ) {
		return it->second;
	} else {
		return NULL;
	}
}

bool TranscriptClusterDomain::IsTranscriptCluster( const string& i_TranscriptClusterId ) const
{
	unsigned long transcriptClusterIdVal = strtoul( i_TranscriptClusterId.c_str(), NULL, 10 );
	hash_map< unsigned long, TranscriptCluster*>::const_iterator it = this->m_TranscriptClusters.find( transcriptClusterIdVal );
	if ( it != this->m_TranscriptClusters.end() ) {
		return true;
	} else {
		return false;
	}
}

Probeset* TranscriptClusterDomain::GetRandomLowConfProbeset( gsl_vector_int* i_pVec )
{
	int numProbes = static_cast<int>(i_pVec->size);
	vector<int> gcContent;
	// A problem if we have fewer than 4 probes!!!  
	for ( int i =0; i < numProbes; ++i ) {
		gcContent.push_back( gsl_vector_int_get(i_pVec, i ));
	}

	int size = static_cast<int>(this->m_LowConfProbesets.count( gcContent ));
	if ( size == 0 ) {
		return NULL;
	}
	int rand = Parameters::GetInstance()->RandomNumber( size );
	multimap< vector<int>, Probeset*>::iterator it = this->m_LowConfProbesets.lower_bound( gcContent );
	// it += rand;
	for (int i =0; i < rand; ++i ) {
		it++;
	}
	return it->second;
}

int TranscriptClusterDomain::GetNumGenomicBgdProbes( int i_Key ) const
{
	int len = static_cast<int>(this->m_GenomicBgdProbes.count( i_Key ));
	return len;
}
int TranscriptClusterDomain::GetNumAntigenomicBgdProbes( int i_Key ) const
{
	return static_cast<int>( this->m_AntigenomicBgdProbes.count( i_Key ) );
}
pair <multimap<int, Probe*>::iterator, multimap<int, Probe*>::iterator> TranscriptClusterDomain::GetGenomicBgdIterators( int i_Key )
{
	return make_pair( this->m_GenomicBgdProbes.lower_bound(i_Key), this->m_GenomicBgdProbes.upper_bound(i_Key) );
}
pair <multimap<int, Probe*>::iterator, multimap<int, Probe*>::iterator> TranscriptClusterDomain::GetAntigenomicBgdIterators( int i_Key )
{
	return make_pair( this->m_AntigenomicBgdProbes.lower_bound(i_Key), this->m_AntigenomicBgdProbes.upper_bound(i_Key) );
}
pair< multimap<int, Probe*>::iterator, multimap<int, Probe*>::iterator > TranscriptClusterDomain::GetGenomicBgdIterators()
{
	return make_pair( this->m_GenomicBgdProbes.begin(), this->m_GenomicBgdProbes.end() );
}
pair< multimap<int, Probe*>::iterator, multimap<int, Probe*>::iterator > TranscriptClusterDomain::GetAntigenomicBgdIterators()
{
	return make_pair( this->m_AntigenomicBgdProbes.begin(), this->m_AntigenomicBgdProbes.end() );
}
int TranscriptClusterDomain::GetNumGenomicBgdProbes() const
{
	return static_cast<int>( this->m_GenomicBgdProbes.size() );
}
int TranscriptClusterDomain::GetNumAntigenomicBgdProbes() const
{
	return static_cast<int>( this->m_AntigenomicBgdProbes.size() );
}

bool TranscriptClusterDomain::OutputCoreProbes( const string& i_File )
{
	ofstream outfile;
	outfile.open( i_File.c_str() );
	if ( outfile.is_open() == false ) {
		return false;
	}
	hash_map<unsigned long, TranscriptCluster*>::iterator itTcBegin = m_TranscriptClusters.begin();
	hash_map<unsigned long, TranscriptCluster*>::iterator itTcEnd = m_TranscriptClusters.end();
	for (; itTcBegin != itTcEnd; ++itTcBegin) {
		TranscriptCluster* pTranscriptCluster = itTcBegin->second;
		unsigned long transcriptClusterId = pTranscriptCluster->GetId();
		pair<multimap<AnnotationConf, Probeset*>::iterator, multimap<AnnotationConf, Probeset*>::iterator > itPsPair = pTranscriptCluster-> GetProbesetIterators( CORE );
		multimap<AnnotationConf, Probeset*>::iterator itPsBegin = itPsPair.first;
		multimap<AnnotationConf, Probeset*>::iterator itPsEnd = itPsPair.second;
		for (; itPsBegin != itPsEnd; ++itPsBegin) {
			Probeset* pProbeset = itPsBegin->second;
			unsigned long probesetId = pProbeset->GetProbesetId();
			pair<map<unsigned long, Probe*>::iterator, map<unsigned long, Probe*>::iterator> itPrPair = pProbeset->GetProbeIterators();
			map< unsigned long, Probe*>::iterator itPrBegin = itPrPair.first;
			map< unsigned long, Probe*>::iterator itPrEnd = itPrPair.second;
			for (; itPrBegin != itPrEnd; ++itPrBegin) {
				Probe* pProbe = itPrBegin->second;
				pair<int, int> posPair = pProbe->GetPosition();
				outfile << ">" << posPair.first << "_" << posPair.second << "_" << pProbe->GetProbeId() << "_" 
					<< probesetId << "_" << transcriptClusterId << endl;
				outfile << pProbe->GetSequence() << endl;
			}
		}
	}


	return true;
}

bool TranscriptClusterDomain::OutputExtendedProbes( const string& i_File )
{
	ofstream outfile;
	outfile.open( i_File.c_str() );
	if ( outfile.is_open() == false ) {
		return false;
	}

	hash_map<unsigned long, TranscriptCluster*>::iterator itTcBegin = m_TranscriptClusters.begin();
	hash_map<unsigned long, TranscriptCluster*>::iterator itTcEnd = m_TranscriptClusters.end();
	for (; itTcBegin != itTcEnd; ++itTcBegin) {
		TranscriptCluster* pTranscriptCluster = itTcBegin->second;
		unsigned long transcriptClusterId = pTranscriptCluster->GetId();
		pair<multimap<AnnotationConf, Probeset*>::iterator, multimap<AnnotationConf, Probeset*>::iterator > itPsPair = pTranscriptCluster-> GetProbesetIterators( EXTENDED );
		multimap<AnnotationConf, Probeset*>::iterator itPsBegin = itPsPair.first;
		multimap<AnnotationConf, Probeset*>::iterator itPsEnd = itPsPair.second;
		for (; itPsBegin != itPsEnd; ++itPsBegin) {
			Probeset* pProbeset = itPsBegin->second;
			unsigned long probesetId = pProbeset->GetProbesetId();
			pair<map<unsigned long, Probe*>::iterator, map<unsigned long, Probe*>::iterator> itPrPair = pProbeset->GetProbeIterators();
			map< unsigned long, Probe*>::iterator itPrBegin = itPrPair.first;
			map< unsigned long, Probe*>::iterator itPrEnd = itPrPair.second;
			for (; itPrBegin != itPrEnd; ++itPrBegin) {
				Probe* pProbe = itPrBegin->second;
				pair<int, int> posPair = pProbe->GetPosition();
				outfile << ">" << posPair.first << "_" << posPair.second << "_" << pProbe->GetProbeId() << "_" 
					<< probesetId << "_" << transcriptClusterId << endl;
				outfile << pProbe->GetSequence() << endl;
			}
		}
	}

	return true;
}

bool TranscriptClusterDomain::OutputFullProbes( const string& i_File )
{
	ofstream outfile;
	outfile.open( i_File.c_str() );
	if ( outfile.is_open() == false ) {
		return false;
	}
	hash_map<unsigned long, TranscriptCluster*>::iterator itTcBegin = m_TranscriptClusters.begin();
	hash_map<unsigned long, TranscriptCluster*>::iterator itTcEnd = m_TranscriptClusters.end();
	for (; itTcBegin != itTcEnd; ++itTcBegin) {
		TranscriptCluster* pTranscriptCluster = itTcBegin->second;
		unsigned long transcriptClusterId = pTranscriptCluster->GetId();
		pair<multimap<AnnotationConf, Probeset*>::iterator, multimap<AnnotationConf, Probeset*>::iterator > itPsPair = pTranscriptCluster-> GetProbesetIterators( FULL );
		multimap<AnnotationConf, Probeset*>::iterator itPsBegin = itPsPair.first;
		multimap<AnnotationConf, Probeset*>::iterator itPsEnd = itPsPair.second;
		for (; itPsBegin != itPsEnd; ++itPsBegin) {
			Probeset* pProbeset = itPsBegin->second;
			unsigned long probesetId = pProbeset->GetProbesetId();
			pair<map<unsigned long, Probe*>::iterator, map<unsigned long, Probe*>::iterator> itPrPair = pProbeset->GetProbeIterators();
			map< unsigned long, Probe*>::iterator itPrBegin = itPrPair.first;
			map< unsigned long, Probe*>::iterator itPrEnd = itPrPair.second;
			for (; itPrBegin != itPrEnd; ++itPrBegin) {
				Probe* pProbe = itPrBegin->second;
				pair<int, int> posPair = pProbe->GetPosition();
				outfile << ">" << posPair.first << "_" << posPair.second << "_" << pProbe->GetProbeId() << "_" 
					<< probesetId << "_" << transcriptClusterId << endl;
				outfile << pProbe->GetSequence() << endl;
			}
		}
	}

	return true;
}


bool TranscriptClusterDomain::OutputBackgroundProbes( const string& i_File )
{
	ofstream outfile;
	outfile.open( i_File.c_str() );
	if ( outfile.is_open() == false ) {
		return false;
	}
	pair< multimap<int, Probe*>::iterator, multimap<int, Probe*>::iterator > itAntiGenPair = this-> GetAntigenomicBgdIterators();
	multimap<int, Probe*>::iterator itAntiGenBegin = itAntiGenPair.first;
	multimap<int, Probe*>::iterator itAntiGenEnd = itAntiGenPair.second;
	for (; itAntiGenBegin != itAntiGenEnd; ++itAntiGenBegin) {
		Probe* pProbe = itAntiGenBegin->second;
		pair<int, int> posPair = pProbe->GetPosition();
		outfile << ">" << posPair.first << "_" << posPair.second << "_" << pProbe->GetProbeId() << endl;
		outfile << pProbe->GetSequence() << endl;
	}

	pair< multimap<int, Probe*>::iterator, multimap<int, Probe*>::iterator > itGenPair = this-> GetGenomicBgdIterators();
	multimap<int, Probe*>::iterator itGenBegin = itGenPair.first;
	multimap<int, Probe*>::iterator itGenEnd = itGenPair.second;
	for (; itGenBegin != itGenEnd; ++itGenBegin) {
		Probe* pProbe = itGenBegin->second;
		pair<int, int> posPair = pProbe->GetPosition();
		outfile << ">" << posPair.first << "_" << posPair.second << "_" << pProbe->GetProbeId() << endl;
		outfile << pProbe->GetSequence() << endl;
	}

	return true;
}


bool TranscriptClusterDomain::OutputRandomProbeIntensities( const string& i_ProbeFile, const string& i_Outfile )
{
	Parameters* pParameters = Parameters::GetInstance();
	map< pair<int, int>, vector<double>* > randProbes;

	ifstream probefile;
	probefile.open(i_ProbeFile.c_str() );
	if ( probefile.is_open() == false ) {		
		Parameters::LOGFILE << "Error:  Cannot open random probe file.  " << endl;
		return false;
	}
	string inputline;
	getline( probefile, inputline ); // capture the header line
	while ( !probefile.eof() ) 
	{
		getline( probefile, inputline );
		StringTokenizer tokline = StringTokenizer( inputline, "\t" );
		string probeId = tokline.nextToken();
		unsigned long probesetId = tokline.nextUnsignedLongIntToken();
		hash_map<unsigned long, Probeset*>::iterator itFind = this->m_Probesets.find( probesetId );
		if ( itFind == this->m_Probesets.end() ) {
			Parameters::LOGFILE << "Missing probeset: "<< probesetId << endl;
			continue; 
		}
		Probeset* pProbeset	= itFind->second;
		Probe* pProbe = pProbeset->GetProbe(probeId);
		if ( pProbe == NULL ) {
			Parameters::LOGFILE << "Missing probe: " << probeId << endl;
			continue;
		}
		pair<int, int> posPair = pProbe->GetPosition();
		vector<double>* pVec = new vector<double>;
		randProbes.insert( make_pair(posPair, pVec) );
	}
	
	// tell logfile how many probes were found
	Parameters::LOGFILE << "Success:  " << randProbes.size() << " random probes loaded. " << endl;
	
	// for each cel file
	// get the probe intensity
	pair< set<string>::iterator, set<string>::iterator > itChipPair = pParameters->GetCelFilesIterators();
	set<string>::iterator itChipBegin = itChipPair.first;
	set<string>::iterator itChipEnd = itChipPair.second;
	for (; itChipBegin != itChipEnd; ++itChipBegin) {
		// Read in the chip's data
		string filename = pParameters->GetArrayCelFolder() + (*itChipBegin);
		ifstream datafile;
		datafile.open(filename.c_str());
		
		if (datafile.is_open())
		{
			string::size_type pos = filename.rfind( Parameters::SEPARATOR );
			int len = static_cast<int>(filename.length());
			string smallname = filename.substr( (pos + 1), (len - 1) );
			CHIP_DATA* pChip = new CHIP_DATA( smallname, pParameters->GetNumCells(), pParameters->GetCellDim()  );				
			bool readOk = pChip->ReadBinaryCel(filename.c_str());
			if (readOk == false) {
				delete pChip;
				continue;
			}
			
			map< pair<int, int>, vector<double>* >::iterator itB = randProbes.begin();
			map< pair<int, int>, vector<double>* >::iterator itE = randProbes.end();
			for (; itB != itE; ++itB ) {
				pair<int, int> loc = itB->first;
				double inten = pChip->GetUnnormalizedIntensity(loc.first,loc.second);
				vector<double>* pVec = itB->second;
				pVec->push_back(inten);
			}
		} else {
			Parameters::LOGFILE << "ERROR:  " << "Cannot open CEL file " << filename << endl;
		}
	}
	
	// Now output the values to the outfile
	ofstream outfile;
	outfile.open(i_Outfile.c_str() );
	map< pair<int, int>, vector<double>* >::iterator itB = randProbes.begin();
	map< pair<int, int>, vector<double>* >::iterator itE = randProbes.end();
	for (; itB != itE; ++itB ) {
		vector<double>* pVec = itB->second;
		vector<double>::iterator itValBegin = pVec->begin();
		vector<double>::iterator itValEnd = pVec->end();
		for (; itValBegin != itValEnd; ++itValBegin) {
			outfile << *itValBegin << ",";
		}
		outfile << endl;
		delete pVec;
		pVec = NULL;
	}
	outfile.close();
	
	return true;
}






















