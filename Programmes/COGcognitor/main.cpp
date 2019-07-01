#include "cognitor.h"
#include "os.h"
#include "cognitorglob.h"
#include "enum.h"

using namespace std;

int main(int argc, char* argv[]) // -p -d -o -a -c -n -t -e -m -g -w -l -r
{
	double	cogOverlapThreshold = 0.25;
	double	overlapThreshold	= 0.5;
	int		numThreshold		= 3; // maximum number of genome-specific hits in each cog; when numThreshold = 1 'classic' BeT algorithm applied
	int		totalNumThreshold	= 15; // maximum total number of hits in each cog (was 15, 100 is for "infinity")
	double	expectThreshold		= 1;
	int		orgThreshold		= 1; // minimum number of genomes into which a given query has hits in given COG; orgThreshold is 2 for 'classic' BeT alghrithm
	int		gapThreshold		= 60; // usually 60
	double	weightExp			= 4;  
	double	lengthExp			= 4; 
	double	rankExp				= 2; 
	double	rankOrg				= 0;
	bool quiet = true;

    string dirpath = ".";
	dirpath = dirpath + SPRTR + "conv";
	string data = "domdat.csv";
	string fout ="out.csv";
	string q2o = "prot2org.csv";
//	string cogs = "cogs.csv";

	string inFile;

	char	c;
	string	str;
	double	cT;
	double	oT;
	int		nT;
	int		tT;
	double	eT;
	int		orT;
	int		gapT;
	double	wE;
	double	lE;
	double	rE;
	double  zE;

	while(--argc > 0 && (*++argv)[0] == '-')
	{
		c = *++argv[0];		
		switch (c)
		{
			case 'i': // input directory (processed BLAST) + hash.csv
				str = argv[0];
				dirpath = str.substr(str.find('=') + 1);
				break;
			case 'q': // queries2orgs.csv
				str = argv[0];
				q2o = str.substr(str.find('=') + 1);
				break;
			case 't': // dom.dat
				str = argv[0];
				data = str.substr(str.find('=') + 1);
				break;
			case 'o': // output file
				str = argv[0];
				fout = str.substr(str.find('=') + 1);
				break;
//			case 'x': // cogs.csv
//				str = argv[0];
//				cogs = str.substr(str.find('=') + 1);
//				break;
			case 'e':
				str = argv[0];
				eT = atof(str.substr(str.find('=') + 1).c_str());
				expectThreshold = eT > 0 ? eT : expectThreshold;
				break;
			case 'a':
				str = argv[0];
				oT = atof(str.substr(str.find('=') + 1).c_str());
				overlapThreshold = oT > 0 && oT <= 1 ? oT : overlapThreshold;
				break;
			case 'b':
				str = argv[0];
				cT = atof(str.substr(str.find('=') + 1).c_str());
				cogOverlapThreshold = cT > 0 && cT <= 1 ? cT : cogOverlapThreshold;
				break;
			case 'g':
				str = argv[0];
				nT = atoi(str.substr(str.find('=') + 1).c_str());
				numThreshold = nT > 0 ? nT : numThreshold;
				break;
			case 'c':
				str = argv[0];
				tT = atoi(str.substr(str.find('=') + 1).c_str());
				totalNumThreshold = tT > 0 ? tT : totalNumThreshold;
				break;
			case 'm':
				str = argv[0];
				orT = atoi(str.substr(str.find('=') + 1).c_str());
				orgThreshold = orT > 0 ? orT : orgThreshold;
				break;
			case 'f':
				str = argv[0];
				gapT = atoi(str.substr(str.find('=') + 1).c_str());
				gapThreshold = gapT > 0 ? gapT : gapThreshold;
				break;
			case 'w':
				str = argv[0];
				wE = atof(str.substr(str.find('=') + 1).c_str());
				weightExp = wE > 0 ? wE : weightExp;
				break;
			case 'l':	
				str = argv[0];
				lE = atof(str.substr(str.find('=') + 1).c_str());
				lengthExp = lE > 0 ? lE : lengthExp;
				break;
			case 'r':
				str = argv[0];
				rE = atof(str.substr(str.find('=') + 1).c_str());
				rankExp = rE > 0 ? rE : rankExp;
				break;
			case 'z':
				str = argv[0];
				zE = atof(str.substr(str.find('=') + 1).c_str());
				rankOrg = zE > 0 ? zE : rankOrg;
			case 'v':
				quiet = false;
				break;
			default:
				cerr << endl;
				cerr << "Usage: COGcognitor [Options]" << endl;
				cerr << "Options:" << endl << endl;
				cerr << "-i=dconv	input directory (with hash.csv, hits.csv etc.; default ./conv)" << endl << endl;
				cerr << "-t=ftarget	file with target COG data (default domdat.csv" << endl << endl;
				cerr << "-q=fquery	file with query data (default prot2org.csv)" << endl << endl;
				cerr << "-o=fout	output file (default out.csv)" << endl << endl;
//				cerr << "-x=cogdef	COG definition file (default cogs.csv)" << endl << endl;				
//				cerr << "-a=overlap	hit coverage threshold (default 0.5)" << endl << endl; //overlapThreshold
//				cerr << "-b=overlap	COG-block overlap treshold (default 0.25)" << endl  << endl; //cogOverlapThreshold
				cerr << "-a=x		block assemble overlap (default 0.5)" << endl << endl; //overlapThreshold
				cerr << "-b=x		block conflict overlap (default 0.25)" << endl  << endl; //cogOverlapThreshold
				cerr << "-g=n		maximum number of hits into tha same genome taken for a given query (default 3)" << endl << endl; //numThreshold; if 1, classic BeT
				cerr << "-c=n		maximum number of COG-specific hits taken for a give query (default 15)" << endl << endl; //totalNumThreshold 				
				cerr << "-e=c		e-value threshold for an HSP (default 1.0)" << endl << endl; //expectThreshold				
				cerr << "-m=n		minimum number of different genomes for COG assignment (default 1)" << endl << endl; //orgThreshold; if 3, classic triangle		
				cerr << "-f=n		minimum length of free domain (default 60)" << endl << endl; //gapThreshold			
				cerr << "-w=x		exponent coefficient for weight of hit (default 4.0)" << endl << endl; //weightExp
				cerr << "-l=x		exponent coefficient for length of hit (default 4.0)" << endl  << endl; //lengthExp
				cerr << "-r=x		exponent coefficient for rank of hit (default 2.0)" << endl << endl; //rankExp
				cerr << "-z=x		exponent coefficient for rank of genome (default 0.0)" << endl << endl; //rankOrg
				cerr << "-v		verbose mode (mostly debugging output to STDOUT)" << endl << endl;	
				exit(1);
		}
	}

	string curPathL = "cognitor.log";
	ofstream log(curPathL.c_str());
	if (!log)
	{
		cerr << "Error opening output file: cognitor.log" << endl;
		exit(1);
	}

	ofstream outFileConfl("conflict.txt");
	if (!outFileConfl)
	{
		cerr << "Error opening output file: conflict.txt" << endl;
		exit(1);
	}

	if (quiet == true)
	{
		streambuf *psbuf;
		psbuf = log.rdbuf();
		cout.rdbuf(psbuf);
	}

	string curRec; // contains the current processing record
	string buffer;

	map<string, int> domHash;
	unsigned int indx;
	string subj;
	int numID= 0;

	map<string, int> orgHash;
	map<string, int> cogHash;
	map<int, string> cogHashR;

	Enumerator enCog(1);
	Enumerator enOrg(1);

	map<int, int>			dom2org;
	map<int, int>			dom2cog;
	map<int, int>			cog2dom;
	map<string, int>		prot2org;
	map<string, long>		offsets; // Presuming all hits for any query are successive
	map<string, string>		selfHits;
	map<int, string>		cogDescript;
	map<int, int>			dom2len;
	map<int, int>			cogMedian;
	map<int, int>			genIdNew2taxidOld; // temporary


	ofstream fOut(fout.c_str());
    if (!fOut)
    {
    	cerr << "Error opening output file: fout" << endl;
    	exit(1);
    }

	inFile = dirpath + SPRTR;

	ifstream inDomHash((inFile + "hash.csv").c_str());
	if (!inDomHash)
	{
		cout << "File " << inFile + "hash.csv not found" << endl;
		cout << "Native target IDs must be numeric" << endl;
		cout << "Trying to proceed..." << endl << endl;
	}
	else
	{
		while(getline(inDomHash, str))
		{
			indx = str.find_last_of(',');
			if (indx != string::npos)
			{
				subj = str.substr(indx + 1);
				numID = atoi((str.substr(0, indx)).c_str());
				domHash.insert(pair<string, int>(subj, numID));
			}
		}
	}

	ifstream inQ2S((inFile + "query2subject.csv").c_str());
	if (!inQ2S)
	{
		cout << "File " << inFile + "query2subect.csv" << " not found" << endl;
	}
	else
	{
		makeOffsets(&offsets, inQ2S);
		if (offsets.size() == 0)
		{
			cout << "File query2subject.csv is empty" << endl;
		}
	}

	ifstream inHitFile((inFile + "hits.csv").c_str());
	if (!inHitFile)
	{
		cerr << "Error opening input file: " << inFile + "hits.csv" << endl;
		exit(1);
	}

	ifstream inSelf((inFile + "self.csv").c_str());
	if (!inSelf)
	{
		cout << "File " << inFile + "self.csv" << " not found; BeT algorithm will be used" << endl;
		cout << "All querys' length set to 100" << endl;
		numThreshold = 1; // maximum number of hits in each genome; when numThreshold = 1 'classic' BeT algorithm applied
		orgThreshold = 2; //
	}
	else
	{
		populate(&selfHits, inSelf);
		if (selfHits.size() == 0)
		{
			cout << "File self.csv empty; BeT algorithm will be used" << endl;
			cout << "All querys' length set to 100" << endl;		
			numThreshold = 1; // maximum number of hits in each genome; when numThreshold = 1 'classic' BeT algorithm applied
			orgThreshold = 2; 
		}
	}

	ifstream inDat(data.c_str());
	if (!inDat)
	{
		cout << "Error opening input file: " << data << endl;
		exit(1);
	}

	ifstream inP2O(q2o.c_str()); // queries
	if (!inP2O)
	{
		cout << "Error opening input file: " << q2o << endl;
		cout << "All quieries presumed from the same organism" << endl;
		exit(1);
	}
	else
	{
		populate(&prot2org, &orgHash, inP2O);
	}

	populate(&dom2org, &dom2len, &dom2cog, &domHash, &orgHash, &enCog, &cogDescript, inDat);

	map<string, int>::iterator m_p = cogHash.begin();
	while (m_p != cogHash.end())
	{
		cogHashR.insert(pair<int, string>((*m_p).second, (*m_p).first));
		m_p++;
	}

/*	ifstream inCogs(cogs.c_str());
	if (!inCogs)
	{
		cout << "Error opening file: " << cogs << endl;
		cout << "Output will not include cog descriptions" << endl;
		populate(&cogDescript, &enCog);
	}
	else
	{
		populate(&cogDescript, &enCog, inCogs);
	}*/

	ifstream inTaxid("genIdNew2taxIdOld.csv");
	if (!inTaxid)
	{
//		cerr << "Error opening file: genIdNew2taxIdOld.csv" << endl;
	}
	else
	{
		populate(&genIdNew2taxidOld, inTaxid);
	}

	// Change directory here
	myOpenDir(dirpath.c_str());

	HitSet hs(
		overlapThreshold, 
		cogOverlapThreshold,
		numThreshold, 
		totalNumThreshold,
		expectThreshold, 
		orgThreshold,
		gapThreshold,
		weightExp,
		lengthExp,
		rankExp,
		rankOrg,
		&domHash,
		&selfHits, 
		&dom2org, 
		&dom2cog, 
		&prot2org,
		&offsets, 
		&cogDescript,
		&dom2len,
		&genIdNew2taxidOld);

	while(getline(inHitFile, buffer))
	{
		curRec = buffer;
		if (hs.insert(curRec) == false)
		{
			hs.print(&fOut, &outFileConfl);
			hs.clear();
			hs.insert(curRec);
		}
	}
	hs.assign();
	hs.print(&fOut, &outFileConfl);

	int length = 0;
	string self;
	map<string, string>::iterator shp = selfHits.begin();

	int id = 0;
	while (shp != selfHits.end()) // some proteins may not have self hit (too short for BLAST)
	{
		if ((*shp).second != "")
		{
			self = (*shp).second;
			length = atoi(self.substr(0, self.find(',')).c_str());	
			id = atoi((*shp).first.c_str());
			if (hs.getName(id) != "")
			{
				fOut << hs.getName(id) << ',' << length << ',' << '1' << ',' << length << ',' << 0 << ',' << -1 << endl;
			}
		}
		shp++;
	}

	return 0;
}
