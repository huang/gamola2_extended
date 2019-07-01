//#define DEBUG

#include "cognitor.h"
#include "math.h"

extern int split(const string str, const string sepStr, vector<string>* v);

void Genomes::myInsert(const int org) 
{
	genomes.insert(org);
	genomesUnique.insert(org);
}

int Genomes::count(const int org)
{
	return genomes.count(org); 
}

int Genomes::size()
{
	return genomes.size();
}

int Genomes::countUnique()
{
	return genomesUnique.size();
}

void CogAssignments::setExp(const double expWeight, const double expLen, const double expRank, const double expOrgRank)
{
	exp.weight	= expWeight;
	exp.length	= expLen;
	exp.rank	= expRank;
	exp.orgRank = expOrgRank;
}


double CogAssignments::killOverlaps(const string& query, ofstream* outFileConfl) // recursion
{
	unsigned int i;
	double maxRank = 0.0;
	int bestSector = -1;
	int start;
	int end;
	float overlap;
	int len;
	int iLen;
	int minLen;
	double rank;

	//*** find sector with the maximum rank from survived sectors ***
	for (i = 0; i < asnts.size(); i++)
	{
		if (maxRank < asnts[i].rank && asnts[i].good == true)
		{
			maxRank = asnts[i].rank;
			bestSector = i;

		}
	}
	if (bestSector == -1)
		return 0.0;

	//*** mark overlapping with 'the best' (including 'the best' itself) sectors as 'false' ***
	if (maxRank != 0.0)
	{
		//////////////////////////////////
		/*start = asnts[bestSector].startExp; // overlap based on extended coordinates
		end = asnts[bestSector].endExp;
		rank = asnts[bestSector].rank;
		len = end - start;*/

		start = asnts[bestSector].start; // or overlap based on core coordinates
		end = asnts[bestSector].end;
		rank = asnts[bestSector].rank;
		len = end - start;
		////////////////////////////////////

//		final.insert(pair<int, Assignment>(start, asnts[bestSector])); // save a copy of 'the best'

		for (i = 0; i < asnts.size(); i++)
		{
			if ((signed)i == bestSector || asnts[i].good == false) // keep C++ quiet
				continue;

			//////////////////////////////////////////////

			/*iLen = asnts[i].endExp - asnts[i].startExp; 
			minLen = iLen < len ? iLen : len;

			if (asnts[i].startExp < start)
			{
				overlap = (float)(asnts[i].endExp - start)/minLen;
			}
			else if (asnts[i].endExp > end)
			{
				overlap = (float)(end - asnts[i].startExp)/minLen;
			}
			else
			{
				overlap = 1;
			}*/

			iLen = asnts[i].end - asnts[i].start; 
			minLen = iLen < len ? iLen : len;

			if (asnts[i].start < start)
			{
				overlap = (float)(asnts[i].end - start)/minLen;
			}
			else if (asnts[i].end > end)
			{
				overlap = (float)(end - asnts[i].start)/minLen;
			}
			else
			{
				overlap = 1;
			}
			/////////////////////////////////////////////
			
			if (overlap > 0.25) //!!!
			{
				if (asnts[i].rank/rank > 0.3 && asnts[bestSector].cog != asnts[i].cog)
				{
					*outFileConfl << query << "," << asnts[bestSector].cog << "," << start << "," << end << "," <<  asnts[bestSector].rank << "," << asnts[i].cog << "," << asnts[i].start << "," << asnts[i].end << "," << asnts[i].rank << endl;
				}
				
				if (asnts[bestSector].cog == asnts[i].cog && (asnts[bestSector].end >= asnts[i].start || asnts[bestSector].start >= asnts[i].end))
				{
					asnts[bestSector].startExp = asnts[bestSector].startExp <= asnts[i].startExp ? asnts[bestSector].startExp : asnts[i].startExp;
					asnts[bestSector].endExp = asnts[bestSector].endExp >= asnts[i].endExp ? asnts[bestSector].endExp : asnts[i].endExp;
					asnts[bestSector].start = asnts[bestSector].start <= asnts[i].start ? asnts[bestSector].start : asnts[i].start;
					asnts[bestSector].end = asnts[bestSector].end >= asnts[i].end ? asnts[bestSector].end : asnts[i].end;

					asnts[i].good = false;
					return maxRank;
				}
				
				asnts[i].good = false;
			}

		}
		final.insert(pair<int, Assignment>(asnts[bestSector].start, asnts[bestSector])); // save a copy of 'the best'
		asnts[bestSector].good = false;
	}
	return maxRank;
}

void CogAssignments::combineHits(const double cogOverlap, const int length)
{
	multimap<int, int> order; // reorder by start
	multimap<double, int> blockByR; // sort inside COG-sector by weights
	unsigned int i;
	int k;
	int l;
	int n;
	int curCog = 0;
	int cog;
	int minStart = 0;
	int start;
	int maxEnd = 0;
	int end;
	int overlap;
	int curLen;
	int len;
	int minLen;
	double rank;
	int count;
	unsigned int last;
	int curStart = 0;
	int curEnd = 0;
	double cumRank;
	int bet;
	int starti;
	set<int> head;
	set<int> tail;
	int startExp;
	int endExp;
	int mT;
	int mH;
	
	for (i = 0; i < asnts.size(); i++)	// presorted by cog, start
	{
		cog = asnts[i].cog;
		starti = asnts[i].start;

		if (curCog == 0) 
		{
			curCog = cog;
		}
		
		if (curCog == cog)
		{
			order.insert(pair<int, int>(starti, i)); // collect all HSPs for the cog and sort them by start
		}

		if (curCog != cog || i == asnts.size() - 1) // next cog found; process all HSPs for the cog
		{
			curCog = cog;

			last = 0;
			curStart = 0;
			curEnd = 0;
			multimap<int, int>::iterator p_m = order.begin();
			while (p_m != order.end())
			{
				if (order.size() == 1)
				{
					p_m++;
					continue;
				}

				last++;
				k = (*p_m).second;
				start = asnts[k].start;
				end = asnts[k].end;
				len = end - start + 1;
				curStart = curStart == 0 ? start : curStart;
				curEnd = curEnd == 0 ? end : curEnd;
				curLen = curEnd - curStart + 1;
				minLen = len < curLen ? len : curLen;
				overlap = curEnd - start + 1; // ok for now
				curEnd = curEnd < end ? end : curEnd;
				rank = asnts[k].rank;


				if ((float)overlap/minLen >= cogOverlap || order.size() == 1) // first hit or hit in same block: accumulate all HSPs for the block and sort them by rank
				{					
					blockByR.insert(pair<double, int>(rank, k));
				}
				if ((float)overlap/minLen < cogOverlap || (last == order.size() && last != 1))  // next block found OR the current block is the last one in cog
				{
					count = 1;
					curStart = 0;
					curEnd = 0;

					if (blockByR.size() == 1)
					{
						p_m++;
						blockByR.clear();
						continue;
					}
					
					multimap<double, int>::reverse_iterator pr_mm = blockByR.rbegin();
					while (pr_mm != blockByR.rend())
					{
						n = (*pr_mm).second;
						//asnts[n].rank = asnts[n].rank/pow(count, exp.rank);
						count++;
						pr_mm++;
					}
					
					multimap<double, int>::reverse_iterator pr1_mm = blockByR.rbegin();					
					bet = (*pr1_mm).second;
					minStart = asnts[bet].start;
					maxEnd = asnts[bet].end;
					cumRank = asnts[bet].rank;
					head.insert(asnts[bet].start - asnts[bet].startExp);
					tail.insert(asnts[bet].endExp - asnts[bet].end);
					pr1_mm++;

					int counter = 0;
					while (pr1_mm != blockByR.rend())
					{
						counter++;
						n = (*pr1_mm).second;
						if (++counter > 15)
						{
							asnts[n].good = false;
							pr1_mm++;
							continue;
						}
						start = asnts[n].start;
						end = asnts[n].end;
						startExp = asnts[n].startExp;
						endExp = asnts[n].endExp;
						head.insert(start - startExp);
						tail.insert(endExp - end);
						minStart = minStart > start ? start : minStart;
						maxEnd = maxEnd < end ? end : maxEnd;
						asnts[bet].start = minStart;
						asnts[bet].end = maxEnd;

//						rank = asnts[n].rank;
						cumRank += asnts[n].rank;
						asnts[bet].rank = cumRank;
						asnts[n].good = false;
						pr1_mm++;

					}


					if (asnts[bet].cog != -1)
					{

						mH = head.size()/2 + 1;
						mT = tail.size()/2 + 1;

						set<int>::iterator psh = head.begin();
						for (l = 1; l != mH; l++)
						{
							psh++;
						}
						asnts[bet].startExp = (asnts[bet].start - *psh < 1 ? 1 : asnts[bet].start - *psh);


						set<int>::iterator pst = tail.begin();
						for (l = 1; l != mT; l++)
						{
							pst++;
						}
						asnts[bet].endExp = (asnts[bet].end + *pst > length ? length : asnts[bet].end + *pst);
					}
					head.clear();
					tail.clear();
					blockByR.clear();
					//blockByR.insert(pair<double, int>(rank, k));

				}
				p_m++;
			}// while
			order.clear();
			order.insert(pair<int, int>(starti, i));

		} // if
	} // for
}

int CogAssignments::print(const string& query, const int length, const int gap, const double cogOverlap, map<int, string>* cogDescript, map<int, string>* keys, ofstream* fOut, ofstream* outFileConfl)
{
	float minRank = 0.0;
	int curStart = 1;
	int curEnd = 1;
	int curStartCore = 1;
	int curEndCore = 1;
	int curCog = 0;
	int start;
	int end;
	int startCore;
	int endCore;
	int cog;
	double curRank = 0;
	double rank;
	int last = 0;
	int i = 0;

	combineHits(cogOverlap, length); 

	if (asnts.size() == 0)
	{
		printSector(query, length, 1, length, 0, -1, cogDescript, keys, fOut);
		return -2;
	}


	while (0.0 != killOverlaps(query, outFileConfl));

	last = final.size();

	multimap<int, Assignment>::iterator p_m = final.begin();

	while(p_m != final.end())
	{
		startCore = (*p_m).second.start;
		endCore = (*p_m).second.end;
		start = (*p_m).second.startExp;
		end = (*p_m).second.endExp;
		cog = (*p_m).second.cog;
		rank = (*p_m).second.rank;

		if (rank < minRank)
		{
			p_m++;
			continue;
		}

		if (i == 0) // first relevant domain
		{
			curCog = cog;
			curStart = start < gap || curCog == -1 ? 1 : start;
//			curStart = 1; // no heads
			curStartCore = startCore < gap ? 1 : startCore;

			if (curStart >= gap)
			{ // print head
				printSector(query, length, 1, curStart - 1, 0, -1, cogDescript, keys, fOut);
			}

			curEnd = end;
			curEndCore = endCore;
			curRank = rank;
			i++;
		}
		
		else if (curCog == cog)// successive sectors of the same COG
		{
			if (start - curEnd > gap) // make new domain from gap
			{
				printSector(query, length, curStart, curEnd, curRank, curCog, cogDescript, keys, fOut);
				printSector(query, length, curEnd + 1, start - 1, 0, -1, cogDescript, keys, fOut);
				curRank = rank;
				curStart = start;
				curEnd = end;
				curStartCore = startCore;
				curEndCore = endCore;
			}
			else
			{
				curEnd = end;
				curEndCore = endCore;
				curRank += rank;				
			}
		}
		else // new COG
		{
			if (start - curEnd == 1) // no gap, no overlap
			{
				printSector(query, length, curStart, curEnd, curRank, curCog, cogDescript, keys, fOut/*, gap, i, last*/);
				curCog = cog;
				curRank = rank;
				curStart = start;
				curEnd = end;
				curStartCore = startCore;
				curEndCore = endCore;
			}
			else
			{
				if (start - curEnd > gap) // make a new free domain from gap
				{
					printSector(query, length, curStart, curEnd, curRank, curCog, cogDescript, keys, fOut);
					printSector(query, length, curEnd + 1, start - 1, 0, -1, cogDescript, keys, fOut);
					curCog = cog;
					curRank = rank;
					curStart = start;
					curEnd = end;
					curStartCore = startCore;
					curEndCore = endCore;
				}
				else if (curEnd < startCore && curEndCore < start) // gap exists and has to be included in domains or overlap exists
				{
					curEnd += (start - curEnd)/2;
					printSector(query, length, curStart, curEnd, curRank, curCog, cogDescript, keys, fOut/*, gap, i, last*/);
					curCog = cog;
					curRank = rank;
					curStart = curEnd + 1;
					curEnd = end;
					curStartCore = startCore;
					curEndCore = endCore;
//					if (curStart > curEnd)
//						cerr << "lazha: " << query << ' ' << start << ' ' << end << endl;
				}
				else // cores overlapped
				{
					curEnd = curEndCore + (startCore - curEndCore)/2;
					printSector(query, length, curStart, curEnd, curRank, curCog, cogDescript, keys, fOut);
					curCog = cog;
					curRank = rank;
					curStart = curEnd + 1;
					curEnd = end;
					curStartCore = startCore;
					curEndCore = endCore;
				}
			}
		}
		p_m++;
	}
	if (curRank >= minRank)
	{
		curEnd = length - curEnd < gap ? length : curEnd;
		printSector(query, length, curStart, curEnd, curRank, curCog, cogDescript, keys, fOut);
	}
	
	if (length - curEnd >= gap)
	{ // print tail or the whole protein if not assigned
		if (curEnd == 1)
		{
			curEnd = 0; //correction for not asigned proteins
		}
		printSector(query, length, curEnd + 1, length, 0, -1, cogDescript, keys, fOut);
	}

	return curCog;
}

void CogAssignments::printSector(const string& query, const int length, const int curStart, int curEnd, const double curRank, const int curCog, map<int, string>* cogDescript, map<int, string>* keys, ofstream* fOut/*, const int gap, const int i, const const last*/)
{
	string prot;
	map<int, string>::iterator mp;

	mp = (*keys).find(atoi(query.c_str()));
	if (mp != (*keys).end())
	{
		prot = (*mp).second;
		*fOut << prot << ',' << length << ',' << curStart << ',' << curEnd << ',' << curRank << ',';
	}

	if ((*cogDescript).find(curCog) == (*cogDescript).end())
	{
		*fOut << curCog << endl;
	}
	else
	{
		*fOut << (*(*cogDescript).find(curCog)).second << endl;
	}
}

void CogAssignments::clear()
{
	asnts.clear();
	final.clear();
}

void CogAssignments::insert(const int cog, const int start, const int startExp, const int end, const int endExp, const double weight, bool isGood)
{
		Assignment* as = new Assignment;
		as->cog = cog;
		as->start = start;
		as->startExp = startExp;
		as->end = end;
		as->endExp = endExp;
		as->rank = weight;
		as->good = isGood;
		asnts.push_back(*as);
		delete as;
}

double CogAssignments::credit(const double scoreW, const double lenW, const double rank, const double rankOrg)
{
	return (9*pow(scoreW, exp.weight) + 1)*(pow(lenW, 1/exp.length))/((pow(rank, exp.rank)*(pow(rankOrg, exp.orgRank))));
}


void Filter::setup(const string& query, map<string, long>* offsets)
{
	long offset;
	string prot;
	string str;
	string subject;
	
	static ifstream inQ2S("query2subject.csv");
	if (!inQ2S)
	{
		cout << query << ": File query2subject.csv not found for setting the mask; all hits will be processed" << endl;
	}	

	mask.clear();
	
	if ((*offsets).find(query) != (*offsets).end())
	{
		offset = (*(*offsets).find(query)).second;

		inQ2S.clear();
		inQ2S.seekg(offset);
		while (getline(inQ2S, str))
		{
			prot = str.substr(0, str.find_first_of(','));
			if (prot == query)
			{
				subject = str.substr(str.find_first_of(',') + 1);
				mask.insert(atoi(subject.c_str()));
			}
			else
			{
				break;
			}
		}
	}
}

bool Filter::isRelevant(const int subj)
{
	if (mask.find(subj) != mask.end()) // || mask.size() == 0) if there is no mask still process
	{
		return true;
	}
	return false;
}

HitSet::HitSet(	const double			overlap,
				const double			cogOverlap,
				const int				num, 
				const int				totalNum,
				const double			expect, 
				const int				orgs,
				const int				gap,
				const double			expWeight,
				const double			expLen,
				const double			expRank,
				const double			expOrgRank,
				map<string, int>*		domHash,
				map<string, string>*	selfHits, 
				map<int, int>*			dom2org, 
				map<int, int>*			dom2cog, 
				map<string, int>*		prot2org,
				map<string, long>*		offsets, 
				map<int, string>*		cogDescript,
				map<int, int>*			dom2len,
				map<int, int>*			genIdNew2taxidOld)
		:organism(0), 
		length(0)
{
	thresholds.overlap		= (overlap > 1 || overlap < 0) ? 1 : overlap;
	thresholds.cogOverlap	= (cogOverlap > 1 || cogOverlap < 0) ? 1 : cogOverlap;
	thresholds.num			= num;
	thresholds.totalNum		= totalNum;
	thresholds.expect		= expect;
	thresholds.orgs			= orgs;
	thresholds.gap			= gap;
	db.domHash				= domHash;
	db.selfHits				= selfHits;
	db.dom2org				= dom2org;
	db.dom2cog				= dom2cog;
	db.prot2org				= prot2org;
	db.offsets				= offsets;
	db.cogDescript			= cogDescript;
	db.dom2len				= dom2len;
	db.genIdNew2taxidOld	= genIdNew2taxidOld;

	cogs.setExp(expWeight, expLen, expRank, expOrgRank);

	map<string, int>::iterator mp = (*db.domHash).begin();
	while (mp != (*db.domHash).end())
	{
		keys.insert(pair<int, string>((*mp).second, (*mp).first));
		mp++;
	}
	(*db.domHash).clear();
}

bool HitSet::insert(const string& curRec) // returns false if trying insert record for next query
{
	string curQuery;
	int curSubject;
	int curOrg;
	int curCog;
	int start;
	int end;
	int tStart;
	int tEnd;
	int startExp;
	int endExp;
	static int subject;
	static multimap<int, Hsp> pluralHit;
	string self;
	vector<string> v;
	static multiset<int> total; // to count nonscpecific hits in each cog

	split(curRec, ",", &v); // format: qEntId bigint,tEntId bigint,qStart int,qEnd int, tStart int, tEnd int, evalue float,score float

	if (v.size() < 8)
	{
		return true; // skip incorrect lines
	}

	curQuery = v[0];
//	if ((*db.dom2cog).find(atoi(curQuery.c_str())) != (*db.dom2cog).end())
//	{
//		return true; // skip queries which have a cog (when 'all on all' blasts)
//	}

	curSubject = atoi(v[1].c_str());

	if ((*db.dom2org).find(curSubject) != (*db.dom2org).end())
	{
		curOrg = (*(*db.dom2org).find(curSubject)).second;
	}
	else
	{
		curOrg = -1;
	}

	if (orgRanks.find(curOrg) == orgRanks.end())
	{	
		orgRanks.insert(pair<int, int>(curOrg, orgRanks.size() + 1));
	}

	if ((*db.dom2cog).find(curSubject) != (*db.dom2cog).end())
	{
		curCog = (*(*db.dom2cog).find(curSubject)).second;
	}
/*	else	
	{
		curCog = -1;
		return true; //!!!!!!!!!! ONLY FOR PROCESSING 'FREE PROTEINS' NOT NEW GENOMES!!!!!!!!!!!!
	}*/

/*	if (atoi(query.c_str()) == curSubject)
	{
		return true; // for selfevaluation only
	}*/

	if (query == "") // first record
	{
		query = curQuery;

		subject = curSubject;
	
		if ((*db.prot2org).find(query) != (*db.prot2org).end())
		{
			organism = (*(*db.prot2org).find(query)).second;
		}
		else
		{
			organism = -1;
		}

		setGenomes();

		filter.setup(query, db.offsets);

		if ((*db.selfHits).find(query) != (*db.selfHits).end())
		{
			self = (*(*db.selfHits).find(query)).second;
			length = atoi(self.substr(0, self.find(',')).c_str());
			score = atof(self.substr(self.find(',') + 1).c_str());
		}
		else
		{
			length = 100;
			score = 1;
			cout << "self hit for " << query << " not found" << endl;
		}
	}


	if (query == curQuery)
	{
		if (relevantGenomes.find(curOrg) == relevantGenomes.end()) 
		{
//			return true; // TURN IT ON IF YOU MUST NOT PROCESS ALL TARGET GENOMES
		}

		if (howMuch.find(curCog) == howMuch.end()) // first hit in this cog
		{
			Genomes* g = new Genomes;
			(*g).myInsert(curOrg);
			howMuch.insert(pair<int, Genomes>(curCog, *g));
			total.insert(curCog);


			delete g;
		}

		else
		{
			(*howMuch.find(curCog)).second.myInsert(curOrg);
			total.insert(curCog);

		}
		if (total.count(curCog) > thresholds.totalNum || (*howMuch.find(curCog)).second.count(curOrg) > thresholds.num)
		{
			return true;
		}

		if (subject != curSubject && pluralHit.size() > 0) // pluralHit is full
		{
			{
				shrinkAndInsert(&pluralHit);
			}
			pluralHit.clear();
		}

		if(filter.isRelevant(curSubject) == true)
		{
			Hsp* hsp = new Hsp;
			hsp->subject = curSubject;
			hsp->org = curOrg;
			hsp->cog = curCog;//(*(*db.dom2cog).find(curSubject)).second;
			hsp->target = atoi(query.c_str());//atoi(v[0].c_str()); // rename to hsp->query, not subject
			start = atoi(v[2].c_str());
			hsp->start = start;
			end = atoi(v[3].c_str());
			hsp->end = end;
			tStart = atoi(v[4].c_str());
			hsp->tStart = tStart;
			tEnd = atoi(v[5].c_str());
			hsp->tEnd = tEnd;
			hsp->evalue = atof(v[6].c_str());
			hsp->score = atof(v[7].c_str());
			startExp = start - getPartLen(curSubject, tStart, tEnd, 'h');
			hsp->startExp = startExp < 1? 1 : startExp;
			endExp = end + getPartLen(curSubject, tStart, tEnd, 't');
			hsp->endExp = endExp > length? length : endExp;

			if (hsp->evalue < thresholds.expect)
			{
				pluralHit.insert(pair<int, Hsp>(hsp->start, *hsp));

				
			}
			delete hsp;

		}
		subject = curSubject;
		return true;
	} // if (query == curQuery)

	else // HitSet completed
	{
		shrinkAndInsert(&pluralHit);
		assign();

		relevantGenomes.clear();
		pluralHit.clear();
		howMuch.clear();
		total.clear();
		subject = 0;

		return false;
	}
}

void HitSet::assign()
{
	int i;
	int j;
	int cog;
	int curCog = 0;
	int len;
	int start;
	int startExp;
	double iscore;
//	int curStart = 0;
	int end;
	int endExp;
	double scoreW;
	float lenW;
	int iorg;
	int lenT;
	int rank;
	int orgRank;

	list<int> hitsByCogStart;
	map<int, int> hits2rank;
//	map<int, int>::iterator p_r;
	multiset<int> orgs;

	multimap<int, int>::iterator p_H = hitsByCog.begin(); 

	multimap<int, int> cogHits; // start, i	or iscore, i
	multimap<int, int>::iterator p_m;
	multimap<int, int>::reverse_iterator pr_m;


	//*** resort hits by cog, start ***
	while(p_H != hitsByCog.end()) 
	{
		i = (*p_H).second;
		cog = hits[i].cog;
		start = hits[i].start;


		if (curCog != cog)
		{
			if (cogHits.size() != 0)
			{
				p_m = cogHits.begin();
				while(p_m != cogHits.end())
				{
					hitsByCogStart.push_back((*p_m).second);
					p_m++;
				}
				cogHits.clear();
			}
			curCog = cog;
		}
		cogHits.insert(pair<int, int>(start, i));
		p_H++;
	}

	p_m = cogHits.begin();
	while(p_m != cogHits.end())
	{
		hitsByCogStart.push_back((*p_m).second);
		p_m++;
	}
	cogHits.clear();

	//*** calculate ranks for hits ***
	curCog = 0;
	p_H = hitsByCog.begin(); 
	while(p_H != hitsByCog.end()) // order hits by cog, score
	{
		i = (*p_H).second;
		cog = hits[i].cog;
		iscore = hits[i].score;

		if (curCog != cog)
		{
			if (cogHits.size() != 0)
			{
				pr_m = cogHits.rbegin();
				while(pr_m != cogHits.rend())
				{
					iorg = hits[(*pr_m).second].org;
					orgs.insert(iorg);
					j = orgs.count(iorg);
					hits2rank.insert(pair<int, int>((*pr_m).second, j));
					pr_m++;
				}
				cogHits.clear();
				orgs.clear();
			}
			curCog = cog;
		}
		cogHits.insert(pair<int, int>(int(iscore), i)); //safe int(double) to keep C++ quiet
		p_H++;
	}

	p_m = cogHits.begin();
	while(p_m != cogHits.end())
	{
		iorg = hits[(*p_m).second].org;
		orgs.insert(iorg);
		j = orgs.count(iorg);
		hits2rank.insert(pair<int, int>((*p_m).second, j));
		p_m++;
	}
	cogHits.clear();
	orgs.clear();

//*** Assign weights to hits ***
	list<int>::iterator p_HCS = hitsByCogStart.begin(); 
	while(p_HCS != hitsByCogStart.end())
	{
		i = (*p_HCS);
		cog = hits[i].cog;

		if ((*howMuch.find(cog)).second.countUnique() >=  thresholds.orgs)
		{
			start = hits[i].start;
			startExp = hits[i].startExp;
			end = hits[i].end;		
			endExp = hits[i].endExp;
 			len = end - start;
//			org = hits[i].org;
			rank = (*hits2rank.find(i)).second;
			scoreW = (hits[i].score/len)/(score/length);
	  		scoreW = scoreW > 1 ? 1 : scoreW;

			if ((*db.dom2len).find(hits[i].subject) != (*db.dom2len).end())
			{
				lenT = (*(*db.dom2len).find(hits[i].subject)).second;
				lenW = (float)len/lenT;
				lenW = lenW > 1 ? 1 : lenW;
			}
			else
			{
				cout << "Not found length for target domain " << hits[i].subject << endl;
				lenW = 0.5;
			}

			orgRank = (*orgRanks.find(hits[i].org)).second;

			double cred = hits[i].score * cogs.credit(scoreW, lenW, rank, orgRank);
			cogs.insert(cog, start, startExp, end, endExp, cred, true);
		}

		p_HCS++;
	}

}

void HitSet::print(ofstream* fOut, ofstream* outFileConfl)
{
	map<string, string>::iterator shp;
//	HitSet& hs = returnObject();
	cogs.print(query, length, thresholds.gap, thresholds.cogOverlap, db.cogDescript, &keys, fOut, outFileConfl);
	shp = (*db.selfHits).find(query);
	if (shp != (*db.selfHits).end())
	{
		(*shp).second = "";
	}

}

HitSet& HitSet::returnObject()
{
	return *this;
}

string	HitSet::getName(const int id)
{
	if (keys.find(id) != keys.end())
		return (*(keys.find(id))).second;
	
	return "";
}

void HitSet::clear()
{
	query = "";
	organism = 0;
	length = 0;
	hits.clear();
	hitsByCog.clear();
	cogs.clear();
	orgRanks.clear();
}

void HitSet::shrinkAndInsert(multimap<int, Hsp>* pluralHit)
{
	int start = 0;
	int nextStart = 0;
	int end = 0;
	int nextEnd = 0;
	double score = 0;
	double nextScore = 0;
	int overlap = 0;
	int minLen = 0;
	int i = 0; // index of inserted HSP
	int cog;

	multimap<int, Hsp>::iterator p_pH = (*pluralHit).begin();
	while(p_pH != (*pluralHit).end())
	{
		if (start == 0) // first HSP
		{
			start = (*p_pH).first;
			end = (*p_pH).second.end;
			score = (*p_pH).second.score;
			hits.push_back((*p_pH).second); // i is still 0

			i = hits.size() - 1;
			cog = (*p_pH).second.cog;
			hitsByCog.insert(pair<int, int>(cog, i));
			
			p_pH++;
			continue;
		}
		
		nextStart = (*p_pH).first;
		nextEnd = (*p_pH).second.end;
		nextScore = (*p_pH).second.score;

		overlap = end - nextStart + 1;
		minLen = (end - start < nextEnd - nextStart) ? end - start : nextEnd - nextStart;

		if ((float)overlap/minLen < thresholds.overlap)
		{
			hits.push_back((*p_pH).second);
			
			i = hits.size() - 1;
			cog = (*p_pH).second.cog;
			
			hitsByCog.insert(pair<int, int>(cog, i));

			start = nextStart;
			end = nextEnd;
		}
		else if (nextScore > score)
		{
			hits.pop_back();
			hits.push_back((*p_pH).second);
			
			i = hits.size() - 1;
			cog = (*p_pH).second.cog;
			
			hitsByCog.insert(pair<int, int>(cog, i));

			start = nextStart;
			end = nextEnd;
		}
		p_pH++;
	}
}

void HitSet::setGenomes()
{
	map<int, int>::iterator p_m = (*db.genIdNew2taxidOld).find(organism);

//*** For update ***

/*	if (p_m != (*db.genIdNew2taxidOld).end())
	{
		relevantGenomes.insert((*p_m).second);
	}*/


//*** For cognitorising ***
/*	{
//		cerr << "Correspondent organism for current query " << query << " not found" << endl;
		p_m = (*db.genIdNew2taxidOld).begin();
		while (p_m != (*db.genIdNew2taxidOld).end())
		{
			if ((*p_m).second == 1012057) // exclude a genome
			{
				p_m++;
				continue;
			}
			relevantGenomes.insert((*p_m).second);
			p_m++;
		}
//		relevantGenomes.insert(-1); // free proteins are relevant
	}*/
}

int HitSet::getPartLen(const int target, const int tStart, const int tEnd, const int part) // part == h(head) || b(hit's body) || t(tail)
{
	int len;
	int tLen;

		switch (part)
		{
		case 'b':
			len = tEnd - tStart + 1;
			break;
		case 'h':
			len = tStart - 1;
			break;
		case 't':
			if ((*db.dom2len).find(target) != (*db.dom2len).end())
			{
				tLen = (*(*db.dom2len).find(target)).second;
			}
			else
			{
				cerr << "target " << target << " not found in dom2len" << endl;
				exit(1);
			}
			len = tLen - tEnd;
			len = len < 0 ? 0 : len;
			break;
		default:
			cerr << "Function getPartLen needs 'h', 'b' or 't' (head, body or tail) as a third parameter";
			exit(1);
		}
	
	return len;	
}
