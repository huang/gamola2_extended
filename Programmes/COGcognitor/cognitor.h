// COGNITOR.H
//#define DEBUG

#ifndef COGNITOR_H
#define COGNITOR_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <list>

#include "enum.h"

using namespace std; 

struct Hsp
{
	int		target;
	int		subject;
	int		org;
	int		cog;
	int		start; // in query
	int		end; //in query
	int		tStart;
	int		tEnd;
	int		startExp;
	int		endExp;
	double	score;
	double	evalue;
};

struct Assignment
{
	int		cog;
	int		start;
	int		end;
	int		startExp; //expanded start and end
	int		endExp;
	double	rank;
	bool	good;
};

struct Thresholds
{
	double	overlap; // minimum overlap between 2 sectors to delete one of them
	double	cogOverlap; // minimum overlap between 2 sectors of the same COG to be joined
	int		num; // how many organism-specific hits for each cog to be taken into account
	unsigned int		totalNum; // how many hits for each cog to be taken into account
	double	expect; // e-value cut-off
	int		orgs; // how many hits in different organisms have to be found for each cog
	int		gap; // maximum length of gap, head or tale to be included in COG sector
};

struct Exp
{
	double weight;
	double length;
	double rank;
	double orgRank;
};

struct DB
{
	map<string, int>*		domHash;
	map<string, string>*	selfHits;
	map<int, int>*			dom2org;
	map<int, int>*			dom2cog;
	map<string, int>*		prot2org;
	map<string, long>*		offsets;
	map<int, string>*		cogDescript;
	map<int, int>*			dom2len;
	map<int, int>*			genIdNew2taxidOld;
	map<int, string>*		cogHashR;
	Enumerator*				enCog;
};

class Genomes
{
private:
	set<int> genomesUnique;
	multiset<int> genomes;
public:
	void myInsert(const int org);
	int count(const int org);
	int countUnique();
	int size();
};

class CogAssignments
{
friend class HitSet;

public:
		void setExp		(const double expWeight,
						const double expLen,
						const double expRank,
						const double expOrgRank);

		int print	(const string& query, 
						const int length,
						const int gap,
						const double cogOverlap,
						map<int, string>* cogDescript,
						map<int, string>* keys,
						ofstream* fOut,
						ofstream* outFileConfl);
			
		void printSector(const	string&				query, 
						const	int					length, 
						const	int					curStart, 
								int					curEnd, 
						const	double				curRank, 
						const	int					curCog,
								map<int, string>*	cogDescript,
								map<int, string>*	keys,
								ofstream*			fOut);
			
		void clear		();
private:
	void	insert		(const int cog, const int start, const int startExp, const int end, const int endExp, const double weight, bool isGood);
	double  credit		(const double scoreW, const double lenW, const double rank, const double rankOrg);
	double	killOverlaps(const string& query, ofstream* outFileConfl);
	void	combineHits	(const double cogOverlap, const int length);


	Exp							exp;
	vector<Assignment>			asnts;
	multimap<int, Assignment>	final; // order by start
};

class Filter
{
public:
	void	setup		(const string& query, map<string, long>* offset);
	bool	isRelevant	(const int subject);
private:
	set<int> mask; // filtered hits for a given query
};

class HitSet // includes all or some HSPs for a given query
{
public:
			HitSet	(const double			overlap,
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
					map<int, int>*			genIdNew2taxidOld);
					

		bool	insert		(const string& curRec);
		void	assign		();
		void	print		(ofstream* fOut, ofstream* outFileConfl);
		void	clear		();	
		string	getName		(const int id);

	
private:
	void shrinkAndInsert	(multimap<int, Hsp>* pluralHit); // makes best minimal set of HSPs to specific subject
	void setGenomes			();
	HitSet&	returnObject	();
	int	getPartLen			(const int target, const int tStart, const int tEnd, const int part); // int == h(head) || b(body) || t(tail)


	Thresholds			thresholds;
	DB					db;
	map<int, string>	keys;
	string				query;
	int					organism;	// query org
	int					length;	// query lenght	
	double				score;	// score of selfHit
	Filter				filter;
	CogAssignments		cogs;
	vector<Hsp>			hits;
	set<int>			relevantGenomes; // hits in what genomes are relevant for the current query
	map<int, Genomes>	howMuch; // cog, genome
	map<int, int>		orgRanks; // genome, rank
	multimap<int, int>	hitsByCog;
};
#endif
