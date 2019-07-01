// global functions
#ifndef GLOBALS_H
#define GLOBALS_H

#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include "enum.h"

int		split		(const string str, const string sepStr, vector<string>* v);
void	makeOffsets	(map<string, long>* offsets, ifstream& inFile);
void	populate	(map<int, int>* m1, map<int, int>* m2, map<int, int>* m3, ifstream& inFile);
void	populate	(map<int, int>* m1, map<int, int>* m2, map<int, int>* m3, map<string, int>* domHash, map<string, int>* orgHash, map<string, int>* cogHash, ifstream& inFile);
void	populate	(map<int, int>* m1, map<int, int>* m2, map<int, int>* m3, map<string, int>* domHash, map<string, int>* orgHash, Enumerator* enCog, map<int, string>* cogDescript, ifstream& inFile);
void	populate	(map<int, int>* p2ent, ifstream& inFile);
void	populate	(map<string, string>* selfHits, ifstream& inFile);
void	populate	(multimap<string, int>* q2t, ifstream& inFile);
void	populate	(map<int, string>* cogs, ifstream& inFile);
void	populate	(map<string, int>* p2o, map<string, int>* orgHash, ifstream& inFile);

int split(const string str, const string sepStr, vector<string>* v)
{
    string	token;
    int		nToken = 0;
    size_t	index;
    size_t	posSep;
    size_t	strLen;
    size_t	tokenLen;

	(*v).clear();

    strLen = str.length();

    index = str.find_first_not_of(sepStr);// skip leading separators

    do
    {
	posSep = str.find_first_of(sepStr, index);
	if (posSep == string::npos)
	{// separator not found
	    tokenLen = strLen - index;
	}
	else
	{
	    tokenLen = posSep - index;
	}
	token = str.substr(index, tokenLen);
	(*v).push_back(token);
	nToken++;
	//index = str.find_first_not_of(sepStr, posSep + 1);//skip all contin. separators
	index = posSep + 1;
    }
	while(posSep != string::npos && index != string::npos);
    
    return nToken;
}


void makeOffsets(map<string, long>* offsets, ifstream& inFile)
{
	long offset;
	string str;
	string curProt;
	string prot;
	string curSubj;
	int indx;
	
	offset = inFile.tellg();
	while (getline(inFile, str))
	{
		indx = str.find_first_of(",");
		prot = str.substr(0, indx);
		curSubj = str.substr(indx + 1);

		if (prot != "")
		{
			if (curProt != prot)
			{
				(*offsets).insert(pair<string, long>(prot, offset));
				curProt = prot;
			}
		}
		offset = inFile.tellg();
	}
}

void populate(map<int, int>* m1, map<int, int>* m2, map<int, int>* m3, ifstream& inFile)
{
	vector<string> v;
	string str;
	int dom;
	int org;
	int len;
	int cog;

	while (getline(inFile, str))
	{
		split(str, ",", &v);
		if (v.size() < 7)
		{
			cerr << "incorrect line: " << endl << str << endl;
			continue;
		}
		dom = atoi(v[0].c_str());
		org = atoi(v[1].c_str());
		len = atoi(v[3].c_str());
		cog = atoi(v[6].c_str());

		(*m1).insert(pair<int, int>(dom, org));
		(*m2).insert(pair<int, int>(dom, len));
		(*m3).insert(pair<int, int>(dom, cog));
	}
}

void populate(map<int, int>* m1, map<int, int>* m2, map<int, int>* m3, map<string, int>* domHash, map<string, int>* orgHash, map<string, int>* cogHash, ifstream& inFile)
{
	vector<string> v;
	string str;
	string dom;
	int domID;
	string org;
	int orgID;
	int orgCount = 1;
	string cog;
	int cogID;
	int cogCount = 1;
	int len;

	while (getline(inFile, str))
	{
		split(str, ",", &v);
		if (v.size() < 7)
		{
			cerr << "incorrect line: " << endl << str << endl;
			continue;
		}

		org = v[1];
		if (((*orgHash).insert(pair<string, int>(org, orgCount))).second == true)
		{
			orgID = orgCount;
			orgCount++;
		}
		else
		{
			orgID = (*(*orgHash).find(org)).second;
		}

		len = atoi(v[3].c_str());
		
		cog = v[6];
		if (((*cogHash).insert(pair<string, int>(cog, cogCount))).second == true)
		{
			cogID = cogCount;
			cogCount++;
		}
		else
		{
			cogID = (*(*cogHash).find(cog)).second;
		}

		dom = v[0];
		if((*domHash).find(dom) != (*domHash).end())
		{
			domID = (*(*domHash).find(dom)).second;	
			(*m1).insert(pair<int, int>(domID, orgID));
			(*m2).insert(pair<int, int>(domID, len));
			(*m3).insert(pair<int, int>(domID, cogID));
		}
	}
}

void populate(map<int, int>* p2ent, ifstream& inFile)
{
	vector<string> v;
	string str;
	int prot;	
	int ent;
	
	while (getline(inFile, str))
	{
		split(str, ",", &v);
		prot = atoi(v[0].c_str());
		ent = atoi(v[1].c_str());
		
		(*p2ent).insert(pair<int, int>(prot, ent));
	}
}


void populate(map<string, string>* selfHits, ifstream& inFile)
{
	string str;
	string query;
	string hit;

	while (getline(inFile, str))
	{
		query = str.substr(0, str.find_first_of(','));
		hit = str.substr(str.find_first_of(',') + 1);
		(*selfHits).insert(pair<string, string>(query, hit));
	}
}

void populate(multimap<string, int>* q2t, ifstream& inFile)
{
	vector<string> v;
	string str;
	string query;	
	int target;
	
	while (getline(inFile, str))
	{
		split(str, ",", &v);
		query = v[0];
		target = atoi(v[1].c_str());
		
		(*q2t).insert(pair<string, int>(query, target));
	}
}

void populate(map<string, int>* p2o, ifstream& inFile)
{
	vector<string> v;
	string str;
	string query;	
	int org;
	
	while (getline(inFile, str))
	{
		split(str, ",", &v);
		query = v[0];
		org = atoi(v[1].c_str());
		
		(*p2o).insert(pair<string, int>(query, org));
	}
}

void populate(map<string, int>* p2o, map<string, int>* orgHash, ifstream& inFile)
{
	vector<string> v;
	string str;
	string query;	
	string org;
	int orgID;
	int orgCount = (*orgHash).size() + 1;
	
	while (getline(inFile, str))
	{
		split(str, ",", &v);
		query = v[0];
		if (v.size() == 1)
			org = "unknown";
		else
			org = v[1];

		if (((*orgHash).insert(pair<string, int>(org, orgCount))).second == true)
		{
			orgID = orgCount;
			orgCount++;
		}
		else
		{
			orgID = (*(*orgHash).find(org)).second;
		}
		
		(*p2o).insert(pair<string, int>(query, orgID));
	}
}

void populate(map<int, string>* cogs, ifstream& inFile)
{
	string str;
	string descript;
	int id;

	while (getline(inFile, str))
	{
		id = atoi(str.substr(0, str.find_first_of(',')).c_str());
		descript = str.substr(str.find_first_of(',') + 1);
		(*cogs).insert(pair<int, string>(id, descript));
	}

}

void populate(map<int, string>* cogs, Enumerator* enCog, ifstream& inFile)
{
	string str;
	string descript;
	string name;
	int id;

	while (getline(inFile, str))
	{
		name = str.substr(0, str.find_first_of(','));
		id = (*enCog).insert(name);
		descript = str.substr(str.find_first_of(',') + 1);
		(*cogs).insert(pair<int, string>(id, name + ',' + descript));
	}
}

void populate(map<int, string>*cogs, Enumerator* enCog)
{
	;
}

void populate(map<int, int>* m1, map<int, int>* m2, map<int, int>* m3, map<string, int>* domHash, map<string, int>* orgHash, Enumerator* enCog, map<int, string>* cogDescript, ifstream& inFile)
{
	vector<string> v;
	string str;
	string dom;
	int domID = 0;
	string org;
	int orgID;
	int orgCount = 1;
	string cog;
	int cogID;
	int len;

	while (getline(inFile, str))
	{
		split(str, ",", &v);
		if (v.size() < 7)
		{
			cerr << "incorrect line: " << endl << str << endl;
			continue;
		}

		org = v[1];
		if (((*orgHash).insert(pair<string, int>(org, orgCount))).second == true)
		{
			orgID = orgCount;
			orgCount++;
		}
		else
		{
			orgID = (*(*orgHash).find(org)).second;
		}

		len = atoi(v[3].c_str());
		
		cog = v[6];

		cogID = (*enCog).insert(cog);

		(*cogDescript).insert(pair<int, string>(cogID, cog));

		dom = v[0];
		if ((*domHash).size() == 0)
		{
			domID = atoi(dom.c_str());	
		}
		else if((*domHash).find(dom) != (*domHash).end())
		{
			domID = (*(*domHash).find(dom)).second;	
		}
		if (domID > 0)
		{
			(*m1).insert(pair<int, int>(domID, orgID));
			(*m2).insert(pair<int, int>(domID, len));
			(*m3).insert(pair<int, int>(domID, cogID));
		}
	}
}

#endif
