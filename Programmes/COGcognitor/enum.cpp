#include "enum.h"
using namespace std;

Enumerator::Enumerator(const int i)
{
	counter = i;
	stringKey.insert(pair<string, int>("NULL", -1));
	stringKey.insert(pair<string, int>("null", -1));
	stringKey.insert(pair<string, int>("Null", -1));
	stringKey.insert(pair<string, int>("0", -1));
	stringKey.insert(pair<string, int>("-1", -1));
	stringKey.insert(pair<string, int>("", -1));
}

int Enumerator::insert(const string name)
{
	if ((stringKey.insert(pair<string, int>(name, counter))).second == true)
	{ //first instance of the name
		intKey.insert(pair<int, string>(counter, name));
		counter++;
		return(counter - 1);
	}
	else
	{ //the name is already in Enumerator
		return findNumber(name);
	}
}

int Enumerator::findNumber(const string name)
{
	if	(stringKey.find(name) != stringKey.end())
		return (*stringKey.find(name)).second;
	else
		return(-1);
}

string Enumerator::findName(const int number)
{
	if	(intKey.find(number) != intKey.end())
		return (*intKey.find(number)).second;
	else
		return "";
}

Printer::Printer(Enumerator* en1, Enumerator* en2, ofstream* ofile)
{
	enCog = en1;
	enOrg = en2;
	outFile = ofile;
}

string Printer::print(vector<string>* sec)
{
	string cogID;

	vector<string>::iterator p_v;
	p_v = (*sec).begin();
	while (p_v != (*sec).end())
	{
		str = (*p_v);
		cogID = str.substr(str.find_last_of(','));
		str = str.substr(0, str.find_last_of(',')) + (*enCog).findName(atoi(cogID.c_str()));
		cout << str << endl;
		fullstr += str;
		p_v++;
	}
	return fullstr;
}
