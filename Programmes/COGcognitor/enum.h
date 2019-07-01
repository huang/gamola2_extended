#ifndef ENUM_H
#define ENUM_H

#include <cstdlib>
#include <string>
#include <map>
#include <vector>
#include <iostream>

using namespace std;

class Enumerator
{
private:
	int counter;
	map<string, int> stringKey;
	map<int, string> intKey;
public:
	Enumerator(const int i);
	int insert(const string);
	int findNumber(const string);
	string findName(const int);
};

class Printer
{
private:
	Enumerator* enCog;
	Enumerator* enOrg;
	ofstream* outFile;
	string str;
	string fullstr;
public:
	Printer(Enumerator* en1, Enumerator* en2, ofstream* ofile);
	string print(vector<string>* sec);
};

#endif
