#ifndef OGRAPH
#define OGRAPH

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <cstdlib>
#include <unordered_map>
#include <stdlib.h>
#include <stdio.h>
//#include "bitarray.h"
#include <cstring>

using namespace std;

#pragma pack(2)
typedef struct KmerHashTable{
    unsigned long long int kmer;
    long int frequency;
    KmerHashTable(){
        kmer = 0;
        frequency = 0;
    }
}KmerHashTable;
#pragma pack ()

#pragma pack(2)
typedef struct KmerSetHashTable{
    unsigned long long int kmer;
    long int * frequency;
    KmerSetHashTable(){
        kmer = 0;
        frequency = NULL;
    }
}KmerSetHashTable;
#pragma pack ()

#pragma pack(2)
typedef struct KmerSetHashTableHead{
    KmerSetHashTable * kmerSetHashTable;
    unsigned long long int kmerSetHashTableCount;
    double * avgKmerFrequency;
    long int setNumber;
    long int kmerLength;
    KmerSetHashTableHead(){
        kmerSetHashTable = NULL;
        kmerSetHashTableCount = 0;
        avgKmerFrequency = NULL;
        setNumber = 0;
        kmerLength = 0;
    }
}KmerSetHashTableHead;
#pragma pack ()

#pragma pack(2)
typedef struct KmerSet{
    KmerHashTable * kmerHashTable;
    unsigned long long int kmerNumber;
    unsigned int kmerLength;
    unsigned long long int kmerHashTableCount;
    char * address;
    int addressNumber;
    KmerSet(){
        kmerHashTable = NULL;
        kmerNumber = 0;
        kmerLength = 0;
        kmerHashTableCount = 0;
        address = NULL;
        addressNumber = 0;
    }
}KmerSet;
#pragma pack ()


//--------------------------------

void SetBit(unsigned long long int * bit_array, unsigned int bit_number, char value);

void SetBit(unsigned long long int * bit_array, unsigned int bit_number, int len, char * value);

char GetBit(unsigned long long int * bit_array, unsigned int bit_number);

void GetBit(unsigned long long int * bit_array, unsigned int bit_number, int len, char * value);

void SetBit(char bit_array[], unsigned int bit_number, char value);

char GetBit(char bit_array[], unsigned int bit_number);

void SetBit(char bit_array[], unsigned int bit_start_number, unsigned int bit_length, char * value);

void SetBit(char bit_array[], unsigned int bit_start_number, unsigned int bit_length, char * value, int index);

char * GetBit(char bit_array[], unsigned int bit_start_number, unsigned int bit_length);

void GetBit(char bit_array[], unsigned int bit_start_number, unsigned int bit_length, char * value);

//-----------------------------------------

void GetKmerSetFromBcalm(const char * bcalmFileName,long int kmerLength);

unsigned long long int KmerHashTableSearch(char * value, KmerHashTable * kmerHashTable, unsigned int kmerLength, unsigned long int kmerHashTableCount);

unsigned long long int KmerHashTableSearch(char * value, unsigned int kmerLength);

void KmerHashTableInsert(char * kmer, long int frequency, unsigned int kmerLength, KmerHashTable * kmerHashTable, unsigned long int kmerHashTableCount, unsigned long long int * kmerNumber);

//-------------------------------------------

void create_hash_function_from_m_mers(int m);
void count_m_mers(string str, int m, int k);
void init_m_mers_table(int m);

typedef unordered_map<string,int> HashMap;

HashMap build_hash_map(int len);

int shash(const string& s, int& previous_shash, unsigned int start_pos = 0, int length = -1);

string inverse_shash (int num, int len);

int minimiserrc(const string &node,const int &minimisersize);

int minimiserv(const string &node,const int &minimisersize);

int minbutbiggerthan(int leftmin, int rightmin, const string &namebucket);

string reversecompletment(const string& str );

bool adjacent (const string& node1,const  string& node2,int k);

string readn(ifstream *file,uint64_t n);

int chartoint(char c);

string minimalsub(const string &w, const int &p,const int &k);

string minimalsub2(const string &w, const int &p,const int &k);

unsigned long int HashKmer(char * str, unsigned int len, unsigned long int max);





class neighbour
{
	public:
		array<pair<uint64_t,unsigned char>,8> list;
		//~ neighbour()
		//~ {
			//~ for(int )
			//~ list[i]=make_pair(0,0);
		//~ }
		uint64_t nbtype(unsigned char c);
		uint64_t gtype(unsigned char c);

		void add(uint64_t p,unsigned char b);
		unsigned char remove(uint64_t v);
		unsigned char removep(uint64_t v,unsigned char c);
		unsigned char removetype(unsigned char c);

};

class graph
{
	public:
		uint64_t n;
		int k;
		vector<string> nodes;
		vector<int> leftmins;
		vector<int> rightmins;
		unordered_multimap<uint64_t,uint64_t> map;
		//no more consider maprev
		unordered_multimap<uint64_t,uint64_t> maprev;
		vector<neighbour> neighbor;

		graph(const int ni)
		{
			k=ni;
			n=1;
			nodes.push_back("");
			leftmins.push_back(-1);
			rightmins.push_back(-1);
		}

		uint64_t getkey(string str);
		uint64_t getkeyrevc(string str);
		uint64_t becompacted(uint64_t nodeindice, int min, unsigned char *);
		int weight();
		void addvertex(const string str);
        void addleftmin(int mini);
        void addrightmin(int mini);
		void debruijn();
		void compressh(int min=-1);
		void compress();
		void importg(const char *name);
		void print(const char *name);
		void printedges(const char *name);
		void compact(uint64_t nodeindice,uint64_t with, unsigned char type);
		void reverse(int64_t with);
		void look(const uint64_t nodeindice, const string min);

};

#endif
