#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

//#include "pretest.h"

#include <pthread.h>
#include <time.h>
#include "bitarray.h"
#include "readset.h"
#include "kmerset.h"
#include "graph.h"
#include "contigmerge.h"
#include "constructcontigset.h"
#include "scaffolding.h"
#include "fillgap.h"


using namespace std;
 
int main(int argc, char *argv[])
{   
    
    int i = 0;
    int j = 0;
    long int setNumber = (argc-6)/5;
    long int kmerLength = atoi(argv[argc-2]);
    long int threadNumber = atoi(argv[argc-1]);
    unsigned long long int kmerSetHashTableCount = atoi(argv[argc-3]);
    char * kmerHashAddress = argv[argc-4];
    ReadSet * readSet = new ReadSet[setNumber];
	double * percentOfVar = new double[setNumber];
    
    for(i = 0; i<setNumber;i++){  
        readSet[i].leftAddress = argv[i*5+1];
        readSet[i].rightAddress = argv[i*5+2];
        readSet[i].insertSize = atoi(argv[i*5+3]);
		percentOfVar[i] = atof(argv[i*5+4]);
        readSet[i].var = percentOfVar[i];
        readSet[i].orientation = atoi(argv[i*5+5]);		
    }
        
	InitReadSet(readSet, setNumber, threadNumber);

 
	DBGraphHead * deBruijnGraphHead = new DBGraphHead;
    GetDBGraphFromAddress(deBruijnGraphHead,argv[argc-5]);
 
	KmerSetHashTableHead * kmerSetHashTableHead = new KmerSetHashTableHead;
	kmerSetHashTableHead->kmerLength = kmerLength;
	kmerSetHashTableHead->kmerSetHashTableCount = kmerSetHashTableCount;
	kmerSetHashTableHead->kmerSetHashTable = new KmerSetHashTable[kmerSetHashTableHead->kmerSetHashTableCount];
	GetKmerSetHashTableFromAddress(kmerSetHashTableHead->kmerSetHashTable, setNumber, kmerHashAddress);

    
    for(i = 0; i<setNumber; i++){
        GetAvgKmerNumberAndGapProblity(deBruijnGraphHead->deBruijnGraph, readSet, i, kmerSetHashTableHead->kmerSetHashTable, kmerSetHashTableHead->kmerSetHashTableCount, kmerLength);
    }

    
    ContigSet * contigSetHead = new ContigSet;
    ConstructContigSet(deBruijnGraphHead, readSet, contigSetHead, setNumber, kmerSetHashTableHead->kmerSetHashTable, kmerSetHashTableHead->kmerSetHashTableCount, kmerLength, threadNumber, extendCutOff);
    
    
    ContigSet * temp = contigSetHead->next;
    temp = contigSetHead->next;

    contigSetHead->next = ContigMerging(contigSetHead->next, deBruijnGraphHead->deBruijnGraph, kmerLength, readSet, setNumber, str1,1);
    
    temp = contigSetHead->next;
    
    contigSetHead->next = ContigMerging(contigSetHead->next, deBruijnGraphHead->deBruijnGraph, kmerLength, readSet, setNumber, str2,2);
    
    char * str = new char[20];
    char * str1 = new char[20];
    strcpy(str,"contigSet.fa");
    strcpy(str1,"contigSetLong.fa");
    
    temp = contigSetHead->next;
    WriteContigSet(temp, str);
    WriteContigSetLong(temp, str1);
    
    
    ScaffoldSetHead * scaffoldSetHead = ScaffoldingContigSet(temp, readSet, setNumber, 3*kmerLength, kmerLength, threadNumber);
    
    FillGap(scaffoldSetHead, readSet, setNumber, kmerSetHashTableHead, kmerLength, threadNumber);
    
    char * address0 = new char[30];
    strcpy(address0, "scaffold0.fa");
    WriteScaffoldSet(scaffoldSetHead, address0);
    contigSetHead = GetContigSetFromScaffoldSetHead(scaffoldSetHead);
    DeleteScaffoldSetHead(scaffoldSetHead);
    
    remove(address0);
    
    temp = contigSetHead->next;
    
    scaffoldSetHead = ScaffoldingContigSet(temp, readSet, setNumber, kmerLength, kmerLength, threadNumber);    
    
    FillGap(scaffoldSetHead, readSet, setNumber, kmerSetHashTableHead, kmerLength, threadNumber); 
    
    
}
