#ifndef KMERSET_HEAD
#define KMERSET_HEAD

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "bitarray.h"
#include <cstring>



using namespace std;


unsigned long int HashKmer(char * str, unsigned int len, unsigned long int max)  
{  
   unsigned int hash = 0;  
   unsigned int i = 0;  
  
   for(i = 0; i < len; str++, i++) {  
      hash = (*str) + (hash << 6) + (hash << 16) - hash;  
   }  
  
   return hash % max;  
} 


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

static KmerHashTable * kmerHashTableStatic = NULL;
static long int kmerHashTableCountStatic = 0;

unsigned long long int KmerHashTableSearch(char * value, KmerHashTable * kmerHashTable, unsigned int kmerLength, unsigned long int kmerHashTableCount){
    
    unsigned long int i = 0;
    unsigned long int j = 0;
    int p = 0;
    i = HashKmer(value, kmerLength, kmerHashTableCount);
    unsigned long long int tempValue = 0;
    SetBit(&tempValue, 0, kmerLength, value);
    tempValue++;
    
    while(kmerHashTable[i].kmer!=0){   
        if(kmerHashTable[i].kmer == tempValue){
            return (i + 1);
        }
        i = (i+1)%kmerHashTableCount; 
    }
       
    return 0;
      
}

//overload
unsigned long long int KmerHashTableSearch(char * value, unsigned int kmerLength){
    
    unsigned long int i = 0;
    unsigned long int j = 0;
    int p = 0;
    i = HashKmer(value, kmerLength, kmerHashTableCountStatic);
    unsigned long long int tempValue = 0;
    SetBit(&tempValue, 0, kmerLength, value);
    tempValue++;
    
    while(kmerHashTableStatic[i].kmer!=0){   
        if(kmerHashTableStatic[i].kmer == tempValue){
            return (i + 1);
        }
        i = (i+1)%kmerHashTableCountStatic; 
    }
       
    return 0;
      
}


void KmerHashTableInsert(char * kmer, long int frequency, unsigned int kmerLength, KmerHashTable * kmerHashTable, unsigned long int kmerHashTableCount, unsigned long long int * kmerNumber){
     
    int j = 0;

    unsigned long int i = HashKmer(kmer, kmerLength, kmerHashTableCount);
    long int previousIndex = i;
    
    unsigned long long int tempKmer = 0;
    SetBit(&tempKmer, 0, kmerLength, kmer);
    tempKmer++;
    unsigned long long int cur = 0;
    do{
        previousIndex = i;
        i = (i + 1) % (kmerHashTableCount);
        cur = __sync_val_compare_and_swap(&kmerHashTable[previousIndex].kmer, 0, tempKmer);
        
    }while(cur!=0&&cur!=tempKmer);
    
    if(cur == 0){
        __sync_fetch_and_add(kmerNumber, 1);
    }
    
    __sync_fetch_and_add(&kmerHashTable[previousIndex].frequency, frequency);
      
}



void GetKmerSetFromBcalm(const char * bcalmFileName,long int kmerLength){
    
    long int i = 0;
    long int j = 0;
    
    //long int kmerLength = 22;
    unsigned long long int kmerNumber = 0;
    kmerHashTableCountStatic = 250000000;
    kmerHashTableStatic = new KmerHashTable[kmerHashTableCountStatic];
    
    char * kmer = new char[kmerLength+1];
    kmer[kmerLength] = '\0';
    
    ifstream icin;
    icin.open(bcalmFileName);
    
    char * contig = new char[150000];
    long int contigLength = 0;
    while(icin.getline(contig, 150000)){
        
        contigLength = strlen(contig);
        for(i = 0; i< contigLength-1-kmerLength; i++){
            strncpy(kmer, contig+i, kmerLength);
            for(j=0;j<kmerLength;j++){
                kmer[j]=toupper(kmer[j]);
            }
            //cout<<kmer<<endl;
            KmerHashTableInsert(kmer, 1, kmerLength, kmerHashTableStatic, kmerHashTableCountStatic, & kmerNumber);
        }    
    }
    
    /*
    ofstream ocout;
    char * outPut = new char[100];
    strcpy(outPut, "bcalmKmer.fa");
    ocout.open(outPut);
    
    char * temp11 = new char[kmerLength + 1];
    unsigned long long int temp = 0;
    
    for(i = 0; i < kmerHashTableCount; i++){
        if(kmerHashTable[i].frequency>1){
            cout<<temp<<"--"<<temp11<<","<<kmerHashTable[i].frequency<<endl;
            exit(0);
        }
        
        if(kmerHashTable[i].kmer!=0){                 
            temp = kmerHashTable[i].kmer - 1;               
            GetBit(&temp,0,kmerLength, temp11);
            ocout<<temp<<"--"<<temp11<<","<<kmerHashTable[i].frequency<<endl;

        }
        
    }
    ocout.close();
    */
    
}


#endif
