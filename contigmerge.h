#ifndef CONTIGMERGE_H_INCLUDED 
#define CONTIGMERGE_H_INCLUDED 

#include <cstring>
#include <iostream>
#include "stdio.h"
#include "stdlib.h"
#include <math.h>

#include "readset.h"
#include "constructcontigset.h"
#include "contigoverlap.h"

ContigSet * GetContigFromContigSet(ContigSet * head, long int index){
    long int i = 0;
    while(head!=NULL){
        if(i == index){
            return head;
        }
        i++;
        head = head->next;
    }
    return NULL;
}

typedef struct GetRemoveSubContigsP{
    ContigSet * contigSet;
    int threadIndex;
    int totalThreadNumber;
    bool * index;
    bool * singleIndex;
    long int contigNumber;
    char * outFileName;
}GetRemoveSubContigsP;

pthread_mutex_t mutexRemoveSubContigs = PTHREAD_MUTEX_INITIALIZER;

void * RemoveSubContigsThread(void * arg){
    GetRemoveSubContigsP * getRemoveSubContigsP = (GetRemoveSubContigsP *)arg;
    unsigned long int i = 0;
    unsigned long int p = 0;
    int j = 0;
	
	ContigSet * contigSet = getRemoveSubContigsP->contigSet;
	p = getRemoveSubContigsP->threadIndex;
	int totalThreadNumber = getRemoveSubContigsP->totalThreadNumber;
	bool * index = getRemoveSubContigsP->index;
	bool * singleIndex = getRemoveSubContigsP->singleIndex;
	long int contigNumber = getRemoveSubContigsP->contigNumber;
	char * outFileName = getRemoveSubContigsP->outFileName;
	ofstream ocout;
    ocout.open(outFileName,ios::app);
	while(contigSet!=NULL){
     
        if(i%totalThreadNumber==p){
            ContigSet * tempContigSet = getRemoveSubContigsP->contigSet;
            long int contigLength = strlen(contigSet->contig);
            cout<<contigLength +1 -42<<endl;
            char * tempContig = new char[contigLength +1 -42];
            strncpy(tempContig, contigSet->contig+21, contigLength -42);
            tempContig[contigLength -42] = '\0';
            char * tempContigRS = new char[contigLength+1 -42];
            ReverseComplement(tempContig,tempContigRS);
            j = 0;
            //cout<<"aa--"<<endl;
            while(tempContigSet!=NULL){
                if(j!=i){
                    char * temp = strstr(tempContigSet->contig,tempContig);
                    if(temp!=NULL){
                        pthread_mutex_lock(&mutexRemoveSubContigs);
                        ocout<<j<<","<<i<<","<<"1"<<endl;
                        pthread_mutex_unlock(&mutexRemoveSubContigs);
                    }
                    temp = strstr(tempContigSet->contig,tempContigRS);
                    if(temp!=NULL){
                        pthread_mutex_lock(&mutexRemoveSubContigs);
                        ocout<<j<<","<<i<<","<<"0"<<endl;
                        pthread_mutex_unlock(&mutexRemoveSubContigs);
                    }
                }
                tempContigSet=tempContigSet->next;
                j++;
            }
            //cout<<"cc"<<endl;
            delete [] tempContigRS;
            tempContigRS = NULL;
            delete [] tempContig;
            tempContig = NULL;
            //cout<<"bb"<<endl;
        }
        
        i++; 
        contigSet = contigSet->next;
    }
	ocout.close();
	
}

void RemoveSubContigs(ContigSet * contigSet, int totalThreadNumber){
    
    pthread_t tid[totalThreadNumber];
    
    ContigSet * tempContigSet = contigSet;
    long int contigNumber = 0;
    while(tempContigSet!=NULL){
        if(tempContigSet->contig!=NULL){
            contigNumber++;
        }
        tempContigSet = tempContigSet->next;
    }
    long int allIndexNumber = contigNumber*contigNumber;
    bool * index = new bool[allIndexNumber];
    bool * singleIndex = new bool[allIndexNumber];
    
    for(long int p = 0; p<allIndexNumber; p++){
        index[p] = false;
    }
    for(long int p = 0; p<contigNumber; p++){
        singleIndex[p] = false;
    }
    
    char * outFileName = new char[30];
	strcpy(outFileName,"SubContig.fa");
    
    remove(outFileName);
    
    int i = 0;
    GetRemoveSubContigsP * getRemoveSubContigsP = new GetRemoveSubContigsP[totalThreadNumber];
    //cout<<"fisrt-------"<<endl;
    for(i = 0; i< totalThreadNumber; i++){
        getRemoveSubContigsP[i].contigSet = contigSet;
        getRemoveSubContigsP[i].threadIndex = i;
        getRemoveSubContigsP[i].totalThreadNumber = totalThreadNumber;
        getRemoveSubContigsP[i].index = index;
        getRemoveSubContigsP[i].singleIndex = singleIndex;
        getRemoveSubContigsP[i].contigNumber = contigNumber;
        getRemoveSubContigsP[i].outFileName = outFileName;
             
        if(pthread_create(&tid[i], NULL, RemoveSubContigsThread, (void *)&getRemoveSubContigsP[i])){
            cout<<"create thread wrong!"<<endl;
            exit(1);
        }      
    }
    //cout<<"second-------"<<endl;
    for(i = 0; i < totalThreadNumber; i++){
        pthread_join(tid[i], NULL);
    }
    //cout<<"third-------"<<endl;
    
    
	ifstream icin;
	icin.open(outFileName);
	char * overlap = new char[100];
	long int contigOverlap[1][2];
	char * tempNumber = new char[10];
	while(icin.getline(overlap,100)){
        int length = strlen(overlap);
        int j = 0;
        int p = 0;
        for(i=0;i<length;i++){
            if(overlap[i]==','){
                tempNumber[j]='\0';
                j=0;
                contigOverlap[0][p] = atoi(tempNumber);
                p++;
                if(p==2){
                    break;
                }
                continue;
            }
            tempNumber[j] = overlap[i];
            j++;
        }
        if(index[contigOverlap[0][0]*contigNumber+contigOverlap[0][1]]!=true 
            && index[contigOverlap[0][1]*contigNumber+contigOverlap[0][0]]!=true){
            //cout<<contigOverlap[0][0]<<","<<contigOverlap[0][1]<<endl;
            ContigSet * tempContigSet = GetContigFromContigSet(contigSet,contigOverlap[0][1]);
            delete [] tempContigSet->contig;
            tempContigSet->contig = NULL;
        }
        index[contigOverlap[0][0]*contigNumber+contigOverlap[0][1]] = true;
        index[contigOverlap[0][1]*contigNumber+contigOverlap[0][0]] = true;
    }
    icin.close();

}


typedef struct GetCompactContigsP{
    ContigSet * contigSet;
    int threadIndex;
    int totalThreadNumber;
    long int contigNumber;
    long int kmerLength;
    bool * compactIndex;
}GetCompactContigsP;

pthread_mutex_t mutexCompactContigs = PTHREAD_MUTEX_INITIALIZER;

void * CompactContigsThread(void * arg){
    GetCompactContigsP * getCompactContigsP = (GetCompactContigsP *)arg;
    unsigned long int i = 0;
    unsigned long int thread = 0;
    int j = 0;
	
	ContigSet * contigSet = getCompactContigsP->contigSet;
	thread = getCompactContigsP->threadIndex;
	int totalThreadNumber = getCompactContigsP->totalThreadNumber;
	long int contigNumber = getCompactContigsP->contigNumber;
	long int kmerLength = getCompactContigsP->kmerLength;
	bool * compactIndex = getCompactContigsP->compactIndex;
	
	long int misMatchCutOff = 3;
	
	char * outFileName = new char[30];
	strcpy(outFileName,"CompactContig.fa");
	ofstream ocout;
    ocout.open(outFileName,ios::app);
	while(contigSet!=NULL){
     
        if(i%totalThreadNumber==thread){
            ContigSet * rightContigSet = getCompactContigsP->contigSet;
            char * leftContig = contigSet->contig;
            long int leftLength = strlen(leftContig);
            j = 0;
            while(rightContigSet!=NULL){
                char * rightContig = rightContigSet->contig;
                long int rightLength = strlen(rightContig);
                if(j!=i && compactIndex[i*contigNumber+j] != true && compactIndex[j*contigNumber+i] != true){
                    long int t = 0;
                    long int p = 0;
                    long int index = 0;
                    long int misMatch = 0;
                    for(t=0;t<leftLength;t++){
                        index = 0;
                        misMatch = 0;
                        for(p=0;p<rightLength&&t+p<leftLength;p++){
                            if(leftContig[t+p]!=rightContig[p]){
                                misMatch++;
                                if(misMatch > misMatchCutOff){
                                    index = 1;
                                    break;
                                } 
                            }
                        }
                        if(index==0){
                            break;
                        }
                    }
                    pthread_mutex_lock(&mutexCompactContigs);
                    if((t+p)==leftLength&&p>kmerLength && compactIndex[i*contigNumber+j] != true){                       
                        ocout<<i<<"	"<<j<<"	"<<p<<"	+	+"<<endl;      
                        compactIndex[i*contigNumber+j] = true;
                        compactIndex[j*contigNumber+i] = true;
                        pthread_mutex_unlock(&mutexCompactContigs);
                    }
                    pthread_mutex_unlock(&mutexCompactContigs);
                    index = 0;
                    char * leftContigRS = new char[leftLength+1];
                    ReverseComplement(leftContig,leftContigRS);
                    for(t=0;t<leftLength;t++){
                        index = 0;
                        misMatch = 0;
                        for(p=0;p<rightLength&&t+p<leftLength;p++){
                            if(leftContigRS[t+p]!=rightContig[p]){
                                misMatch++;
                                if(misMatch > misMatchCutOff){
                                    index = 1;
                                    break;
                                } 
                            }
                        }
                        if(index==0){
                            break;
                        }
                    }
                    delete [] leftContigRS;
                    leftContigRS = NULL;
                    pthread_mutex_lock(&mutexCompactContigs);
                    if((t+p)==leftLength&&p>kmerLength && compactIndex[i*contigNumber+j] != true){                     
                        ocout<<i<<"	"<<j<<"	"<<p<<"	-	+"<<endl;
                        compactIndex[i*contigNumber+j] = true;
                        compactIndex[j*contigNumber+i] = true;
                        pthread_mutex_unlock(&mutexCompactContigs);
                    }
                    pthread_mutex_unlock(&mutexCompactContigs);
                    
                    index = 0;
                    char * rightContigRS = new char[rightLength+1];
                    ReverseComplement(rightContig, rightContigRS);
                    for(t=0;t<leftLength;t++){
                        index = 0;
                        misMatch = 0;
                        for(p=0;p<rightLength&&t+p<leftLength;p++){
                            if(leftContig[t+p]!=rightContigRS[p]){
                                misMatch++;
                                if(misMatch > misMatchCutOff){
                                    index = 1;
                                    break;
                                } 
                            }
                        }
                        if(index==0){
                            break;
                        }
                    }
                    delete [] rightContigRS;
                    rightContigRS = NULL;
                    pthread_mutex_lock(&mutexCompactContigs);
                    if((t+p)==leftLength&&p>kmerLength && compactIndex[i*contigNumber+j] != true){                     
                        ocout<<i<<"	"<<j<<"	"<<p<<"	+	-"<<endl;
                        compactIndex[i*contigNumber+j] = true;
                        compactIndex[j*contigNumber+i] = true;
                        pthread_mutex_unlock(&mutexCompactContigs);
                    }
                    compactIndex[i*contigNumber+j] = true;
                    compactIndex[j*contigNumber+i] = true;
                    pthread_mutex_unlock(&mutexCompactContigs);
                }
                rightContigSet=rightContigSet->next;
                j++;
            }

        }
        
        i++; 
        contigSet = contigSet->next;
    }
	ocout.close();
	
}

void CompactContigs(ContigSet * contigSet, long int kmerLength, int totalThreadNumber){
    
    pthread_t tid[totalThreadNumber];
    
    ContigSet * tempContigSet = contigSet;
    long int contigNumber = 0;
    while(tempContigSet!=NULL){
        if(tempContigSet->contig!=NULL){
            contigNumber++;
        }
        tempContigSet = tempContigSet->next;
    }
    
    long int allIndexNumber = contigNumber*contigNumber;
    bool * index = new bool[allIndexNumber];   
    for(long int p = 0; p<allIndexNumber; p++){
        index[p] = false;
    }
    
    int i = 0;
    GetCompactContigsP * getCompactContigsP = new GetCompactContigsP[totalThreadNumber];
    //cout<<"fisrt-------"<<endl;
    for(i = 0; i< totalThreadNumber; i++){
        getCompactContigsP[i].contigSet = contigSet;
        getCompactContigsP[i].threadIndex = i;
        getCompactContigsP[i].totalThreadNumber = totalThreadNumber;
        getCompactContigsP[i].contigNumber = contigNumber;
        getCompactContigsP[i].kmerLength = kmerLength;
        getCompactContigsP[i].compactIndex = index;
             
        if(pthread_create(&tid[i], NULL, CompactContigsThread, (void *)&getCompactContigsP[i])){
            cout<<"create thread wrong!"<<endl;
            exit(1);
        }      
    }
    //cout<<"second-------"<<endl;
    for(i = 0; i < totalThreadNumber; i++){
        pthread_join(tid[i], NULL);
    }
    //cout<<"third-------"<<endl;
}



double WeightContigMerge(char * contigLeft, char * contigRight, ReadSet * readSet, long int readSetIndex,long int kmerLength, long int overlap, long int length){
    
    long int i = 0;
    long int j = 0;
    long int readLength = readSet[readSetIndex].readLength;
    long int insertSize = readSet[readSetIndex].insertSize;
    long int var = readSet[readSetIndex].var;
    
    long int nullNumber = 0;
    long int positionMatchNumber = 0;
    long int positionNonExist = 0;
    
    long int min = 0;
    long int index1 = 0;

    char * tempContigLeft = contigLeft;
    long int tempLengthLeft = strlen(contigLeft);
    if(tempLengthLeft>=insertSize + lambda*var){
        min = insertSize + lambda*var;
        tempContigLeft = new char[min+1];
        SubContig(tempContigLeft,contigLeft,tempLengthLeft-min,tempLengthLeft);
        tempLengthLeft = min;
    }

    char * tempContigRight = contigRight;
    long int tempLengthRight = strlen(contigRight);
    if(tempLengthRight>=insertSize + lambda*var){
        min = insertSize + lambda*var;
        tempContigRight = new char[min+1];
        SubContig(tempContigRight,contigRight,0,min);
        tempLengthRight = min;
    }

    char * temp = new char[readLength + 1];
    char * tempRC = new char[readLength + 1];
    
    ReadMate * temp1;
    
    long int * gapIndex = new long int[tempLengthRight];
    long int * nullIndex = new long int[tempLengthRight];
    long int * matchIndex = new long int[tempLengthRight];
    
    for(j=overlap-kmerLength;j<tempLengthRight - readLength - length && j<insertSize - lambda*var - readLength && tempLengthLeft > insertSize;j++){
        index1 = 0;
        matchIndex[j] = 0;
        
        SubContig(temp, tempContigRight, j, j + readLength);
        
        ReadMate * tempReadMate = SearchLeftReadMate(temp, readSet+readSetIndex);
        long int a = SearchRead(temp, readSet + readSetIndex);
        
        if(tempReadMate==NULL && a ==-1){
            positionNonExist++;
            gapIndex[j] = 1;
            index1 = 1;
        }else{
            gapIndex[j] = 0;
        }
        
        if(tempReadMate!=NULL){
            temp1 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete temp1;
        }   
        
        if(tempReadMate==NULL&&index1!=1){
            nullIndex[j] = 1;
            nullNumber++;
        }else{
            nullIndex[j] = 0;
        }
        
        bool token = false;

        while(tempReadMate!=NULL){
            index1 = KMPIndexOfContigOfMisMatch(tempContigLeft,tempReadMate->readMate);
            long int tempDD = j + tempLengthLeft - index1 + readLength - overlap;
            if(index1!=-1 && tempDD<insertSize+lambda*var && tempDD>insertSize-lambda*var){
                if(token!=true){
                    positionMatchNumber++;
                    matchIndex[j] = 1;
                    token = true;
                    break;
                }
            }
            temp1 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete [] temp1->readMate;
            delete temp1;       
        }
        if(token!=true){
            matchIndex[j] = 0;
        } 
        
        i++;
        if(i >= length){
            if(positionMatchNumber<=1){
                
                if(tempLengthLeft == insertSize + lambda*var){
                    delete [] tempContigLeft;
                }
                if(tempLengthRight == insertSize + lambda*var){
                    delete [] tempContigRight;
                }
                delete [] temp;
                delete [] tempRC;
                delete [] gapIndex;
                delete [] nullIndex;
                delete [] matchIndex;
                
                return 0;
            }
            double a = 0;
            if(length - nullNumber - positionNonExist>1){
                a = (double)(positionMatchNumber)/(double)(length - nullNumber - positionNonExist);
            }else{
                a = 0.2;
            }
            
            if(a<0.2){
                
                if(tempLengthLeft == insertSize + lambda*var){
                    delete [] tempContigLeft;
                }
                if(tempLengthRight == insertSize + lambda*var){
                    delete [] tempContigRight;
                }
                delete [] temp;
                delete [] tempRC;
                delete [] gapIndex;
                delete [] nullIndex;
                delete [] matchIndex;
                
                return 0;
            }
            
            nullNumber = 0;
            positionMatchNumber = 0;
            positionNonExist = 0;
            i = 0;
        }
                       
    }

    i = 0;
    nullNumber = 0;
    positionMatchNumber = 0;
    positionNonExist = 0;
    
    for(j=overlap - kmerLength;j<tempLengthLeft - readLength - length && j<insertSize - lambda*var - readLength && tempLengthRight>insertSize;j++){
        index1 = 0;
        matchIndex[j] = 0;
        SubContig(temp, tempContigLeft, tempLengthLeft - j - readLength, tempLengthLeft - j);

        ReadMate * tempReadMate = SearchRightReadMate(temp, readSet+readSetIndex);
        ReverseComplement(temp, tempRC);
        long int a = SearchRead(tempRC, readSet + readSetIndex);
        
        if(tempReadMate==NULL && a ==-1){
            positionNonExist++;
            gapIndex[j] = 1;
            index1 = 1;
        }else{
            gapIndex[j] = 0;
        }
        
        if(tempReadMate!=NULL){
            temp1 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete temp1;
        }   
        
        if(tempReadMate==NULL&&index1!=1){
            nullNumber++;
            nullIndex[j] = 1;
        }else{
            nullIndex[j] = 0;
        }
        
        bool token = false;

        while(tempReadMate!=NULL){
            index1 = KMPIndexOfContigOfMisMatch(tempContigRight,tempReadMate->readMate);
            long int tempDD = index1 + 2*readLength + j - overlap;
            if(index1!=-1 && tempDD<insertSize+lambda*var && tempDD>insertSize-lambda*var){
                if(token!=true){
                    positionMatchNumber++;
                    matchIndex[j] = 1;
                    token = true;
                    break;
                }
            }
            temp1 = tempReadMate;
            tempReadMate = tempReadMate->next;
            delete [] temp1->readMate;
            delete temp1;       
        } 
        
        if(token!=true){
            matchIndex[j] = 0;
        }
        
        i++;
        if(i >= length){
            
            if(positionMatchNumber<=1){
            
                
                if(tempLengthLeft == insertSize + lambda*var){
                    delete [] tempContigLeft;
                }
                if(tempLengthRight == insertSize + lambda*var){
                    delete [] tempContigRight;
                }
                delete [] temp;
                delete [] tempRC;
                delete [] gapIndex;
                delete [] nullIndex;
                delete [] matchIndex;
                
                
                return 0;
            }
            
            double a = 0;
            if(length - nullNumber - positionNonExist>1){
                a = (double)(positionMatchNumber)/(double)(length - nullNumber - positionNonExist);
            }else{
                a = 0.2;
            }
            
            if(a<0.2){
                
                
                if(tempLengthLeft == insertSize + lambda*var){
                    delete [] tempContigLeft;
                }
                if(tempLengthRight == insertSize + lambda*var){
                    delete [] tempContigRight;
                }
                delete [] temp;
                delete [] tempRC;
                delete [] gapIndex;
                delete [] nullIndex;
                delete [] matchIndex;
                
                return 0;
            }
            
            
            nullNumber = 0;
            positionMatchNumber = 0;
            positionNonExist = 0;
            i = 0;
        }
                       
    }


    delete [] temp;
    delete [] tempRC;
    delete [] gapIndex;
    delete [] matchIndex;
    delete [] nullIndex;
    if(tempLengthRight==insertSize+ lambda*var){
        delete [] tempContigRight;
    }
    if(tempLengthLeft==insertSize+ lambda*var){
        delete [] tempContigLeft;
    }
      
    return 1;
    
}


int MergeSubContig(ContigSet * contigSet, long int kmerLength){
    long int i = 0;
    long int j = 0;
    int count = 0;
    ContigSet * temp;
    
    ContigSet * first = contigSet;

    while(contigSet!=NULL){
        if(contigSet->contig==NULL){
            contigSet = contigSet->next;
            i++;
            continue;
        }
        temp = first;
        j = 0;
        while(temp!=NULL){
            if(temp->contig==NULL || i==j){
                j++;
                temp = temp->next;
                continue;
            }
            int len = strlen(temp->contig);
            long int * next = new long int[len+1-kmerLength];
            char * pattern = new char[len+1-kmerLength];
            SubContig(pattern, temp->contig, kmerLength, len);
            long int p = KMPIndexOfContig(contigSet->contig,pattern,next);
            if(p!=-1){
                
                delete []temp->contig;
                temp->contig = NULL;
                delete []next;
                next = NULL;
                delete [] pattern;
                pattern = NULL;
                j++;
                temp = temp->next;
                continue;
            }
            
            SubContig(pattern, temp->contig, 0, len-kmerLength);
            p = KMPIndexOfContig(contigSet->contig,pattern,next);
            if(p!=-1){
                
                delete []temp->contig;
                temp->contig = NULL;
                delete []next;
                next = NULL;
                delete [] pattern;
                pattern = NULL;
                j++;
                temp = temp->next;
                continue;
            }
                       
      
            char * tempReverseContig = new char[len+1];
            ReverseComplement(temp->contig,tempReverseContig);
            SubContig(pattern, tempReverseContig, kmerLength, len);
            p = KMPIndexOfContig(contigSet->contig,pattern,next);
            if(p!=-1){
                
                delete []temp->contig;
                temp->contig = NULL;
                delete []next;
                next = NULL;
                delete []tempReverseContig;
                tempReverseContig = NULL;
                delete [] pattern;
                pattern = NULL;
                j++;
                temp = temp->next;
                continue;
            }
            
            SubContig(pattern, tempReverseContig, 0, len-kmerLength);
            p = KMPIndexOfContig(contigSet->contig,pattern,next);
            if(p!=-1){
                
                delete []temp->contig;
                temp->contig = NULL;
                
            }
            
            
            delete []tempReverseContig;
            delete []next;
            next = NULL; 
            delete [] pattern;
            pattern = NULL;
            temp = temp->next; 
            j++;          
        }
        i++;
        contigSet = contigSet->next;
    }
    return 1;
}


char * LengthOfOverlapBetweenContig(char * left, char * right, ReadSet * readSet, int readSetCount, int kmerLength){
    long int i = 0;
    long int j = 0;
    long int t = 0;
    long int index = 0;
    long int leftLength = (long int)strlen(left);
    long int rightLength = (long int)strlen(right);
    
    
    for(i=0;i<leftLength;i++){
        index = 0;
        for(j=0;j<rightLength&&i+j<leftLength;j++){
            if(left[i+j]!=right[j]){
                index = 1;
                break;
            }
        }
        if(index==0){
            break;
        }
    }
    
    if((i+j)==leftLength&&j>kmerLength){ 
        for(int p = 0;p<readSetCount;p++){
            
            if(readSet[p].insertSize - lambda*readSet[p].var<j||leftLength<readSet[p].insertSize- 2*readSet[p].readLength||rightLength<readSet[p].insertSize- 2*readSet[p].readLength){
                continue;
            }
            double temp = WeightContigMerge(left, right, readSet, p, kmerLength, j, kmerLength);
            
            if(temp!=1){
                return NULL;
            }    
        } 
        char * contig = new char[leftLength + rightLength - j + 1];
        AppendRight(contig,left,right,j+1);
        return contig;
    }
    return NULL;
}


typedef struct OverlapInformation{
    char * infor;
    OverlapInformation(){
        infor = new char[10];
    }
}OverlapInformation;

typedef struct ContigMerge{
    long int contigIndex;
    bool orientation;
    long int overlapLength;
    bool merge;
    ContigMerge * next;
    ContigMerge(){
        contigIndex = -1;
        orientation = 1;
        overlapLength = 0;
        next = NULL;
        merge = false;
    }
}ContigMerge;

typedef struct ContigMergeSet{
    ContigMerge * contigMergeSet;
    ContigMergeSet * next;
    ContigMergeSet(){
        contigMergeSet = NULL;
        next = NULL;
    }
}ContigMergeSet;
/*
void CompactContigSet(ContigSet * contigSet, char * file){
    
    long int i = 0;
    long int j = 0;
    long int p = 0;
    
    long int contigNumber = 0;
    ContigSet * tempContigSet = contigSet;
    while(tempContigSet !=NULL){
        contigNumber ++;
        tempContigSet = tempContigSet->next;
    }
    
    ifstream icin;
    icin.open(file);
    
    int * leftOrientation = new int[contigNumber];
    int * rightOrientation = new int[contigNumber];
    int * leftNumber = new int[contigNumber];
    int * rightNumber = new int[contigNumber];
    long int * leftIndex = new long int[contigNumber];
    long int * rightIndex = new long int[contigNumber];
    long int * overlapLength = new long int[contigNumber];
    bool * visited = new bool[contigNumber];
    for(i=0;i<contigNumber;i++){
        leftOrientation[i] = -1;
        rightOrientation[i] = -1;
        leftNumber[i] = 0;
        rightNumber[i] = 0;
        overlapLength[i] = 0;
        visited[i] = false;
    }

    char * overlapInfor = new char[1000];
    
    while(icin.getline(overlapInfor,1000)){
        long int start = 0;
        long int previous = 0;
        long int len = strlen(overlapInfor);
        OverlapInformation * overlap = new OverlapInformation[5];
        i = 0;
        j = 0;
        for(start = 0;start<len;start++){   
            if(overlapInfor[start]=='	'){
                overlap[i].infor[j] = '\0';
                i++;
                j = 0;
                continue;
            }
            overlap[i].infor[j] = overlapInfor[start];
            j++;
        }
        leftIndex[p] = atoi(overlap[0].infor);
        rightIndex[p] = atoi(overlap[1].infor);
        overlapLength[p] = atoi(overlap[2].infor);
        leftOrientation[p] = 1;
        if(overlap[3].infor[0]=='-'){
            leftOrientation[p] = 0;
        }
        rightOrientation[p] = 1;
        if(overlap[4].infor[0]=='-'){
            rightOrientation[p] = 0;
        }   
        p++;      
    }
    
    for(i=0;i<p;i++){
        
    }
    
    
}
*/
ContigSet * ContigMerging(ContigSet * contigSet, DBGraph * deBruijnGraph, int kmerLength, ReadSet * readSet, int readSetCount, char * contigSetFile, int token){
    long int i = 0;
    long int j = 0;
    int count = 0;
    long int p = 0;
    
    
    char * contigSetNonIncludeFile = new char[20];
    char * overlapInformationFile = new char[20];
    strcpy(contigSetNonIncludeFile,"nonInclude.fa");
    strcpy(overlapInformationFile,"CompactContig.fa");
    remove(overlapInformationFile);
    
    //cout<<"sssssssssssssssssssssssssss0000000000000000000000"<<endl;
    
    RemoveSubContigs(contigSet,16);
    WriteContigSetLong(contigSet, contigSetNonIncludeFile);
    contigSet = GetContigSetLong(contigSetNonIncludeFile);
    CompactContigs(contigSet,kmerLength,16);

    //cout<<"ffffffffffffffffffffffffff"<<endl;
    //mapmerge(kmerLength, contigSet, contigSetNonIncludeFile, overlapInformationFile, 1,token);
    
    //cout<<"sssssssssssssssssssssssssss111111111111111111111"<<endl;
    
    
    //contigSet = GetContigSetLong(contigSetNonIncludeFile);
    
    long int overlapNumber = 0;
    ifstream icin1;
    icin1.open(overlapInformationFile);
    char * line = new char[1000];
    while(icin1.getline(line,1000)){
        overlapNumber++;
    }
    icin1.close();
    
    
    
    long int contigNumber = 0;
    ContigSet * tempContigSet = contigSet;
    while(tempContigSet !=NULL){
        contigNumber ++;
        tempContigSet = tempContigSet->next;
    }
    
    //cout<<"overlapNumber:"<<overlapNumber<<endl;
    //cout<<"contigNumber:"<<contigNumber<<endl;
    
    int * leftOrientation = new int[overlapNumber];
    int * rightOrientation = new int[overlapNumber];
    int * leftNumber = new int[contigNumber];
    int * rightNumber = new int[contigNumber];
    long int * leftIndex = new long int[overlapNumber];
    long int * rightIndex = new long int[overlapNumber];
    long int * overlapLength = new long int[overlapNumber];
    long int * leftMaxOverlapLength = new long int[contigNumber];
    long int * rightMaxOverlapLength = new long int[contigNumber];
    long int * leftMaxIndex = new long int[contigNumber];
    long int * rightMaxIndex = new long int[contigNumber];
    bool * visited = new bool[contigNumber];
    for(i=0;i<contigNumber;i++){
        leftNumber[i] = 0;
        rightNumber[i] = 0;
        leftMaxIndex[i] = -1;
        rightMaxIndex[i] = -1;
        leftMaxOverlapLength[i] = 0;
        rightMaxOverlapLength[i] = 0;
        visited[i] = false;
    }
    
    for(i=0;i<overlapNumber;i++){
        leftOrientation[i] = -1;
        rightOrientation[i] = -1;
        leftIndex[i] = -1;
        rightIndex[i] = -1;
        overlapLength[i] = 0;
    }
    
    ifstream icin;
    icin.open(overlapInformationFile);
    char * overlapInfor = new char[1000];
    while(icin.getline(overlapInfor,1000)){
        long int start = 0;
        long int previous = 0;
        long int len = strlen(overlapInfor);
        OverlapInformation * overlap = new OverlapInformation[5];
        i = 0;
        j = 0;
        for(start = 0;start<len;start++){   
            if(overlapInfor[start]=='	'){
                overlap[i].infor[j] = '\0';
                i++;
                j = 0;
                continue;
            }
            overlap[i].infor[j] = overlapInfor[start];
            j++;
        }
        leftIndex[p] = atoi(overlap[0].infor);
        rightIndex[p] = atoi(overlap[1].infor);
        overlapLength[p] = atoi(overlap[2].infor);
        leftOrientation[p] = 1;
        if(overlap[3].infor[0]=='-'){
            leftOrientation[p] = 0;
        }
        rightOrientation[p] = 1;
        if(overlap[4].infor[0]=='-'){
            rightOrientation[p] = 0;
        }   
        p++;      
    }
    for(i=0;i<p;i++){
        
        //cout<<"id--"<<i<<"--------------"<<leftIndex[i]<<"--"<<rightIndex[i]<<"--"<<leftOrientation[i]<<"--"<<rightOrientation[i]<<endl;
        
        ContigSet * left = GetContigFromContigSet(contigSet, leftIndex[i]);
        ContigSet * right = GetContigFromContigSet(contigSet, rightIndex[i]);
        
        if(left == NULL || right == NULL){
            //cout<<"NULLLLL"<<endl;
        }
        
        long int leftLength = strlen(left->contig);
        long int rightLength = strlen(right->contig);
        
        //cout<<leftLength<<"--"<<rightLength<<endl;
        
        char * leftContig = left->contig;
        char * rightContig = right->contig;
        
        if(leftOrientation[i]==0){
            leftContig = new char[leftLength+1];
            ReverseComplement(left->contig,leftContig);
        }
        if(rightOrientation[i]==0){
            rightContig = new char[rightLength+1];
            ReverseComplement(right->contig,rightContig);
        }
        double temp11 = 0;
        for(int t = 0;t<readSetCount;t++){
            
            if(readSet[t].insertSize - lambda*readSet[t].var<overlapLength[i]||leftLength<readSet[t].insertSize- 2*readSet[t].readLength||rightLength<readSet[t].insertSize- 2*readSet[t].readLength){
                temp11 = 1;
                continue;
            }
        
            temp11 = WeightContigMerge(leftContig, rightContig, readSet, t, kmerLength, overlapLength[i], kmerLength);
            
            if(temp11!=1){
                break;
            } 
        } 
        
        if(temp11 == 1){
            if(leftOrientation[i] == 1){
                rightNumber[leftIndex[i]]++;
                if(rightMaxOverlapLength[leftIndex[i]] < overlapLength[i]){
                    rightMaxIndex[leftIndex[i]] = rightIndex[i];
                    rightMaxOverlapLength[leftIndex[i]] = overlapLength[i];
                }
            }else{
                leftNumber[leftIndex[i]]++;
                if(leftMaxOverlapLength[leftIndex[i]] < overlapLength[i]){
                    leftMaxIndex[leftIndex[i]] = rightIndex[i];
                    leftMaxOverlapLength[leftIndex[i]] = overlapLength[i];
                }
            }
            if(rightOrientation[i] == 1){
                leftNumber[rightIndex[i]]++;
                if(leftMaxOverlapLength[rightIndex[i]] < overlapLength[i]){
                    leftMaxIndex[rightIndex[i]] = leftIndex[i];
                    leftMaxOverlapLength[rightIndex[i]] = overlapLength[i];
                }
            }else{
                rightNumber[rightIndex[i]]++;
                if(rightMaxOverlapLength[rightIndex[i]] < overlapLength[i]){
                    rightMaxIndex[rightIndex[i]] = leftIndex[i];
                    rightMaxOverlapLength[rightIndex[i]] = overlapLength[i];
                }
            }
        }else{
            leftIndex[i] = -1;
            rightIndex[i] = -1;
        }
        
        if(leftOrientation[i]==0){
            delete [] leftContig;
        }
        if(rightOrientation[i]==0){
            delete [] rightContig;
        }
        //cout<<"id--end"<<endl;
    }
    //cout<<"ffffffffffffffffffffffffffff"<<endl;
    ofstream ocout;
    char * file = new char[10];
    strcpy(file, "overlap11.fa");
    remove(file);
    ocout.open(file);
    
    
    
    for(i=0;i<p;i++){
        
        if((rightNumber[leftIndex[i]]>1 && leftOrientation[i]==1 && rightMaxIndex[leftIndex[i]] != rightIndex[i])||(leftNumber[leftIndex[i]]>1 && leftOrientation[i]==0 && leftMaxIndex[leftIndex[i]] != rightIndex[i])){
            leftIndex[i] = -1;
            rightIndex[i] = -1;
        }
        if((leftNumber[rightIndex[i]]>1 && rightOrientation[i]==1 && leftMaxIndex[rightIndex[i]] != leftIndex[i])||(rightNumber[rightIndex[i]]>1 && rightOrientation[i]==0 && rightMaxIndex[rightIndex[i]] != leftIndex[i])){
            leftIndex[i] = -1;
            rightIndex[i] = -1;
        }
        
        ocout<<leftIndex[i]<<","<<rightIndex[i]<<","<<overlapLength[i]<<","<<leftOrientation[i]<<","<<rightOrientation[i]<<endl;
    }
    
    ContigMergeSet * head = new ContigMergeSet;
    ContigMergeSet * head1 = head;
    
    ocout<<endl;
    for(i=0;i<p;i++){
        
        if(leftIndex[i] == -1 ){
            continue;
        }
        if(visited[leftIndex[i]]!=false){
            continue;
        }
        
        ContigMergeSet * temp = new ContigMergeSet;
        ContigMerge * tempContigMerge = new ContigMerge;
        tempContigMerge->contigIndex = leftIndex[i];
        tempContigMerge->orientation = leftOrientation[i];
        tempContigMerge->overlapLength = overlapLength[i];
        temp->contigMergeSet = tempContigMerge;
        head1->next = temp;
        head1 = temp;
        
        bool end = 0;
        ocout<<leftIndex[i]<<",";
        visited[leftIndex[i]] = true;
        long int next = rightIndex[i];
        int ori = rightOrientation[i];
        long int start = i;
        
        ContigMerge * temp2 = new ContigMerge;
        temp2->contigIndex = next;
        temp2->orientation = ori;
        tempContigMerge->next = temp2;
        tempContigMerge = temp2;
        
        while(!end && visited[next]!=true){
            ocout<<next<<",";
            visited[next] = true;
            for(j=0;j<p;j++){
                if(leftIndex[j] == -1 || start==j){
                    continue;
                }
                if(leftIndex[j] == next && leftOrientation[j] == ori){
                    next = rightIndex[j];
                    ori = rightOrientation[j];
                    start = j;
                    
                    ContigMerge * temp2 = new ContigMerge;
                    temp2->contigIndex = next;
                    temp2->orientation = ori;
                    tempContigMerge->overlapLength = overlapLength[j];
                    tempContigMerge->next = temp2;
                    tempContigMerge = temp2;
                    
                    break;
                }
                if(rightIndex[j] == next && rightOrientation[j] == !ori){
                    next = leftIndex[j];
                    if(leftOrientation[j]==1){
                        ori = 0;
                    }else{
                        ori = 1;
                    }
                    start = j;
                   
                    ContigMerge * temp2 = new ContigMerge;
                    temp2->contigIndex = next;
                    temp2->orientation = ori;
                    tempContigMerge->overlapLength = overlapLength[j];
                    tempContigMerge->next = temp2;
                    tempContigMerge = temp2;
                   
                    break;
                }       
            }
            if(j==p){
                end = 1;
            }
        }
        
        for(j=i+1;j<p;j++){
            if(leftIndex[j]==leftIndex[i] || rightIndex[j] == leftIndex[i]){
                break;
            }
        }
        
        if(j==p){
            //cout<<endl;
            continue;
        }
        
        
        ocout<<"****";
        end = 0;
        next = leftIndex[i];
        visited[next] = false;
        start = i;
        ori = !leftOrientation[i];
        
        
        
        while(!end && visited[next]!=true){
            ocout<<next<<",";
            visited[next] = true;
            for(j=0;j<p;j++){
                if(leftIndex[j] == -1 || start==j){
                    continue;
                }
                if(leftIndex[j] == next && leftOrientation[j] == ori){
                    next = rightIndex[j];
                    if(rightOrientation[j] == 1){
                        ori = 1;
                    }else{
                        ori = 0;
                    }
                    
                    ContigMerge * temp2 = new ContigMerge;
                    temp2->contigIndex = next;
                    temp2->orientation = !rightOrientation[j];
                    temp2->overlapLength = overlapLength[j];
                    temp2->next = head1->contigMergeSet;
                    head1->contigMergeSet = temp2;
                    
                    
                    start = j;
                    break;
                }
                if(rightIndex[j] == next && rightOrientation[j] == !ori){
                    next = leftIndex[j];
                    ori = !leftOrientation[j];
                    start = j;
                   
                    ContigMerge * temp2 = new ContigMerge;
                    temp2->contigIndex = next;
                    temp2->orientation = leftOrientation[j];
                    temp2->overlapLength = overlapLength[j];
                    temp2->next = head1->contigMergeSet;
                    head1->contigMergeSet = temp2;
                   
                    break;
                }       
            }
            if(j==p){
                end = 1;
            }
        }
        
        
        ocout<<endl;
        
    }
    /*
    head1 = head->next;
    ocout<<endl;
    cout<<"ffffffffffffffffffffff"<<endl;
    while(head1!=NULL){
        cout<<"ttttttttttttttttt"<<endl;
        ContigMerge * temp = head1->contigMergeSet;
        cout<<"ttttttttttttttttt0000000000000000"<<endl;
        if(temp == NULL){
            cout<<"error!!!!!!!!!!!!!!!"<<endl;
        }
        while(temp->next!=NULL){
            cout<<temp->contigIndex<<","<<temp->orientation<<","<<temp->overlapLength<<"---";
            double temp11 = 0;
            for(int p = 0;p<readSetCount;p++){
            
                ContigSet * left = GetContigFromContigSet(contigSet, temp->contigIndex);
                ContigSet * right = GetContigFromContigSet(contigSet, temp->next->contigIndex);
                
                long int leftLength = strlen(left->contig);
                long int rightLength = strlen(right->contig);
                
                char * rightContig = right->contig;
                
                
                if(readSet[p].insertSize - lambda*readSet[p].var<temp->overlapLength||leftLength<readSet[p].insertSize- 2*readSet[p].readLength||rightLength<readSet[p].insertSize- 2*readSet[p].readLength){
                    temp11 = 1;
                    continue;
                }
                cout<<"dddd"<<endl;
                temp11 = WeightContigMerge(left->contig, rightContig, readSet, p, kmerLength, temp->overlapLength, kmerLength);
                cout<<"dddd000"<<endl;
                if(temp11!=1){
                    break;
                } 
            }   
            if(temp11==1){
                temp->merge = true;
            }
            //char * contig = new char[leftLength + rightLength - j + 1];
            //AppendRight(contig,left,right,j+1);
            //return contig;   
                   
            temp = temp->next;
        }
        ocout<<endl;
        head1 = head1->next;
    }
    */
    //cout<<"first"<<endl;
    head1 = head->next;
    char * contig = NULL;
    ContigSet * first = NULL;
    while(head1!=NULL){
        ContigMerge * temp = head1->contigMergeSet;
        first = NULL;
        contig = NULL;
        while(temp->next!=NULL){
            if(temp->merge == true || temp->merge==false){           
                ContigSet * right = GetContigFromContigSet(contigSet, temp->next->contigIndex);
                long int rightLength = strlen(right->contig);
                char * rightContig = right->contig;
                if(temp->next->orientation==0){
                    rightContig = new char[rightLength+1];
                    ReverseComplement(right->contig,rightContig);
                }
                
                if(contig==NULL){
                    ContigSet * left = GetContigFromContigSet(contigSet, temp->contigIndex);
                    long int leftLength = strlen(left->contig);
                    char * leftContig = left->contig;
                    
                    if(temp->orientation==0){
                        leftContig = new char[leftLength+1];
                        ReverseComplement(left->contig,leftContig);
                    }
                    
                    first = left;
                    contig = new char[leftLength + rightLength - temp->overlapLength + 1];
                    AppendRight(contig,leftContig,rightContig,temp->overlapLength+1);
                    if(temp->orientation==0){
                        delete [] leftContig;
                    }
                    if(temp->next->orientation==0){
                        delete [] rightContig;
                    }
                    leftContig = NULL;
                    rightContig = NULL;
                    delete [] left->contig;
                    delete [] right->contig;
                    left->contig = NULL;
                    right->contig = NULL;
                }else{
                    long int leftLength = strlen(contig);
                    char * contig1 = new char[leftLength + rightLength - temp->overlapLength + 1];
                    AppendRight(contig1,contig,rightContig,temp->overlapLength+1);
                    if(temp->next->orientation==0){
                        delete [] rightContig;
                    }
                    rightContig = NULL;
                    delete [] contig;
                    contig = contig1;
                    contig1 = NULL;
                    delete [] right->contig;
                    right->contig = NULL;
                }
            }else{
                if(first!=NULL){
                    first->contig = contig;
                    contig = NULL;
                    first = NULL;
                } 
            }  
                   
            temp = temp->next;
        }
        
        if(first!=NULL){
            first->contig = contig;
            contig = NULL;
            first = NULL;
        }
        
        ocout<<endl;
        head1 = head1->next;
    }
    //cout<<"second"<<endl;
    
    
    
    //exit(0);
    /*
    ContigSet * temp;
    ContigSet * first = contigSet;
    ContigSet * previous = NULL;
    
    MergeSubContig(contigSet, kmerLength);
    cout<<"ss000000"<<endl;
    
    char * fileName = new char[20];
    strcpy(fileName, "contigLong00.fa");
    WriteContigSetLong(contigSet, fileName);
    
    while(contigSet!=NULL){
        
        if(contigSet->contig==NULL){
            contigSet = contigSet->next;
            i ++;
            continue;
        }
        temp = first;
        j = 0;
        while(temp!=NULL){
            if(temp->contig==NULL || i==j){
                temp = temp->next;
                j++;
                continue;
            }
            char * tempContig = NULL;
            
            tempContig = LengthOfOverlapBetweenContig(temp->contig,contigSet->contig,readSet,readSetCount,kmerLength);
            
            if(tempContig!=NULL){
                
                delete []contigSet->contig;
                contigSet->contig = tempContig;
                delete []temp->contig;
                temp->contig = NULL;
                temp = temp->next;
                j++;
                continue;
            }
            char * tempReverseLeft = new char[strlen(temp->contig)+1];
            ReverseComplement(temp->contig,tempReverseLeft);
            
            tempContig = LengthOfOverlapBetweenContig(tempReverseLeft,contigSet->contig,readSet,readSetCount,kmerLength);
            
            if(tempContig!=NULL){
                
                delete []contigSet->contig;
                contigSet->contig = tempContig;
                delete []temp->contig;
                temp->contig = NULL;
                delete []tempReverseLeft;
                tempReverseLeft = NULL;
                j++;
                temp = temp->next;
                continue;
                
            }
            delete []tempReverseLeft;
            tempReverseLeft = NULL;
            char * tempReverseRight = new char[strlen(contigSet->contig)+1];
            ReverseComplement(contigSet->contig,tempReverseRight);
            
            tempContig = LengthOfOverlapBetweenContig(temp->contig,tempReverseRight,readSet,readSetCount,kmerLength);
            
            if(tempContig!=NULL){
                
                delete []contigSet->contig;
                contigSet->contig = tempContig;
                delete []temp->contig;
                temp->contig = NULL;
                delete []tempReverseRight;
                tempReverseRight = NULL;
                temp = temp->next;
                j++;
                continue;
                
            }
            temp = temp->next;
            j++;
            delete []tempReverseRight;
            tempReverseRight = NULL;
        }
        i++;
        contigSet = contigSet->next;
    }
    */
    
    
    
    ContigSet * temp;
    ContigSet * previous = NULL;
    //previous = NULL;
    first = contigSet;
    temp = first;
    //cout<<"ss1111111111"<<endl;
    
    while(first!=NULL){
        if(first->contig==NULL){
            if(previous != NULL){
                previous->next = first->next;
                delete first;
            }else{
                ContigSet * temp11 = first;
                first = first->next;
                temp11->next = NULL;
                temp = first;
                delete temp11;
                continue;
            }           
        }else{
            previous = first;
        }
        first = previous->next;
    }
    
    
    //remove(overlapInformationFile);
    //remove(file);
    return temp;
}










#endif
