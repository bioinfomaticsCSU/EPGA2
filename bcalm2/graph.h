#ifndef GRAPH_H_INCLUDED
#define GRAPH_H_INCLUDED

#include "fstream"
#include <cstring>
#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <utility>      // std::swap
#include <time.h>
#include <string>
#include "stdio.h"
#include "stdlib.h"





using namespace std;

typedef struct GraphNode{
    long int index;
    struct GraphNode * next;
    GraphNode(){
        index = -1;
        next = NULL;
    }
}GraphNode;


#pragma pack(2)
typedef struct DBGraph{
    char * contig;
    GraphNode * outNode;
    GraphNode * inNode;
    DBGraph(){
        contig = NULL;
        outNode = NULL;
        inNode = NULL;
    }
}DBGraph;
#pragma pack ()

#pragma pack(2)
typedef struct DBGraphHead{
    DBGraph * deBruijnGraph;
    unsigned long int nodeNumber;
    DBGraphHead(){
        deBruijnGraph = NULL;
        nodeNumber = 0;
    }
}DBGraphHead;
#pragma pack ()

//current system time
const string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}

DBGraph * readDBG(char* contig){
      DBGraph* graph = new DBGraph();
      graph->contig = new char[strlen((const char*)contig)+1];
      if((contig[0] == '\n') || (contig[0] == '\r')){
          strcpy(graph->contig,contig+1);
      }else{
          strcpy(graph->contig,contig);
      }
      graph->contig[strlen((const char*)contig)]='\0';
      return graph;  
}

//find the lefemost k-mer of dbg
char * LKmer(DBGraph * dbg, int k, char* l_kmer){
     strncpy(l_kmer, dbg->contig,k);
     l_kmer[k] = '\0';
     return l_kmer;    
}

//find the rightmost k-mer of dbg
char * RKmer(const DBGraph * dbg, int k){
     unsigned int contig_size = strlen(dbg->contig);
     return  dbg->contig+(contig_size-k);
}

void printInOut(vector<DBGraph*> vec,  int index, ofstream& out){
     if((vec.size()==0) || (index <0) || (index >=vec.size())){
            return;       
     }
     GraphNode* currentNode = vec[index]->inNode;
     //cout << "--------------------------"<<endl;
     //cout << index << ": " << vec[index]->contig << endl;
	 out << index <<endl;
     out << "InNode:" << endl;
	 
     while(currentNode != NULL){
          out << currentNode->index<<endl;
			  //<< "::" <<vec[currentNode->index]->contig<<endl;
          currentNode = currentNode->next;              
     }
     
     currentNode = vec[index]->outNode;
     out << "OutNode:" <<endl;
     while(currentNode != NULL){
          out << currentNode->index<<endl;
			  //<< "::" <<vec[currentNode->index]->contig<<endl;
          currentNode = currentNode->next;              
     }
	 out<<"*****"<<endl;
}

//compare tail string of g1 and head string of g2, both of length cmplen
int comBinaryDBG(const DBGraph * g1, const DBGraph * g2, int cmplen)
{ 
    int g1_size = strlen(g1->contig);
    char * g1_tail = g1->contig+(g1_size-cmplen);
    //cout <<"g1_tail:" << g1_tail << endl;
    //cout <<"g2:     "<<g2->contig << endl;
    return strncmp(g1_tail,g2->contig, cmplen);
}

//compare if two DBGraph contain an identical contig
bool compDBG(const DBGraph * g1, const DBGraph * g2) {
        if(g1 == NULL || g2 == NULL){
              return 0;      
        }
        if(strcmp(g1->contig,g2->contig) > 0){
            return false;                                                                                            
        }else{
            return true;
        }
}

//search the position of contig whose head has overlap of length overlap_lengh with vec[index]
//-1 returned if not found

int BinarySearch(vector<DBGraph*>& vec, int index, int overlap_length)  
{   
    int len = vec.size();
    if ( len <= 0)  
        return -1;  
    
    int low = 0;  
    int high = len - 1;  
    while (low <= high)  
    {  
        int mid = low + (high - low) / 2;
        int flag = comBinaryDBG(vec[index],vec[mid],overlap_length);
        if ( flag == 0)  
            return mid;  
        else if (flag < 0)  
            high = mid - 1;  
        else  
            low = mid + 1;  
    }  
  
    return -1;  
}

//form a connection from vec[in] to vec[out] and possiblely neighbor contigs vec[potential_out]
void connectDBG(vector<DBGraph*>& vec,  int in, int out,int overlap_len){
     /*if((in < 0) || (in >=vec.size()) || (out <0) || (out >= vec.size()) || (overlap_len <= 0)){
            return -1;
     }*/
     //int count = 0;//#of overlaped contigs with vec[in]
	 //int low_bound = out - 3;
     int up_bound = out + 3;
	 //if(low_bound < 0) low_bound = 0;
     if(up_bound > vec.size()-1) up_bound = vec.size()-1;
     //cout << "low_bound: " << low_bound <<endl;
     //cout << "up_bound: " << up_bound <<endl;
     for( int potential_out = out; potential_out <= up_bound; potential_out++){        
         if((in != potential_out) && comBinaryDBG(vec[in],vec[potential_out],overlap_len) ==0)//confirm and edge from vec[in] to vec[out]
         {
             //cout << "Start to connect " << in << " to " << potential_out <<endl;                                                    
             GraphNode* inNode = new GraphNode();
             GraphNode* outNode = new GraphNode();
             inNode->index = in;
             outNode->index = potential_out;
             //find the last outNode
             GraphNode* currentOutNode = vec[in]->outNode;
             if(currentOutNode == NULL)
             {
                   vec[in]->outNode = outNode;                    
             }
             else
             {
                   while(currentOutNode->next != NULL)
                   {
                       currentOutNode = currentOutNode->next;                                          
                   }
                   currentOutNode->next = outNode;
             }
             //find the last inNode
             GraphNode* currentInNode = vec[potential_out]->inNode;
             if(currentInNode == NULL)
             {
                 vec[potential_out]->inNode = inNode;                      
             }
             else
             {
                 while(currentInNode->next != NULL)
                 {
                     currentInNode = currentInNode->next;
                 }
                 currentInNode->next = inNode;     
             }          
          }
          else{  //no more edges from vec[in] to vec[out] 
               break;   
          }                        
     }  
}

//form connections to all pairs of DBGraphs with overlap of length k in vec
/*void buildDBG(vector<DBGraph*> vec,int k){      
     for(unsigned int in = 0;in < vec.size();++in){
             if(in%10000 == 0){
                       cout << currentDateTime() << endl;
                       cout << "---------BuildingDBG: Loop " << in << "-----------" <<endl;
                       
             }
             int out =  BinarySearch(vec, in, k-1);
             
             if(out != -1){
                   //cout << in <<" found out " << out << endl; 
                   connectDBG(vec,in,out,k-1);        
             }
     }    
}*/


void buildDBG(vector<DBGraph*>& vec, vector<int>& tail_array, int k){
     int i = 0, j = 0;
     int flag = 0;
     while((i <= vec.size()-1) && (j <= tail_array.size()-1))
     {
              /*
			  if(j%10000 == 1){
                       cout << currentDateTime() << endl;
                 
					   cout << "---------BuildingDBG: Loop " << j << "-----------" <<endl;
                       
                }
			   */
              //flag = strncmp(tail,head,k-1);
              flag = comBinaryDBG(vec[tail_array[j]],vec[i],k-1);
              if(flag < 0){ //tail is less than head
                   j++;  
              }
              else if(flag > 0) //head is less than tail
              {
                   i++;
              }
              else          //head and tail are identical
              {
                   connectDBG(vec,tail_array[j],i,k-1);
                   j++;     
              }   
     }     
}



#endif
