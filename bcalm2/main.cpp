#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <cstring>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include "lm.h"
#include "ograph.h"
#include "debug.h"


using namespace std;

//whether use (k+1)mer file to solidly determin the adjacency, default:false
bool isBySolid = false;

void menu();

// this code crashes now; this is the old behavior when no arguments was specified, it runs some test
void test()
{
    srand(time(NULL));
    bool b(true);
    int n(1000);
    while(b)
    {

        int sys(0);
        int k(20);
        int m(2);

        cout<<"Test "<<"n:"<<n<<" k:"<<k<<" m:"<<2*m<<endl;
        createinputlm(n,k,"randomgenome");
        sys+=system("cat randomgenome | sort | uniq > input.dot");
        cout<<"GO!"<<endl;
        graph G(k),G2(k);
        G.importg("input.dot");
        G.debruijn();
        G.compress();
        G.print("output.dot");
        cout<<"GO!"<<endl;
        createoutfile("input.dot","outputlowmemory.dot",k,m);
        if(!checkfile("output.dot","outputlowmemory.dot",k)){
            cout<<"Errors occurred !"<<endl;
            b=false;
        }
        else{
            cout<<"Success !"<<endl;
            n*=10;
            cout<<n<<endl;
        }
        b=false;
    }

}







int main(int argc, char ** argv)
{

	if(argc<=1 || argc>8)
	{
		menu();
		exit(1);
	}
	string input="";
	string solidin="";
	string output("compacted.dot");
	int m = 5;
	if(argc>=2)
	{
		input = string(argv[1]);
	}
	
	//identify parameters
	for(int i = 2;i < argc;i++)
	{
		if(!strcmp("-l",argv[i]) && i+1<argc)
		{
			m = atoi(argv[i+1])/2;
			if(2*m >10)
			{
				cout << "-l minimizer value exceeds 10" << endl;
				exit(1);
			}
		}
		else if(!strcmp("-s",argv[i]) && i+1<argc)
		{
			solidin = string(argv[i+1]);
			ifstream fin(solidin);
			if(!fin)
			{
				cout << solidin << " not exist!" << endl;
				exit(1);
			}
			else
				isBySolid = true;
				fin.close();
		}
		else if(!strcmp("-o",argv[i]) && i+1<argc)
		{
			output = string(argv[i+1]);
		}
	}

	//check the largest # of open file
	if(testulimit(pow(4,m)+50))
	{
		int k(detectk(input));
		cout <<"k="<<k <<endl;
		if(k<=2*m){
			cout<<"k too low"<<endl;
		}
		else{
			GetKmerSetFromBcalm(solidin.c_str(),k+1);//get all k+1 mer
			createoutfile(input.c_str(),output.c_str(),k,m);
		}
	}
	else
	{
		cout<<"ulimit too low"<<endl;
	}

	return 0;
}

void menu()
{
	cout << "./bcalm: get all maximal simple paths from kmers\n\n";
	cout << "usage:\n";
	cout << "./bcalm input_kmer_file [-s solid_(k+1)mer_file] [-l minimiser_length] [-o output_file]\n\n";
	cout << "note: input_kmer_file can only be .dot format\n\n";
	cout << "details:\n";
	cout << "[-s solid_(k+1)mer_file] (k+1)mer file to solidfy the result\n";
	cout << "[-l minimizer_length] minimizer length, i.e. # of nucleotides of minimizer\n";
	cout << "default: 10\n";
	cout << "[-o output_file] output_file, i.e. maximal simple paths, in .dot format\n";
	cout << "default: compacted.dot\n";
}
