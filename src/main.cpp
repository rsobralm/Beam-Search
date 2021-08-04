#include "readData.h"
//#include "CustoIn.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
//#include "infoSeq.h"
#include <sys/timeb.h>
#include <sys/resource.h>
#include <fstream>
#include <tuple>

using namespace std;

//ofstream out("arquivo.txt");
double ** mJobs; // matriz de adjacencia
double ** mSetupTimes; // matriz reorganizada;
int n; // quantidade total de vertices

double timeSequenceTime = 0;

int binarySearch(vector<int> &vec, int posInit, int posEnd, int compare);
void beamSearch(pair< vector<int>, double> &bestSol, int totalBranchesPerNode, int totalBranchesPerLevel);
void printSolution(vector<int> solucao, double **mJobs, double ** mSetupTimes);
bool compBSLWT(tuple <pair<double, double>, vector<int>, vector<int>, double > a, tuple <pair<double, double>, vector<int>, vector<int>, double> b);
bool compBSLCT(tuple <pair<double, double>, vector<int>, vector<int>, double > a, tuple <pair<double, double>, vector<int>, vector<int>, double > b);
double idleTime(vector<int> s, double cTime, int node, double ** mJobs, double **mSetupTimes);
double cpuTime();
double sequenceTime(vector<int> s, double prevCompTime, vector<int> remainingVertices, double ** mJobs, double **mSetupTimes);
int genRandomInteger(int min, int max);
bool compMakespan(tuple <pair<double, double>, vector<int>, vector<int>, double > a, tuple <pair<double, double>, vector<int>, vector<int>, double > b);


int main(int argc, char** argv) {
    //Lê os dados da instância
    readData(argc, argv, &n, &mJobs, &mSetupTimes);

    //Declara variáveis
    pair< vector<int>, double> bestSol, sol;
    bestSol.second = numeric_limits<double>::max();
    sol.first = {};
    sol.second = 0;
    std::string arg1(argv[1]);

    unsigned seed = time(0);
        //cout << "\nseed: " << seed << endl;
    srand(seed);

    double execTime = cpuTime();
    //Chama Beam Search e mostra solução
    beamSearch(bestSol, 5, 500);
    //beamSearch(bestSol, 5, 500);

    execTime = cpuTime() - execTime;
    //printSolution(bestSol.first, mJobs, mSetupTimes);
    cout << arg1 << ";" << bestSol.second << ";" << execTime << endl;

    //out.close();

    return 0;
}

//Algorítimos
void beamSearch(pair< vector<int>, double> &bestSol, int totalBranchesPerNode, int totalBranchesPerLevel){
    int removePosition = 0;

    vector<int> listOfNodesRemaining;
    for (int i = 1; i <= n; i++){
        listOfNodesRemaining.push_back(i);
    }

    vector <tuple<pair<double, double>, vector<int>, vector<int>, double >> listOfSols;
    vector <tuple<pair<double, double>, vector<int>, vector<int>, double >> nodeListOfSols;
    vector <tuple<pair<double, double>, vector<int>, vector<int>, double >> levelListOfSols;
    vector<int> empty;

    vector <tuple<pair<double, double>, vector<int>, vector<int>, double >> emptSols;
    tuple <pair<double, double>, vector<int>, vector<int>, double > sol({0, 0}, empty, listOfNodesRemaining, 0.0);
    listOfSols.push_back(sol);                                                  

    for (int i = 0; i < n-1; i++){
        //cout << "Nivel: " << i << endl;

        //Limpa a lista de nós por nível para nova análise
        levelListOfSols.clear();

        for (int j = 0; j < listOfSols.size(); j++){
            //Limpa a lista de nós para uma nova análise
            nodeListOfSols.clear();
            
            //Obtem todas as possíveis expações a partir do nó J 
            for (int k = 0; k < get<2>(listOfSols[j]).size(); k++){
                nodeListOfSols.push_back(listOfSols[j]); 
                get<0>(nodeListOfSols[k]).first = idleTime(get<1>(nodeListOfSols[k]), get<3>(nodeListOfSols[k]), 
                                                           get<2>(listOfSols[j])[k], mJobs, mSetupTimes);
                get<3>(nodeListOfSols[k]) += get<0>(nodeListOfSols[k]).first + mJobs[get<2>(listOfSols[j])[k]][2]; // att comp time                                       
                get<1>(nodeListOfSols[k]).push_back(get<2>(listOfSols[j])[k]);
                removePosition = binarySearch(get<2>(nodeListOfSols[k]), 0, get<2>(nodeListOfSols[k]).size()-1,
                                              get<1>(nodeListOfSols[k])[ get<1>(nodeListOfSols[k]).size()-1 ]);
                get<2>(nodeListOfSols[k]).erase( get<2>(nodeListOfSols[k]).begin() + removePosition);
            }

            //Ordena pelo Waiting Time a lista de nós obtidos
            sort(nodeListOfSols.begin(), nodeListOfSols.end(), compBSLWT);

            //Escolhe os W melhores nós para expadir
            
            for (int k = 0; k < nodeListOfSols.size(); k++){
                levelListOfSols.push_back(nodeListOfSols[k]);  
                if(k == totalBranchesPerNode-1) break; 
            } 
            
            
        }  

        //Caso a quantidade de nós para expansão ultrapasse o limite por nível, ordena pelo LBMakespan
        if(levelListOfSols.size() > totalBranchesPerLevel){
            for(int k = 0; k < levelListOfSols.size(); k++){
                get<0>(levelListOfSols[k]).second = sequenceTime(get<1>(levelListOfSols[k]), get<3>(levelListOfSols[k]), get<2>(levelListOfSols[k]), 
                                                                mJobs, mSetupTimes);
            }
            sort(levelListOfSols.begin(), levelListOfSols.end(), compBSLCT);
        }

        //Limpa a lista de nós para expanção para nova análise  
        listOfSols.clear();

       
        //Expande os N melhores nós do nível
        for (int j = 0; j < levelListOfSols.size(); j++){
            listOfSols.push_back(levelListOfSols[j]);
            if(j == totalBranchesPerLevel) break;
        }
        

        //cout << "Best: " << get<0>(listOfSols[0]).second << endl << endl;
        //cout << "C_time: " << get<3>(listOfSols[0]) << endl;
    }

    for (int j = 0; j < listOfSols.size(); j++){
        get<0>(listOfSols[j]).first = idleTime(get<1>(listOfSols[j]), get<3>(listOfSols[j]), 
                                                           get<2>(listOfSols[j])[0], mJobs, mSetupTimes);
        get<3>(listOfSols[j]) += get<0>(listOfSols[j]).first + mJobs[get<2>(listOfSols[j])[0]][2];
        get<1>(listOfSols[j]).push_back(get<2>(listOfSols[j])[0]);     
    }


    sort(listOfSols.begin(), listOfSols.end(), compMakespan);
    bestSol.first = get<1>(listOfSols[0]);
    bestSol.second = get<3>(listOfSols[0]);

    //Obtem melhor solução ao final do algorítimo
    /*bestSol.first = get<1>(listOfSols[0]);
    bestSol.second = get<0>(listOfSols[0]).second;*/
}

int binarySearch(vector<int> &vec, int posInit, int posEnd, int compare) { 
    if (posEnd >= posInit) { 
        int midlePos = posInit + (posEnd - posInit) / 2; 
  
        if (vec[midlePos] == compare) 
            return midlePos; 
  
        if (vec[midlePos] > compare) 
            return binarySearch(vec, posInit, midlePos - 1, compare); 
  
        return binarySearch(vec, midlePos + 1, posEnd, compare); 
    } 
  
    return -1; 
} 

//Funções de ordenação
bool compBSLWT(tuple <pair<double, double>, vector<int>, vector<int>, double > a, tuple <pair<double, double>, vector<int>, vector<int>, double > b) {
    //Ordena pelo Waiting Time

    //Esta comparação apresenta melhoras em algumas instâncias
    //return get<0>(a).first*get<0>(a).second < get<0>(b).first*get<0>(b).second; 
    //Comparação do artigo
    return get<0>(a).first < get<0>(b).first;
}

bool compBSLCT(tuple <pair<double, double>, vector<int>, vector<int>, double > a, tuple <pair<double, double>, vector<int>, vector<int>, double > b) {
    //Ordena pelo LB Completion Time

   return get<0>(a).second < get<0>(b).second;
}

//Funçoes de cálculos
double idleTime(vector<int> s, double cTime, int node, double ** mJobs, double **mSetupTimes){
    double idleTime;
    int lastNode = (s.size() > 0) ? s[s.size()-1] : 0;
    double setupTime = mSetupTimes[lastNode][node];

    if(mJobs[node][1] <= cTime)
        idleTime = 0;
    else
        idleTime = mJobs[node][1] - cTime;

    idleTime += setupTime;

    return idleTime;
}

double sequenceTime(vector<int> s, double prevCompTime, vector<int> remainingVertices, double ** mJobs, double **mSetupTimes){
    double LBMakespan;
    double cTime = prevCompTime;
    double minimunReleaseDate = numeric_limits<double>::max();
    double remainingProcessingTime = 0;
    double remainingLBSetupTime = 0;
    double minimumSetupTime = numeric_limits<double>::max();
    vector<int> unionVertices = remainingVertices;
    unionVertices.push_back(s[s.size()-1]);

    for(int i = 0; i < remainingVertices.size(); i++){

        remainingProcessingTime += mJobs[remainingVertices[i]][2]; //Calcula o processing time total restante

        if(mJobs[remainingVertices[i]][1] < minimunReleaseDate) minimunReleaseDate = mJobs[remainingVertices[i]][1]; //Calcula menor release date dos nós ainda não inseridos

        for(int j = 0; j < unionVertices.size(); j++){ //Calcula o lower bound dos setup times
            if(unionVertices[j] == remainingVertices[i]) continue;
            if(mSetupTimes[unionVertices[j]][remainingVertices[i]] < minimumSetupTime)
                minimumSetupTime = mSetupTimes[unionVertices[j]][remainingVertices[i]];
        }
        remainingLBSetupTime += minimumSetupTime;
        minimumSetupTime = numeric_limits<double>::max();
    }

    if( minimunReleaseDate == numeric_limits<double>::max()) minimunReleaseDate = 0;

    //Calcula o lower bound do makespan 
    LBMakespan = (cTime > minimunReleaseDate) ? cTime : minimunReleaseDate;
    LBMakespan += remainingProcessingTime;
    LBMakespan += remainingLBSetupTime;

    return LBMakespan;
}


//Funções de utilidades
void printSolution(vector<int> solucao, double **mJobs, double ** mSetupTimes){
    cout << "Solution: [ ";
    for(int i = 0; i < solucao.size(); i++){
        cout << solucao[i] << " ";
    }
    cout << " ]" << endl;
    /*for(int i = 0; i < n; i++){
        cout << mJobs[solucao[i]][1] << " ";
    }*/

    double cTime = mSetupTimes[0][solucao[0]] + mJobs[solucao[0]][2] + mJobs[solucao[0]][1];
    double totalWT = mJobs[solucao[0]][1];
    /*for(int i = 0; i < s.size(); i++){
        cout << s[i] << " ";
    }*/
    cout << totalWT;
    for(int i = 0, j = 1; j < solucao.size(); i++ ,j++){
        //cout<< i << " " << j << endl;
        if(cTime >= mJobs[solucao[j]][1]){
            cTime += mSetupTimes[solucao[i]][solucao[j]] + mJobs[solucao[j]][2];
            cout << " " << 0;
        }
        else{
            totalWT = mJobs[solucao[j]][1] - cTime;
            cout << " " << totalWT;
            cTime += mSetupTimes[solucao[i]][solucao[j]] + mJobs[solucao[j]][2] + (mJobs[solucao[j]][1] - cTime);
        }
    }

    cout << endl << "TotalWT: " <<  totalWT << endl;
    cout << "CTime: " <<  cTime << endl;
}

double cpuTime() {
	static struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);
	return ((double)usage.ru_utime.tv_sec)+(((double)usage.ru_utime.tv_usec)/((double)1000000));
}

int genRandomInteger(int min, int max){
   return min + (rand() % (max - min + 1));
}

bool compMakespan(tuple <pair<double, double>, vector<int>, vector<int>, double > a, tuple <pair<double, double>, vector<int>, vector<int>, double > b) {
    //Ordena pelo LB Completion Time
    return get<3>(a) < get<3>(b);
}
