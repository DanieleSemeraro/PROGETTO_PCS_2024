#include "utilis.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include"Eigen/Eigen"

using namespace std;
using namespace Eigen;

ostream& operator<<(ostream& os, const vector<double> a)
{
    for (size_t i = 0; i < a.size(); i++) {
        os<< a[i]<< " ";
    }

    return os;
}

ostream& operator<<(ostream& os, const vector<MatrixXd> a)
{
    for (size_t i = 0; i < a.size(); i++) {
        os<<"Id frattura: "<<i<<endl<<"Matrice vertici: "<<endl<< a[i]<<endl;
        cout<<" "<<endl;
    }
    return os;
}

vector<double> ImportDFN(string filename,int n,vector<double> &FractureId,vector<double> &NumVertices,vector<MatrixXd> &ListVertices)
{
    ifstream fin(filename);
    string line;
    int c=0;//serve come contatore per scegliere cosa memorizzare e dove
    vector<double> p;//serve a memorizzare la posizione del punto e virgola
    int a=0;//altro contatore utile alla memorizzazione dei vertici
    vector<double> pv;//serve a memorizzare la posizione del punto e virgola nei vertici
    MatrixXd Vertices(3,4);

    while (getline(fin,line)) {
        if (line.empty() || n==atoi(line.c_str())) {
            continue;
        }
        else if(line[0]=='#' && line[2]=='N'){
            continue;
        }
        else if(line[0]=='#' && line[2]=='F'){
            c=1;
            continue;
        }
        else if(line[0]=='#' && line[2]=='V'){
            c=2;
            continue;
        }

        if(c==0){
            continue;
        }
        else if(c==1){
            for (unsigned i = 0; i < line.size(); i++) {
                if(line[i] == ';'){
                    p.push_back(i);
                }
            }
            FractureId.push_back(stod(line.substr(0,p[0])));
            NumVertices.push_back(stod(line.substr(p[0]+1)));
            c=0;
            p.clear();
        }
        else if(c==2){
            for (unsigned i = 0; i < line.size(); i++) {
                if(line[i] == ';'){
                    pv.push_back(i);
                }
            }
            Vertices(a,0)=stod(line.substr(0,pv[0]));
            for (int j = 1; j < NumVertices.back()-1; j++) {
                Vertices(a,j)=stod(line.substr(pv[j-1]+1,pv[j]));
            }

            int b=NumVertices.back()-1;
            Vertices(a,b)=stod(line.substr(pv[b-1]+1));

            a=a+1;
            pv.clear();
            if(a==3){
                ListVertices.push_back(Vertices);
                a=0;
                c=0;
            }
        }
    }

    fin.close();

    return FractureId;

}








