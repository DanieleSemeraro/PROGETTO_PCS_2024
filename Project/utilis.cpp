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

vector<double> ImportDFN(string filename,int n,vector<double> &FractureId,vector<double> &NumVertices,MatrixXd &Vertices)
{
    ifstream fin(filename);
    string line;
    int c=0;
    int p;

    while (getline(fin,line)) {
        if (line.empty() || n==stoi(line)) {
            continue;
        }
        else if(line[0]=='#' && line[2]=='N'){
            continue;
        }
        else if(line[0]=='#' && line[2]=='F'){
            c=1;
        }
        else if(line[0]=='#' && line[2]=='V'){
            c=2;
        }

        if(c==0){
            continue;
        }
        else if(c==1){
            for (unsigned i = 0; i < line.size(); i++) {
                if(line[i] == ';'){
                    p=i;
                }
            }
            FractureId.push_back(stod(line.substr(0,p)));
            NumVertices.push_back(stod(line.substr(p+1)));
            c=0;
        }
        else if(c==2){
            c=0;

        }







    }






    fin.close();

    return FractureId;

}








