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

MatrixXd SystemSolve(vector<MatrixXd> &ListVertices,Vector3d P0,Vector3d P1,Vector3d P2,Vector3d P3,MatrixXd A,Vector3d b,int i,int j,int &c,vector<double> NumVertices){
    Vector2d sol;
    double norma;
    MatrixXd punti(3,2);
    Vector3d P;

    for (int k = 0; k < NumVertices[i]; ++k) {
        if(k==3){
            P0=ListVertices[i].col(k);
            P1=ListVertices[i].col(0);
        }
        else{
            P0=ListVertices[i].col(k);
            P1=ListVertices[i].col(k+1);
        }
        for (int z = 0; z < NumVertices[j]; ++z) {
            if(z==3){
                P2=ListVertices[j].col(z);
                P3=ListVertices[j].col(0);
            }
            else{
                P2=ListVertices[j].col(z);
                P3=ListVertices[j].col(z+1);
            }
            A.col(0)=P1-P0;
            A.col(1)=P3-P2;
            b=P2-P0;
            norma=((P1-P0).cross(P3-P2)).norm();
            if(norma>0)
            {
                sol=A.colPivHouseholderQr().solve(b);
                if(A*sol==b)
                {
                    P=P0+sol(0)*(P1-P0);
                    //cout<<"punto "<<P<<endl;
                    if((P0(0)<=P(0) && P(0)<=P1(0)) || (P1(0)<=P(0) && P(0)<=P0(0)))
                    {
                        if((P0(1)<=P(1) && P(1)<=P1(1)) || (P1(1)<=P(1) && P(1)<=P0(1)))
                        {
                            if((P0(2)<=P(2) && P(2)<=P1(2)) || (P1(2)<=P(2) && P(2)<=P0(2)))
                            {
                                //cout<<"punto p "<<P<<endl;
                                punti.col(c)=P;
                                c=c+1;
                            }
                        }


                    }


                }

            }


        }


    }

    return punti;

}

vector<double>CalcoloTracce(int n,vector<double> &FractureId,vector<double> &NumVertices,vector<MatrixXd> &ListVertices,vector<MatrixXd> &ListCord){
    vector<double> TraceId;
    Vector3d P0;
    Vector3d P1;
    Vector3d P2;
    Vector3d P3;
    MatrixXd A(3,2);
    Vector3d b;
    MatrixXd punti(3,2);

    int c=0;//serve a verificare se abbiamo ottenuto un inntersezione passante per entrambi i bordi

    for (int i = 0; i < n-1; ++i) {
        for (int j = i+1; j < n; ++j) {
            punti=SystemSolve(ListVertices,P0,P1,P2,P3,A,b,i,j,c,NumVertices);
            if(c==2){
                ListCord.push_back(punti);
                TraceId.push_back(ListCord.size()-1);
                cout<<"TraceId: "<<TraceId<<" FractureId: "<<FractureId[i]<<" "<<FractureId[j]<<endl<<" Cordinates: "<<punti.col(0).transpose()<<"; "<<punti.col(1).transpose()<<endl;
            }
            c=0;
        }

    }

    return TraceId;

}








