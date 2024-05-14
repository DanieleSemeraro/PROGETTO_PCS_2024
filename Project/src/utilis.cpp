#include "utilis.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include"Eigen/Eigen"
#include <cmath>

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

vector<double>CalcoloTraccePassanti(int n,vector<double> &FractureId,vector<double> &NumVertices,vector<MatrixXd> &ListVertices,vector<MatrixXd> &ListCord){
    vector<double> TraceId;
    Vector3d P0;
    Vector3d P1;
    Vector3d P2;
    Vector3d P3;
    MatrixXd A(3,2);
    Vector3d b;
    MatrixXd punti(3,2);//matrici che immagazinano le due cordinate

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

Vector3d SoluzioneSistema(Vector3d &P,vector<MatrixXd> &ListVertices,Vector3d P0,Vector3d P1,Vector3d P2,Vector3d P3,MatrixXd A,Vector3d b,int i,int j,vector<double> NumVertices){
    Vector3d u;
    Vector3d v;
    Vector3d n1;
    Vector3d n2;
    Vector3d t;//direzione retta di intersezione
    double d1;
    double d2;
    MatrixXd alpha(3,2);
    Vector2d sol;
    Vector3d p;
    Vector3d sis2;
    Vector3d sis1;
    //MatrixXd punti(3,2);
    int c1=0;//contatore
    int c2=0;//contatore

    P0=ListVertices[i].col(0);
    P1=ListVertices[i].col(1);
    P2=ListVertices[i].col(2);
    u=P2-P0;
    v=P1-P0;
    n1=(u.cross(v))/(u.norm()*v.norm());
    d1=n1.dot(P0);

    P0=ListVertices[j].col(0);
    P1=ListVertices[j].col(1);
    P2=ListVertices[j].col(2);
    u=P2-P0;
    v=P1-P0;
    n2=(u.cross(v))/(u.norm()*v.norm());
    d2=n2.dot(P0);

    if((n1.cross(n2)).norm()>0){
        t=n1.cross(n2);
        A.row(0)=n1;
        A.row(1)=n2;
        A.row(2)=t;
        b<<d1,d2,0;
        if(A.determinant()!=0){
            P=A.colPivHouseholderQr().solve(b);
            sis1=A*P;
            for (int z = 0; z < 3; ++z) {
                b(z)=round(b(z)*pow(10,6))/pow(10,6);
                sis1(z)=round(sis1(z)*pow(10,6))/pow(10,6);
            }
            if((sis1(0)<=b(0)+0.000001 && sis1(0)>=b(0)-0.000001) && (sis1(1)<=b(1)+0.000001 && sis1(1)>=b(1)-0.000001) && (sis1(2)<=b(2)+0.000001 && sis1(2)>=b(2)-0.000001)){
                for (int k = 0; k < NumVertices[i]; ++k) {
                    if(k==3){
                        alpha.col(0)=t;
                        alpha.col(1)=ListVertices[i].col(0)-ListVertices[i].col(k);
                        P0=ListVertices[i].col(k);
                        P1=ListVertices[i].col(0);
                        b=ListVertices[i].col(k)-P;

                    }
                    else{
                        alpha.col(0)=t;
                        alpha.col(1)=ListVertices[i].col(k+1)-ListVertices[i].col(k);
                        P0=ListVertices[i].col(k);
                        P1=ListVertices[i].col(k+1);
                        b=ListVertices[i].col(k)-P;
                    }
                    if((t.cross(P1-P0)).norm()>0){
                        sol=alpha.colPivHouseholderQr().solve(b);
                        sis2=alpha*sol;
                        for (int z = 0; z < sis2.size(); ++z) {
                            if(sis2(z)==-0){
                                sis2(z)=0;
                            }
                        }
                        for (int z = 0; z < b.size(); ++z) {
                            if(b(z)==-0){
                                b(z)=0;
                            }
                        }
                        for (int z = 0; z < 3; ++z) {
                            sis2(z)=round(sis2(z)*pow(10,6))/pow(10,6);
                            b(z)=round(b(z)*pow(10,6))/pow(10,6);
                        }
                        if((sis2(0)<=b(0)+0.000001 && sis2(0)>=b(0)-0.000001) && (sis2(1)<=b(1)+0.000001 && sis2(1)>=b(1)-0.000001) && (sis2(2)<=b(2)+0.000001 && sis2(2)>=b(2)-0.000001)){
                            p=P+sol(0)*t;
                            //cout<<"p"<<p<<endl;
                            for (int z = 0; z < 3; ++z) {
                                p(z)=round(p(z)*pow(10,6))/pow(10,6);
                            }
                            //cout<<"gigi "<<p<<endl;
                            for (int z = 0; z < 3; ++z) {
                                P0(z)=round(P0(z)*pow(10,6))/pow(10,6);
                            }
                            for (int z = 0; z < 3; ++z) {
                                P1(z)=round(P1(z)*pow(10,6))/pow(10,6);
                            }
                            if((P0(0)<=p(0) && p(0)<=P1(0)) || (P1(0)<=p(0) && p(0)<=P0(0)))
                            {
                                if((P0(1)<=p(1) && p(1)<=P1(1)) || (P1(1)<=p(1) && p(1)<=P0(1)))
                                {
                                    if((P0(2)<=p(2) && p(2)<=P1(2)) || (P1(2)<=p(2) && p(2)<=P0(2)))
                                    {
                                        //cout<<"intersez: "<<c1<<endl<<p<<endl;
                                        c1=c1+1;
                                    }
                                }
                            }
                        }
                    }
                }
                if(c1==2){
                    for (int k = 0; k < NumVertices[j]; ++k) {
                        if(k==3){
                            alpha.col(0)=t;
                            alpha.col(1)=ListVertices[j].col(0)-ListVertices[j].col(k);
                            P0=ListVertices[j].col(k);
                            P1=ListVertices[j].col(0);
                            b=ListVertices[j].col(k)-P;

                        }
                        else{
                            alpha.col(0)=t;
                            alpha.col(1)=ListVertices[j].col(k+1)-ListVertices[j].col(k);
                            P0=ListVertices[j].col(k);
                            P1=ListVertices[j].col(k+1);
                            b=ListVertices[j].col(k)-P;
                        }
                        if((t.cross(P1-P0)).norm()>0){
                            sol=alpha.colPivHouseholderQr().solve(b);
                            sis2=alpha*sol;
                            for (int z = 0; z < sis2.size(); ++z) {
                                if(sis2(z)==-0){
                                    sis2(z)=0;
                                }
                            }
                            for (int z = 0; z < b.size(); ++z) {
                                if(b(z)==-0){
                                    b(z)=0;
                                }
                            }
                            for (int z = 0; z < 3; ++z) {
                                sis2(z)=round(sis2(z)*pow(10,6))/pow(10,6);
                                b(z)=round(b(z)*pow(10,6))/pow(10,6);
                            }
                            if((sis2(0)<=b(0)+0.000001 && sis2(0)>=b(0)-0.000001) && (sis2(1)<=b(1)+0.000001 && sis2(1)>=b(1)-0.000001) && (sis2(2)<=b(2)+0.000001 && sis2(2)>=b(2)-0.000001)){
                                p=P+sol(0)*t;
                                for (int z = 0; z < 3; ++z) {
                                    p(z)=round(p(z)*pow(10,6))/pow(10,6);
                                }
                                //cout<<"gigi "<<p<<endl;
                                for (int z = 0; z < 3; ++z) {
                                    P0(z)=round(P0(z)*pow(10,6))/pow(10,6);
                                }
                                for (int z = 0; z < 3; ++z) {
                                    P1(z)=round(P1(z)*pow(10,6))/pow(10,6);
                                }
                                if((P0(0)<=p(0) && p(0)<=P1(0)) || (P1(0)<=p(0) && p(0)<=P0(0)))
                                {
                                    if((P0(1)<=p(1) && p(1)<=P1(1)) || (P1(1)<=p(1) && p(1)<=P0(1)))
                                    {
                                        if((P0(2)<=p(2) && p(2)<=P1(2)) || (P1(2)<=p(2) && p(2)<=P0(2)))
                                        {
                                            //cout<<"intersez: "<<c2<<endl<<p<<endl;
                                            c2=c2+1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                if(c1==2 && c2==2){
                    cout<<"P "<<P.transpose()<<endl<<"t "<<t.transpose()<<endl;
                }

            }
        }



    }

    return t;
}

vector<double>CalcoloTracce(vector<vector<double>> &IDs,int n,vector<double> &FractureId,vector<double> &NumVertices,vector<MatrixXd> &ListVertices,vector<MatrixXd> &ListCord){
    vector<double> TraceId;
    Vector3d P0;
    Vector3d P1;
    Vector3d P2;
    Vector3d P3;
    MatrixXd A(3,3);
    Vector3d b;
    MatrixXd cord(3,2);
    Vector3d P;// punto piano,retta
    Vector3d t;// direzione retta
    for (int i = 0; i < n-1; ++i) {
        for (int j = i+1; j < n; ++j) {
            t=SoluzioneSistema(P,ListVertices,P0,P1,P2,P3,A,b,i,j,NumVertices);


        }

    }

    return TraceId;

}









