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

vector<VectorXd>CalcoloDirezioneTracce(int &NumberOfTraces,vector<VectorXd> &IDs,int n,vector<double> &FractureId,vector<double> &NumVertices,vector<MatrixXd> &ListVertices,vector<MatrixXd> &ListCord){
    Vector3d u;
    Vector3d P0;
    Vector3d P1;
    Vector3d P2;
    MatrixXd A(3,3);
    Vector3d b;
    Vector3d P;
    Vector3d v;
    Vector2d vet;
    MatrixXd mat(3,2);
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

    for (int i = 0; i < n-1; ++i) {
        for (int j = i+1; j < n; ++j) {
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
                                        //cout<<"point "<<p<<endl;
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
                            //cout<<"P "<<P.transpose()<<endl<<"t "<<t.transpose()<<endl;
                            vet<<FractureId[i],FractureId[j];
                            IDs.push_back(vet);
                            for (int z = 0; z < 3; ++z) {
                                P(z)=round(P(z)*pow(10,6))/pow(10,6);
                            }
                            for (int z = 0; z < 3; ++z) {
                                t(z)=round(t(z)*pow(10,6))/pow(10,6);
                            }
                            mat.col(0)=P;
                            mat.col(1)=t;
                            ListCord.push_back(mat);
                            NumberOfTraces=NumberOfTraces+1;
                            //cout<<"not"<<NumberOfTraces<<endl;
                        }
                        c1=0;
                        c2=0;

                    }
                }



            }



        }

    }

    return IDs;

}

vector<double> CalcoloEstremi(int &NumberOfTraces,vector<VectorXd> &IDs,vector<double> &NumVertices,vector<MatrixXd> &ListVertices,vector<MatrixXd> &ListCord){
    vector<vector<double>> estremi;
    Vector3d P;
    Vector3d t;
    Vector3d P0;
    Vector3d P1;
    MatrixXd a(3,2);
    Vector3d b;
    Vector2d sol;
    Vector3d sis;
    Vector3d p;
    Vector3d A;
    Vector3d B;
    Vector3d C;
    Vector3d D;
    vector<double> estr;
    int c=0;
    for (int k = 0; k < NumberOfTraces; ++k) {
        for (int i = 0; i < NumVertices[IDs[k](0)]; ++i) {
            if(i==3){
                P=ListCord[k].col(0);
                t=ListCord[k].col(1);
                P0=ListVertices[IDs[k](0)].col(i);
                P1=ListVertices[IDs[k](0)].col(0);
            }
            else{
                P=ListCord[k].col(0);
                t=ListCord[k].col(1);
                P0=ListVertices[IDs[k](0)].col(i);
                P1=ListVertices[IDs[k](0)].col(i+1);
            }
            a.col(0)=t;
            a.col(1)=P1-P0;
            b=P0-P;
            if((t.cross(P1-P0)).norm()>0){
                sol=a.colPivHouseholderQr().solve(b);
                sis=a*sol;
                for (int z = 0; z < 3; ++z) {
                    sis(z)=round(sis(z)*pow(10,6))/pow(10,6);
                    b(z)=round(b(z)*pow(10,6))/pow(10,6);
                }
                if((sis(0)<=b(0)+0.000001 && sis(0)>=b(0)-0.000001) && (sis(1)<=b(1)+0.000001 && sis(1)>=b(1)-0.000001) && (sis(2)<=b(2)+0.000001 && sis(2)>=b(2)-0.000001)){
                    p=P+sol(0)*t;
                    for (int z = 0; z < 3; ++z) {
                        p(z)=round(p(z)*pow(10,6))/pow(10,6);
                    }
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
                                //cout<<"sas "<<i<<endl;
                                //cout<<"punti di intersez "<<p<<endl;
                                A=ListVertices[IDs[k](1)].col(0);
                                B=ListVertices[IDs[k](1)].col(1);
                                C=ListVertices[IDs[k](1)].col(2);
                                D=ListVertices[IDs[k](1)].col(3);
                                for (int z = 0; z < 3; ++z) {
                                    A(z)=round(A(z)*pow(10,6))/pow(10,6);
                                    B(z)=round(B(z)*pow(10,6))/pow(10,6);
                                    C(z)=round(C(z)*pow(10,6))/pow(10,6);
                                    D(z)=round(D(z)*pow(10,6))/pow(10,6);
                                }
                                if(puntoInRettangolo(p,A,B,C,D)){
                                    for (int z = 0; z < 3; ++z) {
                                        estr.push_back(p(z));
                                    }
                                    estremi.push_back(estr);
                                    c=c+1;
                                    estr.clear();
                                }
                            }
                        }
                    }

                }
            }
        }
        if(c!=2){
            for (int j = 0; j < NumVertices[IDs[k](1)]; ++j) {
                if(j==3){
                    P=ListCord[k].col(0);
                    t=ListCord[k].col(1);
                    P0=ListVertices[IDs[k](1)].col(j);
                    P1=ListVertices[IDs[k](1)].col(0);
                }
                else{
                    P=ListCord[k].col(0);
                    t=ListCord[k].col(1);
                    P0=ListVertices[IDs[k](1)].col(j);
                    P1=ListVertices[IDs[k](1)].col(j+1);
                }
                a.col(0)=t;
                a.col(1)=P1-P0;
                b=P0-P;
                if((t.cross(P1-P0)).norm()>0){
                    sol=a.colPivHouseholderQr().solve(b);
                    sis=a*sol;
                    for (int z = 0; z < 3; ++z) {
                        sis(z)=round(sis(z)*pow(10,6))/pow(10,6);
                        b(z)=round(b(z)*pow(10,6))/pow(10,6);
                    }
                    if((sis(0)<=b(0)+0.000001 && sis(0)>=b(0)-0.000001) && (sis(1)<=b(1)+0.000001 && sis(1)>=b(1)-0.000001) && (sis(2)<=b(2)+0.000001 && sis(2)>=b(2)-0.000001)){
                        p=P+sol(0)*t;
                        for (int z = 0; z < 3; ++z) {
                            p(z)=round(p(z)*pow(10,6))/pow(10,6);
                        }
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
                                    //cout<<"six "<<j<<endl;
                                    //cout<<"punti di intersez "<<p<<endl;
                                    A=ListVertices[IDs[k](0)].col(0);
                                    B=ListVertices[IDs[k](0)].col(1);
                                    C=ListVertices[IDs[k](0)].col(2);
                                    D=ListVertices[IDs[k](0)].col(3);
                                    for (int z = 0; z < 3; ++z) {
                                        A(z)=round(A(z)*pow(10,6))/pow(10,6);
                                        B(z)=round(B(z)*pow(10,6))/pow(10,6);
                                        C(z)=round(C(z)*pow(10,6))/pow(10,6);
                                        D(z)=round(D(z)*pow(10,6))/pow(10,6);
                                    }
                                    if(puntoInRettangolo(p,A,B,C,D)){
                                        for (int z = 0; z < 3; ++z) {
                                            estr.push_back(p(z));
                                        }
                                        estremi.push_back(estr);
                                        c=c+1;
                                        estr.clear();
                                    }
                                }
                            }
                        }

                    }
                }

            }
        }
        cout<<"Per la traccia numero "<<k<<" abbiamo queste intersezioni:"<<endl;
        for (int z = 0; z < c; ++z) {
            cout<<estremi[z]<<endl;
        }
        c=0;
        estremi.clear();

    }

    return NumVertices;

}


bool puntoInRettangolo(Vector3d &p, Vector3d& A, Vector3d& B,Vector3d& C,Vector3d& D) {
    Vector3d AB = B - A;
    Vector3d AD = D - A;
    Vector3d AP = p - A;

    double dotABAP = AB.dot(AP);
    double dotABAB = AB.dot(AB);
    double dotADAP = AD.dot(AP);
    double dotADAD = AD.dot(AD);

    dotABAP=round(dotABAP*pow(10,6))/pow(10,6);
    dotABAB=round(dotABAB*pow(10,6))/pow(10,6);
    dotADAP=round(dotADAP*pow(10,6))/pow(10,6);
    dotADAD=round(dotADAD*pow(10,6))/pow(10,6);

    return 0 <= dotABAP && dotABAP <= dotABAB &&
           0 <= dotADAP && dotADAP <= dotADAD;
}






