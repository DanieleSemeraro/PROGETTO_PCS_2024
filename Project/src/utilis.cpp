#include "utilis.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include"Eigen/Eigen"
#include <iomanip>
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
    double tol=0.00000001;

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
                    if((sis1(0)<=b(0)+tol && sis1(0)>=b(0)-tol) && (sis1(1)<=b(1)+tol && sis1(1)>=b(1)-tol) && (sis1(2)<=b(2)+tol && sis1(2)>=b(2)-tol)){
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
                                if((sis2(0)<=b(0)+tol && sis2(0)>=b(0)-tol) && (sis2(1)<=b(1)+tol && sis2(1)>=b(1)-tol) && (sis2(2)<=b(2)+tol && sis2(2)>=b(2)-tol)){
                                    p=P+sol(0)*t;

                                    if(((P0(0)<p(0)||(p(0)<=P0(0)+tol && p(0)>=P0(0)-tol)) && (p(0)<P1(0)||(p(0)<=P1(0)+tol && p(0)>=P1(0)-tol))) || ((P0(0)>p(0)||(p(0)<=P0(0)+tol && p(0)>=P0(0)-tol)) && (p(0)>P1(0)||(p(0)<=P1(0)+tol && p(0)>=P1(0)-tol))))
                                    {
                                        if(((P0(1)<p(1)||(p(1)<=P0(1)+tol && p(1)>=P0(1)-tol)) && (p(1)<P1(1)||(p(1)<=P1(1)+tol && p(1)>=P1(1)-tol))) || ((P0(1)>p(1)||(p(1)<=P0(1)+tol && p(1)>=P0(1)-tol)) && (p(1)>P1(1)||(p(1)<=P1(1)+tol && p(1)>=P1(1)-tol))))
                                        {
                                            if(((P0(2)<p(2)||(p(2)<=P0(2)+tol && p(2)>=P0(2)-tol)) && (p(2)<P1(2)||(p(2)<=P1(2)+tol && p(2)>=P1(2)-tol))) || ((P0(2)>p(2)||(p(2)<=P0(2)+tol && p(2)>=P0(2)-tol)) && (p(2)>P1(2)||(p(2)<=P1(2)+tol && p(2)>=P1(2)-tol))))
                                            {

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
                                    if((sis2(0)<=b(0)+tol && sis2(0)>=b(0)-tol) && (sis2(1)<=b(1)+tol && sis2(1)>=b(1)-tol) && (sis2(2)<=b(2)+tol && sis2(2)>=b(2)-tol)){
                                        p=P+sol(0)*t;

                                        if(((P0(0)<p(0)||(p(0)<=P0(0)+tol && p(0)>=P0(0)-tol)) && (p(0)<P1(0)||(p(0)<=P1(0)+tol && p(0)>=P1(0)-tol))) || ((P0(0)>p(0)||(p(0)<=P0(0)+tol && p(0)>=P0(0)-tol)) && (p(0)>P1(0)||(p(0)<=P1(0)+tol && p(0)>=P1(0)-tol))))
                                        {
                                            if(((P0(1)<p(1)||(p(1)<=P0(1)+tol && p(1)>=P0(1)-tol)) && (p(1)<P1(1)||(p(1)<=P1(1)+tol && p(1)>=P1(1)-tol))) || ((P0(1)>p(1)||(p(1)<=P0(1)+tol && p(1)>=P0(1)-tol)) && (p(1)>P1(1)||(p(1)<=P1(1)+tol && p(1)>=P1(1)-tol))))
                                            {
                                                if(((P0(2)<p(2)||(p(2)<=P0(2)+tol && p(2)>=P0(2)-tol)) && (p(2)<P1(2)||(p(2)<=P1(2)+tol && p(2)>=P1(2)-tol))) || ((P0(2)>p(2)||(p(2)<=P0(2)+tol && p(2)>=P0(2)-tol)) && (p(2)>P1(2)||(p(2)<=P1(2)+tol && p(2)>=P1(2)-tol))))
                                                {

                                                    c2=c2+1;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        if(c1==2 && c2==2){

                            vet<<FractureId[i],FractureId[j];
                            IDs.push_back(vet);
                            mat.col(0)=P;
                            mat.col(1)=t;
                            ListCord.push_back(mat);
                            NumberOfTraces=NumberOfTraces+1;

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

vector<MatrixXd> CalcoloEstremi(int &NumberOfTraces,vector<VectorXd> &IDs,vector<double> &NumVertices,vector<MatrixXd> &ListVertices,vector<MatrixXd> &ListCord,vector<double> &pass){
    Vector3d P;
    Vector3d t;
    Vector3d P0;
    Vector3d P1;
    MatrixXd a(3,2);
    Vector3d b;
    Vector2d sol;
    Vector3d sis;
    Vector3d p;
    double tol=0.00000001;
    MatrixXd intersez(3,4);
    MatrixXd estremi(3,2);
    int c=0;
    int c1=0;
    vector<MatrixXd> cordinate;


    ofstream Outfile("Foglio1.txt");
    Outfile<<"# Number of traces"<<endl;
    Outfile<<NumberOfTraces<<endl;
    Outfile<<"# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2"<<endl;

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
                if((sis(0)<=b(0)+tol && sis(0)>=b(0)-tol) && (sis(1)<=b(1)+tol && sis(1)>=b(1)-tol) && (sis(2)<=b(2)+tol && sis(2)>=b(2)-tol)){
                    p=P+sol(0)*t;
                    if(((P0(0)<p(0)||(p(0)<=P0(0)+tol && p(0)>=P0(0)-tol)) && (p(0)<P1(0)||(p(0)<=P1(0)+tol && p(0)>=P1(0)-tol))) || ((P0(0)>p(0)||(p(0)<=P0(0)+tol && p(0)>=P0(0)-tol)) && (p(0)>P1(0)||(p(0)<=P1(0)+tol && p(0)>=P1(0)-tol))))
                    {
                        if(((P0(1)<p(1)||(p(1)<=P0(1)+tol && p(1)>=P0(1)-tol)) && (p(1)<P1(1)||(p(1)<=P1(1)+tol && p(1)>=P1(1)-tol))) || ((P0(1)>p(1)||(p(1)<=P0(1)+tol && p(1)>=P0(1)-tol)) && (p(1)>P1(1)||(p(1)<=P1(1)+tol && p(1)>=P1(1)-tol))))
                        {
                            if(((P0(2)<p(2)||(p(2)<=P0(2)+tol && p(2)>=P0(2)-tol)) && (p(2)<P1(2)||(p(2)<=P1(2)+tol && p(2)>=P1(2)-tol))) || ((P0(2)>p(2)||(p(2)<=P0(2)+tol && p(2)>=P0(2)-tol)) && (p(2)>P1(2)||(p(2)<=P1(2)+tol && p(2)>=P1(2)-tol))))
                            {
                                intersez.col(c)=p;
                                c=c+1;
                            }
                        }
                    }

                }
            }
        }
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
                if((sis(0)<=b(0)+tol && sis(0)>=b(0)-tol) && (sis(1)<=b(1)+tol && sis(1)>=b(1)-tol) && (sis(2)<=b(2)+tol && sis(2)>=b(2)-tol)){
                    p=P+sol(0)*t;
                    if(((P0(0)<p(0)||(p(0)<=P0(0)+tol && p(0)>=P0(0)-tol)) && (p(0)<P1(0)||(p(0)<=P1(0)+tol && p(0)>=P1(0)-tol))) || ((P0(0)>p(0)||(p(0)<=P0(0)+tol && p(0)>=P0(0)-tol)) && (p(0)>P1(0)||(p(0)<=P1(0)+tol && p(0)>=P1(0)-tol))))
                    {
                        if(((P0(1)<p(1)||(p(1)<=P0(1)+tol && p(1)>=P0(1)-tol)) && (p(1)<P1(1)||(p(1)<=P1(1)+tol && p(1)>=P1(1)-tol))) || ((P0(1)>p(1)||(p(1)<=P0(1)+tol && p(1)>=P0(1)-tol)) && (p(1)>P1(1)||(p(1)<=P1(1)+tol && p(1)>=P1(1)-tol))))
                        {
                            if(((P0(2)<p(2)||(p(2)<=P0(2)+tol && p(2)>=P0(2)-tol)) && (p(2)<P1(2)||(p(2)<=P1(2)+tol && p(2)>=P1(2)-tol))) || ((P0(2)>p(2)||(p(2)<=P0(2)+tol && p(2)>=P0(2)-tol)) && (p(2)>P1(2)||(p(2)<=P1(2)+tol && p(2)>=P1(2)-tol))))
                            {
                                intersez.col(c)=p;
                                c=c+1;
                            }
                        }
                    }
                }
            }

        }
        if( ( ((intersez(0,0)<=intersez(0,2)+tol && intersez(0,0)>=intersez(0,2)-tol)&&(intersez(1,0)<=intersez(1,2)+tol && intersez(1,0)>=intersez(1,2)-tol)&&(intersez(2,0)<=intersez(2,2)+tol && intersez(2,0)>=intersez(2,2)-tol)) || ((intersez(0,0)<=intersez(0,3)+tol && intersez(0,0)>=intersez(0,3)-tol)&&(intersez(1,0)<=intersez(1,3)+tol && intersez(1,0)>=intersez(1,3)-tol)&&(intersez(2,0)<=intersez(2,3)+tol && intersez(2,0)>=intersez(2,3)-tol)) ) && ( ((intersez(0,1)<=intersez(0,2)+tol && intersez(0,1)>=intersez(0,2)-tol)&&(intersez(1,1)<=intersez(1,2)+tol && intersez(1,1)>=intersez(1,2)-tol)&&(intersez(2,1)<=intersez(2,2)+tol && intersez(2,1)>=intersez(2,2)-tol)) || ((intersez(0,1)<=intersez(0,3)+tol && intersez(0,1)>=intersez(0,3)-tol)&&(intersez(1,1)<=intersez(1,3)+tol && intersez(1,1)>=intersez(1,3)-tol)&&(intersez(2,1)<=intersez(2,3)+tol && intersez(2,1)>=intersez(2,3)-tol)) ) ){
            pass.push_back(0);
        }
        else{
            pass.push_back(1);
        }
        c=0;
        for (int z = 0; z < 4; ++z) {
            if((z==0 || z==1) && c1<2){
                if(((intersez(0,2)<intersez(0,z)||(intersez(0,z)<=intersez(0,2)+tol && intersez(0,z)>=intersez(0,2)-tol)) && (intersez(0,z)<intersez(0,3)||(intersez(0,z)<=intersez(0,3)+tol && intersez(0,z)>=intersez(0,3)-tol))) || ((intersez(0,2)>intersez(0,z)||(intersez(0,z)<=intersez(0,2)+tol && intersez(0,z)>=intersez(0,2)-tol)) && (intersez(0,z)>intersez(0,3)||(intersez(0,z)<=intersez(0,3)+tol && intersez(0,z)>=intersez(0,3)-tol))))
                {
                    if(((intersez(1,2)<intersez(1,z)||(intersez(1,z)<=intersez(1,2)+tol && intersez(1,z)>=intersez(1,2)-tol)) && (intersez(1,z)<intersez(1,3)||(intersez(1,z)<=intersez(1,3)+tol && intersez(1,z)>=intersez(1,3)-tol))) || ((intersez(1,2)>intersez(1,z)||(intersez(1,z)<=intersez(1,2)+tol && intersez(1,z)>=intersez(1,2)-tol)) && (intersez(1,z)>intersez(1,3)||(intersez(1,z)<=intersez(1,3)+tol && intersez(1,z)>=intersez(1,3)-tol))))
                    {
                        if(((intersez(2,2)<intersez(2,z)||(intersez(2,z)<=intersez(2,2)+tol && intersez(2,z)>=intersez(2,2)-tol)) && (intersez(2,z)<intersez(2,3)||(intersez(2,z)<=intersez(2,3)+tol && intersez(2,z)>=intersez(2,3)-tol))) || ((intersez(2,2)>intersez(2,z)||(intersez(2,z)<=intersez(2,2)+tol && intersez(2,z)>=intersez(2,2)-tol)) && (intersez(2,z)>intersez(2,3)||(intersez(2,z)<=intersez(2,3)+tol && intersez(2,z)>=intersez(2,3)-tol))))
                        {
                            estremi.col(c1)=intersez.col(z);
                            c1=c1+1;
                        }
                    }
                }
            }
            else if((z==2 || z==3) && c1<2){

                if(((intersez(0,0)<intersez(0,z)||(intersez(0,z)<=intersez(0,0)+tol && intersez(0,z)>=intersez(0,0)-tol)) && (intersez(0,z)<intersez(0,1)||(intersez(0,z)<=intersez(0,1)+tol && intersez(0,z)>=intersez(0,1)-tol))) || ((intersez(0,0)>intersez(0,z)||(intersez(0,z)<=intersez(0,0)+tol && intersez(0,z)>=intersez(0,0)-tol)) && (intersez(0,z)>intersez(0,1)||(intersez(0,z)<=intersez(0,1)+tol && intersez(0,z)>=intersez(0,1)-tol))))
                {
                    if(((intersez(1,0)<intersez(1,z)||(intersez(1,z)<=intersez(1,0)+tol && intersez(1,z)>=intersez(1,0)-tol)) && (intersez(1,z)<intersez(1,1)||(intersez(1,z)<=intersez(1,1)+tol && intersez(1,z)>=intersez(1,1)-tol))) || ((intersez(1,0)>intersez(1,z)||(intersez(1,z)<=intersez(1,0)+tol && intersez(1,z)>=intersez(1,0)-tol)) && (intersez(1,z)>intersez(1,1)||(intersez(1,z)<=intersez(1,1)+tol && intersez(1,z)>=intersez(1,1)-tol))))
                    {
                        if(((intersez(2,0)<intersez(2,z)||(intersez(2,z)<=intersez(2,0)+tol && intersez(2,z)>=intersez(2,0)-tol)) && (intersez(2,z)<intersez(2,1)||(intersez(2,z)<=intersez(2,1)+tol && intersez(2,z)>=intersez(2,1)-tol))) || ((intersez(2,0)>intersez(2,z)||(intersez(2,z)<=intersez(2,0)+tol && intersez(2,z)>=intersez(2,0)-tol)) && (intersez(2,z)>intersez(2,1)||(intersez(2,z)<=intersez(2,1)+tol && intersez(2,z)>=intersez(2,1)-tol))))
                        {
                            estremi.col(c1)=intersez.col(z);
                            c1=c1+1;
                        }
                    }
                }

            }

        }
        cordinate.push_back(estremi);
        c1=0;
    }
    for (int z = 0; z < NumberOfTraces; ++z) {
        Outfile<<z<<"; "<<IDs[z](0)<<"; "<<IDs[z](1)<<"; "<<cordinate[z].col(0).transpose()<<"; "<<cordinate[z].col(1).transpose()<<endl;
    }
    Outfile.close();

    return cordinate;

}

vector<MatrixXd> Ordinamento(vector<double> FractureId,int &NumberOfTraces,vector<VectorXd> &IDs,vector<double> &NumVertices,vector<MatrixXd> &ListVertices,vector<MatrixXd> &ListCord, vector<MatrixXd> &cordinate,vector<double> &pass){
    int Ntraces=0;
    bool Tips ;
    double length;

    ofstream Outfile("Foglio2.txt");

    for (unsigned int k = 0; k < FractureId.size(); ++k) {
        for (unsigned int i = 0; i < IDs.size(); ++i) {
            if(k==IDs[i](0) || k==IDs[i](1)){
                Ntraces=Ntraces+1;
            }
        }
        Outfile<<"# FractureId; NumTraces"<<endl;
        Outfile<<k<<"; "<<Ntraces<<endl;
        Ntraces=0;
        Outfile<<"# TraceId; Tips; Length"<<endl;
        for (unsigned int i = 0; i < IDs.size(); ++i) {
            if(k==IDs[i](0) || k==IDs[i](1)){
                if(pass[i]==1){
                    Tips=true;
                }
                else{
                    Tips=false;
                }
                length=sqrt( (cordinate[i](0,0)-cordinate[i](0,1))*(cordinate[i](0,0)-cordinate[i](0,1)) + (cordinate[i](1,0)-cordinate[i](1,1))*(cordinate[i](1,0)-cordinate[i](1,1)) + (cordinate[i](2,0)-cordinate[i](2,1))*(cordinate[i](2,0)-cordinate[i](2,1)) );
                Outfile<<i<<"; "<<Tips<<"; "<<setprecision(15)<<length<<endl; //Tips stampa 1(true) se la traccia è non passante e 0(false) se la traccia è passante
            }
        }
    }
    return cordinate;

}

void BubbleSort(vector<double>& data)
{
    size_t rem_size = data.size();
    size_t last_seen = rem_size;
    bool swapped = true;

    while (swapped) {
        swapped = false;
        for (size_t i = 1; i < rem_size; i++) {
            if (data[i-1] > data[i]) {
                swap(data[i-1], data[i]);
                swapped = true;
                last_seen = i;
            }
        }
        //        rem_size = rem_size - 1;
        rem_size = last_seen;
    }
}










