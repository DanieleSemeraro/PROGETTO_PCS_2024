#include "utilis.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include"Eigen/Eigen"
#include <iomanip>
#include <cmath>
#include <algorithm>

using namespace std;
using namespace Eigen;

namespace DFNLibrary{


void ImportDFN(const string &filename,int n,Fractures& fractures)
{
    ifstream fin(filename);
    string line;
    int c=0;//serve come contatore per scegliere cosa memorizzare e dove
    vector<double> p;//serve a memorizzare la posizione del punto e virgola
    int a=0;//altro contatore utile alla memorizzazione dei vertici
    MatrixXd Vertices;

    while (getline(fin,line )) {
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
            fractures.FractureId.push_back(stoi(line.substr(0,p[0])));
            fractures.NumVertices.push_back(stoi(line.substr(p[0]+1)));
            c=0;
            p.clear();
            Vertices.resize(3,fractures.NumVertices[fractures.NumVertices.size()-1]);
        }
        else if(c==2){
            for (unsigned i = 0; i < line.size(); i++) {
                if(line[i] == ';'){
                    p.push_back(i);
                }
            }
            Vertices(a,0)=stod(line.substr(0,p[0]));
            for (int j = 1; j < fractures.NumVertices.back()-1; j++) {
                Vertices(a,j)=stod(line.substr(p[j-1]+1,p[j]));
            }

            int b=fractures.NumVertices.back()-1;
            Vertices(a,b)=stod(line.substr(p[b-1]+1));

            a=a+1;
            p.clear();
            if(a==3){
                fractures.ListVertices.push_back(Vertices);
                a=0;
                c=0;
            }
        }
    }

    fin.close();
}

void Fractures::sfere(int n,vector<Vector2i> &fratturescluse){
    vector<Vector3d> baricentri;//contine i baricentri delle fratture
    int numver=0;//numero vertici frattura
    Vector3d A;//vertici
    Vector3d B;
    Vector3d C;
    double area;//area dei triangoli
    Vector3d AB;//lati triangoli per calcolo area
    Vector3d AC;
    Vector3d centroidtriangle;//baricentro triangolo
    double areatotale=0;
    Vector3d bari;//baricentro della frattura
    bari<<0,0,0;
    vector<double> raggi;//contiene i raggi di tutte le fratture
    double R=0;//raggio singolo
    double distanza;//distanza baricentro vertice
    Vector2i id;
    for (int i = 0; i < n; ++i) {
        numver=NumVertices[i];
        A=ListVertices[i].col(0);
        for (int j = 1; j < numver-1; ++j) {
            B=ListVertices[i].col(j);
            C=ListVertices[i].col(j+1);
            AB<<B[0]-A[0],B[1]-A[1],B[2]-A[2];
            AC<<C[0]-A[0],C[1]-A[1],C[2]-A[2];
            area=0.5*((AB.cross(AC)).norm());
            centroidtriangle[0]=(A[0]+B[0]+C[0])/3.0;
            centroidtriangle[1]=(A[1]+B[1]+C[1])/3.0;
            centroidtriangle[2]=(A[2]+B[2]+C[2])/3.0;
            bari[0]+=area*centroidtriangle[0];
            bari[1]+=area*centroidtriangle[1];
            bari[2]+=area*centroidtriangle[2];
            areatotale+=area;
        }
        bari[0]/=areatotale;
        bari[1]/=areatotale;
        bari[2]/=areatotale;

        baricentri.push_back(bari);

        areatotale=0;
        bari[0]=0;
        bari[1]=0;
        bari[2]=0;

        //cout<<"bari n: "<<i<<" "<<scientific<<setprecision(16)<<baricentri[i].transpose()<<endl;

        for (int j = 0; j < numver; ++j) {
            A=ListVertices[i].col(j);
            distanza=sqrt((A[0]-baricentri[i](0))*(A[0]-baricentri[i](0)) + (A[1]-baricentri[i](1))*(A[1]-baricentri[i](1)) + (A[2]-baricentri[i](2))*(A[2]-baricentri[i](2)));
            if(distanza>R){
                R=distanza;
            }
        }
        raggi.push_back(R);
        R=0;
        //cout<<"raggio n: "<<i<<" "<<scientific<<setprecision(16)<<raggi[i]<<endl;
    }

    for (int i = 0; i < n-1; ++i) {
        for (int j = i+1; j < n; ++j) {
            distanza=sqrt((baricentri[i](0)-baricentri[j](0))*(baricentri[i](0)-baricentri[j](0)) + (baricentri[i](1)-baricentri[j](1))*(baricentri[i](1)-baricentri[j](1)) + (baricentri[i](2)-baricentri[j](2))*(baricentri[i](2)-baricentri[j](2)));
            if(distanza > (raggi[i]+raggi[j])){
                id[0]=FractureId[i];
                id[1]=FractureId[j];
                fratturescluse.push_back(id);
                //cout<<"non si intersecano le sfere delle fratture "<<scientific<<setprecision(16)<<id<<endl;
            }
            else{
                //cout<<"si intersecano le sfere delle fratture "<<scientific<<setprecision(16)<<i<<" "<<j<<endl;
            }
        }
    }
}

void Traces::CalcoloDirezioneTracce(int &NumberOfTraces,Fractures& fractures,int n,vector<Vector2i> &fratturescluse){//serve per ottenere la retta su cui somo presenti le tracce, restituisce punto e direzione di ogni retta utile con tracce
    Vector3d u;//vettore 1
    Vector3d P0;//primo punto del piano
    Vector3d P1;//secondo punto del piano
    Vector3d P2;//terzo punto del piano
    MatrixXd A(3,3);//matrice del sistema lineare per trovare il punto
    Vector3d b;//vettore del sistema con matrice A. usato anche per il secondo sist linear
    Vector3d P;//punto soluzione del sistema A*P=b.
    Vector3d v;//vettore 2
    Vector2i vet;//var usata per memorizzare i due id delle fratture con traccia
    MatrixXd mat(3,2);//var usata per memor il punto e la direzione della retta
    Vector3d n1;//direzione piano 1 normalizzata
    Vector3d n2;//direzione piano 2 normalizzata
    Vector3d t;//direzione retta di intersezione
    double d1;//parametro d del primo piano. equazione del piano ax+by+cz=d
    double d2;//parametro d del secondo piano. equazione del piano ax+by+cz=d
    MatrixXd alpha(3,2);// matrice per il sistema lineare dell intersezione tra due rette
    Vector2d sol;//soluzione del secondo sist lineare
    Vector3d p;// intersezione tra retta della traccia e retta tra vertici
    Vector3d sis2;//variabile creata per verificare che effettivamente alpha*p sia uguale a b.
    Vector3d sis1;//variabile creata per verificare che effettivamente A*P sia uguale a b.
    //MatrixXd punti(3,2);
    int c1=0;//contatore
    int c2=0;//contatore
    double tol = 1e-8; // Tolleranza di 10^-8
    Vector2i id;//trova coppie di id da verificare
    int cont=0;//contatore per capire se una coppia di fratture si trova o meno in quelle d aescludere

    for (int i = 0; i < n-1; ++i) {//doppio ciclo per confrontare ogni rettangolo con gli altri per trovare eventuai tracce
        for (int j = i+1; j < n; ++j) {
            id[0]=fractures.FractureId[i];
            id[1]=fractures.FractureId[j];
            for (const auto& item : fratturescluse) {
                if(id == item){
                    cont=1;
                }
            }

            if(cont==1){
                cont=0;
                continue;
            }
            else{
                cont=0;
                P0=fractures.ListVertices[i].col(0);
                P1=fractures.ListVertices[i].col(1);
                P2=fractures.ListVertices[i].col(2);
                u=P2-P0;
                v=P1-P0;
                n1=(u.cross(v)).normalized();
                d1=n1.dot(P0);

                P0=fractures.ListVertices[j].col(0);
                P1=fractures.ListVertices[j].col(1);
                P2=fractures.ListVertices[j].col(2);
                u=P2-P0;
                v=P1-P0;
                n2=(u.cross(v)).normalized();
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
                        if((sis1(0)<=b(0)+tol && sis1(0)>=b(0)-tol) && (sis1(1)<=b(1)+tol && sis1(1)>=b(1)-tol) && (sis1(2)<=b(2)+tol && sis1(2)>=b(2)-tol)){// verifico che sis1 sia uguale a b rispetto una tolleranza tol
                            for (int k = 0; k < fractures.NumVertices[i]; ++k) {//ora devo vedere se la retta trovata sia effetivamente una frattura oppure no. devo trovare le intersezioni tra la retta e la prima frattura e verificare che siano comprese tra i vertici del bordo
                                if(k==3){
                                    alpha.col(0)=t;
                                    alpha.col(1)=fractures.ListVertices[i].col(0)-fractures.ListVertices[i].col(k);
                                    P0=fractures.ListVertices[i].col(k);
                                    P1=fractures.ListVertices[i].col(0);
                                    b=fractures.ListVertices[i].col(k)-P;
                                }
                                else{
                                    alpha.col(0)=t;
                                    alpha.col(1)=fractures.ListVertices[i].col(k+1)-fractures.ListVertices[i].col(k);
                                    P0=fractures.ListVertices[i].col(k);
                                    P1=fractures.ListVertices[i].col(k+1);
                                    b=fractures.ListVertices[i].col(k)-P;
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
                                                    //con questi if verifico che il punto p sia compreso tra due vertici
                                                    c1=c1+1;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            if(c1==2){//nel caso affermativo che la retta intersechi due volte la prima frattura verifico che lo faccia anche con la seconda, in caso contrario ho risparmiato dei calcoli
                                for (int k = 0; k < fractures.NumVertices[j]; ++k) {
                                    if(k==3){
                                        alpha.col(0)=t;
                                        alpha.col(1)=fractures.ListVertices[j].col(0)-fractures.ListVertices[j].col(k);
                                        P0=fractures.ListVertices[j].col(k);
                                        P1=fractures.ListVertices[j].col(0);
                                        b=fractures.ListVertices[j].col(k)-P;

                                    }
                                    else{
                                        alpha.col(0)=t;
                                        alpha.col(1)=fractures.ListVertices[j].col(k+1)-fractures.ListVertices[j].col(k);
                                        P0=fractures.ListVertices[j].col(k);
                                        P1=fractures.ListVertices[j].col(k+1);
                                        b=fractures.ListVertices[j].col(k)-P;
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

                            if(c1==2 && c2==2){//nel caso in cui la retta intersechi due volte sia un rettang che l altro vuol dire che è quella di un frattura fra essi quindi procedo a memorizzarmi le info utili

                                vet<<fractures.FractureId[i],fractures.FractureId[j];
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
    }



}
void Traces::CalcoloEstremi(int &NumberOfTraces,Fractures &fractures){//funzione per trovare gli stremi delle fratture

    Vector3d P;//punto della retta i prima
    Vector3d t;//direzione della retta
    Vector3d P0;//primo vertice
    Vector3d P1;//secondo vertice
    MatrixXd a(3,2);//mat del sistema lineare di intersez tra due rette
    Vector3d b;// vet del sist lineare
    Vector2d sol;//soluzione del sistema lineare
    Vector3d sis;//var per il confronto che a*sol=b
    Vector3d p;//punto di intersezione
    double tol=0.00000001;// solita tol
    MatrixXd intersez(3,4);//mat che memorizza le 4 intersezioni (due del primo rettangolo e due del secondo rettangolo)
    MatrixXd estremi(3,2);//mat dei due estremi della frattura
    int c=0;//contatore che permette di memorizzare le intersez nella matrice intersez al posto giusto
    int c1=0;// contatore che aiuta a memorizzare gli estremi della frattura nella mat estremi

    ofstream Outfile("Foglio1.txt");//inizio a stampare sul foglio
    Outfile<<"# Number of traces"<<endl;
    Outfile<<NumberOfTraces<<endl;
    Outfile<<"# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2"<<endl;

    for (int k = 0; k < NumberOfTraces; ++k) {//per ogni traccia devo trovare gli estremi
        for (int i = 0; i < fractures.NumVertices[IDs[k](0)]; ++i) {//IDs[k](0) significa che k sarebbe il traceID , quindi nella traccia numero k so che è data da una frattura a e una frattura b e ora trovo le intersez della a (per questo lo 0) dopo mettero 1
            if(i==3){
                P=ListCord[k].col(0);
                t=ListCord[k].col(1);
                P0=fractures.ListVertices[IDs[k](0)].col(i);
                P1=fractures.ListVertices[IDs[k](0)].col(0);
            }
            else{
                P=ListCord[k].col(0);
                t=ListCord[k].col(1);
                P0=fractures.ListVertices[IDs[k](0)].col(i);
                P1=fractures.ListVertices[IDs[k](0)].col(i+1);
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
                                //solita verifica che il punto sia compreso negli estremi
                                intersez.col(c)=p;
                                c=c+1;
                            }
                        }
                    }

                }
            }
        }
        for (int j = 0; j < fractures.NumVertices[IDs[k](1)]; ++j) {//verifico le stesse cose ma con il secondo rettangolo impiegato nella traccia
            if(j==3){
                P=ListCord[k].col(0);
                t=ListCord[k].col(1);
                P0=fractures.ListVertices[IDs[k](1)].col(j);
                P1=fractures.ListVertices[IDs[k](1)].col(0);
            }
            else{
                P=ListCord[k].col(0);
                t=ListCord[k].col(1);
                P0=fractures.ListVertices[IDs[k](1)].col(j);
                P1=fractures.ListVertices[IDs[k](1)].col(j+1);
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
        c=0;
        for (int z = 0; z < 4; ++z) {//ora trovo finalmente i due estremi della frattura, immagino di sovrappore due segmenti e di trovare gli estremi della loro intersezione
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
    for (int z = 0; z < NumberOfTraces; ++z) {// stampo sul foglio le rispettive informazioni
        Outfile<<z<<"; "<<IDs[z](0)<<"; "<<IDs[z](1)<<"; "<<scientific<<setprecision(16)<<cordinate[z].col(0).transpose()<<"; "<<scientific<<setprecision(16)<<cordinate[z].col(1).transpose()<<endl;
    }

    Outfile.close();



}


void Traces::Ordinamento(Fractures& fractures){//ultima funzione che permette di calcolare il numero di tracce presenti su ogni frattura, la loro lunghezza e ordinarle in maniera decrescente

    int Ntraces=0;//num tracce totali di ogni singola frattura
    bool Tips ;// ottiene valore true se è non passante, false se è passante
    double length;//lunghezza traccia
    vector<double> lung;//vettore contenente le lunghezzze delle tracce
    double l;//var di supporto lunghezza
    double tol = 1e-8; // Tolleranza di 10^-8
    vector<double>  t;//memorizza i traceid , serve per evitare di stampare più volte la stessa riga

    ofstream Outfile("Foglio2.txt");

    for (int k = 0; k < static_cast<int>(fractures.FractureId.size()); ++k) {// per ogni frattura devo visualizzare delle cose
        for (unsigned int i = 0; i < IDs.size(); ++i) {//calcolo il num di tracce  totali su di essa
            if(k==IDs[i](0) || k==IDs[i](1)){
                Ntraces=Ntraces+1;
            }
        }
        Outfile<<"# FractureId; NumTraces"<<endl;//inizio a stampare
        Outfile<<k<<"; "<<Ntraces<<endl;
        if(Ntraces!=0){
            Outfile<<"# TraceId; Tips; Length"<<endl;// non la visualizzo sul foglio se non ha tracce
        }
        Ntraces=0;
        for (unsigned int i = 0; i < IDs.size(); ++i) {
            if(k==IDs[i](0) || k==IDs[i](1)){
                length=sqrt( (cordinate[i](0,0)-cordinate[i](0,1))*(cordinate[i](0,0)-cordinate[i](0,1)) + (cordinate[i](1,0)-cordinate[i](1,1))*(cordinate[i](1,0)-cordinate[i](1,1)) + (cordinate[i](2,0)-cordinate[i](2,1))*(cordinate[i](2,0)-cordinate[i](2,1)) );
                lung.push_back(length);
            }
        }
        //le ordino in maniera decrescente con l'algoritmo spiegato in classe bubblesort
        BubbleSort(lung);
        for (unsigned int j = 0; j < lung.size(); ++j) {//adesso devo stampare queste lunghezze nel nuovo ordine ottenuto ma con le loro rispettive informazioni
            l=lung[j];
            for (unsigned int i = 0; i < IDs.size(); ++i) {//ripete finchè non trova la lunghezza corretta
                auto it=find(t.begin(),t.end(),i);//serve per evitare di stampare più volte le stesse righe nel caso siano presenti lunghezze uguali
                if(k==IDs[i](0) || k==IDs[i](1)){
                    CalcoloPassante(k,i,fractures);
                    if(pass==0){
                        Tips=false;
                        length=sqrt( (cordinate[i](0,0)-cordinate[i](0,1))*(cordinate[i](0,0)-cordinate[i](0,1)) + (cordinate[i](1,0)-cordinate[i](1,1))*(cordinate[i](1,0)-cordinate[i](1,1)) + (cordinate[i](2,0)-cordinate[i](2,1))*(cordinate[i](2,0)-cordinate[i](2,1)) );
                        if(l-tol<=length && length<=l+tol){
                            if(it==t.end()){
                                t.push_back(i);//memorizza gli traceID sempre per il fatto di non stampare più volte le stesse righe
                                Outfile<<i<<"; "<<Tips<<"; "<<scientific<<setprecision(16)<<length<<endl;
                            }
                        }
                    }
                    else{
                        Tips=true;
                        length=sqrt( (cordinate[i](0,0)-cordinate[i](0,1))*(cordinate[i](0,0)-cordinate[i](0,1)) + (cordinate[i](1,0)-cordinate[i](1,1))*(cordinate[i](1,0)-cordinate[i](1,1)) + (cordinate[i](2,0)-cordinate[i](2,1))*(cordinate[i](2,0)-cordinate[i](2,1)) );
                        if(l-tol<=length && length<=l+tol){
                            if(it==t.end()){
                                t.push_back(i);//memorizza gli traceID sempre per il fatto di non stampare più volte le stesse righe
                                Outfile<<i<<"; "<<Tips<<"; "<<scientific<<setprecision(16)<<length<<endl;
                            }
                        }
                    }

                }
            }

        }

        t.clear();
        lung.clear();
    }

    Outfile.close();
}

void Traces::CalcoloPassante(const int k,const int i,Fractures& fractures){
    int numvertici=fractures.NumVertices[k];
    int c=0;// contatore
    double tol = 1e-8;
    Vector3d AB;
    Vector3d AC;
    Vector3d BC;
    double crossproduct=0;
    for (int j = 0; j < numvertici; ++j) {//verifico se il primo estremo della traccia sia sul bordo della frattura
        AB=fractures.ListVertices[k].col((j+1)%numvertici)-fractures.ListVertices[k].col(j);
        AC=cordinate[i].col(0)-fractures.ListVertices[k].col(j);
        BC=cordinate[i].col(0)-fractures.ListVertices[k].col(j+1);
        crossproduct=AB.cross(AC).norm();
        if(crossproduct > tol){//verifica che i punti siano in linea
            continue;
        }
        else if(AB.dot(AC)<0 || AB.dot(AC) > AB.dot(AB)){// C si trova tra A e B se 0 <= dot(AB, AC) <= dot(AB, AB)
            continue;
        }
        else{
            c=c+1;
        }
    }

    if(c==1){//nel caso in cui il primo estremo della traccia sia sul bordo proseguo e faccio lo stesso con il secondo
        for (int j = 0; j < numvertici; ++j) {
            AB=fractures.ListVertices[k].col((j+1)%numvertici)-fractures.ListVertices[k].col(j);
            AC=cordinate[i].col(1)-fractures.ListVertices[k].col(j);
            BC=cordinate[i].col(1)-fractures.ListVertices[k].col((j+1)%numvertici);
            crossproduct=AB.cross(AC).norm();
            if(crossproduct > tol){
                continue;
            }
            else if(AB.dot(AC)<0 || AB.dot(AC) > AB.dot(AB)){
                continue;
            }
            else{
                c=c+1;
            }
        }
    }

    if(c==2){//se entrambi sono sul bordo è passante
        pass=0;//traccia passante

    }
    else{
        pass=1;//traccia non passante
    }
}

void BubbleSort(vector<double>& data)//algoritmo bubblesort spiegato in classe, dato un vet in input lo restituisce riordinato
{
    size_t rem_size = data.size();
    size_t last_seen = rem_size;
    bool swapped = true;

    while (swapped) {
        swapped = false;
        for (size_t i = 1; i < rem_size; i++) {
            if (data[i-1] < data[i]) {
                swap(data[i-1], data[i]);
                swapped = true;
                last_seen = i;
            }
        }
        rem_size = last_seen;
    }
}

}

ostream& operator<<(ostream& os, const vector<int> a)
{
    for (size_t i = 0; i < a.size(); i++) {
        os<<scientific<<setprecision(16)<< a[i]<< " ";
    }

    return os;
}

ostream& operator<<(ostream& os, const vector<MatrixXd> a)
{
    for (size_t i = 0; i < a.size(); i++) {
        os<<"Id frattura: "<<i<<endl<<"Matrice vertici: "<<endl<<scientific<<setprecision(16)<< a[i]<<endl;
        cout<<" "<<endl;
    }
    return os;
}
