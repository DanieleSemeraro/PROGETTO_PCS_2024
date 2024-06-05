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

ostream& operator<<(ostream& os, const vector<double> a)
{
    for (size_t i = 0; i < a.size(); i++) {
        os<<setprecision(15)<< a[i]<< " ";
    }

    return os;
}

ostream& operator<<(ostream& os, const vector<MatrixXd> a)
{
    for (size_t i = 0; i < a.size(); i++) {
        os<<"Id frattura: "<<i<<endl<<"Matrice vertici: "<<endl<<setprecision(15)<< a[i]<<endl;
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
    MatrixXd Vertices;

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
            Vertices.resize(3,NumVertices[NumVertices.size()-1]);
        }
        else if(c==2){
            for (unsigned i = 0; i < line.size(); i++) {
                if(line[i] == ';'){
                    p.push_back(i);
                }
            }
            Vertices(a,0)=stod(line.substr(0,p[0]));
            for (int j = 1; j < NumVertices.back()-1; j++) {
                Vertices(a,j)=stod(line.substr(p[j-1]+1,p[j]));
            }

            int b=NumVertices.back()-1;
            Vertices(a,b)=stod(line.substr(p[b-1]+1));

            a=a+1;
            p.clear();
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

vector<VectorXd>CalcoloDirezioneTracce(int &NumberOfTraces,vector<VectorXd> &IDs,int n,vector<double> &FractureId,vector<double> &NumVertices,vector<MatrixXd> &ListVertices,vector<MatrixXd> &ListCord){//serve per ottenere la retta su cui somo presenti le tracce, restituisce punto e direzione di ogni retta utile con tracce
    Vector3d u;//vettore 1
    Vector3d P0;//primo punto del piano
    Vector3d P1;//secondo punto del piano
    Vector3d P2;//terzo punto del piano
    MatrixXd A(3,3);//matrice del sistema lineare per trovare il punto
    Vector3d b;//vettore del sistema con matrice A. usato anche per il secondo sist linear
    Vector3d P;//punto soluzione del sistema A*P=b.
    Vector3d v;//vettore 2
    Vector2d vet;//var usata per memorizzare i due id delle fratture con traccia
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
    double tol=0.00000001;//tolleranza di 10^-8

    for (int i = 0; i < n-1; ++i) {//doppio ciclo per confrontare ogni rettangolo con gli altri per trovare eventuai tracce
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
                    if((sis1(0)<=b(0)+tol && sis1(0)>=b(0)-tol) && (sis1(1)<=b(1)+tol && sis1(1)>=b(1)-tol) && (sis1(2)<=b(2)+tol && sis1(2)>=b(2)-tol)){// verifico che sis1 sia uguale a b rispetto una tolleranza tol
                        for (int k = 0; k < NumVertices[i]; ++k) {//ora devo vedere se la retta trovata sia effetivamente una frattura oppure no. devo trovare le intersezioni tra la retta e la prima frattura e verificare che siano comprese tra i vertici del bordo
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
                                                //con questi if verifico che il punto p sia compreso tra due vertici
                                                c1=c1+1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        if(c1==2){//nel caso affermativo che la retta intersechi due volte la prima frattura verifico che lo faccia anche con la seconda, in caso contrario ho risparmiato dei calcoli
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

                        if(c1==2 && c2==2){//nel caso in cui la retta intersechi due volte sia un rettang che l altro vuol dire che è quella di un frattura fra essi quindi procedo a memorizzarmi le info utili

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

vector<MatrixXd> CalcoloEstremi(int &NumberOfTraces,vector<VectorXd> &IDs,vector<double> &NumVertices,vector<MatrixXd> &ListVertices,vector<MatrixXd> &ListCord,vector<double> &pass){//funzione per trovare gli stremi delle fratture
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
    vector<MatrixXd> cordinate;//vettore di matrici contenenti tutti gli estremi delle fratture

    ofstream Outfile("Foglio1.txt");//inizio a stampare sul foglio
    Outfile<<"# Number of traces"<<endl;
    Outfile<<NumberOfTraces<<endl;
    Outfile<<"# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2"<<endl;

    for (int k = 0; k < NumberOfTraces; ++k) {//per ogni traccia devo trovare gli estremi
        for (int i = 0; i < NumVertices[IDs[k](0)]; ++i) {//IDs[k](0) significa che k sarebbe il traceID , quindi nella traccia numero k so che è data da una frattura a e una frattura b e ora trovo le intersez della a (per questo lo 0) dopo mettero 1
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
                                //solita verifica che il punto sia compreso negli estremi
                                intersez.col(c)=p;
                                c=c+1;
                            }
                        }
                    }

                }
            }
        }
        for (int j = 0; j < NumVertices[IDs[k](1)]; ++j) {//verifico le stesse cose ma con il secondo rettangolo impiegato nella traccia
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
        //questo servirà poi dopo, comunque calcola se la traccia è passante o non passante.
        if( ( ((intersez(0,0)<=intersez(0,2)+tol && intersez(0,0)>=intersez(0,2)-tol)&&(intersez(1,0)<=intersez(1,2)+tol && intersez(1,0)>=intersez(1,2)-tol)&&(intersez(2,0)<=intersez(2,2)+tol && intersez(2,0)>=intersez(2,2)-tol)) || ((intersez(0,0)<=intersez(0,3)+tol && intersez(0,0)>=intersez(0,3)-tol)&&(intersez(1,0)<=intersez(1,3)+tol && intersez(1,0)>=intersez(1,3)-tol)&&(intersez(2,0)<=intersez(2,3)+tol && intersez(2,0)>=intersez(2,3)-tol)) ) && ( ((intersez(0,1)<=intersez(0,2)+tol && intersez(0,1)>=intersez(0,2)-tol)&&(intersez(1,1)<=intersez(1,2)+tol && intersez(1,1)>=intersez(1,2)-tol)&&(intersez(2,1)<=intersez(2,2)+tol && intersez(2,1)>=intersez(2,2)-tol)) || ((intersez(0,1)<=intersez(0,3)+tol && intersez(0,1)>=intersez(0,3)-tol)&&(intersez(1,1)<=intersez(1,3)+tol && intersez(1,1)>=intersez(1,3)-tol)&&(intersez(2,1)<=intersez(2,3)+tol && intersez(2,1)>=intersez(2,3)-tol)) ) ){
            pass.push_back(0);
        }
        else{
            pass.push_back(1);
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
        Outfile<<z<<"; "<<IDs[z](0)<<"; "<<IDs[z](1)<<"; "<<setprecision(15)<<cordinate[z].col(0).transpose()<<"; "<<setprecision(15)<<cordinate[z].col(1).transpose()<<endl;
    }
    Outfile.close();

    return cordinate;

}

vector<MatrixXd> Ordinamento(vector<double> FractureId,vector<VectorXd> &IDs, vector<MatrixXd> &cordinate,vector<double> &pass){//ultima funzione che permette di calcolare il numero di tracce presenti su ogni frattura, la loro lunghezza e ordinarle in maniera decrescente
    int Ntraces=0;//num tracce totali di ogni singola frattura
    bool Tips ;// ottiene valore true se è non passante, false se è passante
    double length;//lunghezza traccia
    vector<double> lungP;//vettore contenente le lunghezzze delle traccie passanti
    vector<double> lungNP;//vettore contenente le lunghezzze delle traccie non passanti
    double l;//var di supporto lunghezza
    double tol=0.00000001;// solita tol
    vector<double>  t;//memorizza i traceid , serve per evitare di stampare più volte la stessa riga

    ofstream Outfile("Foglio2.txt");

    for (unsigned int k = 0; k < FractureId.size(); ++k) {// per ogni frattura devo visualizzare delle cose
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
                if(pass[i]==1){//memorizzo tutte le lunghezze delle tracce separandole in pass e non pass
                    length=sqrt( (cordinate[i](0,0)-cordinate[i](0,1))*(cordinate[i](0,0)-cordinate[i](0,1)) + (cordinate[i](1,0)-cordinate[i](1,1))*(cordinate[i](1,0)-cordinate[i](1,1)) + (cordinate[i](2,0)-cordinate[i](2,1))*(cordinate[i](2,0)-cordinate[i](2,1)) );
                    lungNP.push_back(length);
                }
                else{
                    length=sqrt( (cordinate[i](0,0)-cordinate[i](0,1))*(cordinate[i](0,0)-cordinate[i](0,1)) + (cordinate[i](1,0)-cordinate[i](1,1))*(cordinate[i](1,0)-cordinate[i](1,1)) + (cordinate[i](2,0)-cordinate[i](2,1))*(cordinate[i](2,0)-cordinate[i](2,1)) );
                    lungP.push_back(length);
                }
            }
        }
        BubbleSort(lungNP);//le ordino in maniera decrescente con l'algoritmo spiegato in classe bubblesort
        BubbleSort(lungP);
        for (unsigned int j = 0; j < lungP.size(); ++j) {//adesso devo stampare queste lunghezze nel nuovo ordine ottenuto ma con le loro rispettive informazioni
            l=lungP[j];
            for (unsigned int i = 0; i < IDs.size(); ++i) {//ripete finchè non trova la lunghezza corretta (siamo nel caso di tracce passanti)
                auto it=find(t.begin(),t.end(),i);//serve per evitare di stampare più volte le stesse righe nel caso siano presenti lunghezze uguali
                if(k==IDs[i](0) || k==IDs[i](1)){
                    if(pass[i]==0){
                        Tips=false;
                        length=sqrt( (cordinate[i](0,0)-cordinate[i](0,1))*(cordinate[i](0,0)-cordinate[i](0,1)) + (cordinate[i](1,0)-cordinate[i](1,1))*(cordinate[i](1,0)-cordinate[i](1,1)) + (cordinate[i](2,0)-cordinate[i](2,1))*(cordinate[i](2,0)-cordinate[i](2,1)) );
                        if(l-tol<=length && length<=l+tol){
                            if(it==t.end()){
                                t.push_back(i);//memorizza gli traceID sempre per il fatto di non stampare più volte le stesse righe
                                Outfile<<i<<"; "<<Tips<<"; "<<setprecision(15)<<length<<endl;
                            }
                        }
                    }

                }
            }

        }

        t.clear();

        for (unsigned int j = 0; j < lungNP.size(); ++j) {// esegue uguale a sopra, solo che lo fa con le tracce non passanti
            l=lungNP[j];
            for (unsigned int i = 0; i < IDs.size(); ++i) {
                auto it=find(t.begin(),t.end(),i);
                if(k==IDs[i](0) || k==IDs[i](1)){
                    if(pass[i]==1){
                        Tips=true;
                        length=sqrt( (cordinate[i](0,0)-cordinate[i](0,1))*(cordinate[i](0,0)-cordinate[i](0,1)) + (cordinate[i](1,0)-cordinate[i](1,1))*(cordinate[i](1,0)-cordinate[i](1,1)) + (cordinate[i](2,0)-cordinate[i](2,1))*(cordinate[i](2,0)-cordinate[i](2,1)) );
                        if(l-tol<=length && length<=l+tol){
                            if(it==t.end()){
                                t.push_back(i);
                                Outfile<<i<<"; "<<Tips<<"; "<<setprecision(15)<<length<<endl;
                            }
                        }
                    }

                }
            }
        }
        t.clear();
        lungNP.clear();
        lungP.clear();
    }
    Outfile.close();
    return cordinate;

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










