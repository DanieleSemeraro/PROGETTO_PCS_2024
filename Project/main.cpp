#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include"Eigen/Eigen"
#include"utilis.hpp"

using namespace std;
using namespace Eigen;


int main() {
    int n;//numero inserito dall utente per decidere quante fratture visualizzare
    string filename;//nome file
    vector<double> FractureId;// id delle fratture
    vector<double> NumVertices;// numero vertici delle fratture
    vector<MatrixXd>ListVertices;//matrici con le cordinate dei vertici
    vector<MatrixXd>ListCord;//contien un punto p e la direzione della retta della traccia t
    vector<VectorXd>IDs;//vettori binari corrispondenti agli id delle fratture per ogni singola traccia
    int NumberOfTraces=0;// numero totale di tracce su ogni rettangolo
    vector<MatrixXd> cordinate;// vettore di matrici 3x2 corrispondenti agli estremi delle fratture
    vector<double> pass;//vettore che in ogni posizione ha 1(non passante) o 0(passante) in base a se è passante o non passante

    cout<<"Inserire il numero di DFN da analizzare (3 10 50 82 200 362):"<<endl;
    cin>>n;
    switch (n) {//switch per decidere
    case 3:
        filename="DFN/FR3_data.txt";
        ImportDFN(filename,n,FractureId,NumVertices,ListVertices);
        cout<<"Id delle fratture: "<<FractureId<<endl;
        cout<<"Numero di vertici di ogni frattura: "<<NumVertices<<endl;
        cout<<setprecision(15)<<ListVertices;
        CalcoloDirezioneTracce(NumberOfTraces,IDs, n,FractureId,NumVertices,ListVertices,ListCord);
        cordinate=CalcoloEstremi(NumberOfTraces,IDs,NumVertices,ListVertices,ListCord,pass);
        Ordinamento(FractureId,IDs,cordinate,pass);
        break;
    case 10:
        filename="DFN/FR10_data.txt";
        ImportDFN(filename,n,FractureId,NumVertices,ListVertices);
        cout<<"Id delle fratture: "<<FractureId<<endl;
        cout<<"Numero di vertici di ogni frattura: "<<NumVertices<<endl;
        cout<<setprecision(15)<<ListVertices;
        CalcoloDirezioneTracce(NumberOfTraces,IDs, n,FractureId,NumVertices,ListVertices,ListCord);
        cordinate=CalcoloEstremi(NumberOfTraces,IDs,NumVertices,ListVertices,ListCord,pass);
        Ordinamento(FractureId,IDs,cordinate,pass);
        break;
    case 50:
        filename="DFN/FR50_data.txt";
        ImportDFN(filename,n,FractureId,NumVertices,ListVertices);
        cout<<"Id delle fratture: "<<FractureId<<endl;
        cout<<"Numero di vertici di ogni frattura: "<<NumVertices<<endl;
        cout<<setprecision(15)<<ListVertices;
        CalcoloDirezioneTracce(NumberOfTraces,IDs, n,FractureId,NumVertices,ListVertices,ListCord);
        cordinate=CalcoloEstremi(NumberOfTraces,IDs,NumVertices,ListVertices,ListCord,pass);
        Ordinamento(FractureId,IDs,cordinate,pass);
        break;
    case 82:
        filename="DFN/FR82_data.txt";
        ImportDFN(filename,n,FractureId,NumVertices,ListVertices);
        cout<<"Id delle fratture: "<<FractureId<<endl;
        cout<<"Numero di vertici di ogni frattura: "<<NumVertices<<endl;
        cout<<setprecision(15)<<ListVertices;
        CalcoloDirezioneTracce(NumberOfTraces,IDs, n,FractureId,NumVertices,ListVertices,ListCord);
        cordinate=CalcoloEstremi(NumberOfTraces,IDs,NumVertices,ListVertices,ListCord,pass);
        Ordinamento(FractureId,IDs,cordinate,pass);
        break;
    case 200:
        filename="DFN/FR200_data.txt";
        ImportDFN(filename,n,FractureId,NumVertices,ListVertices);
        cout<<"Id delle fratture: "<<FractureId<<endl;
        cout<<"Numero di vertici di ogni frattura: "<<NumVertices<<endl;
        cout<<setprecision(15)<<ListVertices;
        CalcoloDirezioneTracce(NumberOfTraces,IDs, n,FractureId,NumVertices,ListVertices,ListCord);
        cordinate=CalcoloEstremi(NumberOfTraces,IDs,NumVertices,ListVertices,ListCord,pass);
        Ordinamento(FractureId,IDs,cordinate,pass);
        break;
    case 362:
        filename="DFN/FR362_data.txt";
        ImportDFN(filename,n,FractureId,NumVertices,ListVertices);
        cout<<"Id delle fratture: "<<FractureId<<endl;
        cout<<"Numero di vertici di ogni frattura: "<<NumVertices<<endl;
        cout<<setprecision(15)<<ListVertices;
        CalcoloDirezioneTracce(NumberOfTraces,IDs, n,FractureId,NumVertices,ListVertices,ListCord);
        cordinate=CalcoloEstremi(NumberOfTraces,IDs,NumVertices,ListVertices,ListCord,pass);
        Ordinamento(FractureId,IDs,cordinate,pass);
        break;
    default:
        cout<<"Data DFN non a sistema" << endl;

    }

    return 0;
}
