#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <iomanip>
#include"Eigen/Eigen"
#include"utilis.hpp"

using namespace std;
using namespace Eigen;


int main() {
    int n;
    string filename;
    vector<double> FractureId;
    vector<double> NumVertices;
    vector<MatrixXd>ListVertices;
    vector<MatrixXd>ListCord;//contien un ounto p e la direzione della retta della traccia
    vector<VectorXd>IDs;//vettori binari corrispondenti agli id delle fratture per ogni singola traccia
    int NumberOfTraces=0;
    vector<MatrixXd> cordinate;
    vector<double> pass;//vettore che in ogni posizione ha 1 o 0 in base a se è passante o non passante

    cout<<"Inserire il numero di DFN da analizzare (3 10 50 82 200 362):"<<endl;
    cin>>n;
    switch (n) {
    case 3:
        filename="DFN/FR3_data.txt";
        ImportDFN(filename,n,FractureId,NumVertices,ListVertices);
        cout<<"Id delle fratture: "<<FractureId<<endl;
        cout<<"Numero di vertici di ogni frattura: "<<NumVertices<<endl;
        cout<<setprecision(15)<<ListVertices;
        break;
    case 10:
        filename="DFN/FR10_data.txt";
        ImportDFN(filename,n,FractureId,NumVertices,ListVertices);
        cout<<"Id delle fratture: "<<FractureId<<endl;
        cout<<"Numero di vertici di ogni frattura: "<<NumVertices<<endl;
        cout<<setprecision(15)<<ListVertices;
        break;
    case 50:
        filename="DFN/FR50_data.txt";
        ImportDFN(filename,n,FractureId,NumVertices,ListVertices);
        cout<<"Id delle fratture: "<<FractureId<<endl;
        cout<<"Numero di vertici di ogni frattura: "<<NumVertices<<endl;
        cout<<setprecision(15)<<ListVertices;
        break;
    case 82:
        filename="DFN/FR82_data.txt";
        ImportDFN(filename,n,FractureId,NumVertices,ListVertices);
        cout<<"Id delle fratture: "<<FractureId<<endl;
        cout<<"Numero di vertici di ogni frattura: "<<NumVertices<<endl;
        cout<<setprecision(15)<<ListVertices;
        break;
    case 200:
        filename="DFN/FR200_data.txt";
        ImportDFN(filename,n,FractureId,NumVertices,ListVertices);
        cout<<"Id delle fratture: "<<FractureId<<endl;
        cout<<"Numero di vertici di ogni frattura: "<<NumVertices<<endl;
        cout<<setprecision(15)<<ListVertices;
        break;
    case 362:
        filename="DFN/FR362_data.txt";
        ImportDFN(filename,n,FractureId,NumVertices,ListVertices);
        cout<<"Id delle fratture: "<<FractureId<<endl;
        cout<<"Numero di vertici di ogni frattura: "<<NumVertices<<endl;
        cout<<setprecision(15)<<ListVertices;
        break;
    default:
        cout<<"Data DFN non a sistema" << endl;

    }
    CalcoloDirezioneTracce(NumberOfTraces,IDs, n,FractureId,NumVertices,ListVertices,ListCord);

    cordinate=CalcoloEstremi(NumberOfTraces,IDs,NumVertices,ListVertices,ListCord,pass);

    Ordinamento(FractureId,IDs,cordinate,pass);

    return 0;
}
