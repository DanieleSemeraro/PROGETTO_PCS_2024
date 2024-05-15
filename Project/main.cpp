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
    vector<MatrixXd>ListCord;//cordinate tracce
    vector<VectorXd>IDs;//vettori binari corrispondenti agli id delle fratture per ogni singola traccia
    int NumberOfTraces=0;

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
    //cout<<"Il numero di tracce: "<<NumberOfTraces<<endl;
    //for (unsigned int i = 0; i < IDs.size(); ++i) {
        //cout<<"TraceID: "<<i<<" Fracture IDs: "<<IDs[i].transpose()<<endl;
        //cout<<"Punto P: "<<setprecision(15)<<ListCord[i].col(0).transpose()<<endl<<"Direzione t: "<<setprecision(15)<<ListCord[i].col(1).transpose()<<endl;

    //}
    CalcoloEstremi(NumberOfTraces,IDs,NumVertices,ListVertices,ListCord);

















    return 0;
}
