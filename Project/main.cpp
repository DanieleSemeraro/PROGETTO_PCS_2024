#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include"Eigen/Eigen"
#include"utilis.hpp"

using namespace std;
using namespace Eigen;


int main() {
    unsigned int n; //numero inserito dall utente per decidere quante fratture visualizzare
    string filename; //nome file
    vector<Vector2i> fratturescluse;//conterra tutte le coppie di fratture che sicuramente non creeranno tracce
    int NumberOfTraces=0; // numero totale di tracce su ogni rettangolo
    DFNLibrary::Fractures fractures;//chiamo la struct Fractures
    DFNLibrary::Traces traces;//chiamo la struct Traces
    cout<<"Inserire il numero di DFN da analizzare (3 10 50 82 200 362): "<<endl;
    cin>>n;
    switch (n) { //switch per decidere

    case 3:
        filename="DFN/FR3_data.txt";
        DFNLibrary::ImportDFN(filename,n, fractures);
        cout<<"Id delle fratture: "<<fractures.FractureId<<endl;
        cout<<"Numero di vertici di ogni frattura: "<<fractures.NumVertices<<endl;
        cout<<setprecision(15)<<fractures.ListVertices;
        fractures.sfere(n,fratturescluse);
        traces.CalcoloDirezioneTracce(NumberOfTraces,fractures, n,fratturescluse);
        traces.CalcoloEstremi(NumberOfTraces,fractures);
        traces.Ordinamento(fractures);
        break;
    case 10:
        filename="DFN/FR10_data.txt";
        DFNLibrary::ImportDFN(filename,n, fractures);
        cout<<"Id delle fratture: "<<fractures.FractureId<<endl;
        cout<<"Numero di vertici di ogni frattura: "<<fractures.NumVertices<<endl;
        cout<<setprecision(15)<<fractures.ListVertices;
        fractures.sfere(n,fratturescluse);
        traces.CalcoloDirezioneTracce(NumberOfTraces,fractures, n,fratturescluse);
        traces.CalcoloEstremi(NumberOfTraces,fractures);
        traces.Ordinamento(fractures);
        break;
    case 50:
        filename="DFN/FR50_data.txt";
        DFNLibrary::ImportDFN(filename,n, fractures);
        cout<<"Id delle fratture: "<<fractures.FractureId<<endl;
        cout<<"Numero di vertici di ogni frattura: "<<fractures.NumVertices<<endl;
        cout<<setprecision(15)<<fractures.ListVertices;
        fractures.sfere(n,fratturescluse);
        traces.CalcoloDirezioneTracce(NumberOfTraces,fractures, n,fratturescluse);
        traces.CalcoloEstremi(NumberOfTraces,fractures);
        traces.Ordinamento(fractures);
        break;
    case 82:
        filename="DFN/FR82_data.txt";
        DFNLibrary::ImportDFN(filename,n, fractures);
        cout<<"Id delle fratture: "<<fractures.FractureId<<endl;
        cout<<"Numero di vertici di ogni frattura: "<<fractures.NumVertices<<endl;
        cout<<setprecision(15)<<fractures.ListVertices;
        fractures.sfere(n,fratturescluse);
        traces.CalcoloDirezioneTracce(NumberOfTraces,fractures, n,fratturescluse);
        traces.CalcoloEstremi(NumberOfTraces,fractures);
        traces.Ordinamento(fractures);
        break;
    case 200:
        filename="DFN/FR200_data.txt";
        DFNLibrary::ImportDFN(filename,n, fractures);
        cout<<"Id delle fratture: "<<fractures.FractureId<<endl;
        cout<<"Numero di vertici di ogni frattura: "<<fractures.NumVertices<<endl;
        cout<<setprecision(15)<<fractures.ListVertices;
        fractures.sfere(n,fratturescluse);
        traces.CalcoloDirezioneTracce(NumberOfTraces,fractures, n,fratturescluse);
        traces.CalcoloEstremi(NumberOfTraces,fractures);
        traces.Ordinamento(fractures);
        break;
    case 362:
        filename="DFN/FR362_data.txt";
        DFNLibrary::ImportDFN(filename,n, fractures);
        cout<<"Id delle fratture: "<<fractures.FractureId<<endl;
        cout<<"Numero di vertici di ogni frattura: "<<fractures.NumVertices<<endl;
        cout<<setprecision(15)<<fractures.ListVertices;
        fractures.sfere(n,fratturescluse);
        traces.CalcoloDirezioneTracce(NumberOfTraces,fractures, n,fratturescluse);
        traces.CalcoloEstremi(NumberOfTraces,fractures);
        traces.Ordinamento(fractures);
        break;
    default:
        cout<<"Data DFN non a sistema" << endl;

    }

    return 0;
}

//Andrea Rostagno 295706
//Jacopo Ferraris 294292
//Daniele Semeraro 284238
