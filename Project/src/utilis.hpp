#ifndef UTILIS_HPP
#define UTILIS_HPP

#include <iostream>
#include <vector>
#include"Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace DFNLibrary {

struct Fractures{

    vector<int> FractureId;// id delle fratture
    vector<int> NumVertices;// numero vertici delle fratture
    vector<MatrixXd> ListVertices;//matrici con le cordinate dei vertici

    Fractures() = default;
    Fractures(const vector<int> &FractureId,
              const vector<int> &NumVertices,
              const vector<MatrixXd> &ListVertices):
        FractureId(FractureId),// Inizializza il membro FractureId con il parametro FractureId
        NumVertices(NumVertices),//
        ListVertices(ListVertices)//
    {}

};
struct Traces{

    vector<MatrixXd> ListCord;//contien un punto p e la direzione della retta della traccia t
    vector<VectorXi> IDs;//vettori binari corrispondenti agli id delle fratture per ogni singola traccia
    VectorXi pass; // vettore che in ogni posizione ha 1(non passante) o 0(passante) in base a se è passante o non passante
    vector<MatrixXd> cordinate;//vettore di matrici contenenti tutti gli estremi delle fratture
    Traces() = default;
    Traces(const vector<MatrixXd> &ListCord,
              const vector<VectorXi> &IDs,
              const vector<MatrixXd> &cordinate,
              const VectorXi& pass):
        ListCord(ListCord),// Inizializza il membro ListCord con il parametro FractureId
        IDs(IDs),//
        pass(pass),//
        cordinate(cordinate)//
    {}

    void CalcoloEstremi(int &NumberOfTraces,Fractures& fractures,Traces& traces);//funzione per trovare gli estremi delle fratture
    void Ordinamento(Fractures& fractures);//ultima funzione che permette di calcolare il numero di tracce presenti su ogni frattura, la loro lunghezza e ordinarle in maniera decrescente
    void CalcoloDirezioneTracce(int &NumberOfTraces,Fractures& fractures,int n,Traces& traces,vector<Vector2i> &fratturescluse);//serve per ottenere la retta su cui somo presenti le tracce, restituisce punto e direzione di ogni retta utile con tracce

};

void ImportDFN(const string& filename,int n,Fractures& fractures);//funzione che importa i dati dai file DFN
void sfere(Fractures& fractures, int n,vector<Vector2i> &fratturescluse);

}
ostream& operator<<(ostream& os,const vector<int> a);

ostream& operator<<(ostream& os,const vector<MatrixXd> a);

void BubbleSort(vector<double>& data);

#endif

