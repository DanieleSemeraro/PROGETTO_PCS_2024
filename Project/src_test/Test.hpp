#ifndef __TEST__H
#define __TEST_H

#include <gtest/gtest.h>
#include "utilis.hpp"
#include "Eigen/Eigen"
#include <iostream>

using namespace std;
using namespace Eigen;

bool VectorsAreEqual(const VectorXd& vec1, const VectorXd& vec2) {
    if (vec1.size() != vec2.size()) {
        return false;
    }
    for (Eigen::Index i = 0; i < vec1.size(); ++i) {
        if (vec1[i] != vec2[i]) {
            return false;
        }
    }
    return true;
}


TEST(FRACTURETEST, TestCalcoloDirezioneTracce){ // test sul caso più semplice (3 fratture) per assicurarsi che il codice funzioni

    int n = 3;
    string filename = "DFN/FR3_data.txt";

    vector<double> FractureId;
    vector<double> NumVertices;
    vector<MatrixXd> ListVertices;//matrici con le cordinate dei vertici
    vector<MatrixXd>ListCord; //vettore di matrici: ogni matrice contiene un punto p e la direzione della retta della traccia t
    vector<VectorXd>IDs;//vettori binari corrispondenti agli id delle fratture per ogni singola traccia
    int NumberOfTraces=0; // numero totale di tracce su ogni rettangolo

    ImportDFN(filename,n,FractureId,NumVertices,ListVertices); // la funzione restituisce FractureId e aggiorna NumVertices e ListVertices
    // un test su take funzione non è necessario in quanto stampa a terminale
    CalcoloDirezioneTracce(NumberOfTraces,IDs,n,FractureId,NumVertices,ListVertices,ListCord);
    // la funzione restituisce IDs, aggiorna NumberOfTraces, ListCord

    EXPECT_EQ(NumberOfTraces, 2);  // check su numero di tracce
    EXPECT_EQ(IDs.size(), 2);
    EXPECT_EQ(ListCord.size(), 2);

    VectorXd first_trace(2); // la prima traccia salvata in IDs è quella
    first_trace << 0,1; // tra la frattura 0 e la frattura 1

    VectorXd second_trace(2);
    second_trace << 0,2;

    EXPECT_TRUE(VectorsAreEqual(IDs[0], first_trace)); //controllo che IDs salvi gli elementi in modo corretto
    EXPECT_TRUE(VectorsAreEqual(IDs[1], second_trace));

    Vector3d DirFirstTrace; // il vettore direzionale della prima traccia è (0,-1,0)
    DirFirstTrace << 0,-1,0;
    Vector3d t1 = ListCord[0].col(1);

    Vector3d DirSecondTrace;
    DirSecondTrace << 1,0,0;
    Vector3d t2 = ListCord[1].col(1);

    EXPECT_TRUE(VectorsAreEqual(DirFirstTrace, t1));
    EXPECT_TRUE(VectorsAreEqual(DirSecondTrace, t2));
}


TEST(FRACTURETEST, TestCalcoloEstremi) {
    int n = 3;
    string filename = "DFN/FR3_data.txt";

    vector<double> FractureId;
    vector<double> NumVertices;
    vector<MatrixXd> ListVertices;//matrici con le cordinate dei vertici
    vector<MatrixXd>ListCord; //vettore di matrici: ogni matrice contiene un punto p e la direzione della retta della traccia t
    vector<VectorXd>IDs;//vettori binari corrispondenti agli id delle fratture per ogni singola traccia
    int NumberOfTraces=0; // numero totale di tracce su ogni rettangolo
    vector<MatrixXd> cordinate;// vettore di matrici 3x2 corrispondenti agli estremi delle fratture
    VectorXd pass;//vettore che in ogni posizione ha 1(non passante) o 0(passante) in base a se è passante o non passante

    ImportDFN(filename,n,FractureId,NumVertices,ListVertices); // la funzione restituisce FractureId e aggiorna NumVertices e ListVertices
    CalcoloDirezioneTracce(NumberOfTraces,IDs,n,FractureId,NumVertices,ListVertices,ListCord);
    // la funzione restituisce IDs, aggiorna NumberOfTraces, ListCord
    cordinate=CalcoloEstremi(NumberOfTraces,IDs,NumVertices,ListVertices,ListCord,pass); // all'interno della funzione viene anche aggiornato il vettore pass

    EXPECT_EQ(pass.size(),2);
    EXPECT_EQ(cordinate.size(), 2);

    MatrixXd ExtremesFirstTrace(3,2);
    ExtremesFirstTrace <<   0.8, 0.8,
                            0, 1,
                            0, 0;

    EXPECT_TRUE(cordinate[0].isApprox(ExtremesFirstTrace, 1e-8));

    MatrixXd ExtremesSecondTrace(3,2);
    ExtremesSecondTrace <<  0, 0.3161837,
                            0.5, 0.5,
                            0, 0;

    EXPECT_TRUE(cordinate[1].isApprox(ExtremesSecondTrace, 1e-8));

    Vector2d TraccePassanti;
    TraccePassanti << 0, 1; // la prima traccia è passante, la seconda no

    EXPECT_TRUE(VectorsAreEqual(TraccePassanti, pass));
}


#endif
