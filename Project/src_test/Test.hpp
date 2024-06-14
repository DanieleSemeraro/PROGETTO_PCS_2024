#ifndef __TEST__H
#define __TEST__H

#include <gtest/gtest.h>
#include "utilis.hpp"
#include "Eigen/Eigen"


using namespace std;
using namespace Eigen;
namespace DFNLibrary {

bool VectorsAreEquali(const VectorXi& vec1, const VectorXi& vec2) {
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

bool VectorsAreEquald (const VectorXd& vec1, const VectorXd& vec2, double tol =1e-8){// test sui double con tolleranza
    if (vec1.size() != vec2.size()) {
        return false;
    }
    for (Eigen::Index i = 0; i < vec1.size(); ++i) {
        if (std::abs(vec1[i] - vec2[i]) > tol) {
            return false;
        }
    }
    return true;
}

TEST(TRACESTEST, TestCalcoloDirezioneTracce){ // test sul caso più semplice (3 fratture) per assicurarsi che il codice funzioni

    int n=3; //numero inserito dall utente per decidere quante fratture visualizzare
    string filename="DFN/FR3_data.txt"; //nome file
    vector<Vector2i> fratturescluse;//lista delle coppie di fratture che sicuramente non avranno tracce
    int NumberOfTraces=0; // numero totale di tracce su ogni rettangolo
    DFNLibrary::Fractures fractures;//chiamo la struct Fractures
    DFNLibrary::Traces traces;
    DFNLibrary::ImportDFN(filename,n, fractures); // la funzione restituisce FractureId e aggiorna NumVertices e ListVertices
    // un test su take funzione non è necessario in quanto stampa a terminale
    fractures.sfere(n,fratturescluse);
    traces.CalcoloDirezioneTracce(NumberOfTraces,fractures, n,fratturescluse);
    // la funzione restituisce IDs, aggiorna NumberOfTraces, ListCord

    EXPECT_EQ(NumberOfTraces, 2);  // check su numero di tracce
    EXPECT_EQ(traces.IDs.size(), 2);
    EXPECT_EQ(traces.ListCord.size(), 2);

    VectorXi first_trace(2); // la prima traccia salvata in IDs è quella
    first_trace << 0,1; // tra la frattura 0 e la frattura 1

    VectorXi second_trace(2);
    second_trace << 0,2;

    EXPECT_TRUE(VectorsAreEquali(traces.IDs[0], first_trace)); //controllo che IDs salvi gli elementi in modo corretto
    EXPECT_TRUE(VectorsAreEquali(traces.IDs[1], second_trace));

    Vector3d DirFirstTrace; // il vettore direzionale della prima traccia è (0,-1,0)
    DirFirstTrace << 0,-1,0;
    Vector3d t1 = traces.ListCord[0].col(1);

    Vector3d DirSecondTrace;
    DirSecondTrace << 1,0,0;
    Vector3d t2 =traces.ListCord[1].col(1);

    EXPECT_TRUE(VectorsAreEquald(DirFirstTrace, t1));
    EXPECT_TRUE(VectorsAreEquald(DirSecondTrace, t2));
}


TEST(TRACESTEST, TestCalcoloEstremi) {

    int n=3; //numero inserito dall utente per decidere quante fratture visualizzare
    string filename="DFN/FR3_data.txt"; //nome file
    int NumberOfTraces=0; // numero totale di tracce su ogni rettangolo
    vector<Vector2i> fratturescluse;
    DFNLibrary::Fractures fractures;//chiamo la struct Fractures
    DFNLibrary::Traces traces;
    DFNLibrary::ImportDFN(filename,n, fractures); // la funzione restituisce FractureId e aggiorna NumVertices e ListVertices
    fractures.sfere(n,fratturescluse);
    traces.CalcoloDirezioneTracce(NumberOfTraces,fractures, n,fratturescluse);
    // la funzione restituisce IDs, aggiorna NumberOfTraces, ListCord
    traces.CalcoloEstremi(NumberOfTraces,fractures); // all'interno della funzione viene anche aggiornato il vettore pass

    EXPECT_EQ(traces.pass.size(),2);
    EXPECT_EQ(traces.cordinate.size(), 2);

    MatrixXd ExtremesFirstTrace(3,2);
    ExtremesFirstTrace <<   0.8, 0.8,
                            0, 1,
                            0, 0;

    EXPECT_TRUE(traces.cordinate[0].isApprox(ExtremesFirstTrace, 1e-8));

    MatrixXd ExtremesSecondTrace(3,2);
    ExtremesSecondTrace <<  0, 0.3161837,
                            0.5, 0.5,
                            0, 0;

    EXPECT_TRUE(traces.cordinate[1].isApprox(ExtremesSecondTrace, 1e-8));

    Vector2i TraccePassanti;
    TraccePassanti << 0, 1; // la prima traccia è passante, la seconda no

    EXPECT_TRUE(VectorsAreEquali(TraccePassanti, traces.pass));
}
}

#endif
