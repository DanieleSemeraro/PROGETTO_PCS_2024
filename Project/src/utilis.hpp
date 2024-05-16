#ifndef UTILIS_HPP
#define UTILIS_HPP

#include <iostream>
#include <vector>
#include"Eigen/Eigen"

using namespace std;
using namespace Eigen;


vector<double>ImportDFN(string filename,int n,vector<double> &FractureId,vector<double> &NumVertices,vector<MatrixXd> &ListVertices);

ostream& operator<<(ostream& os,const vector<double> a);

ostream& operator<<(ostream& os,const vector<MatrixXd> a);

vector<VectorXd> CalcoloDirezioneTracce(int &NumberOfTraces,vector<VectorXd> &IDs,int n,vector<double> &FractureId,vector<double> &NumVertices,vector<MatrixXd> &ListVertices,vector<MatrixXd> &ListCord);

vector<MatrixXd> CalcoloEstremi(int &NumberOfTraces,vector<VectorXd> &IDs,vector<double> &NumVertices,vector<MatrixXd> &ListVertices,vector<MatrixXd> &ListCord,vector<double> &pass);

vector<MatrixXd> Ordinamento(vector<double> FractureId,int &NumberOfTraces,vector<VectorXd> &IDs,vector<double> &NumVertices,vector<MatrixXd> &ListVertices,vector<MatrixXd> &ListCord, vector<MatrixXd> &cordinate,vector<double> &pass);

void BubbleSort(vector<double>& data);

#endif
