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

vector<double>CalcoloTraccePassanti(int n,vector<double> &FractureId,vector<double> &NumVertices,vector<MatrixXd> &ListVertices,vector<MatrixXd> &ListCord);

MatrixXd SystemSolve(vector<MatrixXd> &ListVertices,Vector3d P0,Vector3d P1,Vector3d P2,Vector3d P3,MatrixXd A,Vector3d b,int i,int j,int &c,vector<double> NumVertices);

vector<double>CalcoloTracce(vector<vector<double>> &IDs,int n,vector<double> &FractureId,vector<double> &NumVertices,vector<MatrixXd> &ListVertices,vector<MatrixXd> &ListCord);

Vector3d SoluzioneSistema(Vector3d &P,vector<MatrixXd> &ListVertices,Vector3d P0,Vector3d P1,Vector3d P2,Vector3d P3,MatrixXd A,Vector3d b,int i,int j,vector<double> NumVertices);


#endif
