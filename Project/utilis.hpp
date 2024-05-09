#ifndef UTILIS_HPP
#define UTILIS_HPP

#include <iostream>
#include <vector>
#include"Eigen/Eigen"

using namespace std;
using namespace Eigen;


vector<double>ImportDFN(string filename,int n,vector<double> &FractureId,vector<double> &NumVertices,MatrixXd &Vertices);

ostream& operator<<(ostream& os,const vector<double> a);




#endif
