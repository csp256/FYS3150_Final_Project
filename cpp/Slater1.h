#include <iostream>
#include <new>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
//#define ARMA_NO_DEBUG
#include <armadillo>

using namespace std;
using namespace arma;

#ifndef DET_H
#define DET_H
	class Det {
		public:
		    double  m_omega, m_alpha;
		    mat m_r;

		    Det(mat& r, double omega, double alpha);
		    double det2();
		    double det6();
	};
#endif 

double Slater(mat&, int, int, double, double);