#include <iostream>
#include <new>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
//#define ARMA_NO_DEBUG
#include <armadillo>
#include "lib.h"
#include <vector>
#include "Random.h"

using namespace std;
using namespace arma;

void mc_sampling(double, double, double, int, int, int,
                 int, double, double,
                 int, double,
                 vec&, mat&,
                 vec&, vec&, mat&, mat&,
                 vec&, vec&, mat&, mat&, int&, vector<Random*> &randoms);

void initialise(double &, double&, double&, int&, double &, double &, int&, int&, int&, int&, double&) ;
double  local_energy(double&, double&, double, mat&, int, double, double , double, int, int, int, int, mat&);
void distance(mat&, mat&, int, int );