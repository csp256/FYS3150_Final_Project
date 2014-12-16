#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
//#define ARMA_NO_DEBUG
#include <armadillo>
#include <stdio.h>
#include <stdlib.h>
#include "E_pot.h"

using namespace std;
using namespace arma;

double E_pot::oscillator() {
    double m_e_potential_harm = 0;
    for (int i = 0; i < m_number_particles; i++) {
        m_r_single_particle = 0;
        for (int j = 0; j < m_dimension; j++) {
            double temp = m_r(i,j);
            m_r_single_particle += temp * temp;
        }
        m_e_potential_harm += m_r_single_particle;
    }
    return m_e_potential_harm*m_omega*m_omega * 0.5;
}

double E_pot::atom() {
    double m_e_potential_atom = 0;
    // contribution from electron-proton potential
    for (int i = 0; i < m_number_particles; i++) {
      m_r_single_particle = 0;
      for (int j = 0; j < m_dimension; j++) {
        double temp = m_r(i,j);
        m_r_single_particle += temp * temp;
      }
      m_e_potential_atom -= m_charge/sqrt(m_r_single_particle);
    }
    return m_e_potential_atom;
}

double E_pot::electronelectron() {
    // repulsion electron-electron
    double m_e_potential_electron = 0;
    for (int i = 0; i < m_number_particles; i++) {
        for (int j = i+1; j < m_number_particles; j++) {
            m_e_potential_electron += 1.0/m_r_12(i,j);
        }
    }
    return m_e_potential_electron;
}

E_pot::E_pot(double omega, mat& r, int dimension, int number_particles, int charge, mat& r_12, double r_single_particle) {
    m_omega = omega;
    m_r     = r;
    m_r_12  = r_12;
    m_r_single_particle = r_single_particle;
    m_dimension         = dimension;
    m_number_particles  = number_particles;
    m_charge = charge;
}