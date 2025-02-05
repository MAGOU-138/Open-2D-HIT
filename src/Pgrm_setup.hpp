#ifndef SIMULATION_CONFIG_HPP
#define SIMULATION_CONFIG_HPP

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <hdf5.h>
#include <complex>
#include <fftw3.h>
#include <mpi.h>
#include <iostream>
#include <filesystem>
using Complex = std::complex<double>;
namespace fs = std::filesystem;

namespace set
{
    // These are what you can modify
    constexpr int NX = 512;
    constexpr int NY = 512;
    constexpr int NP = 8;
    constexpr double dt = -5.e-5;
    constexpr int Nsave = 100;
    constexpr int NT = 100000;
    constexpr double nu = 1.e-4;

    // The rest can not be modified
    constexpr int NXH = NX / 2;
    constexpr int NYH = NY / 2;
    constexpr int NX2 = NX * 3 / 2;
    constexpr int NY2 = NY * 3 / 2;
    constexpr int NX2NP = NX2 / NP;
    constexpr int NY2NP = NY2 / NP;
    constexpr int NXHNP = NXH / NP;

    constexpr Complex imag_unit = Complex(0.0, 1.0);

    extern double kx_local[NXHNP];
    extern double ky_local[NY];
    extern Complex u[NXHNP][NY];
    extern Complex v[NXHNP][NY];
    extern Complex vor[NXHNP][NY], vor0[NXHNP][NY], vor1[NXHNP][NY], vor2[NXHNP][NY];
    extern Complex jac[NXHNP][NY], jac0[NXHNP][NY], jac1[NXHNP][NY], jac2[NXHNP][NY];

    extern int id_step;

    extern fftw_plan plan_xf, plan_xb, plan_yf, plan_yb;
    extern fftw_complex *fft_temp_x1;
    extern fftw_complex *fft_temp_x2;
    extern fftw_complex *fft_temp_y1;
    extern fftw_complex *fft_temp_y2;

    extern int rank;
    extern int numProcs;

    void initField();

    constexpr double pi = 3.141592653589793238462643383279502884197;

    // Init wave numbers
}

// program

void checkCase();
void initCase();
void timeIntegration();
void solvePoisson_negsource(const Complex (&source_)[set::NXHNP][set::NY], Complex (&var_)[set::NXHNP][set::NY]);
void solveVelocity(const Complex (&phi_)[set::NXHNP][set::NY]);
void EXCHANGE_Y2X(const Complex (&var_y)[set::NXHNP][set::NY2], Complex (&var_x)[set::NXH][set::NY2NP]);
void EXCHANGE_X2Y(const Complex (&var_x)[set::NXH][set::NY2NP], Complex (&var_y)[set::NXHNP][set::NY2]);

void FOURIER_S2P_DEALIAS(const Complex (&var_)[set::NXHNP][set::NY], double (&var_p)[set::NX2][set::NY2NP]);
void FOURIER_P2S_DEALIAS(const double (&var_p)[set::NX2][set::NY2NP], Complex (&var_)[set::NXHNP][set::NY]);

void nonlinear(const Complex (&u_)[set::NXHNP][set::NY], const Complex (&v_)[set::NXHNP][set::NY], Complex (&res_)[set::NXHNP][set::NY]);

hid_t FILE_CREATE(const std::string &IO_FILENAME);
void WRITE_2D(hid_t file_id, const std::string &dataset_name, const double (&data)[set::NXHNP][set::NY]);
void save_Instantaneous(const int id_step_);

#endif