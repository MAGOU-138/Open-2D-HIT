#include "Pgrm_setup.hpp"

void FOURIER_S2P_DEALIAS(const Complex (&var_)[set::NXHNP][set::NY], double (&var_p)[set::NX2][set::NY2NP])
{
    Complex var_x1y2_py[set::NXHNP][set::NY2];
    // zero pad in ky
    for (size_t i = 0; i < set::NXHNP; i++)
    {
        for (size_t j = 0; j < set::NYH; j++)
        {
            var_x1y2_py[i][j] = var_[i][j];
            var_x1y2_py[i][j + set::NYH] = 0.0;
            var_x1y2_py[i][j + set::NY] = var_[i][j + set::NYH];
        }
    }

    // IFFT in y
    for (size_t i = 0; i < set::NXHNP; i++)
    {
        for (size_t j = 0; j < set::NY2; j++)
        {
            set::fft_temp_y1[j][0] = var_x1y2_py[i][j].real();
            set::fft_temp_y1[j][1] = var_x1y2_py[i][j].imag();
        }
        set::fft_temp_y1[set::NY][0] = 0; // set Nyquist mode zero
        set::fft_temp_y1[set::NY][1] = 0; // set Nyquist mode zero

        fftw_execute(set::plan_yb); // IFFT

        for (size_t j = 0; j < set::NY2; j++)
        {
            var_x1y2_py[i][j] = Complex(set::fft_temp_y2[j][0], set::fft_temp_y2[j][1]);
        }
    }

    // exchange parallel direction
    Complex var_x1y2_px[set::NXH][set::NY2NP];
    EXCHANGE_Y2X(var_x1y2_py, var_x1y2_px);

    // zero pad in kx
    Complex var_x2y2_px[set::NX2][set::NY2NP];
    for (size_t j = 0; j < set::NY2NP; j++)
    {
        for (size_t i = 0; i < set::NXH; i++)
        {
            var_x2y2_px[i][j] = var_x1y2_px[i][j];
            var_x2y2_px[i + set::NXH][j] = 0.0;
            var_x2y2_px[i + set::NX][j] = std::conj(var_x1y2_px[set::NXH - i][j]);
        }
    }

    // IFFT in x
    for (size_t j = 0; j < set::NY2NP; j++)
    {
        for (size_t i = 0; i < set::NX2; i++)
        {
            set::fft_temp_x1[i][0] = var_x2y2_px[i][j].real();
            set::fft_temp_x1[i][1] = var_x2y2_px[i][j].imag();
        }
        set::fft_temp_x1[set::NX][0] = 0;
        set::fft_temp_x1[set::NX][1] = 0;

        fftw_execute(set::plan_xb);

        for (size_t i = 0; i < set::NX2; i++)
        {
            var_x2y2_px[i][j] = Complex(set::fft_temp_x2[i][0], set::fft_temp_x2[i][1]);
        }
    }

    // variable in physical space (refined in space disceritization essentially)
    for (size_t i = 0; i < set::NX2; i++)
    {
        const Complex *row_var = var_x2y2_px[i];
        double *row_var_p = var_p[i];
        for (size_t j = 0; j < set::NY2NP; j++)
        {
            row_var_p[j] = row_var[j].real();
        }
    }
}

void FOURIER_P2S_DEALIAS(const double (&var_p)[set::NX2][set::NY2NP], Complex (&var_)[set::NXHNP][set::NY])
{
    // FFT in x
    Complex var_x2y2_px[set::NX2][set::NY2NP];
    for (size_t j = 0; j < set::NY2NP; j++)
    {
        for (size_t i = 0; i < set::NX2; i++)
        {
            set::fft_temp_x1[i][0] = var_p[i][j];
            set::fft_temp_x1[i][1] = 0.0;
        }

        fftw_execute(set::plan_xf);

        for (size_t i = 0; i < set::NX2; i++)
        {
            var_x2y2_px[i][j] = Complex(set::fft_temp_x2[i][0], set::fft_temp_x2[i][1]) / (double)set::NX2;
        }
    }

    // zero pad in kx
    Complex var_x1y2_px[set::NXH][set::NY2NP];
    for (size_t j = 0; j < set::NY2NP; j++)
    {
        for (size_t i = 0; i < set::NXH; i++)
        {
            var_x1y2_px[i][j] = var_x2y2_px[i][j];
        }
    }

    // exchange paralle direction
    Complex var_x1y2_py[set::NXHNP][set::NY2];
    EXCHANGE_X2Y(var_x1y2_px, var_x1y2_py);

    // FFT in y
    for (size_t i = 0; i < set::NXHNP; i++)
    {
        for (size_t j = 0; j < set::NY2; j++)
        {
            set::fft_temp_y1[j][0] = var_x1y2_py[i][j].real();
            set::fft_temp_y1[j][1] = var_x1y2_py[i][j].imag();
        }

        fftw_execute(set::plan_yf);

        for (size_t j = 0; j < set::NY2; j++)
        {
            var_x1y2_py[i][j] = Complex(set::fft_temp_y2[j][0], set::fft_temp_y2[j][1]) / (double)set::NY2;
        }
    }

    // zero pad in ky
    for (size_t i = 0; i < set::NXHNP; i++)
    {
        for (size_t j = 0; j < set::NYH; j++)
        {
            var_[i][j] = var_x1y2_py[i][j];
            var_[i][j + set::NYH] = var_x1y2_py[i][j + set::NY];
        }
    }
}