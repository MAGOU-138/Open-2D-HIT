#include "Pgrm_setup.hpp"

fftw_complex *set::fft_temp_x1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * set::NX2);
fftw_complex *set::fft_temp_x2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * set::NX2);
fftw_complex *set::fft_temp_y1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * set::NY2);
fftw_complex *set::fft_temp_y2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * set::NY2);

fftw_plan set::plan_xf = fftw_plan_dft_1d(set::NX2, set::fft_temp_x1, set::fft_temp_x2, FFTW_FORWARD, FFTW_ESTIMATE);
fftw_plan set::plan_yf = fftw_plan_dft_1d(set::NY2, set::fft_temp_y1, set::fft_temp_y2, FFTW_FORWARD, FFTW_ESTIMATE);
fftw_plan set::plan_xb = fftw_plan_dft_1d(set::NX2, set::fft_temp_x1, set::fft_temp_x2, FFTW_BACKWARD, FFTW_ESTIMATE);
fftw_plan set::plan_yb = fftw_plan_dft_1d(set::NY2, set::fft_temp_y1, set::fft_temp_y2, FFTW_BACKWARD, FFTW_ESTIMATE);

double set::kx_local[set::NXHNP] = {};
double set::ky_local[set::NY] = {};
Complex set::u[set::NXHNP][set::NY] = {};
Complex set::v[set::NXHNP][set::NY] = {};
Complex set::vor[set::NXHNP][set::NY] = {};
Complex set::vor0[set::NXHNP][set::NY] = {};
Complex set::vor1[set::NXHNP][set::NY] = {};
Complex set::vor2[set::NXHNP][set::NY] = {};
Complex set::jac[set::NXHNP][set::NY] = {};
Complex set::jac0[set::NXHNP][set::NY] = {};
Complex set::jac1[set::NXHNP][set::NY] = {};
Complex set::jac2[set::NXHNP][set::NY] = {};