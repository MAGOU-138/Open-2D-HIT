

#include <iostream>
#include <stdexcept>
#include <string>
#include "Pgrm_setup.hpp"

/*
2d turbulence
program by Chutian Wu
*/

namespace set
{
    int numProcs = 1;
    int rank = 0;
    int id_step = 0;
}

int main(int argc, char **argv)
{

    // Init MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &set::numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &set::rank);

    checkCase();
    initCase();

    // init field
    set::initField();

    for (size_t i = 0; i < set::NT; i++)
    {
        if (set::rank == 0)
        {
            std::cout << "time step " << set::id_step << std::endl;
        }

        if (set::id_step % set::Nsave == 0)
        {
            save_Instantaneous(set::id_step);
        }
        timeIntegration();
        set::id_step += 1;
    }

    // Free FFT
    fftw_destroy_plan(set::plan_xf);
    fftw_destroy_plan(set::plan_xb);
    fftw_destroy_plan(set::plan_yf);
    fftw_destroy_plan(set::plan_yb);
    fftw_free(set::fft_temp_x1);
    fftw_free(set::fft_temp_x2);
    fftw_free(set::fft_temp_y1);
    fftw_free(set::fft_temp_y2);
    // Finalize MPI
    MPI_Finalize();

    return 0;
}
