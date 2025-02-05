#include "Pgrm_setup.hpp"

void EXCHANGE_Y2X(const Complex (&var_y)[set::NXHNP][set::NY2], Complex (&var_x)[set::NXH][set::NY2NP])
{
    // This function changes parallel direction from y to x

    Complex send[set::NP][set::NY2NP][set::NXHNP];
    Complex recv[set::NP][set::NY2NP][set::NXHNP];
    constexpr int all2all_size = set::NXHNP * set::NY2NP;

    // The send data package
    for (size_t i = 0; i < set::NXHNP; ++i)
    {
        for (size_t j = 0; j < set::NY2NP; ++j)
        {
            for (size_t p = 0; p < set::NP; ++p)
            {
                size_t global_j = j + p * set::NY2NP; // Global y-index
                send[p][j][i] = var_y[i][global_j];
            }
        }
    }

    MPI_Alltoall(send, all2all_size, MPI_DOUBLE_COMPLEX, recv, all2all_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);

    // Unpack the revc buffer into var_x
    for (size_t i = 0; i < set::NXHNP; i++)
    {
        for (size_t j = 0; j < set::NY2NP; j++)
        {
            for (size_t p = 0; p < set::NP; p++)
            {
                size_t global_i = i + p * set::NXHNP;
                var_x[global_i][j] = recv[p][j][i];
            }
        }
    }
}

void EXCHANGE_X2Y(const Complex (&var_x)[set::NXH][set::NY2NP], Complex (&var_y)[set::NXHNP][set::NY2])
{
    // This function changes parallel direction from x to y

    Complex send[set::NP][set::NY2NP][set::NXHNP];
    Complex recv[set::NP][set::NY2NP][set::NXHNP];
    constexpr int all2all_size = set::NXHNP * set::NY2NP;

    // The send data package
    for (size_t i = 0; i < set::NXHNP; ++i)
    {
        for (size_t j = 0; j < set::NY2NP; ++j)
        {
            for (size_t p = 0; p < set::NP; ++p)
            {
                size_t global_i = i + p * set::NXHNP; // Global x-index
                send[p][j][i] = var_x[global_i][j];
            }
        }
    }

    MPI_Alltoall(send, all2all_size, MPI_DOUBLE_COMPLEX, recv, all2all_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);

    // Unpack the revc buffer into var_x
    for (size_t i = 0; i < set::NXHNP; i++)
    {
        for (size_t j = 0; j < set::NY2NP; j++)
        {
            for (size_t p = 0; p < set::NP; p++)
            {
                size_t global_j = j + p * set::NY2NP;
                var_y[i][global_j] = recv[p][j][i];
            }
        }
    }
}