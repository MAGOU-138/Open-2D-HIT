#include "Pgrm_setup.hpp"

void checkCase()
{
    if (set::rank == 0)
    {
        if (set::numProcs != set::NP)
        {
            std::cout << "error: number of process not right" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        if (set::NX % set::numProcs != 0)
        {
            std::cout << "error: Nx cannot be divided by NP" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        if (set::NXH % set::numProcs != 0)
        {
            std::cout << "error: number of process cannot divide 0.5NX" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        if (set::NY2 % set::numProcs != 0)
        {
            std::cout << "error: number of process cannot divide 3/2NY" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }
}

void initCase()
{
    if (set::rank == 0)
    {
        if (!fs::exists("DATA"))
        {
            fs::create_directory("DATA");
        }
    }
}

void set::initField()
{
    // This function initilize flow field

    // The wavenumbers
    for (size_t j = 0; j < set::NYH; j++)
    {
        set::ky_local[j] = (double)j;
        set::ky_local[j + set::NYH] = (double)j - (double)set::NYH;
    }

    for (size_t i = 0; i < set::NXHNP; i++)
    {
        set::kx_local[i] = (double)(set::rank * set::NXHNP + i);
    }

    // Initialize velocity with single wave
    // if (set::rank == 0)
    // {
    //     set::u[1][1] = Complex(0.0, -0.25);
    //     set::u[1][63] = Complex(0.0, 0.25);

    //     set::v[1][1] = Complex(0.0, 0.25);
    //     set::v[1][63] = Complex(0.0, 0.25);
    // }

    // initialize velocity with specific energy spetcrum and random phase, the total kinetic energy is 0.5
    for (size_t i = 0; i < set::NXHNP; i++)
    {
        for (size_t j = 0; j < set::NY; j++)
        {
            double k2 = set::kx_local[i] * set::kx_local[i] + set::ky_local[j] * set::ky_local[j];
            double k = sqrt(k2);
            double E = 0.0;
            if (k > 0.0)
            {
                E = 1.0 / (k * (1.0 + powf(k, 4) / powf(6.0, 4)));
                double theta = static_cast<double>(std::rand()) / RAND_MAX * (2 * M_PI);
                set::u[i][j] = set::ky_local[j] / k * sqrt(2 * E) * Complex(cos(theta), sin(theta));
                set::v[i][j] = -set::ky_local[j] / k * sqrt(2 * E) * Complex(cos(theta), sin(theta));
                // std::cout << set::kx_local[i] << ' ' << set::ky_local[j] << ' ' << set::u[i][j] << ' ' << set::v[i][j] << std::endl;
            }
        }
    }

    // calculate vorticity by velocity
    for (size_t i = 0; i < set::NXHNP; i++)
    {
        for (size_t j = 0; j < set::NY; j++)
        {
            set::vor[i][j] = set::u[i][j] * set::imag_unit * set::ky_local[j] - set::v[i][j] * set::imag_unit * set::kx_local[i];
        }
    }

    // get the nonlinear part
    Complex dvordx[set::NXHNP][set::NY];
    Complex dvordy[set::NXHNP][set::NY];
    Complex jac_x[set::NXHNP][set::NY];
    Complex jac_y[set::NXHNP][set::NY];

    for (size_t i = 0; i < set::NXHNP; i++)
    {
        for (size_t j = 0; j < set::NY; j++)
        {
            dvordx[i][j] = set::vor[i][j] * set::imag_unit * set::kx_local[i];
            dvordy[i][j] = set::vor[i][j] * set::imag_unit * set::ky_local[j];
        }
    }
    nonlinear(set::u, dvordx, jac_x);
    nonlinear(set::v, dvordy, jac_y);
    for (size_t i = 0; i < set::NXHNP; i++)
    {
        for (size_t j = 0; j < set::NY; j++)
        {
            set::jac[i][j] = jac_x[i][j] + jac_y[i][j];
        }
    }

    // luanch initialize
    for (size_t i = 0; i < set::NXHNP; i++)
    {
        for (size_t j = 0; j < set::NY; j++)
        {
            set::vor2[i][j] = set::vor[i][j];
            set::vor1[i][j] = set::vor[i][j];
            set::vor0[i][j] = set::vor[i][j];

            set::jac2[i][j] = set::jac[i][j];
            set::jac1[i][j] = set::jac[i][j];
            set::jac0[i][j] = set::jac[i][j];
        }
    }
}

void solvePoisson_negsource(const Complex (&source_)[set::NXHNP][set::NY], Complex (&var_)[set::NXHNP][set::NY])
{
    // This function solves the Poisson equation: [p_{xx}+p_{yy}]var_=-source_
    for (size_t i = 0; i < set::NXHNP; i++)
    {
        for (size_t j = 0; j < set::NY; j++)
        {
            double k2 = set::kx_local[i] * set::kx_local[i] + set::ky_local[j] * set::ky_local[j];
            if (k2 == 0.0)
            {
                var_[i][j] = 0.0;
            }
            else
            {
                var_[i][j] = source_[i][j] / k2;
            }
        }
    }
}

void solveVelocity(const Complex (&phi_)[set::NXHNP][set::NY])
{
    // u_= d(phi_)/dy
    for (size_t i = 0; i < set::NXHNP; i++)
    {
        for (size_t j = 0; j < set::NY; j++)
        {
            set::u[i][j] = phi_[i][j] * set::imag_unit * set::ky_local[j];
            set::v[i][j] = -phi_[i][j] * set::imag_unit * set::kx_local[i];
        }
    }
}

void nonlinear(const Complex (&u_)[set::NXHNP][set::NY], const Complex (&v_)[set::NXHNP][set::NY], Complex (&res_)[set::NXHNP][set::NY])
{
    // This function calculates Fourier modes of u*v

    // get u,v in physical space (refined)
    double up[set::NX2][set::NY2NP];
    double vp[set::NX2][set::NY2NP];
    double res_p[set::NX2][set::NY2NP];

    FOURIER_S2P_DEALIAS(u_, up);
    FOURIER_S2P_DEALIAS(v_, vp);

    // multiplication in physical space
    for (size_t i = 0; i < set::NX2; ++i)
    {
        for (size_t j = 0; j < set::NY2NP; ++j)
        {
            res_p[i][j] = up[i][j] * vp[i][j];
        }
    }

    // return to Fourier space
    FOURIER_P2S_DEALIAS(res_p, res_);
}

void timeIntegration()
{
    for (size_t i = 0; i < set::NXHNP; i++)
    {
        for (size_t j = 0; j < set::NY; j++)
        {
            Complex residual = 3.0 * set::vor0[i][j] - 1.5 * set::vor1[i][j] + 1.0 / 3.0 * set::vor2[i][j] + set::dt * (3.0 * set::jac0[i][j] - 3.0 * set::jac1[i][j] + set::jac2[i][j]);
            double k2 = set::kx_local[i] * set::kx_local[i] + set::ky_local[j] * set::ky_local[j];
            set::vor[i][j] = residual / (11.0 / 6.0 + set::dt * (set::nu * k2));
        }
    }

    // solve streamwise function
    Complex phi[set::NXHNP][set::NY];
    Complex dvordx[set::NXHNP][set::NY];
    Complex dvordy[set::NXHNP][set::NY];
    Complex jac_x[set::NXHNP][set::NY];
    Complex jac_y[set::NXHNP][set::NY];
    solvePoisson_negsource(set::vor, phi);
    solveVelocity(phi);
    for (size_t i = 0; i < set::NXHNP; i++)
    {
        for (size_t j = 0; j < set::NY; j++)
        {
            dvordx[i][j] = set::vor[i][j] * set::imag_unit * set::kx_local[i];
            dvordy[i][j] = set::vor[i][j] * set::imag_unit * set::ky_local[j];
        }
    }
    nonlinear(set::u, dvordx, jac_x);
    nonlinear(set::v, dvordy, jac_y);
    for (size_t i = 0; i < set::NXHNP; i++)
    {
        for (size_t j = 0; j < set::NY; j++)
        {
            set::jac[i][j] = jac_x[i][j] + jac_y[i][j];
        }
    }

    for (size_t i = 0; i < set::NXHNP; i++)
    {
        for (size_t j = 0; j < set::NY; j++)
        {
            set::vor2[i][j] = set::vor1[i][j];
            set::vor1[i][j] = set::vor0[i][j];
            set::vor0[i][j] = set::vor[i][j];

            set::jac2[i][j] = set::jac1[i][j];
            set::jac1[i][j] = set::jac0[i][j];
            set::jac0[i][j] = set::jac[i][j];
        }
    }
}