#include <iomanip>
#include <sstream>
#include <iostream>
#include <hdf5.h>
#include "Pgrm_setup.hpp"

hid_t FILE_CREATE(const std::string &IO_FILENAME)
{
    // HDF5 variables
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    if (plist_id < 0)
    {
        return -1;
    }

    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

    // Create the file collectively
    hid_t file_id = H5Fcreate(IO_FILENAME.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

    // Close property list
    H5Pclose(plist_id);

    return file_id;
}

void WRITE_2D(hid_t file_id, const std::string &dataset_name, const double (&data)[set::NXHNP][set::NY])
{
    int rank = 2;
    hsize_t dims[2] = {(hsize_t)set::NXH, (hsize_t)set::NY};
    hsize_t count[2] = {(hsize_t)set::NXHNP, (hsize_t)set::NY};
    hsize_t offset[2] = {(hsize_t)set::NXHNP * set::rank, 0};

    hid_t file_space_id = H5Screate_simple(rank, dims, NULL);
    hid_t dset_id = H5Dcreate(file_id, dataset_name.c_str(), H5T_NATIVE_DOUBLE, file_space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(file_space_id);

    hid_t mem_space_id = H5Screate_simple(rank, count, NULL);
    file_space_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, offset, NULL, count, NULL);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, mem_space_id, file_space_id, plist_id, data[0]);

    H5Pclose(plist_id);
    H5Sclose(mem_space_id);
    H5Sclose(file_space_id);
    H5Dclose(dset_id);
}

void save_Instantaneous(const int id_step_)
{
    std::stringstream filename_stream;
    filename_stream << "./DATA/F-" << std::setw(8) << std::setfill('0') << id_step_ << ".H5";
    std::string filename = filename_stream.str();
    if (set::rank == 0)
    {
        std::cout << "step " << set::id_step << ", saving data to file: " << filename << std::endl;
    }

    // IO init
    H5open();

    hid_t file_id = FILE_CREATE(filename);
    if (file_id < 0)
    {
        return;
    }

    double var_real[set::NXHNP][set::NY];
    double var_imag[set::NXHNP][set::NY];

    // save velocity u
    for (size_t i = 0; i < set::NXHNP; i++)
    {
        for (size_t j = 0; j < set::NY; j++)
        {
            var_real[i][j] = set::u[i][j].real();
            var_imag[i][j] = set::u[i][j].imag();
        }
    }
    WRITE_2D(file_id, "u-real", var_real);
    WRITE_2D(file_id, "u-imag", var_imag);

    // save velocity v
    for (size_t i = 0; i < set::NXHNP; i++)
    {
        for (size_t j = 0; j < set::NY; j++)
        {
            var_real[i][j] = set::v[i][j].real();
            var_imag[i][j] = set::v[i][j].imag();
        }
    }
    WRITE_2D(file_id, "v-real", var_real);
    WRITE_2D(file_id, "v-imag", var_imag);

    // save vorticity
    for (size_t i = 0; i < set::NXHNP; i++)
    {
        for (size_t j = 0; j < set::NY; j++)
        {
            var_real[i][j] = set::vor[i][j].real();
            var_imag[i][j] = set::vor[i][j].imag();
        }
    }
    WRITE_2D(file_id, "omega-real", var_real);
    WRITE_2D(file_id, "omega-imag", var_imag);

    // IO finish
    H5Fclose(file_id);
    H5close();
}