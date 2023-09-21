/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Write a global array from multiple processors.
 *
 * A global array is an N-dimensional array. A process can write a sub-array
 * into the global array by stating the N-dimensional offset and the size of
 * the sub-array. At reading, one can read back any portion of the array
 * regardless of how many processors wrote that data.
 *
 * Processes are NOT required
 * - to stay in the boundaries of the global dimensions. However, one will not
 * be able to read back data outside of the boundaries.
 * - to fill the whole global array, i.e. one can leave holes in it. At reading,
 * one will get the fill-value set for the array for those coordinates that
 * are not written by any process.
 *
 * The global dimensions of a global array MUST NOT change over time.
 * If they are, then the array should be handled as a local array. Of course, if
 * only a single output step is written to a file, that still shows up at
 * reading as a global array.
 *
 * The decomposition of the array across the processes, however, can change
 * between output steps.
 *
 * Created on: Jun 2, 2017
 *      Author: pnorbert
 */
#include <getopt.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>

#include <adios2.h>
#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

void read1D(int rank, const std::string &filename, const int NSTEPS, adios2::ADIOS &adios, std::vector<size_t> &start, std::vector<size_t> &count);
/* test 2 and 3
 * read 1 or many 1D variables
 */
void read3D(int rank, const std::string &filename, const int NSTEPS, adios2::ADIOS &adios, std::vector<size_t> &start, std::vector<size_t> &count);
/* test 5
 * A 3D subset from 3D variable
 */
void read1D(int rank, const std::string &filename, const int NSTEPS, adios2::ADIOS &adios, size_t start, size_t count)
{
    unsigned int startX = start;
    unsigned int countX = count;
    std::vector<double> data1D(countX);

    try
    {
        std::chrono::time_point<std::chrono::system_clock> start;
        std::chrono::time_point<std::chrono::system_clock> end;
        start = std::chrono::system_clock::now();
        adios2::IO io = adios.DeclareIO("Input");
        io.SetEngine("BP4");

        adios2::Engine reader = io.Open(filename, adios2::Mode::Read);

        std::vector<std::string> variables;
        /* get variables with 1D shape */
        variables.emplace_back("var1");

        for (size_t step = 0; step < NSTEPS; step++) {
            reader.BeginStep();
            for (auto const &name: variables) {
                adios2::Variable<double> var =
                        io.InquireVariable<double>(name);

                var.SetSelection(adios2::Box<adios2::Dims>({startX},
                                                           {countX}));
                reader.Get<double>(var, data1D[0]);
            }

            reader.EndStep();
        }

        reader.Close();
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << "finished reading at " << elapsed_seconds.count() << "s\n";
    }
    catch (std::invalid_argument &e)
    {
        if (rank == 0)
        {
            std::cout << "Invalid argument exception, STOPPING PROGRAM\n";
            std::cout << e.what() << "\n";
        }
    }
    catch (std::ios_base::failure &e)
    {
        if (rank == 0)
        {
            std::cout << "System exception, STOPPING PROGRAM\n";
            std::cout << e.what() << "\n";
        }
    }
    catch (std::exception &e)
    {
        if (rank == 0)
        {
            std::cout << "Exception, STOPPING PROGRAM\n";
            std::cout << e.what() << "\n";
        }
    }
}

void read3D(int rank, const std::string &filename, const int NSTEPS, adios2::ADIOS &adios, std::vector<size_t> &start, std::vector<size_t> &count)
{

    unsigned int startX = start[0];
    unsigned int startY = start[1];
    unsigned int startZ = start[2];
    unsigned int countX = count[0];
    unsigned int countY = count[1];
    unsigned int countZ = count[2];
    std::vector<std::vector<std::vector<double>>> data3D(countX, std::vector<std::vector<double>>(countY, std::vector<double>(countZ)));

    try
    {
        std::chrono::time_point<std::chrono::system_clock> start_time;
        std::chrono::time_point<std::chrono::system_clock> end_time;
        start_time = std::chrono::system_clock::now();
        // Get io settings from the config file or
        // create one with default settings here
        adios2::IO io = adios.DeclareIO("Output");
        io.SetEngine("BP4");

        /*
         * Define global array: type, name, global dimensions
         * The local process' part (start, count) can be defined now or later
         * before Write().
         */

        adios2::Engine reader = io.Open(filename, adios2::Mode::Read);

        std::vector<std::string> variables;
        variables.emplace_back("var1");
        for (size_t step = 0; step < NSTEPS; step++)
        {
            reader.BeginStep();
            for (auto const& name : variables)
            {
                adios2::Variable<double> var =
                        io.InquireVariable<double>(name);

                var.SetSelection(adios2::Box<adios2::Dims>({start[0], start[1], start[2]},
                                                            {count[0], count[1], count[2]}));
                reader.Get<double>(var, data3D[0][0][0]);
            }

            reader.EndStep();
        }
        reader.Close();
        end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_time - start_time;
        std::cout << "finished reading at " << elapsed_seconds.count() << "s\n";
    }
    catch (std::invalid_argument &e)
    {
        if (rank == 0)
        {
            std::cout << "Invalid argument exception, STOPPING PROGRAM\n";
            std::cout << e.what() << "\n";
        }
    }
    catch (std::ios_base::failure &e)
    {
        if (rank == 0)
        {
            std::cout << "System exception, STOPPING PROGRAM\n";
            std::cout << e.what() << "\n";
        }
    }
    catch (std::exception &e)
    {
        if (rank == 0)
        {
            std::cout << "Exception, STOPPING PROGRAM\n";
            std::cout << e.what() << "\n";
        }
    }
}

int main(int argc, char *argv[])
{

    int rank = 0;
    int nproc;
    int nprocx = 2;
    int nprocy = 2;
    int nprocz = 2;
#if ADIOS2_USE_MPI
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#endif
    const int NSTEPS = 1;
    const int GLOBAL_SIZE = 128;
    const int LOCAL_SIZE = 64;

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    // Application variables for output get parameters from bp file
    const unsigned int Nx = 128;
    const unsigned int Ny = 128;
    const unsigned int Nz = 128;
    const unsigned int countX = Nx/nprocx;
    const unsigned int countY = Ny/nprocy;
    const unsigned int countZ = Nz/nprocz;
    const std::string filename = "remote_reading.bp";

    //should be adjusted
    auto start = std::vector<size_t>(3);
    auto count = std::vector<size_t>(3);


    option longopts[] = {
                {"help", no_argument, NULL, 'h'},
                {"case", required_argument, NULL, 'c'},
                {"filename", optional_argument, NULL, 'f'}, {0}};

        while (1) {
            const int opt = getopt_long(argc, argv, "hcf::", longopts, 0);

            if (opt == -1) {
                break;
            }

            switch (opt) {
                case 'h':
                    std::cout << "help" << std::endl;
                    break;
                case 'c':
                    read1D(rank, filename, NSTEPS, adios, start[0], count[0]);
                    read3D(rank, filename, NSTEPS, adios, start, count);
                    break;
                default:
                    break;
            }
        }

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}


