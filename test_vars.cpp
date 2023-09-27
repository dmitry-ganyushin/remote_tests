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
 * Created 9/27/2023
 *      Author: Dmitry Ganyushin ganyushin@gmail.com
 */
#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>

#include <adios2.h>
#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

enum test_cases
{
    DIM1,
    DIM3
};
void read1D(int nproc, int rank, const std::string &filename, const int NSTEPS, adios2::ADIOS &adios, std::vector<std::string> &variables, std::vector<size_t> &start, std::vector<size_t> &count);
/* test 2 and 3
 * read 1 or many 1D variables
 */
void read3D(int nproc, int rank, const std::string &filename, const int NSTEPS, adios2::ADIOS &adios, std::vector<std::string> &variables, std::vector<size_t> &start, std::vector<size_t> &count);
/* test 5
 * A 3D subset from 3D variable
 */
void read1D(int nproc, int rank, const std::string &filename, const int NSTEPS, adios2::ADIOS &adios, std::vector<std::string> &variables_in, size_t start, size_t count)
{
    unsigned int startX = start;
    unsigned int countX = count;

    try
    {
        std::chrono::time_point<std::chrono::system_clock> start;
        std::chrono::time_point<std::chrono::system_clock> end;
        start = std::chrono::system_clock::now();
        adios2::IO io = adios.DeclareIO("Input");
       // io.SetEngine("BP4");

        adios2::Engine reader = io.Open(filename, adios2::Mode::Read);

        /* get variables with 1D shape */


        for (size_t step = 0; step < NSTEPS; step++) {
            reader.BeginStep();
            auto variables = io.AvailableVariables(true);

            for (auto const &var_name: variables) {
                adios2::Variable<double> var =
                        io.InquireVariable<double>(var_name.first);
                if (var.Shape().size() == 1 ){
                    auto globalSize = var.Shape()[0];
                    auto localSize = globalSize / nproc;
                    std::cout << var.Shape()[0] << std::endl;
                    startX = localSize * rank;
                    countX = localSize;
                    std::vector<double> data1D(countX);
                    var.SetSelection(adios2::Box<adios2::Dims>({startX},
                                                               {countX}));
                    reader.Get<double>(var, data1D[0]);
                }
            }

            reader.EndStep();
        }

        reader.Close();
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << elapsed_seconds.count() << std::endl;
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

void read3D(int nproc, int rank, const std::string &filename, const int NSTEPS, adios2::ADIOS &adios, std::vector<std::string> &variables, std::vector<size_t> &start, std::vector<size_t> &count)
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
        //io.SetEngine("BP4");

        /*
         * Define global array: type, name, global dimensions
         * The local process' part (start, count) can be defined now or later
         * before Write().
         */

        adios2::Engine reader = io.Open(filename, adios2::Mode::Read);
        auto variables = io.AvailableVariables(true);

        for (size_t step = 0; step < NSTEPS; step++)
        {
            reader.BeginStep();
            for (auto const& var_name : variables)
            {
                adios2::Variable<double> var =
                        io.InquireVariable<double>(var_name.first);
                if (var.Shape().size() == 3 ) {
                    auto globalSizeX = var.Shape()[0];
                    auto globalSizeY = var.Shape()[1];
                    auto globalSizeZ = var.Shape()[2];
                    auto globalSize = globalSizeX * globalSizeY * globalSizeZ;
                    auto localSize = globalSize / nproc;
                    std::cout << var.Shape()[0] << std::endl;

                    countX = globalSizeX / nproc;
                    startX = countX * rank;
                    if (rank == nproc - 1)
                    {
                        // last process need to read all the rest of slices
                        countX = globalSizeX - countX * (nproc - 1);
                    }

                    var.SetSelection(adios2::Box<adios2::Dims>(
                            {startX, 0, 0}, {countX, globalSizeY, globalSizeZ}));
                    reader.Get<double>(var, data3D[0][0][0]);
                }
            }

            reader.EndStep();
        }
        reader.Close();
        end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_time - start_time;
        std::cout << elapsed_seconds.count() << std::endl;
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
    const unsigned int Nx = 0;
    const unsigned int Ny = 0;
    const unsigned int Nz = 0;
    const unsigned int countX = 0;
    const unsigned int countY = 0;
    const unsigned int countZ = 0;
    //options
    std::string filename;
    int mode;

    //should be adjusted for getopt
    auto start = std::vector<size_t>(3);
    auto count = std::vector<size_t>(3);
    std::vector<std::string> variables;


    option longopts[] = {
                {"help", no_argument, NULL, 'h'},
                {"case", required_argument, NULL, 'c'},
                {"filename", required_argument, NULL, 'f'},
                {0,0,0,0}
    };

        while (1) {
            int option_index = 0;
            const int opt = getopt_long(argc, argv, "hc:f:", longopts, &option_index);

            if (opt == -1) {
                break;
            }

            switch (opt) {
                case 0:
                    /* If this option set a flag, do nothing else now. */
                    if (longopts[option_index].flag != 0)
                        break;
                    printf ("option %s", longopts[option_index].name);
                    if (optarg)
                        printf (" with arg %s", optarg);
                    printf ("\n");
                    break;
                case 'h':
                    std::cout << "help" << std::endl;
                    break;
                case 'c':
                    if(strcmp("1D", optarg) == 0) {
                        mode = DIM1;
                    }
                    if(strcmp("3D", optarg) == 0) {
                        mode = DIM3;
                    }
                case 'f':
                    if (strlen(optarg) > 0){
                        filename = optarg;
                    }

                    break;
                default:
                    break;
            }
        }

    switch (mode) {
            case (DIM1):
                read1D(nproc, rank, filename, NSTEPS, adios, variables, start[0], count[0]);
            break;
        case (DIM3):
            read3D(nproc, rank, filename, NSTEPS, adios, variables, start, count);
            break;
            default:
                break;
    }




#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}


