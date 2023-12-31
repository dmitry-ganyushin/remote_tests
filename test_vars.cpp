/*
 *  Test program for remote reading
 * Created 9/27/2023
 *      Author: Dmitry Ganyushin ganyushin@gmail.com
 */
#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>

#include <adios2.h>

#if ADIOS2_USE_MPI
#include <mpi.h>
#endif
bool DEBUG = false;
enum test_cases
{
    DIM1,
    DIM3X,
    DIM3Y,
    DIM3Z,
    DIM3PLANEYZ,
    DIM3PLANEXY,
    DIM3PLANEXZ
};
void read1D(int nproc, int rank, const std::string &filename, const int NSTEPS, adios2::IO &io, std::vector<std::string> &variables, int output_length);
/* test 2 and 3
 * read 1 or many 1D variables
 */
void read3D(int nproc, int rank, const std::string &filename, const int NSTEPS, adios2::IO &io, std::vector<std::string> &variables, int direction, double ratio, int output_line_length);
/* test 5
 * A 3D subset from 3D variable
 */

void read3DPlane(int nproc, int rank, const std::string &filename, const int NSTEPS, adios2::IO &io, std::vector<std::string> &variables, int direction, double ratio, int output_line_length);

std::string &getOutputString(int rank, int output_line_length, std::string &out, std::vector<double> &data);

/* test 5
 * A 3D subset from 3D variable
 */

void read1D(int nproc, int rank, const std::string &filename, const int NSTEPS, adios2::IO &io, std::vector<std::string> &variables_in, double ratio, int output_line_length)
{
    unsigned int startX;
    unsigned int countX;

    try
    {
        std::chrono::time_point<std::chrono::system_clock> start;
        std::chrono::time_point<std::chrono::system_clock> end;
        start = std::chrono::system_clock::now();


        adios2::Engine reader = io.Open(filename, adios2::Mode::Read);

        /* get variables with 1D shape */

        std::string out;
        for (size_t step = 0; step < NSTEPS; step++) {
            reader.BeginStep();
            auto variables = io.AvailableVariables(true);

            for (auto const &var_name: variables) {
                if(!variables_in.empty() && var_name.first != variables_in[0] ) continue;
                adios2::Variable<double> var =
                        io.InquireVariable<double>(var_name.first);
                if (var.Shape().size() == 1 ){
                    auto globalSize = var.Shape()[0];
                    auto localSize = globalSize / nproc;
                    startX = localSize * rank;
                    countX = localSize;
                    if (rank == nproc - 1) {
                        // last process need to read all the rest of slices
                        countX = globalSize - countX * (nproc - 1);
                    }
                    std::vector<double> data1D(countX);
                    var.SetSelection(adios2::Box<adios2::Dims>({startX},
                                                               {countX}));
                    reader.Get<double>(var, data1D.data(), adios2::Mode::Sync);

                    if (DEBUG)
                        out = getOutputString(rank, output_line_length, out, data1D);
                }
            }

            reader.EndStep();
        }

        reader.Close();
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        out += "Duration :" + std::to_string(elapsed_seconds.count()) + "\n";
        std::ofstream out_file;
        out_file.open ("output1D-" + std::to_string(rank) + ".log");
        out_file << out;
        out_file.close();
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

std::string &getOutputString(int rank, int output_line_length, std::string &out, std::vector<double> &data) {
    {
        out += "rank " + std::to_string(rank) + ":" + "\n";
        int counter = 0;

        for (auto v: data)
        {
            counter++;
            out += std::to_string(v) + " " ;
            if (counter == output_line_length){
                out += "\n";
                counter = 0;
            }
        }
    }
    return out;
}

void read3D(int nproc, int rank, const std::string &filename, const int NSTEPS, adios2::IO &io, std::vector<std::string> &variables_in, int direction, double ratio, int output_line_length)
{
    unsigned int startX;
    unsigned int startY;
    unsigned int startZ;
    unsigned int countX;
    unsigned int countY;
    unsigned int countZ;

    try
    {
        std::chrono::time_point<std::chrono::system_clock> start_time;
        std::chrono::time_point<std::chrono::system_clock> end_time;
        start_time = std::chrono::system_clock::now();
        // Get io settings from the config file or
        // create one with default settings here

        /*
         * Define global array: type, name, global dimensions
         * The local process' part (start, count) can be defined now or later
         * before Write().
         */

        adios2::Engine reader = io.Open(filename, adios2::Mode::Read);

        std::string out;
        for (size_t step = 0; step < NSTEPS; step++)
        {
            reader.BeginStep();
            auto variables = io.AvailableVariables(true);
            for (auto const& var_name : variables)
            {
                if(!variables_in.empty() && var_name.first != variables_in[0] ) continue;
                adios2::Variable<double> var =
                        io.InquireVariable<double>(var_name.first);
                if (var.Shape().size() == 3 ) {
                    auto globalSizeX = var.Shape()[0];
                    auto globalSizeY = var.Shape()[1];
                    auto globalSizeZ = var.Shape()[2];

                    if (ratio == 1.0){
                        startX = 0;
                        countX = globalSizeX;

                        startY = 0;
                        countY = globalSizeY;

                        startZ = 0;
                        countZ = globalSizeZ;
                    }
                    else{
                        startX = (1 - ratio)/2.0 * globalSizeX;
                        countX = globalSizeX * ratio;

                        startY = (1 - ratio)/2.0 * globalSizeY;
                        countY = globalSizeY * ratio;

                        startZ = (1 - ratio)/2.0 * globalSizeZ;
                        countZ = globalSizeZ * ratio;
                    }

                   switch(direction) {
                       case DIM3X:
                           countX = globalSizeX / nproc;
                           startX = countX * rank;

                           if (rank == nproc - 1) {
                               // last process need to read all the rest of slices
                               countX = globalSizeX - countX * (nproc - 1);
                           }
                           break;
                       case DIM3Y:
                           countY = globalSizeY / nproc;
                           startY = countY * rank;

                           if (rank == nproc - 1) {
                               // last process need to read all the rest of slices
                               countY = globalSizeY - countY * (nproc - 1);
                           }
                           break;
                       case DIM3Z:
                           countZ = globalSizeZ / nproc;
                           startZ = countZ * rank;

                           if (rank == nproc - 1) {
                               // last process need to read all the rest of slices
                               countZ = globalSizeZ - countZ * (nproc - 1);
                           }
                           break;
                   }
                    var.SetSelection(adios2::Box<adios2::Dims>(
                            {startX, startY, startZ}, {countX, countY, countZ}));
                    size_t elementsSize = var.SelectionSize();
                    std::vector<double> data3D(elementsSize);
                    reader.Get<double>(var, data3D.data(), adios2::Mode::Sync);
                    if (DEBUG)
                        out = getOutputString(rank, output_line_length, out, data3D);
                }
            }

            reader.EndStep();
        }
        reader.Close();
        end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_time - start_time;
        out += "Duration: " + std::to_string(elapsed_seconds.count()) + "\n";
        std::ofstream out_file;
        out_file.open ("output3D-" + std::to_string(rank) + ".log");
        out_file << out;
        out_file.close();
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

void read3DPlane(int nproc, int rank, const std::string &filename, const int NSTEPS, adios2::IO &io, std::vector<std::string> &variables_in, int direction, double ratio, int output_line_length)
{

    unsigned int startX;
    unsigned int startY;
    unsigned int startZ;
    unsigned int countX;
    unsigned int countY;
    unsigned int countZ;

    try
    {
        std::chrono::time_point<std::chrono::system_clock> start_time;
        std::chrono::time_point<std::chrono::system_clock> end_time;
        start_time = std::chrono::system_clock::now();

        adios2::Engine reader = io.Open(filename, adios2::Mode::Read);
        std::string out;
        for (size_t step = 0; step < NSTEPS; step++)
        {
            reader.BeginStep();
            auto variables = io.AvailableVariables(true);
            for (auto const& var_name : variables)
            {
                if(!variables_in.empty() && var_name.first != variables_in[0] ) continue;
                adios2::Variable<double> var =
                        io.InquireVariable<double>(var_name.first);
                if (var.Shape().size() == 3 ) {
                    auto globalSizeX = var.Shape()[0];
                    auto globalSizeY = var.Shape()[1];
                    auto globalSizeZ = var.Shape()[2];
                    //read some plane in the middle
                    switch(direction) {
                        case DIM3PLANEXZ:
                            countX = globalSizeX / nproc;
                            startX = countX * rank;

                            if (rank == nproc - 1) {
                                // last process need to read all the rest of slices
                                countX = globalSizeX - countX * (nproc - 1);
                            }
                            countY = 1;
                            startY = globalSizeY/2;

                            if (ratio == 1.0){
                                startZ = 0;
                                countZ = globalSizeZ;
                            }
                            else{
                                startZ = (1 - ratio)/2.0 * globalSizeZ;
                                countZ = globalSizeZ * ratio;
                            }



                            break;
                        case DIM3PLANEYZ:
                            countY = globalSizeY / nproc;
                            startY = countY * rank;
                            if (rank == nproc - 1) {
                                // last process need to read all the rest of slices
                                countY = globalSizeY - countY * (nproc - 1);
                            }
                            countX = 1;
                            startX = globalSizeX/2;

                            if (ratio == 1.0){
                                startZ = 0;
                                countZ = globalSizeZ;
                            }
                            else{
                                startZ = (1 - ratio)/2.0 * globalSizeZ;
                                countZ = globalSizeZ * ratio;
                            }

                            break;

                        case DIM3PLANEXY:
                            countY = globalSizeY / nproc;
                            startY = countY * rank;

                            if (rank == nproc - 1) {
                                // last process need to read all the rest of slices
                                countY = globalSizeY - countY * (nproc - 1);
                            }

                            countZ = 1;
                            startZ = globalSizeZ / 2;

                            if (ratio == 1.0){
                                startX = 0;
                                countX = globalSizeX;
                            }
                            else{
                                startX = (1 - ratio)/2.0 * globalSizeX;
                                countX = globalSizeX * ratio;
                            }

                            break;
                        default:
                            break;
                    }
                    var.SetSelection(adios2::Box<adios2::Dims>(
                            {startX, startY, startZ}, {countX, countY, countZ}));

                    size_t elementsSize = var.SelectionSize();
                    std::vector<double> data2D(elementsSize);
                    reader.Get<double>(var, data2D.data(), adios2::Mode::Sync);
                }
            }

            reader.EndStep();
        }
        reader.Close();
        end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_time - start_time;
        out += "Duration: " + std::to_string(elapsed_seconds.count()) + "\n";
        std::ofstream out_file;
        out_file.open ("output3DPlane-" + std::to_string(rank) + ".log");
        out_file << out;
        out_file.close();
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
#if ADIOS2_USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#endif
    int NSTEPS = 1;

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    //options
    std::string filename;
    std::string engine = "BP5";
    std::string  transport = "filesystem";
    int mode = -1;
    double ratio = 1.0;
    int output_line_length = 10;

    //should be adjusted for getopt
    auto start = std::vector<size_t>(3);
    auto count = std::vector<size_t>(3);
    std::vector<std::string> variables;


    option longopts[] = {
                {"help", no_argument, NULL, 'h'},
                {"case", required_argument, NULL, 'c'},
                {"filename", required_argument, NULL, 'f'},
                {"engine", required_argument, NULL, 'e'},
                {"transport", required_argument, NULL, 't'},
                {"variables", required_argument, NULL, 'v'},
                {"ratio", required_argument, NULL, 'r'},
                {"debug", no_argument, NULL, 'd'},
                {"length", required_argument, NULL, 'l'},
                {"nsteps", required_argument, NULL, 'n'},
                {0,0,0,0}
    };

        while (1) {
            int option_index = 0;
            const int opt = getopt_long(argc, argv, "hc:f:e:t:v:dl:n:", longopts, &option_index);

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
                    std::cout << "Help: --case: 1D, 3DX, 3DY, 3DZ, 3DYZ, 3DXY, 3DXZ" << std::endl;
                    std::cout << "Usage: mpirun -n 8 ./test_vars  --case 3DX --variables var1 --filename /absolute/path/on/remote/machine/remote.bp" << std::endl;
                    break;
                case 'd':
                    DEBUG = true;
                    break;
                case 'c':
                    if(strcmp("1D", optarg) == 0) {
                        mode = DIM1;
                    }
                    if(strcmp("3DX", optarg) == 0) {
                        mode = DIM3X;
                    }
                    if(strcmp("3DY", optarg) == 0) {
                        mode = DIM3Y;
                    }
                    if(strcmp("3DZ", optarg) == 0) {
                        mode = DIM3Z;
                    }

                    if(strcmp("3DYZ", optarg) == 0) {
                        mode = DIM3PLANEYZ;
                    }
                    if(strcmp("3DXZ", optarg) == 0) {
                        mode = DIM3PLANEXZ;
                    }
                    if(strcmp("3DXY", optarg) == 0) {
                        mode = DIM3PLANEXY;
                    }
                    break;
                case 'f':
                    if (strlen(optarg) > 0){
                        filename = optarg;
                    }
                    break;
                case 'e':
                    if (strlen(optarg) > 0){
                        engine = optarg;
                    }
                    break;
                case 't':
                    if (strlen(optarg) > 0){
                        transport = optarg;
                    }
                    break;
                case 'v':
                    if (strlen(optarg) > 0){
                        std::string variables_string = optarg;
                        //split and put into vector
                        variables.push_back(variables_string);
                    }
                    break;
                case 'r':
                    if (strlen(optarg) > 0){
                        ratio = std::stod(optarg);
                    }
                    break;
                case 'l':
                    if (strlen(optarg) > 0){
                        output_line_length = std::stoi(optarg);
                    }
                    break;
                case 'n':
                    if (strlen(optarg) > 0){
                        NSTEPS = std::stoi(optarg);
                    }
                    break;
                default:
                    break;
            }
        }

    adios2::IO io = adios.DeclareIO("Input");
    if (engine == "BP5")
        io.SetEngine("BP5");
    else if (engine == "BP4")
        io.SetEngine("BP4");
    if (transport == "http" || transport == "HTTP" ){
        std::map<std::string, std::string> parameters;
        parameters["Library"] = "HTTP";
        io.AddTransport("File", parameters);
    }
    else if (transport == "daos" || transport == "DAOS"){
        std::map<std::string, std::string> parameters;
        parameters["Library"] = "DAOS";
        io.AddTransport("File", parameters);
    }


    switch (mode) {
        case DIM1:
            read1D(nproc, rank, filename, NSTEPS, io, variables, ratio, output_line_length);
            break;
        case DIM3X:
            read3D(nproc, rank, filename, NSTEPS, io, variables, DIM3X, ratio, output_line_length);
            break;
        case DIM3Y:
            read3D(nproc, rank, filename, NSTEPS, io, variables, DIM3Y, ratio, output_line_length);
            break;
        case DIM3Z:
            read3D(nproc, rank, filename, NSTEPS, io, variables, DIM3Z, ratio, output_line_length);
            break;
        case DIM3PLANEYZ:
            read3DPlane(nproc, rank, filename, NSTEPS, io, variables, DIM3PLANEYZ, ratio, output_line_length);
            break;
        case DIM3PLANEXY:
            read3DPlane(nproc, rank, filename, NSTEPS, io, variables, DIM3PLANEXY, ratio, output_line_length);
            break;
        case DIM3PLANEXZ:
            read3DPlane(nproc, rank, filename, NSTEPS, io, variables, DIM3PLANEXZ, ratio, output_line_length);
            break;

        default:
            break;
    }

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}


