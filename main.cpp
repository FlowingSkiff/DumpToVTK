#include <iostream>
#include <string>
#include "properties.hpp"
#include "outputvtks.hpp"

Properties Properties::m_instance;

static void PrintHelp()
{
    using std::cout;
    cout << "Usage is: LammpsToVTK [OPTION] -i [INPUTFILE] -o [OUTPUTFILE]\n"
            "Converts a lammps dump file to a vtk format and writes a seperate file for each\n"
            "timestep. \n"
            "\n"
            "Options include:\n"
            "  -h,    Displays this help message.\n"
            "  -i,    Specifies the next argument is the input file.\n"
            "  -l,    Specifies that lines should be generated from the molecule file.\n"
            "  -m,    Specifies the molecule filename from this directory. The molecule\n"
            "         file must have bond or angle definitions in order for the mapping\n"
            "         to work correctly. It is best to use the input file specified to\n"
            "         LAMMPS for the specific output file. Bonds are treated as lines \n"
            "         and angles are treated as triangle polygons.\n"
            "  -o,    Specifies that the next argument is the outputfile string.\n"
            "  -t,    Enables triangle generation from the molecule file.\n"
            "  -v,    Enables verbose messaging to the console.\n"
            "  -n,    Specifies the maximum number of threads used, default: 1\n"
            " \n"
            "Exit status is always 0 when everything runs fine,  and 1 in the event of\n"
            "an error.\n";
}

/**
 * @brief ParseArguments takes the command line input and
 * ensures that the options are correct before starting to 
 * parse the files. Encountered errors are printed to the 
 * console.
 * 
 * @param argc Number of input arguments
 * @param argv string array of input arguments
 * @return true Should exit, encountered error
 * @return false Argument parsed successfully
 */
static bool ParseArguments(const int& argc, char** argv)
{
    bool shouldExit = false;
    // Minimum is be 5 arguments, 
    // 0            1   2       3   4
    //LammpsToVtk   -i  input   -o  output
    if (argc < 5)
    {
        std::cout << "Not enough input arguments\n";
        PrintHelp();
        return true;
    }
    for (int i = 1; i < argc; ++i)
    {
        std::string tmp(argv[i]);
        if (tmp[0] == '-')
        {
            for (const char& c : tmp.substr(1))
            {
                switch (c)
                {
                    case 'v': // verbose
                        Properties::Get().SetVerbose(true);
                        break;
                    case 'i': // input file name
                        if (i == argc - 1)
                        {
                            shouldExit = true;
                            break;
                        }
                        Properties::Get().SetInputFile(argv[++i]);
                        break;
                    case 'o': // ouptut file string
                        if (i == argc - 1)
                        {
                            shouldExit = true;
                            break;
                        }
                        Properties::Get().SetOutputFile(argv[++i]);
                        break;
                    case 'm': // molecule file name
                        if (i == argc - 1)
                        {
                            shouldExit = true;
                            break;
                        }
                        Properties::Get().SetMoleculeFile(argv[++i]);
                        break;
                    case 't': // generate triangles
                        Properties::Get().SetGenerateTriangles(true);
                        break;
                    case 'l': // generate lines
                        Properties::Get().SetGenerateLines(true);
                        break;
                    case 'h': // print help
                        PrintHelp();
                        shouldExit = true;
                        break;
                    case 'n': // number of threads
                        {
                            if (i == argc - 1)
                            {
                                shouldExit = true;
                                break;
                            }
                            if (Properties::Get().SetMultithreaded(true, std::stoi(argv[++i])))
                            {
                                shouldExit = true;
                            }
                            break;
                        }
                    default:
                        shouldExit = true;
                        break;
                }
            }
        }
        else shouldExit = true;
    }
    if (Properties::Get().Check())
        shouldExit = true;
    return shouldExit;
}

int main(int argc, char ** argv){
    if (ParseArguments(argc, argv)) 
        return 1;
    OutputVTKs();
    return 0;
}
