#include "outputvtks.hpp"
#include "properties.hpp"
#include "dumpfile.hpp"
#include "geofile.hpp"
#include <fstream>
#include <thread>
#include <memory>
#include <vector>
#include <regex>
#include <numeric>
#include <sstream>
#include "wrapped.hpp"
#include <deque>
#include <list>

// Creates a filename string using the base string and current timestep
// hard coded to use 8 digits to represent the timestep
static std::string CreateFilename(const std::string& base, const std::string& timestep)
{
    return base + std::string(8-timestep.size(), '0') + timestep + ".vtk";
}

/**
 * @brief Reads in the input file from LAMMPS and writes out
 * the data in vtk format
 * 
 * @return true Always return true right now
 */
bool OutputVTKs()
{
    Properties& properties = Properties::Get();
    std::deque<std::thread> threads;
    Wrapped<uint32_t> threadIndex(0u, properties.GetNumberOfThreads());

    // Geofile handling
    std::shared_ptr<Geofile> geofile(nullptr);
    std::thread* geofilethread = nullptr;
    bool geofileHasJoined = (properties.GetMultithreaded()) ? false : true;
    if (properties.GetProcessGeoFile())
    {
        if (properties.GetVerbose())
        {
            std::cout << "Reading geography file.\n";
        }
        geofile = std::shared_ptr<Geofile>(new Geofile);
        if (properties.GetMultithreaded())
        {
            threads.emplace_back(&Geofile::Read, std::ref(*geofile),  properties.GetMoleculeFile());
            geofilethread = &threads.back();
            ++threadIndex;
        }
        else
        {
            geofile->Read(properties.GetMoleculeFile());
        }
    }

    // load in and start working with the dump file
    std::list<Dumpfile> dumpfiles; 
    std::ifstream inputfile(Properties::Get().GetInputFile());
    const std::regex itemExpression("ITEM: ([A-Z\\sa-z]*)");
    std::string line, timestep;
    uint32_t numAtoms, numProperties;
    std::array<std::array<double, 2>, 3> bounds;
    std::vector<std::string> atomProperties;
    while (std::getline(inputfile, line))
    {
        std::smatch match;
        if (std::regex_match(line, match, itemExpression))
        {
            if (match[1] == "TIMESTEP")
            {
                inputfile >> timestep;
                if (properties.GetVerbose())
                {
                    std::cout << "Reading timestep: " << timestep << '\n';
                }
            }
            else if (match[1] == "NUMBER OF ATOMS")
            {
                inputfile >> numAtoms;
                if (properties.GetVerbose())
                {
                    std::cout << "Found " << numAtoms << '\n';
                }
            }
            else if (match[1].str().find("BOX BOUNDS") != std::string::npos)
            {
                if (properties.GetVerbose())
                {
                    std::cout << "Bounds are: ";
                }
                for (auto& i : bounds)
                {
                    for (auto& v : i)
                    {
                        inputfile >> v;
                        if (properties.GetVerbose())
                            std::cout << v << ' ';
                    }
                }
                if (properties.GetVerbose())
                {
                    std::cout << '\n';
                }
            }
            else if (match[1].str().find("ATOMS") != std::string::npos)
            {
                // match is: 'ATOMS id mol type x y z vx vy vz'
                // columns is the same number of spaces
                if (properties.GetVerbose())
                {
                    std::cout << "Reading atoms: \n";
                }
                if (atomProperties.size() == 0)
                {
                    std::istringstream stream(match[1].str());
                    std::string tmpproperty;
                    stream >> tmpproperty;
                    while (stream >> tmpproperty)
                    {
                        atomProperties.emplace_back(tmpproperty);
                    }
                    numProperties = atomProperties.size();
                    if (properties.GetVerbose())
                    {
                        std::cout << "Found " << numProperties << ": ";
                        for (const auto& v : atomProperties)
                            std::cout << v << ' ';
                        std::cout << '\n';
                    }
                }
                if (geofile)
                {
                    if (!geofileHasJoined && geofilethread)
                    {
                        if (geofilethread->joinable())
                            geofilethread->join();
                        geofileHasJoined = true;
                    }
                    dumpfiles.emplace_back( CreateFilename(properties.GetOutputFile(), timestep),
                                            numAtoms, numProperties, bounds, atomProperties, geofile);
                }
                else
                {
                    dumpfiles.emplace_back( CreateFilename(properties.GetOutputFile(), timestep),
                                            numAtoms, numProperties, bounds, atomProperties);
                }
                if (properties.GetVerbose())
                {
                    std::cout << "Reading atom information: \n";
                }
                inputfile >> dumpfiles.back();
                if (properties.GetMultithreaded())
                {
                    if (threads.size() == properties.GetNumberOfThreads())
                    {
                        threads.at(threadIndex.GetValue()).join();
                        threads.at(threadIndex.GetValue()) = std::move(std::thread(&Dumpfile::Write, std::ref(dumpfiles.back())));
                        ++threadIndex;
                    }
                    else
                    {
                        threads.emplace_back(&Dumpfile::Write, std::ref(dumpfiles.back()));
                        ++threadIndex;
                    }
                }
                else
                {
                    // if not using multi-threaded, the memory can be freed right away
                    dumpfiles.back().Write();
                    dumpfiles.pop_back();
                }
            } // ATOM regex
        } // REGEX match
    } // Line in file

    if (properties.GetMultithreaded())
    {
        for (auto& t : threads)
            t.join();
    }
    // TODO: Add in error handling
    return true;
}