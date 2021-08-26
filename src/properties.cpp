#include "properties.hpp"
#include <iostream>
#include <fstream>
#include <thread>
Properties::Properties():
    m_verbose(false), m_moleculefile(),
    m_inputfile(), m_outputfile(), m_generatetriangles(false),
    m_generatelines(false), m_processgeofile(false), m_multithreaded(false),
    m_numthreads(1)
{ 
    m_errorcodes = {
        "No error.\n",
        "Received bad input file.\n",
        "Received bad output file.\n",
        "Received bad molecule file.\n"
    };
}

Properties& Properties::Get()
{
    return m_instance;
}

std::string Properties::GetInputFile() const
{
    return m_inputfile;
}
std::string Properties::GetOutputFile() const
{
    return m_outputfile;
}
std::string Properties::GetMoleculeFile() const
{
    return m_moleculefile;
}
void Properties::SetInputFile(const std::string& str)
{
    m_inputfile = str;
}
void Properties::SetOutputFile(const std::string& str)
{
    m_outputfile = str;
}
void Properties::SetMoleculeFile(const std::string& str)
{
    m_moleculefile = str;
}
void Properties::SetVerbose(const bool& b)
{
    m_verbose = b;
}
void Properties::SetGenerateTriangles(const bool& b)
{
    m_generatetriangles = b;
}
void Properties::SetGenerateLines(const bool& b)
{
    m_generatelines = b;
}

int Properties::Check() const
{
    int status = OKAY;
    std::ifstream i_file(m_inputfile, std::ios::in);
    if (!i_file.is_open()) 
    {
        status |= BAD_INPUT;
        std::cout << m_errorcodes[1];
    }
    i_file.close();
    
    i_file.open(m_outputfile, std::ios::in);
    if (i_file.is_open())
    {
        status |= BAD_OUTPUT;
        std::cout << m_errorcodes[2];
    }
    i_file.close();

    if (m_generatelines || m_generatetriangles)
    {
        Get().SetProcessGeoFile(true);
        i_file.open(m_moleculefile, std::ios::in);
        if (!i_file.is_open())
        {
            status |= BAD_MOLECULE;
            std::cout << m_errorcodes[3];
            Get().SetProcessGeoFile(false);
        }
        i_file.close();
    }
    return status != OKAY;
}

void Properties::SetProcessGeoFile(const bool& b)
{
    m_processgeofile = b;
}

bool Properties::GetProcessGeoFile() const
{
    return m_generatelines || m_generatetriangles;
}

bool Properties::GetVerbose() const
{
    return m_verbose;
}

bool Properties::GetMultithreaded() const
{
    return m_multithreaded;
}

bool Properties::SetMultithreaded(const bool& val, const int& numThreads)
{
    m_multithreaded = val;
    m_numthreads = static_cast<uint32_t>(numThreads);
    return !(numThreads > 0);
}

uint32_t Properties::GetNumberOfThreads() const
{
    return m_numthreads;
}


bool Properties::GenerateLines() const
{
    return m_generatelines;
}
bool Properties::GenerateTris() const
{
    return m_generatetriangles;
}