#pragma once
#include <string>
#include <vector>

/**
 * @brief The properties class should be a static class
 * which handles the command line options passed to the 
 * program. 
 * 
 */
class Properties
{
public:
    Properties(const Properties&) = delete;
    static Properties& Get();
public:
    std::string GetInputFile() const;
    std::string GetOutputFile() const;
    std::string GetMoleculeFile() const;
    void SetInputFile(const std::string&);
    void SetOutputFile(const std::string&);
    void SetMoleculeFile(const std::string&);
    void SetVerbose(const bool&);
    void SetGenerateTriangles(const bool&);
    void SetGenerateLines(const bool&);
    bool SetMultithreaded(const bool&, const int&);
    // Checks the parameters which are stored in the 
    // properties class
    int Check() const;
    bool GetProcessGeoFile() const;
    bool GetVerbose() const;
    bool GetMultithreaded() const;
    uint32_t GetNumberOfThreads() const;
    bool GenerateLines() const;
    bool GenerateTris() const;
public:
    // Flag class for status of input options
    enum{ OKAY, BAD_INPUT = 0x1, BAD_OUTPUT = 0x10, BAD_MOLECULE = 0x100};
private:
    Properties();
    void SetProcessGeoFile(const bool&);
private:
    static Properties m_instance;
    bool m_verbose;
    std::string m_moleculefile;
    std::string m_inputfile;
    std::string m_outputfile;
    bool m_generatetriangles;
    bool m_generatelines;
    bool m_processgeofile;
    bool m_multithreaded;
    uint32_t m_numthreads;
    std::vector<std::string> m_errorcodes;
};