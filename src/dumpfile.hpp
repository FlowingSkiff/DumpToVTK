#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <memory>
#include "geofile.hpp"
#include <unordered_map>
/**
 * @brief Dumpfile class stores the data from a LAMMPS dump
 * including atom properties, timestep, etc.
 * 
 */
class Dumpfile
{
public:
    // Constructor when there is no geofile
    Dumpfile(   const std::string& , uint32_t, uint32_t, 
                const std::array<std::array<double, 2>, 3>&,
                const std::vector<std::string>&); 
    // Constuctor with a geofile
    Dumpfile(   const std::string& , uint32_t, uint32_t, 
                const std::array<std::array<double, 2>, 3>&,
                const std::vector<std::string>&,
                const std::shared_ptr<const Geofile>&);
    Dumpfile() = delete;
    Dumpfile(const Dumpfile&) = delete;
    void SetSize(uint32_t, uint32_t); // Set the size of m_data
    void SetBoundary(const std::array<std::array<double, 2>, 3>&);
public: 
    friend std::istream& operator>>(std::istream&, Dumpfile&); // reading in data
    bool Write(); // can't be const, checks molecules outside of the domain before writing
private:
    void Wrap();
    void ParseProperties(const std::vector<std::string>&);
    void Scale();
    inline double Difference(const uint32_t&, const uint32_t&, const uint32_t&) const;
    inline double MaxDifference(const uint32_t&, const uint32_t&, const uint32_t&, const uint32_t&) const;
private:
    std::string m_name; // output filename
    std::vector<std::vector<double>> m_data; // atom data
    std::array<std::array<double, 2>, 3> m_bounds;
    std::shared_ptr<const Geofile> m_geofile;
private:
    struct AtomProperties
    {
        std::array<uint32_t, 3> xyz;
        std::unordered_map<std::string, std::array<uint32_t, 3>> vectors;
        std::unordered_map<std::string, uint32_t> scalars;
        bool isScaled;
    };
    AtomProperties m_atomProperties;
};