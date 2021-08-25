#include "dumpfile.hpp"
#include "assert.h"
#include <fstream>
#include <algorithm>
#include "properties.hpp"
#include <cmath>
#include <limits>
Dumpfile::Dumpfile( const std::string& str, uint32_t n_rows, uint32_t n_columns,
                    const std::array<std::array<double, 2>, 3>& bounds, 
                    const std::vector<std::string>& atomProperties):
                        m_name(str), 
                        m_data(n_rows, std::vector<double>(n_columns, 0.00)), 
                        m_bounds(bounds), 
                        m_geofile(nullptr) 
{ 
    ParseProperties(atomProperties);
    if (Properties::Get().GetVerbose())
    {
        std::cout << "Dump file atom scalars are: \n";
        for (const auto& p : m_atomProperties.scalars)
        {
            std::cout << p.first << ' ';
        }
        std::cout << '\n' << "Dump file atom vectors are: \n";
        for (const auto& p : m_atomProperties.vectors)
            std::cout << p.first << ' ';
        std::cout << '\n';
    }
}

Dumpfile::Dumpfile( const std::string& str, uint32_t n_rows, uint32_t n_columns, 
                    const std::array<std::array<double, 2>, 3>& bounds,
                    const std::vector<std::string>& atomProperties, 
                    const std::shared_ptr<const Geofile>& geofileptr): 
                        m_name(str), 
                        m_data(n_rows, std::vector<double>(n_columns, 0.00)),
                        m_bounds(bounds), 
                        m_geofile(geofileptr) 
{ 
    ParseProperties(atomProperties);
}


// Reading in the atom information from the input file
// The size of m_data needs to be specified before calling
// this function
std::istream& operator>>(std::istream& input, Dumpfile& dmp)
{
    if (dmp.m_data.empty() || dmp.m_data[0].empty()) 
    {
        std::cout << "Error handling dump file. \n";
        assert(false);
    }
    for (auto& v : dmp.m_data)
        for (double& i : v)
            input >> i;
    return input;
}


// Resizes m_data for atom storage
// rows is the number of atoms in the simulation
// columns are the number of properties per atom in the 
//   dump file
void Dumpfile::SetSize(uint32_t n_rows, uint32_t n_cols)
{
    m_data.resize(n_rows, std::vector<double>(n_cols, 0));
}


// Sets the domain boundary, it should be taken in as a 
// ctor parameter but here for legacy
void Dumpfile::SetBoundary(const std::array<std::array<double, 2>, 3>& boundary)
{
    m_bounds = std::move(boundary);
}


// Writes the dumpfile into VTK format
// First wraps and scales the atoms in the domain
bool Dumpfile::Write()
{
    if (Properties::Get().GetVerbose())
    {
        std::cout << "Writing vtk: " << m_name << '\n';
    }
    if (m_atomProperties.isScaled)
        Scale();
    if (m_geofile.get())
        Wrap();
    std::ofstream outfile(m_name, std::ios::out);
    if (!outfile.is_open())
    {
        std::cout << "Could not open file.\n";
        return false;
    }
    outfile << "# vtk DataFile Version 2.0\n";
    outfile << "From Lammps 2 vtk code\n";
    outfile << "ASCII\n";
    outfile << "DATASET POLYDATA\n";
    outfile << "POINTS " << m_data.size() << " FLOAT\n";
    for (const auto& atom : m_data)
    {
        for (const auto& x : m_atomProperties.xyz)
        {
            outfile << atom[x] << ' ';
        }
        outfile << '\n';
    }
    
    outfile << "VERTICES " << m_data.size() << ' ' << m_data.size() * 2 << '\n';
    for (uint32_t i = 0; i < m_data.size(); ++i)
        outfile << "1 " << i << '\n';
    if (Properties::Get().GenerateLines())
    {
        const auto& bonds = m_geofile->GetBonds();
        uint32_t numLines = bonds.size();
        outfile << "LINES " << numLines << " " << 3 * numLines <<'\n';
        for (const auto& b : bonds)
        {
            outfile << "2 " << b[0]-1 << " " << b[1]-1 << '\n';
        }
    }
    if (Properties::Get().GenerateTris())
    {
        const auto& tris = m_geofile->GetTriangles();
        uint32_t numTris = tris.size();
        outfile << "POLYGONS " << numTris << " " << 4*numTris << '\n';
        for (const auto& t : tris)
        {
            outfile << "3 " << t[0] << ' ' << t[1] << ' ' << t[2] << '\n';
        }

    }
    if (m_atomProperties.scalars.size() > 0 || m_atomProperties.vectors.size() > 0)
        outfile << "POINT_DATA " << m_data.size() << '\n';
    for (const auto& p : m_atomProperties.scalars)
    {
        outfile << "SCALARS " << p.first << " float 1\n";
        outfile << "LOOKUP_TABLE default\n";
        for (const auto& atom : m_data)
        {
            outfile << atom[p.second] << '\n';
        }
    }
    for (const auto& p : m_atomProperties.vectors)
    {
        outfile << "VECTORS " << p.first << " float\n";
        for (const auto& atom : m_data)
        {
            for (const auto& i : p.second)
            {
                outfile << atom[i] << ' ';
            }
            outfile << '\n';
        }
    }
    if (Properties::Get().GetVerbose())
    {
        std::cout << "Finished Writing vtk: " << m_name << '\n';
    }
    return true;
}


// Parses the atom properties inside the dump file. Scalars
// are generated by finding properties which do not have 
// x, y, or z at the end of the property and vice versa.
void Dumpfile::ParseProperties(const std::vector<std::string>& props)
{
    for (uint32_t i = 0; i < props.size(); ++i)
    {
        if (props[i] == "x" || props[i] == "xs")
            m_atomProperties.xyz[0] = i;
        else if (props[i] == "y" || props[i] == "ys")
            m_atomProperties.xyz[1] = i;
        else if (props[i] == "z" || props[i] == "zs")
            m_atomProperties.xyz[2] = i;
    }
    m_atomProperties.isScaled = props[m_atomProperties.xyz[0]] == "xs";

    // search through the rest of the properties
    // vectors will end in x, y, z and scalars do not
    for (uint32_t i = 0; i < props.size(); ++i)
    {
        // Check to make sure index is not x, y, z
        if (std::find(std::begin(m_atomProperties.xyz), std::end(m_atomProperties.xyz) , i) == std::end(m_atomProperties.xyz))
        {
            const std::string& prop = props[i];
            if (prop.back() == 'x' || prop.back() == 'y' || prop.back() == 'z') // vector
            {
                uint32_t tmpindex = prop.back() - 'x';
                std::string tmpvalue = prop.substr(0, prop.size() - 1);
                if (m_atomProperties.vectors.find(tmpvalue) == m_atomProperties.vectors.end())
                    m_atomProperties.vectors[tmpvalue] = {0u, 0u, 0u};
                m_atomProperties.vectors[tmpvalue][tmpindex] = i;
                if (Properties::Get().GetVerbose())
                {
                    std::cout << "Found vector: " << prop.substr(0, prop.size() - 1) << '\n';
                }
            }
            else  // scalar
            {
                m_atomProperties.scalars[prop] = i;
                if (Properties::Get().GetVerbose())
                {
                    std::cout << "Found scalar: " << prop << '\n';
                }
            }
        }
    }
}

// Scales the position of the atom to domain scale rather than 
// only having values between 0 and 1. This is done using the 
// box bounds passed to the dumpfile
void Dumpfile::Scale()
{
    // compute the multiplier using the domain bounds
    std::array<double,3> scalefactor { 0.00, 0.00, 0.00};
    for (uint16_t i = 0; i < 3u; ++i)
    {
        scalefactor.at(i) = m_bounds.at(i).at(1) - m_bounds.at(i).at(0);
    }
    // reshape each atom using the multiplier and adding the lower bound
    for (auto& atom : m_data)
    {
        for (uint16_t i = 0; i < 3u; ++i)
        {
            const auto& x = m_atomProperties.xyz.at(i);
            atom.at(x) = atom.at(x) * scalefactor.at(i) + m_bounds.at(i).at(0);
        }
    }
}

void Dumpfile::Wrap()
{
    if (Properties::Get().GetVerbose())
    {
        std::cout << "Wraping atoms.\n";
    }
    bool bondFlag =  Properties::Get().GenerateLines();
    bool angleFlag = Properties::Get().GenerateTris();
    
    // For each molecule, a map for each bond/tri holding each
    // atom id is needed. A molecule is wrapped to the +x, +y, +z
    // side if the length between two atoms in a bond/tri is 
    // greater than half the total domain.

    const auto& bondMap = m_geofile->GetMoleculeBondMap();
    const auto& angleMap = m_geofile->GetMoleculeAngleMap();
    if (bondFlag && bondMap.size() == 0)
        bondFlag = false;
    if (angleFlag && angleMap.size() == 0)
        angleFlag = false;
    if (!bondFlag && !angleFlag)
        return;

    std::array<double,3> shift { 0.0, 0.0, 0.0}, maxLength{ 0, 0, 0};
    for (uint32_t i = 0; i < 3u; ++i)
    {
        shift[i] = m_bounds[i][1] - m_bounds[i][0];
        maxLength[i] = shift[i] / 2.0;
    }
    if (bondFlag)
    {
        for (const auto& mol : bondMap)
        {
            if (mol.first == 0) continue;
            std::array<bool, 3> shouldShift{{false, false, false}};
            bool shouldBreak = false;
            for (const auto& bond : mol.second)
            {
                for (uint32_t x = 0u; x < 3u; ++x)
                {
                    uint32_t index = m_atomProperties.xyz.at(x);
                    uint32_t atom0 = bond.at(0) - 1; // Atom id starts at 1
                    uint32_t atom1 = bond.at(1) - 1; // atom id starts at 1
                    if (Difference(atom0, atom1, index) > maxLength.at(x))
                    {
                        shouldShift.at(x) = true;
                        if (std::all_of(shouldShift.begin(), shouldShift.end(), [](const auto& b){return b;}))
                        {
                            shouldBreak = true;
                        }
                    }
                    if (shouldBreak) break;
                }
            }
            if (std::any_of(shouldShift.begin(), shouldShift.end(), [](const auto& b){ return b;}))
            {
                for (const auto& bond : mol.second)
                {
                    for (const auto& atom : bond)
                    {
                        for (uint32_t x = 0u; x < 3u && shouldShift.at(x); ++x)
                        {
                            uint32_t index = m_atomProperties.xyz.at(x);
                            if (m_data.at(atom - 1).at(index) < shift.at(x))
                            {
                                m_data.at(atom - 1).at(index) += shift.at(x);
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        for (const auto& mol : angleMap)
        {
            if (mol.first == 0) continue;
            std::array<bool, 3> shouldShift{{false, false, false}};
            bool shouldBreak = false;
            for (const auto& angle : mol.second)
            {
                for (uint32_t x = 0u; x < 3u; ++x)
                {
                    uint32_t index = m_atomProperties.xyz.at(x);
                    uint32_t atom0 = angle.at(0) - 1; // atom id starts at 1
                    uint32_t atom1 = angle.at(1) - 1; // atom id starts at 1
                    uint32_t atom2 = angle.at(2) - 1; // atom id starts at 1
                    if (MaxDifference(atom0, atom1, atom2, index) > maxLength.at(x))
                    {
                        shouldShift.at(x) = true;
                        if (std::all_of(shouldShift.begin(), shouldShift.end(), [](const auto& b){return b;}))
                        {
                            shouldBreak = true;
                        }
                    }
                    if (shouldBreak) break;
                }
            }
            if (std::any_of(shouldShift.begin(), shouldShift.end(), [](const auto& b){ return b;}))
            {
                for (const auto& angle : mol.second)
                {
                    for (const auto& atom : angle)
                    {
                        for (uint32_t x = 0u; x < 3u && shouldShift.at(x); ++x)
                        {
                            uint32_t index = m_atomProperties.xyz.at(x);
                            if (m_data.at(atom - 1).at(index) < shift.at(x))
                            {
                                m_data.at(atom - 1).at(index) += shift.at(x);
                            }
                        }
                    }
                }
            }
        }
    }
    if (Properties::Get().GetVerbose())
    {
        std::cout << "End wraping atoms.\n";
    }
}

inline double Dumpfile::Difference(const uint32_t& a0, const uint32_t& a1, const uint32_t& index) const
{
    return std::abs(m_data.at(a0).at(index) - m_data.at(a1).at(index));
}

inline double Dumpfile::MaxDifference(  const uint32_t& a0, 
                                        const uint32_t& a1, 
                                        const uint32_t& a2, 
                                        const uint32_t& index) const
{
    return std::max({Difference(a0, a1, index), Difference(a0, a2, index), Difference(a1, a2, index)});
}