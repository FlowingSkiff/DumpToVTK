#include "geofile.hpp"
#include <fstream>
#include <regex>
#include "properties.hpp"
Geofile::Geofile(const std::string& filename):
    m_status(false)
{
    Read(filename);
}

Geofile::Geofile(): m_status(false) { }

bool Geofile::Read(const std::string& filename)
{
    std::ifstream input(filename, std::ifstream::in);
    if (!input.is_open())
    {
        m_status = false;
        return false;
    }
    const std::regex expression("^([0-9]*) ([a-zA-Z\\s]*)");
    std::string line;

    uint32_t num_triangles = 0;
    uint32_t num_bonds = 0;
    uint32_t num_atoms = 0;
    std::vector<std::vector<double>> atomData;
    while (std::getline(input, line))
    {
        std::smatch stringmatch;
        if (std::regex_match(line, stringmatch, expression))
        {
            if (stringmatch[2] == "atoms")
            {
                num_atoms = std::stoi(stringmatch[1]);
                atomData.resize(num_atoms, std::vector<double>(6, 0));
            }
            else if (stringmatch[2] == "bonds")
            {
                num_bonds = std::stoi(stringmatch[1]);
            }
            else if (stringmatch[2] == "angles")
            {
                num_triangles = std::stoi(stringmatch[1]);
            }
        }
        else if (line == "Atoms")
        {
            if (num_atoms == 0)
            {
                m_status = false;
                return false;
            }
            for (uint32_t atom = 0; atom < num_atoms; ++atom)
            {
                for (auto& column : atomData.at(atom))
                {
                    input >> column;
                }
                m_map[static_cast<uint32_t>(atomData.at(atom).at(Attribute::MOL))].push_back(static_cast<uint32_t>(atomData.at(atom).at(Attribute::ID)));
            }
        }
        else if (line == "Bonds")
        {
            if (num_bonds == 0) 
            {
                m_status = false;
                return false;
            }
            m_bonds.reserve(num_bonds);
            for (uint32_t bond = 0; bond < num_bonds; ++bond)
            {
                uint32_t id, type, firstID, secondID;
                input >> id >> type >> firstID >> secondID;
                m_bonds.push_back({firstID, secondID});
            }
        }
        else if (line == "Angles")
        {
            if (num_triangles == 0)
            {   
                m_status = false;
                return false;
            }
            m_triangles.reserve(num_triangles);
            for (uint32_t angle = 0; angle < num_triangles; ++angle)
            {
                uint32_t id, type, firstID, secondID, thirdID;
                input >> id >> type >> firstID >> secondID >> thirdID;
                m_triangles.push_back({firstID, secondID, thirdID});
            }
        }
    }
    if (Properties::Get().GetVerbose())
    {
        std::cout << "Geography file analysis: \n";
        std::cout << "Found " << m_map.size() << " molecules\n";
        std::cout << "Found " << m_bonds.size() << " bonds\n";
        std::cout << "Found " << m_triangles.size() << " angles\n";
    }
    for (const auto& bond : m_bonds)
    {
        uint32_t mol = static_cast<uint32_t>(atomData.at(bond.at(0)).at(Attribute::MOL));
        m_molBondMap[mol].emplace_back(bond);
    }
    for (const auto& angle : m_triangles)
    {
        uint32_t mol = static_cast<uint32_t>(atomData.at(angle.at(0)).at(Attribute::MOL));
        m_molAngleMap[mol].emplace_back(angle);
    }
    return m_status = true;;
}

bool Geofile::IsOkay() const
{
    return m_status;
}

// returns the bond data, should be returned as a const reference
const std::vector<Geofile::Bond>& Geofile::GetBonds() const
{
    return m_bonds;
}
// returns the triangle data, should be returned as a const reference
const std::vector<Geofile::Triangle>& Geofile::GetTriangles() const
{
    return m_triangles;
}


// wraps the atom data using the map generated when reading in
const std::unordered_map<uint32_t, std::vector<uint32_t>>& Geofile::GetMoleculeMap() const
{
    return m_map;
}

// gets the molecule map based on bond information
const std::unordered_map<uint32_t, std::list<Geofile::Bond>>& Geofile::GetMoleculeBondMap() const
{
    return m_molBondMap;
}

// gets the molecule triangle map 
const std::unordered_map<uint32_t, std::list<Geofile::Triangle>>& Geofile::GetMoleculeAngleMap() const
{
    return m_molAngleMap;
}