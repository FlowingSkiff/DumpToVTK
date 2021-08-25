#pragma once
#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>
#include <array>
#include <list>
/**
 * @brief Gofile class stores the geography information 
 * from the input molecule file. Only one of these should be 
 * created/needed and a shared ptr for the data should be 
 * used to share the Geofile in multi-threaded applications.
 * 
 */
class Geofile
{ 
private:
    using Bond = std::array<uint32_t, 2>;
    using Triangle = std::array<uint32_t, 3>;
    enum Attribute{ID, MOL, TYPE, X, Y, Z};
public:
    Geofile(const std::string&); // filename;
    Geofile(); // user must call read
    Geofile(const Geofile& copy) = delete;
    Geofile(const Geofile&& move) = delete;
    friend std::istream& operator>>(std::istream& , Geofile&);
    bool Read(const std::string&);
    bool IsOkay() const;
    const std::vector<Bond>& GetBonds() const;
    const std::vector<Triangle>& GetTriangles() const;
    const std::unordered_map<uint32_t, std::vector<uint32_t>>& GetMoleculeMap() const;
    const std::unordered_map<uint32_t, std::list<Bond>>& GetMoleculeBondMap() const;
    const std::unordered_map<uint32_t, std::list<Triangle>>& GetMoleculeAngleMap() const;
private:
    // Molecule map for each molecule id
    // Molecule is a vector of atom id's 
    bool m_status;
    std::unordered_map<uint32_t, std::vector<uint32_t>> m_map;
    std::vector<Bond> m_bonds;
    std::vector<Triangle> m_triangles;
    std::unordered_map<uint32_t, std::list<Bond>> m_molBondMap;
    std::unordered_map<uint32_t, std::list<Triangle>> m_molAngleMap;
};