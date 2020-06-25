#pragma once
#include<vector>
#include<string>
namespace D{
class Dump{
    public:
        Dump(std::string, std::string);
        Dump(std::string, std::string, std::string);
        ~Dump();
        void WriteVTKs();
    private:
        std::string generateFileName(int);
        std::string filename;
        std::string outputfilename;
        bool geoflag, scaleflag;
        std::string geofilename;
        int numsteps, numcols, numatoms, formatsize, numlines;
        std::vector<std::string>cols;
        std::vector<int> coli;
        int x[3];
        char boundaryConditions[3];
        int geonumlines, geonumtris;
        int ** geolinedata;
        int ** geotridata;
        void GenerateGeoData();
        double boxbounds[3][2];
        void setGeometryFile(std::string);
        void parseFile();
        void VtkOut(double** , std::string);
        bool IsPeriodic();
        void AdjustForPeriodicity(double **&);
        void WriteScalarToFile(FILE *&, double **&, int, std::string);
        void WriteVectorToFile(FILE *&, double **&, int *, std::string);
		void BuildMolData(double **&);
};
};
