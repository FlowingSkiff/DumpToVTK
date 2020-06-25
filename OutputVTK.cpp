/*
  Written by Michael Roeing-Donna


*/

#pragma once
#include "OutputVTK.h"
#include<string>
#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<cmath>
namespace D{
  Dump::Dump(std::string filename_, std::string outputfilename_){
      filename = filename_;
      outputfilename = outputfilename_;
      geoflag = false;
      parseFile();
  }
  Dump::Dump(std::string filename_, std::string outputfilename_, std::string geofilename_){
      filename = filename_;
      outputfilename = outputfilename_;
      geofilename = geofilename_;
      geoflag = true;
      parseFile();
  }
  Dump::~Dump(){
      if (geoflag){
          for (int i = 0; i < geonumlines; i++)
              delete [] geolinedata[i];
          delete [] geolinedata;
      }
  }
  void Dump::setGeometryFile(std::string geofilename_){
      geoflag = true;
      geofilename = geofilename_;
  }
  std::string Dump::generateFileName(int step){
      char buffer [100];
      int formatsize = 8;
      sprintf(buffer, "%0*d", formatsize, step);
      std::string bufferstring(buffer);
      return outputfilename+bufferstring+std::string(".vtk");
  }
  void Dump::parseFile(){
      cols = {"id","mol","type","x ","y ","z ","xs","ys","zs","vx","vy","vz"};
      coli = std::vector<int> (cols.size(),0);
      scaleflag = false;
      std::ifstream dmp;
      dmp.open(filename);
      std::string line;
      int c(0), col(1), t(0);
      if (dmp.is_open()){
          while(getline(dmp, line)){
              c++;
              if (col==1){
                  if(line.find("ITEM: ATOMS")!=std::string::npos){
                      int ncol(0), ncolu(0);
                      for (std::string::iterator it = line.begin(); it!=line.end(); ++it){
                          if (*it==' ') ncol++;
                      }
                      ncol-=1;
                      if (line.back() == ' ') ncol-=1;
                      std::cout << "Located "<<ncol<<" per atom attributes\n";
                      for (unsigned int i=0; i< coli.size(); i++){
                          if (line.find(cols[i])!=std::string::npos){
                            ncolu++;
                            coli[i] = line.find(cols[i]);
                          }
                      }
                    numcols = ncolu;
                    col = 0;
                  }
            }
            if(line.find("ITEM: TIMESTEP")!=std::string::npos){
                t++;
            }
            if(line.find("ITEM: NUMBER OF ATOMS")!=std::string::npos){
                c++;
                getline(dmp,line);
                numatoms = stoi(line);
            }
          }
          dmp.close();
          numsteps = t;
          numlines = c;
      }
      std::vector<int> tempcoli = coli;
      c = 1;
      std::sort(tempcoli.begin(), tempcoli.end());
      for (std::vector<int>::iterator i = coli.begin(); i!=coli.end(); ++i){
          for (std::vector<int>::iterator j = tempcoli.begin(); j!=tempcoli.end(); ++j){
              if (*i==*j && *i!=0) *i = c++;
          }
      }

      if (coli[6] > 0){
          scaleflag = true;
          for (int i = 0; i < 3; i++){
              x[i] = coli[6]+i-1;
          }
      } else {
          for (int i = 0; i < 3; i++) x[i] = coli[3]+i-1;
      }
  }
  void Dump::WriteVTKs(){
      double **thisstep;
      thisstep = new double *[numatoms];
      for (int i = 0; i < numatoms; i++)
          thisstep[i] = new double[numcols];
      int thistimestep = 0;
      std::ifstream dumpfile;
      dumpfile.open(filename);
      std::string line;
      if (geoflag) GenerateGeoData();
      if (dumpfile.is_open()){
          while(getline(dumpfile, line)){
              if(line.find("ITEM: BOX")!=std::string::npos){
                  boundaryConditions[0] = line.at(17);
                  boundaryConditions[1] = line.at(20);
                  boundaryConditions[2] = line.at(23);
                  for(int i = 0; i < 3; i++){
                      getline(dumpfile, line);
                      sscanf(line.c_str(),"%le %le",&boxbounds[i][0],&boxbounds[i][1]);
                  }
              }
              else if (line.find("TIMESTEP")!=std::string::npos){
                  getline(dumpfile, line);
                  thistimestep = stoi(line);
              }
              else if (line.find("ITEM: ATOMS")!=std::string::npos){
                  std::cout << "num atoms = " << numatoms << ": num cols = " << numcols << std::endl;
                  for(int i = 0; i< numatoms; i++){
                      for( int j = 0; j < numcols; j++){
                          dumpfile >> thisstep[i][j];
                      }
                  }
                  VtkOut(thisstep, generateFileName(thistimestep));
              }
          }
      }
      dumpfile.close();
      for (int i = 0; i < numatoms; i++)
          delete [] thisstep[i];
      delete [] thisstep;
  }
  void Dump::WriteScalarToFile(FILE *&outfile, double **&data, int colIndex, std::string fieldName){
          fprintf(outfile, "SCALARS %s int 1\n", fieldName.c_str());
          fprintf(outfile, "LOOKUP_TABLE default\n");
          fprintf(outfile, "%i", int(data[0][colIndex]));
          for (int i = 1; i < numatoms; i++)
              fprintf(outfile, " %i", int(data[i][colIndex]));
          fprintf(outfile, "\n");
  } 
  void Dump::WriteVectorToFile(FILE *&outfile, double **&data, int * colIndex, std::string fieldName){
          fprintf(outfile, "VECTORS %s float\n",fieldName.c_str());
          fprintf(outfile, "%f %f %f", data[0][colIndex[0]],
                  data[0][colIndex[1]], data[0][colIndex[2]]);
          for (int i = 1; i < numatoms; i++){
              fprintf(outfile, "%f %f %f", data[i][colIndex[0]],
                      data[i][colIndex[1]], data[i][colIndex[2]]);
          }
  } 
  void Dump::VtkOut(double **data, std::string fname){
      FILE * outfile;
      outfile = fopen(fname.c_str(),"w");

      fprintf(outfile, "# vtk DataFile Version 2.0\n");
      fprintf(outfile, "From Lammps 2 vtk code\n");
      fprintf(outfile, "ASCII\n");
      fprintf(outfile, "DATASET POLYDATA\n");
      fprintf(outfile, "POINTS %i float\n",numatoms);
      AdjustForPeriodicity(data);
      for(int i = 0; i < numatoms; i++){
          if (!scaleflag){
              double tmpx(data[i][coli[3]-1]), tmpy(data[i][coli[4]-1]), tmpz(data[i][coli[5]-1]);
              fprintf(outfile,"%f %f %f\n",tmpx,tmpy,tmpz);
          } else {
              fprintf(outfile, "%f %f %f\n", 
                      data[i][coli[6]-1]*(boxbounds[0][1]-boxbounds[0][0])+boxbounds[0][0],
                      data[i][coli[7]-1]*(boxbounds[1][1]-boxbounds[1][0])+boxbounds[1][0],
                      data[i][coli[8]-1]*(boxbounds[2][1]-boxbounds[2][0])+boxbounds[2][0]);
          }
      }
      if (geoflag){
          if (geonumlines > 0){
              fprintf(outfile, "LINES %i %i\n",geonumlines, geonumlines*3);
              for (int i = 0; i < geonumlines; i++){
                  fprintf(outfile, "%i %i %i\n",2,geolinedata[i][2]-1,geolinedata[i][3]-1);
              }
          }
          if (geonumtris > 0){
              fprintf(outfile, "POLYGONS %i %i\n",geonumtris,geonumtris*4);
              for (int i = 0; i < geonumtris; i++){
                  fprintf(outfile, "%i %i %i %i\n",3,geotridata[i][2]-1,geotridata[i][3]-1,geotridata[i][4]-1);
              }
          }
      }
      fprintf(outfile, "VERTICES %i %i\n", numatoms, 2*numatoms);
      for (int i = 0; i < numatoms; i++)
          fprintf(outfile, "%i %i\n",1,i);

      fprintf(outfile, "POINT_DATA %i\n", numatoms);
      if (coli[2]){
          fprintf(outfile, "SCALARS atom_type int 1\n");
          fprintf(outfile, "LOOKUP_TABLE default\n");
          fprintf(outfile, "%i", int(data[0][coli[2]-1]));
          for (int i = 1; i < numatoms; i++)
              fprintf(outfile, " %i", int(data[i][coli[2]-1]));
          fprintf(outfile, "\n");
      }
      if (coli[1]){
          fprintf(outfile, "SCALARS mol int 1\n");
          fprintf(outfile, "LOOKUP_TABLE default\n");
          fprintf(outfile, "%i", int(data[0][coli[1]-1]));
          for (int i = 1; i < numatoms; i++)
              fprintf(outfile, " %i", int(data[i][coli[1]-1]));
          fprintf(outfile, "\n");
      }
      if (coli[0]){
          fprintf(outfile, "SCALARS id int 1\n");
          fprintf(outfile, "LOOKUP_TABLE default\n");
          fprintf(outfile, "%i", int(data[0][coli[0]-1]));
          for(int i = 1; i < numatoms; i++)
              fprintf(outfile, " %i", int(data[i][coli[0]-1]));
          fprintf(outfile, "\n");
      }
      if (coli[9]){
          fprintf(outfile, "VECTORS vel float\n");
          fprintf(outfile, "%f %f %f", data[0][coli[9]-1],
                  data[0][coli[10]-1], data[0][coli[11]-1]);
          for (int i = 1; i < numatoms; i++){
              fprintf(outfile, " %f %f %f", data[i][coli[9]-1],
                      data[i][coli[10]-1], data[i][coli[11]-1]);
          }
      }
      fclose(outfile);
      }

 void Dump::GenerateGeoData(){
  if (!geoflag) return;
  geonumlines = 0;
  geonumtris = 0;
  std::ifstream geofile;
  std::string line;
  geofile.open(geofilename);
  if(geofile.is_open()){
    while(getline(geofile,line)){
      if (line.find("bonds")!=std::string::npos){
        geonumlines = stoi(line.substr(0,line.size()-6));
        geolinedata = new int * [geonumlines];
        for (int i = 0; i < geonumlines; i++)
            geolinedata[i] = new int [4];
      }
      if (line.find("angles")!=std::string::npos){
          geonumtris = stoi(line.substr(0,line.size()-7));
          geotridata = new int * [geonumtris];
          for (int i = 0; i < geonumtris; i++)
              geotridata[i] = new int [5];
      }
      if (line.find("Bonds")!=std::string::npos){
          getline(geofile, line);
          for (int i = 0; i < geonumlines; i++)
              for (int j = 0; j < 4; j++)
                  geofile >> geolinedata[i][j];
      }
      if (line.find("Angles")!=std::string::npos){
          getline(geofile, line);
          for (int i = 0; i < geonumtris; i++)
              for (int j = 0; j < 5; j++)
                  geofile >> geotridata[i][j];

      }
    }
    geofile.close();
    std::cout << "GeoLines = " << geonumlines << std::endl;
    std::cout << "GeoTris = " << geonumtris << std::endl;
  }
  
 }
	void Dump::BuildMolData( double **&data){

	}
  void Dump::AdjustForPeriodicity(double **&data){
      if (!geoflag || !IsPeriodic()) return;
      if (coli[0]==0) std::cout << "Error: Atom id must be in dump file if geometry file is used\n";
      double periodShift[3];
      periodShift[0] = boxbounds[0][1]-boxbounds[0][0];
      periodShift[1] = boxbounds[1][1]-boxbounds[1][0];
      periodShift[2] = boxbounds[2][1]-boxbounds[2][0];
      if (scaleflag){
          periodShift[0]=periodShift[1]=periodShift[2]=1;
      }
      for (int i = 0; i < geonumlines; i++){ // Loop over all bonds chekcing for periodicity
          int point0(geolinedata[i][2]-1), point1(geolinedata[i][3]-1);
          for (int j = 0; j < 3; j++){
              double component0(data[point0][x[j]]), component1(data[point1][x[j]]);
              if ((component1-component0)>periodShift[j]/2){
                  data[point1][x[j]]-=periodShift[j];
                  i = -1;
              }
              else if ((component0-component1)>periodShift[j]/2) {
                  data[point0][x[j]]-=periodShift[j];
                  i = -1;
              }
          }
      }
  } 
  bool Dump::IsPeriodic(){
      for (int i = 0; i < 3; i++){
          if (boundaryConditions[i] == 'p') return true;
      }
      return false;
  }  
};
