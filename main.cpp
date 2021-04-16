#include<iostream>
#include"OutputVTK.cpp"
#include<string.h>
using namespace std;
using namespace D;

void DisplayHelp(int reason){
    if (reason == 0){
        cout << "--------------------------------------------------------------------\n";
        cout << "Help with LammpsToVtk\n";
        cout << "Writen by Michael R-D\n\n";
        cout << "Usage:   ./LammpsToVTK inputfile outputfilestring geometryfile\n\n";
        cout << "Note:    Geometry file is optional.\n";
        cout << "Note:    Dump file should be sorted by id for use with geometry file\n";
        cout << "Note:    Can handle scaled and unscale dump files\n";
        cout << "Note:    Supports position, velocity, mol, id, and atom type\n";
        cout << "Note:    Timestep and vtk suffix is automatically appended\n";
        cout << "Note:    For outputfilestring = dumpvtk, file on timestep 0 \n";
        cout << "         would be dumpvtk00000000.vtk\n";
        cout << "--------------------------------------------------------------------\n";
    } else if (reason == 1){
        cout << "Wrong number of input arguments. 2 or 3 required\n";
        cout << "use --help for more information\n";
    }
}


int main(int argc, char ** argv){
	if (argc < 2){
		DisplayHelp(1);
		return 0;
	} else if (strcmp(argv[1],"--help")==0){
        DisplayHelp(0);
        return 0;
    } else if (argc < 3){
        DisplayHelp(1);
        return 0;
    }
    else if (argc == 3) {
        Dump dmpfile(argv[1], argv[2]);
        dmpfile.WriteVTKs();
    } else if (argc == 4){
        Dump dmpfile(argv[1], argv[2], argv[3]);
        dmpfile.WriteVTKs();
    } else {
        DisplayHelp(1);
        return 0;
    }
    return 0;
}
