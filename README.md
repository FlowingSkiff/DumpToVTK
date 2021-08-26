# LammpsToVTK
LammpsToVTK converts a LAMMPS dump file to vtk format. It can take in a 2D or 3D file with an arbitrary amount of other dumped atom properties. 
The read atom properties are assumed to follow the standard LAMMPS format where _%x, %y, %z_ is treated as a vector quantity and all other properties are scalars. 
Reading is done single threaded, but multi-threaded writing is supported. 

Lines and triangles can optionally be created by passing an input file where atoms, bonds, and angles are defined. 
Bonds are used to create lines and angles are used to create triangles. 
When a bond or angle has atoms which are wrapped about a periodic boundary, atom positioned on the lower boundary are shifted above the upper boundary. 
_This unwrapping functionality only occurs when a geography file is used._

## Installation
```bash
git clone https://github.com/FlowingSkiff/LammpsToVTK.git
cd LammpsToVTK
mkdir build
cd build
cmake ..
make
```

## Use

Usage is: LammpsToVTK [OPTION] -i [INPUTFILE] -o [OUTPUTFILE]

Options include:
| Arg | Next Arg Type | Description |
| --- | --- | --- |
| -i | __STRING_PATH__ | Specifies the next argument is the input file |
| -o | __STRING_PATH_PREFIX__ | Specifies that the next argument is the outputfile string |
| -m | __STRING_PATH__ | Specifies the molecule filename from this directory. The molecule file must have bond or angle definitions in order for the mapping to work correctly. It is best to use the input file specified to LAMMPS for the specific output file. Bonds are treated as lines and angles are treated as triangle polygons |
| -l | __NONE__ | Specifies that lines should be generated from the molecule file |
| -t | __NONE__ | Enables triangle generation from the molecule file |
| -n | __INTEGER__ | Specifies the maximum number of threads used, defaults to single-threaded |
| -v | __NONE__ | Enables verbose messaging to the console |
| -h | __NONE__ | Displays the help message |
 
Exit status is always 0 when everything runs fine, and 1 in the event of an error.

Example for running the test files:
>.\LammpsToVTK -i output.dump -o vtks\OutputVTK -m atomDefinition.in -n 8 -l -t
