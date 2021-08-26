# LammpsToVTK
LammpsToVTK converts a LAMMPS dump file to vtk format. It can take in a 2D or 3D file with an arbitrary amount of other dumped atom properties. 
The read atom properties are assumed to follow the standard LAMMPS format where %x, %y, %z is a vector and all other properties are scalars. 
Reading is done single threaded, but multi-threaded writing is supported. 

Lines and triangles can optionally be created by passing an input file where atoms, bonds, and angles are defined. 
Bonds are used to create lines and angles are used to create triangles. 

Usage is: LammpsToVTK [OPTION] -i [INPUTFILE] -o [OUTPUTFILE]

Options include:

  -h,    Displays the help message.
  
  -i,    Specifies the next argument is the input file.
  
  -l,    Specifies that lines should be generated from the molecule file.
  
  -m,    Specifies the molecule filename from this directory. The molecule
         file must have bond or angle definitions in order for the mapping
         to work correctly. It is best to use the input file specified to
         LAMMPS for the specific output file. Bonds are treated as lines 
         and angles are treated as triangle polygons.
         
  -o,    Specifies that the next argument is the outputfile string.
  
  -t,    Enables triangle generation from the molecule file.
  
  -v,    Enables verbose messaging to the console.
  
  -n,    Specifies the maximum number of threads used, defaults to single-threaded
 
Exit status is always 0 when everything runs fine,  and 1 in the event of
an error.

Example for running the test files:
LammpsToVTK -i output.dump -o vtks\OutputVTK -m atomDefinition.in -n 8 -l -t
