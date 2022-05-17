# SLME-fortran90
SLME implementation in Fortran90 

Code compilation: ifort pce-code.f90 -o pce.x
Code Running: ./pce.x dirgap indgap thickmaxm < abs_coef > pce.out

dirgap -> direct gap value in eV
indgap -> indirect gap in eV
thickmaxm -> maximum crystal thickness in meters
abs_coef -> file with the absorption coefficient data
