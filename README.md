# SLME-fortran90
A Fortran 90 implementation of the Spectroscopic Limited Maximum Efficiency (SLME) analysis of solar absorbers.

############### CITATIONS AND FURTHER REFERENCES ######################### 
The SLME method was created by Liping Yu and Alex Zunger. Please cite their work: L. Yu, A. Zunger, Phys. Rev. Lett. 108, 068701 (2012). https://doi.org/10.1103/PhysRevLett.108.068701

For more analysis of the method, see the book chapter by M. Bercx, R. Saniz, B. Partoens, D. Lamoen called "Exceeding the Shockley–Queisser Limit Within the Detailed Balance Framework": Bercx M., Saniz R., Partoens B., Lamoen D. (2018) Exceeding the Shockley–Queisser Limit Within the Detailed Balance Framework. In: Angilella G., Amovilli C. (eds) Many-body Approaches at Different Scales. Springer, Cham https://doi.org/10.1007/978-3-319-72374-7_15

also at arxiv: https://arxiv.org/pdf/1705.07762.pdf

############ INSTRUCTIONS FOR USE ######################################

Code compilation: ifort pce-code.f90 -o pce.x

Code Running: ./pce.x dirgap indgap thickmaxm < abs_coef > pce.out

dirgap -> direct gap value in eV

indgap -> indirect gap in eV

thickmaxm -> maximum crystal thickness in meters

abs_coef -> file with the absorption coefficient data
