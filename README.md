> fus-evap-code

This Python code calculates some physical quantites and all possible residues relaated to a fusion-evaporation nuclear reaction with semi-classical calculations. This code is part of an undergraduate thesis of the author at Nawaz Sharif Universiy of Agriculture.

> Details

The structure of this code is based in a study of fudion-evaporation nuclear reactions, from a semi-classical perspective. Quantites such as the cross section of compound nucleus formation and various evaporation residues after it's formation, as well as their cross sections (proportional to the events number)
, are the main objectives of the calculations. The code splits the compound nucleus formation process and it's subsequent decay into several residue nuclei, which ocurs as a sequential
particle emission. In order to priotizea first several residue nuclei, which occurs as a sequential particle emission. In order to priortize a first promximation theory(muchsimpler to
understandfor new students interested in this topic), different nuclear models, with semi classical and statistical origin, related to projectile-target fusion, light particle evaporation
(n, p, a) and fussion, were implemented here.

Then, given a target and projectile nucleus (at a certain energy), this code is capable of is calculate cross section of each exit channel of this type of reaction. The theoretical framework used 
focuses in the Bass model and Wong formula (classical formula is also available) to calculate fusion cross section, Monte Carlo algorithm to select a random decay mode. Weisskopf-Ewing formalism
to calculate evaporation decay widths and Bohr-Wheeler model to calculate fission decay width.

> content

This respository conatins:

1. A copy of the thesis on which is based this code(in spanish), in the file tHESIS (in spanish).pdf.
2. The list of all routines used in the main code, in the file routines_fus_evap_code.py.
3. The main code main_fus_evap_code.py which makes the calculations and show all results in console.
4. A list of the information of all isotopesfound in NuDat Database (extracted Dec-09-2024) in isotope_data_NuDat.txt
