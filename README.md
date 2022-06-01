# PublicationData
Data + programs needed for reproduction of published results

All (well, most) of the data needed to reproduce all results in 'Barone et. al: Properties of AgBiI4 using 
high through-put DFT and machine learning methods'.  Missing data includes the 12870 POS / CONT CAR files, 
the eigenvalue files that are input into the average effective mass code, and the Brillouin zone
interpolations output from the effective mass code.  

Most importantly, the random forest model code is 'graphRfCompare.py'  Other than the usual machine learning
packages, you'll need the files 'headerRndmFrst.py' and 'headerPoscar.py'.  It also relies on the data in
the directory "agbii4-cubic-cryst-111-gga-rfCompare//' to run.  
