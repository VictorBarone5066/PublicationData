import os
import shutil

TOP_100_LOC = "agbii4BestSpecialOutput.csv"
ALL_FILE_LOC = "cubic-special//all//"
NEW_FILE_LOC = "cubic-struct//all//"

with open(TOP_100_LOC, 'r') as infile:
    for n, line in enumerate(infile):
        if(n == 0):
            continue
        wDir = line.split(',')[1]
        shutil.copytree(ALL_FILE_LOC + wDir, NEW_FILE_LOC + wDir)
        os.remove(NEW_FILE_LOC + wDir + "//POSCAR")
        os.rename(NEW_FILE_LOC + wDir + "//CONTCAR", NEW_FILE_LOC + wDir + "//POSCAR")
    infile.close()
