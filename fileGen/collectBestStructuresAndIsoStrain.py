import os
import shutil
import headerPoscar as hp

TOP_100_LOC = "agbii4BestSpecialOutput.csv"
ALL_FILE_LOC = "cubic-special//all//"
NEW_FILE_LOC = "cubic-struct//all//"

with open(TOP_100_LOC, 'r') as infile:
    for n, line in enumerate(infile):
        if(n == 0):
            continue
        wDir = line.split(',')[1]
        shutil.copytree(ALL_FILE_LOC + wDir, NEW_FILE_LOC + wDir)
        #os.remove(NEW_FILE_LOC + wDir + "//POSCAR")
        os.rename(NEW_FILE_LOC + wDir + "//CONTCAR", NEW_FILE_LOC + wDir + "//POSCAR")

        for eps in [0.94, 0.97, 1.00, 1.03, 1.6]:
            name = str(eps)
            os.makedirs(NEW_FILE_LOC + wDir + "//" + name)

            pos = hp.Poscar(NEW_FILE_LOC + wDir + "//POSCAR")
            pos.superCellVecA[0] = pos.superCellVecA[0] * (eps)
            pos.superCellVecB[1] = pos.superCellVecB[1] * (eps)
            pos.superCellVecC[2] = pos.superCellVecC[2] * (eps)
            pos.Write(NEW_FILE_LOC + wDir + "//" + name + "//POSCAR")

        os.remove(NEW_FILE_LOC + wDir + "//POSCAR0")
        os.remove(NEW_FILE_LOC + wDir + "//POSCAR")

    infile.close()
