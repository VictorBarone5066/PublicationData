import pymatgen as mg
import pymatgen.symmetry.analyzer as sa
import headerBoltzmann as hb

POSCAR_LOC = "agbii4-cubic-cryst-111-gga-poscars//"
CONTCAR_LOC = "agbii4-cubic-cryst-111-gga-contcars//"
BOLTZMANN_LOC = "agbii4-cubic-cryst-444-gga-boltzmann//output.csv"
OUT_LOC = "symms.csv"
TEMPERATURE = 600 #kelvin

outputDat = hb.GetOutputData(BOLTZMANN_LOC)
hb.SetAllProbs(datLis=outputDat, temp=TEMPERATURE)

with open(OUT_LOC, 'w') as outfile:
    outfile.write("genNum,pSpaceGroupNum,pSystem,pSpaceGroupSym,cSpaceGroupNum,cSystem,cSpaceGroupSym,\n")
    #POSCARS
    for x in range(0, 12870):
        if(x%100 == 0):
            print(x)

        s = sa.Structure.from_file(POSCAR_LOC + 'POSCAR0-' + str(x))
        a = sa.SpacegroupAnalyzer(s, symprec=0.1, angle_tolerance=0.1)
        n = a.get_space_group_number()
        k = a.get_crystal_system()
        r = a.get_space_group_symbol()
    #CONTCARS
        s_ = sa.Structure.from_file(CONTCAR_LOC + 'CONTCAR-' + str(x))
        a_ = sa.SpacegroupAnalyzer(s_, symprec=0.01, angle_tolerance=0.1)
        n_ = a_.get_space_group_number()
        k_ = a_.get_crystal_system()
        r_ = a_.get_space_group_symbol()

        outfile.write(str(x) + ',' + str(n) + ',' + str(k) + ',' + str(r) + \
                      ',' + str(n_) + ',' + str(k_) + ',' + str(r_) + "\n")

    outfile.close()
