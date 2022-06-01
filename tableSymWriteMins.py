import pymatgen as mg
import pymatgen.symmetry.analyzer as sa
import headerBoltzmann as hb

POSCAR_LOC = "agbii4-cubic-cryst-111-gga-poscars//"
CONTCAR_LOC = "agbii4-cubic-cryst-444-gga-contcars//"
BOLTZMANN_LOC = "agbii4-cubic-cryst-444-gga-boltzmann//output.csv"
OUT_LOC = "symms.csv"
TEMPERATURE = 600 #kelvin

outputDat = hb.GetOutputData(BOLTZMANN_LOC)
hb.SetAllProbs(datLis=outputDat, temp=TEMPERATURE)

with open(OUT_LOC, 'w') as outfile:
    outfile.write("genNum,nrg,pSpaceGroupNum,pSystem,pSpaceGroupSym,cSpaceGroupNum,cSystem,cSpaceGroupSym,\n")
    #POSCARS
    for x in outputDat:
        s = sa.Structure.from_file(POSCAR_LOC + 'POSCAR0-' + str(x.wDirectory))
        a = sa.SpacegroupAnalyzer(s, symprec=0.1, angle_tolerance=0.1)
        n = a.get_space_group_number()
        k = a.get_crystal_system()
        r = a.get_space_group_symbol()
    #CONTCARS
        s_ = sa.Structure.from_file(CONTCAR_LOC + 'CONTCAR-' + str(x.wDirectory))
        a_ = sa.SpacegroupAnalyzer(s_, symprec=0.01, angle_tolerance=0.1)
        n_ = a_.get_space_group_number()
        k_ = a_.get_crystal_system()
        r_ = a_.get_space_group_symbol()

        outfile.write(str(x.wDirectory) + ',' + str(x.nrg) + ',' + str(n) + ',' + str(k) + ',' + str(r) + \
                      ',' + str(n_) + ',' + str(k_) + ',' + str(r_) + "\n")

    outfile.close()
