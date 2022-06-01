#!/usr/bin/python3.6

TOL = 0.01
beginLineNum = 8

vbm, cbm = 100, -100
with open("EIGENVAL", 'r') as infile:
	for nu, li in enumerate(infile):
		if(nu > beginLineNum):
			tl = li.split()
			nrg, occ = float(tl[1]), float(tl[2])
			
			#We read up from the v band, so always update that
			if(abs(1.0 - occ) < TOL):
				vbm = nrg
			if(abs(0.0 - occ) < TOL):
				cbm = nrg
				break
	infile.close()

with open("bandgapOut", 'w') as outfile:
	outfile.write("BANDGAP " + str(cbm - vbm) + "\n")
	outfile.close()
