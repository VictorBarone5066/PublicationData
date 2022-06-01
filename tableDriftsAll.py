import headerPoscar as h

POSCARS_LOC = "agbii4-cubic-cryst-111-gga-poscars//"
CONTCARS_LOC = "agbii4-cubic-cryst-111-gga-contcars//"
ORIG_FILE_NAME = "POSCAR0-"
FIN_FILE_NAME = "CONTCAR-"

OUTPUT_FILE_LOC = "driftsAll.csv"

DIFF = 0.75 #fractional deviation to see if an atom went to an equiv mirror image
            #if |final position - init position| > DIFF, need to reset atom location

def Dist(st, fi):
    return ((fi.a-st.a)**2. + (fi.b-st.b)**(2.) + (fi.c-st.c)**(2.))**(1./2.)

avgDevAg, avgDevBi, avgDevI = 0.0, 0.0, 0.0
counterAg, counterBi, counterI = 0, 0, 0
maxDevAg, maxDevBi, maxDevI = 0.0, 0.0, 0.0 ##angstroms
#This is going to take a while
for i in range(0, 12870):
    pos = h.Poscar(POSCARS_LOC + ORIG_FILE_NAME + str(i))
    con = h.Poscar(CONTCARS_LOC + FIN_FILE_NAME + str(i))

    pos.ConvertToDirect()
    con.ConvertToDirect()

    #Make sure atoms didn't take a trip from one side of the cell to the other
    for atomIndex in range(0, len(pos.atoms)):
        while(abs(con.atoms[atomIndex].a - pos.atoms[atomIndex].a) > DIFF):
            if(con.atoms[atomIndex].a - pos.atoms[atomIndex].a > DIFF): ###contcar location > poscar location
                con.atoms[atomIndex].a = con.atoms[atomIndex].a - 1.0
            elif(con.atoms[atomIndex].a - pos.atoms[atomIndex].a < DIFF): ###contcar location < poscar location
                con.atoms[atomIndex].a = con.atoms[atomIndex].a + 1.0
        while(abs(con.atoms[atomIndex].b - pos.atoms[atomIndex].b) > DIFF):
            if(con.atoms[atomIndex].b - pos.atoms[atomIndex].b > DIFF): ###contcar location > poscar location
                con.atoms[atomIndex].b = con.atoms[atomIndex].b - 1.0
            elif(con.atoms[atomIndex].b - pos.atoms[atomIndex].b < DIFF): ###contcar location < poscar location
                con.atoms[atomIndex].b = con.atoms[atomIndex].b + 1.0
        while(abs(con.atoms[atomIndex].c - pos.atoms[atomIndex].c) > DIFF):
            if(con.atoms[atomIndex].c - pos.atoms[atomIndex].c > DIFF): ###contcar location > poscar location
                con.atoms[atomIndex].c = con.atoms[atomIndex].c - 1.0
            elif(con.atoms[atomIndex].c - pos.atoms[atomIndex].c < DIFF): ###contcar location < poscar location
                con.atoms[atomIndex].c = con.atoms[atomIndex].c + 1.0

    pos.ConvertToCartesian()
    con.ConvertToCartesian()


    for atomIndex in range(0, len(pos.atoms)):
        ##Ag
        if(pos.atoms[atomIndex].atomType == "Ag"):
            dev = Dist(pos.atoms[atomIndex], con.atoms[atomIndex])
            counterAg += 1
            avgDevAg += dev
            if(dev > maxDevAg):
                maxDevAg = dev
        ##Bi
        if(pos.atoms[atomIndex].atomType == "Bi"):
            dev = Dist(pos.atoms[atomIndex], con.atoms[atomIndex])
            counterBi += 1
            avgDevBi += dev
            if(dev > maxDevBi):
                maxDevBi = dev
        ##I
        if(pos.atoms[atomIndex].atomType == "I"):
            dev = Dist(pos.atoms[atomIndex], con.atoms[atomIndex])
            counterI += 1
            avgDevI += dev
            if(dev > maxDevI):
                maxDevI = dev

    if(i%100 == 0):
        print(i)

with open(OUTPUT_FILE_LOC, 'w') as outfile:
    outfile.write("elem,avgDev(angst),maxDev(angst)\n")
    outfile.write("Ag," + str(avgDevAg/float(counterAg)) + ',' + str(maxDevAg) + "\n")
    outfile.write("Bi," + str(avgDevBi/float(counterBi)) + ',' + str(maxDevBi) + "\n")
    outfile.write("I," + str(avgDevI/float(counterI)) + ',' + str(maxDevI) + "\n")

    outfile.close()
