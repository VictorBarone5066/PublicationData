#!/usr/bin/python3.9
import headerPoscar as Poscar
import os
import sys

genNum = -1
seedPoscarLoc = ""
newPoscarLoc = ""
bin = ""
for i in range(0, len(sys.argv)):
    if(sys.argv[i][0] == '-' and (sys.argv[i][1] == 's' or sys.argv[i][1] == 'S')):
        seedPoscarLoc = sys.argv[i + 1]
    if(sys.argv[i][0] == '-' and (sys.argv[i][1] == 'n' or sys.argv[i][1] == 'N')):
        newPoscarLoc = sys.argv[i + 1]
    if(sys.argv[i][0] == '-' and (sys.argv[i][1] == 'g' or sys.argv[i][1] == 'G')):
        genNum = int(sys.argv[i + 1])
    if (sys.argv[i][0] == '-' and (sys.argv[i][1] == 'b' or sys.argv[i][1] == 'B')):
        bin = str(sys.argv[i + 1]) ##string to avoid losing prefixed zeros

origPoscar = Poscar.Poscar(seedPoscarLoc)

         #0     1
elems = ['Ag', 'Bi']
idsToSwap = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

#Create a list of split digits (and exclude the newline character at the end)
arr = [char for char in bin][:]

if(len(arr) != 16):
    print("Uh oh...")
    exit(1)

arr = [int(i) for i in arr]
#Create a new POSCAR instance
newPoscar = Poscar.deepcopy(origPoscar)
newPoscar.comment = str(genNum) + ' ' + str(bin)
#The index of the split digits array corresponds to the ID of the atom to swap,
#and the value at that location in the array tells what atom to swap it with
for i in range (0, len(arr)):
    newPoscar.atoms[idsToSwap[i] - 1].atomType = elems[arr[i]]
#Refresh the atom counts and the atom locations
newAtoms = []
for i in range(0, len(newPoscar.atoms)):
    if(newPoscar.atoms[i].atomType != "Vac"):
        newAtoms.append(newPoscar.atoms[i])
newPoscar.atoms = newAtoms
newPoscar.Refresh()
#Make sure the elements are in order Ag, Bi, I
newPoscar.ChangeAtomOrder(['Ag', 'Bi', 'I'])
#The sum of the first half of the intigers in the array determine the number of swaps,
#so I can decide where to place the files based off of this
newPoscar.Write(str(newPoscarLoc))
