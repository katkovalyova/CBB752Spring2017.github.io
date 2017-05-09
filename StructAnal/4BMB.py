#!/usr/bin/python

__author__ = "Yekaterina Kovalyova"
__copyright__ = "Copyright 2017"
__credits__ = ["Yekaterina Kovalyova"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Yekaterina Kovalyova"
__email__ = "yekaterina.kovalyova@yale.edu"

### Usage:      python 4BMB.py -p1 <pdb prot 1> -p2 <pdb prot 2> -u1 <U_RLJ 1> -u2 <U_RLJ 2> 
### Example:    python 4BMB.py -p1 4bmb_aligned.pdb -p2 4bme_aligned.pdb -u1 4BMB_F19.txt -u2 4BMB_F19Y.txt
### Note:       Structural comparison of protein to its mutant.

import argparse
import numpy as np
import math
import mpl_toolkits
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Protein Structure Analysis')
parser.add_argument('-p1', '--p1', help='protein structure 1', required=True)
parser.add_argument('-p2', '--p2', help='protein structure 2', required=True)
parser.add_argument('-u1', '--u1', help='U_RLJs of protein 1', required=True)
parser.add_argument('-u2', '--u2', help='U_RLJs of protein 2', required=True)
args = parser.parse_args()

# Calculates root-mean-square deviation of positions
# of heavy atoms in backbones of two proteins
def RMSD(pdb1, pdb2):
    (res1, atom1, r1) = processPDB(pdb1) #residue #, atom names, & coordinates from pdb file1
    (res2, atom2, r2) = processPDB(pdb2) #residue #, atom names, & coordinates from pdb file2
    
    n = 0 #total number of atoms compared for RMSD calculation    
    rmsd = 0 #RMSD value
    
    bb1 = [0, 0, 0, 0] #backbone 1 [C, N, O, CA]
    bb2 = [0, 0, 0, 0] #backbone 2 [C, N, O, CA]

    #For each residue, get RMSD of the backbone atoms
    #Not optimal, since some residues may be skipped
    for i in range (res1[0], res1[len(res1)-1]):
        #Get coordinates of each backbone atom
        bb1[0] = getCoord(i, 'C  ', [res1, atom1, r1])
        bb2[0] = getCoord(i, 'C  ', [res2, atom2, r2])
        bb1[1] = getCoord(i, 'N  ', [res1, atom1, r1])
        bb2[1] = getCoord(i, 'N  ', [res2, atom2, r2])
        bb1[2] = getCoord(i, 'O  ', [res1, atom1, r1])
        bb2[2] = getCoord(i, 'O  ', [res2, atom2, r2])
        bb1[3] = getCoord(i, 'CA ', [res1, atom1, r1])
        bb2[3] = getCoord(i, 'CA ', [res2, atom2, r2])
        
        n += 4 #Increment # atoms compared (each loop is 4 backbone atoms)

        #Calculate RMSD for each of current backbone atoms
        for j in range (0, 4):
            if not math.isnan(bb1[0][0]):
                rmsd = rmsd + (np.linalg.norm(bb1[j] - bb2[j], 2, 0))**2
            else: #residue number was skipped
                break
    rmsd = (rmsd / n) ** 0.5
    print "RMSD = " + str(rmsd)

# Extracts residue #s, atoms, and coordinates from pdb file
def processPDB(pdb):
    res = [] #residues
    atom = [] #atoms
    r = [] #x, y, z coordinates
    f = open(pdb, 'r')
    for line in f:
        if line[0:4] == "ATOM": #exclude 'TER' and 'END'
            res.append(int(line[23:26]))
            atom.append(line[13:16])
            r.append(np.array((float(line[31:38]), float(line[39:46]), float(line[47:54]))))
    return (res, atom, r)    

# Returns coordinates of atom at specified residue
def getCoord(resNum, atomType, pdb):
    #Check that residue actually exists
    #If not, return not-a-number for each coordinate
    try:
        pdb[0].index(resNum)
    except:
        return [np.nan, np.nan, np.nan]
    
    #Check that next residue exists
    #If not, next residue is one after that
    #Doesn't check if more than one residue is skipped
    try:        
        if resNum == pdb[0][len(pdb[0])-1]:
            begin = pdb[0].index(resNum)-1
            end = len(pdb[0])
        elif resNum == pdb[0][0]:
            begin = 0
            end = pdb[0].index(resNum+1)                
        else:
            begin = pdb[0].index(resNum)-1
            end = pdb[0].index(resNum+1)
    
    except ValueError:
        if resNum == pdb[0][len(pdb[0])-1]:
            begin = pdb[0].index(resNum)-1
            end = len(pdb[0])
        elif resNum == pdb[0][0]:
            begin = 0
            end = pdb[0].index(resNum+2)                
        else:
            begin = pdb[0].index(resNum)-1
            end = pdb[0].index(resNum+2)
    
    #Return coordinates
    for i in range (begin, end):
        if pdb[1][i] == atomType:
            return pdb[2][i]

# Plot repulsive Lennard Jones potential of residues as fxn of Chi_1 and Chi_2
# Hard-coded for the given files
def U_rlj(u1, u2):
    #Process the given files with chi1, chi2, and energies
    (wt_chi1, wt_chi2, wt_u) = processU(u1) #wildtype residue
    (mut_chi1, mut_chi2, mut_u) = processU(u2) #mutant residue

    #######################################################
    ###For wireframe plots, make 2D arrays of each value###
    #######################################################
    
    WT_chi1 = [] #2D array for WT chi1
    WT_chi2 = [] #2D array for WT chi2
    WT_u = [] #2D array for WT energy
    MUT_chi1 = [] #2D array for mutant chi1
    MUT_chi2 = [] #2D array for mutant chi2
    MUT_u = [] #2D array for mutant energy
    count = 0

    #Actually make the arrays 2D and populate
    #Hard-coded 72x72, for 72 chi1 angles x 72 chi2 angles
    for i in range(0, 72):
        WT_chi1.append([])
        WT_chi2.append([])
        WT_u.append([])
        
        MUT_chi1.append([])
        MUT_chi2.append([])
        MUT_u.append([])
        
        for j in range(0, 72):
            WT_chi1[i].append(wt_chi1[count])
            WT_chi2[i].append(wt_chi2[count])
            WT_u[i].append(wt_u[count])
             
            MUT_chi1[i].append(mut_chi1[count])
            MUT_chi2[i].append(mut_chi2[count])
            MUT_u[i].append(mut_u[count])

            count += 1
    
    ###############################################
    ###Plot the energies as fxn of chi1 and chi2###
    ###############################################
    
    fig = plt.figure()
    plot = fig.add_subplot(111, projection='3d')

    #blue = WT, red = mutant
    plot.plot_wireframe(WT_chi1, WT_chi2, WT_u, rstride=1, cstride=1, color="blue")
    plot.plot_wireframe(MUT_chi1, MUT_chi2, MUT_u, rstride=1, cstride=1, color="red")

    plot.set_xlim3d(155, 185)
    plot.set_ylim3d(40,100)
    plot.set_zlim3d(0, 1)    
    plot.set_xlabel('Chi_1')
    plot.set_ylabel('Chi_2')
    plot.set_zlabel('U_RLJ')    
    
    plt.show()

# Get the chi1, chi2, and energies in separate arrays
def processU(Ufile):
    chi1 = []
    chi2 = []
    U = []
    f = open(Ufile, 'r')
    f.readline()    
    for line in f:
        tmp = line.split(" ")
        chi1.append(int(tmp[0]))
        chi2.append(int(tmp[1]))
        U.append(float(tmp[2]))
    f.close()
    return (chi1, chi2, U)  
           
# Run programs
RMSD(args.p1, args.p2)
U_rlj(args.u1, args.u2)