#!/usr/bin/python

__author__ = "Yekaterina Kovalyova"
__copyright__ = "Copyright 2017"
__credits__ = ["Yekaterina Kovalyova"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Yekaterina Kovalyova"
__email__ = "yekaterina.kovalyova@yale.edu"

### Usage:      python crisprTargets.py -i <input DNA file> -g <guide RNAs file>
### Example:    python crisprTargets.py -i DNAseq.txt -g gRNA.txt
### Note:       Tool for finding CRSIPR guide RNA targets and off-targets.

import argparse
import numpy as np

parser = argparse.ArgumentParser(description='gRNA targets and off-targets')
parser.add_argument('-i', '--input', help='input DNA file', required=True)
parser.add_argument('-g', '--grna', help='gRNA file', required=True)
args = parser.parse_args()

# Main method
# Defines alignment parameters and PAM sequence(s)
# Calls on other fxns to get alignments and fund targets
def gRNATargets(DNAFile, gRNAFile):
    #################################
    ###Define Alignment Parameters###
    #################################
    match = 1 #score for matched pair in alignment
    mismatch = 0 #score for mismatched pair in alignment
    initGapPen = 1 #initial gap penalty
    extGapPen = 10 #extension gap penalty; some arbitrary large number
    maxMis = 5 #maximum mismatches/gaps allowed for gRNA to induce effect
    params = [match, mismatch, initGapPen, extGapPen, maxMis] #compile the params together
    
    #Define PAM sequence(s) to search for
    #(may be changed, or have multiple, depending on type of Cas9 used)
    PAM = ["N", "G", "G"]  
    
    ################################################################
    ###Get DNA, gRNA, and write the targets/off-targets to output###
    ################################################################
    output = open('gRNAtargets.txt', 'w') #write results to gRNAtargets.txt
    #Get input DNA out of file
    seq = open(DNAFile, 'r')
    dnaSeq = []
    for char in seq.readline():
        if char != '\n' and char != '\r':
            dnaSeq.append(char)
    seq.close()
    #Get guide RNAs out of file
    rna = open(gRNAFile, 'r')
    for line in rna:
        gRNA = []
        for char in line:
            #Convert U to T for comparison to DNA
            if char == "U":
                gRNA.append("T")           
            if char != '\n' and char != '\r':
                gRNA.append(char)
        ### For every gRNA, find the targets/off-targets in given DNA
        output.write("---Forward---\n \n")
        getTargets(dnaSeq, gRNA, PAM, params, output)
        ### And check for potential targets/off-targets in reverse complement of DNA
        output.write("---Reverse Complement---\n \n")
        getTargets(getRevComp(dnaSeq), gRNA, PAM, params, output)
    rna.close()
    output.close()

# Produces reverse complement of given (DNA) sequence
# Does not account for presence of uracil
def getRevComp(seq):   
    comp = []
    for i in (len(seq)-1, -1, -1):
        if seq[i] == "A":
            comp.append("T")
        elif seq[i] == "T":
            comp.append("A")
        elif seq[i] == "C":
            comp.append("G")
        elif seq[i] == "G":
            comp.append("C")
    return comp 

# Smith-Waterman algorithm for alignment
# Takes in DNA sequence, gRNA, and alignment parameters defined in main method
def align(DNA, gRNA, maxMis, initGapPen, extGapPen, match, mismatch): 
    arr = np.zeros((len(DNA)+1,len(gRNA)+1)) #array with scores
    gaps = np.zeros((len(DNA)+1,len(gRNA)+1)) #corresponding gap distances
    prev = np.zeros((len(DNA)+1,len(gRNA)+1)) #pointers for traceback    
    
    for i in range(1,arr.shape[0]):
        for j in range(1,arr.shape[1]):
            if DNA[i-1] == gRNA[j-1]:
                s = match
            else:
                s = mismatch
                
            #Consider each possible score for current member of matrix
            diag = arr[i-1,j-1]+s #score if looking at i-1,j-1
            
            if gaps[i-1,j] > 0: #score if looking at i-1,j, consider length of gap
                left = arr[i-1,j]-extGapPen
            else:
                left = arr[i-1,j]-initGapPen
            
            if gaps[i,j-1] > 0: #score if looking at i,j-1, consider length of gap
                up = arr[i,j-1]-extGapPen
            else:
                up = arr[i,j-1]-initGapPen
            
            #Find max value for member, populate arr, gaps, and prev accordingly
            if max(diag, left, up, 0) == diag:
                arr[i,j] = diag
                gaps[i,j] = 0
                prev[i,j] = 1 #value of 1 corresponds to pointer to diagonal member
            elif max(diag, left, up, 0) == left:
                arr[i,j] = left
                gaps[i,j] = gaps[i-1,j] + 1
                prev[i,j] = 2 #value of 2 corresponds to pointer to left member
            elif max(diag, left, up, 0) == up:
                arr[i,j] = up
                gaps[i,j] = gaps[i,j-1] + 1
                prev[i,j] = 3 #value of 3 corresponds to pointer to up member
            else:
                arr[i,j] = 0
                gaps[i,j] = 0
                prev[i,j] = 0 #value of 0 corresponds to no pointer
    
    #Perform traceback only if alignment has at most max number of mismatches
    if np.amax(arr) > (len(gRNA) - maxMis - 1):
        #For traceback, populate the arrays for proper alignment in reverse
        alDNARev = [] #seq A, for alignment, in reverse
        alLine = [] #corresponding array for lines to indicate matching segments, in reverse
        algRNARev = [] #seq B, for alignment, in reverse
 
        i = arr.shape[0] - 1
        j = arr.shape[1] - 1
        #Perform alignment (traceback) based on scores in arr and pointers from prev
        currScore = arr[i,j]        
        while currScore > 0: #perform traceback until reach score 0
            if prev[i,j] == 1: #go back diagonally
                alDNARev.append(DNA[i-1])
                algRNARev.append(gRNA[j-1])
                #add line for matching segments
                if DNA[i-1] == gRNA[j-1]:
                    alLine.append('|')
                else:
                    alLine.append(' ')
                i = i-1
                j = j-1            
            elif prev[i,j] == 2: #go left        
                alDNARev.append(DNA[i-1])
                algRNARev.append('-')
                alLine.append(' ')
                i = i-1
            elif prev[i,j] == 3: #go up
                alDNARev.append('-')
                algRNARev.append(gRNA[j-1])
                alLine.append(' ')
                j = j-1
            currScore = arr[i,j]
        #Populate the rest of DNA and gRNA in reverse
        while i>0 and j>0: #if both DNA and gRNA left to populate
            alDNARev.append(DNA[i-1])
            algRNARev.append(gRNA[j-1])
            if DNA[i-1] == gRNA[j-1]:
                    alLine.append('|')
            else:
                alLine.append(' ')
            i = i-1
            j = j-1
        while i>0: #if only DNA left
            alDNARev.append(DNA[i-1])
            algRNARev.append(' ')
            alLine.append(' ')
            i -= 1
        while j>0: #if only gRNA left
            alDNARev.append(' ')
            algRNARev.append(gRNA[j-1])
            alLine.append(' ')
            j -= 1

        return (alDNARev, algRNARev, alLine, np.amax(arr))
    #if alignment yields too low sequence similarity
    else:
        #don't even bother aligning if too many mismatches
        return (0, 0, 0, 0)

# Determine the gRNA targets and off-targets in the given dnaSeq
# Output file provided to write results to
def getTargets(dnaSeq, gRNA, PAM, params, output):
    #Unpack the alignment parameters from params
    #Params provided to pass to align fxn
    match = params[0]
    mismatch = params[1]
    initGapPen = params[2]
    extGapPen = params[3]
    maxMis = params[4]

    alDNARev = [] #DNA, for alignment, in reverse
    alLine = [] #corresponding array for lines to indicate matching segments, in reverse
    algRNARev = [] #gRNA, for alignment, in reverse    
    
    #Iterate backward
    for i in range(len(dnaSeq)-1, len(gRNA) + 3, -1):
        #Find PAM sequence in DNA
        if dnaSeq[i] == PAM[2] or PAM[2] == "N":            
            if dnaSeq[i-1] == PAM[1] or PAM[1] == "N":               
                if dnaSeq[i-2] == PAM[0] or PAM[0] == "N":
                    #Check alignment between DNA right before PAM sequence, and gRNA
                    (subDNA, subRNA, subLine, score) = align(dnaSeq[i-2-len(gRNA): i-2], gRNA, 
                                                             maxMis, initGapPen, extGapPen, 1, 0)                   
                    if not subDNA == 0: #if alignment is good enough
                        ####################################################################
                        ###Populate and format reverse arrays to properly write to output###
                        ####################################################################
                        curDNARev = list(alDNARev)
                        curLine = list(alLine)
                        curRNARev = list(algRNARev)
                        
                        curDNARev.extend([dnaSeq[i],dnaSeq[i-1], dnaSeq[i-2], ' ) '])
                        curDNARev.extend(subDNA)
                        curDNARev.extend([' ( '])
                        
                        curLine.extend([' ', ' ', ' ', '   '])
                        curLine.extend(subLine)
                        curLine.extend(['   '])
                        
                        curRNARev.extend([' ', ' ', ' ', '   '])
                        curRNARev.extend(subRNA)
                        curRNARev.extend(['   '])
                        
                        for j in range(i-4-len(gRNA), -1, -1):
                            curDNARev.append(dnaSeq[j])
                            curLine.append(' ')
                            curRNARev.append(' ')
                        output.write("Guide RNA: ")
                        for j in range(0, len(gRNA)):
                            output.write(gRNA[j])
                        output.write("   ")
                        if score == len(gRNA):
                            output.write("***Target.*** \n")
                        else:
                            output.write("Off-target. \n")
                        for j in range(len(curDNARev)-1, -1, -1):
                            output.write(curDNARev[j])
                        output.write('\n')
                        for j in range(len(curLine)-1, -1, -1):
                            output.write(curLine[j])
                        output.write('\n')
                        for j in range(len(curRNARev)-1, -1, -1):
                            output.write(curRNARev[j])
                        output.write('\n \n')
        alDNARev.append(dnaSeq[i])
        alLine.append(' ')
        algRNARev.append(' ') 

#Run main fxn to find targets and off-targets of gRNA   
gRNATargets(args.input, args.grna)
