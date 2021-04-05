#!/usr/bin/python
__author__ = "A Ram"
__copyright__ = "Copyright 2021"
__license__ = "GPL"
__version__ = "1.0.0"

### Usage: python pset1.py -i <input file> -s <score file>
### Example: python pset1.py -i input.txt -s blosum62.txt
### Note: Smith-Waterman Algorithm

import argparse
import numpy as np
import pandas as pd

### Read in the arguments from command line

parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False,default=-2)
parser.add_argument('-e', '--extgap', help='extension gap', required=False,default=-1)

args = parser.parse_args()

### Implement your Smith-Waterman Algorithm

def runSW(inputFile, scoreFile, openGap, extGap):

	#convert gap opening/extension penalties to ints since argparse passes them in as strings
	gapOpen = int(openGap)
	gapExtend = int(extGap)

	#####################################
	#grab input sequences from inputFile
	#####################################

	with open(inputFile, 'r') as infile: 
		seqs = infile.read().splitlines()
		seq1 = seqs[0]
		seq2 = seqs[1]

	infile.close()

	##########################################
	#grab ordering of scores from blosum file
	##########################################

	'''
	The current blosum file has weird spacing that messes things up for me when I want to turn it into a df
	What I'm doing is just taking that file, stripping all the weird spacing and making a new one with uniform spacing that's easier to work with
	I'm creating a new blosum file to work with here but I still use the old one as input, so this should satisfy the assignment requirements
	From there I can convert nicely into a df
	'''

	blosumold = open(scoreFile, 'r')
	blosumnew = open("blosum_new.txt", "w")

	for line in blosumold: 
		#adapted from https://pythonexamples.org/python-replace-multiple-spaces-with-single-space-in-text-file/
		#getting rid of the inconsistent spacing in the blosum file-now each char is single space-separated
		blosumnew.write(' '.join(line.split())+"\n")

	blosumnew.close()
	blosumold.close()	

	#we can use pd.read_csv to convert the text file into a df
	#see here https://stackoverflow.com/questions/21546739/load-data-from-txt-with-pandas
	blosumdf = pd.read_csv("blosum_new.txt", delimiter = " ")
	
	###########################################
	#initialize scoring and traceback matrices
	###########################################

	#create arrays for the matrix indices so we can index into the blosum df -- they're just the letters of the input sequences
	colNames = [char for char in seq1] 
	rowNames = [char for char in seq2] 

	#add an empty string to the beginning of both arrays because it's supposed to be a m+1 x n+1 matrix 	
	colNames.insert(0, "")
	rowNames.insert(0,"")

	#create m+1xn+1 matrix for scoring
	scoringmat = np.zeros((len(rowNames),len(colNames)))	

	#create a traceback 2d list [I tried using a numpy matrix but it doesn't like to have lists as values]
	#populate it with the coords of each cell, we can add to it as we do the traceback
	tracebackmat = [[[[i,j]] for j in range(len(colNames))] for i in range(len(rowNames))]

	#variables for best score/location of this best score--we'll use this for the traceback
	bestscore = 0
	bestloc = [0,0]

	###################
	#calculate scores
	###################

	#to prevent an index out of bounds error we'll start our counting at 1
	for i in range(1,len(rowNames)):
		for j in range(1,len(colNames)):
	
			#get row and column names so we can look up blosum scores
			row = rowNames[i]
			col = colNames[j] 		

			#score of upper left neighbor + score of aligning the two residues of interest
			#if i don't use the .values[0] it'll cause issues
			#i'm indexing into the df based on the current two residues
			#values idea from here https://stackoverflow.com/questions/49277640/how-to-extract-values-from-a-pandas-dataframe-rather-than-a-series-without-ref
			#i have to write [0][0] because values[0] gives me something like "[3]" and I need to work with an int
			diag = scoringmat[i-1][j-1]+ int(blosumdf.loc[[row],[col]].values[0][0])

			#score of neighbor to the left [vertical gap] #previously set to horiz gap
			#takes the form max_k(uk+v) where u=extension penalty, k=gap size, v=opening penalty			
			#we calculate all possible values of k and find the value that produces the largest score
			vertGapVals = np.array([gapExtend*(i-k-1)+gapOpen for k in range(i)])
			vertGapArray = scoringmat[0:i,j]+vertGapVals
			vertGap = max(vertGapArray)

			#score of neighbor above [horizontal gap], calculated same as above
			horizGapVals = np.array([gapExtend*(j-k-1)+gapOpen for k in range(j)])
			horizGapArray = scoringmat[i,0:j]+horizGapVals
			horizGap = max(horizGapArray)

			#calculate the new score and add to scoringdf
			#we add 0 into the newScore calculation because the matrix can't have negative values, so 0 has to be the smallest
			newScore = int(max(diag, horizGap, vertGap, 0))
			scoringmat[i][j]=int(newScore)

			#if this is the highest score set bestscore and bestloc to its value and coords respectively
			if newScore>bestscore:
				bestscore = newScore
				bestloc = [i,j]

			###################################
			#add traceback to traceback matrix
			###################################

			#includes the path used to reach a given cell
			#each iteration you add the neighbor [diagonal, horizontal, vertical] that produced the highest score
			
			#highest score came from diagonal neighbor
			if newScore == diag: 
				tracebackmat[i][j].append([i-1,j-1])
			
			#highest score came from horizontal gap
			elif newScore == horizGap: 
				#find all the values that produced the max horizontal gap score
				kvals = np.where(horizGapArray==newScore)[0]
				#there may be a few values that produced the max score so we'll keep them all for now	
				for k in kvals:
					tracebackmat[i][j].append([i,k])
			
			#highest score came from vertical gap--same idea as horiz gap
			elif newScore == vertGap: 
				#find all the values that produced the max horizontal gap score
				kvals = np.where(vertGapArray==newScore)[0]
				#there could be a few of them so keep them all for now
				for k in kvals:
					tracebackmat[i][j].append([k,j])

	print("Best score: " + str(bestscore))
	print(scoringmat)


	###############################################
	#Run the traceback and find the best alignment
	###############################################

	#extract the row/col producing the best score
	bestrow = bestloc[0]
	bestcol = bestloc[1]

	#initialize arrays for the sequences to print
	#we have to do this and then convert to string because strings aren't mutable
	topvals = ["(",")"] #first sequence corresponding to the matrix rows
	matching = [" "," "] #the "|"'s for matching values
	bottomvals = ["(",")"] #second sequence corresponding to the matrix columns

	#vars to hold current col/row, just so i don't overwrite bestrow and bestcol
	row = bestrow
	col = bestcol

	#row and col are currently at the end of the alignment so we can add everything after it to the alignment arrays
	topvals.insert(1,seq2[bestrow-1])
	bottomvals.insert(1,seq1[bestcol-1])

	#run the traceback -- until you hit a cell with score 0 keep going
	while scoringmat[row][col] != 0: 
	
		#the [0] value will be that of the current cell coords so choose the one after it
		#I know this is inefficient, it's just how I wrote the traceback matrix
		[r,c] = tracebackmat[row][col][1]

		#diagonal
		if r == row-1 and c == col-1: 
			
			#extract the cell producing the [r,c] score
			#we're working right->left on the alignment so each new char gets placed immediately after the "("
			topvals.insert(1,seq2[r])
			bottomvals.insert(1,seq1[c])

			#check to see if these values are equivalent, if so make a note in the matching array
			if seq2[r] == seq1[c]: 
				matching.insert(1,"|")	
			else: 
				matching.insert(1," ")
		
		#horizontal gap 
		elif r == row and c < col:
			
			#if i'm at cell [1,5] and the traceback produces cell [1,3] i need to indicate that i have a gap of length 2
			#I tried range(col,c) and (col-1,c-1) worked ... weirdly
			#we have to count backwards because we're starting from the end of the alignment--hence the -1
			for i in range(col-1,c-1,-1): 
				topvals.insert(1,'-')
				bottomvals.insert(1,seq1[i])
				matching.insert(1," ")

		#vertical gap--same idea as horiz gap
		elif r < row and c == col: 
			#again we have to count backwards
			for j in range(row-1,r-1,-1):
				topvals.insert(1,seq2[j])
				bottomvals.insert(1,'-')
				matching.insert(1," ")

		row = r
		col = c

	###############################################################
	#process/format the traceback and produce the strings to print
	###############################################################

	#get rid of the last string in the alignment which has a score of 0
	#the -2 is because the alignment is of the form "[alignment0]" so to pop that 0 and produce "[alignment]" we need the -2 index i.e. second from last
	topvals.pop(-2)
	bottomvals.pop(-2)

	#now, let's add on the rest of the alignment

	#portions of the string before the alignment
	topvalStart = seq2[:row]
	bottomvalStart = seq1[:col]

	#portions of the strings after the alignment
	topvalEnd = seq2[bestrow:]
	bottomvalEnd = seq1[bestcol:]	

	#the strings have different length so we need to figure out how many spaces to append to the shorter one
	#this just takes the length of the two pre-alignment sequences and sees which is bigger
	toAppend = len(max([topvalStart,bottomvalStart],key=len))

	#add toAppend-start zeroes to the beginning of the shorter sequence so the alignments show up properly
	topLength = len(topvalStart)
	bottomLength = len(bottomvalStart)

	topGap = toAppend-topLength
	bottomGap = toAppend-bottomLength

	#actually append the zeros here
	topStart = topGap*" " +topvalStart
	bottomStart = bottomGap*" " + bottomvalStart
	#the matching array has to be the same size as the longest string since we've only made it the same size as the alignment
	matchingStart = toAppend*" "

	#combine this all to create the top, bottom and matching alignment strings 
	topString = topStart+"".join(topvals)+topvalEnd
	matchingString = matchingStart+"".join(matching)
	bottomString = bottomStart+"".join(bottomvals)+bottomvalEnd

	#this is a stupid step i have to take for diff -- add some extra " "'s at the end otherwise it says the alignments are "different" when they're not actually
	#I'm just finding the longest string length and appending " "'s to any strings shorter than it
	#It is ugly, redundant code but it is also an ugly, useless step to have to take 
	maxStringLength = len(max([topString, matchingString,bottomString],key=len))
	if len(topString) < maxStringLength: 
		lengthDiff = maxStringLength - len(topString)
		topString = topString + " "*lengthDiff	
	if len(matchingString) < maxStringLength: 
		lengthDiff = maxStringLength - len(matchingString)
		matchingString = matchingString + " "*lengthDiff	
	if len(bottomString) < maxStringLength: 
		lengthDiff = maxStringLength - len(bottomString)
		bottomString = bottomString + " "*lengthDiff	

	#just so yall can see the output I print it nicely in shell
	print(topString)
	print(matchingString)
	print(bottomString)	

	#############################
	#write output to output file
	#############################

	outputfile = open('output.txt','w')
	outputfile.write("-----------\n")
	outputfile.write("|Sequences|\n")
	outputfile.write("-----------\n")
	outputfile.write("sequence1\n")
	outputfile.write(seq1+"\n")
	outputfile.write("sequence2\n")
	outputfile.write(seq2+"\n")
	outputfile.write("--------------\n")
	outputfile.write("|Score Matrix|\n")
	outputfile.write("--------------\n")

	##########################
	#write the scoring matrix
	##########################
	
	#for some reason these sample output files have an extra "\t" at the beginning so, for the purposes of diff, here you go
	outputfile.write("\t")
	for char in colNames: 
		outputfile.write(char+"\t")
	outputfile.write("\n")
	
	for i in range(len(rowNames)):
		outputfile.write(rowNames[i]+"\t")
		for j in range(len(colNames)):
			outputfile.write(str(int(scoringmat[i][j]))+"\t")
		outputfile.write("\n")		

	outputfile.write("----------------------\n")
	outputfile.write("|Best Local Alignment|\n")
	outputfile.write("----------------------\n")
	outputfile.write("Alignment Score:"+str(bestscore)+"\n")
	outputfile.write("Alignment Results:\n")
	
	#################
	#write traceback
	#################

	outputfile.write(bottomString+"\n")
	outputfile.write(matchingString+"\n")
	outputfile.write(topString+"\n")

	outputfile.close()	

###################################
#Run your Smith-Waterman Algorithm
###################################

runSW(args.input, args.score, args.opengap, args.extgap)

