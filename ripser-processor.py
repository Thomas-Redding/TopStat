# processor.py
# a python script to parse ripser data

from os import path
import math

'''
NOTE:
This script will only work for data that is formatted by the 
ripser topological data processing script.  
'''

# Simple implementation of insertion sort
def Insertion(bar, barList):
	barList.insert(0,bar)
	for i in range(0,5):
		if barList[i][1] > barList[i+1][1]:
			temp = barList[i]
			barList[i] = barList[i+1]
			barList[i+1] = temp
	barList.pop(0)
	return barList

# Determines length of a given interval
def ProcessInterval(interval):
	interval = interval.rstrip()
	interval = interval.split(",")
	if(interval[1] == " )"):
		return float("inf")
	else:
		lowBound = float(interval[0][1:])
		highBound = float(interval[1][0:-2])
		return math.fabs(highBound - lowBound)

def main():
	dataFile = raw_input("Please enter name of data file: ")

	data = open(dataFile, 'r')

	# Initialize dimension/ list of significant bars
	dimension = -1
	bars = [[-1,-1], [-1,-1], [-1,-1], [-1,-1], [-1,-1]]
	# Process all intervals in file
	for interval in data:
		if interval[0] == "p":
			if dimension != -1:
				print "Significant bars for dimension " + str(dimension) + ":"
				for bar in bars:
					print bar
				dimension += 1
				bars = [[-1,-1], [-1,-1], [-1,-1], [-1,-1], [-1,-1]]
			else:
				dimension += 1
		else:
			difference = ProcessInterval(interval)
			newBar = [interval.rstrip(), difference]
			bars = Insertion(newBar, bars)

	print "Significant bars for dimension " + str(dimension) + ":"
	for bar in bars:
		print bar