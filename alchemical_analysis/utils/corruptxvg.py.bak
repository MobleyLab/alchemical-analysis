import os
import numpy as np
import commands as c
from shutil import copy2
###===============================###
# Module written by Caitlin Bannan
# From ManipulateXVG
###==============================###

def findDataStart(lines, delin = ' '):
    """
    Finds the line where the data starts

    input:
        lines = list of strings (probably from a data file)
        delin = optional string, tells how the data would be separated, default is a space (' ')
    
    output:
        i = integer where the data stops being a string (returns -1 if no data was found)
        header = lines that make up the strings at the top of the file
        datalines = lines starting from where the data ends
    """
    
    for i, l in enumerate(lines):
        test = l.split(delin)[0]
        # If it can't be turned into a float is is still a part of the header
        try:
            float(test)
            return i, lines[0:i], lines[i:]
        except:
            continue
    return -1, lines, []

def removeCorruptLines(inputFN, outputFN = 'output.xvg', headerLine = None, fail = 10.0):
    """
    This method takes an .xvg file, removes any corrupted lines of data or data lines that are not sufficiently long. 

    input:
        inputFN - (String) name of input file
        outputFN - (optional string) name of output file, default is output.xvg, you can provide the same name as the input file, but it will be over written. 
        outputFN - (optional, String) name of outputfile, will out output a file with the same name as input if not provided
        headerLine - (optional, int) number of lines at the beginning of the file consisting of strings, if None this will be found, if you know there is no header assign header = 0.
        warn - (optional, float) the percent of data that can be removed before a warning is printed. 
        failing - (optional, float) the percent of data that can be removed and still run the program. If more data is run an exception is raised. 

    output:
        no returns, but creates an output file with all of the non-corrupt data from the input file. 
    """

    # If input file not found raise an exception
    if not os.path.isfile(inputFN):
        fileError = Exception("Input File not found!")
        raise fileError
    
    # Read in file to find header and length of data
    inputFile = open(inputFN,'r')
    lines = inputFile.readlines()
    inputFile.close()

    # Need header lines for numpy.genfromtxt to work
    # If not provided find it:
    if headerLine == None:
        headerLine, header, datalines = findDataStart(lines)
        # Turn list of lines in header into a single string
        header = ''.join(header)

    # If headerLine was provided and it's bigger than zero get header:
    elif headerLine > 0:
        header = lines[0:headerLine]
        # Make header lines into string
        header = ''.join(header)
        
    else: # headerLine == 0:
        header = '' # This is default for numpy.savetxt

    # For writing header purposes we do not want the last character to create new line
    if header[-1] == '\n':
        header = header[:-1]

    # length of original data set:
    dataLength = len(lines) - headerLine

    # Use genfromtxt to get data
    # skip_header skips lines at the beginning of the file
    # invalid_raise means that if a line is found with the wrong number of entries then it will be ignored and a warning will be printed
    data = np.genfromtxt(inputFN, skip_header = headerLine, invalid_raise = False)
    
    
    # Figure out how much data was removed:
    try:
        a = len(data[0])
        difLength = dataLength - len(data)
    except:
        difLength = dataLength - 1    
        data = [data]
    failing = int(fail * dataLength/100.0)
    #print failing

    # Print how many lines were removed or
    # Throw an error if more than 'fail'% of data removed
    if difLength == 0:
        print "No corrupt lines found, no data was removed from %s" %(inputFN)
    elif difLength >= failing:
        DataRemovalError = Exception("%i corrupted lines of %i found, that is more than %.1f percent. Removal failed, examine input file  %s" % (difLength, dataLength, fail, inputFN))
        raise DataRemovalError
    else: # Lines removed, but less than warning limit
        if not os.path.exists('xvg-bak'):
            os.makedirs('xvg-bak')
        copy2(inputFN,'./xvg-bak/'+inputFN)
        print "%i corrupt lines of %i found and removed from %s" % (difLength, dataLength, inputFN) 


    # Use numpy savetxt to write data (with header) to the output file
    # fmt assigns the format of the float, for now it is the longest decimal in examples considered
    # header prints the header from the original file to the top of the output
    # comments adds that character to the beginning of ever line in the header
    np.savetxt(outputFN, data, fmt='%.10f', header = header, comments = '')

    #print "file %s created with corrected data" % outputFN 
