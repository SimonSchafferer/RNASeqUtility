'''
Created on 21.04.2014

@author: schaffrr
'''
import getopt
import sys
import re

class FastaUtils():

    def __init__(self, pathToFasta):
        self._pathToFasta = pathToFasta
  
    def checkValidity(self):
        try:
            f = open(self._pathToFasta, "r")
            print "opening file: " + self._pathToFasta 
            fastaSign = re.compile('^>')
            if fastaSign.match(f.readline()) != None:
                f.close()
                return 1
            else:
                f.close()
                print "no valid fasta file given, each ID line must start with >"
                return -1
            print "Fasta file valid"
            return 1
        except:
            print "cannot open file!"
            return -1
        
    def eliminateNewLineFsa(self, outFile, removeWhiteSpaceFlag, replacementChar, changeToDNA):
        f = open(self._pathToFasta,"r")
        fw = open(outFile,"w")
        print "Removing White Spaces: " + str(removeWhiteSpaceFlag)
        fastaSign = re.compile('>')
        # writing the first line
        firstLine = f.readline()
        toStrip = re.compile("\s+")
        endSign = re.compile("_$")
        
        if removeWhiteSpaceFlag:
            fw.writelines( toStrip.sub("_",firstLine.rstrip()) + "\n" )
        else:
            fw.writelines( firstLine )
        
        for line in f:
            if fastaSign.match(line) != None:
                if removeWhiteSpaceFlag:
                    currl = toStrip.sub(replacementChar,line)
                    currl = endSign.sub("",currl)
                    fw.writelines( "\n" + currl + "\n" )
                else:
                    fw.writelines("\n" + line)
            else:
		if changeToDNA:
                    currL = line.replace("U","T")
                else:
                    currL = line
                fw.writelines(  toStrip.sub("",currL) )
        f.close()
        fw.close()
        return 1

    def fetchFastaEntriesByRegex(self, outFile, regx):
        f = open(self._pathToFasta,"r")
        fw = open(outFile,"w")
        toSearch = re.compile(regx)
        wnext = False
        
        for line in f:
            if wnext:
                if ">" in line:
                    wnext = False
                else:
                    fw.writelines(line)
            if toSearch.search(line) != None:
                fw.writelines(line)
                wnext = True
        
        f.close()
        fw.close()
        return 1      


def main(argv):
    inputfile = ''
    outputfile = ''
    removeWhiteSpaceFlag = False
    replacementChar = ""#default
    changeToDNA = False
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["help","ifile=", "ofile=", "rmSpaceInID=", "replacementChar=", "RNAtoDNA="])
#         opts, args = getopt.getopt(argv, "hi:o", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print 'PreparencRNAfasta.py -i <inputfile> -o <outputfile> --rmSpaceInID <F> --replacementChar <char> --RNAtoDNA <F>'
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', "--help"):
            print 'PreparencRNAfasta.py -i <inputfile> -o <outputfile> --rmSpaceInID <F> --replacementChar <char> --RNAtoDNA <F>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
            print outputfile
        elif opt == "--rmSpaceInID":
            if arg.lower() in ("true", "t"):
                removeWhiteSpaceFlag = True
                print "Removing Whitespace flag set"
        elif opt in ("--replacementChar"):
            replacementChar = arg
            print "Replacing white spaces with " + replacementChar
        elif opt == "--RNAtoDNA":
            if arg.lower() in ("true", "t"):
                changeToDNA = True
                print "Changing U to T"
    if inputfile == "":
        print "input filename must not be empty!"
        sys.exit()
    if outputfile == "":
        print "no output filename, choosing <input>.mod"
        outputfile = inputfile + ".mod"

    fsaUtil = FastaUtils(inputfile)
    if fsaUtil.checkValidity() == 1:
        fsaUtil.eliminateNewLineFsa(outputfile, removeWhiteSpaceFlag, replacementChar, changeToDNA)
    else:
        sys.exit()
        
if __name__ == "__main__":
    main(sys.argv[1:])

