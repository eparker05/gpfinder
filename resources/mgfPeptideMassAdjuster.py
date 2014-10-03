'''/**   save file as "mgfPepMassAdjuster.py"
 * Python MGF file pepmass corrector.
 * Make sure to edit the first three variables
 *
 *  indicate the observed bias with the correct sign. The
 *  the inverse of this indicated error will be applied
 *  to the data
 *
 * This product is the work of Evan Parker @ UC Davis and is distributed
 * in the public domain
**/'''
 
 
inputFileName = "Fet_high_tryptic-filtered.mgf"
outputFile = "Fet_high_tryptic-filtered2.mgf"
observedBiasPPM = -22
 
 
 
 
#setup some stuf
correctionFactor = 1.0/(1.0+(observedBiasPPM*(1e-06)))
file = open(inputFileName, 'r')
outputList = []
 
#main loop
for line in file:
        line = line.strip()
        try:
            splitLine = line.split('=')
            a = float(splitLine[1])
            if splitLine[0] == "PEPMASS":
                print(splitLine[1])
                splitLine[1] = str(float(a)*correctionFactor)
                outputList.append(splitLine[0]+"="+splitLine[1]+"\n")
            else:
                outputList.append(splitLine[0]+"="+splitLine[1]+"\n")
        except:
            outputList.append(line+"\n")
 
#output
with open(outputFile, 'w') as endFile:
    endFile.writelines(outputList)