from os import listdir, getcwd
from os.path import isfile, join
from collections import namedtuple
from re import match
import re
import matplotlib.pyplot as plt
import sys
import math


csv1_fname = "CSV_1_fetuin_dualM-highZslopes.csv"
csv5_fname = "CSV_5_fetuin_dualM-highZslopes.csv"
lowE_mgf = "fetuin_dualM-highZslopes_lower.mgf"
highE_mgf = "fetuin_dualM-highZslopes_higher.mgf"


def read_csv5():
    goodMatches = []
    firstLinePassed = False
    with open(csv5_fname,'r') as handle:
        for line in handle:
            #deal with passing the first line
            if not firstLinePassed:
                if "Peptide Theoretical Mass" in line:
                    firstLinePassed = True
                continue
            
            #parse the products
            splitLine = line.split(",")
            featureNo = int(splitLine[0])
            mass = float(splitLine[9])
            rt = float(splitLine[14])
            
            goodMatches.append({"feature":featureNo, "mass":mass, "rt":rt})
    goodMatches.sort(key = lambda match: match["rt"])
    return goodMatches

def read_CSV1():
    matchTuple = namedtuple("matchTuple", ['featureId', 'matchId',
                                            'rt', 'mz', 'seq', 'glycan',
                                            'xyz', 'abc', 'other', 'pepmass'])

    
    with open(csv1_fname, 'r') as handle:
        
        lineTupleList = []
        
        titlePassed = False
        rt, mz = 0,0
        for line in handle:
            lineList = line.split(",")
            
            if titlePassed and len(lineList)>35:
                #import all lines
                
                
                newRtFound = False
                try:
                    rt2 = float(lineList[5].strip())
                    if rt2 < 0.9:
                        rt2 = rt2*60
                    mz2 = float(lineList[1].strip())
                    rt = rt2
                    mz = mz2
                    newRtFound =True
                except:
                    pass
                
                if newRtFound:
                    pepmass = float(lineList[9])
                    featNo = int(lineList[0])
                    matchNo = int(lineList[7])
                    seq = "".join(lineList[21].split("-"))
                    
                    glycan = "_".join([lineList[25],
                                      lineList[26],
                                      lineList[27],
                                      lineList[28],
                                      lineList[30]])
                    if lineList[32] == "1":
                        glycan = glycan + "_1"
                                      
                    xyz =  [result.strip() for result in lineList[34].split(";") if match('''.\d+''', result.strip())]
                    xyz.sort(key =lambda ionString: int(match('''.\d+''', ionString).group()[1:]))
                    abc = [result.strip() for result in lineList[35].split(";") if match('''.\d+''', result.strip())]
                    abc.sort(key =lambda ionString: int(match('''.\d+''', ionString).group()[1:]))
                    other = [result.strip() for result in lineList[36].split(";")]
                    
                    newMatch = matchTuple(featureId=featNo,
                                          matchId=matchNo,
                                          rt=rt, mz=mz, 
                                          seq=seq,
                                          glycan=glycan,
                                          xyz=xyz,
                                          abc=abc,
                                          other=other,
                                          pepmass=pepmass)
                    lineTupleList.append(newMatch)
                    
            #find title
            
            if lineList[0:2] == ["Feature No.","mz"]:
                titlePassed = True
    return lineTupleList
    
Feature = namedtuple("Feature", ["featureIndex",
                                 "compoundNeutralMass",
                                 "compoundMzValue",
                                 "compoundCharge",
                                 "compoundRtCenter",
                                 "compoundRtTuple",
                                 "compoundIntensity",
                                 "xyzPeakList"])

def parse_mgf(mgfFileName, noFilters = True, verbose=False):

    # Function Description
    #     Read the mgf file, which is the output from MassHunter
    # Function Input:
    #     The mgf file from MassHunter
    
    #setting charges list
        
    inFile = open(mgfFileName, 'r')
    
    
    featureList = []
    
    featureCounter = 0
    
    #used for protonated peptides
    adduct = 1.007276
    
    for line in inFile:
        if verbose:
            print line
        
        line = line.strip()
        #print(line.split())
        if "BEGIN" in line and "IONS" in line:
            featureCounter += 1
            compoundRtCenter = featureCounter  #default value when RTINSECONDS is absent
            print("BEGIN IONS tag seen. Start of Feature " + str(featureCounter))
            xyzPeakList = []
            compoundCharge = 1 # default charge value when CHARGE is absent

                
        elif "PEPMASS" in line:
            massString = line.split('=')[1]
            compoundMzValue = float(massString.split()[0])    
            compoundNeutralMass = compoundMzValue - adduct #assumes charge = 1, this assumption is fixed later with 'CHARGE' line
            
        elif "CHARGE" in line:
            compoundCharge = int("".join([char
                                          for char in line
                                          if char.isdigit() == True]))
            compoundNeutralMass = compoundMzValue * compoundCharge - adduct * compoundCharge
            
        elif "TITLE" in line:
            pass # old version uses agilent specific times - compoundRtCenter = float(line.split()[-2])
            
        elif "END" in line and "IONS" in line:
            #check xyzPeakList
            #print("END IONS tag seen. End of Feature " + str(featureCounter) + "\n")

            newFeature = Feature(featureCounter,
                                 compoundNeutralMass,
                                 compoundMzValue,
                                 compoundCharge,
                                 compoundRtCenter,
                                 None, # compoundRtTuple unknown
                                 None, # compoundIntensity unknown
                                 xyzPeakList )

            featureList.append(newFeature)
            del compoundNeutralMass
            del compoundMzValue
            del compoundCharge
            del compoundRtCenter
            del xyzPeakList
       
            
        elif "RTINSECONDS" in line:
            junk, rtSeconds = line.split("=")
            compoundRtCenter = float(rtSeconds)/60.0
            
        
        elif line[0:2].isdigit() == True and len(line.split()) == 2:
            fragmentNeutralMass, fragmentIntensity = line.split()
            # In CEF file reading,
            #     we use x for fragment ion mass, y for fragment intensity
            # In Mgf file reading,
            #     Change 1:
            #     we still use the same legacy notation of x and y
            #     except that x is fragment neutral mass
            #     where the fragment neutral mass is dealt with later
            #         withOUT subtracting the adducts
            #     Change 2:
            #     there is no fragment charge information
            #     which means there is nothing like 1, 2, 3, 4 charges
            #         for us to traverse
            x = float(fragmentNeutralMass)
            y = float(fragmentIntensity)
            multiplier = 1.0  #this may need to be changed
            xyzTuple = (x, y, 1, multiplier)
            xyzPeakList.append(xyzTuple)

    #used for testing:
    import time
    print("found %d features that have diagnostic ions out of %d features in file"%(len(featureList), featureCounter))
    
    return featureList
    
def plot_xyz_peak_list_with_annotations(xyzpeakList, xyzotherPeakList, matchTuple, type="backbone"):
    adduct  = 1.007276
    mzlist, intenlist, dummylist, dummylist2 = zip(*xyzpeakList)
    mzlistoth, intenlistoth, dummylist, dummylist2 = zip(*xyzotherPeakList)
    
    maxinten = float(max(intenlist))
    intenlist = [100*inten/maxinten for inten in intenlist]
    for i in range(len(intenlist)):
        if intenlist[i] < 50:
            intenlist[i] = intenlist[i]*2
    intenzeros = [0 for inten in intenlist]
    
    maxintenoth = float(max(intenlistoth))
    intenzerosoth = [0 for inten in intenlistoth]
    intenlistoth = [100*inten/maxintenoth for inten in intenlistoth]
    for i in range(len(intenlistoth)):
        if intenlistoth[i] < 50:
            intenlistoth[i] = intenlistoth[i]*2
 
    
    fig = plt.figure()
    fig.set_size_inches(7,10)
    

    def xyzabc_sort_fx(xyzabc_entry):
        try:
            massstring = xyzabc_entry.split(" ")[-2]
            masss = float(massstring)
        except:
            masss = 2000
        return masss
        
    
    ################################setting up some peptide info
    glycanseq = matchTuple.glycan
    peptideseq = matchTuple.seq 
    peptideseq = peptideseq.replace("^carb'ed", "")
    firstN = list(re.finditer("""N.[ST]""", peptideseq))
    peplen = len(peptideseq)
    peptideList = [aa for aa in peptideseq]
    if len(firstN) > 0:
        peptideList[firstN[0].start()] = "N*"
        peptideseq = "".join(peptideList)
    cleavageStatus = [0 for p in range(peplen)]
    
    ############################### evaluating backbone cleavages
    ax = plt.axes([0.05, 0.07, 0.9, 0.36])
    ax.vlines(mzlist, intenzeros, intenlist)
    ax.set_ylim(0, 130)
    ax.set_yticklabels([])
    ax.set_yticks([])
    
    types = {"abc":matchTuple.abc, "xyz":matchTuple.xyz, "other":matchTuple.other}
    
    matchset = []
    matchset.extend(matchTuple.xyz)
    matchset.extend(matchTuple.abc)
    matchset.sort(key= xyzabc_sort_fx)
    bb_count = 0
    max_x = 0
    
    for xyzlabel in matchset:
        xyzlabelsplit = xyzlabel.split(" ")
        mz = float(xyzlabelsplit[-2])
        label = xyzlabelsplit[0].split("-")[0]
        color = 'r' if "z" in label or "y" in label else "b"

        for peaki in range(len(mzlist)):
            
            if mzlist[peaki]*0.99995 <= mz <= mzlist[peaki]*1.00005 and intenlist[peaki] > 0.75:
                print "comparing {} to mzlist={}".format(mz, mzlist[peaki]) 
                print("{} {} {} {}".format(mzlist[peaki], 105, 0, -103+intenlist[peaki]))
                arrowTop = 5*(3-(bb_count%4)) + 105
                arrowBottom = -1*(arrowTop-2)
                ax.arrow(mzlist[peaki], arrowTop, 0, arrowBottom+intenlist[peaki], \
                          head_length=0.5, length_includes_head=True, fc=color, ec=color, ls='dashed') 
                ax.text(mzlist[peaki], arrowTop, label, horizontalalignment='center')
                max_x = mzlist[peaki] + 100
                bb_count += 1
                if color == 'r': #this is a y/z ion:
                    peakPlace = peplen - int(label[1:]) -1
                    cleavageStatus[peakPlace] += 1
                else: #this is an a  b or c ion
                    peakPlace = int(label[1:]) -1
                    cleavageStatus[peakPlace] += 2
                break
    if bb_count > 5:
        max_x = max(max_x, 1000)
    else:
        max_x = 1200
    fig.suptitle("Compound: mz={} at rt={}.\nGlycan composition:{}".format( \
                  str(matchTuple.mz)[0:6], str(matchTuple.rt)[0:5], glycanseq), fontsize=14)  
    ax.set_title("High energy CID")
    ax.set_xlabel("m/z")
    ax.set_ylabel("Intensity")
    ax.set_xlim(0, max_x)
    #plt.show()
   

    ################## evaluating backbone annotation
    ax3 = plt.axes([0.05, 0.83, 0.9, 0.08])
    ax3.set_ylim(0, 10)
    ax3.set_yticklabels([])
    ax3.set_yticks([])
    ax3.set_xticklabels([])
    ax3.set_xticks([])
    voffset = 5
    for aai in range(len(peptideList)):
        maxLen = 10*aai + 20
        ax3.text(10*aai + 10, voffset, peptideList[aai], fontsize=14, horizontalalignment='center', verticalalignment='center')
    
    for fragi in range(len(cleavageStatus)):
        hoffset = 10*fragi + 15
        if cleavageStatus[fragi] == 3:
            ax3.plot([hoffset, hoffset, hoffset+2], [voffset, voffset+3, voffset+5], 'r-')
            ax3.plot([hoffset, hoffset, hoffset-2], [voffset, voffset-3, voffset-5], 'b-')
        elif cleavageStatus[fragi] == 2:
            ax3.plot([hoffset, hoffset, hoffset-2], [voffset, voffset-3, voffset-5], 'b-')
        if cleavageStatus[fragi] == 1:
            ax3.plot([hoffset, hoffset, hoffset+2], [voffset, voffset+3, voffset+5], 'r-')
    print cleavageStatus
    #raw_input()
    ax3.set_ylim(0, 10)
    ax3.set_xlim(0, maxLen)
    #axtext = matchTuple.seq + " " + matchTuple.glycan
    #axtext = axtext.replace("^carb'ed", "")
    #ax.text(50, 126, axtext)
    
    ################## evaluating glycan losses
    ax2 = plt.axes([0.05, 0.49, 0.9, 0.3])
    ax2.vlines(mzlistoth, intenzerosoth, intenlistoth)
    ax2.set_ylim(0, 130)
    ax2.set_yticklabels([])
    ax2.set_yticks([])
    matchset = []
    matchset.extend(matchTuple.other)
    matchset.sort(key= xyzabc_sort_fx)
    PN_mass = (matchTuple.pepmass+204.08)
    matchset.append(  "peptide-HexNAc-Y (5/1): " + str(PN_mass) + " (" + str(PN_mass) +")" )
    matchset.append(  "peptide-HexNAc-Y (5/1): " + str((PN_mass+1.007)/2.0) + " (" + str(PN_mass) +")") 
    bb_count = 0
    max_x = 0
    min_x = 1000
    
    for xyzlabel in matchset:
        xyzlabelsplit = xyzlabel.split(" ")
        mz = float(xyzlabelsplit[-2])
        label = xyzlabelsplit[0].replace("Peptide", "Pep_")
        label = label.replace("peptide", "Pep_")
        label = label.replace("glycoPep", "GP")
        label = label.replace("minus", "minus_")
        label = label.replace("Glyco", "G")
        label = label.replace("Peptide", "P")
        label = label.replace("HexNAc", "N")
        label = label.replace("DeoxyHex", "F")
        label = label.replace("Hex", "H")
        label = label.replace("Neu5Ac", "S")
        label = label.replace("-Y", "")
        label = label.replace("-", "")
        color = 'r' if "z" in label or "y" in label else "b"
        lenmzlist = len(mzlistoth)
        for peaki in range(len(mzlistoth)):
            
            
            if mzlistoth[peaki]*0.99995 <= mz <= mzlistoth[peaki]*1.00005 and intenlistoth[peaki] > 05 or \
               (mzlistoth[peaki]*0.99995 <= mz+0.3333 <= mzlistoth[peaki]*1.00005 and intenlistoth[peaki] and peaki+1 < lenmzlist and\
                mzlistoth[peaki+1]*0.99995 <= mz+0.66667 <= mzlistoth[peaki+1]*1.00005 and intenlistoth[peaki]):
                print "comparing {} to mzlistoth={}".format(mz, mzlistoth[peaki]) 
                print("{} {} {} {}".format(mzlistoth[peaki], 105, 0, -103+intenlistoth[peaki]))
                arrowTop = 5*(3-(bb_count%4)) + 105
                arrowBottom = -1*(arrowTop-2)
                ax2.arrow(mzlistoth[peaki], arrowTop, 0, arrowBottom+intenlistoth[peaki], \
                          head_length=0.5, length_includes_head=True, fc=color, ec=color, ls='dashed') 
                ax2.text(mzlistoth[peaki], arrowTop, label)
                max_x = max(max_x, mzlistoth[peaki] + 200)
                min_x = min(min_x, mzlistoth[peaki])
                bb_count += 1
                break
    if bb_count < 3:
        min_x = 200
        max_x = 2000
    max_x = max(max_x, 1000)
    ax2.set_title("Low energy CID")
    #axtext = matchTuple.seq + " " + matchTuple.glycan
    #axtext = axtext.replace("^carb'ed", "")
    #ax2.set_xlabel("m/z")
    ax2.set_ylabel("Intensity")
    min_x = min_x - 100
    ax2.set_xlim(min_x, max_x)
    figname = "figs/{}rt__{}mz.png".format(str(matchTuple.rt)[0:5], str(matchTuple.mz)[0:6])
    plt.savefig(figname, dpi=250, facecolor='w', edgecolor='w',
            orientation='portrait', papertype=None, format=None,
            transparent=False, bbox_inches=None, pad_inches=0.1,
            frameon=None)
    #plt.show() 
               
    plt.cla() #clear 
    plt.clf()  
    plt.close()
        
        
if __name__ == "__main__":
    lineTupleList = read_CSV1()
    
    print("reading low energy")
    #lowEList = parse_mgf(lowE_mgf, noFilters = True)
    print("reading high energy")
    highEList = parse_mgf(highE_mgf, noFilters = True)
    
    goodmatches =  read_csv5()
    for m in goodmatches:
        print "csv5 match with mass= {}, rt = {} and #{}".format(m["mass"], m['rt'], m["feature"])
        for csv1feat in lineTupleList:
            #each csv1Feat is define as
            #matchTuple = namedtuple("matchTuple", ['featureId', 'matchId',
            #                                'rt', 'mz', 'seq',
            #                                'xyz', 'abc', 'other'])
            if csv1feat.featureId == m["feature"]:
                print "csv1 match with mass= {}, rt = {} and #{}".format(csv1feat.mz, csv1feat.rt, csv1feat.featureId)
                for f in highEList:
                    rtdiff = math.fabs(m["rt"] - f.compoundRtCenter) 
                    mzdiff = math.fabs(m["mass"] - f.compoundMzValue)
                    if mzdiff < 0.2 and rtdiff < 0.02:
                        m["backbonelist"] = f.xyzPeakList
                        print "mgf match with mass={}, rt = {}".format( f.compoundMzValue, f.compoundRtCenter)
                        print " "
                    
                        
                        #plot_xyz_peak_list_with_annotations(f.xyzPeakList, csv1feat, type="backbone")
                        #breaks out of looing at the list
                        break
                #raw_input()
                break
    
    lowEList = parse_mgf(lowE_mgf, noFilters = True)
    for m in goodmatches:
        print "csv5 match with mass= {}, rt = {} and #{}".format(m["mass"], m['rt'], m["feature"])
        for csv1feat in lineTupleList:
            #each csv1Feat is define as
            #matchTuple = namedtuple("matchTuple", ['featureId', 'matchId',
            #                                'rt', 'mz', 'seq',
            #                                'xyz', 'abc', 'other'])
            if csv1feat.featureId == m["feature"]:
                print "csv1 match with mass= {}, rt = {} and #{}".format(csv1feat.mz, csv1feat.rt, csv1feat.featureId)
                for f in lowEList:
                    rtdiff = math.fabs(m["rt"] - f.compoundRtCenter) 
                    mzdiff = math.fabs(m["mass"] - f.compoundMzValue)
                    if mzdiff < 0.2 and rtdiff <= 0.02:
                        
                        print "mgf match with mass={}, rt = {}".format( f.compoundMzValue, f.compoundRtCenter)
                        print " "
                    
                        
                        plot_xyz_peak_list_with_annotations(m["backbonelist"],f.xyzPeakList, csv1feat)
                        #breaks out of looing at the list
                        break
                #raw_input()
                break    
    raw_input("done")
    

    
