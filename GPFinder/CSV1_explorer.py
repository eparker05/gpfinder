from os import listdir, getcwd
from os.path import isfile, join
from collections import namedtuple
from re import match
import sys

matchTuple = namedtuple("matchTuple", ['featureId', 'matchId',
                                        'rt', 'mz', 'seq',
                                        'xyz', 'abc', 'other'])

thisPath = getcwd()
onlyCSVs = [f for f in listdir(thisPath) if isfile(join(thisPath,f)) \
                and f.split(".")[-1] == "csv" and "CSV_1" in f]
                
for i in range(len(onlyCSVs)):
    print("{}: '{}'".format(i+1, onlyCSVs[i][0:70]))
fIndex = -1

while fIndex < 0 or fIndex >= len(onlyCSVs):
    try:
        fIndex = int(raw_input("Select a CSV: ")) - 1
        if fIndex < 0 or fIndex >= len(onlyCSVs):
            raise
    except:
        print("enter a valid number")

with open(onlyCSVs[fIndex], 'r') as handle:
    
    lineTupleList = []
    
    titlePassed = False
    rt, mz = 0,0
    for line in handle:
        lineList = line.split(",")
        
        if titlePassed and len(lineList)>35:
            #import all lines
            featNo = int(lineList[0])
            matchNo = int(lineList[6])
            
            try:
                rt2 = float(lineList[4].strip())
                if rt2 < 0.9:
                    rt2 = rt2*60
                mz2 = float(lineList[1].strip())
                rt = rt2
                mz = mz2
            except:
                pass
            seq = "".join(lineList[20].split("-"))
            
            glycan = "_".join([lineList[24],
                              lineList[25],
                              lineList[26],
                              lineList[27],
                              lineList[29]])
                              
            xyz =  [result.strip() for result in lineList[33].split(";") if match('''.\d+''', result.strip())]
            xyz.sort(key =lambda ionString: int(match('''.\d+''', ionString).group()[1:]))
            abc = [result.strip() for result in lineList[34].split(";") if match('''.\d+''', result.strip())]
            abc.sort(key =lambda ionString: int(match('''.\d+''', ionString).group()[1:]))
            other = [result.strip() for result in lineList[35].split(";")]
            
            newMatch = matchTuple(featureId=featNo,
                                  matchId=matchNo,
                                  rt=rt, mz=mz, 
                                  seq=seq,
                                  xyz=xyz,
                                  abc=abc,
                                  other=other)
            lineTupleList.append(newMatch)
            
        #find title
        
        if lineList[0:2] == ["Feature No.","mz"]:
            titlePassed = True
            
while True:
    try:
        getStr = raw_input("Enter a feature or a feature and a match number (space delimited): ")
        getList = [int(idx.strip()) for idx in getStr.split()]
        printed = False
        if len(getList) == 1:
            for T in lineTupleList:
                if [T.featureId] == getList:
                    print '\n\n'
                    print "Feature mz:{}  at rt:{}".format(T.mz, T.rt)
                    print "Sequence: {}\n".format(T.seq)
                    print "xyz type ions:"
                    for ion in T.xyz:
                        print ion
                    print "\n"
                    print "abc type ions:"
                    for ion in T.abc:
                        print ion
                    print "\n"
                    print "other\oxonium type ions:"
                    for ion in T.other:
                        print ion
                    print "\n"
                    printed = True
                    break
        elif len(getList) == 2:
            for T in lineTupleList:
                if [T.featureId, T.matchId] == getList:
                    print '\n\n'
                    print "Feature mz:{}  at rt:{}".format(T.mz, T.rt)
                    print "Sequence: {}\n".format(T.seq)
                    print "xyz type ions:"
                    for ion in T.xyz:
                        print ion
                    print "\n"
                    print "abc type ions:"
                    for ion in T.abc:
                        print ion
                    print "\n"
                    print "other\oxonium type ions:"
                    for ion in T.other:
                        print ion
                    print "\n"
                    printed = True
                    break
        if not printed:
            print "feature requested not available"
          
    except:
        print("enter a valid number or pair of numbers")
        sys.exit()
            