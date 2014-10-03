from xml.dom import minidom, Node
from operator import attrgetter
from os.path import splitext
import os
import re
import glob
import time
import numpy
import math
import array
import base64
import struct
import bisect
from copy import deepcopy

from gaussian_kde import *
from pyFindMassMatch import *
from MassScoreObject import *



class GetFeatures():
    
    def __init__(self, config):
    
        self.config = config
        
        #self.configFileParse   will set the variables manually set below
        #self.adduct = 1.0078250321
        self.adduct = 1.00727646
        self.featureList = []
        
        fileName = self.config.get("msFile")
        if self.config["msFormat"] == "cef" and "cef" in fileName:
            self.parse_cef_xml()
        elif self.config["msFormat"] == "mgf" and "mgf" in fileName:
            self.parse_mgf()
        elif self.config["msFormat"] == "mzdata":
            self.parse_mzData_file()
        else:
            raise EnvironmentError("unsupported output, check file name and file extension") 
            
        #initialize parameters used to judge rescore effacacy
        self.passingFeaturesBeforeRescore = 0
        self.passingFeaturesAfterRescore = 0
        self.scoreCutoffForFdrAnalysisFirstPass = 0
        self.scoreCutoffForFdrAnalysis = 0
        #populate feature list


        
    def find_mass_matches(self):

        # Function Description:
        #     For every feature (every parent mass)
        #         Search for fragment mass matches
        # Function Input:
        #     The user-specified configuration file
        #     The mass matches
        # Function Output:
        #     The fragment mass matches for every parent mass match

        listOfMatchFinderObjects = []
        listOfProteinFiles = self.config["listOfProteinFiles"]

        # If there is only one protein file,
        # We do the original way
        if len(listOfProteinFiles) == 1:
            matchFinderObject = findMassMatch(self.config)
            listOfMatchFinderObjects.append(matchFinderObject)
        # If there are more than one proteins
        # For keeping the interface of findMassMatch intact,
        # We replace the config file with the same one
        #     that has the protein file entry replaced by every single protein file
        # The new fake config is named fakeConfigForSingleProtein
        elif len(listOfProteinFiles) > 1:
            for proteinFile in listOfProteinFiles:
                fakeConfigForSingleProtein = self.config
                fakeConfigForSingleProtein["proteinFile"] = proteinFile
                matchFinderObject = findMassMatch(fakeConfigForSingleProtein)
                listOfMatchFinderObjects.append(matchFinderObject)
        
        # iterate through features for every protein and add matching attributes
        # then doing fragment neutral mass matching
        
        for feature in self.featureList:

            for matchFinder in listOfMatchFinderObjects:
                
                matchFinder.find_all_mass_matches(feature)
                    
                for match in feature.massMatchList:
                    match.featureIndex = feature.featureIndex
                    massScoreObject = MassScoreObject(self.config, match)
                        
                    [match.listOfFragmentNeutralMassMatches,
                     match.totalScore,
                     match.numberOfPeptideIonFragmentMatches,
                     match.numberOfGlycanIonFragmentMatches,
                     match.numberOfDeoxyHexIonFragmentMatches,
                     match.numberOfNeu5AcIonFragmentMatches,
                     match.numberOfNeu5GcIonFragmentMatches,
                     match.numberOfCommonIonFragmentMatches,
                     match.numberOfBackboneFragmentMatches
                     ] = massScoreObject.score_feature_object(feature)

                massMatchList = feature.massMatchList
                massMatchList.sort(key = attrgetter("totalScore"), reverse = True)
                feature.massMatchList = massMatchList

                if self.config["printMatchDetails"] == True:
                    for match in feature.massMatchList:
                        featureIndex = feature.featureIndex
                        match.printMatchDetailsWithFeatureIndices()
                            
        # Improved Series Scoring Algorithm:
        # Adding series scores for all mass match objects
        # The series score of every mass match equals to
        #     the total number of features (including the feature where the mass match lies in)
        #     that has at least one mass match (where the mass match can match to itself)
        #         with the same glycosylation site index and same glycan composition at the site
        
        '''FDR_Scale = [0, 0.67, 1.33, 2, 3, 4, 5, 6.5, 8, 10, 12]
        FDR_Result = []
        
        for FDRScale in FDR_Scale:
            for featureIndex in range(len(self.featureList)):
                for massMatchIndex in range(len(self.featureList[featureIndex].massMatchList)):
                    numberOfFeaturesThatContainAtLeastOneMatchOfSameSiteAndGlycan = len([
                            feature
                            for feature in self.featureList
                            if len([massMatch
                                    for massMatch in feature.massMatchList
                                    if (massMatch.glycosylationSiteIndex
                                        == self.featureList[featureIndex].massMatchList[
                                            massMatchIndex].glycosylationSiteIndex)
                                    if (massMatch.glycanCompositionMatched
                                        == self.featureList[featureIndex].massMatchList[
                                            massMatchIndex].glycanCompositionMatched)
                                    ]) > 0
                            ])

                    if numberOfFeaturesThatContainAtLeastOneMatchOfSameSiteAndGlycan == 0:
                        raise Exception("Bug here, should at least be 1")
                    elif numberOfFeaturesThatContainAtLeastOneMatchOfSameSiteAndGlycan > 0:
                        seriesScore = numberOfFeaturesThatContainAtLeastOneMatchOfSameSiteAndGlycan
                        
                    seriesScore = min(math.log(seriesScore), 3)*(FDRScale/3.0)
                    self.featureList[featureIndex].massMatchList[
                            massMatchIndex].seriesScore = seriesScore
            
            # calculate FDR and append to FDR_Result 
            totalMatchCount = 0
            decoyCount = 0
            for feature in self.featureList:
                if feature.massMatchList != []:
                    massMatchWithHighestScoreWithinOneFeature = max(feature.massMatchList,
                                                                    key = lambda feature: feature.totalScore + feature.seriesScore)
                    totalMatchCount += 1
                    if massMatchWithHighestScoreWithinOneFeature.isDecoyUsed == True:
                        decoyCount += 1
            FDR_Result.append(float(decoyCount) / max(1.0, float(totalMatchCount - decoyCount)))
        
        # plotting (for Michael and Evan)
        # do NOT ship code with the plotting activated'''
        '''import matplotlib.pyplot as plt
        plt.plot(FDR_Scale, FDR_Result, 'r-')
        plt.plot(FDR_Scale, FDR_Result, 'o')
        plt.show()'''
        
        #choose best scaling factor to minimize FDR
        '''bestFDRScale = FDR_Scale[FDR_Result.index(min(FDR_Result))]
        
        for featureIndex in range(len(self.featureList)):
            for massMatchIndex in range(len(self.featureList[featureIndex].massMatchList)):
                numberOfFeaturesThatContainAtLeastOneMatchOfSameSiteAndGlycan = len([
                        feature
                        for feature in self.featureList
                        if len([massMatch
                                for massMatch in feature.massMatchList
                                if (massMatch.glycosylationSiteIndex
                                        == self.featureList[featureIndex].massMatchList[
                                                massMatchIndex].glycosylationSiteIndex)
                                        if (massMatch.glycanCompositionMatched
                                                == self.featureList[featureIndex].massMatchList[
                                                        massMatchIndex].glycanCompositionMatched)
                                ]) > 0
                        ])

                if numberOfFeaturesThatContainAtLeastOneMatchOfSameSiteAndGlycan == 0:
                    raise Exception("Bug here, should at least be 1")
                elif numberOfFeaturesThatContainAtLeastOneMatchOfSameSiteAndGlycan > 0:
                    seriesScore = numberOfFeaturesThatContainAtLeastOneMatchOfSameSiteAndGlycan
                                    
                seriesScore = min(math.log(seriesScore), 3)*(bestFDRScale/3.0)
                self.featureList[featureIndex].massMatchList[
                        massMatchIndex].seriesScore = seriesScore
                self.featureList[featureIndex].massMatchList[
                        massMatchIndex].totalScore += seriesScore
                            
        '''
        ####end of evan's edits


        
        # Get the difference of the highest total score and the second highest total score
        #     based on all matches within every feature

        for feature in self.featureList:

            massMatchWithTheHighestTotalScore = None
            massMatchWithThe2ndHighestTotalScore = None

            if len(feature.massMatchList) == 0:
                feature.differenceBetweenHighestAnd2ndHighestScore = "No matches"
            elif len(feature.massMatchList) == 1:
                feature.differenceBetweenHighestAnd2ndHighestScore = "Only 1 match"
            else:
                
                for massMatch in feature.massMatchList:
                    
                    if (massMatchWithTheHighestTotalScore == None and
                            massMatchWithThe2ndHighestTotalScore == None):
                        massMatchWithTheHighestTotalScore = massMatch
                    elif (massMatchWithTheHighestTotalScore != None and
                              massMatchWithThe2ndHighestTotalScore == None):
                        if massMatch.totalScore > massMatchWithTheHighestTotalScore.totalScore:
                            tempMassMatch = massMatchWithTheHighestTotalScore
                            massMatchWithTheHighestTotalScore = massMatch
                            massMatchWithThe2ndHighestTotalScore = tempMassMatch
                        else:
                            massMatchWithThe2ndHighestTotalScore = massMatch
                    elif (massMatchWithTheHighestTotalScore != None and
                              massMatchWithThe2ndHighestTotalScore != None):
                        if massMatch.totalScore > massMatchWithTheHighestTotalScore.totalScore:
                            tempMassMatch = massMatchWithTheHighestTotalScore
                            massMatchWithTheHighestTotalScore = massMatch
                            massMatchWithThe2ndHighestTotalScore = tempMassMatch
                        elif (massMatch.totalScore <= massMatchWithTheHighestTotalScore.totalScore and
                                  massMatch.totalScore > massMatchWithThe2ndHighestTotalScore.totalScore):
                            massMatchWithThe2ndHighestTotalScore = massMatch
                    else:
                        raise Exception("Error: 2nd highest scored mass match cannot be assigned before 1st highest scored mass match got assigned.")

                feature.differenceBetweenHighestAnd2ndHighestScore = (massMatchWithTheHighestTotalScore.totalScore
                                                                      - massMatchWithThe2ndHighestTotalScore.totalScore)
            
        # count number of decoy vs. non-decoy matches
        totalMatches = 0
        decoys = 0
        for feature in self.featureList:
            for match in feature.massMatchList:
                totalMatches += 1
                if match.decoyMassToMakeMatch == 11:
                    decoys += 1
                    
        self.totalMatches = totalMatches
        self.decoyMatches = decoys
        
        self.tolerancePPM = str(self.config["tolerancePPM"])
        self.scoringTolerancePPM = str(self.config["scoringTolerancePPM"])

        
        msFile = splitext(self.config["msFile"].split("/")[-1])[0]

            
        if self.config["isSequenceGivenByUser"] == False:
            listOfProteinFiles = self.config["listOfProteinFiles"]
            listOfProteinFileNames = [splitext(proteinFile.split("/")[-1])[0]
                                      for proteinFile in listOfProteinFiles]
            stringProteinFileNames = "_".join(listOfProteinFileNames)
            self.sequenceSource = "Uniprot_" + stringProteinFileNames
        else:
            listOfProteinFiles = ["User entered"]
            listOfProteinFileNames = ["User entered"]
            self.config["listOfProteinNames"].append("none")
            stringProteinFileNames = "_".join(listOfProteinFileNames)
            self.sequenceSource = "User-given"

        self.glycanSource = self.config["glycanSearchMethod"]

        self.digestion = " AND ".join([digestType
                                    for digestType
                                    in self.config["digestTypes"]])

        self.maxMissedCleavages = str(self.config["missedCleavages"])


        self.csvOutputFileName = msFile + "_" + self.config["msFormat"] + "_matching_and_fragment_matching"

        stringProteinFileNamesToPreappend = " AND ".join(listOfProteinFileNames)
        
        listOfProteinNames = self.config["listOfProteinNames"]
        stringProteinNamesToPreappend = " AND ".join(listOfProteinNames)
        
        self.specToPreappendToCsv = ("Sequence source:," + stringProteinFileNamesToPreappend
                                     + "\nProtein name(s):," + stringProteinNamesToPreappend
                                     + "\nGlycan source:," + self.glycanSource
                                     + "\nDigestion types:," + self.digestion
                                     + "\nMax miss cleavages:," + self.maxMissedCleavages
                                     + "\nMatching tolerance:," + self.tolerancePPM
                                     + "\nFragment tolerance:," + self.scoringTolerancePPM
                                     + "\nTotal Matches:," + str(self.totalMatches)
                                     + "\nDecoy Matches:," + str(self.decoyMatches)
                                     + "\n\n")
        


    def find_fragment_mass_matches(self):
    
        print "This function shouldn't be used any more, delete this"
        raise
        '''
        # Function Description
        #     Find fragment mass matches for every parent mass
        
        for feature in self.featureList:
            for match in feature.massMatchList:
                massScoreObject = MassScoreObject(match)
                match.listOfFragmentNeutralMassMatches = massScoreObject.score_feature_object(feature)'''


                
    def parse_mzData_file(self):
        """ this function will parse a mzData file """
        
        #open the file
        try:
            print("Parsing mzData file for XML formatting, this can take long")
            mzDatafile = open(self.config["msFile"], "rb")
            #mzDataDom = minidom.parse(mzDatafile)
        except IOError:
            print("ioError: mzDataFile not found")
        #ExpatError generated when parsing XML fails
        except expat.ExpatError:
            print("GP3Config.xml not in XML format:")
            for entry in sys.exc_info():
                print(entry)
                
        #getting the parent node for spectra
        """spectrumNodeList = mzDataDom.getElementsByTagName("spectrumList")[0]
        spectrumNodeListLength = int(spectrumNodeList.getAttribute('count').strip())
        spectrumNodeList = spectrumNodeList.getElementsByTagName("spectrum")
        assert spectrumNodeListLength == len(spectrumNodeList)"""
        
        featureList = []
        
        featureCount = 1
        addedFeatureCount = 0
        passedFeatureCount = 0
        level1Spectrum = []
        level1MzList = []
        previousSavedParent = -1
        previousUnUsedSpectrum = None
        
        #for spectrumNode in spectrumNodeList:
        for spectrumNode in self.iterate_through_mzData_spectrum_nodes(mzDatafile):
        
            featureCount += 1
            level, feature = self.parse_mzData_spectrum_node(spectrumNode, featureCount)
            if level == 1:
                level1Spectrum = feature
                level1MzList = [mz for mz, inten in level1Spectrum]
            else:
                newSpectrum = feature
                #should spectrum be saved at all?
                passes_diagnostic_check = self.check_for_diagnostic_ions(newSpectrum.xyzPeakList, newSpectrum.compoundNeutralMass)
                if passes_diagnostic_check:
                    passedFeatureCount += 1
                if not passes_diagnostic_check and \
                    math.fabs(previousSavedParent-newSpectrum.isolatedMzValue) > 0.01 :
                    previousUnUsedSpectrum = newSpectrum
                    continue
                
                ## should we just add this to the previous spectrum?
                if math.fabs(previousSavedParent-newSpectrum.isolatedMzValue) <= 0.01:
                    if self.config["applyFragmentNoiseFilter"] == True:
                        newSpectrum = self.plot_spectrum_and_filter_out_noisy_peaks_using_kernel_density_estimation(newSpectrum)
                    featureList[-1].combine_spectrum(newSpectrum.xyzPeakList, collisionEnergy=newSpectrum.collisionEnergy)
                    addedFeatureCount += 1
                    previousSavedParent = -1
                    continue
                else:
                    previousSavedParent = -1
                
                ##should the spectrum be saved as a new spectrum?
                if passes_diagnostic_check:
                    #note the fact it was saved:
                    previousSavedParent = newSpectrum.isolatedMzValue
                    #de noise the data 
                    if self.config["applyFragmentNoiseFilter"] == True:
                        newSpectrum = self.plot_spectrum_and_filter_out_noisy_peaks_using_kernel_density_estimation(newSpectrum)
                    #fix the monoisotopic peak
                    correctionFactor, nm, mz = self.fix_m0_isotope_selection(level1Spectrum, level1MzList, newSpectrum.isolatedMzValue, newSpectrum.compoundCharge, newSpectrum.compoundRtCenter)
                    newSpectrum.compoundNeutralMass = nm
                    newSpectrum.compoundMzValue = mz
                    newSpectrum.isotopeCorrection = correctionFactor
                    #save the spectrum
                    featureList.append(newSpectrum)
                    addedFeatureCount += 1
                    
                    #finally: should a previous unsaved spectrum be added?
                    if previousUnUsedSpectrum and \
                        math.fabs(previousUnUsedSpectrum.isolatedMzValue - newSpectrum.isolatedMzValue) <= 0.01:
                        if self.config["applyFragmentNoiseFilter"] == True:
                            previousUnUsedSpectrum = self.plot_spectrum_and_filter_out_noisy_peaks_using_kernel_density_estimation(previousUnUsedSpectrum)
                        featureList[-1].combine_spectrum(previousUnUsedSpectrum.xyzPeakList)
                        addedFeatureCount += 1
                        previousUnUsedSpectrum = None
                        continue
                    else:
                        previousUnUsedSpectrum = None
        
        #secret feature to write a mgf from the combined features:
        if "secret-write-mgf" in self.config:
            filename = self.config["secret-write-mgf"]
            self.write_new_mgf(featureList, filename)
            
        #deconvolute (or pseudodeconvolute)
        for feature in featureList:
            feature.convert_convoluted_peak_list_to_pseudo_deconvoluted(self.config["maxPossibleFragmentIonCharge"])
                
        print "SUCCESSFULLY RAN PARSE FUNCTION, press enter"
        print "Found {} passing features, added {} total features. ListLength = {}".format(passedFeatureCount, addedFeatureCount, len(featureList))
        time.sleep(2)
        self.featureList = featureList

    def iterate_through_mzData_spectrum_nodes(self, handle):
        spectrumBeginRe = """\A +<spectrum id"""
        spectrumBeginRe = re.compile(spectrumBeginRe)
        spectrumEndRe = """\A +</spectrum>"""
        spectrumEndRe = re.compile(spectrumEndRe)
        
        currentLine = ""
        
        for line in handle:
            
            if spectrumBeginRe.match(line):
                currentLine = line
                continue
                
            if spectrumEndRe.match(line):
                currentLine += line
                docNode = minidom.parseString(currentLine)
                spectrumNode = docNode.getElementsByTagName("spectrum")[0]
                #print currentLine
                #raw_input()
                yield spectrumNode
                continue
                
            currentLine += line
        
    
    def write_new_mgf(self, featureList, fileName ):
        "this simple function writes all feature-type objects back to file"
        
        with open(fileName, 'wb') as outputHandle:
            for feature in featureList:
                outputHandle.write("BEGIN IONS\n")
                if self.config.get("cidHigher"):
                    outputHandle.write("PEPMASS={}\n".format(feature.compoundMzValue))
                    outputHandle.write("CHARGE={}+\n".format(feature.compoundCharge))
                    outputHandle.write("TITLE={}-higheronly\n".format("isotopeCorrection={}".format(feature.isotopeCorrection)))
                    outputHandle.write("RTINSECONDS={}\n".format(feature.compoundRtCenter*60))
                    for xyzPeak in feature.xyzPeakListHigher:
                        outputHandle.write("{}\t{}\n".format(xyzPeak[0],xyzPeak[1]))
                    outputHandle.write("END IONS\n") 
                elif self.config.get("cidLower"):
                    outputHandle.write("PEPMASS={}\n".format(feature.compoundMzValue))
                    outputHandle.write("CHARGE={}+\n".format(feature.compoundCharge))
                    outputHandle.write("TITLE={}-loweronly\n".format("isotopeCorrection={}".format(feature.isotopeCorrection)))
                    outputHandle.write("RTINSECONDS={}\n".format(feature.compoundRtCenter*60))
                    for xyzPeak in feature.xyzPeakListLower:
                        outputHandle.write("{}\t{}\n".format(xyzPeak[0],xyzPeak[1]))
                    outputHandle.write("END IONS\n")  
                else:
                    outputHandle.write("PEPMASS={}\n".format(feature.compoundMzValue))
                    outputHandle.write("CHARGE={}+\n".format(feature.compoundCharge))
                    outputHandle.write("TITLE={}\n".format("isotopeCorrection={}-interlaced".format(feature.isotopeCorrection)))
                    outputHandle.write("RTINSECONDS={}\n".format(feature.compoundRtCenter*60))
                    for xyzPeak in feature.xyzPeakList:
                        outputHandle.write("{}\t{}\n".format(xyzPeak[0],xyzPeak[1]))
                    outputHandle.write("END IONS\n")  
        raise EnvironmentError("disable the output dummy!")                
                
    def fix_m0_isotope_selection(self, level1Spectrum, level1MzList, parentIon, charge, rt):
        """this will fix the selection of m0
           returns corrected c13 offset
           if offset == 0, no correction is required
           valid returns are [0, -1, -2, -3, -4, -5]"""
        
        c13Diff = 1.00335
        
        def binary_search(a, x, lo=0, hi=None, side="Left"):   # can't use a to specify default for hi
            hi = hi if hi is not None else len(a) # hi defaults to len(a)   
            if side=="Left":
                pos = bisect.bisect_left(a,x,lo,hi)          # find insertion position
            else:
                pos = bisect.bisect_right(a,x,lo,hi) 
            return (pos if pos != hi else -1) # don't walk off the end

        
        #the following distributions were calculated using http://prospector.ucsf.edu/prospector/cgi-bin/msform.cgi?form=msisotope
        # each composition is the core N-glycan + averagine, gaps of roughly 500Da were used.
        distributionDict = {500: [100, 20] ,   # this is arbitrary
                            1493:[100,45.6,15.67],    # C38 H63 N3 O27 S0
                            1930:[100,70.51,31.23,10.33],  # C58 H94 N8 O33 S0
                            2480:[98.12,100,58.45,24.95],  # C83 H133 N15 O40 S0
                            3045:[74.98,100,73.51,38.75],  # C108 H172 N22 O48 S0
                            3582:[61.11,100,88.17,55.01],  # C132 H211 N29 O55 S0
                            4149:[48.29,94.38,100,75.25],  # C157 H249 N35 O62 S1
                            4714:[36.37,82.52,100,85.28],  # C182 H288 N42 O70 S1
                            5251:[28.61,73.58,100,94.95,70.42],  # C206 H327 N49 O77 S1
                            5817:[21.77,62.84,95.11,100,81.77],  # C231 H366 N56 O85 S1
                            6366:[16.25,52.01,86.74,100,89.33],  # C256 H405 N63 O92 S1
                            6920:[12.28,43.05,78.85,100,98.39],  # C280 H443 N69 O99 S2
                            7486:[9.07,34.67,68.85,94.35,100],  # C305 H482 N76 O107 S2
                            8035:[6.76,27.95,59.83,88.06,100],  # C330 H521 N83 O114 S2
                            8588:[5.17,22.97,52.61,82.67,100],  # C354 H560 N90 O122 S2
                            9138:[3.8,18.06,44.17,73.92,95.06,100],  # C379 H599 N97 O129 S2
                            9704:[2.73,13.84,36.13,64.62,88.85,100],  # C404 H637 N103 O136 S3
                            10269:[2.04,11,30.42,57.54,83.53,99.12,100],  # C429 H676 N110 O144 S3
                            10807:[1.51,8.58,25.01,49.78,75.93,94.55,100],  # C453 H715 N117 O151 S3
                            11356:[1.12,6.71,20.61,43.13,69.09,90.24,100]}   # C478 H754 N124 O158 S3

        distributionKeyList = [500, 1493, 1930, 2480, 3045, 3582, 
                               4149, 4714, 5251, 5817, 6366,
                               6920, 7486, 8035, 8588, 9138, 
                               9704, 10269, 10807, 11356]

        
        correctionScopeDict = {1500:1, 3000:2, 4500:3}
        correctionScopelist = [ 1500, 3000, 4500]
        
        roughNeutralMass = parentIon*charge
        correctionScope = 4  # default value
        for correctionScopekey in correctionScopelist:
            if roughNeutralMass < correctionScopekey:
                correctionScope = correctionScopeDict[correctionScopekey]
                break
        
        distKeyIdx = binary_search(distributionKeyList, roughNeutralMass)
        distKey = distributionKeyList[distKeyIdx]
        
        isotopePattern = distributionDict[distKey]
        lenIsotopePattern = len(isotopePattern)
        lenSpectrum = len(level1MzList)
        
        ppmError = min( parentIon*(20.0/1000000), 0.3*(c13Diff/charge) )
        
        correctionFactorDeltaList = []
        for correctionFactor in [-1*L for L in range(correctionScope)]:
            #this makes the miniList
            minMz = (parentIon+((c13Diff*correctionFactor)/charge))-ppmError
            maxMz = (parentIon+((c13Diff*correctionFactor)/charge)+((c13Diff*lenIsotopePattern)/charge))+ppmError
            minMzIdx = binary_search(level1MzList, minMz, lo=0)
            maxMzIdx = binary_search(level1MzList, maxMz, lo=0, side="Right")
            miniList = level1Spectrum[minMzIdx:maxMzIdx]
            miniListMz = [mz for mz,inten in miniList]
            if len(miniList) == 0:
                neutralMass = parentIon*charge - self.adduct*charge
                massToCharge = parentIon
                return 0, neutralMass, massToCharge
                
            #this makes the intensity ajusted isotopePattern
            miniListMaxInten = max([inten for mz, inten in miniList])
            maxPattern = max(isotopePattern)
            miniList = [(minMz,miniInt*(maxPattern/miniListMaxInten)) for minMz, miniInt in miniList]
            
            diffList = []
            for isotope in range(lenIsotopePattern):
                #make the microlist:
                isoMass = parentIon+(c13Diff*(isotope+correctionFactor))/charge
                microListIdxL = binary_search(miniListMz, isoMass-ppmError, lo=0)
                microListIdxR = binary_search(miniListMz, isoMass+ppmError, lo=0, side="Right")
                if microListIdxR == -1:
                    microListIdxR = len(miniList)
                if microListIdxL == microListIdxR or microListIdxL == -1:
                    microList = [(isoMass, 0)]
                else:
                    microList = miniList[microListIdxL:microListIdxR]
                
                isoInten = isotopePattern[isotope]
                scoreAndRankFunction = lambda mzValInten: math.sqrt((0.5*maxPattern*(mzValInten[0]-isoMass))**2 + (mzValInten[1]-isoInten)**2)
                matchingIon = min(microList,key=scoreAndRankFunction)
                differenceFound = scoreAndRankFunction(matchingIon)
                diffList.append(differenceFound)
            #fix weight of first 2 isotopes then do geometric mean of differences
            diffList = [math.fabs(d) for d in diffList]
            diffList[0] = diffList[0]**2
            diffList[1] = diffList[1]**2
            geomMeanFixer = 3

            geometricDiff = np.prod([diff+0.2 for diff in diffList])**(1.0/(lenIsotopePattern+geomMeanFixer))
            correctionFactorDeltaList.append(geometricDiff)

        #correctionFactorDeltaList = [correctionFactorDeltaList[i]*(1-(i*0.1)) for i in range(len(correctionFactorDeltaList))]
        if min(correctionFactorDeltaList) < 5:
            bestIsotopeCorrection = -1* correctionFactorDeltaList.index(min(correctionFactorDeltaList))
        else:
            bestIsotopeCorrection = 0
        neutralMass = parentIon*charge + (c13Diff*bestIsotopeCorrection) - self.adduct*charge
        massToCharge = parentIon+((c13Diff*bestIsotopeCorrection)/charge)
        #temporary plotting functions
        """ """

        
        if "make_test_isotope_plot" in self.config \
            and self.config["make_test_isotope_plot"] \
            and roughNeutralMass>2500 and bestIsotopeCorrection != 0:
            
            if "fix_m0_isotope_pltFigure" in self.config:
                fig = self.config["fix_m0_isotope_pltFigure"]
                ax = self.config["fix_m0_isotope_ax"]
        
            
            import matplotlib.pyplot as plt
            
            if "fix_m0_isotope_pltFigure" not in self.config:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                self.config["fix_m0_isotope_pltFigure"] = fig
                self.config["fix_m0_isotope_ax"] = ax
                plt.ion()
            
            ax.cla()
            minMz = (parentIon+((c13Diff*bestIsotopeCorrection)/charge))-ppmError
            maxMz = (parentIon+((c13Diff*bestIsotopeCorrection)/charge)+((c13Diff*lenIsotopePattern)/charge))+ppmError
            minMzIdx = binary_search(level1MzList, minMz, lo=0, hi=lenSpectrum)
            maxMzIdx = binary_search(level1MzList, maxMz, lo=0, hi=lenSpectrum)
            maxMzIdx = min(maxMzIdx, lenSpectrum)
            miniList = level1Spectrum[minMzIdx:maxMzIdx]
            if len(miniList) == 0:
                miniList = [(minMz/0.9998,0)]
            miniListMaxInten = max([inteni for mzi, inteni in miniList])
            maxPattern = max(isotopePattern)
            miniList = [(minMz,miniInt) for minMz, miniInt in miniList]
            mz, inten = zip(*miniList)
            zeros = [0 for i in mz]
            collection2 = ax.vlines(mz, zeros, inten, 'k', linewidths=2.5)
            adjustedPatternMass = [parentIon+(c13Diff*(isotope+bestIsotopeCorrection))/charge for isotope in range(len(isotopePattern))]
            adjustedPatternInten = [intenBuilder*(miniListMaxInten/maxPattern) for intenBuilder in isotopePattern]
            adjustedPatternZeros = [0 for isotope in range(len(isotopePattern))]
            collection1 = ax.vlines(adjustedPatternMass, adjustedPatternZeros, adjustedPatternInten, 'DarkSlateGray', linestyles='dotted', linewidths=9)
            if len(mz) > 0:
                pattern_intenParent = inten[binary_search(mz, parentIon-ppmError)]
                pattern_intenFix = inten[binary_search(mz, parentIon-ppmError+(bestIsotopeCorrection/float(charge)))]
            else:
                pattern_intenParent = 0
                pattern_intenFix = 0
            machineAsigned = max(pattern_intenParent, adjustedPatternInten[int(-1*bestIsotopeCorrection)])*1.05
            fixedAsigned = max(pattern_intenFix, adjustedPatternInten[0])*1.1
            dots1 = ax.plot([parentIon], [machineAsigned], 'r.', markersize=14, label="machine assigned")
            dots2 = ax.plot([massToCharge], [fixedAsigned], 'b.', markersize=14, label="corrected")
            l1 = ax.legend([collection1,collection2,dots1[0], dots2[0]], ["Simulated MS","Actual MS","Machine assigned m0","Corrected m0"], loc=1, numpoints=1)
            plt.title('Selected ion ={},  z={}, rt={}, delta={}'.format(parentIon,charge, rt,min(correctionFactorDeltaList)))
            plt.xlabel('m/z')
            plt.ylabel('inten')
            if "fix_m0_isotope_shown" not in self.config:
                raw_input("----- first plot ready ----------")
                plt.show()
                self.config["fix_m0_isotope_shown"] = True
            else:
                raw_input("next plot ready")
                fig.canvas.draw()
                
        """ """
        
        return bestIsotopeCorrection, neutralMass, massToCharge
            
        
    def parse_mzData_spectrum_node(self, spectrumNode, featureCount):
        """this will parse a node containing a spectrum""" 
        
        #set up some function specific abstractions
        def decode_scalar_list(encodedText, endian, length, precision):
            if endian=="little":
                fmtStrA = "<"
            elif endian == "big":
                fmtStrA = ">"
            else:
                fmtStrA = ""
            
            if precision == "64":
                fmtStrC = "d"
                bytes = 8
            elif precision == "32":
                fmtStrC = "f"
                bytes = 4

            dividerCounts = 500*32       # 32 packs 3 structures of 8 byte or two of 4 byte
            dividerLen = int(dividerCounts*(3/4.0)/float(bytes))
            fmtCount, fmtRemainder = divmod(len(encodedText), dividerCounts)
            remainingEntries = int(fmtRemainder*(3/4.0)*(1.0/bytes))
            returnList = []
            for idx in range(fmtCount):
                idx0 = dividerCounts*idx
                idx1 = dividerCounts*(idx+1)
                fmtStr = fmtStrA + str(dividerLen) + fmtStrC
                subEncodedText = encodedText[idx0:idx1]
                binary_data = base64.b64decode(subEncodedText)
                valList = struct.unpack(fmtStr,binary_data)
                returnList.extend(valList)
            
            fmtStr = fmtStrA + str(remainingEntries) + fmtStrC
            base = dividerCounts*fmtCount
            subEncodedText = encodedText[base:base+fmtRemainder]
            binary_data = base64.b64decode(subEncodedText)
            valList = struct.unpack(fmtStr,binary_data)
            returnList.extend(valList)
            
            return returnList
        
        #get level, if level=1, this will be shorter
        spectrumInstrumentNode = spectrumNode.getElementsByTagName("spectrumInstrument")[0]
        msLevel = int(spectrumInstrumentNode.getAttribute('msLevel').strip())
        
        #get all ms data
        mzArrayBinaryDataNode = spectrumNode.getElementsByTagName("mzArrayBinary")[0]
        mzArrayBinaryDataNode = mzArrayBinaryDataNode.getElementsByTagName("data")[0]
        precision = mzArrayBinaryDataNode.getAttribute('precision').strip()
        endian= mzArrayBinaryDataNode.getAttribute('endian').strip()
        length = int(mzArrayBinaryDataNode.getAttribute('length').strip())
        if length > 0:
            encodedText = mzArrayBinaryDataNode.firstChild.data
            mzValueList = decode_scalar_list(encodedText, endian, length, precision)
        else:
            mzValueList = []
        
        intenBinaryDataNode = spectrumNode.getElementsByTagName("intenArrayBinary")[0]
        intenBinaryDataNode = intenBinaryDataNode.getElementsByTagName("data")[0]
        precision= intenBinaryDataNode.getAttribute('precision').strip()
        endian= intenBinaryDataNode.getAttribute('endian').strip()
        length= int(intenBinaryDataNode.getAttribute('length').strip())
        if length > 0:
            encodedText = intenBinaryDataNode.firstChild.data
            intenValueList = decode_scalar_list(encodedText, endian, length, precision)
        else:
            intenValueList = []
        #return data if no more meta data is needed
        if length < 1:
            cvParamsToGather = {"MassToChargeRatio":None,
                            "CollisionEnergy":None,
                            "MassToChargeRatio":None,
                            "ChargeState":None,
                            "TimeInMinutes":None}
            for cvPram in spectrumNode.getElementsByTagName("cvParam"):
                name = cvPram.getAttribute('name').strip() 
                if name in cvParamsToGather:
                    cvParamsToGather[name] = float(cvPram.getAttribute('value').strip())
            newFeature = Feature(featureCount,
                         -1*featureCount, #compound neutral mass
                         -1*featureCount, # m/z
                         1,
                         cvParamsToGather["TimeInMinutes"], #rt
                         None, # compoundRtTuple unknown
                         None, # compoundIntensity unknown
                         collisionEnergy=0 )
            newFeature.xyzPeakList = []
            return 2, newFeature
        
        if msLevel == 1:
            return msLevel, zip(mzValueList, intenValueList)
        
        
        cvParamsToGather = {"MassToChargeRatio":None,
                            "CollisionEnergy":None,
                            "MassToChargeRatio":None,
                            "ChargeState":None,
                            "TimeInMinutes":None}
                           
        for cvPram in spectrumNode.getElementsByTagName("cvParam"):
            name = cvPram.getAttribute('name').strip() 
            if name in cvParamsToGather:
                if name == "ChargeState":
                    cvParamsToGather[name] = int(cvPram.getAttribute('value').strip())
                else: 
                    cvParamsToGather[name] = float(cvPram.getAttribute('value').strip())

        parentCharge = cvParamsToGather["ChargeState"] if cvParamsToGather["ChargeState"] is not None else 1
        compoundNeutralMass = (cvParamsToGather["MassToChargeRatio"] * parentCharge) - parentCharge*self.adduct
        
        xyzPeakList = []
        for ion in zip(mzValueList, intenValueList):
            x = ion[0]
            y = ion[1]
            multiplier = 1.0  #this may need to be changed
            xyzTuple = (x, y, 0, multiplier)
            xyzPeakList.append(xyzTuple)
        
        newFeature = Feature(featureCount,
                             compoundNeutralMass,
                             cvParamsToGather["MassToChargeRatio"] ,
                             parentCharge,
                             cvParamsToGather["TimeInMinutes"],
                             None, # compoundRtTuple unknown
                             None, # compoundIntensity unknown
                             collisionEnergy=cvParamsToGather["CollisionEnergy"] )
        newFeature.xyzPeakList = xyzPeakList
        return msLevel, newFeature
    
    def parse_mgf(self):

        # Function Description
        #     Read the mgf file, which is the output from MassHunter
        # Function Input:
        #     The mgf file from MassHunter
        
        #setting charges list
        maximumFragmentCharge = self.config["maxPossibleFragmentIonCharge"]
        if maximumFragmentCharge > 0:
            charges = [charge + 1 for charge in range(maximumFragmentCharge)]
        else:
            charges = [0]
            
        inFile = open(self.config["msFile"], 'r')
        print("Reading " + self.config["msFile"] + "...")
        
        featureList = []
        
        featureCounter = 0
        
        for line in inFile:
            
            line = line.strip()
            #print(line.split())
            if "BEGIN" in line and "IONS" in line:
                isotopeCorrection = 0
                featureCounter += 1
                compoundRtCenter = featureCounter  #default value when RTINSECONDS is absent
                print("BEGIN IONS tag seen. Start of Feature " + str(featureCounter))
                xyzPeakList = []
                compoundCharge = 1 # default charge value when CHARGE is absent
                if "SKIP" not in line:
                    skipSpectrum = False
                if "SKIP" in line:
                    skipSpectrum = True
                    
            elif "PEPMASS" in line:
                massString = line.split('=')[1]
                compoundMzValue = float(massString.split()[0])    
                compoundNeutralMass = compoundMzValue - self.adduct #assumes charge = 1, this assumption is fixed later with 'CHARGE' line
                
            elif "CHARGE" in line:
                compoundCharge = int("".join([char
                                              for char in line
                                              if char.isdigit() == True]))
                compoundNeutralMass = compoundMzValue * compoundCharge - self.adduct * compoundCharge
                
            elif "TITLE" in line:
                if "isotopeCorrection" in line:
                    for potentialDigit in line:
                        if potentialDigit.isdigit():
                            isotopeCorrection = -1*int(potentialDigit)
                else:
                    pass # old version uses agilent specific times - compoundRtCenter = float(line.split()[-2])
                
            elif "END" in line and "IONS" in line:
                #check xyzPeakList
                validByDiagnosticPeakChecker = self.check_for_diagnostic_ions(xyzPeakList, compoundNeutralMass)
                print("END IONS tag seen. End of Feature " + str(featureCounter) + "\n")
                
                if skipSpectrum == True:
                    print("Spectrum skipped, not added to the list of features.\n")
                    
                newFeature = Feature(featureCounter,
                                     compoundNeutralMass,
                                     compoundMzValue,
                                     compoundCharge,
                                     compoundRtCenter,
                                     None, # compoundRtTuple unknown
                                     None # compoundIntensity unknown
                                     )

                newFeature.xyzPeakList = xyzPeakList
                newFeature.isotopeCorrection = isotopeCorrection

                if self.config["applyFragmentNoiseFilter"] == True and validByDiagnosticPeakChecker:
                    newFeature = self.plot_spectrum_and_filter_out_noisy_peaks_using_kernel_density_estimation(
                                                                newFeature)
                
                if skipSpectrum == False and validByDiagnosticPeakChecker:
                    featureList.append(newFeature)
                
                
                del compoundNeutralMass
                del compoundMzValue
                del compoundCharge
                del compoundRtCenter
                del xyzPeakList
                skipSpectrum = False
                
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
                for charge in charges:
                    xyzTuple = (x, y, charge, multiplier)
                    xyzPeakList.append(xyzTuple)

        self.featureList = featureList
        
        #used for testing:
        print("found %d features that have diagnostic ions out of %d features in file"%(len(featureList), featureCounter))
        time.sleep(2)

    
    def check_for_diagnostic_ions(self, xyzPeakList, compoundNeutralMass):
        #given an xyzPeakList, this function will return true if the diagnostic peaks are present
        #  diagnostic peaks are found within tolerances used for scoring. If the sum of the signal
        #  of the diagnostic peaks is greater than the maximum peak multiplied by some adjustment factor 
        #  (usually a fraction), then this function will evaluate as 'True'
        #  In addition to checking intensity, a simple check is made such that at least 2 diagnostic peaks
        #  must be present for this routine to evaluate as True
        
        #first, check if we need to run this algorithm, return True if user has not enabled it
        if self.config["applyDiagnosticIonFilter"] == False:
            return True
        
        #next, ignore empty lists or lists with too few ions
        if len(xyzPeakList) <= max(0,self.config["diagnosticNumberRequired"]):
            return False
            
        #setting diagnostic peaks, this depends on charge of input files:
        #   deconvoluted input files need different diagnostic ions
        maximumFragmentCharge = self.config["maxPossibleFragmentIonCharge"]
        if maximumFragmentCharge > 0:
            diagnosticPeakList = [# Hex-minus-3H2O-plus-H+:
                                  # C6H12O6 - 3H2O + H = C6H7O3
                                  127.0395190913,
                                  # Hex-minus-H2O-plus-H+:
                                  # C6H12O6 - H2O + H = C6H11O5
                                  163.0606484641,
                                  # HexNAc-minus-3H2O-plus-H+:
                                  # C8H15NO6 - 3H2O + H = C8H10NO3
                                  168.066068193,
                                  # HexNAc-minus-2H2O-plus-H+:
                                  # C8H15NO6 - 2H2O + H = C8H12NO4
                                  186.0766328794,
                                  # Hex-minus-3H2O-plus-H+:
                                  # C6H12O6 - 3H2O + H = C6H7O3
                                  204.0871975658,
                                  # Neu5Ac-minus-2H2O-plus-H+:
                                  # C11H19NO9 - 2H2O + H = C11H16NO7
                                  274.0926768743,
                                  # Neu5Ac-minus-H2O-plus-H+:
                                  # C11H19NO9 - H2O + H = C11H18NO8
                                  292.1032415607,
                                  # Hex-HexNAc-minus-H2O-plus-H+:
                                  # C6H12O6 + C8H15NO6 - 2H2O + H = C14H24NO10
                                  366.1400209978,
                                  compoundNeutralMass - 162.052823432 + 1.0078250321,
                                  compoundNeutralMass - 203.0793725337 + 1.0078250321,
                                  compoundNeutralMass - 291.0954165286 + 1.0078250321,
                                  compoundNeutralMass - 365.1321959657 + 1.0078250321,
                                  compoundNeutralMass - 453.1482399606 + 1.0078250321,
                                  (compoundNeutralMass - 162.052823432 + 2*1.0078250321)/2.0,
                                  (compoundNeutralMass - 203.0793725337 + 2*1.0078250321)/2.0,
                                  (compoundNeutralMass - 291.0954165286 + 2*1.0078250321)/2.0,
                                  (compoundNeutralMass - 365.1321959657 + 2*1.0078250321)/2.0,
                                  (compoundNeutralMass - 453.1482399606 + 2*1.0078250321)/2.0]
        else: 
            # same as masses from above, reduced by mass of one proton
            diagnosticPeakList = [126.0316940592,
                                  162.052823432,
                                  167.0582431609,
                                  185.0688078473,
                                  203.0793725337,
                                  273.0848518422,
                                  291.0954165286,
                                  365.1321959657,
                                  #glycopeptide-minus-Hex
                                  compoundNeutralMass - 162.052823432,
                                  #glycopeptide-minus-HexNAc
                                  compoundNeutralMass - 203.0793725337,
                                  #glycopeptide-minus-SialicAcid,
                                  compoundNeutralMass - 291.0954165286,
                                  #glycopeptide-minus-Hex-HexNAc,
                                  compoundNeutralMass - 365.1321959657,
                                  #glycopeptide-minus-Hex-SialicAcid,
                                  compoundNeutralMass - 453.1482399606]
        
        
        maxPeak = max(xyzPeakList, key=lambda xyzPeakList: xyzPeakList[1])
        maxPeakIntensity = maxPeak[1]
        scoreTolerance = 1.5*self.config["scoringTolerance"]
        scoreToleranceAbs= 300.0 * scoreTolerance  #a rough tolerance used for low mass fragments
        
        
        scorePeakIntensity = 0
        scorePeakCount = 0

        for diagnosticPeak in diagnosticPeakList:
            for xyzPeak in xyzPeakList:
                realPeakLower = xyzPeak[0] - scoreToleranceAbs
                realPeakHigher = xyzPeak[0] + scoreToleranceAbs
                if realPeakLower <=  diagnosticPeak <= realPeakHigher:
                    scorePeakIntensity += xyzPeak[1]
                    scorePeakCount += 1

        if (scorePeakIntensity >= maxPeakIntensity * self.config["diagnosticIonCutoff"]
                and scorePeakCount >= self.config["diagnosticNumberRequired"]):
            return True
        else:
            print("Not enough diagnostic ions found, spectrum will not be appended")
            return False
        
        
    
    def parse_cef_xml(self):

        # Function Description:
        #     Read the cef file, which is the output from MassHunter Qual 5
        #         (in the wet chemistry lab)
        # Function Input:
        #     The cef file from MassHunter Qual 5
        
        inFile = open(self.config["msFile"], 'r')
        print("Reading " + self.config["msFile"] + "...")
        
        document = minidom.parse(inFile)
        
        # find compounds in CEF file
        # take care of base case; given empty file or nonstandard format.
        compoundNodes = document.getElementsByTagName("Compound")
        if not compoundNodes:
            compoundNodes = document.getElementsByTagName("compound")
            if not compoundNodes:
                raise IOError("no compounds found in CEF file")

        featureCounter = 0
        
        #iterate through compounds
        for compound in compoundNodes:
            
            #get polarity of MS spectrum polarity is either (1) for positive or (-1) for negative
            details = compound.getElementsByTagName("MSDetails")[0]
            polarity = details.getAttribute('p').strip()
            if polarity == '+':
                polarity = 1
            else:
                polarity = -1
                
            #get charge Charge of parent ion
            MzOfInterestNode = compound.getElementsByTagName("MzOfInterest")[0]
            charge = MzOfInterestNode.getAttribute("z").strip()
            if charge == '':
                compoundCharge = 1*polarity
            else:
                compoundCharge = abs(int(charge))*polarity
            
            #get mz, neutral mass
            mzNodes = MzOfInterestNode.getElementsByTagName("mz")
            mzList = []
            for mz in mzNodes:
                mzValue = float(mz.firstChild.data)
                mzList.append(mzValue)
            compoundMzValue = float(sum(mzList)) / len(mzList)
            compoundNeutralMass = compoundMzValue*compoundCharge - self.adduct*compoundCharge
            
            # get retention time (rt)
            rtNode = compound.getElementsByTagName("RTRange")[0]
            rtMin = float(rtNode.getAttribute("min").strip())
            rtMax = float(rtNode.getAttribute("max").strip())
            compoundRtTuple = (rtMin, rtMax)
            
            # get intensity and center of RT from cef file
            locationNode = compound.getElementsByTagName("Location")[0]
            compoundIntensity = int(locationNode.getAttribute('y').strip())
            compoundRtCenter = float(locationNode.getAttribute('rt').strip())

            featureCounter = featureCounter + 1
            
            # make new feature and append
            newFeature = Feature(featureCounter,
                                 compoundNeutralMass,
                                 compoundMzValue,
                                 compoundCharge,
                                 compoundRtCenter,
                                 compoundRtTuple,
                                 compoundIntensity
                                 )
            
            featureObject = self.get_compound_mass_spec(compound, newFeature)
            
            self.featureList.append(featureObject)
                                
            
        
    def get_compound_mass_spec(self,
                               compoundNode,
                               featureObject):

        # Function Description:
        #     Read and get the compound parent ion mass and charge
        # Function Input:
        #     The CEF file from MassHunter Qual 5
        
        peakList = []
        
        msPeaks = compoundNode.getElementsByTagName("MSPeaks")
        xyData = compoundNode.getElementsByTagName("XYData")
        
        if msPeaks:
            msPeaks = msPeaks[0]
            peakNodes = msPeaks.getElementsByTagName("p")
            numberOfOriginalPeaks = len(peakNodes)
            originalIntensityList = []
            
            for peak in peakNodes:
                mass = float(peak.getAttribute('x').strip())
                height = float(peak.getAttribute('y').strip())
                z = peak.getAttribute('z').strip()
                originalIntensityList.append(height)
                
                if z:
                    charge = int(z)
                    multiplier = 1.0
                    xyzPeak = (mass, height, charge, multiplier)
                    xyzPeakList = [xyzPeak]
                else:
                    maxPossibleFragmentIonCharge = self.config["maxPossibleFragmentIonCharge"]
                    charges = [charge
                               for charge in range(1, maxPossibleFragmentIonCharge + 1)]
                    multiplier = 1.0 / maxPossibleFragmentIonCharge
                    xyzPeakList = []
                    for charge in charges:
                        xyzPeakSingle = (mass, height, charge, multiplier)
                        xyzPeakList.append(xyzPeakSingle)
                
                peakList.extend(xyzPeakList)

            featureObject.numberOfOriginalPeaks = numberOfOriginalPeaks
            featureObject.originalIntensityList = originalIntensityList
            featureObject.xyzPeakList = peakList
            
            if self.config["applyFragmentNoiseFilter"] == True:
                featureObject = self.plot_spectrum_and_filter_out_noisy_peaks_using_kernel_density_estimation(
                                                                featureObject)
                
            #featureObject.xyzPeakList.sort()
            
        elif xyData:
            print("xyData not supported (yet)!!!!!")
            
        else:
            print("file missing tandem data or tandem data type not supported")
            
        return featureObject

    

    def plot_spectrum_and_filter_out_noisy_peaks_using_kernel_density_estimation(
                                                                self,
                                                                featureObject):

        # Function Description:
        #     Do a Kernel Density Estimation of all peak intensities within every feature,
        #     Get the probability density function based on all peak intensities
        #     Get the maximum value of the probability density function
        #     Discard all peak intensities that falls below alpha times the maximum value
        #     where alpha still needs to be defined by the user
        #     Kernel Density Estimation:
        #     http://en.wikipedia.org/wiki/Kernel_density_estimation
        # Function Input:
        #     One feature object
        # Function Output:
        #     A selective subset of list of peaks

        xyzPeakList = featureObject.xyzPeakList
        textDisplay = ("MS/MS for Feature " + str(featureObject.featureIndex)
                            + " for parent ion mz = " + str(featureObject.compoundMzValue)
                            + " at RT = " + str(featureObject.compoundRtCenter))
        newXyzPeakList = []
        successfulNoiseFiltering = True # assume success, failure will negate this
        lenXyzPeakList = len(xyzPeakList)
        noiseCutoffIntensity = None
        
        if lenXyzPeakList > 20:
        
            groups = 5
            print("Performing kernel density estimation on " + textDisplay)
            segmentLength = int(len(xyzPeakList)/groups)
            segCounter = 0
            
            while segCounter < groups:
                try:
                    #set up intensity List
                    indexZero = segCounter*segmentLength
                    if segCounter < (groups - 1):
                        indexOne = (segCounter+1)*segmentLength
                    else:
                        indexOne = lenXyzPeakList
                    xyzPeakListShortCopy = xyzPeakList[indexZero:indexOne]
                    mzList = [xyzPeak[0] for xyzPeak in xyzPeakListShortCopy]
                    intensityList = [xyzPeak[1] for xyzPeak in xyzPeakListShortCopy]
                    
                    
                    #running kernel density estimation
                    kde_function = gaussian_kde(intensityList)
                    maxIntensity = max(intensityList)
                    
                    
                    #used for finding true max precision = 0.5 ppt
                    xPointsKDE = np.linspace(0, maxIntensity, 2000)
                    yPointsKDE = kde_function(xPointsKDE)
                    maxKDE = numpy.max(yPointsKDE)
                    max_inten_KDE = xPointsKDE[yPointsKDE.tolist().index(maxKDE)]
                    
                    # filter out xyzPeaks
                    # that have a intensity being less than the cutoff intensity
                    noiseCutoffIntensity = max_inten_KDE * self.config["noiseCutoffThreshold"]
                    
                    newXyzPeakListEndSegment = [xyzPeak
                                               for xyzPeak in xyzPeakListShortCopy
                                               if xyzPeak[1] >= noiseCutoffIntensity]
                    
                    newXyzPeakList.extend(newXyzPeakListEndSegment)
                    
                    
                    # PLOTTING FUNCTION
                    #  uncomment the code block below to plot KDE
                    """
                    import matplotlib.pyplot as plt
                    xPointsKDE_plot = np.linspace(0, maxIntensity, 200)
                    yPointsKDE_plot = kde_function(xPointsKDE_plot)
                    
                    left, width = 0.1, 0.57
                    bottom, height = 0.1, 0.75
                    left_h = left + width + 0.065

                    rect_scatter = [left, bottom, width, height]
                    rect_histy = [left_h, bottom, 0.2, height]

                    # start with a rectangular Figure
                    f = plt.figure(num=1, figsize=(8,6))
                    f.text(0.5, 0.95, textDisplay, ha="center")
                    axScatter = plt.axes(rect_scatter)
                    axScatter.set_title("Mass spectrum")
                    axHist = plt.axes(rect_histy)
                    
                    # now determine nice limits for plotting and binning:
                    ymax = np.max(intensityList)
                    binwidth = 50
                    binwidth = ymax/40
                    binLimit = ( int(ymax/binwidth) + 1) * binwidth
                    
                    # determine how to plot the cutoff line
                    cutoffVlinesX = (0, max(mzList))
                    cutoffVlinesY = (max_inten_KDE, max_inten_KDE)
                        
                    cutoffHistogramX = (max_inten_KDE, max_inten_KDE)
                    #print cutoffHistogramX
                    cutoffHistogramY = (0, maxKDE)
                    #print cutoffHistogramY
                    
                    #set up histogram and probability density function derived from KDE
                    bins = numpy.arange(0, binLimit + binwidth, binwidth)
                    axHist.hist(intensityList,
                                bins=bins,
                                orientation='horizontal',
                                normed=True,
                                alpha=0.25)
                    axHist.plot(yPointsKDE_plot, xPointsKDE_plot, "b-")
                    axHist.plot(cutoffHistogramY, cutoffHistogramX, 'b-', lw=2)
                    cutoffHistogramX = [noiseCutoffIntensity,noiseCutoffIntensity]
                    axHist.plot(cutoffHistogramY, cutoffHistogramX, 'r--', lw=2)
                    axHist.set_ylim( (0,ymax) )
                    axHist.set_title("PDF")
                    xticks = np.linspace(0,axHist.get_xlim()[1],3)
                    axHist.set_xticks(xticks)
                    axHist.get_yaxis().set_visible(False)
                    
                    # the vlines plot:
                    intenZeros = np.zeros(len(mzList))
                    axScatter.vlines(mzList ,intenZeros, intensityList)
                    axScatter.plot(cutoffVlinesX, cutoffVlinesY, 'b-', lw=2)
                    cutoffVlinesY = [noiseCutoffIntensity,noiseCutoffIntensity]
                    axScatter.plot(cutoffVlinesX, cutoffVlinesY, 'r--', lw=2)
                    axScatter.set_xlim( (0, 2500) )
                    axScatter.set_ylim( (0, ymax) )
                    
                    plt.show()"""
                    

                # Sometimes doing gaussian_kde over a list of values gives an error:
                # numpy.linalg.linalg.LinAlgError: Singular matrix
                # So I catch this error to deal with it
                # When the program sees a singular matrix
                # which is often due to a single fragment intensity value in a feature
                # we skip this KDE and cutoff procedure to keep the peak list intact
                except numpy.linalg.linalg.LinAlgError as err:
                    if "Singular matrix" in err.message:
                        if self.config["testEnabled"] == True:
                            print("\nIntensity List: " + str(intensityList))
                            print("Avoided throwing numpy.linalg.linalg.LinAlgError")
                        print("Cutoff intensity: " + str(intensityList[0]))
                        print("Number of xyzPeaks originally:      " + str(len(xyzPeakList)))
                        print("Unchanged.\n")
                        successfulNoiseFiltering = False
                        
                    else:
                        raise
                
                #incriment segment counter no matter what
                segCounter += 1
                
        else:
            successfulNoiseFiltering = False
            print("No cuttoff intensity.")
            print("Number of xyzPeaks                  " + str(len(xyzPeakList)))
            print("No filtering performed.")
        
        
        if self.config["testEnabled"] == True:
            print("Cutoff intensity: " + str(noiseCutoffIntensity))
            print("Number of xyzPeaks originally:      " + str(lenXyzPeakList))
            print("Number of xyzPeaks after filtering: " + str(len(newXyzPeakList)) + "\n")
        
        
        if successfulNoiseFiltering:
            featureObject.xyzPeakList = newXyzPeakList
        
        return featureObject
            

        
    def write_matches_and_fragment_matches_to_csv(self):

        # Function Description:
        #     Write the parent mass match and fragment mass match information
        #         into a CSV file
        # Function Input:
        #     All the match details for parent mass match and fragment mass match
        # Function Output:
        #     A CSV file containing those
        
        csvOutputFileName = self.csvOutputFileName
        
        if self.config.get("filePrefix"):
            filePrefix = self.config.get("filePrefix")
        else:
            filePrefix = ""
            
        csvOutputFileName = filePrefix + "CSV_1_" + csvOutputFileName
        
        baseOutput, extensionOutput = splitext(csvOutputFileName)
        csvOutputFileName = (baseOutput
                             + time.strftime("_%I_%M_%S_%p_%b_%d_%Y")
                             + extensionOutput) #time: %I:%M:%S_%p
        
        listOfRepetitiveFileNames = [fileName
                                     for fileName in glob.glob("*.csv")
                                     if fileName == csvOutputFileName]

        while listOfRepetitiveFileNames != []:
            baseOutput, extensionOutput = splitext(csvOutputFileName)
            csvOutputFileName = baseOutput + "_1" + extensionOutput
    
            listOfRepetitiveFileNames = [fileName
                                         for fileName in glob.glob("*.csv")
                                         if fileName == csvOutputFileName]
    
        outMatchesToCsv = open(csvOutputFileName + ".csv", 'w')

        outMatchesToCsv.write(self.specToPreappendToCsv)

        outMatchesToCsv.write("(1st ordering),,,,,,,,,,,,(2nd ordering)\n")
        
        outMatchesToCsv.write("Feature No.,mz,isotopeCorrection,z,Neutral Mass,"
                              + "RT,Inten.,Match No.,Protein,Peptide Theoretical Mass,Glycan Theoretical Mass,GP Theoretical Mass,ppm Error,"
                              + "Total Fragm. Score,Peptide Matches (Score 5),Backbone(score 1), Glycan Matches (Score 4),"
                              + "DeoxyHex Matches (Score 3),Neu5Ac Matches (Score 3),"
                              + "Neu5Gc Matches (Score 3),Common Matches (Score 1),"
                              + "Peptide,Glycosylation Site Position,is Decoy?,Peptide Mass Minus H2O,Hex,HexNAc,DeoxyHex,"
                              + "Neu5Ac,Neu5Gc,Pentose,KDN,HexA,xyz_abc counts,xyzBackboneFragments,abcBackboneFragments,peptide_plus_glycan_fragments"
                              + "Fragment Ion Match Details -- Fragment Name (Score Contribution): Fragment Mz (Fragment Neutral Mass)\n")

        for featureIndex in range(len(self.featureList)):
            massMatchList = self.featureList[featureIndex].massMatchList
            massMatchList = sorted(massMatchList,
                                   key = attrgetter("totalScore"),
                                   reverse = True)
            self.featureList[featureIndex].massMatchList = massMatchList

        isFeaturePrinted = False
        
        for feature in self.featureList:
            for massMatch in feature.massMatchList:
                if isFeaturePrinted == False:
                    outMatchesToCsv.write(  str(feature.featureIndex) + ","
                                          + str(feature.compoundMzValue) + ","
                                          + str(feature.isotopeCorrection) + ","
                                          + str(feature.compoundCharge) + ","
                                          + str(feature.compoundNeutralMass) + ","
                                          + str(feature.compoundRtCenter) + ","
                                          + str(feature.compoundIntensity) + ",")
                    isFeaturePrinted = True
                else:
                    outMatchesToCsv.write(str(feature.featureIndex) + ",,,,,,")

                xyz_abc = [0, 0]   # format = [xyz, abc]
                fragmentIonMatchDetails = ""
                xyzIonMatchDetails = ""
                abcIonMatchDetails = ""
                peptideIonMatchDetails = ""
                
                xyzRegExpression = re.compile(r"""[xyz]([0-9]{1,3})-""")
                abcRegExpression = re.compile(r"""[abc]([0-9]{1,3})-""")
                peptideRegExpression = re.compile(r"""peptide""")
                
                for fragmentNeutralMassMatch in massMatch.listOfFragmentNeutralMassMatches:
                    multiplicity = int(1.0 / fragmentNeutralMassMatch.matchedPeakTuple[3])
                    fragmentName = fragmentNeutralMassMatch.matchedFragmentName
                    if xyzRegExpression.search(fragmentName):
                        xyz_abc[0] += 1
                        xyzIonMatchDetails = (xyzIonMatchDetails
                                                   + fragmentNeutralMassMatch.matchedFragmentName
                                                   + " ("
                                                   + str(int(fragmentNeutralMassMatch.scoreContributed * multiplicity))
                                                   + "/"
                                                   + str(multiplicity)
                                                   + "): "
                                                   + str(fragmentNeutralMassMatch.matchedActualFragmentMz)
                                                   + " ("
                                                   + str(fragmentNeutralMassMatch.matchedActualFragmentNeutralMass)
                                                   + "); ")
                    elif abcRegExpression.search(fragmentName):
                        xyz_abc[1] += 1
                        abcIonMatchDetails = (abcIonMatchDetails
                                                   + fragmentNeutralMassMatch.matchedFragmentName
                                                   + " ("
                                                   + str(int(fragmentNeutralMassMatch.scoreContributed * multiplicity))
                                                   + "/"
                                                   + str(multiplicity)
                                                   + "): "
                                                   + str(fragmentNeutralMassMatch.matchedActualFragmentMz)
                                                   + " ("
                                                   + str(fragmentNeutralMassMatch.matchedActualFragmentNeutralMass)
                                                   + "); ")
                    elif peptideRegExpression.search(fragmentName):
                        peptideIonMatchDetails = (peptideIonMatchDetails
                                                   + fragmentNeutralMassMatch.matchedFragmentName
                                                   + " ("
                                                   + str(int(fragmentNeutralMassMatch.scoreContributed * multiplicity))
                                                   + "/"
                                                   + str(multiplicity)
                                                   + "): "
                                                   + str(fragmentNeutralMassMatch.matchedActualFragmentMz)
                                                   + " ("
                                                   + str(fragmentNeutralMassMatch.matchedActualFragmentNeutralMass)
                                                   + "); ")
                    else:
                        fragmentIonMatchDetails = (fragmentIonMatchDetails
                                                   + fragmentNeutralMassMatch.matchedFragmentName
                                                   + " ("
                                                   + str(int(fragmentNeutralMassMatch.scoreContributed * multiplicity))
                                                   + "/"
                                                   + str(multiplicity)
                                                   + "): "
                                                   + str(fragmentNeutralMassMatch.matchedActualFragmentMz)
                                                   + " ("
                                                   + str(fragmentNeutralMassMatch.matchedActualFragmentNeutralMass)
                                                   + "); ")

                fragmentIonMatchDetails = fragmentIonMatchDetails[:-2]
                xyzIonMatchDetails = xyzIonMatchDetails[:-2]
                abcIonMatchDetails = abcIonMatchDetails[:-2]
                peptideIonMatchDetails = peptideIonMatchDetails[:-2]
                xyz_abc = "_".join([str(num) for num in xyz_abc])
                
                if type(massMatch.peptideMass) != type(None):
                    peptideMass = massMatch.peptideMass + 18.0105646864
                else:
                    peptideMass = ""
                    
                isDecoyUsed = ""
                if massMatch.isDecoyUsed == True:
                    isDecoyUsed = "1"
                    ppmError = massMatch.deviationFromActualToDecoyedTheoreticalMass
                else:
                    isDecoyUsed = "0"
                    ppmError = massMatch.deviationFromActualToTheoreticalMass
                '''
                if "." in str(massMatch.totalScore):
                    if len(str(massMatch.totalScore).split(".")[-1]) == 1:
                        totalScore = str(massMatch.totalScore) + "0"
                    elif len(str(massMatch.totalScore).split(".")[-1]) == 2:
                        totalScore = str(massMatch.totalScore)
                else:
                    totalScore = str(massMatch.totalScore) + ".00"
                '''

                    
                outMatchesToCsv.write(str(massMatch.matchIndex) + ","
                                      + massMatch.proteinName + ","
                                      + str(peptideMass) + ","
                                      + str(massMatch.glycanMass) + ","
                                      + str(massMatch.totalTheoreticalMass) + ","
                                      + ppmError + ","
                                      + str(massMatch.totalScore) + ","
                                      + str(massMatch.numberOfPeptideIonFragmentMatches) + ","
                                      + str(massMatch.numberOfBackboneFragmentMatches) + ","
                                      + str(massMatch.numberOfGlycanIonFragmentMatches) + ","
                                      + str(massMatch.numberOfDeoxyHexIonFragmentMatches) + ","
                                      + str(massMatch.numberOfNeu5AcIonFragmentMatches) + ","
                                      + str(massMatch.numberOfNeu5GcIonFragmentMatches) + ","
                                      + str(massMatch.numberOfCommonIonFragmentMatches) + ","
                                      + str(massMatch.peptideStringMatched) + ","
                                      + str(massMatch.glycosylationSiteIndex) + ","
                                      + isDecoyUsed + ","
                                      + str(massMatch.peptideMass) + ","
                                      + str(massMatch.glycanCompositionMatched["Hex"]) + ","
                                      + str(massMatch.glycanCompositionMatched["HexNAc"]) + ","
                                      + str(massMatch.glycanCompositionMatched["DeoxyHex"]) + ","
                                      + str(massMatch.glycanCompositionMatched["Neu5Ac"]) + ","
                                      + str(massMatch.glycanCompositionMatched["Neu5Gc"]) + ","
                                      + str(massMatch.glycanCompositionMatched["Pentose"]) + ","
                                      + str(massMatch.glycanCompositionMatched["KDN"]) + ","
                                      + str(massMatch.glycanCompositionMatched["HexA"]) + ","
                                      + xyz_abc + ","
                                      + xyzIonMatchDetails + ","
                                      + abcIonMatchDetails + ","
                                      + peptideIonMatchDetails + ","
                                      + fragmentIonMatchDetails + "\n")
                                            
            isFeaturePrinted = False                           
                
        outMatchesToCsv.close()

        print("Matches and Fragment Matches writing to CSV done.")
        print("File name: " + csvOutputFileName + "\n")
        


    def write_matches_and_fragment_matches_to_csv_based_on_glycan_site_order(self):

        # Function Description:
        #     Write the parent mass match and fragment mass match information
        #         ordered by the glycosylation site position
        #         into a CSV file
        # Function Input:
        #     All the match details for parent mass match and fragment mass match
        # Function Output:
        #     A CSV file containing those
        
        csvOutputFileName = self.csvOutputFileName

        if self.config.get("filePrefix"):
            filePrefix = self.config.get("filePrefix")
        else:
            filePrefix = ""
            
        csvOutputFileName = filePrefix + "CSV_2_site_specific_" + self.csvOutputFileName
        
        baseOutput, extensionOutput = splitext(csvOutputFileName)
        csvOutputFileName = (baseOutput
                             + time.strftime("_%I_%M_%S_%p_%b_%d_%Y")
                             + extensionOutput) #time: %I:%M:%S_%p
        
        listOfRepetitiveFileNames = [fileName
                                     for fileName in glob.glob("*.csv")
                                     if fileName == csvOutputFileName]

        while listOfRepetitiveFileNames != []:
            baseOutput, extensionOutput = splitext(csvOutputFileName)
            csvOutputFileName = baseOutput + "_1" + extensionOutput
    
            listOfRepetitiveFileNames = [fileName
                                         for fileName in glob.glob("*.csv")
                                         if fileName == csvOutputFileName]
    
        outMatchesToCsv = open(csvOutputFileName + ".csv", 'w')

        outMatchesToCsv.write(self.specToPreappendToCsv)

        outMatchesToCsv.write("(3rd ordering),,,,,,,,,(1st ordering),(2nd ordering)\n")
        outMatchesToCsv.write("Feature No.,Mass Match No.,Protein,Peptide Theoretical Mass,Glycan Theoretical Mass,GP Theoretical Mass,GP Actual Mass,ppm Error,GP Compound MZ,Glycosylation Site Position,"
                              + "Mass Match Score,RT,Nonstandard Glycan Components,Peptide Fragment Matches (Score 5),Peptide Sequence,"
                              + "Hex,HexNAc,DeoxyHex,Sialic Acid (Neu5Ac and Neu5Gc),is Decoy?,Mass Match Intensity\n")
        
        totalMassMatchList = [massMatch
                              for feature in self.featureList
                              for massMatch in feature.massMatchList]

        totalMassMatchList = sorted(totalMassMatchList,
                                    key = attrgetter("featureIndex"),
                                    reverse = False)

        totalMassMatchList = sorted(totalMassMatchList,
                                    key = attrgetter("totalScore"),
                                    reverse = True)
        
        totalMassMatchList = sorted(totalMassMatchList,
                                    key = attrgetter("glycosylationSiteIndex"),
                                    reverse = False)
        
        for massMatch in totalMassMatchList:
            nonstandardGlycanComponents = ""
            
            for glycanComponent in massMatch.glycanCompositionMatched.keys():
                if (glycanComponent == "Pentose" or
                        glycanComponent == "KDN" or
                        glycanComponent == "HexA"):
                    if massMatch.glycanCompositionMatched[glycanComponent] != 0:
                        nonstandardGlycanComponents = (nonstandardGlycanComponents
                                + glycanComponent + ": "
                                + str(massMatch.glycanCompositionMatched[glycanComponent]) + "; ")

            nonstandardGlycanComponents = nonstandardGlycanComponents[:-2]
            '''
            if "." in str(massMatch.totalScore):
                if len(str(massMatch.totalScore).split(".")[-1]) == 1:
                    totalScore = str(massMatch.totalScore) + "0"
                elif len(str(massMatch.totalScore).split(".")[-1]) == 2:
                    totalScore = str(massMatch.totalScore)
            else:
                    totalScore = str(massMatch.totalScore) + ".00"
            '''

            if type(massMatch.peptideMass) != type(None):
                peptideMass = massMatch.peptideMass + 18.0105646864
            else:
                peptideMass = ""
                    
            if nonstandardGlycanComponents == "":
                nonstandardGlycanComponents = "not any"
                
            isDecoyUsed = ""
            if massMatch.isDecoyUsed == True:
                isDecoyUsed = "1"
                ppmError = massMatch.deviationFromActualToDecoyedTheoreticalMass
            else:
                isDecoyUsed = "0"
                ppmError = massMatch.deviationFromActualToTheoreticalMass
                    
            outMatchesToCsv.write(str(massMatch.featureIndex) + ","
                                  + str(massMatch.matchIndex) + ","
                                  + massMatch.proteinName + ","
                                  + str(peptideMass) + ","
                                  + str(massMatch.glycanMass) + ","
                                  + str(massMatch.totalTheoreticalMass) + ","
                                  + str(massMatch.inputMassToMatch) + ","
                                  + ppmError + ","
                                  + str(massMatch.compoundMzValue) + ","
                                  + str(massMatch.glycosylationSiteIndex) + ","
                                  + str(massMatch.totalScore) + ","
                                  + str(massMatch.compoundRtCenter) + ","
                                  + nonstandardGlycanComponents + ","
                                  + str(massMatch.numberOfPeptideIonFragmentMatches) + ","
                                  + str(massMatch.peptideStringMatched) + ","
                                  + str(massMatch.glycanCompositionMatched["Hex"]) + ","
                                  + str(massMatch.glycanCompositionMatched["HexNAc"]) + ","
                                  + str(massMatch.glycanCompositionMatched["DeoxyHex"]) + ","
                                  + str(massMatch.glycanCompositionMatched["Neu5Ac"]
                                        + massMatch.glycanCompositionMatched["Neu5Gc"]) + ","
                                  + isDecoyUsed + ","
                                  + str(massMatch.compoundIntensity) + "\n")

        print("Matches and Fragment Matches writing to CSV done.")
        print("File name: " + csvOutputFileName + "\n")



    def write_highest_score_only_matches_and_fragment_matches_to_csv_based_on_glycan_site_order(self, doCSV5=False):

        # Function Description:
        #     Write the parent mass match and fragment mass match information
        #         ordered by the glycosylation site position
        #         into a CSV file
        # Function Input:
        #     All the match details for parent mass match and fragment mass match
        # Function Output:
        #     A CSV file containing those
        
        csvOutputFileName = self.csvOutputFileName
        
        if self.config.get("filePrefix"):
            filePrefix = self.config.get("filePrefix")
        else:
            filePrefix = ""
        
        if doCSV5:
            csvOutputFileName = filePrefix + "CSV_5_FDR_passing_only_" + self.csvOutputFileName
        else:
            csvOutputFileName = filePrefix + "CSV_3_highest_score_only_site_specific_" + self.csvOutputFileName
        
        baseOutput, extensionOutput = splitext(csvOutputFileName)
        csvOutputFileName = (baseOutput
                             + time.strftime("_%I_%M_%S_%p_%b_%d_%Y")
                             + extensionOutput) #time: %I:%M:%S_%p
        
        listOfRepetitiveFileNames = [fileName
                                     for fileName in glob.glob("*.csv")
                                     if fileName == csvOutputFileName]

        while listOfRepetitiveFileNames != []:
            baseOutput, extensionOutput = splitext(csvOutputFileName)
            csvOutputFileName = baseOutput + "_1" + extensionOutput
    
            listOfRepetitiveFileNames = [fileName
                                         for fileName in glob.glob("*.csv")
                                         if fileName == csvOutputFileName]
    
        outMatchesToCsv = open(csvOutputFileName + ".csv", 'w')

        outMatchesToCsv.write(self.specToPreappendToCsv)
        
        if doCSV5:
            scoreCutoff = self.scoreCutoffForFdrAnalysis
            #rescore parameters
            rescoreParams = ("pre Rescore passes," + str(self.passingFeaturesBeforeRescore)
                             + "\npost Rescore passes," + str(self.passingFeaturesAfterRescore)
                             + "\nfirst threshold," + str(self.scoreCutoffForFdrAnalysisFirstPass)
                             + "\nscore threshold," + str(scoreCutoff)
                             + "\n\n")
            outMatchesToCsv.write(rescoreParams)

        outMatchesToCsv.write("(3rd ordering),,,,,,,,,(1st ordering),(2nd ordering)\n")
        
        outMatchesToCsv.write("Feature No.,Mass Match No.,Protein,Peptide Theoretical Mass,Glycan Theoretical Mass,GP Theoretical Mass,isotopeError,GP Actual Mass,ppm Error,GP Compound MZ,Glycosylation Site Position,"
                              + "Mass Match Score,Highest Score - 2nd Highest Score in Feature,Number of Score Occurrences,RT,Nonstandard Glycan Components,Peptide Matches (Score 5), Backbone Matches,Peptide Sequence,"
                              + "Hex,HexNAc,DeoxyHex,Sialic Acid (Neu5Ac and Neu5Gc),is Decoy?,Mass Match Intensity\n")

        totalMassMatchList = []
        ''' previous: to print all mass matches
        totalMassMatchList = [massMatch
                              for feature in self.featureList
                              for massMatch in feature.massMatchList]

        totalMassMatchList = sorted(totalMassMatchList,
                                    key = attrgetter("featureIndex"),
                                    reverse = False)
        '''

        ''' current: to only print the highest scored mass match within every feature '''
        passingPeptides = 0
        for feature in self.featureList:
            if feature.massMatchList != []:
                massMatchWithHighestScoreWithinOneFeature = max(feature.massMatchList,
                                                                key = attrgetter("totalScore"))
                
                highestScoreWithinOneFeature = massMatchWithHighestScoreWithinOneFeature.totalScore
                numberOfHighestScoreOccurrences = len([massMatch
                                                       for massMatch in feature.massMatchList
                                                       if massMatch.totalScore == highestScoreWithinOneFeature])
                massMatchWithHighestScoreWithinOneFeature.numberOfHighestScoreOccurrences = numberOfHighestScoreOccurrences
                totalMassMatchList.append(massMatchWithHighestScoreWithinOneFeature)
                
        totalMassMatchList = sorted(totalMassMatchList,
                                    key = attrgetter("totalScore"),
                                    reverse = True)
        
        totalMassMatchList = sorted(totalMassMatchList,
                                    key = attrgetter("glycosylationSiteIndex"),
                                    reverse = False)
        
        for massMatch in totalMassMatchList:
            nonstandardGlycanComponents = ""
            
            for glycanComponent in massMatch.glycanCompositionMatched.keys():
                if (glycanComponent == "Pentose" or
                        glycanComponent == "KDN" or
                        glycanComponent == "HexA"):
                    if massMatch.glycanCompositionMatched[glycanComponent] != 0:
                        nonstandardGlycanComponents = (nonstandardGlycanComponents
                                + glycanComponent + ": "
                                + str(massMatch.glycanCompositionMatched[glycanComponent]) + "; ")

            nonstandardGlycanComponents = nonstandardGlycanComponents[:-2]
            '''
            if "." in str(massMatch.totalScore):
                if len(str(massMatch.totalScore).split(".")[-1]) == 1:
                    totalScore = str(massMatch.totalScore) + "0"
                elif len(str(massMatch.totalScore).split(".")[-1]) == 2:
                    totalScore = str(massMatch.totalScore)
            else:
                    totalScore = str(massMatch.totalScore) + ".00"
            '''
            if type(massMatch.peptideMass) != type(None):
                peptideMass = massMatch.peptideMass + 18.0105646864
            else:
                peptideMass = ""
                    
            if nonstandardGlycanComponents == "":
                nonstandardGlycanComponents = "not any"
                
            isDecoyUsed = ""
            if massMatch.isDecoyUsed == True:
                isDecoyUsed = "1"
                ppmError = massMatch.deviationFromActualToDecoyedTheoreticalMass
            else:
                isDecoyUsed = "0"
                ppmError = massMatch.deviationFromActualToTheoreticalMass
                
            featureObject = [featureObject
                             for featureObject in self.featureList
                             if massMatch in featureObject.massMatchList][0]
            if not doCSV5:
                outMatchesToCsv.write(str(massMatch.featureIndex) + ","
                                      + str(massMatch.matchIndex) + ","
                                      + str(massMatch.proteinName) + ","
                                      + str(peptideMass) + ","
                                      + str(massMatch.glycanMass) + ","
                                      + str(massMatch.totalTheoreticalMass) + ","
                                      + str(featureObject.isotopeCorrection) + ","
                                      + str(massMatch.inputMassToMatch) + ","
                                      + ppmError + ","
                                      + str(massMatch.compoundMzValue) + ","
                                      + str(massMatch.glycosylationSiteIndex) + ","
                                      + str(massMatch.totalScore) + ","
                                      + str(featureObject.differenceBetweenHighestAnd2ndHighestScore) + ","
                                      + str(massMatch.numberOfHighestScoreOccurrences) + ","
                                      + str(massMatch.compoundRtCenter) + ","
                                      + nonstandardGlycanComponents + ","
                                      + str(massMatch.numberOfPeptideIonFragmentMatches) + ","
                                      + str(massMatch.peptideStringMatched) + ","
                                      + str(massMatch.glycanCompositionMatched["Hex"]) + ","
                                      + str(massMatch.glycanCompositionMatched["HexNAc"]) + ","
                                      + str(massMatch.glycanCompositionMatched["DeoxyHex"]) + ","
                                      + str(massMatch.glycanCompositionMatched["Neu5Ac"]
                                            + massMatch.glycanCompositionMatched["Neu5Gc"]) + ","
                                      + isDecoyUsed + ","
                                      + str(massMatch.compoundIntensity) + "\n")
            else:
                if massMatch.totalReScore >= scoreCutoff:
                    passingPeptides += 1
                    outMatchesToCsv.write(str(massMatch.featureIndex) + ","
                                      + str(massMatch.matchIndex) + ","
                                      + str(massMatch.proteinName) + ","
                                      + str(peptideMass) + ","
                                      + str(massMatch.glycanMass) + ","
                                      + str(massMatch.totalTheoreticalMass) + ","
                                      + str(featureObject.isotopeCorrection) + ","
                                      + str(massMatch.inputMassToMatch) + ","
                                      + ppmError + ","
                                      + str(massMatch.compoundMzValue) + ","
                                      + str(massMatch.glycosylationSiteIndex) + ","
                                      + str(massMatch.totalReScore) + ","
                                      + str(featureObject.differenceBetweenHighestAnd2ndHighestScore) + ","
                                      + str(massMatch.numberOfHighestScoreOccurrences) + ","
                                      + str(massMatch.compoundRtCenter) + ","
                                      + nonstandardGlycanComponents + ","
                                      + str(massMatch.numberOfPeptideIonFragmentMatches) + ","
                                      + str(massMatch.numberOfBackboneFragmentMatches) + ","
                                      + str(massMatch.peptideStringMatched) + ","
                                      + str(massMatch.glycanCompositionMatched["Hex"]) + ","
                                      + str(massMatch.glycanCompositionMatched["HexNAc"]) + ","
                                      + str(massMatch.glycanCompositionMatched["DeoxyHex"]) + ","
                                      + str(massMatch.glycanCompositionMatched["Neu5Ac"]
                                      + massMatch.glycanCompositionMatched["Neu5Gc"]) + ","
                                      + isDecoyUsed + ","
                                      + str(massMatch.compoundIntensity) + "\n")
        print("Matches and Fragment Matches writing to CSV done.")
        print("File name: " + csvOutputFileName + "\n")
        if doCSV5:
            return (scoreCutoff, self.scoreCutoffForFdrAnalysisFirstPass, passingPeptides)



    def write_all_highest_score_only_matches_and_fragment_matches_to_csv_based_on_glycan_site_order(self):

        # Function Description:
        #     Write the parent mass match and fragment mass match information
        #         ordered by the glycosylation site position
        #         into a CSV file
        # Function Input:
        #     All the match details for parent mass match and fragment mass match
        # Function Output:
        #     A CSV file containing those
        
        csvOutputFileName = self.csvOutputFileName
        
        csvOutputFileName = "CSV_4_all_highest_score_only_site_specific_" + self.csvOutputFileName
        
        baseOutput, extensionOutput = splitext(csvOutputFileName)
        csvOutputFileName = (baseOutput
                             + time.strftime("_%I_%M_%S_%p_%b_%d_%Y")
                             + extensionOutput) #time: %I:%M:%S_%p
        
        listOfRepetitiveFileNames = [fileName
                                     for fileName in glob.glob("*.csv")
                                     if fileName == csvOutputFileName]

        while listOfRepetitiveFileNames != []:
            baseOutput, extensionOutput = splitext(csvOutputFileName)
            csvOutputFileName = baseOutput + "_1" + extensionOutput
    
            listOfRepetitiveFileNames = [fileName
                                         for fileName in glob.glob("*.csv")
                                         if fileName == csvOutputFileName]
    
        outMatchesToCsv = open(csvOutputFileName + ".csv", 'w')

        outMatchesToCsv.write(self.specToPreappendToCsv)

        outMatchesToCsv.write("(3rd ordering),,,,,,,,,(1st ordering),(2nd ordering)\n")
        
        outMatchesToCsv.write("Feature No.,Mass Match No.,Protein,Peptide Theoretical Mass,Glycan Theoretical Mass,GP Theoretical Mass,GP Actual Mass,IsotopeError,ppm Error,GP Compound MZ,Glycosylation Site Position,"
                              + "Mass Match Score,Highest Score - 2nd Highest Score in Feature,Number of Score Occurrences,RT,Nonstandard Glycan Components,Peptide Matches (Score 5),Peptide Sequence,"
                              + "Hex,HexNAc,DeoxyHex,Sialic Acid (Neu5Ac and Neu5Gc),is Decoy?,Mass Match Intensity\n")

        totalMassMatchList = []
        ''' previous: to print all mass matches
        totalMassMatchList = [massMatch
                              for feature in self.featureList
                              for massMatch in feature.massMatchList]

        totalMassMatchList = sorted(totalMassMatchList,
                                    key = attrgetter("featureIndex"),
                                    reverse = False)
        '''

        ''' current: to only print the highest scored mass match within every feature '''
        
        for feature in self.featureList:
            if feature.massMatchList != []:
                massMatchWithHighestScoreWithinOneFeature = max(feature.massMatchList,
                                                                key = attrgetter("totalScore"))
                
                highestScoreWithinOneFeature = massMatchWithHighestScoreWithinOneFeature.totalScore
                massMatchesWithHighestScoreWithinOneFeature = [massMatch
                                                               for massMatch in feature.massMatchList
                                                               if massMatch.totalScore == highestScoreWithinOneFeature]
                
                numberOfHighestScoreOccurrences = len(massMatchesWithHighestScoreWithinOneFeature)

                for massMatch in massMatchesWithHighestScoreWithinOneFeature:
                    massMatch.numberOfHighestScoreOccurrences = numberOfHighestScoreOccurrences
                
                totalMassMatchList.extend(massMatchesWithHighestScoreWithinOneFeature)
                
        totalMassMatchList = sorted(totalMassMatchList,
                                    key = attrgetter("totalScore"),
                                    reverse = True)
        
        totalMassMatchList = sorted(totalMassMatchList,
                                    key = attrgetter("glycosylationSiteIndex"),
                                    reverse = False)
        
        for massMatch in totalMassMatchList:
            nonstandardGlycanComponents = ""
            
            for glycanComponent in massMatch.glycanCompositionMatched.keys():
                if (glycanComponent == "Pentose" or
                        glycanComponent == "KDN" or
                        glycanComponent == "HexA"):
                    if massMatch.glycanCompositionMatched[glycanComponent] != 0:
                        nonstandardGlycanComponents = (nonstandardGlycanComponents
                                + glycanComponent + ": "
                                + str(massMatch.glycanCompositionMatched[glycanComponent]) + "; ")

            nonstandardGlycanComponents = nonstandardGlycanComponents[:-2]
            '''
            if "." in str(massMatch.totalScore):
                if len(str(massMatch.totalScore).split(".")[-1]) == 1:
                    totalScore = str(massMatch.totalScore) + "0"
                elif len(str(massMatch.totalScore).split(".")[-1]) == 2:
                    totalScore = str(massMatch.totalScore)
            else:
                    totalScore = str(massMatch.totalScore) + ".00"
            '''
            if type(massMatch.peptideMass) != type(None):
                peptideMass = massMatch.peptideMass + 18.0105646864
            else:
                peptideMass = ""
                    
            if nonstandardGlycanComponents == "":
                nonstandardGlycanComponents = "not any"
                
            isDecoyUsed = ""
            if massMatch.isDecoyUsed == True:
                isDecoyUsed = "1"
                ppmError = massMatch.deviationFromActualToDecoyedTheoreticalMass
            else:
                isDecoyUsed = "0"
                ppmError = massMatch.deviationFromActualToTheoreticalMass
                
            featureObject = [featureObject
                             for featureObject in self.featureList
                             if massMatch in featureObject.massMatchList][0]

            outMatchesToCsv.write(str(massMatch.featureIndex) + ","
                                  + str(massMatch.matchIndex) + ","
                                  + str(massMatch.proteinName) + ","
                                  + str(peptideMass) + ","
                                  + str(massMatch.glycanMass) + ","
                                  + str(massMatch.totalTheoreticalMass) + ","
                                  + str(massMatch.inputMassToMatch) + ","
                                  + str(featureObject.isotopeCorrection) + ","
                                  + ppmError + ","
                                  + str(massMatch.compoundMzValue) + ","
                                  + str(massMatch.glycosylationSiteIndex) + ","
                                  + str(massMatch.totalScore) + ","
                                  + str(featureObject.differenceBetweenHighestAnd2ndHighestScore) + ","
                                  + str(massMatch.numberOfHighestScoreOccurrences) + ","
                                  + str(massMatch.compoundRtCenter) + ","
                                  + nonstandardGlycanComponents + ","
                                  + str(massMatch.numberOfPeptideIonFragmentMatches) + ","
                                  + str(massMatch.peptideStringMatched) + ","
                                  + str(massMatch.glycanCompositionMatched["Hex"]) + ","
                                  + str(massMatch.glycanCompositionMatched["HexNAc"]) + ","
                                  + str(massMatch.glycanCompositionMatched["DeoxyHex"]) + ","
                                  + str(massMatch.glycanCompositionMatched["Neu5Ac"]
                                        + massMatch.glycanCompositionMatched["Neu5Gc"]) + ","
                                  + isDecoyUsed + ","
                                  + str(massMatch.compoundIntensity) + "\n")

        print("Matches and Fragment Matches writing to CSV done.")
        print("File name: " + csvOutputFileName + "\n")

    
    def rescore_based_on_crpytic_proteolytic_specificity_andor_ppm_distr(self, scoreIncreasePPC=30, pValue=0.05, useCrypticSpecificity=True, useErrorDistribution=False):
        
        #
        #
        #  the following code generates a plot of score occurrence verses count for both
        #     normal and decoy results. Less overlap = better
        #
        #
        
        #generate dataset for score distribution figure
        massMatchIntScores = [i for i in range(200)]
        massMatchCount = [0 for i in range(200)]
        decoyMatchIntScores = [i for i in range(200)]
        decoyMatchCount = [0 for i in range(200)]
        
        #empty lists for score vs. ppm error plot
        ppmVsScoreScoreList = []
        ppmVsScorePPMList = []
        ppmVsScoreScoreListDecoy = []
        ppmVsScorePPMListDecoy = []
        
        # Loop used to accumulate all values for lists above
        for feature in self.featureList:
            for massMatch in feature.massMatchList:
                totalScore = massMatch.totalScore
                intScore = int(totalScore)
                
                #discrete score distributions
                if massMatch.isDecoyUsed == 0:
                    if intScore in massMatchIntScores:
                        index = massMatchIntScores.index(intScore)
                        massMatchCount[index] += 1
                else:
                    if intScore in decoyMatchIntScores:
                        index = decoyMatchIntScores.index(intScore)
                        decoyMatchCount[index] += 1
                
                #get ppm error and store in lists
                if massMatch.isDecoyUsed == True:
                    ppmVsScoreScoreListDecoy.append(totalScore)
                    ppmErr = massMatch.deviationFromActualToDecoyedTheoreticalMass.split(' ')[0]
                    ppmErr = float(ppmErr)
                    ppmVsScorePPMListDecoy.append(ppmErr)
                else:
                    ppmVsScoreScoreList.append(totalScore)
                    ppmErr = massMatch.deviationFromActualToTheoreticalMass.split(' ')[0]
                    ppmErr = float(ppmErr)
                    ppmVsScorePPMList.append(ppmErr)


        
        occupiedIndexList = [i for i,count in enumerate(massMatchCount) if count > 0]
        occupiedIndexList.extend([i for i,count in enumerate(decoyMatchCount) if count > 0])
        maxIndexForPlotting = min(max(occupiedIndexList) + 5, 199)
        
        import matplotlib.pyplot as plt
        from scipy.interpolate import interp1d
        MIP = maxIndexForPlotting
        
        #perform polynomial interpolation
        xnew = np.linspace(0,MIP,300)
        massMatch_smooth_fx = interp1d(massMatchIntScores[0:MIP],massMatchCount[0:MIP], kind='cubic', bounds_error=False, fill_value=0)
        massMatch_smooth = [max(0,y) for y in massMatch_smooth_fx(xnew)]
        decoyMatch_smooth_fx = interp1d(decoyMatchIntScores[0:MIP],decoyMatchCount[0:MIP], kind='cubic', bounds_error=False, fill_value=0)
        decoyMatch_smooth = [max(0,y) for y in decoyMatch_smooth_fx(xnew)]
        
        # plot ppm error vs score boxes
        if self.config.get("dontMakePlots"):
            pass
        else:
            plt.plot(ppmVsScoreScoreListDecoy , ppmVsScorePPMListDecoy, 'sr')
            plt.plot(ppmVsScoreScoreList , ppmVsScorePPMList, 'sb')
            plt.title('PPM Error vs. Score')
            plt.xlabel('Fragmentation Score')
            plt.ylabel('PPM Error')
            plt.show()

        #
        #  BEGIN Cryptic specificity
        #
        #find cutoff score for p= pValue
        #
        if useCrypticSpecificity:
            cutoffScore = 401
            for score in range(0,400):
                areaDecoy = sum(decoyMatchCount[score:])
                areaReal = sum(massMatchCount[score:])
                totalArea = areaReal+areaDecoy
                if totalArea == 0:
                    decoyRate = 1
                else:
                    decoyRate = (2.0*areaDecoy)/float(areaReal+areaDecoy)
                if decoyRate <= pValue:
                    cutoffScore = score
                    break
            if cutoffScore == 401:
                raise NameError("Dataset and configuration does not support p=%f rescore cutoff"%(pValue))
            if self.config.get("dontMakePlots"):
                pass
            else:
                plt.plot(xnew,decoyMatch_smooth, 'r-')
                plt.plot(xnew,massMatch_smooth, 'b-')
                plt.plot(massMatchIntScores[0:MIP] , massMatchCount[0:MIP], 'ob')
                plt.plot(decoyMatchIntScores[0:MIP] , decoyMatchCount[0:MIP], 'xr')
                plt.plot([cutoffScore, cutoffScore], [0, 10], '-k', lw=2.8)
                plt.title('Decoy and Library Score Intensities')
                plt.xlabel('Total Score')
                plt.ylabel('Score Intensity')
                plt.show()
            
            #export cutoff score for later output
            self.scoreCutoffForFdrAnalysisFirstPass = cutoffScore

            #find list of peptides above score threshold
            peptideDict = {}
            for feature in self.featureList:
                for massMatch in feature.massMatchList:
                    if massMatch.totalScore >= cutoffScore:
                        self.passingFeaturesBeforeRescore += 1
                        peptide = massMatch.peptideStringMatched
                        if peptide in peptideDict:
                            peptideDict[peptide] += 1
                        else:
                            peptideDict[peptide] = 1
            peptideCountList = [(count, peptide) for (peptide, count) in peptideDict.iteritems()]
            peptideCountList.sort()
            countList, peptideList = zip(*peptideCountList)
            meanCount = sum(countList)/len(countList)    
            maxCount = max(countList)
            
            maxIncrease = 1.01 + (scoreIncreasePPC/100.0)
            middleIncrease = (maxIncrease - 1.01)*0.5 + 1.01
            adjustment = [1.01, middleIncrease, maxIncrease]
            countAxis = [0.99, meanCount, maxCount+0.1]
            scoreAdjustmentFunction = interp1d(countAxis, adjustment, kind='linear', bounds_error=False, fill_value=1)
            
            print(peptideDict)
            
            #perform re-scoring
            for feature in self.featureList:
                for massMatch in feature.massMatchList:
                    peptideMatched = massMatch.peptideStringMatched
                    if peptideMatched in peptideDict:
                        count = peptideDict[peptideMatched]
                        scoreAdjustmentFactor = float(scoreAdjustmentFunction(count))
                        print("count, adjustF = %i, %f"%(count, scoreAdjustmentFactor))
                        massMatch.totalReScore = massMatch.totalScore*scoreAdjustmentFactor
                    else:
                        massMatch.totalReScore = massMatch.totalScore
                        
            
                        
        if useErrorDistribution:
            #make error distribution
            #adjust scores
            pass
                    
                    
    #generate dataset for score distribution figure
        massMatchIntScores = [i for i in range(200)]
        massMatchCount = [0 for i in range(200)]
        decoyMatchIntScores = [i for i in range(200)]
        decoyMatchCount = [0 for i in range(200)]
        
        #empty lists for score vs. ppm error plot
        ppmVsScoreScoreList = []
        ppmVsScorePPMList = []
        ppmVsScoreScoreListDecoy = []
        ppmVsScorePPMListDecoy = []
        
        # Loop used to accumulate all values for lists above
        for feature in self.featureList:
            for massMatch in feature.massMatchList:
                totalScore = massMatch.totalReScore
                intScore = int(totalScore)
                
                #discrete score distributions
                if massMatch.isDecoyUsed == 0:
                    if intScore in massMatchIntScores:
                        index = massMatchIntScores.index(intScore)
                        massMatchCount[index] += 1
                else:
                    if intScore in decoyMatchIntScores:
                        index = decoyMatchIntScores.index(intScore)
                        decoyMatchCount[index] += 1
                
                #get ppm error and store in lists
                if massMatch.isDecoyUsed == True:
                    ppmVsScoreScoreListDecoy.append(totalScore)
                    ppmErr = massMatch.deviationFromActualToDecoyedTheoreticalMass.split(' ')[0]
                    ppmErr = float(ppmErr)
                    ppmVsScorePPMListDecoy.append(ppmErr)
                else:
                    ppmVsScoreScoreList.append(totalScore)
                    ppmErr = massMatch.deviationFromActualToTheoreticalMass.split(' ')[0]
                    ppmErr = float(ppmErr)
                    ppmVsScorePPMList.append(ppmErr)


        occupiedIndexList = [1]
        occupiedIndexList.extend([i for i,count in enumerate(massMatchCount) if count > 0])
        occupiedIndexList.extend([i for i,count in enumerate(decoyMatchCount) if count > 0])
        maxIndexForPlotting = min(max(occupiedIndexList) + 5, 199)
        
        import matplotlib.pyplot as plt
        from scipy.interpolate import interp1d
        MIP = maxIndexForPlotting
        
        #perform polynomial interpolation
        xnew = np.linspace(0,MIP,300)
        massMatch_smooth_fx = interp1d(massMatchIntScores[0:MIP],massMatchCount[0:MIP], kind='linear', bounds_error=False, fill_value=0)
        massMatch_smooth = [max(0,y) for y in massMatch_smooth_fx(xnew)]
        decoyMatch_smooth_fx = interp1d(decoyMatchIntScores[0:MIP],decoyMatchCount[0:MIP], kind='linear', bounds_error=False, fill_value=0)
        decoyMatch_smooth = [max(0,y) for y in decoyMatch_smooth_fx(xnew)]
        
        # plot ppm error vs score boxes
        if self.config.get("dontMakePlots"):
            pass
        else:
            plt.plot(ppmVsScoreScoreListDecoy , ppmVsScorePPMListDecoy, 'sr')
            plt.plot(ppmVsScoreScoreList , ppmVsScorePPMList, 'sb')
            plt.title('PPM Error vs. Score; pass 2')
            plt.xlabel('Fragmentation Score')
            plt.ylabel('PPM Error')
            plt.show()
        
        #report on updated scoring parameters and export 
        if useCrypticSpecificity:
            cutoffScore = 201
            for score in range(0,100):
                areaDecoy = sum(decoyMatchCount[score:])
                areaReal = sum(massMatchCount[score:])
                totalArea = areaReal+areaDecoy
                if totalArea == 0:
                    decoyRate = 1
                else:
                    decoyRate = (2.0*areaDecoy)/float(totalArea)
                if decoyRate <= pValue:
                    cutoffScore = score
                    break
            if cutoffScore == 202:
                raise NameError("Rescored dataset does not support p=%f rescore cutoff"%(pValue))
            if self.config.get("dontMakePlots"):
                pass
            else:
                plt.plot(xnew,decoyMatch_smooth, 'r-')
                plt.plot(xnew,massMatch_smooth, 'b-')
                plt.plot([cutoffScore, cutoffScore], [0, 10], '-k', lw=2.8)
                plt.plot(massMatchIntScores[0:MIP] , massMatchCount[0:MIP], 'ob')
                plt.plot(decoyMatchIntScores[0:MIP] , decoyMatchCount[0:MIP], 'xr')
                plt.title('Decoy and Library Score Intensities; pass 2')
                plt.xlabel('Total Score')
                plt.ylabel('Score Intensity')
                plt.show()
            
            peptideDict = {}
            for feature in self.featureList:
                for massMatch in feature.massMatchList:
                    if massMatch.totalReScore >= cutoffScore:
                        self.passingFeaturesAfterRescore += 1
                        peptide = massMatch.peptideStringMatched
                        if peptide in peptideDict:
                            peptideDict[peptide] += 1
                        else:
                            peptideDict[peptide] = 1
            
            self.scoreCutoffForFdrAnalysis = cutoffScore
            print(peptideDict)
        
    
    def write_matches_and_fragment_matches_to_xml(self):

        xmlOutputFileName = self.csvOutputFileName

        xmlOutputFileName = "XML_1_" + xmlOutputFileName

        baseOutput, extensionOutput = splitext(xmlOutputFileName)
        xmlOutputFileName = (baseOutput
                             + time.strftime("_%I_%M_%S_%p_%b_%d_%Y")
                             + extensionOutput) #time: %I:%M:%S_%p
        
        listOfRepetitiveFileNames = [fileName
                                     for fileName in glob.glob("*.csv")
                                     if fileName == xmlOutputFileName]

        while listOfRepetitiveFileNames != []:
            baseOutput, extensionOutput = splitext(xmlOutputFileName)
            xmlOutputFileName = baseOutput + "_1" + extensionOutput
    
            listOfRepetitiveFileNames = [fileName
                                         for fileName in glob.glob("*.csv")
                                         if fileName == xmlOutputFileName]

        outMatchesToXml = open(xmlOutputFileName + ".xml", 'w')
        
        doc = minidom.Document()

        xmlNode = doc.createElement("xml")
        doc.appendChild(xmlNode)
        
        gpFinderConfigNode = doc.createElement("gpfinderconfig")
        xmlNode.appendChild(gpFinderConfigNode)
        
        if self.config["isSequenceGivenByUser"] == True:
            proteinSequenceNode = doc.createElement("protein")
            proteinSequence = self.config["sequence"]
            proteinSequenceNode.setAttribute("sequence", proteinSequence)
            gpFinderConfigNode.appendChild(proteinSequenceNode)
        else:
            proteinFileNode = doc.createElement("protein")
            proteinFileName = self.config["proteinFile"]
            proteinFileNode.setAttribute("filename", proteinFileName)
            gpFinderConfigNode.appendChild(proteinFileNode)

        if self.config["msFormat"] == "cef":
            cefFileNode = doc.createElement("msFile")
            cefFileName = self.config["msFile"]
            cefFileNode.setAttribute("filename", cefFileName)
            gpFinderConfigNode.appendChild(cefFileNode)
        elif self.config["msFormat"] == "mgf":
            mgfFileNode = doc.createElement("msFile")
            mgfFileName = self.config["msFile"]
            mgfFileNode.setAttribute("filename", mgfFileName)
            gpFinderConfigNode.appendChild(mgfFileNode)

        glycanNode = doc.createElement("glycan")
        gpFinderConfigNode.appendChild(glycanNode)
        
        glycanSearchMethodNode = doc.createElement("method")
        glycanSearchMethod = self.config["glycanSearchMethod"]
        glycanSearchMethodNode.setAttribute("glycansearchmethod", glycanSearchMethod)
        glycanNode.appendChild(glycanSearchMethodNode)

        if glycanSearchMethod == "combinatorial":
            combinatorialNode = doc.createElement("combinatorial")
            glycanComponentMinMaxDict = self.config["glycanComponentMinMaxDict"]
            combinatorialDictKeys = ['Hex', 'HexNAc', 'DeoxyHex', 'Neu5Ac',
                                     'Neu5Gc', 'Pentose', 'KDN', 'HexA']
            for key in combinatorialDictKeys:
                componentName = key
                componentMinMaxNode = doc.createElement(key.lower())
                componentMin = glycanComponentMinMaxDict[key]["min"]
                componentMax = glycanComponentMinMaxDict[key]["max"]
                componentMinMaxNode.setAttribute("min", str(componentMin))
                componentMinMaxNode.setAttribute("max", str(componentMax))
                combinatorialNode.appendChild(componentMinMaxNode)

        toleranceNode = doc.createElement("tolerance")
        gpFinderConfigNode.appendChild(toleranceNode)
        
        parentNode = doc.createElement("parent")
        toleranceNode.appendChild(parentNode)
        
        parentPPMNode = doc.createElement("ppm")
        parentPPMTolerance = self.config["tolerancePPM"]
        parentPPMNode.setAttribute("value", str(parentPPMTolerance))
        parentNode.appendChild(parentPPMNode)
        
        fragmentNode = doc.createElement("fragment")
        toleranceNode.appendChild(fragmentNode)
        
        fragmentPPMNode = doc.createElement("ppm")
        fragmentPPMTolerance = self.config["scoringTolerancePPM"]
        fragmentPPMNode.setAttribute("value", str(fragmentPPMTolerance))
        fragmentNode.appendChild(fragmentPPMNode)

        fragmentAbsoluteNode = doc.createElement("absolute")
        fragmentAbsoluteTolerance = self.config["scoringToleranceMass"]
        fragmentAbsoluteNode.setAttribute("value", str(fragmentAbsoluteTolerance))
        fragmentNode.appendChild(fragmentAbsoluteNode)
                       
        gpFinderResultNode = doc.createElement("gpfinderresult")
        xmlNode.appendChild(gpFinderResultNode)
        
        for feature in self.featureList:
            
            featureIndexNode = doc.createElement("feature")
            featureIndexNode.setAttribute("index", str(feature.featureIndex))
            xmlNode.appendChild(featureIndexNode)

            mzNode = doc.createElement("mz")
            mzNode.setAttribute("value", str(feature.compoundMzValue))
            featureIndexNode.appendChild(mzNode)

            zNode = doc.createElement("z")
            zNode.setAttribute("value", str(feature.compoundCharge))
            featureIndexNode.appendChild(zNode)

            neutralMassNode = doc.createElement("neutralmass")
            neutralMassNode.setAttribute("value", str(feature.compoundNeutralMass))
            featureIndexNode.appendChild(neutralMassNode)
            
            rtNode = doc.createElement("rt")
            rtNode.setAttribute("value", str(feature.compoundRtCenter))
            featureIndexNode.appendChild(rtNode)

            intensityNode = doc.createElement("intensity")
            intensityNode.setAttribute("value", str(feature.compoundIntensity))
            featureIndexNode.appendChild(intensityNode)

            fragmentMzListNode = doc.createElement("fragmentmzlist")
            fragmentMzList = [xyzPeak[0] for xyzPeak in feature.xyzPeakList]
            # Following is the encoding scheme
            # list of floats -> double array -> double array string -> base64 string (precision 64)
            fragmentMzListArray = array.array('d', fragmentMzList)
            fragmentMzListString = fragmentMzListArray.tostring()
            fragmentMzListBase64 = base64.standard_b64encode(fragmentMzListString)
            # In order to get the list of floats back:
            # base64 string (precision 64) -> double array string -> double array -> list of floats
            # We need to do:
            '''
            fragmentMzListString = base64.standard_b64decode(fragmentMzListBase64)
            fragmentMzListArray = array.array('d')
            fragmentMzListArray.fromstring(fragmentMzListString)
            fragmentMzList = list(fragmentMzListArray)
            '''
            fragmentMzListTextNode = doc.createTextNode(fragmentMzListBase64)
            fragmentMzListNode.appendChild(fragmentMzListTextNode)
            featureIndexNode.appendChild(fragmentMzListNode)

            # Following is the encoding scheme
            # list of floats -> float array -> float array string -> base64 string (precision 32)
            fragmentIntensityListNode = doc.createElement("fragmentintensitylist")
            fragmentIntensityList = [xyzPeak[1] for xyzPeak in feature.xyzPeakList]
            fragmentIntensityListArray = array.array('f', fragmentIntensityList)
            fragmentIntensityListString = fragmentIntensityListArray.tostring()
            fragmentIntensityListBase64 = base64.standard_b64encode(fragmentIntensityListString)
            # In order to get the list of floats back:
            # base64 string (precision 32) -> float array string -> float array -> list of floats
            # We need to do:
            '''
            fragmentIntensityListString = base64.standard_b64decode(fragmentIntensityListBase64)
            fragmentIntensityListArray = array.array('d')
            fragmentIntensityListArray.fromstring("fragmentIntensityListString)
            fragmentIntensityList = list(fragmentIntensityListArray)
            '''
            fragmentIntensityListTextNode = doc.createTextNode(fragmentIntensityListBase64)
            fragmentIntensityListNode.appendChild(fragmentIntensityListTextNode)
            featureIndexNode.appendChild(fragmentIntensityListNode)
            
            massMatchListNode = doc.createElement("massmatchlist")
            featureIndexNode.appendChild(massMatchListNode)
            
            for massMatch in feature.massMatchList:

                massMatchNode = doc.createElement("massmatch")
                massMatchNode.setAttribute("index", str(massMatch.matchIndex))
                massMatchListNode.appendChild(massMatchNode)

                massMatchProteinFileNameNode = doc.createElement("proteinfilename")
                massMatchProteinFileNameNode.setAttribute("value", str(massMatch.proteinFileName))
                massMatchNode.appendChild(massMatchProteinFileNameNode)
                
                massMatchProteinNameNode = doc.createElement("proteinname")
                massMatchProteinNameNode.setAttribute("value", massMatch.proteinName)
                massMatchNode.appendChild(massMatchProteinNameNode)
                
                massMatchTheoreticalMassNode = doc.createElement("theoreticalmass")
                massMatchTheoreticalMassNode.setAttribute("value", str(massMatch.totalTheoreticalMass))
                massMatchNode.appendChild(massMatchTheoreticalMassNode)

                massMatchTotalScoreNode = doc.createElement("totalscore")
                massMatchTotalScoreNode.setAttribute("value", str(massMatch.totalScore))
                massMatchNode.appendChild(massMatchTotalScoreNode)

                massMatchFragmentMatchesNode = doc.createElement("fragmentmatches")
                massMatchNode.appendChild(massMatchFragmentMatchesNode)

                massMatchPeptideIonFragmentMatchesNode = doc.createElement("peptideionmatches")
                massMatchPeptideIonFragmentMatchesNode.setAttribute("numberoffragmentmatches", str(massMatch.numberOfPeptideIonFragmentMatches))
                massMatchPeptideIonFragmentMatchesNode.setAttribute("scoreperfragmentmatch", "5")
                massMatchFragmentMatchesNode.appendChild(massMatchPeptideIonFragmentMatchesNode)

                massMatchGlycanIonFragmentMatchesNode = doc.createElement("glycanionmatches")
                massMatchGlycanIonFragmentMatchesNode.setAttribute("numberoffragmentmatches", str(massMatch.numberOfGlycanIonFragmentMatches))
                massMatchGlycanIonFragmentMatchesNode.setAttribute("scoreperfragmentmatch", "4")
                massMatchFragmentMatchesNode.appendChild(massMatchGlycanIonFragmentMatchesNode)

                massMatchDeoxyHexIonFragmentMatchesNode = doc.createElement("deoxyhexionmatches")
                massMatchDeoxyHexIonFragmentMatchesNode.setAttribute("numberoffragmentmatches", str(massMatch.numberOfDeoxyHexIonFragmentMatches))
                massMatchDeoxyHexIonFragmentMatchesNode.setAttribute("scoreperfragmentmatch", "3")
                massMatchFragmentMatchesNode.appendChild(massMatchDeoxyHexIonFragmentMatchesNode)

                massMatchNeu5AcIonFragmentMatchesNode = doc.createElement("neu5acionmatches")
                massMatchNeu5AcIonFragmentMatchesNode.setAttribute("numberoffragmentmatches", str(massMatch.numberOfNeu5AcIonFragmentMatches))
                massMatchNeu5AcIonFragmentMatchesNode.setAttribute("scoreperfragmentmatch", "3")
                massMatchFragmentMatchesNode.appendChild(massMatchNeu5AcIonFragmentMatchesNode)

                massMatchNeu5GcIonFragmentMatchesNode = doc.createElement("neu5gcionmatches")
                massMatchNeu5GcIonFragmentMatchesNode.setAttribute("numberoffragmentmatches", str(massMatch.numberOfNeu5GcIonFragmentMatches))
                massMatchNeu5GcIonFragmentMatchesNode.setAttribute("scoreperfragmentmatch", "3")
                massMatchFragmentMatchesNode.appendChild(massMatchNeu5GcIonFragmentMatchesNode)

                massMatchCommonIonFragmentMatchesNode = doc.createElement("commonionmatches")
                massMatchCommonIonFragmentMatchesNode.setAttribute("numberoffragmentmatches", str(massMatch.numberOfCommonIonFragmentMatches))
                massMatchCommonIonFragmentMatchesNode.setAttribute("scoreperfragmentmatch", "1")
                massMatchFragmentMatchesNode.appendChild(massMatchCommonIonFragmentMatchesNode)

                massMatchPeptideStringNode = doc.createElement("peptidestringmatched")
                massMatchPeptideStringNode.setAttribute("value", str(massMatch.peptideStringMatched))
                massMatchNode.appendChild(massMatchPeptideStringNode)

                massMatchGlycosylationSiteIndexNode = doc.createElement("glycosylationsiteindex")
                massMatchGlycosylationSiteIndexNode.setAttribute("value", str(massMatch.glycosylationSiteIndex))
                massMatchNode.appendChild(massMatchGlycosylationSiteIndexNode)

                isDecoyUsed = ""
                if massMatch.isDecoyUsed == True:
                    isDecoyUsed = "1"
                else:
                    isDecoyUsed = "0"

                massMatchIsDecoyUsedNode = doc.createElement("isdecoyused")
                massMatchIsDecoyUsedNode.setAttribute("value", isDecoyUsed)
                massMatchNode.appendChild(massMatchIsDecoyUsedNode)

                massMatchPeptideMassMinusWaterNode = doc.createElement("peptidemassminuswater")
                massMatchPeptideMassMinusWaterNode.setAttribute("value", str(massMatch.peptideMass))
                massMatchNode.appendChild(massMatchPeptideMassMinusWaterNode)

                massMatchParentCompositionNode = doc.createElement("parentcomposition")
                massMatchNode.appendChild(massMatchParentCompositionNode)

                massMatchHexComponentNode = doc.createElement("hex")
                massMatchHexComponentNode.setAttribute("numberofcomponents", str(massMatch.glycanCompositionMatched["Hex"]))
                massMatchParentCompositionNode.appendChild(massMatchHexComponentNode)

                massMatchHexNAcComponentNode = doc.createElement("hexnac")
                massMatchHexNAcComponentNode.setAttribute("numberofcomponents", str(massMatch.glycanCompositionMatched["HexNAc"]))
                massMatchParentCompositionNode.appendChild(massMatchHexNAcComponentNode)

                massMatchDeoxyHexComponentNode = doc.createElement("deoxyhex")
                massMatchDeoxyHexComponentNode.setAttribute("numberofcomponents", str(massMatch.glycanCompositionMatched["DeoxyHex"]))
                massMatchParentCompositionNode.appendChild(massMatchDeoxyHexComponentNode)

                massMatchNeu5AcComponentNode = doc.createElement("neu5ac")
                massMatchNeu5AcComponentNode.setAttribute("numberofcomponents", str(massMatch.glycanCompositionMatched["Neu5Ac"]))
                massMatchParentCompositionNode.appendChild(massMatchNeu5AcComponentNode)

                massMatchNeu5GcComponentNode = doc.createElement("neu5gc")
                massMatchNeu5GcComponentNode.setAttribute("numberofcomponents", str(massMatch.glycanCompositionMatched["Neu5Gc"]))
                massMatchParentCompositionNode.appendChild(massMatchNeu5GcComponentNode)
                
                massMatchPentoseComponentNode = doc.createElement("pentose")
                massMatchPentoseComponentNode.setAttribute("numberofcomponents", str(massMatch.glycanCompositionMatched["Pentose"]))
                massMatchParentCompositionNode.appendChild(massMatchPentoseComponentNode)

                massMatchKDNComponentNode = doc.createElement("kdn")
                massMatchKDNComponentNode.setAttribute("numberofcomponents", str(massMatch.glycanCompositionMatched["KDN"]))
                massMatchParentCompositionNode.appendChild(massMatchKDNComponentNode)

                massMatchHexAComponentNode = doc.createElement("hexa")
                massMatchHexAComponentNode.setAttribute("numberofcomponents", str(massMatch.glycanCompositionMatched["HexA"]))
                massMatchParentCompositionNode.appendChild(massMatchHexAComponentNode)

                massMatchFragmentIonMatchDetailsNode = doc.createElement("fragmentionmatchdetails")
                massMatchNode.appendChild(massMatchFragmentIonMatchDetailsNode)

                fragmentCounterWithinOneParent = 0
                
                for fragmentNeutralMassMatch in massMatch.listOfFragmentNeutralMassMatches:

                    fragmentNeutralMassMatchNode = doc.createElement("fragmentmatch")
                    fragmentCounterWithinOneParent += 1
                    fragmentNeutralMassMatchNode.setAttribute("index", str(fragmentCounterWithinOneParent))
                    massMatchFragmentIonMatchDetailsNode.appendChild(fragmentNeutralMassMatchNode)

                    fragmentNeutralMassMatchFragmentNameNode = doc.createElement("fragmentname")
                    fragmentNeutralMassMatchFragmentNameNode.setAttribute("value", fragmentNeutralMassMatch.matchedFragmentName)
                    fragmentNeutralMassMatchNode.appendChild(fragmentNeutralMassMatchFragmentNameNode)

                    fragmentNeutralMassMatchScoreWithoutMultiplicityNode = doc.createElement("fragmentscorewithoutmultiplicity")
                    fragmentNeutralMassMatchScoreWithoutMultiplicityNode.setAttribute("value", str(fragmentNeutralMassMatch.scoreContributed * fragmentNeutralMassMatch.matchedPeakTuple[3]))
                    fragmentNeutralMassMatchNode.appendChild(fragmentNeutralMassMatchScoreWithoutMultiplicityNode)
                    
                    fragmentNeutralMassMatchScoreMultiplierNode = doc.createElement("fragmentscoremultiplier")
                    fragmentNeutralMassMatchScoreMultiplierNode.setAttribute("value", str(fragmentNeutralMassMatch.matchedPeakTuple[3]))
                    fragmentNeutralMassMatchNode.appendChild(fragmentNeutralMassMatchScoreMultiplierNode)

                    fragmentNeutralMassMatchScoreContributedNode = doc.createElement("fragmentscorecontributed")
                    fragmentNeutralMassMatchScoreContributedNode.setAttribute("value", str(fragmentNeutralMassMatch.scoreContributed))
                    fragmentNeutralMassMatchNode.appendChild(fragmentNeutralMassMatchScoreContributedNode)

                    fragmentNeutralMassMatchActualFragmentMzNode = doc.createElement("fragmentactualmz")
                    fragmentNeutralMassMatchActualFragmentMzNode.setAttribute("value", str(fragmentNeutralMassMatch.matchedActualFragmentMz))
                    fragmentNeutralMassMatchNode.appendChild(fragmentNeutralMassMatchActualFragmentMzNode)

                    fragmentNeutralMassMatchActualFragmentNeutralMassNode = doc.createElement("fragmentactualneutralmass")
                    fragmentNeutralMassMatchActualFragmentNeutralMassNode.setAttribute("value", str(fragmentNeutralMassMatch.matchedActualFragmentNeutralMass))
                    fragmentNeutralMassMatchNode.appendChild(fragmentNeutralMassMatchActualFragmentNeutralMassNode)

        outMatchesToXml.write(doc.toprettyxml(indent = "  "))

        print("Done with xml writing for CSV 1.")
        
        outMatchesToXml.close()

        
                
class Feature():
    
    def __init__(self,
                 featureIndex,
                 compoundNeutralMass,
                 compoundMzValue,
                 compoundCharge,
                 compoundRtCenter,
                 compoundRtTuple,
                 compoundIntensity,
                 collisionEnergy=0):
        
        # deprecated ---
        #a rundant implementation of the input as a dictionary

        self.featureIndex = featureIndex
        self.compoundNeutralMass = compoundNeutralMass
        self.compoundMzValue = compoundMzValue
        self.compoundCharge = compoundCharge
        self.compoundRtCenter = compoundRtCenter
        self.compoundRtTuple = compoundRtTuple
        self.compoundIntensity = compoundIntensity
        self.isolatedMzValue = compoundMzValue
        self.collisionEnergy = collisionEnergy
        self.isotopeCorrection = 0
                    
        #initializing empty list for fragmentation spectrum
        self.xyzPeakList = []
        self.xyzPeakListLower = []
        self.xyzPeakListHigher = []
        
        self.massMatchList = []
    
    def convert_convoluted_peak_list_to_pseudo_deconvoluted(self, maxFragmentCharge):
        if maxFragmentCharge < 0:
            raise InputError("negative charges not supported")
        if maxFragmentCharge == 0:
            pass
        else:
            charges = [charge + 1 for charge in range(maxFragmentCharge)]
        newPeakList = []
        for ion in self.xyzPeakList:
            x = ion[0]
            y = ion[1]
            multiplier = 1.0  #this may need to be changed
            for charge in charges:
                xyzTuple = (x, y, charge, multiplier)
                newPeakList.append(xyzTuple)
        self.xyzPeakList = newPeakList

    #this function adds the x/y system of time and intensity, it also integrates
    def add_inten_coords(self, timeCoordinates, intenCoordinates):
        
        self.timeCoordinates = timeCoordinates
        self.intenCoordinates = intenCoordinates
        
        # base case: only 1 scan contains the target compounds
        #  intensity assumes a pseudo-chromatographic length of 1 second
        if len(timeCoordinates) == 1:
            newInten = intenCoordinates[0]*(1.0/60.0)
            self.intensity = newInten
        
        #sum area under curve using squares and triangles
        else:    
            newInten = 0
            for index in range(len(timeCoordinates - 1)):
                squareHeight = min(intenCoordinates[index],intenCoordinates[index+1])
                triangleHeight = max(intenCoordinates[index],intenCoordinates[index+1]) - squareHeight
                width = timeCoordinates[index+1] - timeCoordinates[index]
                newInten += width*(squareHeight + 0.5*triangleHeight)
            self.intensity = newInten
            

    
    def add_match(self, matchObject):
        
        match = (matchObject, 0)
        self.listOfMassMatchO.append(match)
    
    def combine_spectrum(self, newXyzPeakList, collisionEnergy=0, ppmDifference = 15, absoluteDifference = 0.005):
        """this function will accept a new xyzPeakList and combine it with the current list"""
        if collisionEnergy > self.collisionEnergy:
            self.xyzPeakListLower = deepcopy(self.xyzPeakList)
            self.xyzPeakListHigher = deepcopy(newXyzPeakList)
        else: 
            self.xyzPeakListLower = deepcopy(newXyzPeakList)
            self.xyzPeakListHigher = deepcopy(self.xyzPeakList)
        
        combinedXyzPeakList = []
        currentXyzPeakList = deepcopy(self.xyzPeakList)
        print "combining spectra of length {} and {}".format(len(currentXyzPeakList),len(newXyzPeakList))
        
        currentPeak = -1
        currentPeakMz = -1
        newPeak = -1
        newPeakMz = -2
        while len(currentXyzPeakList) > 0 or len(newXyzPeakList) > 0:
            if newPeakMz == -2 and len(newXyzPeakList) > 0:
                newPeak = newXyzPeakList.pop()
                newPeakMz = newPeak[0]
            if currentPeakMz == -1 and len(currentXyzPeakList) > 0:
                currentPeak = currentXyzPeakList.pop()
                currentPeakMz = currentPeak[0]
            
            if newPeakMz > 0 and currentPeakMz > 0:
                newOldDiff = math.fabs(currentPeakMz - newPeakMz)
                newOldDiffPpm = 1000000*newOldDiff/newPeakMz
                #if less than 15 ppm different or 0.005 Thompsons
                if newOldDiffPpm < ppmDifference or newOldDiff < absoluteDifference:
                    realMz = (newPeakMz+currentPeakMz)/2.0
                    realInten = (newPeak[1] + currentPeak[1])/2.0
                    newXyzList = (realMz, realInten, newPeak[2], newPeak[3])
                    combinedXyzPeakList.append(newXyzList)
                    #reset the current tabs (since they have been added)
                    newPeakMz = -2
                    currentPeakMz = -1
                    continue
            
            #since they are not roughly equal, append the larger one
            if newPeakMz != -2 and newPeakMz > currentPeakMz:
                combinedXyzPeakList.append(newPeak)
                newPeakMz = -2
            elif  currentPeakMz != -1: 
                combinedXyzPeakList.append(currentPeak)
                currentPeakMz = -1
        
        combinedXyzPeakList.reverse()
        combinedXyzPeakList.sort()
        self.xyzPeakList = combinedXyzPeakList
        print "combined spectra are of length {}".format(len(combinedXyzPeakList))
    
    #def score_compound_matches
        #run scorring algorthm, sort via score
    
    #def output_xml
                
#--------------- test script below --------------            
'''
for feature in ss2features.featureList:
    print(' ')
    for entry in feature.dictionary:
        print(str(entry) + ": " + str(feature.dictionary[entry])) '''
