from os.path import splitext
from xml.dom import minidom
from xml.dom.minidom import parse
from xml.dom.minidom import parseString
from xml.dom.minidom import Document
import time
import math

# Python 2 only
 
class findMassMatch:

    def __init__(self, config):
        
        # set program configuration dictionary as instance variable
        self.config = config

        self.isSequenceGivenByUser = config["isSequenceGivenByUser"]
        # The dictionary of all smaller components names and masses (sugars mostly)
        # We can expand the list to include more smaller components if we want
        
        # ATTENTION: glycan component masses here are dehydrated already
        self.glycanComponentMassDict = {}

        self.waterMass = 18.0105646864
        self.hydrogenMass = 0

        self.AMINO_ACID_NAME_LENGTH_ADDITION = 8
        # The Min and Max of all smaller components (sugars mostly)
        #     that can be chosen
        #     are all set to 0 to be the default value
        #     And we only assign them to non-negative values when
        #     we are in "combinatorial" method
        
        # Setting up all member variables

        '''self.decoy = 0
        self.decoyMin = 0
        self.decoyMax = 0'''
      
        # If we call "findMassMatch.to_pretty_xml(self)"
        #     The new name is having the same extension
        #     But it has the base name with "_pretty" appended
        # Else if we do not call that
        #     The new name is an empty string
        self.newInputXmlFileName = ""
        self.inputMassToMatch = 0

        # All Post-Translational Modifications (PTM's)
        self.listOfPhosphorylationTypes = []
        self.listOfCarbamidomethylationTypes = []
        self.listOfDeamidationTypes = []
        self.listOfOxidationTypes = []
        
        self.listOfPhosphorylationSitePositions = []
        self.listOfCarbamidomethylationSitePositions = []
        self.listOfDeamidationSitePositions = []
        self.listOfOxidationSitePositions = []
        
        self.phosphorylationSitePositionDict = {}
        self.carbamidomethylationSitePositionDict = {}
        self.deamidationSitePositionDict = {}
        self.oxidationSitePositionDict = {}
        
        self.listOfGlycosylationSitePositions = []
        self.aminoAcidSequence = ""
        
        #Data structure for aminoAcidSequenceDict:
        #    aminoAcidSequenceDict = {1:'S', 2:'K', 3:"PT", 4:'V', 5:'M', ...}

        self.aminoAcidSequenceDict = {}
            
        self.aminoAcidSequenceLength = 0
        self.glycanDict = {}
        '''self.listOfMassIndexPairs = []'''
        self.numberOfGlycansInLibrary = 0
        self.aminoAcidDict = {}

        self.listOfFragmentLeftRightIndexPairs = []
        
        self.globalCounterForMatches = 0

        '''self.allMassMatches = []'''                                               

        # Scoring tolerance in ppm
        self.scoringTolerance = 0
        
        # Comment out this line to skip the XML file formatting
        # Although minidom does NOT preserve the order of attributes
        # It does NOT hurt our desired functionalities
        # to_pretty_xml() is not needed in the functionarlity
        # It only helps to give a pretty layout of the xml file
        # self.to_pretty_xml()

        self.read_amino_acid_monoisotopic_masses_csv()

        self.read_glycan_component_monoisotopic_masses_csv()

        self.read_user_input_xml()

        self.read_hydrogen_mass_for_reduced_glycan()

        if self.config["glycanSearchMethod"]  == "scottslibrary":
            self.read_scotts_library_csv()
        elif self.config["glycanSearchMethod"]  == "genericlibrary":
            self.read_generic_glycan_library_csv()
        
        # If the reduced attribute is NOT checked,
        #     we do the normal glycopeptide mass matching
        if not self.config['free']:

            if self.config["uniprot"] == False:
                self.aminoAcidSequenceDict = config["sequenceDict"]
                self.proteinName = "Sequence " + self.config["sequence"]
                self.proteinFileName = "Sequence " + self.config["sequence"]
                self.find_potential_phosphorylation_site_positions_from_user_input()

            if self.config["uniprot"] == True:
                self.read_protein_xml()

            self.find_potential_carbamidomethylation_site_positions_from_user_input()
            self.find_potential_deamidation_site_positions_from_user_input()
            self.find_potential_oxidation_site_positions_from_user_input()
            
            if ("trypsin" in self.config["digestTypes"]
                    and "lysc" in self.config["digestTypes"]
                    and "gluc" in self.config["digestTypes"]):
                self.trypsin_and_lysc_and_gluc_cleavages()
            elif ("trypsin" in self.config["digestTypes"]
                      and "lysc" in self.config["digestTypes"]):
                self.trypsin_and_lysc_cleavages()
            elif ("trypsin" in self.config["digestTypes"]
                      and "elastase" in self.config["digestTypes"]):
                self.trypsin_and_elastase()
            elif ("trypsin" in self.config["digestTypes"]
                    and "elastase" in self.config["digestTypes"]
                    and "chymotrypsin" in self.config["digestTypes"]):
                self.trypsin_and_elastase_and_chymotrypsin()
            elif ("chymotrypsin" in self.config["digestTypes"]
                    and "elastase" in self.config["digestTypes"]):
                self.chymotrypsin_and_elastase()
            elif "pepsin" in self.config["digestTypes"]:
                self.pepsin_cleavages()
            elif "elastase" in self.config["digestTypes"]:
                self.elastase_cleavages()
            elif ("trypsin" in self.config["digestTypes"]
                      and "gluc" in self.config["digestTypes"]):
                self.trypsin_and_gluc_cleavages()
            elif ("lysc" in self.config["digestTypes"]
                      and "gluc" in self.config["digestTypes"]):
                self.lysc_and_gluc_cleavages()
            elif "trypsin" in self.config["digestTypes"]:
                self.trypsin_cleavages()
            elif "lysc" in self.config["digestTypes"]:
                self.lysc_cleavages()
            elif "gluc" in self.config["digestTypes"]:
                self.gluc_cleavages()
            elif "argc" in self.config["digestTypes"]:
                self.argc_cleavages()
            elif "chymotrypsin" in self.config["digestTypes"]:
                self.chymotrypsin_cleavages()
                
            
        # If the reduced attribute is checked,
        #     we do the glycan mass matching without the peptide chain
        else:
            self.proteinFileName = "Glycan only"
            self.proteinName = "Glycan only"
        


    """
    def to_pretty_xml(self):

        # Function Description: 
        #     Read the content in the xml file
        #     Add pretty indentations
        #     Trim redundant spaces at the end of each line
        #     Write to pretty content to a new xml file that ends with "_pretty.xml"
        # Function Input:
        #     Original XML file
        # Function Output:
        #     New XML file with nice formatting

        # ------ Pass member variable values to local variables ------
        # ------ It boosts performance, as suggested by Evan Parker ------
        # ------ See link below: ------
        # ------ http://wiki.python.org/moin/PythonSpeed/PerformanceTips#Local_Variables------
        inputXmlFileName = self.config["proteinFile"]
        oldInputXmlFileName = inputXmlFileName
        
        print("Converting the file " + inputXmlFileName +
              " to a pretty indented look without redundant spaces at the end of each line...")
 
        # ------ Open XML file and parse it into an XML tree of nodes ------
        inUglyXml = open(inputXmlFileName, 'r')
        uglyXmlDocument = parse(inUglyXml)
        inUglyXml.close()

        # ------ And "_pretty.xml" to the new XML file name ------
        newInputXmlFileName = splitext(inputXmlFileName)[0] + "_pretty.xml"
        outPrettyXml = open(newInputXmlFileName, 'w')

        # ------ Set up formatting, using 2-space indentation ------
        content = uglyXmlDocument.toprettyxml(indent = "  ",
                                              newl = "",
                                              encoding = "utf-8")

        # Trivial: Make a new line at the <?xml> line
        content = str(content).replace("?>", "?>\n", 1)

        # outPrettyXml.write(str(content))
        # ------ Write all the lines into the new XML file
        lines = content.split("\n")
        for lineIndex in range(len(lines)):
            lines[lineIndex] = lines[lineIndex].rstrip()
            outPrettyXml.write(lines[lineIndex] + "\n")
            
        outPrettyXml.close()
        
        print("Done.\n\n\n")

        # ------ Pass local variable values back to member variables ------
        self.newInputXmlFileName = newInputXmlFileName
    """
        
        
    def read_protein_xml(self):

        # Function Description:
        #     Read the XML that contains protein information
        #     to get
        #     the list of glycosylation sites,
        #     the list of phosphorylation sites,
        #     the amino acid sequence, and
        #     the amino acid sequence length
        # Function Input:
        #     XML file that contains protein information
        # Function Output:
        #     the list of glycosylation sites,
        #     the list of phosphorylation sites,
        #     the amino acid sequence, and
        #     the amino acid sequence length
        
        # ------ Pass member variable values to local variables ------
        inputXmlFileName = self.config["proteinFile"]
        config = self.config

        proteinFileName = inputXmlFileName.split("/")[-1]
        # all positions in the list start at 1
        listOfGlycosylationSitePositions = self.listOfGlycosylationSitePositions
        phosphorylationSitePositionDict = self.phosphorylationSitePositionDict
        listOfPhosphorylationSitePositions = self.listOfPhosphorylationSitePositions
        aminoAcidSequenceLength = self.aminoAcidSequenceLength
        aminoAcidSequence = self.aminoAcidSequence
        aminoAcidSequenceDict = self.aminoAcidSequenceDict
        aminoAcidDict = self.aminoAcidDict
            
        print("Reading the file " + inputXmlFileName + "...")

        # ------ Open XML file and parse it into an XML tree of nodes ------
        
        inXml = open(inputXmlFileName, 'r')
        xmlDocument = parse(inXml)

        print("Done.\n\n\n")

        # ------ Get the recommended full name of this protein ------
        
        recommendedNameNode = xmlDocument.getElementsByTagName("recommendedName")[0]
        recommendedFullNameNode = recommendedNameNode.getElementsByTagName("fullName")[0]
        recommendedFullName = recommendedFullNameNode.firstChild.data.strip()
        listOfProteinNames = self.config["listOfProteinNames"]
        
        if recommendedFullName not in listOfProteinNames:
            listOfProteinNames.append(recommendedFullName)
            self.config["listOfProteinNames"] = listOfProteinNames

        # ------ Get all nodes with name "feature" ------
        
        featureNodeList = xmlDocument.getElementsByTagName("feature")
        
        print("Obtaining a list of glycosylation site positions...")

        # ------ Get the list of glycosylation sites ------
        
        # For every node:
        #     Get the 1st location node with attribute "glycosylation site"
        #         but withOUT description "glycation"
        #     Get the 1st position node
        #     Store all integer position indices
        #         to the list of glycosylation sites
        #         according to glycosylation types that user configured to search for
        for featureNode in featureNodeList:
            if "glycosylation site" in featureNode.getAttribute("type"):
                if "glycation" not in featureNode.getAttribute("description"):
                    
                    locationNode = featureNode.getElementsByTagName("location")[0]
                    positionNode = locationNode.getElementsByTagName("position")[0]
                    positionAttribute = int(positionNode.getAttribute("position").strip())
                    print(str(positionAttribute))

                    # all positions in the list start at 1
                    listOfGlycosylationSitePositions.append(positionAttribute)
                    
        print("Done.\n")

        if config["searchForNLinked"] == True and config["searchForOLinked"] == False:
            print("Glycosylation types to search: N-Linked.")
        elif config["searchForNLinked"] == False and config["searchForOLinked"] == True:
            print("Glycosylation types to search: O-Linked.")
        else: # both True, because both False has been excluded above already.
            print("Glycosylation types to search: N-Linked and O-Linked (result may be confusing).")
            
        print("The glycosylation site positions (those start at 1) are:")

        # Trivial: Set up a string to display the glycosylation sites
        displayString = ", ".join([str(glycosylationSitePositions)
                                   for glycosylationSitePositions
                                   in listOfGlycosylationSitePositions])

        if len(listOfGlycosylationSitePositions) == 0:
            displayString = ""
            
        print(displayString + "\n\n\n")

        print("Obtaining a dictionary of phosphorylation sites...")
        
        # ------ Get the list of phosphorylation sites ------

        # For every node:
        #     Get the 1st location node that has attribute "modified residue"
        #     Get the 1st position node
        #     Store all integer position indices
        #         to the dictionary of phosphorylation sites
        #         in a format of
        #         phosphorylationSitePositionDict = {"Phosphoserine": [35, 68, 93, 139],
        #                                            "Phosphothreonine: [75, 115, 150],
        #                                            "Phosphotyrosine": [90, 190, 223],
        #                                            ...}
        #         Although at this point,
        #             we only need phosphoserine and phosphothreonine,
        #         We keep all of them here to allow possible expansion
        #             to consider more of them in the future
        for featureNode in featureNodeList:
            if featureNode.getAttribute("type") == "modified residue":
                description = featureNode.getAttribute("description").strip()
                if "Phospho" in description:
                    phosphorylationType = description
                    locationNode = featureNode.getElementsByTagName("location")[0]
                    positionNode = locationNode.getElementsByTagName("position")[0]
                    positionAttribute = int(positionNode.getAttribute("position").strip())
                    listOfPhosphorylationSitePositions.append(positionAttribute)

                    for aminoAcid in aminoAcidDict.keys():
                        aminoAcidName = aminoAcidDict[aminoAcid]["Name"]
                        if (aminoAcidName.lower() in phosphorylationType.lower() and
                                    "Phospho" in aminoAcidName):
                            phosphorylationType = aminoAcidName
                            if phosphorylationType in phosphorylationSitePositionDict.keys():
                                phosphorylationSitePositionDict[phosphorylationType].append(
                                                                        positionAttribute)
                            else:
                                phosphorylationSitePositionDict[phosphorylationType] = [
                                                                        positionAttribute]
                            break
                    
        print("Done.\n")

        print("The phosphorylation site positions (those start at 1) are:")

        # Trivial: Set up a string to display the phosphorylation sites
        for phosphorylationType in phosphorylationSitePositionDict.keys():
            displayString = ", ".join([str(phosphorylationSitePosition)
                                       for phosphorylationSitePosition
                                       in phosphorylationSitePositionDict[
                                               phosphorylationType]])
            print(phosphorylationType + ": " + displayString)

        print("\n")
        
        print("Obtaining the amino acid sequence...")

        # ------ Get the amino acid sequence from the XML file ------

        # Get the 1st node with name "entry"
        # Get the 1st node with name "sequence"
        entryNode = xmlDocument.getElementsByTagName("entry")[0]        
        sequenceNodeList = entryNode.getElementsByTagName("sequence")

        # For every node:
        #     Get the length attribute if the node has length attribute
        #     Get the text data of the sequence
        for sequenceNode in sequenceNodeList:
            if sequenceNode.hasAttribute("length"):
                aminoAcidSequenceLength = int(sequenceNode.getAttribute("length").strip())
                sequenceUnstripped = sequenceNode.firstChild.data.strip()
        
        print("Done.\n")

        print("The amino acid sequence is:")

        # ------ Delete spaces and newlines in the amino acid sequence ------
        aminoAcidSequence = ''.join([char
                                     for char in sequenceUnstripped
                                     if char.isalpha()])
         
        print(aminoAcidSequence + "\n")

        print("The amino acid sequence length is:")

        print(str(aminoAcidSequenceLength) + "\n\n\n")

        # ------ Check whether the amino acid sequence data length ------
        # ------ is the same as the length attribute read ------
        # ------ If not, there must be an error ------
        if len(aminoAcidSequence) != aminoAcidSequenceLength:
            print("ERROR!! Amino Acid Sequence Length does NOT match the length attribute!!")

        # ------ Make a dictionary fo all the amino acids in the sequence ------
        # ------ dictionary key: the string position of the amino acid starting at 1 ------
        # ------ dictionary value: the string value of the amino acid based on the key ------
        for i in range(len(aminoAcidSequence)):
            aminoAcidSequenceDict[i + 1] = aminoAcidSequence[i]

        # ------ Change the amino acid to be phosphorylated based on uniprot XML ------
        for position in listOfPhosphorylationSitePositions:
            aminoAcidSequenceDict[position] = aminoAcidSequenceDict[position] + "^phos'ed"

        if (config["searchForNLinked"] == True and
                config["searchForOLinked"] == False):
            
            listOfGlycosylationSitePositions = [sitePosition
                                                for sitePosition in listOfGlycosylationSitePositions
                                                if (aminoAcidSequenceDict[sitePosition] == "N" or
                                                        aminoAcidSequenceDict[sitePosition] == "N^phos'ed")]

        elif (config["searchForNLinked"] == False and
                  config["searchForOLinked"] == True):
            
            listOfGlycosylationSitePositions = [sitePosition
                                                for sitePosition in listOfGlycosylationSitePositions
                                                if (aminoAcidSequenceDict[sitePosition] != "N" and
                                                        aminoAcidSequenceDict[sitePosition] != "N^phos'ed")]

        elif (config["searchForNLinked"] == True and
                  config["searchForOLinked"] == True):

            pass

        else:
            # both False
            if config["free"] == False:

                raise Exception("searchForNLinked and searchForOLinked cannot be both unchecked except you're doing free glycan searching.")
            
        # ------ Pass local variable values back to member variables ------
        self.proteinFileName = proteinFileName
        self.proteinName = recommendedFullName
        self.listOfGlycosylationSitePositions = listOfGlycosylationSitePositions
        self.phosphorylationSitePositionDict = phosphorylationSitePositionDict 
        self.aminoAcidSequenceLength = aminoAcidSequenceLength
        self.aminoAcidSequence = aminoAcidSequence
        self.aminoAcidSequenceDict = aminoAcidSequenceDict



    def read_scotts_library_csv(self):

        # ATTENTION: The mass values in the csv are not yet dehydrated

        # Function Description:
        #     Read the CSV library file that contains glycan information
        #     Store all the values for every glycan
        #     Into a dictionary of dictionaries
        # Function Input:
        #     The CSV library file that contains glycan information
        # Function Output:
        #     The dictionary of dictionaries that contains all useful information
        #     Format shown below
        
        # ------ Pass member variable values to local variables ------
        glycanDict = self.glycanDict
        '''listOfMassIndexPairs = self.listOfMassIndexPairs'''
        numberOfGlycansInLibrary = self.numberOfGlycansInLibrary
        
        print("Reading Scott's Library...")

        # ------ Open and read the CSV file ------
        inLibrary = open(self.config["glycanLibraryFile"], 'r')
        # Skip the first 2 lines
        line = inLibrary.readline()
        line = inLibrary.readline()

        # ------ From every line in the CSV file ------
        # ------ Store all useful values: ------
        # ------ Index, Molecular Mass, ------
        # ------ Hex, HexNac, DxyHex, NeuAc, ------
        # ------ Carbohydrate Type, Glycan Code, Class. ------
        # ------ into a dictionary "glycanDict" of dictionaries "tableEntryDicts" ------
        # ------ in the format of ------
        # ------ glycanDict = {"1": {"Molecular Mass": 910.32778,
        # ------                     "Hex": 3,
        # ------                     "HexNac": 2,
        # ------                     "DxyHex": 0,
        # ------                     "neu5ac": 0,
        # ------                     "Carbohydrate Type": "Aldehyde",
        # ------                     "Glycan Code": "32000",
        # ------                     "Class": "Core"},
        # ------               "2": {"Molecular Mass": 1056385688,
        # ------                     "Hex": 3,
        # ------                     "HexNac": 2,
        # ------                     "DxyHex": 1,
        # ------                     "neu5ac": 0,
        # ------                     "Carbohydrate Type": "Aldehyde",
        # ------                     "Glycan Code": "32100",
        # ------                     "Class": "High Mannose"},
        # ------               "3": {...},
        # ------               "4": {...},
        # ------               ...,
        # ------               "Last index (331 now but may change)": {...}
        # ------               }
        for line in inLibrary:

            # Some values are converted to float or int necessarily to
            # make the future manipulation of data more intuitive and facilitated
            line = line.strip()
            fields = line.split(",")
            [index,molecularMass,
             Hexose,HexNAc,DeoxyHex,Neu5Ac,
             carbType,glycanCode,glycanClass,
             glycanFormula,glycanMassReference] = fields
                
            # Attach all key-value pairs in the dictionary "tableEntryDict"
            tableEntryDict = {}
            tableEntryDict["Molecular Mass"] = float(molecularMass)
            tableEntryDict["Hex"] = int(Hexose)
            tableEntryDict["HexNAc"] = int(HexNAc)
            tableEntryDict["DeoxyHex"] = int(DeoxyHex)
            tableEntryDict["Neu5Ac"] = int(Neu5Ac)
            tableEntryDict["Neu5Gc"] = 0
            tableEntryDict["Pentose"] = 0
            tableEntryDict["KDN"] = 0
            tableEntryDict["HexA"] = 0
            #tableEntryDict["Sulf"] = 0
            tableEntryDict["Carbohydrate Type"] = carbType
            tableEntryDict["Glycan Code"] = glycanCode
            tableEntryDict["Class"] = glycanClass
            tableEntryDict["Chemical Formula"] = glycanFormula

            self.decoy = 0

            # Add the key(index)-value(dictionary) pair in the dictionary "glycanDict"
            glycanDict[index] = tableEntryDict
            '''listOfMassIndexPairs.append((molecularMass, index))'''

            # Count the number of glycans in the library
            numberOfGlycansInLibrary = numberOfGlycansInLibrary + 1

        '''listOfMassIndexPairs.sort()'''
        
        print("Done.\n\n\n")
        
        # ------ Pass local variable values back to member variables ------               
        self.glycanDict = glycanDict
        '''self.listOfMassIndexPairs = listOfMassIndexPairs'''
        self.numberOfGlycansInLibrary = numberOfGlycansInLibrary


       
    def read_generic_glycan_library_csv(self):
        #
        # A new glycan library reader function, see the previous function (read_scotts_library_csv) for further documentation
        #
        glycanDict = self.glycanDict
        '''listOfMassIndexPairs = self.listOfMassIndexPairs'''
        numberOfGlycansInLibrary = self.numberOfGlycansInLibrary
        
        print("Reading generic Library...")

        # ------ Open and read the CSV file ------
        inLibrary = open(self.config["glycanLibraryFile"], 'r')
        # Skip the first line
        line = inLibrary.readline()
        
        for line in inLibrary:

            # Some values are converted to float or int necessarily to
            # make the future manipulation of data more intuitive and facilitated
            line = line.strip()
            fields = line.split(",")
            if len(fields) >= 9:
                [index,Hexose,HexNAc,DeoxyHex,Neu5Ac,Neu5Gc,
                 Pentose, KDN, HexA] = fields[0:9]
                

                
                # Attach all key-value pairs in the dictionary "tableEntryDict"
                tableEntryDict = {}
                tableEntryDict["Hex"] = int(Hexose)
                tableEntryDict["HexNAc"] = int(HexNAc)
                tableEntryDict["DeoxyHex"] = int(DeoxyHex)
                tableEntryDict["Neu5Ac"] = int(Neu5Ac)
                tableEntryDict["Neu5Gc"] = int(Neu5Gc)
                tableEntryDict["Pentose"] = int(Pentose)
                tableEntryDict["KDN"] = int(KDN)
                tableEntryDict["HexA"] = int(HexA)
                #tableEntryDict["Sulf"] = 0
                tableEntryDict["Carbohydrate Type"] = ""
                tableEntryDict["Glycan Code"] = "".join([Hexose, HexNAc, DeoxyHex, Neu5Ac, "0"])
                tableEntryDict["Class"] = ""
                tableEntryDict["Chemical Formula"] = "Generic"
                tableEntryDict["Molecular Mass"] = self.get_glycan_mass(tableEntryDict)
                
                # Add the key(index)-value(dictionary) pair in the dictionary "glycanDict"
                glycanDict[index] = tableEntryDict
                
                numberOfGlycansInLibrary = numberOfGlycansInLibrary + 1
        
        # ------ Pass local variable values back to member variables ------ 
        if len(glycanDict) == 0:
            raise IOError("generic glycan library contains no glycans!")
        
        self.glycanDict = glycanDict

        '''self.listOfMassIndexPairs = listOfMassIndexPairs'''
        self.numberOfGlycansInLibrary = numberOfGlycansInLibrary
        

        
    def read_user_input_xml(self):

        # Function Description:
        #     Read GP3Config.xml file that is generated by the online form
        # Function Input:
        #     The GP3Config.xml file that is generated by the online form
        # Function Output:
        #     None
        
        # ------ Pass member variable values to local variables ------
        inputXmlFileName = self.config["proteinFile"]
        glycanComponentMassDict = self.glycanComponentMassDict
        glycanComponentMinMaxDict = self.config["glycanComponentMinMaxDict"]
        
        listOfPhosphorylationTypes = []
        listOfCarbamidomethylationTypes = []
        listOfDeamidationTypes = []
        listOfOxidationTypes = []
        
        aminoAcidSequenceDict = {}
        aminoAcidSequence = self.aminoAcidSequence
        aminoAcidSequenceLength = self.aminoAcidSequenceLength
        listOfFragmentLeftRightIndexPairs = self.listOfFragmentLeftRightIndexPairs
        
        print("The maximum tolerance in PPM allowed in the matching is: " +
                                          str(self.config["tolerancePPM"]) + ".\n")
        
        if not self.config['free']:

            digestTypeDisplay = ", ".join([digestType
                                           for digestType
                                           in self.config["digestTypes"]])
            
            print("The digestion method is: " + digestTypeDisplay + ".\n\n\n")

            print("Done.\n")

            if self.config["uniprot"] == False:
                
                print("The amino acid sequence is:")

                print(self.config["sequence"] + "\n")
                
                aminoAcidSequenceLength = len(self.config["sequence"])
                
                for i in range(len(self.config["sequence"])):
                    aminoAcidSequenceDict[i + 1] = self.config["sequence"][i]
                    
                # ------ Pass local variable values back to member variables ------
                
                self.aminoAcidSequenceLength = aminoAcidSequenceLength
                self.aminoAcidSequenceDict = aminoAcidSequenceDict

                if self.config["autoDetectNSites"] == True:
                    self.find_potential_n_linked_glycosylation_site_positions_from_input_string()            
                else:
                    self.listOfGlycosylationSitePositions = self.config["listOfSites"]
                    
                # getting phosphorylation sites
                print("Obtaining the list of phosphorylation types to be considered...")
                for phosphorylationTuple in self.config["listOfPhosphorylationTuples"]:
                    if phosphorylationTuple[1]:
                        listOfPhosphorylationTypes.append(phosphorylationTuple[0])                       
                print("Done.\n")
                print("The types of phosphorylation sites to be considered are:")
                
                displayString = ", ".join([phosphorylationType
                                           for phosphorylationType
                                           in listOfPhosphorylationTypes])
                if len(listOfPhosphorylationTypes) == 0:
                    displayString = "None."
                print(displayString + "\n\n\n")
                
                # ------ Pass local variable values back to member variables ------
                self.listOfPhosphorylationTypes = listOfPhosphorylationTypes
                self.aminoAcidSequence = self.config["sequence"]
                
            elif self.config["uniprot"] == True:
                
                print("\n\nUser input amino acid sequence is empty.\n\n\n")

            # getting carbamidomethylation sites
            print("Obtaining the list of carbamidomethylation types to be considered...")
            for carbamidomethylationTuple in self.config["listOfCarbamidomethylationTuples"]:
                if carbamidomethylationTuple[1] == True:
                    listOfCarbamidomethylationTypes.append(carbamidomethylationTuple[0])
            print("Done.\n")
            print("The types of carbamidomethylation sites to be considered are:")
                    
            displayString = ", ".join([carbamidomethylationType
                                        for carbamidomethylationType
                                        in listOfCarbamidomethylationTypes])
            if len(listOfCarbamidomethylationTypes) == 0:
                displayString = "None."
            print(displayString + "\n\n\n")

            # getting deamidation sites
            print("Obtaining the list of deamidation types to be considered...")
            for deamidationTuple in self.config["listOfDeamidationTuples"]:
                if deamidationTuple[1] == True:
                    listOfDeamidationTypes.append(deamidationTuple[0])
            print("Done.\n")
            print("The types of deamidation sites to be considered are:")

            displayString = ", ".join([deamidationType
                                       for deamidationType
                                       in listOfDeamidationTypes])
            if len(listOfDeamidationTypes) == 0:
                displayString = "None."
            print(displayString + "\n\n\n")

            # getting oxidation sites
            print("Obtaining the list of oxidation types to be considered...")
            for oxidationTuple in self.config["listOfOxidationTuples"]:
                if oxidationTuple[1] == True:
                    listOfOxidationTypes.append(oxidationTuple[0])
            print("Done.\n")
            print("The types of oxidation sites to be considered are:")

            displayString = ", ".join([oxidationType
                                       for oxidationType
                                       in listOfOxidationTypes])
            if len(listOfOxidationTypes) == 0:
                displayString = "None."
            print(displayString + "\n\n\n")

            # ------ Pass local variable values back to member variables ------
            self.listOfCarbamidomethylationTypes = listOfCarbamidomethylationTypes
            self.listOfDeamidationTypes = listOfDeamidationTypes
            self.listOfOxidationTypes = listOfOxidationTypes
                
        else:
            pass
        
        print("Done.\n\n\n")

        

    def read_glycan_component_monoisotopic_masses_csv(self):

        # Function Description:
        #      Read the CSV file that contains glycan monoisotopic mass information
        #      Store the key(glycan name)-value(mass) pairs in a dictionary
        # Function Input:
        #     The CSV file that contains glycan monoisotopic mass information
        # Function Output:
        #     The dictionary of key(glycan name)-value(mass) pairs
        
        print("Reading " + self.config["saccharideFile"] + "...")
        print("Obtaining glycan component monoisotopic masses...")
        
        glycanComponentMassDict = self.glycanComponentMassDict

        # ------ Open and read the csv file ------
        inGlycanComponentMasses = open(self.config["saccharideFile"], 'r')
        inGlycanComponentMasses.readline()

        # Get the glycan component names and glycan component masses
        #     to build a dictionary of all glycan components
        for line in inGlycanComponentMasses:
            
            line = line.strip()
            fields = line.split(',')
            
            glycanComponentName = fields[0]
            glycanComponentMass = float(fields[1])
            
            glycanComponentMassDict[glycanComponentName] = glycanComponentMass

        print("Done.\n\n\n")

        # ------ Pass local variable values back to member variables ------
        self.glycanComponentMassDict = glycanComponentMassDict

        
        
    def read_amino_acid_monoisotopic_masses_csv(self):

        # ATTENTION: The mass values in the csv are dehydrated already

        # Function Description:
        #     Read the CSV file that contains amino acid monoisotopic mass information
        #     Store the key(amino acid name)-value(mass) pairs in a dictionary
        # Function Input:
        #     The CSV file that contains amino acid monoisotopic mass information
        # Function Output:
        #     The dictionary of key(amino acid name)-value(mass) pairs

        # ------ Pass member variable values to local variables ------
        aminoAcidDict = self.aminoAcidDict
        
        print("Reading " + self.config["aaMassFile"] + "...")
        print("Obtaining Amino Acid Molecular Masses...")

        # ------ Open and read the CSV file ------
        inAminoAcidMasses = open(self.config["aaMassFile"], 'r')
        inAminoAcidMasses.readline()

        # ------ Read from the file and
        # ------ stores the key(amino acid name)-value(mass) pairs
        # ------ in a dictionary
        for line in inAminoAcidMasses:

            line = line.strip()
            fields = line.split(",")
            aminoAcidEntryDict = {}
            acronym, abbreviation, name, mass, formula, reference = fields
            aminoAcidEntryDict["Abbreviation"] = abbreviation
            aminoAcidEntryDict["Name"] = name
            aminoAcidEntryDict["Mass"] = float(mass)
            #aminoAcidEntryDict["Chemical Formula"] = formula
            #aminoAcidEntryDict["Mass Reference"] = reference
            aminoAcidDict[acronym] = aminoAcidEntryDict

        print("Done.\n\n\n")
        '''
        # If cysteine carbamidomethylation is applied:
        #     Modified Cys and Phospho-Cys masses
        # Principle:
        #     One peptide can have any two Cysteines (Cys, C) connected
        #         which makes the fragmentation analysis difficult
        #     In order to prevent this S-S bond between two C's,
        #     We use Dithothreitol and Iodoacetamide
        #         to perform Cysteine carbamidomethylation
        #     This is done by replacing the H- from HS- in Cysteine
        #         with NH2COCH2-
        #     The mass is increased by 57.0214637238
        #         which makes C (Cysteine) mass 160.0306482031
        #             instead of 103.0091844793
        #         and makes PC (Phosphocysteine) mass 239.9969786136
        #             instead of 182.9755148898
        if self.config["applyCysteineCarbamidomethylation"] == True:
            aminoAcidDict["C"]["Mass"] = 160.0306482031
            aminoAcidDict["PC"]["Mass"] = 239.9969786136
        '''
            
        # ------ Pass local variable values back to member variables ------
        # ------ also pass back to config so it is globally available -----
        aminoAcidDict[""] = {'Abbreviation': 'None', 'Mass': 0, 'Name': 'None'}
        self.config["aminoAcidMassDict"] = aminoAcidDict
        self.aminoAcidDict = aminoAcidDict
        


    def read_water_mass(self):

        # Function Description:
        #     Read in the water mass
        # Function Input:
        #     The CSV file that contains the water mass
        # Function Output:
        #     None

        inWater = open("water_mass.csv", 'r')

        line = inWater.readline()
        line = inWater.readline()
        
        line = line.strip()
        line = line.split(',')

        self.waterMass = float(line[1])


        
    def read_hydrogen_mass_for_reduced_glycan(self):

        # Function Description:
        #     Read in the hydrogen mass
        # Function Input:
        #     The CSV file that contains the hydrogen mass
        # Function Output:
        #     None
        
        print("mass used for glycan addition: "),
        if self.config['free']:
            if self.config['reduced'] and self.config['deuterated']:
                self.hydrogenMass = 2.0141017778
            elif self.config['reduced'] and not self.config['deuterated']:
                self.hydrogenMass = 1.0078250321
            else:
                self.hydrogenMass = 0
               
        else:
            self.hydrogenMass = 0
        
        print(self.hydrogenMass)


        
    def find_potential_n_linked_glycosylation_site_positions_from_input_string(self):

        # Function Description:
        #     For a list of amino acid sequence
        #     discover potential N-linked glycosylation site positions
        # Function Input:
        #     an amino acid sequence
        # Function Output:
        #     a list of potential N-linked glycosylation site positions

        # The rule:
        #     Any amino acid sequence of N-X-S or N-X-T
        #     indicates a potential N-linked glycosylation site at N in the sequence
        #     where X can be any amino acid except P

        # ------ Pass member variable values to local variables ------
        aminoAcidSequenceLength = self.aminoAcidSequenceLength
        aminoAcidSequenceDict = self.aminoAcidSequenceDict
        listOfGlycosylationSitePositions = self.listOfGlycosylationSitePositions

        print("Calculating for a list of potential glycosylation site positions...")
        
        # If it's a sequence of N-X-S or N-X-T, except when X is P,
        #     Add the index of N in this short sequence
        #     to the list of glycosylation site positions
        for i in range(1, aminoAcidSequenceLength - 1):
            if 'N' in aminoAcidSequenceDict[i]:
                if 'P' not in aminoAcidSequenceDict[i + 1]:
                    if ('S' in aminoAcidSequenceDict[i + 2] or
                            'T' in aminoAcidSequenceDict[i + 2]):
                        listOfGlycosylationSitePositions.append(i)

        print("Done.\n")
        
        print("The potential glycosylation site positions (those start at 1) are:")

        # Trivial: Set up a string to display the glycosylation sites
        displayString = ", ".join([str(glycosylationSitePosition)
                                   for glycosylationSitePosition
                                   in listOfGlycosylationSitePositions])

        if len(listOfGlycosylationSitePositions) == 0:
            displayString = "None."
            
        print(displayString + ".\n\n\n")
        
        # ------ Pass local variable values back to member variables ------
        self.listOfGlycosylationSitePositions = listOfGlycosylationSitePositions

        

    def find_potential_phosphorylation_site_positions_from_user_input(self):

        # Function Description:
        #     For a list of amino acid sequence
        #     based on user input of the list of phosphorylation sites to consider
        #     discover potential phosphorylation site positions
        # Function Input:
        #     an amino acid sequence
        #     the list of phosphorylation sites to consider
        # Function Output:
        #     a dictionary of potential phosphorylation site positions

        # ------ Pass member variable values to local variables ------
        listOfPhosphorylationTypes = self.listOfPhosphorylationTypes
        phosphorylationSitePositionDict = self.phosphorylationSitePositionDict
        listOfPhosphorylationSitePositions = self.listOfPhosphorylationSitePositions
        aminoAcidSequenceLength = self.aminoAcidSequenceLength
        aminoAcidSequenceDict = self.aminoAcidSequenceDict
        aminoAcidDict = self.aminoAcidDict

        print("Obtaining a dictionary of potential phosphorylation site positions...")

        for i in range(1, aminoAcidSequenceLength + 1):
            for phosphorylationType in listOfPhosphorylationTypes:

                # Get the amino acid acronym from position i
                aminoAcid = aminoAcidSequenceDict[i]

                # If an amino acid is a potential phosphorylation site
                # append "^phos'ed" to the amino acid acronym
                # so that "S" -> "S^phos'ed", "T" -> "T^phos'ed", "Y" -> "Y^phos'ed", etc.
                # append the index of this position to the list of phosphorylation sites
                # append this position to the phosphorylation site dictionary for printing
                if phosphorylationType in aminoAcid:
                    phosphorylatedAminoAcid = aminoAcid + "^phos'ed"
                    phosphorylatedAminoAcidName = aminoAcidDict[phosphorylatedAminoAcid]["Name"]
                    listOfPhosphorylationSitePositions.append(i)
                    
                    if phosphorylatedAminoAcidName in phosphorylationSitePositionDict.keys():                       
                        phosphorylationSitePositionDict[phosphorylatedAminoAcidName].append(i)
                    else:
                        phosphorylationSitePositionDict[phosphorylatedAminoAcidName] = [i]
                        
        print("Done.\n")

        if listOfPhosphorylationSitePositions != []:        
            print("The potential phosphorylation site positions (those start at 1) are:")

            # Trivial: Set up a string to display the phosphorylation sites
            displayString = ", ".join([phosphorylationSitePosition
                                       for phosphorylationSitePosition
                                       in listOfPhosphorylationSitePositions])
            
            print(displayString + ".\n")

            if listOfPhosphorylationTypes != []:
                
                print("Within which:\n")
            
                for phosphorylationType in phosphorylationSitePositionDict.keys():
                    displayString = ", ".join([phosphorylationSitePosition
                                               for phosphorylationSitePosition
                                               in phosphorylationSitePositionDict[
                                                       phosphorylationType]])
                    
                    print(phosphorylationType + ": " + displayString + ".\n")
                    
                    displayString = ""
                    
                print("\n")                 
            
        # ------ Pass local variable values back to member variables ------
        self.phosphorylationSitePositionDict = phosphorylationSitePositionDict
        self.listOfPhosphorylationSitePositions = listOfPhosphorylationSitePositions
        


    def find_potential_carbamidomethylation_site_positions_from_user_input(self):

        # Function Description:
        #     For a list of amino acid sequence
        #     based on user input of the list of carbamidomethylation sites to consider
        #     discover potential carbamidomethylation site positions
        # Function Input:
        #     an amino acid sequence
        #     the list of carbamidomethylation sites to consider
        # Function Output:
        #     a dictionary of potential carbamidomethylation site positions

        # ------ Pass member variable values to local variables ------
        listOfCarbamidomethylationTypes = self.listOfCarbamidomethylationTypes
        carbamidomethylationSitePositionDict = self.carbamidomethylationSitePositionDict
        listOfCarbamidomethylationSitePositions = self.listOfCarbamidomethylationSitePositions
        aminoAcidSequenceLength = self.aminoAcidSequenceLength
        aminoAcidSequenceDict = self.aminoAcidSequenceDict
        aminoAcidDict = self.aminoAcidDict

        print("Obtaining a dictionary of potential carbamidomethylation site positions...")

        for i in range(1, aminoAcidSequenceLength + 1):
            for carbamidomethylationType in listOfCarbamidomethylationTypes:

                # Get the amino acid acronym from position i
                aminoAcid = aminoAcidSequenceDict[i]

                # If an amino acid is a potential carbamidomethylation site
                # append "^carb'ed" to the amino acid acronym
                # so that "C" -> "C^carb'ed", etc.
                # append the index of this position to the list of carbamidomethylation sites
                # append this position to the carbamidomethylation site dictionary for printing
                if carbamidomethylationType in aminoAcid:
                    carbamidomethylatedAminoAcid = aminoAcid + "^carb'ed"
                    carbamidomethylatedAminoAcidName = aminoAcidDict[carbamidomethylatedAminoAcid]["Name"]
                    listOfCarbamidomethylationSitePositions.append(i)
                    
                    if carbamidomethylatedAminoAcidName in carbamidomethylationSitePositionDict.keys():                       
                        carbamidomethylationSitePositionDict[carbamidomethylatedAminoAcidName].append(i)
                    else:
                        carbamidomethylationSitePositionDict[carbamidomethylatedAminoAcidName] = [i]
                        
        print("Done.\n")

        if listOfCarbamidomethylationSitePositions != []:        
            print("The potential carbamidomethylation site positions (those start at 1) are:")

            # Trivial: Set up a string to display the carbamidomethylation sites
            displayString = ", ".join([str(carbamidomethylationSitePosition)
                                       for carbamidomethylationSitePosition
                                       in listOfCarbamidomethylationSitePositions])
            
            print(displayString + ".\n")

            if listOfCarbamidomethylationTypes != []:
                
                print("Within which:\n")
            
                for carbamidomethylationType in carbamidomethylationSitePositionDict.keys():
                    displayString = ", ".join([str(carbamidomethylationSitePosition)
                                               for carbamidomethylationSitePosition
                                               in carbamidomethylationSitePositionDict[
                                                       carbamidomethylationType]])
                    
                    print(carbamidomethylationType + ": " + displayString + ".\n")
                    
                    displayString = ""
                    
                print("\n")                 
        
        # ------ Pass local variable values back to member variables ------
        self.carbamidomethylationSitePositionDict = carbamidomethylationSitePositionDict
        self.listOfCarbamidomethylationSitePositions = listOfCarbamidomethylationSitePositions
        


    def find_potential_deamidation_site_positions_from_user_input(self):

        # Function Description:
        #     For a list of amino acid sequence
        #     based on user input of the list of deamidation sites to consider
        #     discover potential deamidation site positions
        # Function Input:
        #     an amino acid sequence
        #     the list of deamidation sites to consider
        # Function Output:
        #     a dictionary of potential deamidation site positions

        # ------ Pass member variable values to local variables ------
        listOfDeamidationTypes = self.listOfDeamidationTypes
        deamidationSitePositionDict = self.deamidationSitePositionDict
        listOfDeamidationSitePositions = self.listOfDeamidationSitePositions
        aminoAcidSequenceLength = self.aminoAcidSequenceLength
        aminoAcidSequenceDict = self.aminoAcidSequenceDict
        aminoAcidDict = self.aminoAcidDict

        print("Obtaining a dictionary of potential deamidation site positions...")

        for i in range(1, aminoAcidSequenceLength + 1):
            for deamidationType in listOfDeamidationTypes:

                # Get the amino acid acronym from position i
                aminoAcid = aminoAcidSequenceDict[i]

                # If an amino acid is a potential deamidation site
                # append "^deam'ed" to the amino acid acronym
                # so that "Q" -> "Q^deam'ed", "Q^phos'ed" -> "Q^phos'ed^deam'ed", etc.
                # append the index of this position to the list of deamidation sites
                # append this position to the deamidation site dictionary for printing
                if deamidationType in aminoAcid:
                    deamidatedAminoAcid = aminoAcid + "^deam'ed"
                    deamidatedAminoAcidName = aminoAcidDict[deamidatedAminoAcid]["Name"]
                    listOfDeamidationSitePositions.append(i)
                    
                    if deamidatedAminoAcidName in deamidationSitePositionDict.keys():                       
                        deamidationSitePositionDict[deamidatedAminoAcidName].append(i)
                    else:
                        deamidationSitePositionDict[deamidatedAminoAcidName] = [i]
                        
        print("Done.\n")

        if listOfDeamidationSitePositions != []:        
            print("The potential deamidation site positions (those start at 1) are:")

            # Trivial: Set up a string to display the deamidation sites
            displayString = ", ".join([str(deamidationSitePosition)
                                       for deamidationSitePosition
                                       in listOfDeamidationSitePositions])
            
            print(displayString + ".\n")

            if listOfDeamidationTypes != []:
                
                print("Within which:\n")
            
                for deamidationType in deamidationSitePositionDict.keys():
                    displayString = ", ".join([str(deamidationSitePosition)
                                               for deamidationSitePosition
                                               in deamidationSitePositionDict[
                                                       deamidationType]])
                    
                    print(deamidationType + ": " + displayString + ".\n")
                    
                    displayString = ""
                    
                print("\n")                 
        
        # ------ Pass local variable values back to member variables ------
        self.deamidationSitePositionDict = deamidationSitePositionDict
        self.listOfDeamidationSitePositions = listOfDeamidationSitePositions



    def find_potential_oxidation_site_positions_from_user_input(self):

        # Function Description:
        #     For a list of amino acid sequence
        #     based on user input of the list of oxidation sites to consider
        #     discover potential oxidation site positions
        # Function Input:
        #     an amino acid sequence
        #     the list of oxidation sites to consider
        # Function Output:
        #     a dictionary of potential oxidation site positions

        # ------ Pass member variable values to local variables ------
        listOfOxidationTypes = self.listOfOxidationTypes
        oxidationSitePositionDict = self.oxidationSitePositionDict
        listOfOxidationSitePositions = self.listOfOxidationSitePositions
        aminoAcidSequenceLength = self.aminoAcidSequenceLength
        aminoAcidSequenceDict = self.aminoAcidSequenceDict
        aminoAcidDict = self.aminoAcidDict

        print("Obtaining a dictionary of potential oxidation site positions...")

        for i in range(1, aminoAcidSequenceLength + 1):
            for oxidationType in listOfOxidationTypes:

                # Get the amino acid acronym from position i
                aminoAcid = aminoAcidSequenceDict[i]

                # If an amino acid is a potential deamidation site
                # append "M^oxid'ed" to the amino acid acronym
                # so that "M" -> "M^oxid'ed", "M^phos'ed" -> "M^phos'ed^oxid'ed", etc.
                # append the index of this position to the list of deamidation sites
                # append this position to the deamidation site dictionary for printing
                if oxidationType in aminoAcid:
                    oxydatedAminoAcid = aminoAcid + "^oxid'ed"
                    oxydatedAminoAcidName = aminoAcidDict[oxydatedAminoAcid]["Name"]
                    listOfOxidationSitePositions.append(i)
                    
                    if oxydatedAminoAcidName in oxidationSitePositionDict.keys():                       
                        oxidationSitePositionDict[oxydatedAminoAcidName].append(i)
                    else:
                        oxidationSitePositionDict[oxydatedAminoAcidName] = [i]
                        
        print("Done.\n")

        if listOfOxidationSitePositions != []:        
            print("The potential oxidation site positions (those start at 1) are:")

            # Trivial: Set up a string to display the deamidation sites
            displayString = ", ".join([str(oxidationSitePosition)
                                       for oxidationSitePosition
                                       in listOfOxidationSitePositions])
            
            print(displayString + ".\n")

            if listOfOxidationTypes != []:
                
                print("Within which:\n")
            
                for oxidationType in oxidationSitePositionDict.keys():
                    displayString = ", ".join([str(oxidationSitePosition)
                                               for oxidationSitePosition
                                               in oxidationSitePositionDict[
                                                       oxidationType]])
                    
                    print(oxidationType + ": " + displayString + ".\n")
                    
                    displayString = ""
                    
                print("\n")                 
        
        # ------ Pass local variable values back to member variables ------
        self.oxidationSitePositionDict = oxidationSitePositionDict
        self.listOfOxidationSitePositions = listOfOxidationSitePositions

    

    def trypsin_cleavages(self):

        # ------ Pass member variable values to local variables ------
        aminoAcidSequenceLength = self.aminoAcidSequenceLength
        aminoAcidSequenceDict = self.aminoAcidSequenceDict

        # Data Structure for listOfFragmentLeftRightIndexPairs:
        # listOfFragmentLeftRightIndexPairs = [(1, 21), (22, 71), ...]
        listOfFragmentLeftRightIndexPairs = self.listOfFragmentLeftRightIndexPairs
        
        # Trypsin cleaves peptide chains
        # mainly at the carboxyl (-COOH) side of
        # the amino acids Lysine (K) or Arginine (R),
        # except when either is followed by Proline (P)

        print("Performing trypsin cleavages...")
        
        fragmentLeftIndexIncluded = 1

        # Append all fragment in a format of left-right pair, except for the last fragment
        for i in range(1, aminoAcidSequenceLength):
            if ('K' in aminoAcidSequenceDict[i]
                    or 'R' in aminoAcidSequenceDict[i]):
                try:
                    nextAA = aminoAcidSequenceDict[i + 1]
                except KeyError:
                    nextAA = ""
                if 'P' not in nextAA:
                    listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, i))
                    fragmentLeftIndexIncluded = i + 1

        # Append the last fragment in a format of left-right pair
        listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded,
                                                  aminoAcidSequenceLength))

        print("Done.\n\n\n")

        # ------ Pass local variable values back to member variables ------
        self.listOfFragmentLeftRightIndexPairs = listOfFragmentLeftRightIndexPairs

        

    def lysc_cleavages(self):

        # ------ Pass member variable values to local variables ------
        aminoAcidSequenceLength = self.aminoAcidSequenceLength
        aminoAcidSequenceDict = self.aminoAcidSequenceDict

        # Data Structure for listOfFragmentLeftRightIndexPairs:
        # listOfFragmentLeftRightIndexPairs = [(1, 21), (22, 71), ...]
        listOfFragmentLeftRightIndexPairs = self.listOfFragmentLeftRightIndexPairs
        
        # Trypsin cleaves peptide chains
        # mainly at the carboxyl (-COOH) side of
        # the amino acids Lysine (K)

        print("Performing Lys-C cleavages...")
        
        fragmentLeftIndexIncluded = 1

        # Append all fragment in a format of left-right pair, except for the last fragment
        for i in range(1, aminoAcidSequenceLength):
            if 'K' in aminoAcidSequenceDict[i]:
                listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, i))
                fragmentLeftIndexIncluded = i + 1

        # Append the last fragment in a format of left-right pair
        listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded,
                                                  aminoAcidSequenceLength))

        print("Done.\n\n\n")

        # ------ Pass local variable values back to member variables ------
        self.listOfFragmentLeftRightIndexPairs = listOfFragmentLeftRightIndexPairs



    def gluc_cleavages(self):

        # ------ Pass member variable values to local variables ------
        aminoAcidSequenceLength = self.aminoAcidSequenceLength
        aminoAcidSequenceDict = self.aminoAcidSequenceDict

        # Data Structure for listOfFragmentLeftRightIndexPairs:
        # listOfFragmentLeftRightIndexPairs = [(1, 21), (22, 71), ...]
        listOfFragmentLeftRightIndexPairs = self.listOfFragmentLeftRightIndexPairs
        
        # Trypsin cleaves peptide chains
        # mainly at the carboxyl (-COOH) side of
        # the amino acids glutemate (E)
        
        print("Performing Glu-C cleavages...")
        
        fragmentLeftIndexIncluded = 1

        # Append all fragment in a format of left-right pair, except for the last fragment
        for i in range(1, aminoAcidSequenceLength):
            if 'E' in aminoAcidSequenceDict[i]:
                listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, i))
                fragmentLeftIndexIncluded = i + 1

        # Append the last fragment in a format of left-right pair
        listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded,
                                                 aminoAcidSequenceLength))

        print("Done.\n\n\n")

        # ------ Pass local variable values back to member variables ------
        self.listOfFragmentLeftRightIndexPairs = listOfFragmentLeftRightIndexPairs
    

    def elastase_cleavages(self):
        aminoAcidSequenceLength = self.aminoAcidSequenceLength
        aminoAcidSequenceDict = self.aminoAcidSequenceDict
        # datastructure: listOfFragmentLeftRightIndexPairs = [(1, 21), (22, 71), ...]
        listOfFragmentLeftRightIndexPairs = self.listOfFragmentLeftRightIndexPairs
        # lastase cleaves at C side of A,V,S,G,L,I
        print("Performing elastase cleavages...")
        fragmentLeftIndexIncluded = 1
        for i in range(1, aminoAcidSequenceLength):
            for aa in ['A','V','S','G','L','I']:
                if aa in aminoAcidSequenceDict[i]:
                    listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, i))
                    fragmentLeftIndexIncluded = i + 1
        listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, aminoAcidSequenceLength))
        self.listOfFragmentLeftRightIndexPairs = listOfFragmentLeftRightIndexPairs
        
    def argc_cleavages(self):
        aminoAcidSequenceLength = self.aminoAcidSequenceLength
        aminoAcidSequenceDict = self.aminoAcidSequenceDict
        # datastructure: listOfFragmentLeftRightIndexPairs = [(1, 21), (22, 71), ...]
        listOfFragmentLeftRightIndexPairs = self.listOfFragmentLeftRightIndexPairs
        # lastase cleaves at C side of A,V,S,G,L,I
        print("Performing argc cleavages...")
        fragmentLeftIndexIncluded = 1
        for i in range(1, aminoAcidSequenceLength):
            for aa in ['R', 'K']:
                if aa in aminoAcidSequenceDict[i]:
                    listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, i))
                    fragmentLeftIndexIncluded = i + 1
        listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, aminoAcidSequenceLength))
        self.listOfFragmentLeftRightIndexPairs = listOfFragmentLeftRightIndexPairs
        
    def pepsin_cleavages(self):
        aminoAcidSequenceLength = self.aminoAcidSequenceLength
        aminoAcidSequenceDict = self.aminoAcidSequenceDict
        # datastructure: listOfFragmentLeftRightIndexPairs = [(1, 21), (22, 71), ...]
        listOfFragmentLeftRightIndexPairs = self.listOfFragmentLeftRightIndexPairs
        # pepsin cleaves at C side of FLYW
        print("Performing pepsin cleavages...")
        fragmentLeftIndexIncluded = 1
        for i in range(1, aminoAcidSequenceLength):
            for aa in ['F', 'L', 'Y', 'W']:
                if aa in aminoAcidSequenceDict[i]:
                    listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, i))
                    fragmentLeftIndexIncluded = i + 1
        listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, aminoAcidSequenceLength))
        self.listOfFragmentLeftRightIndexPairs = listOfFragmentLeftRightIndexPairs
        
    def chymotrypsin_cleavages(self):
        aminoAcidSequenceLength = self.aminoAcidSequenceLength
        aminoAcidSequenceDict = self.aminoAcidSequenceDict
        # datastructure: listOfFragmentLeftRightIndexPairs = [(1, 21), (22, 71), ...]
        listOfFragmentLeftRightIndexPairs = self.listOfFragmentLeftRightIndexPairs
        # chymotrypsin = F,Y,W
        print("Performing chymotryp cleavages...")
        fragmentLeftIndexIncluded = 1
        for i in range(1, aminoAcidSequenceLength):
            for aa in ['F','Y','W']:
                if aa in aminoAcidSequenceDict[i]:
                    listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, i))
                    fragmentLeftIndexIncluded = i + 1
        listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, aminoAcidSequenceLength))
        self.listOfFragmentLeftRightIndexPairs = listOfFragmentLeftRightIndexPairs
    
    def chymotrypsin_and_elastase(self):
        aminoAcidSequenceLength = self.aminoAcidSequenceLength
        aminoAcidSequenceDict = self.aminoAcidSequenceDict
        # datastructure: listOfFragmentLeftRightIndexPairs = [(1, 21), (22, 71), ...]
        listOfFragmentLeftRightIndexPairs = self.listOfFragmentLeftRightIndexPairs
        # trp + elastase cleaves at C side of A,V,S,G,L,I,F,Y,W
        print("Performing elastase+chymotryp cleavages...")
        fragmentLeftIndexIncluded = 1
        for i in range(1, aminoAcidSequenceLength):
            for aa in ['A','V','S','G','L','I','F','Y','W']:
                if aa in aminoAcidSequenceDict[i]:
                    listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, i))
                    fragmentLeftIndexIncluded = i + 1
        listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, aminoAcidSequenceLength))
        self.listOfFragmentLeftRightIndexPairs = listOfFragmentLeftRightIndexPairs
    
    def trypsin_and_elastase_and_chymotrypsin(self):
        aminoAcidSequenceLength = self.aminoAcidSequenceLength
        aminoAcidSequenceDict = self.aminoAcidSequenceDict
        # datastructure: listOfFragmentLeftRightIndexPairs = [(1, 21), (22, 71), ...]
        listOfFragmentLeftRightIndexPairs = self.listOfFragmentLeftRightIndexPairs
        # trp + elastase +chymotryp cleaves at C side of A,V,S,G,L,I,F,Y,W & K,R but not before P
        print("Performing trp+elastase+chymotryp cleavages...")
        fragmentLeftIndexIncluded = 1
        for i in range(1, aminoAcidSequenceLength):
            for aa in ['A','V','S','G','L','I','F','Y','W']:
                if aa in aminoAcidSequenceDict[i]:
                    listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, i))
                    fragmentLeftIndexIncluded = i + 1
                elif ('R' in aminoAcidSequenceDict[i] or 'K' in aminoAcidSequenceDict[i]): 
                    try:
                        nextAA = aminoAcidSequenceDict[i + 1]
                    except KeyError:
                        nextAA = ""
                    if 'P' not in nextAA:
                        listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, i))
                        fragmentLeftIndexIncluded = i + 1
        listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, aminoAcidSequenceLength))
        self.listOfFragmentLeftRightIndexPairs = listOfFragmentLeftRightIndexPairs
        
    def trypsin_and_elastase(self):
        aminoAcidSequenceLength = self.aminoAcidSequenceLength
        aminoAcidSequenceDict = self.aminoAcidSequenceDict
        # datastructure: listOfFragmentLeftRightIndexPairs = [(1, 21), (22, 71), ...]
        listOfFragmentLeftRightIndexPairs = self.listOfFragmentLeftRightIndexPairs
        # trp + elastase cleaves at C side of A,V,S,G,L,I, & K,R but not before P
        print("Performing trp+elastase cleavages...")
        fragmentLeftIndexIncluded = 1
        for i in range(1, aminoAcidSequenceLength):
            for aa in ['A','V','S','G','L','I']:
                if aa in aminoAcidSequenceDict[i]:
                    listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, i))
                    fragmentLeftIndexIncluded = i + 1
                elif ('R' in aminoAcidSequenceDict[i] or 'K' in aminoAcidSequenceDict[i]):
                    try:
                        nextAA = aminoAcidSequenceDict[i + 1]
                    except KeyError:
                        nextAA = ""
                    if 'P' not in nextAA:
                        listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, i))
                        fragmentLeftIndexIncluded = i + 1
        listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, aminoAcidSequenceLength))
        self.listOfFragmentLeftRightIndexPairs = listOfFragmentLeftRightIndexPairs
        
    def trypsin_and_lysc_cleavages(self):

        # ------ Pass member variable values to local variables ------
        aminoAcidSequenceLength = self.aminoAcidSequenceLength
        aminoAcidSequenceDict = self.aminoAcidSequenceDict

        # Data Structure for listOfFragmentLeftRightIndexPairs:
        # listOfFragmentLeftRightIndexPairs = [(1, 21), (22, 71), ...]
        listOfFragmentLeftRightIndexPairs = self.listOfFragmentLeftRightIndexPairs
        
        # Trypsin and Lys-C together cleave peptide chains
        # mainly at the carboxyl (-COOH) side of
        # the amino acids Lysine (K) or Arginine (R),
        # except when Arginine (R) is followed by Proline (P)

        print("Performing trypsin and Lys-C cleavages...")
        
        fragmentLeftIndexIncluded = 1

        # Append all fragment in a format of left-right pair, except for the last fragment
        for i in range(1, aminoAcidSequenceLength):
            if ('R' in aminoAcidSequenceDict[i] or 'K' in aminoAcidSequenceDict[i]):
                try:
                    nextAA = aminoAcidSequenceDict[i + 1]
                except KeyError:
                    nextAA = ""
                if 'P' not in aminoAcidSequenceDict[i + 1]:
                    listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, i))
                    fragmentLeftIndexIncluded = i + 1

        # Append the last fragment in a format of left-right pair
        listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded,
                                                 aminoAcidSequenceLength))

        print("Done.\n\n\n")

        # ------ Pass local variable values back to member variables ------
        self.listOfFragmentLeftRightIndexPairs = listOfFragmentLeftRightIndexPairs



    def trypsin_and_gluc_cleavages(self):

        # ------ Pass member variable values to local variables ------
        aminoAcidSequenceLength = self.aminoAcidSequenceLength
        aminoAcidSequenceDict = self.aminoAcidSequenceDict

        # Data Structure for listOfFragmentLeftRightIndexPairs:
        # listOfFragmentLeftRightIndexPairs = [(1, 21), (22, 71), ...]
        listOfFragmentLeftRightIndexPairs = self.listOfFragmentLeftRightIndexPairs
        
        # Trypsin and Glu-C together cleave peptide chains
        # mainly at the carboxyl (-COOH) side of
        # the amino acids Lysine (K) or Arginine (R) or Glutamic acid (E),
        # except when Lysine (K) or Arginine (R) is followed by Proline (P)
        
        print("Performing trypsin and Glu-C cleavages...")
        
        fragmentLeftIndexIncluded = 1

        # Append all fragment in a format of left-right pair, except for the last fragment
        for i in range(1, aminoAcidSequenceLength):
            if 'E' in aminoAcidSequenceDict[i]:
                listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, i))
                fragmentLeftIndexIncluded = i + 1
            elif ('K' in aminoAcidSequenceDict[i]
                    or 'R' in aminoAcidsequenceDict[i]):
                if 'P' not in aminoAcidSequenceDict[i + 1]:
                    listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, i))
                    fragmentLeftIndexIncluded = i + 1

        # Append the last fragment in a format of left-right pair
        listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded,
                                                 aminoAcidSequenceLength))

        print("Done.\n\n\n")

        # ------ Pass local variable values back to member variables ------
        self.listOfFragmentLeftRightIndexPairs = listOfFragmentLeftRightIndexPairs



    def lysc_and_gluc_cleavages(self):

        # ------ Pass member variable values to local variables ------
        aminoAcidSequenceLength = self.aminoAcidSequenceLength
        aminoAcidSequenceDict = self.aminoAcidSequenceDict

        # Data Structure for listOfFragmentLeftRightIndexPairs:
        # listOfFragmentLeftRightIndexPairs = [(1, 21), (22, 71), ...]
        listOfFragmentLeftRightIndexPairs = self.listOfFragmentLeftRightIndexPairs
        
        # Lys-C and Glu-C together cleave peptide chains
        # mainly at the carboxyl (-COOH) side of
        # the amino acids Lysine (K) or Glutamic acid (E)

        print("Performing Lys-C and Glu-C cleavages...")
        
        fragmentLeftIndexIncluded = 1

        # Append all fragment in a format of left-right pair, except for the last fragment
        for i in range(1, aminoAcidSequenceLength):
            if ('K' in aminoAcidSequenceDict[i]
                    or 'E' in aminoAcidSequenceDict[i]):
                listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, i))
                fragmentLeftIndexIncluded = i + 1

        # Append the last fragment in a format of left-right pair
        listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded,
                                                 aminoAcidSequenceLength))

        print("Done.\n\n\n")

        # ------ Pass local variable values back to member variables ------
        self.listOfFragmentLeftRightIndexPairs = listOfFragmentLeftRightIndexPairs



    def trysin_and_lysc_and_gluc_cleavages(self):

        # ------ Pass member variable values to local variables ------
        aminoAcidSequenceLength = self.aminoAcidSequenceLength
        aminoAcidSequenceDict = self.aminoAcidSequenceDict

        # Data Structure for listOfFragmentLeftRightIndexPairs:
        # listOfFragmentLeftRightIndexPairs = [(1, 21), (22, 71), ...]
        listOfFragmentLeftRightIndexPairs = self.listOfFragmentLeftRightIndexPairs
        
        # Trypsin and Lys-C and Glu-C together cleave peptide chains
        # mainly at the carboxyl (-COOH) side of
        # the amino acids Lysine (K) or Arginine (R) or Glutamic Acid (E)
        # except when Arginine (R) is followed by Proline (P)

        print("Performing trypsin and Lys-C and Glu-C cleavages...")
        
        fragmentLeftIndexIncluded = 1

        # Append all fragment in a format of left-right pair, except for the last fragment
        for i in range(1, aminoAcidSequenceLength):
            if ('K' in aminoAcidSequenceDict[i] or
                        'E' in aminoAcidSequenceDict[i]):
                listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, i))
                fragmentLeftIndexIncluded = i + 1
            elif 'R' in aminoAcidSequenceDict[i]:
                if 'P' not in aminoAcidSequenceDict[i + 1]:
                    listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded, i))
                    fragmentLeftIndexIncluded = i + 1
                    

        # Append the last fragment in a format of left-right pair
        listOfFragmentLeftRightIndexPairs.append((fragmentLeftIndexIncluded,
                                                 aminoAcidSequenceLength))

        print("Done.\n\n\n")

        # ------ Pass local variable values back to member variables ------
        self.listOfFragmentLeftRightIndexPairs = listOfFragmentLeftRightIndexPairs



    def find_all_mass_matches(self, featureObject):

        # Function Description:
        #     For all glycosylation sites, from interval length 1 to max_length (10 now),
        #     check all intervals that contains the site
        #     Find the matches that has a tentative glycopeptide mass
        #     be within the range of mass spectrometer input +- the deviation given
        #     Return a list of massMatch objects
        # Function Output:
        #     massMatch objects for all glycosylation sites

        
        inputMassToMatch = featureObject.compoundNeutralMass
        self.inputMassToMatch = inputMassToMatch 
        # ------ Pass function input variable values to member variables ------
        decisionDeviation = self.config["tolerance"]
        hydrogenMass = self.hydrogenMass
        
        print("Looking for mass matches between a glycopeptide and" +
              " the mass spectrometer output...\n")

        if hasattr(featureObject, "massMatchList") == False:
            massMatchList = []
        else:
            massMatchList = featureObject.massMatchList
        
        # set the upper and lower bounds for a match
        lowerMassMatchLimit = inputMassToMatch / (1 + decisionDeviation)
        upperMassMatchLimit = inputMassToMatch / (1 - decisionDeviation)
        #print(str(lowerMassMatchLimit) + str(upperMassMatchLimit) + "/n")
        
        # If the reduced attribute is NOT checked,
        #     we consider the normal glycopeptide mass matching
        if not self.config['free']:
            
            # ------ Find the mass matches for every glycosylation site ------
            for sitePosition in self.listOfGlycosylationSitePositions:

                # ------ Get an iterator from the generator for one glycosylation site ------
                # ------ the iterator is used to record all massMatches for one glycosylation site ------
                iterator_massMatches = self.generator__find_mass_matches_for_one_glycosylation_site(featureObject,
                                                                                                    sitePosition)

                massMatches = list(iterator_massMatches)
                # Append the list of matches for one glycosylation site to
                #     the list of all detected mass matches
                #     which in the end records all mass matches for all glycosylation sites
                massMatchList.extend(massMatches)
            
            featureObject.massMatchList = massMatchList
               
            print("Done.\n\n\n")

            print(self.proteinName + ":")
            print("There are " + str(self.globalCounterForMatches) + " matches found!\n\n\n")

        # If the reduced attribute is checked,
        #     we only consider the reduced glycan mass matching without the peptide chain
        else:

            for glycanAlone in self.generator__glycan_maker():

                glycanMass = self.get_glycan_mass(glycanAlone)
                reducedGlycanMass = glycanMass + 2 * self.hydrogenMass
                reducedGlycanMass = round(reducedGlycanMass, 10)
                
                isNewMatch = False
                isDecoyUsed = False
                
                # If we do NOT apply decoy mass
                #     If the reduced glycan mass is within input mass from Mass Spectrometer
                #             +- the maximum deviation allowed
                #         One match is found
                # If we apply decoy mass
                #     If the reduced glycan mass is within input mass from Mass Spectrometer
                #             +- the maximum deviation allowed
                #         One match is found
                #     If the reduced glycan mass plus decoy mass is within input mass from Mass Spectrometer
                #             +- the maximum deviation allowed
                #         One match is found
                if self.config["decoyAnalysis"] == False:
                    isNewMatch = (reducedGlycanMass >= lowerMassMatchLimit and
                            reducedGlycanMass <= upperMassMatchLimit)
                    decoyMassToMakeMatch = None
                elif self.config["decoyAnalysis"] == True:
                    if (reducedGlycanMass >= lowerMassMatchLimit and
                            reducedGlycanMass <= upperMassMatchLimit):
                        isNewMatch = True
                        decoyMassToMakeMatch = None
                    elif (reducedGlycanMass + self.config["decoyMass"] >= lowerMassMatchLimit and
                            reducedGlycanMass + self.config["decoyMass"] <= upperMassMatchLimit):
                        isNewMatch = True
                        decoyMassToMakeMatch = self.config["decoyMass"]
                        isDecoyUsed = True                                                        
                """
                print("isNewMatch: " + str(isNewMatch))
                print("reducedGlycanMass: " + str(reducedGlycanMass))
                print("lowerMassMatchLimit: " + str(lowerMassMatchLimit))
                print("upperMassMatchLimit: " + str(upperMassMatchLimit))
                """
                if isNewMatch:
                    print( "works!!!")
                    # Increment the counter by one for result printing purposes
                    self.globalCounterForMatches = self.globalCounterForMatches + 1

                    matchIndex = self.globalCounterForMatches

                    peptideStringMatched = None
                    peptideStringLengthMatched = None

                    sitePosition = None

                    leftIndexIncluded = None
                    rightIndexIncluded = None
                    
                    glycosylationSiteAminoAcid = None
                    
                    glycosylationSiteAminoAcidMass = None

                    peptideLeftSideAnimoAcid = None

                    peptideLeftSideAnimoAcidMass = None

                    peptideRightSideAnimoAcid = None

                    peptideRightSideAnimoAcidMass = None

                    peptideMass = None

                    glycopeptideMass = None

                    totalTheoreticalMass = reducedGlycanMass
                    
                    '''# Get a length-5 tuple, in the order of
                    # (Hex, HexNac, DxyHex, NeuAc, Sulf)'''
                    glycanCodeMatched = glycanAlone["Glycan Code"]
                    '''Hexose = glycanAlone["Hex"]
                    HexNac = glycanAlone["HexNAc"]
                    DeoxyHex = glycanAlone["DeoxyHex"]
                    Neu5Ac = glycanAlone["Neu5Ac"]
                    Sulf = 0
                    glycanComponentsLengthFiveTuple = (Hexose,
                                                       HexNac,
                                                       DeoxyHex,
                                                       Neu5Ac,
                                                       Sulf)'''

                    # Calculate the actual deviaßtion
                    # from the tentative mass to the mass spectrometer input mass
                    # in a unit of ppm (parts per million)
                    if self.config["decoyAnalysis"] == False:
                        
                        deviationFromActualToTheoreticalMass = ((inputMassToMatch - reducedGlycanMass)/reducedGlycanMass) * 1000000
                        deviationFromActualToTheoreticalMass = str(round(deviationFromActualToTheoreticalMass, 2)) + " ppm"

                        deviationFromActualToDecoyedTheoreticalMass = None

                    elif self.config["decoyAnalysis"] == True:

                        deviationFromActualToTheoreticalMass = ((inputMassToMatch - reducedGlycanMass)/reducedGlycanMass) * 1000000
                        deviationFromActualToTheoreticalMass = str(round(deviationFromActualToTheoreticalMass, 2)) + " ppm"               
                        
                        deviationFromActualToDecoyedTheoreticalMass = ((inputMassToMatch - reducedGlycanMass - self.config["decoyMass"])/reducedGlycanMass) * 1000000
                        deviationFromActualToDecoyedTheoreticalMass = str(round(deviationFromActualToDecoyedTheoreticalMass, 2)) + " ppm"        


                    print("\n******************************\n")

                    print("Glycan only:")
                    print("Match " + str(self.globalCounterForMatches) +
                          " found!! With Glycan Code " + glycanCodeMatched + ".\n\n\n")

                    # Yield a massMatch object
                    newMassMatch = massMatch(matchIndex,
                                             self.proteinFileName,
                                             self.proteinName,
                                             peptideStringMatched,
                                             peptideStringLengthMatched,
                                             glycanCodeMatched,
                                             glycanAlone,
                                             sitePosition,
                                             glycosylationSiteAminoAcid,
                                             glycosylationSiteAminoAcidMass,
                                             leftIndexIncluded,
                                             peptideLeftSideAnimoAcid,
                                             peptideLeftSideAnimoAcidMass,
                                             rightIndexIncluded,
                                             peptideRightSideAnimoAcid,
                                             peptideRightSideAnimoAcidMass,
                                             peptideMass,
                                             glycanMass,
                                             reducedGlycanMass,
                                             glycopeptideMass,
                                             totalTheoreticalMass,
                                             isDecoyUsed,
                                             decoyMassToMakeMatch,
                                             inputMassToMatch,
                                             deviationFromActualToTheoreticalMass,
                                             deviationFromActualToDecoyedTheoreticalMass,
                                             featureObject.compoundMzValue,
                                             featureObject.compoundIntensity,
                                             featureObject.compoundRtCenter,
                                             hydrogenMass)                                                     

                    massMatchList.append(newMassMatch)

                isNewMatch = False
                
            featureObject.massMatchList = massMatchList

                
            print("Done.\n\n\n")

            print("Glycan only:")
            print("There are " + str(self.globalCounterForMatches) + " matches found!\n\n\n")

        '''# ------ Print out the details of all mass matchings ------
        for massMatch in allMassMatches:
            massMatch.printMatchDetails()'''

        '''# ------ Pass local variable values back to member variables ------
        self.allMassMatches = massMatchList'''

                
        
    def generator__find_mass_matches_for_one_glycosylation_site(self,
                                                                featureObject,
                                                                sitePosition):

        # Function Description:
        #     Find all mass matches for one glycosylation site
        #     Yield all massMatch objects one by one
        # Function Input:
        #     One glycosylation site position
        # Function Output:
        #     massMatch objects for one glycosylation site
        
        print("Looking for mass matches for Glycosylation Site at Position " + str(sitePosition) + "...")

        # ------ Pass member variable values to local variables ------
        aminoAcidDict = self.aminoAcidDict
        aminoAcidSequenceDict = self.aminoAcidSequenceDict
        decisionDeviation = self.config["tolerance"]
        inputMassToMatch = self.inputMassToMatch
        glycanDict = self.glycanDict
        hydrogenMass = self.hydrogenMass

        '''print("".join([aminoAcidSequenceDict[index]
              for index
              in range(1, len(aminoAcidSequenceDict.keys()) + 1)]))'''
        # ------ Get an interval iterator from the interval generator ------
        # ------ Each interval has minimum length 1 and maximum length max_length (10 given now) ------
        iterator_peptideIntervals = self.generator__find_peptide_intervals(sitePosition)

        # define mass intervals for finding match
        lowerMassMatchLimit = inputMassToMatch * (1 - decisionDeviation)
        upperMassMatchLimit = inputMassToMatch * (1 + decisionDeviation)
        
        print("Done.\n")

        # ------ Find all mass matches for every interval and every glycosylation site ------

        # For all intervals that contains one glycosylation site
        for leftRightPairs in iterator_peptideIntervals:
            # For all glycan masses in the library
            for glycanAttached in self.generator__glycan_maker():
                
                decoyStatusList = [0, 1] if self.config["decoyAnalysis"] else [0]
                for decoyStatus in decoyStatusList:
                    # leftRightPairs = [leftIndexIncluded, rightIndexIncluded]
                    leftIndexIncluded = leftRightPairs[0]
                    rightIndexIncluded = leftRightPairs[1]
                    '''print("".join([aminoAcidSequenceDict[index]
                                   for index
                                   in range(leftIndexIncluded, rightIndexIncluded + 1)]))'''
    ##              print("Checking sequence from " + str(leftIndexIncluded) + " to " +
    ##                             str(rightIndexIncluded) + " for glycosylation site " + str(sitePosition) + "...")
                    glycanMassAttached = glycanAttached["Molecular Mass"]
                    # Calculate the mass for one tentative glycopeptide
                    massTriple = self.get_peptide_and_glycopeptide_molecular_mass(
                                                                leftIndexIncluded,
                                                                rightIndexIncluded,
                                                                glycanMassAttached)
                    # massTriple = [peptideMass (minus one water), glycanMass, glycopeptideMass]
                    glycopeptideMass = massTriple[2]
                    #print(glycopeptideMass)
                    
                    isNewMatch = False
                    isDecoyUsed = False
                    
                    #print(self.config["decoyAnalysis"])
                    #print(lowerMassMatchLimit)
                    #print(upperMassMatchLimit)
                    # If we do NOT apply decoy mass
                    #     If the glycopeptide mass is within input mass from Mass Spectrometer
                    #             +- the maximum deviation allowed
                    #         One match is found
                    # If we apply decoy mass
                    #     If the glycopeptide mass is within input mass from Mass Spectrometer
                    #             +- the maximum deviation allowed
                    #         One match is found
                    #     If the glycopeptide mass plus decoy mass is within input mass from Mass Spectrometer
                    #             +- the maximum deviation allowed
                    #         One match is found
                    if decoyStatus == 0:
                        if (glycopeptideMass >= lowerMassMatchLimit and
                                glycopeptideMass <= upperMassMatchLimit):
                            isNewMatch = True
                            decoyMassToMakeMatch = None
                    else:
                        if (glycopeptideMass + self.config["decoyMass"] >= lowerMassMatchLimit and
                                glycopeptideMass + self.config["decoyMass"] <= upperMassMatchLimit):
                            isNewMatch = True
                            decoyMassToMakeMatch = self.config["decoyMass"]                                                         
                            isDecoyUsed = True

                    """
                    print("isNewMatch: " + str(isNewMatch))
                    print("glycopeptideMass: " + str(glycopeptideMass))
                    print("decoyMass: " + str(self.config["decoyMass"]))
                    print("lowerMassMatchLimit: " + str(lowerMassMatchLimit))
                    print("upperMassMatchLimit: " + str(upperMassMatchLimit))
                    print("")
                    """
                    if isNewMatch:
                        print("match found")
                        # Increment the counter by one for result printing purposes
                        self.globalCounterForMatches = self.globalCounterForMatches + 1

                        # The counter value is the index for the current match
                        #     which starts at 1
                        matchIndex = self.globalCounterForMatches

                        # Calculate the peptide string

                        
                        # We store the string in a format of N-R-S-S-R-T isntead of NRSSRT
                        # So that if it is N-R-S-PS-R-PT
                        # where PS and PT are phosphorylation sites
                        # It can still be displayed correctly with ambiguity
                        # And we can make the string back by
                        '''
                        # sequenceList = "N-R-S-S-R-T".split("-")
                        # sequence = ""
                        # for aminoAcid in sequenceList:
                        #     sequence = sequence + aminoAcid
                        '''
                        
                        peptideStringMatched = "-".join([aminoAcidSequenceDict[index]
                                                         for index
                                                         in range(leftIndexIncluded,
                                                                  rightIndexIncluded + 1)])
                        
                        peptideStringLengthMatched = rightIndexIncluded - leftIndexIncluded + 1
                        
                        # Get the amino acid and amino acid mass at the glycosylation site
                        glycosylationSiteAminoAcid = aminoAcidSequenceDict[sitePosition]

                        aminoAcidMass = aminoAcidDict[glycosylationSiteAminoAcid]["Mass"]
                        glycosylationSiteAminoAcidMass = aminoAcidMass

                        # Get the amino acid and amino acid mass at the left end on the matched peptide
                        peptideLeftSideAnimoAcid = aminoAcidSequenceDict[leftIndexIncluded]

                        aminoAcidMass = aminoAcidDict[peptideLeftSideAnimoAcid]["Mass"]
                        peptideLeftSideAnimoAcidMass = aminoAcidMass
                                 
                        # Get the amino acid and amino acid mass at the right end on the matched peptide
                        peptideRightSideAnimoAcid = aminoAcidSequenceDict[rightIndexIncluded]

                        aminoAcidMass = aminoAcidDict[peptideRightSideAnimoAcid]["Mass"]
                        peptideRightSideAnimoAcidMass = aminoAcidMass

                        # massTriple = [peptideMass (minus one water), glycanMass, glycopeptideMass]
                        peptideMass = massTriple[0]
                        glycanMass = massTriple[1]

                        totalTheoreticalMass = glycopeptideMass
                        '''# Get a length-5 tuple, in the order of
                        # (Hex, HexNac, DxyHex, NeuAc, Sulf)'''
                        glycanCodeMatched = glycanAttached["Glycan Code"]
                        '''Hexose = glycanAttached["Hex"]
                        HexNac = glycanAttached["HexNAc"]
                        DeoxyHex = glycanAttached["DeoxyHex"]
                        Neu5Ac = glycanAttached["Neu5Ac"]
                        Sulf = 0
                        glycanComponentsLengthFiveTuple = (Hexose,
                                                           HexNac,
                                                           DeoxyHex,
                                                           Neu5Ac,
                                                           Sulf)'''

                        # Calculate the actual devition
                        # from the tentative mass to the mass spectrometer input mass
                        # in a unit of ppm (parts per million)
                        if self.config["decoyAnalysis"] == False:
                            
                            deviationFromActualToTheoreticalMass = ((inputMassToMatch - glycopeptideMass)/glycopeptideMass) * 1000000
                            deviationFromActualToTheoreticalMass = str(round(deviationFromActualToTheoreticalMass, 2)) + " ppm"                                                            

                            deviationFromActualToDecoyedTheoreticalMass = None
                            
                        elif self.config["decoyAnalysis"] == True:
                            
                            deviationFromActualToTheoreticalMass = ((inputMassToMatch - glycopeptideMass)/glycopeptideMass) * 1000000
                            deviationFromActualToTheoreticalMass = str(round(deviationFromActualToTheoreticalMass, 2)) + " ppm"      

                            deviationFromActualToDecoyedTheoreticalMass = ((inputMassToMatch - self.config["decoyMass"] - glycopeptideMass)/glycopeptideMass) * 1000000
                            deviationFromActualToDecoyedTheoreticalMass = str(round(deviationFromActualToDecoyedTheoreticalMass, 2)) + " ppm"
                            
                        print("\n******************************\n")
                        
                        print(self.proteinName + ":")
                        print("Match " + str(self.globalCounterForMatches) +
                              " found!! From Position " +
                              str(leftIndexIncluded) +
                              " inclusive to Position " +
                              str(rightIndexIncluded) +
                              " inclusive on Glycosylation Site at Position " +
                              str(sitePosition) +
                              " with Glycan Code " + glycanCodeMatched + ".\n\n\n")

                        reducedGlycanMass = glycanMass + 2 * hydrogenMass

                        
                        # Yield a massMatch object
                        yield massMatch(matchIndex,
                                        self.proteinFileName,
                                        self.proteinName,
                                        peptideStringMatched,
                                        peptideStringLengthMatched,
                                        glycanCodeMatched,
                                        glycanAttached,
                                        sitePosition,
                                        glycosylationSiteAminoAcid,
                                        glycosylationSiteAminoAcidMass,
                                        leftIndexIncluded,
                                        peptideLeftSideAnimoAcid,
                                        peptideLeftSideAnimoAcidMass,
                                        rightIndexIncluded,
                                        peptideRightSideAnimoAcid,
                                        peptideRightSideAnimoAcidMass,
                                        peptideMass,
                                        glycanMass,
                                        reducedGlycanMass, 
                                        glycopeptideMass,
                                        totalTheoreticalMass,
                                        isDecoyUsed,
                                        decoyMassToMakeMatch,
                                        inputMassToMatch,
                                        deviationFromActualToTheoreticalMass,
                                        deviationFromActualToDecoyedTheoreticalMass,
                                        featureObject.compoundMzValue,
                                        featureObject.compoundIntensity,
                                        featureObject.compoundRtCenter,
                                        hydrogenMass)                                                     
                


    def generator__glycan_maker(self):

        # Function Description:
        #     Return the glycans
        #     from either Scott's Library or combinatorial numbers of glycan components
        # Function Input:
        #     Scott's Library of combinatorial min and max numbers of glycan components
        # Function Output:
        #     One glycan, with all detailed information of the glycan included
        
        glycanDict = self.glycanDict
        methodToGetGlycanData = self.config["glycanSearchMethod"] 
        
        if methodToGetGlycanData == "scottslibrary" or methodToGetGlycanData == "genericlibrary":
            
            for key in glycanDict.keys():
                yield glycanDict[key]

        elif methodToGetGlycanData == "combinatorial":
            
            glycanComponentMinMaxDict = self.config["glycanComponentMinMaxDict"]
            
            HexoseMin = glycanComponentMinMaxDict["Hex"]["min"]
            HexoseMax = glycanComponentMinMaxDict["Hex"]["max"]

            HexNAcMin = glycanComponentMinMaxDict["HexNAc"]["min"]
            HexNAcMax = glycanComponentMinMaxDict["HexNAc"]["max"]

            DeoxyHexMin = glycanComponentMinMaxDict["DeoxyHex"]["min"]
            DeoxyHexMax = glycanComponentMinMaxDict["DeoxyHex"]["max"]

            Neu5AcMin = glycanComponentMinMaxDict["Neu5Ac"]["min"]
            Neu5AcMax = glycanComponentMinMaxDict["Neu5Ac"]["max"]

            Neu5GcMin = glycanComponentMinMaxDict["Neu5Gc"]["min"]
            Neu5GcMax = glycanComponentMinMaxDict["Neu5Gc"]["max"]

            PentoseMin = glycanComponentMinMaxDict["Pentose"]["min"]
            PentoseMax = glycanComponentMinMaxDict["Pentose"]["max"]

            KDNMin = glycanComponentMinMaxDict["KDN"]["min"]
            KDNMax = glycanComponentMinMaxDict["KDN"]["max"]

            HexAMin = glycanComponentMinMaxDict["HexA"]["min"]
            HexAMax = glycanComponentMinMaxDict["HexA"]["max"]

            SulfMin = 0
            SulfMax = 0
            
            listOfCircumstances = [(Hexose, HexNAc, DeoxyHex, Neu5Ac,
                                    Neu5Gc, Pentose, KDN, HexA, Sulf)
                                  for Hexose in range(HexoseMin, HexoseMax + 1)
                                  for HexNAc in range(HexNAcMin, HexNAcMax + 1)
                                  for DeoxyHex in range(DeoxyHexMin, DeoxyHexMax + 1)
                                  for Neu5Ac in range(Neu5AcMin, Neu5AcMax + 1)
                                  for Neu5Gc in range(Neu5GcMin, Neu5GcMax + 1)
                                  for Pentose in range(PentoseMin, PentoseMax + 1)
                                  for KDN in range(KDNMin, KDNMax + 1)
                                  for HexA in range(HexAMin, HexAMax + 1)
                                  for Sulf in range(SulfMin, SulfMax + 1)]
                                   
            for circumstance in listOfCircumstances:
                glycan = {}
                
                [glycan["Hex"], 
                 glycan["HexNAc"], 
                 glycan["DeoxyHex"], 
                 glycan["Neu5Ac"], 
                 glycan["Neu5Gc"], 
                 glycan["Pentose"], 
                 glycan["KDN"],
                 glycan["HexA"],
                 glycan["Sulf"]] = circumstance
                                    
                #glycan["decoy"] = decoy     # mass = 11.00000
                glycan["Carbohydrate Type"] = ""
                glycan["Glycan Code"] = (str(glycan["Hex"]) +
                                         str(glycan["HexNAc"]) +
                                         str(glycan["DeoxyHex"]) +
                                         str(glycan["Neu5Ac"]) +
                                         str(glycan["Sulf"]))
                glycan["Class"] = ""
                
                glycan["Molecular Mass"] = self.get_glycan_mass(glycan)

                yield glycan



    def generator__find_peptide_intervals(self,
                                          sitePosition):

        # Function Description:
        #     Generates and returns all peptide intervals,
        #     based on all intervals with a max length or trypsin cleavages
        #     and modify the amino acids into phosphorylated ones if necessary
        # Function Input:
        #     One glycosylation site position
        # Function Output:
        #     The left index included on the peptide chain
        #     The right index included on the peptide chain
        #     Modifications from amino acids to phosphorylated ones and vice versa
        
        # ------ Pass member variable values to local variables ------
        aminoAcidSequenceLength = self.aminoAcidSequenceLength
        
        listOfPhosphorylationSitePositions = self.listOfPhosphorylationSitePositions
        listOfCarbamidomethylationSitePositions = self.listOfCarbamidomethylationSitePositions
        listOfDeamidationSitePositions = self.listOfDeamidationSitePositions
        listOfOxidationSitePositions = self.listOfOxidationSitePositions

        '''print(listOfPhosphorylationSitePositions)
        print(listOfCarbamidomethylationSitePositions)
        print(listOfDeamidationSitePositions)
        print(listOfOxidationSitePositions)'''
        
        listOfFragmentLeftRightIndexPairs = self.listOfFragmentLeftRightIndexPairs
        digestTypes = self.config["digestTypes"]
        isSequenceGivenByUser = self.config["isSequenceGivenByUser"]
                    
        if (len(digestTypes) == 1
                and "nonspecific" in digestTypes
                and isSequenceGivenByUser == False):

            maxPeptideLengthInMatch = int(self.config["missedCleavages"])
            
            if len(listOfCarbamidomethylationSitePositions) == 0:

                if len(listOfDeamidationSitePositions) == 0:

                    if len(listOfOxidationSitePositions) == 0:
                        
                        # Interval enumeration algorithm:
                        # Assume the glycopeptide is at position 20
                        # Then the intervals are (both sides included, in the order enumerated below)
                        # [11, 20],
                        # [12, 20], [12, 21],
                        # [13, 20], [13, 21], [13, 22],
                        # [14, 20], [14, 21], [14, 22], [14, 23],
                        # [15, 20], [15, 21], .................., [15, 24]
                        # [16, 20], [16, 21], ....................., [16, 25]
                        # [17, 20], [17, 21], ........................, [17, 26]
                        # [18, 20], [18, 21], ..........................., [18, 27]
                        # [19, 20], [19, 21], .............................., [19, 28]
                        # [20, 20], [20, 21], ................................., [20, 29]
                        for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                       sitePosition + 1):
                            for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                               maxPeptideLengthInMatch + 1):
                                rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                # If the peptide left-right indices are NOT out of bound,
                                #     yield the left-right index pair
                                if (leftIndexIncluded > 0 and
                                        rightIndexIncluded <= aminoAcidSequenceLength):
                                        
                                    yield (leftIndexIncluded, rightIndexIncluded)

                    else:

                        for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                       sitePosition + 1):
                            for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                               maxPeptideLengthInMatch + 1):
                                rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                # If the peptide left-right indices are NOT out of bound,
                                #     yield the left-right index pair
                                if (leftIndexIncluded > 0 and
                                        rightIndexIncluded <= aminoAcidSequenceLength):

                                    for situation in generator__binary_select_oxidation_sites(
                                            leftIndexIncluded, rightIndexIncluded):
                                        
                                        yield (leftIndexIncluded, rightIndexIncluded)                        

                else:

                    if len(listOfOxidationSitePositions) == 0:
                        
                        for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                       sitePosition + 1):
                            for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                               maxPeptideLengthInMatch + 1):
                                rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                # If the peptide left-right indices are NOT out of bound,
                                #     yield the left-right index pair
                                if (leftIndexIncluded > 0 and
                                        rightIndexIncluded <= aminoAcidSequenceLength):

                                    for situation in generator__binary_select_deamidation_sites(
                                            leftIndexIncluded, rightIndexIncluded):

                                        yield (leftIndexIncluded, rightIndexIncluded)

                    else:

                        for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                       sitePosition + 1):
                            for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                               maxPeptideLengthInMatch + 1):
                                rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                # If the peptide left-right indices are NOT out of bound,
                                #     yield the left-right index pair
                                if (leftIndexIncluded > 0 and
                                        rightIndexIncluded <= aminoAcidSequenceLength):

                                    for situation in generator__binary_select_deamidation_sites(
                                            leftIndexIncluded, rightIndexIncluded):

                                        for situation in generator__binary_select_oxidation_sites(
                                                leftIndexIncluded, rightIndexIncluded):
                                            
                                            yield (leftIndexIncluded, rightIndexIncluded)

            else:

                if len(listOfDeamidationSitePositions) == 0:

                    if len(listOfOxidationSitePositions) == 0:
                        
                        for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                       sitePosition + 1):
                            for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                               maxPeptideLengthInMatch + 1):
                                rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                # If the peptide left-right indices are NOT out of bound,
                                #     yield the left-right index pair
                                if (leftIndexIncluded > 0 and
                                        rightIndexIncluded <= aminoAcidSequenceLength):

                                    for situation in self.generator__binary_select_carbamidomethylation_sites(
                                            leftIndexIncluded, rightIndexIncluded):
                                        
                                        yield (leftIndexIncluded, rightIndexIncluded)

                    else:

                        for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                       sitePosition + 1):
                            for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                               maxPeptideLengthInMatch + 1):
                                rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                # If the peptide left-right indices are NOT out of bound,
                                #     yield the left-right index pair
                                if (leftIndexIncluded > 0 and
                                        rightIndexIncluded <= aminoAcidSequenceLength):

                                    for situation in self.generator__binary_select_carbamidomethylation_sites(
                                            leftIndexIncluded, rightIndexIncluded):

                                        for situation in self.generator__binary_select_oxidation_sites(
                                                leftIndexIncluded, rightIndexIncluded):
                                            
                                            yield (leftIndexIncluded, rightIndexIncluded)

                else:

                    if len(listOfOxidationSitePositions) == 0:
                        
                        for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                       sitePosition + 1):
                            for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                               maxPeptideLengthInMatch + 1):
                                rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                # If the peptide left-right indices are NOT out of bound,
                                #     yield the left-right index pair
                                if (leftIndexIncluded > 0 and
                                        rightIndexIncluded <= aminoAcidSequenceLength):

                                    for situation in self.generator__binary_select_carbamidomethylation_sites(
                                            leftIndexIncluded, rightIndexIncluded):

                                        for situation in self.generator__binary_select_deamidation_sites(
                                                leftIndexIncluded, rightIndexIncluded):
                                            
                                            yield (leftIndexIncluded, rightIndexIncluded)

                    else:

                        for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                       sitePosition + 1):
                            for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                               maxPeptideLengthInMatch + 1):
                                rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                # If the peptide left-right indices are NOT out of bound,
                                #     yield the left-right index pair
                                if (leftIndexIncluded > 0 and
                                        rightIndexIncluded <= aminoAcidSequenceLength):

                                    for situation in self.generator__binary_select_carbamidomethylation_sites(
                                            leftIndexIncluded, rightIndexIncluded):

                                        for situation in self.generator__binary_select_deamidation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in self.generator__binary_select_oxidation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):
                                                
                                                yield (leftIndexIncluded, rightIndexIncluded)

        elif (len(digestTypes) == 1
                  and "nonspecific" in digestTypes
                  and isSequenceGivenByUser == True):

            maxPeptideLengthInMatch = int(self.config["missedCleavages"])
            
            if len(listOfPhosphorylationSitePositions) == 0:

                if len(listOfCarbamidomethylationSitePositions) == 0:

                    if len(listOfDeamidationSitePositions) == 0:

                        if len(listOfOxidationSitePositions) == 0:
                            
                            for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                           sitePosition + 1):
                                for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                                   maxPeptideLengthInMatch + 1):
                                    rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                    # If the peptide left-right indices are NOT out of bound,
                                    #     yield the left-right index pair
                                    if (leftIndexIncluded > 0 and
                                            rightIndexIncluded <= aminoAcidSequenceLength):

                                        yield (leftIndexIncluded, rightIndexIncluded)

                        else:

                            for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                           sitePosition + 1):
                                for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                                   maxPeptideLengthInMatch + 1):
                                    rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                    # If the peptide left-right indices are NOT out of bound,
                                    #     yield the left-right index pair
                                    if (leftIndexIncluded > 0 and
                                            rightIndexIncluded <= aminoAcidSequenceLength):

                                        yield (leftIndexIncluded, rightIndexIncluded)

                    else:

                        if len(listOfOxidationSitePositions) == 0:
                            
                            for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                           sitePosition + 1):
                                for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                                   maxPeptideLengthInMatch + 1):
                                    rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                    # If the peptide left-right indices are NOT out of bound,
                                    #     yield the left-right index pair
                                    if (leftIndexIncluded > 0 and
                                            rightIndexIncluded <= aminoAcidSequenceLength):

                                        for situation in generator__binary_select_deamidation_sites(
                                                leftIndexIncluded, rightIndexIncluded):
                                            
                                            yield (leftIndexIncluded, rightIndexIncluded)

                        else:

                            for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                           sitePosition + 1):
                                for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                                   maxPeptideLengthInMatch + 1):
                                    rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                    # If the peptide left-right indices are NOT out of bound,
                                    #     yield the left-right index pair
                                    if (leftIndexIncluded > 0 and
                                            rightIndexIncluded <= aminoAcidSequenceLength):

                                        for situation in generator__binary_select_deamidation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in generator__binary_select_oxidation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):
                                            
                                                yield (leftIndexIncluded, rightIndexIncluded)


                else:

                    if len(listOfDeamidationSitePositions) == 0:

                        if len(listOfOxidationSitePositions) == 0:
                            
                            for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                           sitePosition + 1):
                                for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                                   maxPeptideLengthInMatch + 1):
                                    rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                    # If the peptide left-right indices are NOT out of bound,
                                    #     yield the left-right index pair
                                    if (leftIndexIncluded > 0 and
                                            rightIndexIncluded <= aminoAcidSequenceLength):

                                        for situation in self.generator__binary_select_carbamidomethylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):
                                        
                                            yield (leftIndexIncluded, rightIndexIncluded)

                        else:

                            for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                           sitePosition + 1):
                                for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                                   maxPeptideLengthInMatch + 1):
                                    rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                    # If the peptide left-right indices are NOT out of bound,
                                    #     yield the left-right index pair
                                    if (leftIndexIncluded > 0 and
                                            rightIndexIncluded <= aminoAcidSequenceLength):

                                        for situation in self.generator__binary_select_carbamidomethylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in self.generator__binary_select_oxidation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):
                                                
                                                yield (leftIndexIncluded, rightIndexIncluded)


                    else:

                        if len(listOfOxidationSitePositions) == 0:
                            
                            for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                           sitePosition + 1):
                                for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                                   maxPeptideLengthInMatch + 1):
                                    rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                    # If the peptide left-right indices are NOT out of bound,
                                    #     yield the left-right index pair
                                    if (leftIndexIncluded > 0 and
                                            rightIndexIncluded <= aminoAcidSequenceLength):

                                        for situation in self.generator__binary_select_carbamidomethylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in self.generator__binary_select_deamidation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):
                                                
                                                yield (leftIndexIncluded, rightIndexIncluded)

                        else:

                            for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                           sitePosition + 1):
                                for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                                   maxPeptideLengthInMatch + 1):
                                    rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                    # If the peptide left-right indices are NOT out of bound,
                                    #     yield the left-right index pair
                                    if (leftIndexIncluded > 0 and
                                            rightIndexIncluded <= aminoAcidSequenceLength):

                                        for situation in self.generator__binary_select_carbamidomethylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in self.generator__binary_select_deamidation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):

                                                for situation in self.generator__binary_select_oxidation_sites(
                                                        leftIndexIncluded, rightIndexIncluded):
                                                    
                                                    yield (leftIndexIncluded, rightIndexIncluded)

            else:

                if len(listOfCarbamidomethylationSitePositions) == 0:

                    if len(listOfDeamidationSitePositions) == 0:

                        if len(listOfOxidationSitePositions) == 0:
                            
                            for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                           sitePosition + 1):
                                for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                                   maxPeptideLengthInMatch + 1):
                                    rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                    # If the peptide left-right indices are NOT out of bound,
                                    #     yield the left-right index pair
                                    if (leftIndexIncluded > 0 and
                                            rightIndexIncluded <= aminoAcidSequenceLength):

                                        for situation in self.generator__binary_select_phosphorylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            yield (leftIndexIncluded, rightIndexIncluded)

                        else:

                            for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                           sitePosition + 1):
                                for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                                   maxPeptideLengthInMatch + 1):
                                    rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                    # If the peptide left-right indices are NOT out of bound,
                                    #     yield the left-right index pair
                                    if (leftIndexIncluded > 0 and
                                            rightIndexIncluded <= aminoAcidSequenceLength):

                                        for situation in self.generator__binary_select_phosphorylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in self.generator__binary_select_oxidation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):
                                                
                                                yield (leftIndexIncluded, rightIndexIncluded)

                    else:

                        if len(listOfOxidationSitePositions) == 0:
                            
                            for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                           sitePosition + 1):
                                for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                                   maxPeptideLengthInMatch + 1):
                                    rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                    # If the peptide left-right indices are NOT out of bound,
                                    #     yield the left-right index pair
                                    if (leftIndexIncluded > 0 and
                                            rightIndexIncluded <= aminoAcidSequenceLength):

                                        for situation in self.generator__binary_select_phosphorylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in self.generator__binary_select_deamidation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):
                                                
                                                yield (leftIndexIncluded, rightIndexIncluded)

                        else:

                            for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                           sitePosition + 1):
                                for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                                   maxPeptideLengthInMatch + 1):
                                    rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                    # If the peptide left-right indices are NOT out of bound,
                                    #     yield the left-right index pair
                                    if (leftIndexIncluded > 0 and
                                            rightIndexIncluded <= aminoAcidSequenceLength):

                                        for situation in self.generator__binary_select_phosphorylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in self.generator__binary_select_deamidation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):

                                                for situation in self.generator__binary_select_oxidation_sites(
                                                        leftIndexIncluded, rightIndexIncluded):
                                                    
                                                    yield (leftIndexIncluded, rightIndexIncluded)

                else:

                    if len(listOfDeamidationSitePositions) == 0:

                        if len(listOfOxidationSitePositions) == 0:
                            
                            for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                           sitePosition + 1):
                                for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                                   maxPeptideLengthInMatch + 1):
                                    rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                    # If the peptide left-right indices are NOT out of bound,
                                    #     yield the left-right index pair
                                    if (leftIndexIncluded > 0 and
                                            rightIndexIncluded <= aminoAcidSequenceLength):

                                        for situation in self.generator__binary_select_phosphorylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in self.generator__binary_select_carbamidomethylation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):
                                                
                                                yield (leftIndexIncluded, rightIndexIncluded)

                        else:

                            for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                           sitePosition + 1):
                                for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                                   maxPeptideLengthInMatch + 1):
                                    rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                    # If the peptide left-right indices are NOT out of bound,
                                    #     yield the left-right index pair
                                    if (leftIndexIncluded > 0 and
                                            rightIndexIncluded <= aminoAcidSequenceLength):

                                        for situation in self.generator__binary_select_phosphorylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in self.generator__binary_select_carbamidomethylation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):

                                                for situation in self.generator__binary_select_oxidation_sites(
                                                        leftIndexIncluded, rightIndexIncluded):
                                                    
                                                    yield (leftIndexIncluded, rightIndexIncluded)

                    else:

                        if len(listOfOxidationSitePositions) == 0:
                        
                            for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                           sitePosition + 1):
                                for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                                   maxPeptideLengthInMatch + 1):
                                    rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                    # If the peptide left-right indices are NOT out of bound,
                                    #     yield the left-right index pair
                                    if (leftIndexIncluded > 0 and
                                            rightIndexIncluded <= aminoAcidSequenceLength):

                                        for situation in self.generator__binary_select_phosphorylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in self.generator__binary_select_carbamidomethylation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):

                                                for situation in self.generator__binary_select_deamidation_sites(
                                                        leftIndexIncluded, rightIndexIncluded):
                                                    
                                                    yield (leftIndexIncluded, rightIndexIncluded)

                        else:
                            
                            for leftIndexIncluded in range(sitePosition - maxPeptideLengthInMatch + 1,
                                                           sitePosition + 1):
                                for peptideIntervalLength in range(sitePosition - leftIndexIncluded + 1,
                                                                   maxPeptideLengthInMatch + 1):
                                    rightIndexIncluded = leftIndexIncluded + peptideIntervalLength - 1
                                    # If the peptide left-right indices are NOT out of bound,
                                    #     yield the left-right index pair
                                    if (leftIndexIncluded > 0 and
                                            rightIndexIncluded <= aminoAcidSequenceLength):

                                        for situation in self.generator__binary_select_phosphorylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in self.generator__binary_select_carbamidomethylation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):

                                                for situation in self.generator__binary_select_deamidation_sites(
                                                        leftIndexIncluded, rightIndexIncluded):

                                                    for situation in self.generator__binary_select_oxidation_sites(
                                                            leftIndexIncluded, rightIndexIncluded):
                                                        
                                                        yield (leftIndexIncluded, rightIndexIncluded)
                   
        elif ("nonspecific" not in digestTypes
                  and isSequenceGivenByUser == False):

            missedCleavages = int(self.config["missedCleavages"])
            
            if len(listOfCarbamidomethylationSitePositions) == 0:

                if len(listOfDeamidationSitePositions) == 0:

                    if len(listOfOxidationSitePositions) == 0:
                        
                        for numberOfMisses in range(missedCleavages + 1):
                            print("Checking all glycopeptides that allow " +
                                          str(numberOfMisses) + " miss cleavages...")
                            for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                           numberOfMisses):
                                leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                # If we are doing trypsin cleavage
                                # and the glycosylation site is not within the fragment,
                                # We skip the checking, because there's no glycopeptide
                                if (leftIndexIncluded <= sitePosition and
                                        sitePosition <= rightIndexIncluded):
                                    
                                    yield (leftIndexIncluded, rightIndexIncluded)

                    else:
                        
                        for numberOfMisses in range(missedCleavages + 1):
                            print("Checking all glycopeptides that allow " +
                                          str(numberOfMisses) + " miss cleavages...")
                            for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                           numberOfMisses):
                                leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                # If we are doing trypsin cleavage
                                # and the glycosylation site is not within the fragment,
                                # We skip the checking, because there's no glycopeptide
                                if (leftIndexIncluded <= sitePosition and
                                        sitePosition <= rightIndexIncluded):

                                    for situation in generator__binary_select_oxidation_sites(
                                            leftIndexIncluded, rightIndexIncluded):
                                        
                                        yield (leftIndexIncluded, rightIndexIncluded)

                else:

                    if len(listOfOxidationSitePositions) == 0:

                        for numberOfMisses in range(missedCleavages + 1):
                            print("Checking all glycopeptides that allow " +
                                          str(numberOfMisses) + " miss cleavages...")
                            for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                           numberOfMisses):
                                leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                # If we are doing trypsin cleavage
                                # and the glycosylation site is not within the fragment,
                                # We skip the checking, because there's no glycopeptide
                                if (leftIndexIncluded <= sitePosition and
                                        sitePosition <= rightIndexIncluded):

                                    for situation in generator__binary_select_deamidation_sites(
                                            leftIndexIncluded, rightIndexIncluded):
                                        
                                        yield (leftIndexIncluded, rightIndexIncluded)

                    else:

                        for numberOfMisses in range(missedCleavages + 1):
                            print("Checking all glycopeptides that allow " +
                                          str(numberOfMisses) + " miss cleavages...")
                            for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                           numberOfMisses):
                                leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                # If we are doing trypsin cleavage
                                # and the glycosylation site is not within the fragment,
                                # We skip the checking, because there's no glycopeptide
                                if (leftIndexIncluded <= sitePosition and
                                        sitePosition <= rightIndexIncluded):

                                    for situation in generator__binary_select_deamidation_sites(
                                            leftIndexIncluded, rightIndexIncluded):

                                        for situation in generator__binary_select_oxidation_sites(
                                                leftIndexIncluded, rightIndexIncluded):
                                            
                                            yield (leftIndexIncluded, rightIndexIncluded)


            else:

                if len(listOfDeamidationSitePositions) == 0:

                    if len(listOfOxidationSitePositions) == 0:
                        
                        for numberOfMisses in range(missedCleavages + 1):
                            print("Checking all glycopeptides that allow " +
                                          str(numberOfMisses) + " miss cleavages...")
                            for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                           numberOfMisses):
                                leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                # If are doing trypsin cleavage
                                # and we the glycosylation site is not within the fragment,
                                # We skip the checking, because there's no glycopeptide
                                if (leftIndexIncluded <= sitePosition and
                                        sitePosition <= rightIndexIncluded):

                                    for situation in self.generator__binary_select_carbamidomethylation_sites(
                                            leftIndexIncluded, rightIndexIncluded):
                                        
                                        yield (leftIndexIncluded, rightIndexIncluded)

                    else:
                        
                        for numberOfMisses in range(missedCleavages + 1):
                            print("Checking all glycopeptides that allow " +
                                          str(numberOfMisses) + " miss cleavages...")
                            for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                           numberOfMisses):
                                leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                # If are doing trypsin cleavage
                                # and we the glycosylation site is not within the fragment,
                                # We skip the checking, because there's no glycopeptide
                                if (leftIndexIncluded <= sitePosition and
                                        sitePosition <= rightIndexIncluded):

                                    for situation in self.generator__binary_select_carbamidomethylation_sites(
                                            leftIndexIncluded, rightIndexIncluded):

                                        for situation in self.generator__binary_select_oxidation_sites(
                                                leftIndexIncluded, rightIndexIncluded):
                                            
                                            yield (leftIndexIncluded, rightIndexIncluded)

                else:

                    if len(listOfOxidationSitePositions) == 0:
                        
                        for numberOfMisses in range(missedCleavages + 1):
                            print("Checking all glycopeptides that allow " +
                                          str(numberOfMisses) + " miss cleavages...")
                            for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                           numberOfMisses):
                                leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                # If are doing trypsin cleavage
                                # and we the glycosylation site is not within the fragment,
                                # We skip the checking, because there's no glycopeptide
                                if (leftIndexIncluded <= sitePosition and
                                        sitePosition <= rightIndexIncluded):

                                    for situation in self.generator__binary_select_carbamidomethylation_sites(
                                            leftIndexIncluded, rightIndexIncluded):

                                        for situation in self.generator__binary_select_deamidation_sites(
                                                leftIndexIncluded, rightIndexIncluded):
                                            
                                            yield (leftIndexIncluded, rightIndexIncluded)

                    else:
                        
                        for numberOfMisses in range(missedCleavages + 1):
                            print("Checking all glycopeptides that allow " +
                                          str(numberOfMisses) + " miss cleavages...")
                            for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                           numberOfMisses):
                                leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                # If are doing trypsin cleavage
                                # and we the glycosylation site is not within the fragment,
                                # We skip the checking, because there's no glycopeptide
                                if (leftIndexIncluded <= sitePosition and
                                        sitePosition <= rightIndexIncluded):

                                    for situation in self.generator__binary_select_carbamidomethylation_sites(
                                            leftIndexIncluded, rightIndexIncluded):

                                        for situation in self.generator__binary_select_deamidation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in self.generator__binary_select_oxidation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):
                                                
                                                yield (leftIndexIncluded, rightIndexIncluded)

        elif ("nonspecific" not in digestTypes
                  and isSequenceGivenByUser == True):

            missedCleavages = int(self.config["missedCleavages"])
            
            if len(listOfPhosphorylationSitePositions) == 0:

                if len(listOfCarbamidomethylationSitePositions) == 0:

                    if len(listOfDeamidationSitePositions) == 0:

                        if len(listOfOxidationSitePositions) == 0:
                            
                            for numberOfMisses in range(missedCleavages + 1):
                                print("Checking all glycopeptides that allow " +
                                              str(numberOfMisses) + " miss cleavages...")
                                for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                               numberOfMisses):
                                    leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                    rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                    # If we are doing trypsin cleavage
                                    # and the glycosylation site is not within the fragment,
                                    # We skip the checking, because there's no glycopeptide
                                    if (leftIndexIncluded <= sitePosition and
                                            sitePosition <= rightIndexIncluded):
                                        
                                        yield (leftIndexIncluded, rightIndexIncluded)

                        else:

                            for numberOfMisses in range(missedCleavages + 1):
                                print("Checking all glycopeptides that allow " +
                                              str(numberOfMisses) + " miss cleavages...")
                                for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                               numberOfMisses):
                                    leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                    rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                    # If we are doing trypsin cleavage
                                    # and the glycosylation site is not within the fragment,
                                    # We skip the checking, because there's no glycopeptide
                                    if (leftIndexIncluded <= sitePosition and
                                            sitePosition <= rightIndexIncluded):

                                        for situation in generator__binary_select_oxidation_sites(
                                                leftIndexIncluded, rightIndexIncluded):
                                            
                                            yield (leftIndexIncluded, rightIndexIncluded)

                    else:

                        if len(listOfOxidationSitePositions) == 0:
                            
                            for numberOfMisses in range(missedCleavages + 1):
                                print("Checking all glycopeptides that allow " +
                                              str(numberOfMisses) + " miss cleavages...")
                                for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                               numberOfMisses):
                                    leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                    rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                    # If we are doing trypsin cleavage
                                    # and the glycosylation site is not within the fragment,
                                    # We skip the checking, because there's no glycopeptide
                                    if (leftIndexIncluded <= sitePosition and
                                            sitePosition <= rightIndexIncluded):

                                        for situation in generator__binary_select_deamidation_sites(
                                                leftIndexIncluded, rightIndexIncluded):
                                            
                                            yield (leftIndexIncluded, rightIndexIncluded)

                        else:

                            for numberOfMisses in range(missedCleavages + 1):
                                print("Checking all glycopeptides that allow " +
                                              str(numberOfMisses) + " miss cleavages...")
                                for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                               numberOfMisses):
                                    leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                    rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                    # If we are doing trypsin cleavage
                                    # and the glycosylation site is not within the fragment,
                                    # We skip the checking, because there's no glycopeptide
                                    if (leftIndexIncluded <= sitePosition and
                                            sitePosition <= rightIndexIncluded):

                                        for situation in generator__binary_select_deamidation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in generator__binary_select_oxidation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):

                                                yield (leftIndexIncluded, rightIndexIncluded)

                else:

                    if len(listOfDeamidationSitePositions) == 0:

                        if len(listOfOxidationSitePositions) == 0:

                            for numberOfMisses in range(missedCleavages + 1):
                                print("Checking all glycopeptides that allow " +
                                              str(numberOfMisses) + " miss cleavages...")
                                for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                               numberOfMisses):
                                    leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                    rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                    # If we are doing trypsin cleavage
                                    # and the glycosylation site is not within the fragment,
                                    # We skip the checking, because there's no glycopeptide
                                    if (leftIndexIncluded <= sitePosition and
                                            sitePosition <= rightIndexIncluded):
                                        
                                        for situation in self.generator__binary_select_carbamidomethylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):
                                        
                                            yield (leftIndexIncluded, rightIndexIncluded)

                        else:

                            for numberOfMisses in range(missedCleavages + 1):
                                print("Checking all glycopeptides that allow " +
                                              str(numberOfMisses) + " miss cleavages...")
                                for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                               numberOfMisses):
                                    leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                    rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                    # If we are doing trypsin cleavage
                                    # and the glycosylation site is not within the fragment,
                                    # We skip the checking, because there's no glycopeptide
                                    if (leftIndexIncluded <= sitePosition and
                                            sitePosition <= rightIndexIncluded):
                                        
                                        for situation in self.generator__binary_select_carbamidomethylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in self.generator__binary_select_oxidation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):
                                                
                                                yield (leftIndexIncluded, rightIndexIncluded)

                    else:

                        if len(listOfOxidationSitePositions) == 0:
                        
                            for numberOfMisses in range(missedCleavages + 1):
                                print("Checking all glycopeptides that allow " +
                                              str(numberOfMisses) + " miss cleavages...")
                                for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                               numberOfMisses):
                                    leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                    rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                    # If we are doing trypsin cleavage
                                    # and the glycosylation site is not within the fragment,
                                    # We skip the checking, because there's no glycopeptide
                                    if (leftIndexIncluded <= sitePosition and
                                            sitePosition <= rightIndexIncluded):
                                        
                                        for situation in self.generator__binary_select_carbamidomethylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in self.generator__binary_select_deamidation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):
                                                
                                                yield (leftIndexIncluded, rightIndexIncluded)

                        else:

                             for numberOfMisses in range(missedCleavages + 1):
                                print("Checking all glycopeptides that allow " +
                                              str(numberOfMisses) + " miss cleavages...")
                                for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                               numberOfMisses):
                                    leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                    rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                    # If we are doing trypsin cleavage
                                    # and the glycosylation site is not within the fragment,
                                    # We skip the checking, because there's no glycopeptide
                                    if (leftIndexIncluded <= sitePosition and
                                            sitePosition <= rightIndexIncluded):
                                        
                                        for situation in self.generator__binary_select_carbamidomethylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in self.generator__binary_select_deamidation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):

                                                for situation in self.generator__binary_select_oxidation_sites(
                                                        leftIndexIncluded, rightIndexIncluded):
                                                    
                                                    yield (leftIndexIncluded, rightIndexIncluded)

            else:

                if len(listOfCarbamidomethylationSitePositions) == 0:

                    if len(listOfDeamidationSitePositions) == 0:

                        if len(listOfOxidationSitePositions) == 0:
                            
                            for numberOfMisses in range(missedCleavages + 1):
                                print("Checking all glycopeptides that allow " +
                                              str(numberOfMisses) + " miss cleavages...")
                                for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                               numberOfMisses):
                                    leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                    rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                    # If we are doing trypsin cleavage
                                    # and the glycosylation site is not within the fragment,
                                    # We skip the checking, because there's no glycopeptide
                                    if (leftIndexIncluded <= sitePosition and
                                            sitePosition <= rightIndexIncluded):
                                        
                                        for situation in self.generator__binary_select_phosphorylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):
                                            
                                            yield (leftIndexIncluded, rightIndexIncluded)

                        else:

                            for numberOfMisses in range(missedCleavages + 1):
                                print("Checking all glycopeptides that allow " +
                                              str(numberOfMisses) + " miss cleavages...")
                                for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                               numberOfMisses):
                                    leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                    rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                    # If we are doing trypsin cleavage
                                    # and the glycosylation site is not within the fragment,
                                    # We skip the checking, because there's no glycopeptide
                                    if (leftIndexIncluded <= sitePosition and
                                            sitePosition <= rightIndexIncluded):
                                        
                                        for situation in self.generator__binary_select_phosphorylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in self.generator__binary_select_oxidation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):
                                                
                                                yield (leftIndexIncluded, rightIndexIncluded)

                    else:

                        if len(listOfOxidationSitePositions) == 0:
                            
                            for numberOfMisses in range(missedCleavages + 1):
                                print("Checking all glycopeptides that allow " +
                                              str(numberOfMisses) + " miss cleavages...")
                                for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                               numberOfMisses):
                                    leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                    rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                    # If we are doing trypsin cleavage
                                    # and the glycosylation site is not within the fragment,
                                    # We skip the checking, because there's no glycopeptide
                                    if (leftIndexIncluded <= sitePosition and
                                            sitePosition <= rightIndexIncluded):
                                        
                                        for situation in self.generator__binary_select_phosphorylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in self.generator__binary_select_deamidation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):
                                                
                                                yield (leftIndexIncluded, rightIndexIncluded)

                        else:

                            for numberOfMisses in range(missedCleavages + 1):
                                print("Checking all glycopeptides that allow " +
                                              str(numberOfMisses) + " miss cleavages...")
                                for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                               numberOfMisses):
                                    leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                    rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                    # If we are doing trypsin cleavage
                                    # and the glycosylation site is not within the fragment,
                                    # We skip the checking, because there's no glycopeptide
                                    if (leftIndexIncluded <= sitePosition and
                                            sitePosition <= rightIndexIncluded):
                                        
                                        for situation in self.generator__binary_select_phosphorylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in self.generator__binary_select_deamidation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):

                                                for situation in self.generator__binary_select_oxidation_sites(
                                                        leftIndexIncluded, rightIndexIncluded):
                                                    
                                                    yield (leftIndexIncluded, rightIndexIncluded)

                else:

                    if len(listOfDeamidationSitePositions) == 0:

                        if len(listOfOxidationSitePositions) == 0:
                            
                            for numberOfMisses in range(missedCleavages + 1):
                                print("Checking all glycopeptides that allow " +
                                              str(numberOfMisses) + " miss cleavages...")
                                for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                               numberOfMisses):
                                    leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                    rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                    # If we are doing trypsin cleavage
                                    # and the glycosylation site is not within the fragment,
                                    # We skip the checking, because there's no glycopeptide
                                    if (leftIndexIncluded <= sitePosition and
                                            sitePosition <= rightIndexIncluded):
                                        
                                        for situation in self.generator__binary_select_phosphorylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in self.generator__binary_select_carbamidomethylation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):
                                                
                                                yield (leftIndexIncluded, rightIndexIncluded)

                        else:

                            for numberOfMisses in range(missedCleavages + 1):
                                print("Checking all glycopeptides that allow " +
                                              str(numberOfMisses) + " miss cleavages...")
                                for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                               numberOfMisses):
                                    leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                    rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                    # If we are doing trypsin cleavage
                                    # and the glycosylation site is not within the fragment,
                                    # We skip the checking, because there's no glycopeptide
                                    if (leftIndexIncluded <= sitePosition and
                                            sitePosition <= rightIndexIncluded):
                                        
                                        for situation in self.generator__binary_select_phosphorylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in self.generator__binary_select_carbamidomethylation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):

                                                for situation in self.generator__binary_select_oxidation_sites(
                                                        leftIndexIncluded, rightIndexIncluded):
                                                    
                                                    yield (leftIndexIncluded, rightIndexIncluded)


                    else:

                        if len(listOfOxidationSitePositions) == 0:
                            
                            for numberOfMisses in range(missedCleavages + 1):
                                print("Checking all glycopeptides that allow " +
                                              str(numberOfMisses) + " miss cleavages...")
                                for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                               numberOfMisses):
                                    leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                    rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                    # If we are doing trypsin cleavage
                                    # and the glycosylation site is not within the fragment,
                                    # We skip the checking, because there's no glycopeptide
                                    if (leftIndexIncluded <= sitePosition and
                                            sitePosition <= rightIndexIncluded):
                                        
                                        for situation in self.generator__binary_select_phosphorylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in self.generator__binary_select_carbamidomethylation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):

                                                for situation in self.generator__binary_select_deamidation_sites(
                                                        leftIndexIncluded, rightIndexIncluded):
                                                    
                                                    yield (leftIndexIncluded, rightIndexIncluded)

                        else:

                            for numberOfMisses in range(missedCleavages + 1):
                                print("Checking all glycopeptides that allow " +
                                              str(numberOfMisses) + " miss cleavages...")
                                for i in range(len(listOfFragmentLeftRightIndexPairs) -
                                                               numberOfMisses):
                                    leftIndexIncluded = listOfFragmentLeftRightIndexPairs[i][0]
                                    rightIndexIncluded = listOfFragmentLeftRightIndexPairs[i + numberOfMisses][1]
                                    # If we are doing trypsin cleavage
                                    # and the glycosylation site is not within the fragment,
                                    # We skip the checking, because there's no glycopeptide
                                    if (leftIndexIncluded <= sitePosition and
                                            sitePosition <= rightIndexIncluded):
                                        
                                        for situation in self.generator__binary_select_phosphorylation_sites(
                                                leftIndexIncluded, rightIndexIncluded):

                                            for situation in self.generator__binary_select_carbamidomethylation_sites(
                                                    leftIndexIncluded, rightIndexIncluded):

                                                for situation in self.generator__binary_select_deamidation_sites(
                                                        leftIndexIncluded, rightIndexIncluded):

                                                    for situation in self.generator__binary_select_oxidation_sites(
                                                            leftIndexIncluded, rightIndexIncluded):
                                                        
                                                        yield (leftIndexIncluded, rightIndexIncluded)

        print("\n\n")



    def generator__binary_select_phosphorylation_sites(self,
                                                       leftIndexIncluded,
                                                       rightIndexIncluded,
                                                       ):

        # Function Description:
        #     Decide a subset of potential phosphorylation sites
        #         that lies within the intervals
        #         specified by leftIndexIncluded and rightIndexIncluded
        #     Run a set of binary selections
        #         on this subset of potential phosphorylation sites
        #         and modify the original amino acids to phosphorylated ones
        #             return to original loop to do mass match
        #             and then restore the original amino acids back
        # Function Input:
        #     The left side index inclusive for the interval
        #     The right side index inclusive for the interval
        # Function Output:
        #     Continuing modification of amino acids
        
        # ------ Pass member variable values to local variables ------
        aminoAcidSequenceDict = self.aminoAcidSequenceDict
        listOfPhosphorylationSitePositions = self.listOfPhosphorylationSitePositions
        listOfPhosphorylationSitesInInterval = []
        
        # Algorithm:
        #     binaryUpperBound gets doubled when a phosporylation site
        #         is found to be in the interval
        #     and UpperBound is then converted to a binary number
        #     then for all binary numbers from UpperBound to UpperBound * 2
        #     we get the string, left strip the "0b1",
        #     and then selectively assign all the phosphorylation sites within the interval
        #         to be phosphorylation sites in a binary behavior
        #     the subset of phosphorylation sites are selected
        #         to be replaced with "S^phos'ed", "T^phos'ed", "Y^phos'ed", etc.
        #     and converted back to "S", "T", "Y", etc. after one mass checking
        #     then we come to the next subset of phosphorylation sites
        #     do the same thing back and forth
        
        # Example:
        #     If there are 5 phosphorylation sites in the interval
        #     UpperBound is assigned to be 32, which is 0b100000 in binary
        #     then we go through all binary numbers from UpperBound to UpperBound * 2
        #     which is 0b100000, 0b100001, 0b100010, 0b100011, ..., 0b111111
        #     and truncate the "0b1",
        #     which makes 00000, 00001, 00010, 00011, ..., 11111
        #     and then, for example, 01001 means the second and fifth
        #         phosphorylation sites are binary selectively chosen to be active
        #     the subset of phosphorylation sites are selected
        #         to be replaced with "S^phos'ed", "T^phos'ed", "Y^phos'ed", etc.
        #     and converted back to "S", "T", "Y", etc. after one mass checking
        #     then we come to the next subset of phosphorylation sites
        #     do the same thing back and forth
        
        UpperBound = 1
        AMINO_ACID_NAME_LENGTH_ADDITION = self.AMINO_ACID_NAME_LENGTH_ADDITION
        
        # Decide the potential phosphorylation sites inside the interval
        for i in range(leftIndexIncluded, rightIndexIncluded + 1):
            if i in listOfPhosphorylationSitePositions:
                UpperBound = UpperBound * 2
                listOfPhosphorylationSitesInInterval.append(i)

        # For all phosphorylation types, run the binary selection
        if len(listOfPhosphorylationSitesInInterval) > 0:
            print("Performing binary phosphorylation on " +
                          str(len(listOfPhosphorylationSitesInInterval)) +
                          " phosphorylation sites from Index " +
                          str(leftIndexIncluded) + " to Index " +
                          str(rightIndexIncluded) + " on peptide chain...")
            
            # Loop through from 0b100001 to 0b111111,
            # which is essentially 00001 to 11111
            for integerChoice in range(UpperBound, UpperBound * 2):
                binaryChoice = bin(integerChoice)
                stringChoice = str(binaryChoice)[3:]

                # Modify all original amino acids to phosphorylated amino acids
                for i in range(len(stringChoice)):
                    if stringChoice[i] == '1':
                        phosphorylatedIndex = listOfPhosphorylationSitesInInterval[i]
                        aminoAcidSequenceDict[phosphorylatedIndex] += "^phos'ed"

                self.aminoAcidSequenceDict = aminoAcidSequenceDict

                # Go back to the original loop to find mass matching
                yield 0

                # Turn all phosphorylated amino acids back into original amino acids
                #     by deleting the string of length AMINO_ACID_NAME_LENGTH_ADDITION
                #         from the very last of the string itself
                # By default, it's deleting the very last 8 characters of the string itself
                for i in range(len(stringChoice)):
                    if stringChoice[i] == '1':
                        phosphorylatedIndex = listOfPhosphorylationSitesInInterval[i]
                        aminoAcidSequenceDict[phosphorylatedIndex] = (
                                aminoAcidSequenceDict[phosphorylatedIndex][:(
                                        - AMINO_ACID_NAME_LENGTH_ADDITION)])

                self.aminoAcidSequenceDict = aminoAcidSequenceDict

        # If there are no phosphorylation type selected in the user input xml file
        #     Do the original for loop
        else:
            yield 0

        # Print out a Done for finishing one set of binary selections
        #     based on one interval.
        # Nothing magic here.
        if len(listOfPhosphorylationSitesInInterval) > 0:
            print("Done.\n")



    def generator__binary_select_carbamidomethylation_sites(self,
                                                            leftIndexIncluded,
                                                            rightIndexIncluded,
                                                            ):

        # Function Description:
        #     Decide a subset of potential carbamidomethylation sites
        #         that lies within the intervals
        #         specified by leftIndexIncluded and rightIndexIncluded
        #     Run a set of binary selections
        #         on this subset of potential carbamidomethylation sites
        #         and modify the original amino acids to carbamidomethylated ones
        #             return to original loop to do mass match
        #             and then restore the original amino acids back
        # Function Input:
        #     The left side index inclusive for the interval
        #     The right side index inclusive for the interval
        # Function Output:
        #     Continuing modification of amino acids
        
        # ------ Pass member variable values to local variables ------
        aminoAcidSequenceDict = self.aminoAcidSequenceDict
        listOfCarbamidomethylationSitePositions = self.listOfCarbamidomethylationSitePositions
        listOfCarbamidomethylationSitesInInterval = []
        
        # Algorithm:
        #     binaryUpperBound gets doubled when a carbamidomethylation site
        #         is found to be in the interval
        #     and UpperBound is then converted to a binary number
        #     then for all binary numbers from UpperBound to UpperBound * 2
        #     we get the string, left strip the "0b1",
        #     and then selectively assign all the carbamidomethylation sites within the interval
        #         to be carbamidomethylation sites in a binary behavior
        #     the subset of carbamidomethylation sites are selected
        #         to be replaced with "C^carb'ed", etc.
        #     and converted back to "C", etc. after one mass checking
        #     then we come to the next subset of carbamidomethylation sites
        #     do the same thing back and forth
        
        # Example:
        #     If there are 5 carbamidomethylation sites in the interval
        #     UpperBound is assigned to be 32, which is 0b100000 in binary
        #     then we go through all binary numbers from UpperBound to UpperBound * 2
        #     which is 0b100000, 0b100001, 0b100010, 0b100011, ..., 0b111111
        #     and truncate the "0b1",
        #     which makes 00000, 00001, 00010, 00011, ..., 11111
        #     and then, for example, 01001 means the second and fifth
        #         carbamidomethylation sites are binary selectively chosen to be active
        #     the subset of carbamidomethylation sites are selected
        #         to be replaced with "C^carb'ed", etc.
        #     and converted back to "C", etc. after one mass checking
        #     then we come to the next subset of carbamidomethylation sites
        #     do the same thing back and forth
        
        UpperBound = 1
        AMINO_ACID_NAME_LENGTH_ADDITION = self.AMINO_ACID_NAME_LENGTH_ADDITION
        config = self.config
        
        # Decide the potential carbamidomethylation sites inside the interval
        for i in range(leftIndexIncluded, rightIndexIncluded + 1):
            if i in listOfCarbamidomethylationSitePositions:
                listOfCarbamidomethylationSitesInInterval.append(i)

        listOfBinarySelectedCarbamidomethylationSitesInInterval = [
                site
                for site in listOfCarbamidomethylationSitesInInterval
                if aminoAcidSequenceDict[site] not in config["listOfCarbamidomethylationsThatApplyToAll"]]
                        
        listOfApplyAllCarbamidomethylationSitesInInterval = [
                site
                for site in listOfCarbamidomethylationSitesInInterval
                if aminoAcidSequenceDict[site] in config["listOfCarbamidomethylationsThatApplyToAll"]]

        # Apply carbamidomethylations for amino acids that marked as "applyToAll"
        for applyAllSite in listOfApplyAllCarbamidomethylationSitesInInterval:
            aminoAcidSequenceDict[applyAllSite] += "^carb'ed"

        
        for i in range(leftIndexIncluded, rightIndexIncluded + 1):
            if i in listOfBinarySelectedCarbamidomethylationSitesInInterval:
                UpperBound = UpperBound * 2

                
        # For all BINARY carbamidomethylation types, run the binary selection
        if len(listOfBinarySelectedCarbamidomethylationSitesInInterval) > 0:

            print("Performing binary carbamidomethylation on " +
                          str(len(listOfBinarySelectedCarbamidomethylationSitesInInterval)) +
                          " carbamidomethylation sites from Index " +
                          str(leftIndexIncluded) + " to Index " +
                          str(rightIndexIncluded) + " on peptide chain...")
            
            # Loop through from 0b100001 to 0b111111,
            # which is essentially 00001 to 11111
            for integerChoice in range(UpperBound, UpperBound * 2):
                binaryChoice = bin(integerChoice)
                stringChoice = str(binaryChoice)[3:]

                # Modify all original amino acids to carbamidomethylated amino acids
                for i in range(len(stringChoice)):
                    if stringChoice[i] == '1':
                        carbamidomethylatedIndex = listOfBinarySelectedCarbamidomethylationSitesInInterval[i]
                        aminoAcidSequenceDict[carbamidomethylatedIndex] += "^carb'ed"

                self.aminoAcidSequenceDict = aminoAcidSequenceDict

                # Go back to the original loop to find mass matching
                yield 0

                # Turn all BINARY carbamidomethylated amino acids back into original amino acids
                #     by deleting the string of length AMINO_ACID_NAME_LENGTH_ADDITION
                #         from the very last of the string itself
                # By default, it's deleting the very last 8 characters of the string itself
                for i in range(len(stringChoice)):
                    if stringChoice[i] == '1':
                        carbamidomethylatedIndex = listOfBinarySelectedCarbamidomethylationSitesInInterval[i]
                        aminoAcidSequenceDict[carbamidomethylatedIndex] = (
                                aminoAcidSequenceDict[carbamidomethylatedIndex][:(
                                        - AMINO_ACID_NAME_LENGTH_ADDITION)])

                self.aminoAcidSequenceDict = aminoAcidSequenceDict

        # If there are no carbamidomethylation type selected in the user input xml file
        #     Do the original for loop
        else:
            yield 0

        # Cancel carbamidomethylations for amino acids that marked as "applyToAll"
        for applyAllSite in listOfApplyAllCarbamidomethylationSitesInInterval:
            aminoAcidSequenceDict[applyAllSite] = (
                    aminoAcidSequenceDict[applyAllSite][:(
                            - AMINO_ACID_NAME_LENGTH_ADDITION)])
                                                    
        # Print out a Done for finishing one set of binary selections
        #     based on one interval.
        # Nothing magic here.
        if len(listOfBinarySelectedCarbamidomethylationSitesInInterval) > 0:
            print("Done.\n")



    def generator__binary_select_deamidation_sites(self,
                                                   leftIndexIncluded,
                                                   rightIndexIncluded,
                                                   ):

        # Function Description:
        #     Decide a subset of potential deamidation sites
        #         that lies within the intervals
        #         specified by leftIndexIncluded and rightIndexIncluded
        #     Run a set of binary selections
        #         on this subset of potential deamidation sites
        #         and modify the original amino acids to deamidated ones
        #             return to original loop to do mass match
        #             and then restore the original amino acids back
        # Function Input:
        #     The left side index inclusive for the interval
        #     The right side index inclusive for the interval
        # Function Output:
        #     Continuing modification of amino acids
        
        # ------ Pass member variable values to local variables ------
        aminoAcidSequenceDict = self.aminoAcidSequenceDict
        listOfDeamidationSitePositions = self.listOfDeamidationSitePositions
        listOfDeamidationSitesInInterval = []
        
        # Algorithm:
        #     binaryUpperBound gets doubled when a deamidation site
        #         is found to be in the interval
        #     and UpperBound is then converted to a binary number
        #     then for all binary numbers from UpperBound to UpperBound * 2
        #     we get the string, left strip the "0b1",
        #     and then selectively assign all the deamidation sites within the interval
        #         to be deamidation sites in a binary behavior
        #     the subset of deamidation sites are selected
        #         to be replaced with "Q^deam'ed", "Q^phos'ed^deam'ed", etc.
        #     and converted back to "Q", "Q^phos'ed", etc. after one mass checking
        #     then we come to the next subset of deamidation sites
        #     do the same thing back and forth
        
        # Example:
        #     If there are 5 deamidation sites in the interval
        #     UpperBound is assigned to be 32, which is 0b100000 in binary
        #     then we go through all binary numbers from UpperBound to UpperBound * 2
        #     which is 0b100000, 0b100001, 0b100010, 0b100011, ..., 0b111111
        #     and truncate the "0b1",
        #     which makes 00000, 00001, 00010, 00011, ..., 11111
        #     and then, for example, 01001 means the second and fifth
        #         deamidation sites are binary selectively chosen to be active
        #     the subset of deamidation sites are selected
        #         to be replaced with "Q^deam'ed", "Q^phos'ed^deam'ed", etc.
        #     and converted back to "Q", "Q^phos'ed", etc. after one mass checking
        #     then we come to the next subset of deamidation sites
        #     do the same thing back and forth
        
        UpperBound = 1
        AMINO_ACID_NAME_LENGTH_ADDITION = self.AMINO_ACID_NAME_LENGTH_ADDITION
        
        # Decide the potential deamidation sites inside the interval
        for i in range(leftIndexIncluded, rightIndexIncluded + 1):
            if i in listOfDeamidationSitePositions:
                UpperBound = UpperBound * 2
                listOfDeamidationSitesInInterval.append(i)

        # For all deamidation types, run the binary selection
        if len(listOfDeamidationSitesInInterval) > 0:
            print("Performing binary deamidation on " +
                          str(len(listOfDeamidationSitesInInterval)) +
                          " deamidation sites from Index " +
                          str(leftIndexIncluded) + " to Index " +
                          str(rightIndexIncluded) + " on peptide chain...")
            
            # Loop through from 0b100001 to 0b111111,
            # which is essentially 00001 to 11111
            for integerChoice in range(UpperBound, UpperBound * 2):
                binaryChoice = bin(integerChoice)
                stringChoice = str(binaryChoice)[3:]

                # Modify all original amino acids to deamidated amino acids
                for i in range(len(stringChoice)):
                    if stringChoice[i] == '1':
                        deamidatedIndex = listOfDeamidationSitesInInterval[i]
                        aminoAcidSequenceDict[deamidatedIndex] += "^deam'ed"

                self.aminoAcidSequenceDict = aminoAcidSequenceDict

                # Go back to the original loop to find mass matching
                yield 0

                # Turn all deamidated amino acids back into original amino acids
                #     by deleting the string of length AMINO_ACID_NAME_LENGTH_ADDITION
                #         from the very last of the string itself
                # By default, it's deleting the very last 8 characters of the string itself
                for i in range(len(stringChoice)):
                    if stringChoice[i] == '1':
                        deamidatedIndex = listOfDeamidationSitesInInterval[i]
                        aminoAcidSequenceDict[deamidatedIndex] = (
                                aminoAcidSequenceDict[deamidatedIndex][:(
                                        - AMINO_ACID_NAME_LENGTH_ADDITION)])

                self.aminoAcidSequenceDict = aminoAcidSequenceDict

        # If there are no deamidation type selected in the user input xml file
        #     Do the original for loop
        else:
            yield 0

        # Print out a Done for finishing one set of binary selections
        #     based on one interval.
        # Nothing magic here.
        if len(listOfDeamidationSitesInInterval) > 0:
            print("Done.\n")

            

    def generator__binary_select_oxidation_sites(self,
                                                 leftIndexIncluded,
                                                 rightIndexIncluded,
                                                 ):

        # Function Description:
        #     Decide a subset of potential oxidation sites
        #         that lies within the intervals
        #         specified by leftIndexIncluded and rightIndexIncluded
        #     Run a set of binary selections
        #         on this subset of potential oxidation sites
        #         and modify the original amino acids to oxidated ones
        #             return to original loop to do mass match
        #             and then restore the original amino acids back
        # Function Input:
        #     The left side index inclusive for the interval
        #     The right side index inclusive for the interval
        # Function Output:
        #     Continuing modification of amino acids
        
        # ------ Pass member variable values to local variables ------
        aminoAcidSequenceDict = self.aminoAcidSequenceDict
        listOfOxidationSitePositions = self.listOfOxidationSitePositions
        listOfOxidationSitesInInterval = []
        
        # Algorithm:
        #     binaryUpperBound gets doubled when a oxidation site
        #         is found to be in the interval
        #     and UpperBound is then converted to a binary number
        #     then for all binary numbers from UpperBound to UpperBound * 2
        #     we get the string, left strip the "0b1",
        #     and then selectively assign all the oxidation sites within the interval
        #         to be oxidation sites in a binary behavior
        #     the subset of oxidation sites are selected
        #         to be replaced with "M^oxid'ed", "M^phos'ed^oxid'ed", etc.
        #     and converted back to "M", "M^phos'ed", etc. after one mass checking
        #     then we come to the next subset of oxidation sites
        #     do the same thing back and forth
        
        # Example:
        #     If there are 5 oxidation sites in the interval
        #     UpperBound is assigned to be 32, which is 0b100000 in binary
        #     then we go through all binary numbers from UpperBound to UpperBound * 2
        #     which is 0b100000, 0b100001, 0b100010, 0b100011, ..., 0b111111
        #     and truncate the "0b1",
        #     which makes 00000, 00001, 00010, 00011, ..., 11111
        #     and then, for example, 01001 means the second and fifth
        #         oxidation sites are binary selectively chosen to be active
        #     the subset of oxidation sites are selected
        #         to be replaced with "M^oxid'ed", "M^phos'ed^oxid'ed", etc.
        #     and converted back to "M", "M^phos'ed", etc. after one mass checking
        #     then we come to the next subset of oxidation sites
        #     do the same thing back and forth
        
        UpperBound = 1
        AMINO_ACID_NAME_LENGTH_ADDITION = self.AMINO_ACID_NAME_LENGTH_ADDITION
        
        # Decide the potential oxidation sites inside the interval
        for i in range(leftIndexIncluded, rightIndexIncluded + 1):
            if i in listOfOxidationSitePositions:
                UpperBound = UpperBound * 2
                listOfOxidationSitesInInterval.append(i)

        # For all deamidation types, run the binary selection
        if len(listOfOxidationSitesInInterval) > 0:
            print("Performing binary oxidation on " +
                          str(len(listOfOxidationSitesInInterval)) +
                          " oxidation sites from Index " +
                          str(leftIndexIncluded) + " to Index " +
                          str(rightIndexIncluded) + " on peptide chain...")
            
            # Loop through from 0b100001 to 0b111111,
            # which is essentially 00001 to 11111
            for integerChoice in range(UpperBound, UpperBound * 2):
                binaryChoice = bin(integerChoice)
                stringChoice = str(binaryChoice)[3:]

                # Modify all original amino acids to oxidated amino acids
                for i in range(len(stringChoice)):
                    if stringChoice[i] == '1':
                        oxidatedIndex = listOfOxidationSitesInInterval[i]
                        aminoAcidSequenceDict[oxidatedIndex] += "^oxid'ed"

                self.aminoAcidSequenceDict = aminoAcidSequenceDict

                # Go back to the original loop to find mass matching
                yield 0

                # Turn all oxidated amino acids back into original amino acids
                #     by deleting the string of length AMINO_ACID_NAME_LENGTH_ADDITION
                #         from the very last of the string itself
                # By default, it's deleting the very last 8 characters of the string itself
                for i in range(len(stringChoice)):
                    if stringChoice[i] == '1':
                        oxidatedIndex = listOfOxidationSitesInInterval[i]
                        aminoAcidSequenceDict[oxidatedIndex] = (
                                aminoAcidSequenceDict[oxidatedIndex][:(
                                        - AMINO_ACID_NAME_LENGTH_ADDITION)])

                self.aminoAcidSequenceDict = aminoAcidSequenceDict

        # If there are no oxidation type selected in the user input xml file
        #     Do the original for loop
        else:
            yield 0

        # Print out a Done for finishing one set of binary selections
        #     based on one interval.
        # Nothing magic here.
        if len(listOfOxidationSitesInInterval) > 0:
            print("Done.\n")

            
                   
    def get_glycan_mass(self, glycanDict):

        # Function Description:
        #     Calculate the glycan mass based on all component masses
        # Function Input:
        #     A glycan dictionary, specifying the number of each component
        # Function Output:
        #     The glycan mass

        # ------ Pass member variable values to local variables ------
        glycanComponentMassDict = self.glycanComponentMassDict
        waterMass = self.waterMass
        
        glycanMass = 0

        # Add up all glycan components mass to make the mass of a specific glycan
        for glycanComponent in glycanComponentMassDict.keys():
            glycanMass = glycanMass + (glycanComponentMassDict[glycanComponent] *
                                               glycanDict[glycanComponent])

        # Add one water mass to restore the water from dehydration with the peptide
        glycanMass = glycanMass + waterMass

        glycanMass = round(glycanMass, 10)

        return glycanMass



    def get_peptide_and_glycopeptide_molecular_mass(self,
                                                    leftIndexIncluded,
                                                    rightIndexIncluded,
                                                    glycanMassAttached):

        # Function Description:
        #     Calculation the glycopeptide mass
        #     Based on the dehydrated (with one water loss) amino acid masses
        #     and the unhydrated (without any water loss) glycan mass
        # Function Input:
        #     The left end index and the right end index
        #     The glycan attached
        # Function Output:
        #     The total mass of the glycopeptide
        
        # ------ Pass member variable values to local variables ------
        aminoAcidSequenceDict = self.aminoAcidSequenceDict
        aminoAcidDict = self.aminoAcidDict
        
        peptideMass = 0

        # Go through all amino acids one by one to get the mass of the peptide
        for aminoAcidIndex in range(leftIndexIncluded, rightIndexIncluded + 1):
            aminoAcidAtThisIndex = aminoAcidSequenceDict[aminoAcidIndex]
            aminoAcidMassAtThisIndex = aminoAcidDict[str(aminoAcidAtThisIndex)]["Mass"]
            peptideMass = peptideMass + aminoAcidMassAtThisIndex
        
        
        
        # Add the mass of the glycan attached to make the glycopeptide mass
        glycopeptideMass = peptideMass + glycanMassAttached

        peptideMass = round(peptideMass, 10)
        glycanMassAttached = round(glycanMassAttached, 10)
        glycopeptideMass = round(glycopeptideMass, 10)
        
        return (peptideMass, glycanMassAttached, glycopeptideMass)



class massMatch:

    def __init__(self,
                 matchIndex,
                 proteinFileName,
                 proteinName,
                 peptideStringMatched,
                 peptideStringLengthMatched,
                 glycanCodeMatched,
                 glycanCompositionMatched,
                 glycosylationSiteIndex,
                 glycosylationSiteAminoAcid,
                 glycosylationSiteAminoAcidMass,
                 peptideLeftIndexIncluded,
                 peptideLeftSideAminoAcid,
                 peptideLeftSideAminoAcidMass,
                 peptideRightIndexIncluded,
                 peptideRightSideAminoAcid,
                 peptideRightSideAminoAcidMass,
                 peptideMass,
                 glycanMass,
                 reducedGlycanMass,
                 glycopeptideMass,
                 totalTheoreticalMass,
                 isDecoyUsed,
                 decoyMassToMakeMatch,
                 inputMassToMatch,
                 deviationFromActualToTheoreticalMass,
                 deviationFromActualToDecoyedTheoreticalMass,
                 compoundMzValue,
                 compoundIntensity,
                 compoundRtCenter,
                 hydrogenMass):

        # Function Description:
        #     Set all variables from outside of the class to class member variables
        self.totalReScore = 0
        
        self.matchIndex = matchIndex

        self.proteinName = proteinName
        self.proteinFileName = proteinFileName
        
        self.peptideStringMatched = peptideStringMatched
        self.peptideStringLengthMatched = peptideStringLengthMatched

        self.glycanCodeMatched = glycanCodeMatched
        self.glycanCompositionMatched = glycanCompositionMatched
        
        self.glycosylationSiteIndex = glycosylationSiteIndex
        self.glycosylationSiteAminoAcid = glycosylationSiteAminoAcid
        self.glycosylationSiteAminoAcidMass = glycosylationSiteAminoAcidMass
        
        self.peptideLeftIndexIncluded = peptideLeftIndexIncluded
        self.peptideLeftSideAminoAcid = peptideLeftSideAminoAcid
        self.peptideLeftSideAminoAcidMass = peptideLeftSideAminoAcidMass

        self.peptideRightIndexIncluded = peptideRightIndexIncluded
        self.peptideRightSideAminoAcid = peptideRightSideAminoAcid
        self.peptideRightSideAminoAcidMass = peptideRightSideAminoAcidMass

        self.peptideMass = peptideMass
        self.glycanMass = glycanMass
        self.reducedGlycanMass = reducedGlycanMass
        self.glycopeptideMass = glycopeptideMass
        self.totalTheoreticalMass = totalTheoreticalMass

        self.isDecoyUsed = isDecoyUsed
        self.decoyMassToMakeMatch = decoyMassToMakeMatch
        
        self.inputMassToMatch = inputMassToMatch
        self.deviationFromActualToTheoreticalMass = deviationFromActualToTheoreticalMass
        self.deviationFromActualToDecoyedTheoreticalMass = deviationFromActualToDecoyedTheoreticalMass

        self.compoundMzValue = compoundMzValue
        self.compoundIntensity = compoundIntensity
        self.compoundRtCenter = compoundRtCenter
        
        self.hydrogenMass = hydrogenMass
        
        self.listOfFragmentNeutralMassMatches = []

        self.totalScore = 0

        

    def printMatchDetails(self):

        # Function Description:
        #     Print the details of the glycopeptide parent mass match
        
        print("******************************\n")
        
        print("Match " + str(self.matchIndex) + " (index starts at 1):")

        print("The Source Protein File Name: \t\t\t" + self.proteinFileName)

        print("The Protein Name: \t\t\t\t" + self.proteinName)
        
        print("The Peptide String Matched: \t\t\t" +
              str(self.peptideStringMatched))
        print("The Peptide Location: \t\t\t\t" +
              str(self.peptideLeftIndexIncluded) + "-" +
              str(self.peptideRightIndexIncluded))
        #print("The Length of the Peptide String Matched: \t" +
        #      str(self.peptideStringLengthMatched))
        #print("The Glycan Code Matched: \t\t\t" +
        #      self.glycanCodeMatched)
        print("The Composition of the Matched Glycan: \n")

        glycanComponentDict = {}
        
        for glycanComponent in self.glycanCompositionMatched.keys():
            if (glycanComponent != "Glycan Code" and
                    glycanComponent != "Carbohydrate Type" and
                    glycanComponent != "Class" and
                    glycanComponent != "Chemical Formula" and
                    glycanComponent != "Molecular Mass"):
                glycanComponentDict[glycanComponent] = self.glycanCompositionMatched[glycanComponent]
        
        for glycanComponent in glycanComponentDict.keys():
            print("\t" + glycanComponent + ": "
                  + str(glycanComponentDict[glycanComponent]))
        
        '''print("The Glycosylation Site: ")

        print("\tIndex: \t\t\t" +
              str(self.glycosylationSiteIndex))
        print("\tAmino Acid: \t\t" +
              str(self.glycosylationSiteAminoAcid))
        print("\tAmino Acid Mass: \t" +
              str(self.glycosylationSiteAminoAcidMass) + "\n")

        print("The Left End of the matched peptide: ")
        print("\tIndex: \t\t\t" +
              str(self.peptideLeftIndexIncluded))
        print("\tAmino Acid: \t\t" +
              str(self.peptideLeftSideAminoAcid))
        print("\tAmino Acid Mass: \t" +
              str(self.peptideLeftSideAminoAcidMass) + "\n")
              
        print("The Right End of the matched peptide: ")
        print("\tIndex: \t\t\t" +
              str(self.peptideRightIndexIncluded))
        print("\tAmino Acid: \t\t" +
              str(self.peptideRightSideAminoAcid))
        print("\tAmino Acid Mass: \t" +
              str(self.peptideRightSideAminoAcidMass) + "\n")'''
        
        print("The Peptide Mass minus one water mass: \t\t\t" +
              str(self.peptideMass))

        print("The Glycan Mass Matched: \t\t\t\t" +
              str(self.glycanMass))

        print("The Reduced Glycan Mass Matched: \t\t\t" +
              str(self.reducedGlycanMass))
            
        print("The Glycopeptide Mass Matched: \t\t\t\t" +
              str(self.glycopeptideMass))

        if self.isDecoyUsed == True:
            print("Decoy Used!!")
        else:
            print("Decoy unused.")
            
        print("The Decoy Mass to Make The Mass Match: \t\t\t" +
              str(self.decoyMassToMakeMatch))
        
        print("The Input Mass to Match: \t\t\t\t" +
              str(self.inputMassToMatch))

        print("The Deviation from Actual to Theoretical Mass: \t\t" +
              self.deviationFromActualToTheoreticalMass)

        print("The Deviation from Actual to Decoyed Theoretical Mass: \t" +
              str(self.deviationFromActualToDecoyedTheoreticalMass) + "\n")

        print("The Total Score for Fragment Matching: \t" +
              str(self.totalScore) + "\n")
        
        if self.listOfFragmentNeutralMassMatches != []:
            print("Fragment Neutral Mass Matchings Found!!")
            for fragmentNeutralMassMatch in self.listOfFragmentNeutralMassMatches:
                fragmentNeutralMassMatch.printFragmentMatchDetails()

        print("")



    def printMatchDetailsWithFeatureIndices(self):

        print("******************************\n")

        print("Feature " + str(self.featureIndex) + " (index starts at 1):")
        print("Match " + str(self.matchIndex) + " (index starts at 1):")
        
        print("The Source Protein File Name: \t\t\t" + self.proteinFileName)

        print("The Protein Name: \t\t\t\t" + self.proteinName)
        
        print("The Peptide String Matched: \t\t\t" +
              str(self.peptideStringMatched))
        print("The Peptide Location: \t\t\t\t" +
              str(self.peptideLeftIndexIncluded) + "-" +
              str(self.peptideRightIndexIncluded))
        #print("The Length of the Peptide String Matched: \t" +
        #      str(self.peptideStringLengthMatched))
        #print("The Glycan Code Matched: \t\t\t" +
        #      self.glycanCodeMatched)
        print("The Composition of the Matched Glycan:")

        glycanComponentDict = {}
        
        for glycanComponent in self.glycanCompositionMatched.keys():
            if (glycanComponent != "Glycan Code" and
                    glycanComponent != "Carbohydrate Type" and
                    glycanComponent != "Class" and
                    glycanComponent != "Chemical Formula" and
                    glycanComponent != "Molecular Mass"):
                glycanComponentDict[glycanComponent] = self.glycanCompositionMatched[glycanComponent]
        
        for glycanComponent in glycanComponentDict.keys():
            print("\t" + glycanComponent + ": "
                  + str(glycanComponentDict[glycanComponent]))
        
        '''print("The Glycosylation Site: ")

        print("\tIndex: \t\t\t" +
              str(self.glycosylationSiteIndex))
        print("\tAmino Acid: \t\t" +
              str(self.glycosylationSiteAminoAcid))
        print("\tAmino Acid Mass: \t" +
              str(self.glycosylationSiteAminoAcidMass) + "\n")

        print("The Left End of the matched peptide: ")
        print("\tIndex: \t\t\t" +
              str(self.peptideLeftIndexIncluded))
        print("\tAmino Acid: \t\t" +
              str(self.peptideLeftSideAminoAcid))
        print("\tAmino Acid Mass: \t" +
              str(self.peptideLeftSideAminoAcidMass) + "\n")
              
        print("The Right End of the matched peptide: ")
        print("\tIndex: \t\t\t" +
              str(self.peptideRightIndexIncluded))
        print("\tAmino Acid: \t\t" +
              str(self.peptideRightSideAminoAcid))
        print("\tAmino Acid Mass: \t" +
              str(self.peptideRightSideAminoAcidMass) + "\n")'''
        
        print("The Peptide Mass minus one water mass: \t\t\t" +
              str(self.peptideMass))

        print("The Glycan Mass Matched: \t\t\t\t" +
              str(self.glycanMass))

        if self.hydrogenMass != 0:
            print("The Reduced Glycan Mass Matched: \t\t\t" +
                  str(self.reducedGlycanMass))
            
        print("The Glycopeptide Mass Matched: \t\t\t\t" +
              str(self.glycopeptideMass))

        print("The Decoy Mass to Make The Mass Match: \t\t\t" +
              str(self.decoyMassToMakeMatch))
        
        print("The Input Mass to Match: \t\t\t\t" +
              str(self.inputMassToMatch))

        print("The Deviation from Actual to Theoretical Mass: \t\t" +
              self.deviationFromActualToTheoreticalMass)

        print("The Deviation from Actual to Decoyed Theoretical Mass: \t" +
              self.deviationFromActualToDecoyedTheoreticalMass + "\n")

        print("The Total Score for Fragment Matching: \t" +
              str(self.totalScore) + "\n")
        
        if self.listOfFragmentNeutralMassMatches != []:
            print("Fragment Neutral Mass Matchings Found!!")
            for fragmentNeutralMassMatch in self.listOfFragmentNeutralMassMatches:
                fragmentNeutralMassMatch.printFragmentMatchDetails()

        print("")
