class MassScoreObject():

    def __init__(self,
                 config,
                 massMatch):

        # Function Description:
        #     Set the user-defined configuration and mass match object
        #         to member variables
        #     Set the pool of fragments
        #         which defines all the fragments
        #         that is taken into account for scoring algorithm
        #     The pool of fragment differs from different configuration
        #         See details below
        
        self.config = config

        self.massMatch = massMatch
        self.matchedPeptideMass = self.massMatch.peptideMass
        self.matchedGlycanMass = self.massMatch.glycanMass
        
        self.totalScore = 0

        waterMass = 20.1565006428
        hydrogenMass = 1.0078250321
        
        # Hex:
        # C6H12O6
        hexMass = 180.0633881184
        # HexNAc:
        # C8H15NO6
        hexNAcMass = 221.0899372201
        
        # Hex-HexNAc:
        # C6H12O6 + C8H15NO6 - H2O = C14H25NO11
        HexHexNAcMinusOneWaterMass = 383.1427606521
        # HexNAc-HexNAc:
        # 2C8H15NO6 - H2O = C16H28N2O11
        HexNAcHexNAcMinusOneWaterMass = 424.1693097538
        # HexNAc-HexNAc-Hex:
        # 2C8H15NO6 + C6H12O6 - 2H2O = C22H38N2O16
        HexNAcHexNAcHexMinusTwoWaterMass = 586.2221331857
        # HexNAc-HexNAc-Hex-Hex:
        # 2C8H15NO6 + 2C6H12O6 - 3H2O = C28H48N2O21
        HexNAcHexNAcHexHexMinusThreeWaterMass = 748.2749566177
        
        # Hex-Neu5Ac-minus-H2O:
        # C6H12O6 + C11H19NO9 - 2H2O = C17H27NO13
        HexNeu5AcMinusTwoWaterMass = 453.1482399606
        # Hex-Neu5Ac:
        # C6H12O6 + C11H19NO9 - H2O = C17H29NO14
        HexNeu5AcMinusOneWaterMass = 471.158804647
        # Hex-HexNAc-Neu5Ac-minus-H2O:
        # C6H12O6 + C8H15NO6 + C11H19NO9 - 3H2O = C25H40N2O18
        HexHexNAcNeu5AcMinusThreeWaterMass = 656.2276124942
        # Hex-HexNAc-Neu5Ac:
        # C6H12O6 + C8H15NO6 + C11H19NO9 - 2H2O = C25H42N2O19
        HexHexNAcNeu5AcMinusTwoWaterMass = 674.2381771806
        
        # Hex-Neu5Gc-minus-H2O:
        # C6H12O6 + C11H19NO10 - 2H2O = C17H27NO14
        HexNeu5GcMinusTwoWaterMass = 469.1431545827
        # Hex-Neu5Gc-H2O:
        # C6H12O6 + C11H19NO10 - H2O = C17H29NO15
        HexNeu5GcMinusOneWaterMass = 487.1537192691
        # Hex-HexNAc-Neu5Gc-minus-H2O:
        # C6H12O6 + C8H15NO6 + C11H19NO10 - 3H2O = C25H40N2O19
        HexHexNAcNeu5GcMinusThreeWaterMass = 672.2225271164
        # Hex-HexNAc-Neu5Gc:
        # C6H12O6 + C8H15NO6 + C11H19NO10 - 2H2O = C25H42N2O20
        HexHexNAcNeu5GcMinusTwoWaterMass = 690.2330918028
        
        # Glycan core: HexNAc-HexNAc-Hex(Hex)2:
        # 2C8H15NO6 + 3C6H12O6 - 4H2O = C34H58N2O26
        glycanCoreMass = 910.3277800497

        self.fragmentPoolDict = {
            #-------------------- SCORE 4: --------------------#

            # Score 4 is for the attached glycan as a B-ion
            #     which happens due to the short time of enzyme digestion
            #     for a long time digestion,
            #     this will be really rare.
            
            # attached glycan:
            "glycan-minus-H2O": {"mass": self.matchedGlycanMass - waterMass,
                                 "score": 4},
            
            #-------------------- SCORE 1: --------------------#

            # Score 1 is for the simple glycan component
            #     including Hex, HexNAc.
            #     These species happen for all glycans
            #     So we give a low score for this
            #     In the future, high mannose glycan scoring
            #         may make these score-1 species higher scores
            #     But due to lack of experience in fragmenting high-mannose glycan
            #         We temporarily leave it here the way it is
            
            # All species

            # Hex minus one water:
            # C6H12O6 - H2O = C6H10O5
            "Hex-minus-H2O": {"mass": 162.052823432, 
                              "score": 1},
            # HexNAc minus one water:
            # C8H15NO6 - H2O = C8H13NO5
            "HexNAc-minus-H2O": {"mass": 203.0793725337,
                                 "score": 1},
            # Hex-Hex minus one water:
            # 2 C6H12O6 - 2 H2O = C12H20O10
            "Hex-Hex-minus-H2O": {"mass": 324.1056468639,
                                  "score": 1},
            # Neu5Ac - 2 water:
            "Neu5Ac-minus-2H2O": {"mass": 273.0848518422,
                                  "score": 1},
            # Neu5Ac - 1 water:
            "Neu5Ac-minus-H2O": {"mass": 291.0954165286,
                                  "score": 1},
            # Hex-HexNAc minus one water:
            # C6H12O6 + C8H15NO6 - 2 H2O = C14H23NO10
            "Hex-HexNAc-minus-H2O": {"mass": 365.1321959657,
                                     "score": 1},
            # Hex-Hex-HexNAc minus one water:
            # 2 C6H12O6 + C8H15NO6 - 3 H2O = C20H33NO15
            "Hex-Hex-HexNAc-minus-H2O": {"mass": 527.1850193976,
                                         "score": 1},
            # Hex-Hex-HexNAc-HexNAc minus one water:
            # 2 C6H12O6 + 2 C8H15NO6 - 4 H2O = C28H46NO20
            "Hex-Hex-HexNAc-HexNAc-minus-H2O": {"mass": 730.2643919313,
                                                "score": 1},
            # Hex-Hex-HexNAc-HexNAc:
            # 2 C6H12O6 + 2 C8H15NO6 - 3 H2O = C28H48NO21
            "Hex-Hex-HexNAc-HexNAc": {"mass": 748.2749566177,
                                      "score": 1}
            }

        if self.config["reduced"] == False:
            
            #-------------------- SCORE 5: --------------------#

            # If we see a peptide in the fragment
            #      that is highly indicating
            # So we give a high score for peptide-related fragment
            
            # protonated peptide chain:
            self.fragmentPoolDict["peptide"] = {
                "mass": self.matchedPeptideMass + waterMass,
                "score": 5}

            self.fragmentPoolDict["peptide-HexNAc-Y"] = {
                "mass": self.matchedPeptideMass + hexNAcMass,
                "score": 5}
                
            self.fragmentPoolDict["Peptide-HexNac-DeoxyHex-Y"] = {
                "mass": self.matchedPeptideMass + hexNAcMass + 146.0579088,
                "score":5}

            self.fragmentPoolDict["peptide-HexNAc-HexNAc-Y"] = {
                "mass": self.matchedPeptideMass + HexNAcHexNAcMinusOneWaterMass,
                "score": 5}
                
            self.fragmentPoolDict["peptide-HexNAc-HexNAc-DeoxyHex-Y"] = {
                "mass": self.matchedPeptideMass + HexNAcHexNAcMinusOneWaterMass + 146.0579088,
                "score": 5}    

            self.fragmentPoolDict["peptide-HexNAc-HexNAc-Hex-Y"] = {
                "mass": self.matchedPeptideMass + HexNAcHexNAcHexMinusTwoWaterMass,
                "score": 5}
            
            self.fragmentPoolDict["peptide-HexNAc-HexNAc-Hex-DeoxyHex-Y"] = {
                "mass": self.matchedPeptideMass + HexNAcHexNAcHexMinusTwoWaterMass + 146.0579088,
                "score": 5}

            self.fragmentPoolDict["peptide-HexNAc-HexNAc-Hex-Hex-Y"] = {
                "mass": self.matchedPeptideMass + HexNAcHexNAcHexHexMinusThreeWaterMass,
                "score": 5}
                
            self.fragmentPoolDict["peptide-HexNAc-HexNAc-Hex-Hex-DeoxyHex-Y"] = {
                "mass": self.matchedPeptideMass + HexNAcHexNAcHexHexMinusThreeWaterMass + 146.0579088,
                "score": 5}

            self.fragmentPoolDict["peptide-HexNAc-HexNAc-Hex(Hex)2-Y"] = {
                "mass": self.matchedPeptideMass + glycanCoreMass,
                "score": 5}
                
            self.fragmentPoolDict["peptide-HexNAc-HexNAc-Hex(Hex)2-DeoxyHex-Y"] = {
                "mass": self.matchedPeptideMass + glycanCoreMass + 146.0579088,
                "score": 5}
            
            #sialated n-linked
            if (massMatch.glycanCompositionMatched["Neu5Ac"] >= 1 
                    and massMatch.glycosylationSiteAminoAcid == "N"):
                
                self.fragmentPoolDict["glycopeptide-minus-Neu5Ac"] = {
                    "mass":self.matchedGlycanMass + self.matchedPeptideMass - 291.0954165286,
                    "score":5}
                    
                self.fragmentPoolDict["glycopeptide-minus-Hex-Neu5Ac"] = {
                    "mass":self.matchedGlycanMass + self.matchedPeptideMass - 453.1482399606,
                    "score":5}
                    
                # complex sialated
                if massMatch.glycanCompositionMatched["Neu5Ac"] >= 1 and massMatch.glycanCompositionMatched["HexNAc"] >= 3:
                        
                    self.fragmentPoolDict["glycopeptide-minus-HexNac-Hex-Neu5Ac"] = {
                        "mass":self.matchedGlycanMass + self.matchedPeptideMass - 656.2276124942,
                        "score":5}
                
                # trisialated or greater
                if massMatch.glycanCompositionMatched["Neu5Ac"] >= 2 and massMatch.glycanCompositionMatched["HexNAc"] >= 3:
                    
                    self.fragmentPoolDict["glycopeptide-minus-2xNeu5Ac"] = {
                        "mass":self.matchedGlycanMass + self.matchedPeptideMass - 291.0954165286 - 291.0954165286,
                        "score":5}
                        
                    self.fragmentPoolDict["glycopeptide-minus-Hex-2xNeu5Ac"] = {
                        "mass":self.matchedGlycanMass + self.matchedPeptideMass - 453.1482399606 - 291.0954165286,
                        "score":5}
                        
                    self.fragmentPoolDict["glycopeptide-minus-HexNac-Hex-2xNeu5Ac"] = {
                        "mass":self.matchedGlycanMass + self.matchedPeptideMass - 656.2276124942 - 291.0954165286,
                        "score":5}
                        
                    self.fragmentPoolDict["glycopeptide-minus-2xHex-2xNeu5Ac"] = {
                        "mass":self.matchedGlycanMass + self.matchedPeptideMass - 453.1482399606 - 453.1482399606,
                        "score":5}
                        
                    self.fragmentPoolDict["glycopeptide-minus-HexNac-2xHex-2xNeu5Ac"] = {
                        "mass":self.matchedGlycanMass + self.matchedPeptideMass - 656.2276124942 - 453.1482399606,
                        "score":5}
                    
            
            # high mannose glycans       
            elif (massMatch.glycanCompositionMatched["Hex"] >= 5
                    and massMatch.glycanCompositionMatched["HexNAc"] == 2 
                    and massMatch.glycosylationSiteAminoAcid == "N"):
                
                self.fragmentPoolDict["glycopeptide-minus-Hex"] = {
                    "mass":self.matchedGlycanMass + self.matchedPeptideMass - 162.052823432,
                    "score":5 }
                    
                if massMatch.glycanCompositionMatched["Hex"] >= 6 and massMatch.glycanCompositionMatched["HexNAc"] == 2:
                    self.fragmentPoolDict["glycopeptide-minus-2Hex"] = {
                        "mass":self.matchedGlycanMass + self.matchedPeptideMass - 324.1056468639,
                        "score":5 }
                
                if massMatch.glycanCompositionMatched["Hex"] >= 8 and massMatch.glycanCompositionMatched["HexNAc"] == 2:
                    
                    self.fragmentPoolDict["glycopeptide-minus-2Hex-Hex"] = {
                        "mass":self.matchedGlycanMass + self.matchedPeptideMass - 324.1056468639 - 162.052823432,
                        "score":5 }
                    
                    self.fragmentPoolDict["glycopeptide-minus-Hex-Hex"] = {
                        "mass":self.matchedGlycanMass + self.matchedPeptideMass - 162.052823432 - 162.052823432,
                        "score":5 }
                    
                    self.fragmentPoolDict["glycopeptide-minus-2Hex-2Hex"] = {
                        "mass":self.matchedGlycanMass + self.matchedPeptideMass - 324.1056468639 - 324.1056468639,
                        "score":5}
            
            #small complex glycans without sialation
            elif (massMatch.glycanCompositionMatched["HexNAc"] >= 4 
                    and massMatch.glycanCompositionMatched["Neu5Ac"] == 0 
                    and massMatch.glycosylationSiteAminoAcid == "N"):
                
                self.fragmentPoolDict["glycopeptide-minus-HexNAc"] = {
                    "mass":self.matchedGlycanMass + self.matchedPeptideMass - 203.0793725337,
                    "score":5}
                
                if massMatch.glycanCompositionMatched["Hex"] >= 4 and massMatch.glycanCompositionMatched["HexNAc"] >= 4 and massMatch.glycanCompositionMatched["Neu5Ac"] == 0:
                    
                    self.fragmentPoolDict["glycopeptide-minus-Hex"] = {
                        "mass":self.matchedGlycanMass + self.matchedPeptideMass - 162.052823432,
                        "score":5}
                        
                    self.fragmentPoolDict["glycopeptide-minus-Hex-HexNAc"] = {
                        "mass":self.matchedGlycanMass + self.matchedPeptideMass - 365.1321959657,
                        "score":5 }
                        
        if massMatch.glycanCompositionMatched["Neu5Ac"] != 0:

            #-------------------- SCORE 3: --------------------#

            # If Neu5Ac is spotted in the composition of the glycan in the matching
            #     We give Neu5Ac-related fragment score 3

            # Sialic Acid (Neu5Ac) species

            # Neu5Ac minus two waters:
            # C11H19NO9 - 2 H2O = C11H15NO7
            self.fragmentPoolDict["Neu5Ac-minus-2(H2O)"] = {
                "mass": 273.0848518422,
                "score": 3}
            # Neu5Ac minus one water:
            # C11H19NO9 - H2O = C11H17NO8
            self.fragmentPoolDict["Neu5Ac-minus-H2O"] = {
                "mass": 291.0954165286,
                "score": 3}
            # Neu5Ac-Hex minus one water:
            # C11H19NO9 + C6H12O6 - 2 H2O = C17H27NO13
            self.fragmentPoolDict["Neu5Ac-Hex-minus-H2O"] = {
                "mass": 453.1482399606,
                "score": 3}
            # Neu5Ac-Hex:
            # C11H19NO9 + C6H12O6 - H2O = C17H29NO14
            self.fragmentPoolDict["Neu5Ac-Hex"] = {
                "mass": 471.158804647,
                "score": 3}
            # Neu5Ac-Hex-HexNAc minus one water:
            # C11H19NO9 + C6H12O6 + C8H15NO6 - 3 H2O = C25H40N2O18
            self.fragmentPoolDict["Neu5Ac-Hex-HexNAc-minus-H2O"] = {
                "mass": 656.2276124942,
                "score": 3}

        if massMatch.glycanCompositionMatched["Neu5Gc"] != 0:

            #-------------------- SCORE 3: --------------------#

            # If Neu5Gc is spotted in the composition of the glycan in the matching
            #     We give Neu5Gc-related fragment score 3
            
            # Neu5Gc species

            # Neu5Gc minus two waters:
            # C11H19NO10 - 2 H2O = C11H15NO8
            self.fragmentPoolDict["Neu5Gc-minus-2(H2O)"] = {
                "mass": 289.0797664643,
                "score": 3}
            # Neu5Gc minus one water:
            # C11H19NO10 - H2O = C11H17NO9
            self.fragmentPoolDict["Neu5Gc-minus-H2O"] = {
                "mass": 307.0903311507,
                "score": 3}
            # Neu5Gc-Hex minus one water:
            # C11H19NO10 + C6H12O6 - 2 H2O = C17H27NO14
            self.fragmentPoolDict["Neu5Gc-Hex-minus-H2O"] = {
                "mass": 469.1431545827,
                "score": 3}
            # Neu5Gc-Hex:
            # C11H19NO10 + C6H12O6 - H2O = C17H29NO15
            self.fragmentPoolDict["Neu5Gc-Hex"] = {
                "mass": 487.1537192691,
                "score": 3}
            # Neu5Gc-Hex-HexNAc minus one water:
            # C11H19NO10 + C6H12O6 + C8H15NO6 - 3 H2O = C25H40N2O19
            self.fragmentPoolDict["Neu5Gc-Hex-HexNAc-minus-H2O"] = {
                "mass": 672.2225271164,
                "score": 3}

        if massMatch.glycanCompositionMatched["DeoxyHex"] != 0:

            #-------------------- SCORE 3: --------------------#

            # If DeoxyHex is spotted in the composition of the glycan in the matching
            #     We give DeoxyHex-related fragment score 3
            # Fucose species
            
            # HexNAc-DeoxyHex minus one water
            # C8H15NO6 + C6H12O5 - 2 H2O = C14H23NO9
            self.fragmentPoolDict["HexNAc-DeoxyHex-minus-H2O"] = {
                "mass": 349.1372813435,
                "score": 3}
            # Hex-HexNAc-DeoxyHex minus one water:
            # C6H12O6 + C8H15NO6 + C6H12O5 - 3 H2O = C20H33NO14
            self.fragmentPoolDict["Hex-HexNAc-DeoxyHex-minus-H2O"] = {
                "mass": 511.1901047755,
                "score": 3}
            # Hex-Hex-HexNAc-HexNAc-DeoxyHex minus one water:
            # 2 C6H12O6 + 2C8H15NO6 + C6H12O5 - 5 H2O = C34H56N2O24
            self.fragmentPoolDict["Hex-Hex-HexNAc-HexNAc-DeoxyHex-minus-H2O"] = {
                "mass": 876.3223007412,
                "score": 3}
            # Hex-Hex-HexNAc-HexNAc-DeoxyHex:
            # 2 C6H12O6 + 2C8H15NO6 + C6H12O5 - 4 H2O = C34H58N2O25
            self.fragmentPoolDict["Hex-Hex-HexNAc-HexNAc-DeoxyHex"] = {
                "mass": 894.3328654276,
                "score": 3}

        if massMatch.glycanCompositionMatched["Pentose"] != 0:

            #-------------------- SCORE 3: --------------------#

            # If Pentose (e.g. Xylose) is spotted in the composition of the glycan in the matching
            #     We give Pentose-related fragment score 3
            # Pentose species
            '''
            # Pentose-Hex:
            # C5H10O5 + C6H12O6 - H2O = C11H20O10
            self.fragmentPoolDict["Pentose-Hex"] = {
                "mass": 312.1056468639,
                "score": 3}
            '''
            # Pentose-Hex minus one water:
            # C5H10O5 + C6H12O6 - 2H2O = C11H18O9dictionary
            self.fragmentPoolDict["Pentose-Hex-minus-H2O"] = {
                "mass": 294.0950821776,
                "score": 3}
            '''
            # Pentose-Hex-HexNAc:
            # C5H10O5 + C6H12O6 + C8H15NO6 - 2H2O = C19H33NO15
            self.fragmentPoolDict["Pentose-Hex-HexNAc"] = {
                "mass": 515.1850193976,
                "score": 3}
            '''
            # Pentose-Hex-HexNAc minus one water:
            # C5H10O5 + C6H12O6 + C8H15NO6 - 3H2O = C19H31NO14
            self.fragmentPoolDict["Pentose-Hex-HexNAc-minus-H2O"] = {
                "mass": 497.1744547112,
                "score": 3}
            '''
            # Pentose-Hex-Hex:
            # C5H10O5 + 2C6H12O6 - 2H2O = C17H30O15
            self.fragmentPoolDict["Pentose-Hex-Hex"] = {
                "mass": 474.1584702959,
                "score": 3}
            '''
            # Pentose-Hex-Hex minus one water:
            # C5H10O5 + 2C6H12O6 - 3H2O = C17H28O14
            self.fragmentPoolDict["Pentose-Hex-Hex-minus-H2O"] = {
                "mass": 456.1479056095,
                "score": 3}
            '''
            # Pentose-Hex(Hex)-HexNAc:
            # C5H10O5 + 2C6H12O6 + C8H15NO6 - 3H2O = C25H43NO20
            self.fragmentPoolDict["Pentose-Hex(Hex)-HexNAc"] = {
                "mass": 677.2378428296,
                "score": 3}
            '''
            # Pentose-Hex(Hex)-HexNAc minus one water:
            # C5H10O5 + 2C6H12O6 + C8H15NO6 - 4H2O = C25H41NO19
            self.fragmentPoolDict["Pentose-Hex(Hex)-HexNAc-minus-H2O"] = {
                "mass": 659.2272781432,
                "score": 3}
            '''
            # Pentose-Hex(Hex)2:
            # C5H10O5 + 3C6H12O6 - 3H2O = C23H40O20
            self.fragmentPoolDict["Pentose-Hex(Hex)2"] = {
                "mass": 636.2112937279,
                "score": 3}
            '''
            # Pentose-Hex(Hex)2 minus one water:
            # C5H10O5 + 3C6H12O6 - 4H2O = C23H38O19
            self.fragmentPoolDict["Pentose-Hex(Hex)2-minus-H2O"] = {
                "mass": 618.2007290415,
                "score": 3}
            '''
            # Pentose-Hex-HexNAc-HexNAc:
            # C5H10O5 + C6H12O6 + 2C8H15NO6 - 3H2O = C27H46N2O20
            self.fragmentPoolDict["Pentose-Hex-HexNAc-HexNAc"] = {
                "mass": 718.2643919313,
                "score": 3}
            '''
            # Pentose-Hex-HexNAc-HexNAc minus water:
            # C5H10O5 + C6H12O6 + 2C8H15NO6 - 4H2O = C27H44N2O19
            self.fragmentPoolDict["Pentose-Hex-HexNAc-HexNAc-minus-H2O"] = {
                "mass": 700.2538272449,
                "score": 3}
            '''
            # Pentose-Hex(Hex)-HexNAc-HexNAc:
            # C5H10O5 + 2C6H12O6 + 2C8H15NO6 - 4H2O = C33H56N2O25
            self.fragmentPoolDict["Pentose-Hex(Hex)-HexNAc-HexNAc"] = {
                "mass": 880.3172153633,
                "score": 3}
            '''
            # Pentose-Hex(Hex)-HexNAc-HexNAc minus one water:
            # C5H10O5 + 2C6H12O6 + 2C8H15NO6 - 5H2O = C33H54N2O24
            self.fragmentPoolDict["Pentose-Hex(Hex)-HexNAc-HexNAc-minus-H2O"] = {
                "mass": 862.3066506769,
                "score": 3}
            '''
            # Pentose-Hex(Hex)2-HexNAc:
            # C5H10O5 + 3C6H12O6 + C8H15NO6 - 4H2O = C31H53NO25
            self.fragmentPoolDict["Pentose-Hex(Hex)2-HexNAc"] = {
                "mass": 839.2906662616,
                "score": 3}
            '''
            # Pentose-Hex(Hex)2-HexNAc minus one water:
            # C5H10O5 + 3C6H12O6 + C8H15NO6 - 5H2O = C31H51NO24
            self.fragmentPoolDict["Pentose-Hex(Hex)2-HexNAc-minus-H2O"] = {
                "mass": 821.2801015752,
                "score": 3}
            '''
            # Pentose-Hex(Hex)2-HexNAc-HexNAc: (NlinkedCore-Xylose)
            # C5H10O5 + 3C6H12O6 + 2C8H15NO6 - 5H2O = C39H66N2O30
            self.fragmentPoolDict["Pentose-Hex(Hex)2-HexNAc-HexNAc"] = {
                "mass": 1042.3700387953,
                "score": 3}
            '''
            # Pentose-Hex(Hex)2-HexNAc-HexNAc minus one water:
            # C5H10O5 + 3C6H12O6 + 2C8H15NO6 - 6H2O = C39H64N2O29
            self.fragmentPoolDict["Pentose-Hex(Hex)2-HexNAc-HexNAc-minus-H2O"] = {
                "mass": 1024.3594741089,
                "score": 3}

        if (self.config["free"] == False):
            if (massMatch.glycosylationSiteAminoAcid != "N"
                    and massMatch.glycosylationSiteAminoAcid != "N^phos'ed"):

                # For O-linked glycosylation sites
                #     We have at least 2 more additional diagnostic ions
                self.fragmentPoolDict["peptide-Hex"] = {
                    "mass": self.matchedPeptideMass + hexMass,
                    "score": 5}

                self.fragmentPoolDict["peptide-Hex-HexNAc"] = {
                    "mass": self.matchedPeptideMass + HexHexNAcMinusOneWaterMass,
                    "score": 5}

                if massMatch.glycanCompositionMatched["Neu5Ac"] != 0:

                    self.fragmentPoolDict["peptide-HexNac-Neu5Ac"] = {
                        "mass": self.matchedPeptideMass + 291.0954165286 + 203.079372533,
                        "score": 5}

                    self.fragmentPoolDict["peptide-Hex-HexNAc-Neu5Ac-minus-H2O"] = {
                        "mass": self.matchedPeptideMass + HexHexNAcNeu5AcMinusThreeWaterMass,
                        "score": 5}

                    self.fragmentPoolDict["peptide-Hex-HexNAc-Neu5Ac"] = {
                        "mass": self.matchedPeptideMass + HexHexNAcNeu5AcMinusTwoWaterMass,
                        "score": 5}

                if massMatch.glycanCompositionMatched["Neu5Gc"] != 0:

                    self.fragmentPoolDict["peptide-HexNac-Neu5Gc"] = {
                        "mass": self.matchedPeptideMass + 307.0903311507 + 203.079372533,
                        "score": 5}

                    self.fragmentPoolDict["peptide-Hex-HexNAc-Neu5Gc-minus-H2O"] = {
                        "mass": self.matchedPeptideMass + HexHexNAcNeu5GcMinusThreeWaterMass,
                        "score": 5}

                    self.fragmentPoolDict["peptide-Hex-HexNAc-Neu5Gc"] = {
                        "mass": self.matchedPeptideMass + HexHexNAcNeu5GcMinusTwoWaterMass,
                        "score": 5}

        #find backbone fragments
        peptideList = massMatch.peptideStringMatched.split("-")
        peptideBegin = massMatch.peptideLeftIndexIncluded
        glycanSite = massMatch.glycosylationSiteIndex
        
        #calculateIonMasses
        aaMassDict = self.config["aminoAcidMassDict"]
        #for aminoAcidIndex in range(leftIndexIncluded, rightIndexIncluded + 1):
        #    aminoAcidAtThisIndex = aminoAcidSequenceDict[aminoAcidIndex]
        #    aminoAcidMassAtThisIndex = aminoAcidDict[str(aminoAcidAtThisIndex)]["Mass"]
        #    peptideMass = peptideMass + aminoAcidMassAtThisIndex
        
        #begin stuff from elseware
        relativePosition = glycanSite - peptideBegin
        #print relativePosition
        ptmDict = {relativePosition:{"name":"glycan", "mass":self.massMatch.glycanMass-18.010565}}   #a list of locations and ptm annotated data

        def peptideBackboneFragmentSequence( sequence, ptmDict, usedIons = "by"):
            """
            this function generates a(name, mass, sequence) tuple for peptide backbone fragmentation
            
            test using: http://db.systemsbiology.net:8080/proteomicsToolkit/FragIonServlet.html
            """
            seqLen = len(sequence)
            
            #the following logic creates a list of modifications mass deltas mapped to AA sequence location
            allMods = [0 for a in range(len(sequence))]
            for i in range(len(sequence)):
                if i in ptmDict:
                    allMods[i] = ptmDict[i]["mass"]
            realSeq = [(sequence[i], allMods[i]) for i in range(seqLen)]
            #print realSeq
            
            #the following adjustments will make a singly charged product
            #y(right; c side) and b(left; n side)
            yAdj = 19.018390 #(+water+H)
            bAdj = 1.0072765 #(+1 proton)
            
            #x(right; c side) and a(left; n side)
            xAdj = 44.997655 #(+COOH)
            aAdj = -26.98709 #(- CO + H)
            
            #z(right; c side) and c(left; n side)
            zAdj = -15.010899 #(-NH)
            cAdj = 17.026549 #(+NH3)
            
            ionTypeDict = {"y":yAdj, "b":bAdj, "a":aAdj, "x":xAdj, "z":zAdj  , "c":cAdj }
            ionDict = {}
            ionNames = []
            for ionType in usedIons.lower():
                if ionType in "abc":
                    for i in range(seqLen-1):
                        ionNames.append(ionType+str(i+1))
                        ionDict[ionType+str(i+1)] = realSeq[0:i+1] + [("",ionTypeDict[ionType])]
                elif ionType in "xyz":
                    for i in range(seqLen-1):
                        ionNames.append(ionType+str(i+1))
                        ionDict[ionType+str(i+1)] = realSeq[seqLen-i-1:seqLen] + [("",ionTypeDict[ionType])]
            
            for name in ionNames:
                aa = ionDict[name]
                seq = "".join([b[0] for b in aa])
                name = name + "-backbone-" + seq
                mass = sum(aaMassDict[b[0]]["Mass"]+b[1] for b in aa) - 1.0072763
                yield (name, mass, seq)
                    
        for fragName, fragMass, sequence in peptideBackboneFragmentSequence(peptideList, ptmDict, usedIons="by"):
            
            self.fragmentPoolDict[fragName] = {
                    "mass": fragMass,
                    "score": 3}
            #print fragName
            #print sequence
            #print self.fragmentPoolDict[fragName]
            #raw_input(" --- waiting -----")
                
        
            
    def score_feature_object(self, feature):

        # Function Description:
        #     Based on the fragment pool
        #         that contains all the diagnostic ion fragments to score
        #     We give a total score for every parent match
        #         where the total score is summed up by all fragment scores
        
        fragmentPoolDict = self.fragmentPoolDict
        
        fragmentDecisionDeviation = self.config["scoringTolerance"]
        fragmentDecisionDeviationFromMassTolerance = self.config["scoringToleranceMass"]
        xyzPeakList = feature.xyzPeakList

        listOfFragmentNeutralMassMatches = []
        
        totalScore = 0

        hydrogenMass = 1.0078250321
        
        '''listOfTentativeFragmentMasses = [mass
                                         for mass in self.fragmentPoolDict.values()]
        listOfTentativeFragmentMasses.sort()'''

        # I avoided both repetitive list comprehension and dictionaries,
        # to make the code faster
        numberOfPeptideIonFragmentMatches = 0
        numberOfGlycanIonFragmentMatches = 0
        numberOfDeoxyHexIonFragmentMatches = 0
        numberOfNeu5AcIonFragmentMatches = 0
        numberOfNeu5GcIonFragmentMatches = 0
        numberOfCommonIonFragmentMatches = 0
        numberOfPepBackboneFragmentMatches = 0
            
        for matchedFragmentName in fragmentPoolDict.keys():
            tentativeMass = fragmentPoolDict[matchedFragmentName]["mass"]
            fragmentLowerLimitMass = tentativeMass * (1 - fragmentDecisionDeviation)
            fragmentUpperLimitMass = tentativeMass * (1 + fragmentDecisionDeviation)
            fragmentLowerLimitMassDueToToleranceMass = tentativeMass - fragmentDecisionDeviationFromMassTolerance
            fragmentUpperLimitMassDueToToleranceMass = tentativeMass + fragmentDecisionDeviationFromMassTolerance
            
            for xyzPeak in xyzPeakList:
                mzFromXyzPeak = xyzPeak[0]
                chargeFromXyzPeak = xyzPeak[2]

                if chargeFromXyzPeak > 0:
                    fragmentActualNeutralMass = (mzFromXyzPeak * chargeFromXyzPeak
                                                 - hydrogenMass * chargeFromXyzPeak)
                else:
                    fragmentActualNeutralMass = mzFromXyzPeak
                    
                if (fragmentActualNeutralMass >= fragmentLowerLimitMass and
                        fragmentActualNeutralMass <= fragmentUpperLimitMass) or \
                   (fragmentActualNeutralMass >= fragmentLowerLimitMassDueToToleranceMass and
                        fragmentActualNeutralMass <= fragmentUpperLimitMassDueToToleranceMass):

                    multiplierFromXyzPeak = xyzPeak[3]
                    matchedTheoreticalFragmentNeutralMass = fragmentPoolDict[matchedFragmentName]["mass"]
                    deviationFromActualToTheoreticalFragmentNeutralMass = ((
                            fragmentActualNeutralMass
                            - tentativeMass) / tentativeMass * 1000000)
                    deviationFromActualToTheoreticalFragmentNeutralMass = str(
                            round(deviationFromActualToTheoreticalFragmentNeutralMass,
                                  2)) + " ppm"
                    scoreContributed = fragmentPoolDict[matchedFragmentName]["score"] * multiplierFromXyzPeak
                    
                    if "times" in fragmentPoolDict[matchedFragmentName].keys():
                        fragmentPoolDict[matchedFragmentName]["times"] += 1
                    else:
                        fragmentPoolDict[matchedFragmentName]["times"] = 1

                    if "peptide" in matchedFragmentName:
                        numberOfPeptideIonFragmentMatches += 1
                    elif "glycan" in matchedFragmentName:
                        numberOfGlycanIonFragmentMatches += 1
                    elif "DeoxyHex" in matchedFragmentName:
                        numberOfDeoxyHexIonFragmentMatches += 1
                    elif "Neu5Ac" in matchedFragmentName:
                        numberOfNeu5AcIonFragmentMatches += 1
                    elif "Neu5Gc" in matchedFragmentName:
                        numberOfNeu5GcIonFragmentMatches += 1
                    elif "backbone" in matchedFragmentName:
                        numberOfPepBackboneFragmentMatches += 1
                    else:
                        numberOfCommonIonFragmentMatches += 1
                        
                    fragmentNeutralMassMatch = FragmentNeutralMassMatch(
                            xyzPeak,
                            matchedFragmentName,
                            matchedTheoreticalFragmentNeutralMass,
                            fragmentActualNeutralMass,
                            mzFromXyzPeak,
                            deviationFromActualToTheoreticalFragmentNeutralMass,
                            scoreContributed)

                    listOfFragmentNeutralMassMatches.append(fragmentNeutralMassMatch)
                    
                    #end the iteration through the fragment spectra if a match is found
                    #not having this break may over-promote some spectra
                    break
                    
                    '''elif (fragmentActualNeutralMass >= fragmentLowerLimitMassDueToToleranceMass and
                        fragmentActualNeutralMass <= fragmentUpperLimitMassDueToToleranceMass):

                    multiplierFromXyzPeak = xyzPeak[3]
                    matchedTheoreticalFragmentNeutralMass = fragmentPoolDict[matchedFragmentName]["mass"]
                    deviationFromActualToTheoreticalFragmentNeutralMass = ((
                            fragmentActualNeutralMass
                            - tentativeMass) / tentativeMass * 1000000)
                    deviationFromActualToTheoreticalFragmentNeutralMass = str(
                            round(deviationFromActualToTheoreticalFragmentNeutralMass,
                                  2)) + " ppm"
                    scoreContributed = fragmentPoolDict[matchedFragmentName]["score"] / multiplierFromXyzPeak
                    
                    if "times" in fragmentPoolDict[matchedFragmentName].keys():
                        fragmentPoolDict[matchedFragmentName]["times"] += 1
                    else:
                        fragmentPoolDict[matchedFragmentName]["times"] = 1

                    if "peptide" in matchedFragmentName:
                        numberOfPeptideIonFragmentMatches += 1
                    elif "glycan" in matchedFragmentName:
                        numberOfGlycanIonFragmentMatches += 1
                    elif "DeoxyHex" in matchedFragmentName:
                        numberOfDeoxyHexIonFragmentMatches += 1
                    elif "Neu5Ac" in matchedFragmentName:
                        numberOfNeu5AcIonFragmentMatches += 1
                    elif "Neu5Gc" in matchedFragmentName:
                        numberOfNeu5GcIonFragmentMatches += 1
                    else:
                        numberOfCommonIonFragmentMatches += 1
                        
                    fragmentNeutralMassMatch = FragmentNeutralMassMatch(
                            xyzPeak,
                            matchedFragmentName,
                            matchedTheoreticalFragmentNeutralMass,
                            fragmentActualNeutralMass,
                            mzFromXyzPeak,
                            deviationFromActualToTheoreticalFragmentNeutralMass,
                            scoreContributed)

                    listOfFragmentNeutralMassMatches.append(fragmentNeutralMassMatch)'''
                        
        '''totalScore = sum([matchableComponentDict["score"]
                          for matchableComponentDict in fragmentPoolDict.values()
                          if "times" in matchableComponentDict])'''

        totalNeu5AcInMatch = self.massMatch.glycanCompositionMatched["Neu5Ac"]
        totalNeu5GcInMatch = self.massMatch.glycanCompositionMatched["Neu5Gc"]          
        totalDeoxyHexInMatch = self.massMatch.glycanCompositionMatched["DeoxyHex"]
        
        
        if totalDeoxyHexInMatch == 0 and numberOfDeoxyHexIonFragmentMatches >= 1:
            totalScore -= 5*numberOfDeoxyHexIonFragmentMatches
        if totalNeu5GcInMatch == 0 and numberOfNeu5GcIonFragmentMatches >= 1:
            totalScore -= 5*numberOfNeu5GcIonFragmentMatches
        if totalNeu5AcInMatch == 0 and numberOfNeu5AcIonFragmentMatches >= 1:
            totalScore -= 5*numberOfNeu5AcIonFragmentMatches
                          
                          
                          
        totalScore += sum([fragmentNeutralMassMatch.scoreContributed
                          for fragmentNeutralMassMatch
                          in listOfFragmentNeutralMassMatches])
        totalScore = max(totalScore, 0)
        
        return (listOfFragmentNeutralMassMatches,
                totalScore,
                numberOfPeptideIonFragmentMatches,
                numberOfGlycanIonFragmentMatches,
                numberOfDeoxyHexIonFragmentMatches,
                numberOfNeu5AcIonFragmentMatches,
                numberOfNeu5GcIonFragmentMatches,
                numberOfCommonIonFragmentMatches,
                numberOfPepBackboneFragmentMatches)


    
# the list of fragmentNeutralMassMatch
# is going to be one massMatch class attribute

class FragmentNeutralMassMatch:

    def __init__(self,
                 matchedPeakTuple,
                 matchedFragmentName,
                 matchedTheoreticalFragmentNeutralMass,
                 matchedActualFragmentNeutralMass,
                 matchedActualFragmentMz,
                 deviationFromActualToTheoreticalFragmentNeutralMass,
                 scoreContributed
                 ):

        # Setting all parameters passed into this class to be member variables
        
        self.matchedPeakTuple = matchedPeakTuple
        self.matchedFragmentName = matchedFragmentName
        self.matchedTheoreticalFragmentNeutralMass = matchedTheoreticalFragmentNeutralMass
        self.matchedActualFragmentNeutralMass = matchedActualFragmentNeutralMass
        self.matchedActualFragmentMz = matchedActualFragmentMz
        self.deviationFromActualToTheoreticalFragmentNeutralMass = deviationFromActualToTheoreticalFragmentNeutralMass
        self.scoreContributed = scoreContributed
        


    def printFragmentMatchDetails(self):

        # Function Description:
        #     Print out the details of the fragment mass matches
        
        print("===============")
        
        print("Fragment matches: \t\t"
              + self.matchedFragmentName)

        print("The Fragment mz in Matched Peak: \t\t\t\t"
              + str(self.matchedPeakTuple[0]))
        print("The Fragment Intensity in Matched Peak: \t\t\t\t"
              + str(self.matchedPeakTuple[1]))
        print("The Fragment Charge in Matched Peak: \t\t\t\t"
              + str(self.matchedPeakTuple[2]))
        
        print("The Matched Theoretical Fragment Neutral Mass: \t\t\t"
              + str(self.matchedTheoreticalFragmentNeutralMass))
        print("The Matched Actual Fragment Neutral Mass: \t\t\t"
              + str(self.matchedActualFragmentNeutralMass))
        print("The Deviation from Actual to Theoretical Fragment Neutral Mass: \t"
              + str(self.deviationFromActualToTheoreticalFragmentNeutralMass))
        print("The Matched Actual Fragment MZ: \t\t\t\t\t"
              + str(self.matchedActualFragmentMz))
        if self.scoreContributed == 4:
            print("The Score Contributed --------- GLYCAN -------- : "
                  + str(self.scoreContributed))
            
        else:
            print("The Score Contributed: "
                  + str(self.scoreContributed))
        '''print("The Score Contributed: \t"
              + str(self.scoreContributed))'''
