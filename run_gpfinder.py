from GPFinder.pyFindMassMatch import *
from GPFinder.readConfigXml import configuration
from GPFinder.FeatureObject import *
from copy import deepcopy





config = configuration()
#config["filePrefix"] = testName + str(score).replace('.','_') # no file prefix by default, str(score) + "_"
#config["peptideScore"] = score      #for Testing
#config["dontMakePlots"] = True
#config["secret-write-mgf"] = "fetuin_dualM-highZslopes_interlaced.mgf"
#config["cidLower"] = True
#config["cidHigher"] = True
#config["make_test_isotope_plot"] = False 

### the following code tests all file-structures before doing the long search
listOfProteinFiles = config["listOfProteinFiles"]
if len(listOfProteinFiles) == 1:
    matchFinderObject = findMassMatch(config)
    del(matchFinderObject)
    # If there are more than one proteins
    # For keeping the interface of findMassMatch intact,
    # We replace the config file with the same one
    #     that has the protein file entry replaced by every single protein file
    # The new fake config is named fakeConfigForSingleProtein
elif len(listOfProteinFiles) > 1:
    for proteinFile in listOfProteinFiles:
        fakeConfigForSingleProtein = deepcopy(config)
        fakeConfigForSingleProtein["proteinFile"] = proteinFile
        matchFinderObject = findMassMatch(fakeConfigForSingleProtein)
        del(matchFinderObject)


# initialize feature object to extract feature data
features = GetFeatures(config)


# use feature object to find and score mass matches
features.find_mass_matches()

######### testFunction
features.rescore_based_on_crpytic_proteolytic_specificity_andor_ppm_distr(scoreIncreasePPC=1, pValue=0.05)

## output
features.write_matches_and_fragment_matches_to_csv()
#features.write_highest_score_only_matches_and_fragment_matches_to_csv_based_on_glycan_site_order()
cutoff, firstPasscutoff, passingPeptides = features.write_highest_score_only_matches_and_fragment_matches_to_csv_based_on_glycan_site_order(doCSV5 = True)
