
"""
INPUT SETTINGS

"""

#Input file name
mgfFileName = "Peroxidase_elastase-1-29.mgf"

## Maximum charge of fragments
# set to 0 for deconvoluted data
fragmentMaxCharges = 2

#output file suffix
outputSuffix = "-filtered"


"""
noise filtering settings

"""

#diagnostic ions required
diagnosticIonCountRequired = 2

#diagnostic ion S/N ratio required for a count
diagnosticSignalToNoise = 5

# noise cutoff threshhold when doing KDE
noiseCutoffThreshold = 3 

## Parts per million error on fragments 
# for speed, ppm error comparisons will be approximate; 
# the absolute error is calculated from 600 m/z
ppmError = 40 



"""
#######
Program
#######
"""
import numpy
from collections import namedtuple
import random

#secret config options:
plotOnly = True

Feature = namedtuple("Feature", ["featureIndex",
                                 "compoundNeutralMass",
                                 "compoundMzValue",
                                 "compoundCharge",
                                 "compoundRtCenter",
                                 "compoundRtTuple",
                                 "compoundIntensity",
                                 "xyzPeakList"])

def parse_mgf(mgfFileName, noFilters = False):

    # Function Description
    #     Read the mgf file, which is the output from MassHunter
    # Function Input:
    #     The mgf file from MassHunter
    
    #setting charges list
    maximumFragmentCharge = fragmentMaxCharges
    if maximumFragmentCharge > 0:
        charges = [charge + 1 for charge in range(maximumFragmentCharge)]
    else:
        charges = [0]
        
    inFile = open(mgfFileName, 'r')
    
    
    featureList = []
    
    featureCounter = 0
    
    #used for protonated peptides
    adduct = 1.007276
    
    for line in inFile:
        
        line = line.strip()
        #print(line.split())
        if "BEGIN" in line and "IONS" in line:
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
            print("END IONS tag seen. End of Feature " + str(featureCounter) + "\n")
            
            if skipSpectrum == True:
                print("Spectrum skipped, not added to the list of features.\n")
            

            newFeature = Feature(featureCounter,
                                 compoundNeutralMass,
                                 compoundMzValue,
                                 compoundCharge,
                                 compoundRtCenter,
                                 None, # compoundRtTuple unknown
                                 None, # compoundIntensity unknown
                                 xyzPeakList )
            if noFilters:
                #don't bother filtering
                featureList.append(newFeature)
            else:
                #do the filtering
                defaultNoiseLevel = 100
                passingFilter, tempFeature, noise = plot_spectrum_and_filter_out_noisy_peaks_using_kernel_density_estimation( newFeature )
                if passingFilter:
                    newFeature = tempFeature
                    noiseLevel = noise
                    validByDiagnosticPeakChecker = check_for_diagnostic_ions(xyzPeakList, compoundNeutralMass, noiseLevel)
                
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
            xyzTuple = (x, y, 1, multiplier)
            xyzPeakList.append(xyzTuple)

    #used for testing:
    import time
    print("found %d features that have diagnostic ions out of %d features in file"%(len(featureList), featureCounter))
    time.sleep(2)
    
    return featureList



def plot_spectrum_and_filter_out_noisy_peaks_using_kernel_density_estimation( featureObject ):

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
    #     if successful:
    #       return Bool(True), featureObject(filtered), noise_level
    #     else:
    #       return Bool(False) None, None

    xyzPeakList = featureObject.xyzPeakList
    textDisplay = ("MS/MS for Feature " + str(featureObject.featureIndex)
                        + " for parent ion mz = " + str(featureObject.compoundMzValue)
                        + " at RT = " + str(featureObject.compoundRtCenter))
    newXyzPeakList = []
    successfulNoiseFiltering = True # assume success, failure will negate this
    lenXyzPeakList = len(xyzPeakList)
    
    
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
                xPointsKDE = numpy.linspace(0, maxIntensity, 2000)
                yPointsKDE = kde_function(xPointsKDE)
                maxKDE = numpy.max(yPointsKDE)
                max_inten_KDE = xPointsKDE[yPointsKDE.tolist().index(maxKDE)]
                
                # filter out xyzPeaks
                # that have a intensity being less than the cutoff intensity
                noiseCutoffIntensity = max_inten_KDE * float(noiseCutoffThreshold)
                
                newXyzPeakListEndSegment = [xyzPeak
                                           for xyzPeak in xyzPeakListShortCopy
                                           if xyzPeak[1] >= noiseCutoffIntensity]
                
                newXyzPeakList.extend(newXyzPeakListEndSegment)
                
                
                # PLOTTING FUNCTION
                #  uncomment the code block below to plot KDE
                """
                import matplotlib.pyplot as plt
                xPointsKDE_plot = numpy.linspace(0, maxIntensity, 200)
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
                ymax = numpy.max(intensityList)
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
                xticks = numpy.linspace(0,axHist.get_xlim()[1],3)
                axHist.set_xticks(xticks)
                axHist.get_yaxis().set_visible(False)
                
                # the vlines plot:
                intenZeros = numpy.zeros(len(mzList))
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
    
    
    if successfulNoiseFiltering:
        featureObject = featureObject._replace(xyzPeakList = newXyzPeakList)
        return successfulNoiseFiltering, featureObject, max_inten_KDE
    else:
        return successfulNoiseFiltering, None, None
        

def check_for_diagnostic_ions(xyzPeakList, compoundNeutralMass, noiseLevel):
    #given an xyzPeakList, this function will return true if the diagnostic peaks are present
    #  diagnostic peaks are found within tolerances used for scoring. If the sum of the signal
    #  of the diagnostic peaks is greater than the maximum peak multiplied by some adjustment factor 
    #  (usually a fraction), then this function will evaluate as 'True'
    #  In addition to checking intensity, a simple check is made such that at least 2 diagnostic peaks
    #  must be present for this routine to evaluate as True
    
        
    #setting diagnostic peaks, this depends on charge of input files:
    #   deconvoluted input files need different diagnostic ions
    maximumFragmentCharge = fragmentMaxCharges
    if maximumFragmentCharge > 0:
        diagnosticPeakList = [# Hex-minus-H2O-plus-H+:
                              # C6H12O6 - H2O + H = C6H11O5
                              163.0606484641,
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
                              (compoundNeutralMass - 162.052823432 + 2*1.0078250321)/2.0,
                              (compoundNeutralMass - 203.0793725337 + 2*1.0078250321)/2.0,
                              (compoundNeutralMass - 291.0954165286 + 2*1.0078250321)/2.0,
                              (compoundNeutralMass - 365.1321959657 + 2*1.0078250321)/2.0]
    else: 
        # same as masses from above, reduced by mass of one proton
        diagnosticPeakList = [162.052823432,
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
                              compoundNeutralMass - 365.1321959657]
    
    maxPeak = max(xyzPeakList, key=lambda xyzPeakList: xyzPeakList[1])
    maxPeakIntensity = maxPeak[1]
    scoreTolerance = 0.50*ppmError*(10**-6)
    scoreToleranceAbs= 600.0 * scoreTolerance  #a rough tolerance used for low mass fragments
    signalToNoiseRatioTest = noiseLevel*float(diagnosticSignalToNoise)
    
    scorePeakIntensity = 0
    scorePeakCount = 0
    scorePeakIndexList = []
    
    for diagnosticPeak in diagnosticPeakList:
        for i in range(len(xyzPeakList)):
            realPeakLower = xyzPeakList[i][0] - scoreToleranceAbs
            realPeakHigher = xyzPeakList[i][0] + scoreToleranceAbs
            if (realPeakLower <=  diagnosticPeak <= realPeakHigher
                    and xyzPeakList[i][1] > signalToNoiseRatioTest
                    and i not in scorePeakIndexList):
                scorePeakIndexList.extend(range(i, i+4))
                scorePeakIntensity += xyzPeakList[i][1]
                scorePeakCount += 1
                

    if scorePeakCount >= diagnosticIonCountRequired:
        return True
    else:
        print("Not enough diagnostic ions found, spectrum will not be appended")
        return False




#
# this modified version of the scipy.stats.kde.gaussian_kde() relies
#   only on numpy. simply using "import numpy as np" will enable any script
#   to use this class.
#
class gaussian_kde(object):
    """
    Representation of a kernel-density estimate using Gaussian kernels.



    Attributes
    ----------
    d : int
        number of dimensions
    n : int
        number of datapoints

    Methods
    -------
    kde.evaluate(points) : array
        evaluate the estimated pdf on a provided set of points
    kde(points) : array
        same as kde.evaluate(points)
    kde.resample(size=None) : array
        randomly sample a dataset from the estimated pdf.
    kde.covariance_factor() : float
        computes the coefficient that multiplies the data covariance matrix to
        obtain the kernel covariance matrix. Set this method to
        kde.scotts_factor or kde.silverman_factor (or subclass to provide your
        own). The default is scotts_factor.
    *For integration methods, see scipy distribution: scipy.stats.kde.gaussian_kde()

    Parameters
    ----------
    dataset : (# of dims, # of data)-array
        datapoints to estimate from

    """

    def __init__(self, dataset, bw_method="scott"):
        self.dataset = numpy.atleast_2d(dataset)
        
        if bw_method == "scott":
            self.covariance_factor = self.scotts_factor
        elif bw_method == "silverman":
            self.covariance_factor = self.silverman_factor
        else:
            self.covariance_factor = self.arbitrary_factor
            self.bw_method = float(bw_method)
            
        self.d, self.n = self.dataset.shape

        self._compute_covariance()
        
        
            
    def evaluate(self, points):
        """Evaluate the estimated pdf on a set of points.

        Parameters
        ----------
        points : (# of dimensions, # of points)-array
            Alternatively, a (# of dimensions,) vector can be passed in and
            treated as a single point.

        Returns
        -------
        values : (# of points,)-array
            The values at each point.

        Raises
        ------
        ValueError if the dimensionality of the input points is different than
        the dimensionality of the KDE.
        """

        points = numpy.atleast_2d(points).astype(self.dataset.dtype)

        d, m = points.shape
        if d != self.d:
            if d == 1 and m == self.d:
                # points was passed in as a row vector
                points = reshape(points, (self.d, 1))
                m = 1
            else:
                msg = "points have dimension %s, dataset has dimension %s" % (d,
                    self.d)
                raise ValueError(msg)

        result = numpy.zeros((m,), points.dtype)

        if m >= self.n:
            # there are more points than data, so loop over data
            for i in range(self.n):
                diff = self.dataset[:,i,numpy.newaxis] - points
                tdiff = numpy.dot(self.inv_cov, diff)
                energy = numpy.sum(diff*tdiff,axis=0)/2.0
                result += numpy.exp(-energy)
        else:
            # loop over points
            for i in range(m):
                diff = self.dataset - points[:,i,numpy.newaxis]
                tdiff = numpy.dot(self.inv_cov, diff)
                energy = numpy.sum(diff*tdiff,axis=0)/2.0
                result[i] = numpy.sum(numpy.exp(-energy),axis=0)

        result /= self._norm_factor

        return result

    __call__ = evaluate

    def resample(self, size=None):
        """Randomly sample a dataset from the estimated pdf.

        Parameters
        ----------
        size : int, optional
            The number of samples to draw.
            If not provided, then the size is the same as the underlying
            dataset.

        Returns
        -------
        dataset : (self.d, size)-array
            sampled dataset
        """

        if size is None:
            size = self.n

        norm = transpose(multivariate_normal(zeros((self.d,), float),
            self.covariance, size=size))
        indices = randint(0, self.n, size=size)
        means = self.dataset[:,indices]

        return means + norm


    def scotts_factor(self):
        return numpy.power(self.n, -1./(self.d+4))

    def silverman_factor(self):
        return power(self.n*(self.d+2.0)/4.0, -1./(self.d+4))
    
    def arbitrary_factor(self):
        return self.bw_method

    # This can be replaced with silverman_factor if one wants to use Silverman's
    # rule for choosing the bandwidth of the kernels.
    

    def _compute_covariance(self):
        """Computes the covariance matrix for each Gaussian kernel using
        covariance_factor
        """
        self.factor = self.covariance_factor()
        self.covariance = numpy.atleast_2d(numpy.cov(self.dataset, rowvar=1, bias=False) *
            self.factor * self.factor)
        self.inv_cov = numpy.linalg.inv(self.covariance)
        self._norm_factor = numpy.sqrt(numpy.linalg.det(2*numpy.pi*self.covariance)) * self.n
        


                                 
def write_new_mgf( featureList, fileName ):
    "this simple function writes all feature-type objects back to file"
    
    with open(fileName, 'wb') as outputHandle:
        for feature in featureList:
            outputHandle.write("BEGIN IONS\n")
            outputHandle.write("PEPMASS={}\n".format(feature.compoundMzValue))
            outputHandle.write("CHARGE={}+\n".format(feature.compoundCharge))
            outputHandle.write("TITLE={}\n".format("title_line_deprecated"))
            outputHandle.write("RTINSECONDS={}\n".format(feature.compoundRtCenter))
            for xyzPeak in feature.xyzPeakList:
                outputHandle.write("{}\t{}\n".format(xyzPeak[0],xyzPeak[1]))
            outputHandle.write("END IONS\n")

            
if __name__ == "__main__":


    outputFileName = mgfFileName.split(".")
    
    if len(outputFileName) > 1:
        outputFileName[-2] += outputSuffix
    else: 
        outputFileName[-1] += outputSuffix
    outputFileName = ".".join(outputFileName)
    
    
    if not plotOnly:
        featureList = parse_mgf( mgfFileName )
        write_new_mgf(featureList, outputFileName)
    else:
        import matplotlib.pyplot as plt
        featureList = parse_mgf( outputFileName, noFilters=True)
        plots = []
        for f in featureList:
            xs, ys, _junk1, _junk2 = zip(*f.xyzPeakList)
            yzeros = [0 for x in range(len(ys))]
            title = "m/z={} and z={}".format(str(f.compoundMzValue)[0:6], f.compoundCharge)
            newPlot = (xs, yzeros, ys, title)
            plots.append(newPlot)
            #plt.vlines(xs, yzeros, ys)
            #plt.title("m/z={} and z={}".format(str(f.compoundMzValue)[0:6], f.compoundCharge))
            #plt.show()
        curr_pos = 0
        def key_event(e):
            global curr_pos

            if e.key == "right" or e.key == "alt+right" :
                curr_pos = curr_pos + 1
            elif e.key == "left" or e.key == "alt+left" :
                curr_pos = curr_pos - 1
            else:
                print e.key
                return
            curr_pos = curr_pos % len(plots)

            ax.cla()
            plt.title(plots[curr_pos][3])
            ax.vlines(plots[curr_pos][0], plots[curr_pos][1], plots[curr_pos][2])
            fig.canvas.draw()

        fig = plt.figure()
        fig.canvas.mpl_connect('key_press_event', key_event)
        ax = fig.add_subplot(111)
        plt.title(plots[0][3])
        ax.vlines(plots[0][0], plots[0][1], plots[0][2])
        plt.show()
    
    


    
    