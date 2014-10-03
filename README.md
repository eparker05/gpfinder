Readme for GPFinder
=======

Glycopeptide Finder 3 main repository

### About
###### What is this
GPFinder 3 is a python script that implements an algorithm to assign
glycopeptide identities (peptide structure + glycan composition) by
their fragmentation mass spectra. 

###### Project Scope
While candidate glycopeptides can be found with a single tandem mass
spectrum, statistical significance of assignments as stated by a false
discovery rate can be determined using decoy analysis. As implemented
in GPF3, decoy analysis requires input from many spectra taken in a single
session, thus, in order to use decoy analysis, GPF3 requires data be
provided as with mass-chromatograms. Other tools exist to work with single
tandem spectra, I recommend glycopep grader made by Carrie Woodin of the
Desaire lab (http://glycopro.chem.ku.edu/GPGHome.php ).


### Setting up an environment for GPF3

Prior to using GPF3 install all dependencies and their subdependencies. 

GP3 depends on Python 2.7, Numpy, and MatPlotLib.

Experienced Python users may be familiar with the process of installing
these dependnecies but for others, I recommend using a scientific python
distribution like Enthought Python Distribution. Others exist that include
these libraries. Importantly, GPF3 does not yet work with python 3.

### GPF 3 quickstart

GP Finder 3 is currently best described as a script, rather than a full
fledged python library. As such, it is used by using git to clone this
repository then copying analysis files into the resources directory.

    ~$ git clone https://github.com/gpfinder/gpfinder.git
    ~$ cd gpfinder
     $ cp ~/mydata.mgf ./resources/mydata.mgf
     $ cp ~/protein_uniprot_record.xml ./resources/myprotein.xml

Now that your data accessable to gpfinder, open the config.xml, edit line
`proteinfile` to contain the name of your protein library and edit the line
`massspecfile` to point to your mgf data file. Now run the analysis.

     $ python run_gpfinder.py

Once complete, the results will be output in the various CSV files described
later.

##### Sample preparation & acquisition

Sample preparation is out of the scope of this guide, but many resources
exist in literature to cover glycopeptide generation and purification.
GPF3 has been tested on glycopeptides derived from tryptic digestion,
pronase digestion, and elastase digestion. Several other proteases are
suported as well as combinations of proteases. New proteases can be added
by copying an existing protease specificity function and implementing logic
to cover behavior of the desired protease.

##### Preparation of spectra

Once glycopeptide spectra have been collected and uploaded into a spectral
viewer, the desired spectra should be exported as a MGF or mzXML file.
NO filtering should be performed before-hand; for the best results use
GPF3's builtin filtering based on kernel density estimation of the noise
level within a spectra.

Once the spectra has been exported, move a copy to the GPF3 resources folder
where it can be used for a search.

##### Configuration of GPF3

GPF3 is configured with the file "GPF3Config.xml". This file should be
in the root directory of the analysis (along side the "run_gpf3.py").
Each configuration option is verbose and for ambigous portions, helper
text has been included as a comment. A future guide will cover this file
in depth but currently rely on the commented configuration file.

In addition to the configuration file, search resources are kept in the
"resources" folder.
 * Protein files supplied as uniprot XML files
 * Glycan libraries are a simple CSV based format containing compositions


### 
