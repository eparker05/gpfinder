Readme for GPFinder
=======

Glycopeptide Finder 3 main repository

### About
###### What is this
GPFinder 3 is a python script that implements an algorithm to assign
glycopeptide identities (peptide structure + glycan composition) by
their fragmentation mass spectra. 

###### Who made this
This software was written by Evan A. Parker and Michael Xin Sun. 
GP Finder 3 is based on GP Finder 2, GP Finder, and Glyco-X.
See http://chemgroups.ucdavis.edu/~lebrilla/ for publications
on these topics.

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

#### Sample preparation & acquisition

GPF3 has been tested on glycopeptides derived from tryptic digestion,
pronase digestion, and elastase digestion. Several other proteases are
suported as well as combinations of proteases. New proteases can be added
by copying an existing protease specificity function and implementing logic
to cover behavior of the desired protease. A detailed guide on sample
preparation is out of the scope of this guide, but many resources exist
in literature that cover glycopeptide production and purification.

#### Preparation of spectra

Once glycopeptide spectra have been collected and uploaded into a spectral
viewer, the desired spectra should be exported as a MGF or mzXML file.
NO filtering should be performed before-hand; for the best results use
GPF3's builtin filtering based on kernel density estimation of the noise
level within a spectra.

Once the spectra has been exported, move a copy to the GPF3 resources folder
where it can be used for a search.

#### Configuration of GPF3

GPF3 is configured with the file "GPF3Config.xml". This file should be
in the root directory of the analysis (along side the "run_gpfinder.py").
Each configuration option is verbose and for ambigous portions, helper
text has been included as a comment. A future guide will cover this file
in depth but currently rely on the commented configuration file. Currently,
the most important configuration options are described here.

In addition to the configuration file, search resources are kept in the
"resources" folder.
 * Protein files supplied as uniprot XML files
 * Glycan libraries are a simple CSV based format containing compositions

###### File resources:

GPF3 imports many resources for the aptly named directory `resources`.
The spectra to be analyzed should reside there along with uniprot xml
files describing the protein and a glycan library. If no protein xml file
exists, a sequence can be added maually later. 

The glycan library is a CSV containing rows describing the composition
of all included glycans. This library should reflect a targeted superset
of the expected glycoforms. A general purpose library is included with
over 300 structures found in mamalian cells. Finally, GPF3 can be run
without a glycan library using the combinatorial glycan option. This is
useful for exploratory searches but we recommend the use of a library.

###### Tags descriptive of MS1

After file tags, there are several tags used to describe the input spectra
to GPF3 at the MS1 level. Parts per million tollerance, maximum expected
charge (use 0 for pre-deconvoluted data), and the noise threshold for
significant data are all included here.

The tag `diagnosticfragmentfilter` is used to instruct GPF3 on the pre-
selection of data. Since kernel density estimation (the noise detection
algorithm) is slow, all spectra are pre-filtered for diagnostic oxonium
ions and neutral losses. Each spectra needs to include the specified count
of oxonium ions that sum to the specified fraction of the base peak.

###### Protein tag group

The next group of tags is bound by the `protein` tag. In this group the
digestion enzyme is specified and 

### Running gpfinder and output

Use your python interpreter to run the script run_gpfinder.py. Doing
this in an interactive shell will allow you to more easily detect errors.
Most configuration and file format errors should stop the run immediately.

During the run, three graphs will be rendered and once the run is
completed, two csv files will be created that contain the results.

CSV_1 contains the unfiltered candidate list for all ions with
glycopeptide candidates.

CSV_3 contains the results for each feature with a compound passing
the false discovery rate threshold of 5%.



