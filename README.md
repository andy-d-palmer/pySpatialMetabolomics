# README #

Reference implementation for the SpatialMetabolomics processing pipeline

### What is this repository for? ###

This is a reference implementation of our pipeline for annotating high-resolution imaging mass spectrometry data. It is scientific code for demonstration purposes only.

### How do I get set up? ###

* Requires pyIMS and pyMS
    * clone from https://github.com/alexandrovteam/pyMS, https://github.com/alexandrovteam/pyIMS
    * add the local directory to the python path os.path.append(my/local//directory/python_code)
* install pyimzml
    * pip install pyimzml
* other dependencies:
    * h5py <for some data type support>
    * numpy
    * scipy
* Clone pySpatialMetabolomics
    * https://github.com/andy-d-palmer/pySpatialMetabolomics
    * add the local directory to the python path

### Processing a dataset ###
* Data should be in .imzml format
* The pipeline currently only supports centroid data
    * using centroided data is **highly** recommended for run time and annotation performance


### Contribution guidelines ###

* Run it, see what happens, ask Andy questions

### Who do I talk to? ###

palmer@embl.de