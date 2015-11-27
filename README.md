# README #

Reference implementation for the SpatialMetabolomics processing pipeline

### What is this repository for? ###

* ipython notebook script
* shows code for a single metabolite
* It's not optimised for run speed **at all**

### How do I get set up? ###
* Requires pyIMS and pyMS
    * clone from https://github.com/alexandrovteam/pyMS, https://github.com/alexandrovteam/pyIMS
    * add the local directory to the python path os.path.append(my/local//directory/python_code)
* install pyimzml
    * pip install pyimzml
* Clone pySpatialMetabolomics
    * https://github.com/andy-d-palmer/pySpatialMetabolomics
    * add the local directory to the python path

### Processing a dataset ###
* Data should be in .imzml format
* The pipeline currently only supports centroid data
    * using centroided data is **highly** recommended for run time and annotation performance
*
### Contribution guidelines ###

* Run it, see what happens, ask Andy questions

### Who do I talk to? ###

* palmer@embl.de