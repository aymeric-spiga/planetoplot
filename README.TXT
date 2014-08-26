
******************************************************************************

You might want to access planetoplot modules and programs from anywhere
This is done through the PYTHONPATH environment variable
We also recommend to add PYTHONPATH to your PATH to use scripts anywhere

In a nutshell:
Add to environment file (e.g. .bashrc) those lines

PYTHONPATH=$PYTHONPATH:you_got_to_adapt_this_to_your_own/planetoplot
export PYTHONPATH
PATH=$PYTHONPATH:$PATH
export PATH

Don't forget to source the environment file. Then enjoy!

******************************************************************************

----------------------------
PLANETOPLOT v2
----------------------------
A. Spiga LMD/UPMC 24/03/2013
----------------------------

Contents
--------
--------
core
----
- ppclass.py --> main class with pp() objects (need ppplot and ppcompute)
- ppplot.py --> plot class (can be used independently, need ppcompute)
- ppcompute.py --> computation class (can be used independently)
scripts
-------
- pp.py --> command line utility to use ppclass [HELP with 'pp.py -h']
- pp_reload.py --> command line utility to load saved plot objects *.ppobj
- example/* --> example scripts using ppclass
settings files
--------------
- set_area.txt --> setting file: predefined areas for plotting (can be omitted)
- set_back.txt --> setting file: predefined backgrounds for plotting (can be omitted)
- set_multiplot.txt --> setting file: predefined coefficients for multiplots (can be omitted)
- set_ppclass.txt --> setting file: predefined variables for x,y,z,t (can be omitted)
- set_var.txt --> setting file: predefined colorbars, format, labels, etc... for variables (can be omitted)
documentation
-------------
- README.TXT --> this README file
data
----
plot objects for a demonstration tour and customizing tests
https://dl.dropbox.com/u/11078310/demo_data.tar.gz

Requirements
------------
------------
python + numpy + matplotlib + netCDF4
- for mapping --> Basemap
- for scientific computations --> scipy
[recommended: Enthought Python Distribution (free for academics)]

Installation
------------
------------
- install required softwares and librairies in requirements
- add planetoplot_v2 in your PYTHONPATH environment variable (and in your PATH to use pp.py)

Take a demo tour
----------------
----------------
pp_reload.py demo_data/*

Improvements compared to v1
---------------------------
---------------------------
* code readability and structuration for future improvements
* modularity (class formulation)
  + easy definition/addition of attributes
* separation between data retrieval and plotting
* versatility
  + command line (pp.py) 
    --> for quick inspection
  + interactive session (ipython) 
    --> for testing and exploring
  + scripts 
    --> for powerful and fully customized results
* performance (overall and on large files)
  + memory consumption (only retrieve what is needed)
* saving/loading plot objects in/from *.ppobj
* plot aesthetics and customizing (see header in ppplot)
* multiplot flexibility with .plotin attribute
* easy definition of presets with set_*.txt files
* function: one field vs. another one
* emulation of + - / * ** operations on fields
  (between two fields or between a field and a int/float)
* computations of min/max in addition to mean
* simple inspection of dimensions+variables in a file (pp.py -f file)
* verbose / non verbose mode

Acknowledgements
----------------
----------------
Thanks to A. Colaitis, T. Navarro, J. Leconte
for feedbacks and contributions on version 1

