**PLANETOPLOT: a cool python-based tool to plot stuff and explore data**
* To get sources through git `git clone https://github.com/aymeric-spiga/planetoplot`
* To get sources through SVN `svn co https://github.com/aymeric-spiga/planetoplot/trunk planetoplot`
* To get a static ZIP file of the current version of the code, [click here](https://github.com/aymeric-spiga/planetoplot/archive/master.zip)
* To install, see [**installation notes**](https://github.com/aymeric-spiga/planetoplot/blob/master/INSTALL.md)
* To look at examples and learn how to use, see [**tutorial**](http://nbviewer.ipython.org/github/aymeric-spiga/planetoplot/blob/master/tutorial/planetoplot_tutorial.ipynb)
* To share experience and examples, write in the [**wiki**](https://github.com/aymeric-spiga/planetoplot/wiki)

**Contents**

* `modules`
 - `ppclass.py`: main class with `pp()` objects 
 - `ppplot.py`: plot class
 - `ppcompute.py`: computation class

* `bin`
 - `pp.py`: main command line utility to use ppclass [HELP with `pp.py -h`]
 - `pp_reload.py`: command line utility to load saved plot objects `*.ppobj`
 - `asciiplot.py`: a simple script to easily plot data in a text file

* `settings`
 - `set_area.txt` --> setting file: predefined areas for plotting (can be omitted)
 - `set_back.txt` --> setting file: predefined backgrounds for plotting (can be omitted)
 - `set_multiplot.txt` --> setting file: predefined coefficients for multiplots (can be omitted)
 - `set_ppclass.txt` --> setting file: predefined variables for x,y,z,t (can be omitted)
 - `set_var.txt` --> setting file: predefined colorbars, format, labels, etc... for variables (can be omitted)

* `example`
 - example scripts using `ppclass`

* `tutorial`
 - tutorial (based on `ipython notebook`)

