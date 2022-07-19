# How to install PLANETOPLOT?
-----------------------------

This install guide is mostly for Linux/UNIX OS users. The following install steps should also work on Mac OSX provided you could access its underlying UNIX system and install `python` librairies. The tool has never been tested on Windows (but should be functional provided `python` requirements are met).

## Important! Getting the environment variables right

Add the `planetoplot/modules` folder to `PYTHONPATH` in your environment file (e.g. `.bashrc`)

	export PYTHONPATH=$PYTHONPATH:adapt_to_your_own/planetoplot/modules

If you plan to use the command line tools, add the `planetoplot/bin` folder to your `PATH`

	export PATH=$PATH:adapt_to_your_own/planetoplot/bin

Do not forget to source the environment file.

## Required librairies

- `python` (2.6 or 2.7)
- `numpy`
- `matplotlib`
- `netCDF4` (not required if `ppplot` only is used)

*Note* : You might want/need to use Python 3+, in this case checkout [python3 branch](https://github.com/aymeric-spiga/planetoplot/tree/python3) but all might not be working correctly. Currently master is set to work on 2.7.

## Recommended librairies

- `basemap` (for mapping)
- `scipy` (for scientific computations)

## Quick guide to install the necessary `python` librairies

A suite of pre-packaged `python` librairies such as [Anaconda](https://store.continuum.io/cshop/anaconda) is probably the best way to get started with a complete python environment.

### Possibility 1: full `anaconda` suite

The full `anaconda` suite is free of charge and easy to install

 1. Go to [this website](https://store.continuum.io/cshop/anaconda/)
 2. Install `anaconda` on your system (Python 2.7 version)
 3. Add the `anaconda` `bin` directory to your `PATH` environment variable
 4. Add `basemap` and `netCDF4` with the shell command `conda install basemap pil netCDF4`

### Possibility 2: add librairies one-by-one with `miniconda`

The only drawback of `anaconda` is that you download a pretty huge amount of librairies. To make the installation lighter -- should you ever want to -- use the `miniconda` software

 1. Go to [this website](https://conda.io/miniconda.html)
 2. Download the appropriate `miniconda` install script (Python 2.7 version)
 3. Install `miniconda` by executing this script
 4. Install the necessary packages using the `conda install` command
~~~
conda install numpy
conda install matplotlib
conda install netCDF4
~~~
(If needed, install the optional packages)

- `basemap` for map projections
- `pillow` for background image (warpimage in `basemap`)
- `scipy` for scientific computations

~~~
conda install -c conda-forge basemap
conda install pil
conda install scipy
~~~

Here is a summary of conda commands using a dedicated environment for planetoplot
~~~
conda create -n planetoplot python=2.7
source activate planetoplot
conda install numpy matplotlib netcdf4
conda install -c conda-forge basemap
conda install pillow
~~~

### Possibility 3: install directly librairies to the built-in `python` in your OS 

If you know what you are doing and have the administration rights to do so.
