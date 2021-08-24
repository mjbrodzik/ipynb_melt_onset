# ipynb_melt_onset
Notebooks for melt onset research.

## To install the python cetbtools on Windows:

1. Download and install Miniconda3 from

http://conda.pydata.org/miniconda.html.

Miniconda is a program that will allow you control what are
called conda "environments," which just means that you can set up
individualized python installations with difference sets of
packages. It does all of the package management for you, so it
will know what combination of packages to install for cetbtools,
and jupyter and all the tool you will want to use. Most
importantly, it keeps this installation completely separate from
any other tools you may have installed that use their own special
version of python.  So you will be playing "nice" and not
stepping on your other tools.

The latest Miniconda version currently is for Python 3.9, this is
fine.

You will download a file called something like:

Miniconda3-py39_4.9.2-Windows-x86_64.exe

Go ahead and run it to install Miniconda.

Accept all defaults during the installation, except for one:

Check to make sure you have a latest conda installed with

conda update -n base -c defaults conda

As of 4/6/2021, the latest version of conda is 4.10.0.
As of 7/13/2021, the latest version of conda is 4.10.3.

2. Get the miniconda shell window by doing:

Windows menu > All Programs > Anaconda Prompt (miniconda3).

3. Set up the list of "channels," which are locations to look
online for the packages we will need, and the priority order that
these locations will be searched. (If you had a previous
.condarc, it might need to be deleted or appended.) The first
conda config command that is run will always include defaults
channel.

Add necessary channels (order is important here):

conda config --add channels conda-forge
conda config --add channels nsidc

You should see the following when checking your channels:

conda config --show channels

channels:
   - nsidc
   - conda-forge
   - defaults

4. Create an environment with the cetbtools and jupyter

First, confirm that you only have the base environment:

conda info --envs

should return:

# conda environments:
base                 * C:<your home directory>\miniconda3

Now install a package called mamba in the base environment that
will speed up a lot of the rest of the process:

conda install mamba

Now use mamba to create a new environment that we will use for our project:

If you are on Windows:

mamba create -n cetb python~=3.7 cetbtools matplotlib scipy jupyter basemap proj4 seaborn

(Including proj4 corrects a problem with basemap 1.2 that doesn't properly
set the PROJ_LIB env variable see this for grizzly details:

https://stackoverflow.com/questions/52295117/basemap-import-error-in-pycharm-keyerror-proj-lib
)

If you are on Mac or linux:

mamba create -n cetb python~=3.7 cetbtools matplotlib scipy jupyter basemap seaborn

To get a very specific cetbtools:

Windows:
mamba create -n cetbTest python~=3.7 cetbtools==1.6.0.rc2 matplotlib scipy jupyter basemap proj4 seaborn
Non-Windows:
mamba create -n cetbTest python~=3.7 cetbtools==1.6.0.rc2 matplotlib scipy jupyter basemap seaborn

N.B.:
1) jupyter is needed for running jupyter notebooks
2) seaborn is a graphics package that does some things much much
better than matplotlib
3) basemap is needed for some of the displays and maps that we are
displaying in the notebooks.

N.B.: confirm that matplotlib and scipy are both coming from conda-forge?

If you run "conda env list" again will will see a second entry,
for cetb. The asterisk will be positioned in the environment that is
"current."

5. Switch to the new conda env that we just created:

conda activate cetb

Note that your command prompt has changed from (base) to (cetb).

This took a little while to run on my machine, be patient. Mine
generated some non-fatal warning messages prepended with
"SafetyError:..."
I think it's OK to ignore them, I'm checking on that.

6. Test the installation.  At the command prompt:

python -c "from cetbtools import algorithms; help(algorithms)"

should return an array of with methods in the algorithms module.

7. Start jupyter notebook. The directory where you start the
jupyter notebook server will be the top directory in the jupyter
notebooks, so I recommend cd'ing to it if it's a long way from
your home directory.

Optional:

cd <directory with your ipython notebooks>

At the command prompt, start the jupyter notebook server:

jupyter notebook

Your browser should start up with jupyter server running in new
tab. 

8. Open a blank notebook, and try doing:

import cetbtools

this should do the import without generating an error message.

help(cetbtools)

this should display the package help message

from cetbtools import ease2conv
help(ease2conv)

this should display the ease2conv help message

Congratulations, you have successfully installed the cetbtools
python package on your local machine.


## Miscellaneous Conda Notes

To list all the conda envs you have created:

conda info --envs

To remove a conda env:

conda remove --name <my_env> --all

