# ipynb_melt_onset
Notebooks for melt onset research.

## To install the python cetbtools:

1. Download and install Miniconda3 from

http://conda.pydata.org/miniconda.html.

cd ~/Downloads

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh



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

As of 8/8/2022, the latest version of conda is 4.13.0.

2.a. Windows: Get the miniconda shell window by doing:

Windows menu > All Programs > Anaconda Prompt (miniconda3).

2.b. non-Windows: Open a terminal on your machine.

3. Set up the list of "channels," which are locations to look
online for the packages we will need, and the priority order that
these locations will be searched. (If you had a previous
.condarc, it might need to be deleted or appended.) The first
conda config command that is run will always include defaults
channel.

Add necessary channels (order is important here):

```
conda config --add channels conda-forge
conda config --add channels nsidc
```

You should see the following when checking your channels:

```
conda config --show channels
```

channels:
   - nsidc
   - conda-forge
   - defaults

4. Create an environment with the cetbtools and jupyter

First, confirm that you only have the base environment:

```
conda info --envs
```

should return:
```
# conda environments:
base                 * C:<your home directory>\miniconda3
```
Now install a package called mamba in the base environment that
will speed up a lot of the rest of the process:

```
conda install mamba
```

Now use mamba to create a new environment that we will use for our project:
```
mamba create -n cetb cetbtools matplotlib scipy jupyter basemap seaborn
```

(Starting summer 2024, pandas is complaining about future need to install
pyarrow for pandas 3.x and later. For now, add pyarrow to the conda env
manually. I think once pandas is pinned to 3.x, the pyarrow will come along with
it.)

Alternatively, some extra packages that are needed for geotiff clipping operations: 
```
mamba create -n cetb cetbtools matplotlib scipy jupyter basemap seaborn geopandas rasterio pycrs
```
This should get cetbtools 1.7.x--this is the latest as of August 2022.

Optional, only if it's needed to get a very specific version of cetbtools:

Windows:
```
mamba create -n cetb cetbtools matplotlib scipy jupyter basemap proj4 seaborn
```

Non-Windows:
```
mamba create -n cetb cetbtools matplotlib scipy jupyter basemap seaborn
```
as of August 2022, apparently mamba won't do the magic on OSX, use conda if working on a Mac.

N.B.:
1) jupyter is needed for running jupyter notebooks
2) seaborn is a graphics package that does some things much much
better than matplotlib
3) basemap is needed for some of the displays and maps that we are
displaying in the notebooks.

N.B.: confirm that matplotlib and scipy are both coming from conda-forge

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

jupyter notebook --ip=0.0.0.0 --no-browser

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

   
## Debugging the ipynb files and scripts

To use the ipynb debugger for understanding what the scripts are doing:

Add the following line to one of the script files in the place you want
to set a breakpoint:

``` python
import pdb; pdb.set_trace()
```

Save the python script.

Then (this is important) go back to the notebook, restart the kernel and clear
all outputs and then recompile all cells up to the point where you are calling
the function in the python script.  This will ensure that you are now importing
the script with the set_trace call in it. 

When you execute this cell, the python debugger will start an ipython debugger
("ipdb") subshell in the jupyter notebook.  To do things in this cell, you can
use the usual 

l : list the code here
l <line> : list the code centered on this <line>
the usual jupter commands to print variables, etc

Remember that this is a subshell to the jupyter notebook, so to execute commands
here you use <ret> (NOT shift-<ret>). 

To exit the ipdb subshell, use 

``` python
exit()
```


