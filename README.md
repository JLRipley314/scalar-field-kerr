# Scalar field dynamics around a Kerr black hole

A C++(14) and python code that solves the equations of 
motion for a real scalar field around a Kerr black hole.

Look under **Releases** for the latest stable version of this code.

Runtime parameters are configured in the `config.py` file.

## Libraries

* FFTW (Fast Fourier Transforms): 
	http://fftw.org

* SHTns (Fast spherical harmonic transforms):
	https://nschaeff.bitbucket.io/shtns/

* OpenMP (optional): 
	https://www.openmp.org/

* googletest (only used to run unit tests under `Test/`): 
	https://github.com/google/googletest

I have successfully compiled the code with 
the GNU C++ compiler g++.

## Visualization

I use pyqtgraph-graph derived software
(see [here](https://github.com/JLRipley314/sci-vis))
to visualize the data, which are saved as csv files. 

## Citation

## Contact

For questions please contact
Justin Ripley: lloydripley [at] gmail [dot] com
