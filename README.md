# Scalar field dynamics around a Kerr black hole

A C++ (version C++11 onwards) code that solves the equations of 
motion for a real scalar field around a Kerr black hole.

Runtime parameters are configured in the `runscript.py` file
(with some preprocessing done in `sim_params.py`).

To compile all the different executables I recommend using
`config.py`.

## Libraries

See `INSTALL.md` for more information on installing some
of the below libraries.

* [FFTW](http://fftw.org) 
	Fast Fourier transforms	

* [SHTns](https://nschaeff.bitbucket.io/shtns/)
	Fast spherical harmonic transforms.

* [OpenMP](https://www.openmp.org/)	

* [googletest](https://github.com/google/googletest)
	This is only used to run unit tests under `Test/`.
	
## Visualization

I've mostly used [paraview](https://www.paraview.org/) to visualize the data.
Contact me for more details about how to, e.g.
read in and visualize the 3d data.
For the moment all data are saved as `.csv` files. 

## Citation

## Contact

For questions please contact
Justin Ripley: lloydripley [at] gmail [dot] com
