# Scalar field dynamics around a Kerr black hole

A C++ (version C++11 onwards) code that solves the equations of 
motion for a real scalar field with a potential around a Kerr black hole.
See the Documenation for more details about the equations of motion.

Runtime parameters are configured in the `runscript.py` file
(with some preprocessing done in `sim_params.py`).

To compile all the different executables I recommend using
`config.py`.

## Libraries

See `INSTALL.md` for more information on installing some
of the below libraries.

* [SHTns](https://nschaeff.bitbucket.io/shtns/)
	Fast spherical harmonic transforms.

* [FFTW](http://fftw.org) 
	Optional (for the pseudospectral radial derivative option).

* [OpenMP](https://www.openmp.org/)
	Optional (and yes, not really a library).

* [googletest](https://github.com/google/googletest)
	Optional (this is only used to run unit tests under `Test/`).
	
## Visualization

I've mostly used [paraview](https://www.paraview.org/) to visualize the data.
Contact me for more details about how to, e.g.
read in and visualize the 3d data.
For the moment all data are saved as `.csv` files. 

## Citation

Exact citation a work in progress,
but please do contact me if you'd like to use/cite this code.

## Contact

For questions please contact
Justin Ripley: ripley [at] illinois [dot] edu
