# Scalar field dynamics around a Kerr black hole

A c++(14) and python code that solves the equations of 
motion for a real scalar field around a Kerr black hole.

Look under **Releases** for the latest stable version of this code.

Runtime parameters are configured in the `config.py` file.

## Libraries

* FFTW: 
	http://fftw.org

* SHTns:
	https://nschaeff.bitbucket.io/shtns/

* OpenMP (optional): 
	https://www.openmp.org/

* googletest (optional; for unit testing; see Test/): 
	https://github.com/google/googletest

I have successfully compiled the code with g++

## Derivation of equations of motion in coordinate form

A Mathematica notebook that contains the equations of motion
(as described in the `code paper` listed under `Citation`) in coordinate
form can be found [here](https://github.com/JLRipley314/2nd-order-teuk-derivations).

## Visualization

I use pyqtgraph-graph derived software
(see [here](https://github.com/JLRipley314/sci-vis))
to visualize the data, which are saved as csv files. 

## Citation

## Contact

For questions please contact
Justin Ripley: lloydripley [at] gmail [dot] com
