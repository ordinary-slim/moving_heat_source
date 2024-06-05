Computational environment for the simulation of convection-diffusion problems with moving heat sources. C++ library with Python interface.

![2D welding animation with coupled and reference method](https://github.com/ordinary-slim/media_moving_heat_source/blob/master/animation.gif)

In particular, it can be used for 3D printing problems. _TODO: add gif_

# Installation

## Docker

Pull the docker image and run an interactive container:
```
docker pull ordinaryslim/mhs_chimera
docker run -ti ordinaryslim/mhs_chimera
```

The container will spawn you in this folder after which you can try out the examples / tests,
e.g.

```
cd examples/3d_ded_unitary_params
python3 main.py --run-coupled --layers=1 --case-name=coupled# Runs Chimera model
```

Use `docker cp` to transfer the results to your machine after which you can open them using `Paraview`

## Source

Dependencies and instructions can be found in the Dockerfile

# Post-processing

I recommend examining the resulting `vtk`s using Paraview.

A Chimera model run produces two datasets, one for the background domain and one for the moving domain.
Opening both datasets on Paraview and thresholding out the inactive elements will produce the results on the computational domain.
Alternatively one can open only the background dataset and threshold out by physicalDomain in order to receover the computational domain.

Similarly, for the examples involving 3D printing, one needs to threshold out the inactive elements.
