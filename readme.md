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

# Citation

If you found this library useful in academic or industry work, we appreciate your support if you consider:

1. Starring the project on Github
2. Citing the relevant paper(s):

[A Chimera method for thermal part-scale metal additive manufacturing simulation](https://doi.org/10.1016/j.finel.2024.104238)

``` bibtex
@article{slimani2024,
title = {A {C}himera method for thermal part-scale metal additive manufacturing simulation},
journal = {Finite Elements in Analysis and Design},
volume = {241},
pages = {104238},
year = {2024},
issn = {0168-874X},
doi = {https://doi.org/10.1016/j.finel.2024.104238},
url = {https://www.sciencedirect.com/science/article/pii/S0168874X2400132X},
author = {Mehdi Slimani and Miguel Cervera and Michele Chiumenti}
}
```

[Substepped and advected subdomain methods for part-scale LPBF modeling (*PREPRINT*)](http://dx.doi.org/10.2139/ssrn.5529518)
