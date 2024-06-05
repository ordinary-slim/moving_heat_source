3-dimensional LPBF example. Results featured in 3rd example of paper.
As the heat source moves to the right,
the material of the elements switches from bulk to powder.
Beware mesh size is around 250K els. To reduce computational time consider
increasing the radius of the heat source or decreasing the dimensions of the domain.
Run as:

```
python3 main.py --run-reference --layers=1 --case-name=ref# Runs reference model and stores post under post_ref
python3 main.py --run-coupled --layers=1 --case-name=coupled# Runs Chimera model
```
