3-dimensional DED example. Fails in reducing the number of global time-steps because of T^n discontinuity issue discussed in the paper.
Run as:

```
python3 main.py --run-reference --layers=1 --case-name=ref# Simulates 1 layers with reference model and stores post under post_ref
python3 main.py --run-coupled --layers=1 --case-name=coupled# Runs Chimera model
```
