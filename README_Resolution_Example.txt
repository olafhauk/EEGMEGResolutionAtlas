Example MNE-Python script to compute resolution metrics for linear inverse
operators.

Resolution_Example.py computes peak localisation error (PLE) and spatial
deviation (SD) for point-spread and cross-talk functions of L2-MNE and sLORETA.
PLE and SD distributions are saved in STC files, as well as the difference
distributions for MNE-sLORETA.

sLORETA has zero PLE, therefore outperforming L2-MNE.
However, L2-MNE outperforms sLORETA with respect to SD.

Resolution_Example_Functions.py contains functions used in
Resolution_Example.py.
Both files should be in the same directory.

Olaf Hauk, May 2019