MNE-Python scripts for the evaluation of EEG/MEG spatial resolution.
These scripts are still under development, and are not intended to
be used out-of-the-box yet. If you would like to use any of these
utilities, please contact olaf.hauk@mrc-cbu.cam.ac.uk.

The following scripts were used to produce the results in
Hauk, O., M. Stenroos, and M. Treder, Towards an objective Evaluation
of EEG/MEG Source Estimation Methods: The Linear Tool Kit. bioRxiv, 2019

WH_Res_ResolutionMetrics_bioRxiv0.py: computes resolution metrics for different methods,
sensor configurations, regularisation parameters, etc.
WH_Res_config_bioRxiv0.py: configuration file with parameters for WH_Res_ResolutionMetrics.py.
WH_Resolution_Function_bioRxiv0.py: Functions to compute resolution metrics etc.
WH_Res_Morph_STCs_bioRxiv0.py: Morph resolution metrics (STC files) to average brain.
WH_Res_AvgSTCs.py: Grand-average morphed STCs.
WH_Res_QSUB_bioRxiv0.py: Submit parallel jobs to computing cluster.

The following example scripts demonstrate the basic functionality,
see README_Resolution_Example.txt:

Resolution_Example.py
Resolution_Example_Functions.py