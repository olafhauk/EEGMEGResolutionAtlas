#!/imaging/local/software/anaconda/latest/x86_64/bin/python
"""
==========================================
Submit qsub jobs for WH Resolution analysis
allow dependencies between jobs
==========================================

OH July 2018
"""

print(__doc__)

import subprocess
from os import path as op

# indices of subjects to process
subjs = range(0,16)

job_list = \
[
    #####
    # COMPUTING and MORPHING Sensitivity Maps
    #####
    # {'N':   'R_SMap',                                       # job name
    #  'Py':  'WH_SensitivityMaps',                         # Python script
    #  'Cf':  'WH_MNE_Resolution_config',                     # configuration script
    #  'Ss':  subjs,                                          # subject indices
    #  'mem': '2GB',                                          # memory for qsub process
    #  'dep': '',
    #  'var': 'SensMap RMS'} # variable string for python script (e.g. for filenames)},                                             # name of preceeding process (optional)
    # {'N':   'R_MphSM',                                      # job name
    #  'Py':  'WH_Res_MorphSTC',                              # Python script
    #  'Cf':  'WH_MNE_Resolution_config',                     # configuration script
    #  'Ss':  subjs,                                          # subject indices
    #  'mem': '2GB',                                          # memory for qsub process
    #  'dep': '',
    #  'var': '"SensitivityMaps SensMap SNR"'},               # additional variables for python script (optional)
    
    # #####
    # # COMPUTING Resolution Metrics
    # #####
    # {'N':   'R_locerr',                                    # job name
    #  'Py':  'WH_Res_ResolutionMetrics',                     # Python script
    #  'Cf':  'WH_MNE_Resolution_config',                     # configuration script
    #  'Ss':  subjs,                                          # subject indices
    #  'mem': '2GB',                                          # memory for qsub process
    #  'var': '"locerr peak"',
    #  'dep': ''},
    
    # {'N':   'R_width',                                    # job name
    #  'Py':  'WH_Res_ResolutionMetrics',                     # Python script
    #  'Cf':  'WH_MNE_Resolution_config',                     # configuration script
    #  'Ss':  subjs,                                          # subject indices
    #  'mem': '2GB',                                          # memory for qsub process
    #  'var': '"width sd"',
    #  'dep': ''},
    
    # {'N':   'R_amp',                                    # job name
    #  'Py':  'WH_Res_ResolutionMetrics',                     # Python script
    #  'Cf':  'WH_MNE_Resolution_config',                     # configuration script
    #  'Ss':  subjs,                                          # subject indices
    #  'mem': '2GB',                                          # memory for qsub process
    #  'var': '"amplitude peak"',
    #  'dep': ''}, 
    
    # #####
    # # MORPHING Resolution Metrics
    # #####
    # {'N':   'R_MphLocErr',                                     # job name
    #  'Py':  'WH_Res_MorphSTC',                              # Python script
    #  'Cf':  'WH_MNE_Resolution_config',                     # configuration script
    #  'Ss':  subjs,                                          # subject indices
    #  'mem': '2GB',                                          # memory for qsub process
    #  'dep': 'R_locerr',
    #  'var': '"ResolutionMetrics locerr_peak"'},
    
    # {'N':   'R_MphWidth',                                     # job name
    #  'Py':  'WH_Res_MorphSTC',                              # Python script
    #  'Cf':  'WH_MNE_Resolution_config',                     # configuration script
    #  'Ss':  subjs,                                          # subject indices
    #  'mem': '2GB',                                          # memory for qsub process
    #  'dep': 'R_width',
    #  'var': '"ResolutionMetrics width_sd"'},
    
    # {'N':   'R_MphAmp',                                     # job name
    #  'Py':  'WH_Res_MorphSTC',                              # Python script
    #  'Cf':  'WH_MNE_Resolution_config',                     # configuration script
    #  'Ss':  subjs,                                          # subject indices
    #  'mem': '2GB',                                          # memory for qsub process
    #  'dep': 'R_amp',
    #  'var': '"ResolutionMetrics amplitude_peak"'}
    
    #####
    # AVERAGING Resolution Metrics
    #####
    ### The following depend on completion of previous jobs for ALL subjects
    {'N':   'R_AvgMphLocErr',                                     # job name
     'Py':  'WH_Res_AvgSTCs',                               # Python script
     'Cf':  'WH_MNE_Resolution_config',                     # configuration script
     'Ss':  [''],                                           # subject indices, '' if across all subjects
     'mem': '2GB',                                          # memory for qsub process
     'dep': '',
     'var': '"ResolutionMetrics locerr_peak"'},
    
    {'N':   'R_AvgMphWidth',                                     # job name
     'Py':  'WH_Res_AvgSTCs',                               # Python script
     'Cf':  'WH_MNE_Resolution_config',                     # configuration script
     'Ss':  [''],                                           # subject indices, '' if across all subjects
     'mem': '2GB',                                          # memory for qsub process
     'dep': '',
     'var': '"ResolutionMetrics width_sd"'},
    
    {'N':   'R_AvgMphAmp',                                     # job name
     'Py':  'WH_Res_AvgSTCs',                               # Python script
     'Cf':  'WH_MNE_Resolution_config',                     # configuration script
     'Ss':  [''],                                           # subject indices, '' if across all subjects
     'mem': '2GB',                                          # memory for qsub process
     'dep': '',
     'var': '"ResolutionMetrics amplitude_peak"'}
]

# directory where python scripts are
dir_py = op.join('/', 'home', 'olaf', 'MEG', 'WakemanHensonEMEG', 'ScriptsResolution')

# directory for qsub output
dir_qsub = op.join('/', 'home', 'olaf', 'MEG', 'WakemanHensonEMEG', 'ScriptsResolution', 'qsub_out')

# wrapper to run python script via qsub
fname_wrap = op.join('/', 'home', 'olaf', 'MEG', 'WakemanHensonEMEG', 'wrapper_qsub_python_v015.sh')


# keep track of qsub Job IDs
Job_IDs = {}

for job in job_list:

    for Ss in job['Ss']:

        Ss = str(Ss) # turn into string for filenames etc.

        N = job['N']
        Py = op.join(dir_py, job['Py'])
        Cf = job['Cf']
        mem = job['mem']

        # files for qsub output
        file_out = op.join(dir_qsub, job['N'] + '_' + Cf + '-%s.out' % str(Ss))
        file_err = op.join(dir_qsub, job['N'] + '_' + Cf + '-%s.err' % str(Ss))

        # if job dependent of previous job, get Job ID and produce command
        if 'dep' in job: # if dependency on previous job specified            
            if job['dep']=='':
                dep_str = ''
            else:
                job_id = Job_IDs[job['dep'], Ss]
                dep_str = '-W depend=afterok:%s' % (job_id)
        else:
            dep_str = ''

        if 'node' in job: # if node constraint present (e.g. Maxfilter)
            node_str = job['node']
        else:
            node_str = ''

        if 'var' in job: # if variables for python script specified
            var_str = job['var']
        else:
            var_str = ''

        # use wrapper to submit python script
        qsub_cmd = 'qsub %s \
                        -N %s%s \
                        -l walltime=24:00:00,mem=%s \
                        -o %s \
                        -e %s \
                        -v pycmd="%s.py %s",subj_idx=%s,var=%s \
                        %s \
                        %s' \
                        % (fname_wrap, N, Ss, mem, file_out, file_err, Py, Cf, Ss, var_str, dep_str, node_str)

        # format string for display
        print_str = qsub_cmd.replace(' ' * 25, '  ')
        print('\n%s\n' % print_str)

        # execute qsub command
        proc = subprocess.Popen(qsub_cmd, stdout=subprocess.PIPE, shell=True)

        # get linux output
        (out, err) = proc.communicate()

        # keep track of Job IDs from qsub
        Job_IDs[N, Ss] = out.split('.')[0]

# Done