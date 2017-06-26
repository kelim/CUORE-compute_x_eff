#!/usr/bin/python
import os

###################################
### Multiplicity cut
n_multiplicity = 1

## in case n_multiplicity >1
TotalEnergy = 500 #keV

### Histogram binsize
h_energy_binsize = 0.5 #keV 

### Determine peak count range (n_peak count range)
n_sigma_for_peak_count_stack = [5, 6]
n_sigma_for_peak_count_index_min = n_sigma_for_peak_count_stack.index(5)
n_sigma_for_peak_count_index_max = n_sigma_for_peak_count_stack.index(5)+1

### Determine background count range, excluding n_sigma_for_peak_count*2
n_sigma_for_bg_count  = 15

###################################

### peak_info
peak_name_stack = ["Te_X-rays","511-keV","K-40","Tl-208","Pt-190","Po-210"]
peak_energy_min_stack = [0, 510-60, 1460-20, 2615-60, 3290-80, 5000] #keV
peak_energy_max_stack = [50, 510+60, 1460+20, 2615+60, 3290+80, 5800] #keV

if len(peak_name_stack)!=len(peak_energy_min_stack):
    print "Error! check the peak info stacks!!"
elif len(peak_energy_max_stack)!=len(peak_energy_min_stack):
    print "Error! check the peak info stacks!!"

peak_index_min = peak_name_stack.index("K-40")
peak_index_max = peak_name_stack.index("K-40")+1

### energy estimator info
e_estimator_stack = ["Energy","EnergyNewOF","EnergyDecorrOF","EnergyWoH","EnergyNewOFWoH","EnergyDecorrOFWoH"]

e_estimator_index_min = e_estimator_stack.index("Energy")
e_estimator_index_max = e_estimator_stack.index("Energy")+1

#print e_estimator_stack[0]

### channel info
ch_stack = [2,3,4,6,7,9,12,14,17,18,19,20,21,24,25,26,31,32,36,41,44,46]

ch_index_min = ch_stack.index(2) 
ch_index_max = ch_stack.index(2)+1
#ch_index_max = len(ch_stack)
#print ch_index_max

### dataset info
ds_stack = ["2079","2085","2088","2091","2097","2100","2103","2109","2124","2130","2139","all"]
ds_index_min = ds_stack.index("2079")
ds_index_max = ds_stack.index("2079")+1

#ds_index_min = ds_stack.index("all")
#ds_index_max = len(ds_stack)

ds_index_all = ds_stack.index("all")

#print ds_index_max

### data type info
# check the list name in the processed data folder
data_type_stack = ["calib","background","physics"]
n_data_type_min = data_type_stack.index("calib")
n_data_type =data_type_stack.index("background")+1 

### Reduced files?
b_reduced = 1

### file names to submit eff computation scripts
script_running_dir = 'job_scripts'
cast_exefile_name = script_running_dir+'/cast_compute_x_eff.C'
cast_readfiles_name = script_running_dir+'/cast_readfiles.C'
cast_qsubfile_name = script_running_dir+'/cast_submit_jobs.pbs'
dir_path = os.path.dirname(os.path.realpath(cast_exefile_name)) + '/'
#print dir_path

### Insert above parameters in the root script
positions = []
job_id = 0
for ik in range (ch_index_min,ch_index_max):
    for jk in range (ds_index_min,ds_index_max):
        for j in range(e_estimator_index_min,e_estimator_index_max):
            for k in range(peak_index_min,peak_index_max):
                for l in range(n_sigma_for_peak_count_index_min,n_sigma_for_peak_count_index_max):
                    for lk in range(n_data_type_min,n_data_type):
                        job_id = job_id+1
                        readfiles_name = 'readfiles_%d.C' %(job_id)
                        exefile_name = 'compute_x_eff_%d.C' %(job_id)
                        qsubfile_name = 'submit_jobs_%d.pbs' %(job_id)
                        positions.append([k, peak_name_stack[k], peak_energy_min_stack[k],peak_energy_max_stack[k], 
                                          j, e_estimator_stack[j], n_multiplicity,
                                          h_energy_binsize, ch_stack[ik],
                                          ds_stack[jk], ds_index_all, jk,
                                          data_type_stack[lk],
                                          n_sigma_for_peak_count_stack[l],
                                          n_sigma_for_bg_count,
                                          TotalEnergy,b_reduced,script_running_dir,readfiles_name,exefile_name,job_id])
                        
print positions 
s_readfiles = os.path.abspath(cast_readfiles_name)+' > '+ dir_path 
#+ readfiles_name
s_exe = os.path.abspath(cast_exefile_name)+' > '+ dir_path 
#+ exefile_name
s_qsub = os.path.abspath(cast_qsubfile_name)+' > '+ dir_path
#+ qsubfile_name

#### fill in the parameters 
for position in positions:
    #common paramters for readfiles and compute_x_eff
    set_parameters = "sed -e \'s/_peak_index_/%d/g; s/_peak_name_/%s/g; s/_energy_min_/%d/g; s/_energy_max_/%d/g; s/_draw_branchIndex_/%d/g; s/_draw_branch_/%s/g; s/_n_multiplicity_/%d/g; s/_h_energy_binsize_/%f/g; s/_ch_number_/%d/g; s/_ds_name_/%s/g; s/_ds_index_all_/%d/g; s/_ds_index_/%d/g; s/_data_type_/%s/g; s/_n_sigma_for_peak_count_/%d/g; s/_n_sigma_for_bg_count_/%d/g; s/_tot_energy_threshold_/%f/g; s/_b_reduced_/%d/g; s/_script_running_dir_/%s/g; s/_readfiles_name_/%s/g; s/_exefile_name_/%s/g\' "
    #cmd1 to read files
    set_readfiles_parameters = set_parameters+s_readfiles+'readfiles_%d.C'
    cmd1 = set_readfiles_parameters %tuple(position)
    os.system(cmd1)
#    print cmd1
    #cmd2 to compute cut eff
    set_exe_parameters = set_parameters+s_exe+'compute_x_eff_%d.C'
    cmd2 = set_exe_parameters %tuple(position)
#    print cmd2
    os.system(cmd2)
    #cmd3 to write qsub file to execute the script
    set_qsub_parameters = set_parameters+s_qsub+'submit_jobs_%d.pbs'
    cmd3 = set_qsub_parameters %tuple(position)
#    print cmd3
    os.system(cmd3)
    #cmd 4 to submit batch jobs contain root -q -b compute_x_eff.C
    cmd4 = 'qsub %s' % (dir_path) + 'submit_jobs_%d.pbs' %(tuple(position)[-1])
#    print cmd4
    os.system(cmd4)

    
