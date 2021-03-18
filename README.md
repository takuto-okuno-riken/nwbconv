# nwbconv
 NWB format convert tool
 
## Requirements: software
* MATLAB R2019a or later

## setup matnwb
~~~
!git clone https://github.com/NeurodataWithoutBorders/matnwb.git
cd matnwb
addpath(genpath(pwd));
generateCore();
~~~

## Command line tool
~~~
>> nwbconv
no input files. please specify nwb registration xlsx files.
usage: nwbconv [options] filename.xlsx ...
  --outpath           output files path (default:".")
  -v, --version       show version number
  -h, --help          show command line help
~~~

## Command line tool Demo
~~~
>> nwbconv ori004_ch1_rp.xlsx
success to load data : ori004_ch1_rp.mat:stack
success to load data : logfile_M36_2p_ori004.mat:TS.stim_log
success to load data : logfile_M36_2p_ori004.mat:TS.start_time
success to load data : logfile_M36_2p_ori004.mat:TS.stop_time
success to load data : logfile_M36_2p_ori004.mat:TS.ori_deg
export NWB file : ./ori004_ch1_rp.nwb
~~~

  