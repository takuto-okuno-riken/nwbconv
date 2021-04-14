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
  --session num       output session range (i.e. "1","1:8")
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

~~~
>> nwbconv --session 1:8 M73.xlsx
success to load data : M73/200508/Tiff/raw-256/M73-200508-1.tif
success to load data : M73/200508/Task-Log/task.mat:task_data(1).Cue,Rew,X,Y,tY,tX,tYs,tXs,pupil_X,pupil_Y,eyelid_X,eyelid_Y,lick,ccdA,ccdB
success to load data : M73/200508/Task-Log/task.mat:task_data(1).Trial_type
converting cue time series with time rate : 30
export NWB file : ./M73_1.nwb
...
success to load data : M73/200508/Tiff/raw-256/M73-200508-8.tif
success to load data : M73/200508/Task-Log/task.mat:task_data(8).Cue,Rew,X,Y,tY,tX,tYs,tXs,pupil_X,pupil_Y,eyelid_X,eyelid_Y,lick,ccdA,ccdB
success to load data : M73/200508/Task-Log/task.mat:task_data(8).Trial_type
converting cue time series with time rate : 30
export NWB file : ./M73_8.nwb
~~~

  