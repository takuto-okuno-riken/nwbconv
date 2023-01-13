[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
# nwbconv
 NWB format convert tool.
 This tool could generate a NWB format file from an excel template file and raw image files.
 
## Requirements: software
* MATLAB R2019a or later

## Installation
1. Run the MATLAB software, and setup the [MatNWB](https://github.com/NeurodataWithoutBorders/matnwb) toolbox. 
~~~
!git clone https://github.com/NeurodataWithoutBorders/matnwb.git
cd matnwb
addpath(genpath(pwd));
generateCore();
~~~
2. Download [nwbconv](https://github.com/takuto-okuno-riken/nwbconv/archive/refs/heads/main.zip) zip file.
3. Extract zip file under your working directory <work_path>.
4. "Add Path" extracted directories (i.e. <work_path>/nwbconv-main).
5. Move to <work_path>/nwbconv-main directory and fill out excel template file with your imaging data.

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
>> nwbconv ori099_ch.xlsx
success to load data : ori099_ch.mat:stack
success to load data : logfile_M99_ori099.mat:TS.stim_log
success to load data : logfile_M99_ori099.mat:TS.start_time
success to load data : logfile_M99_ori099.mat:TS.stop_time
success to load data : logfile_M99_ori099.mat:TS.ori_deg
export NWB file : ./ori099_ch.nwb
~~~

~~~
>> nwbconv --session 1:8 M99.xlsx
success to load data : M99/10099/Tiff/raw-256/M99-1.tif
success to load data : M99/10099/Task-Log/task.mat:task_data(1).Cue,Rew,X,Y,tY,tX,tYs,tXs,pupil_X,pupil_Y,eyelid_X,eyelid_Y,lick,ccdA,ccdB
success to load data : M99/10099/Task-Log/task.mat:task_data(1).Trial_type
converting cue time series with time rate : 30
export NWB file : ./M99_1.nwb
...
success to load data : M99/10099/Tiff/raw-256/M99-8.tif
success to load data : M99/10099/Task-Log/task.mat:task_data(8).Cue,Rew,X,Y,tY,tX,tYs,tXs,pupil_X,pupil_Y,eyelid_X,eyelid_Y,lick,ccdA,ccdB
success to load data : M99/10099/Task-Log/task.mat:task_data(8).Trial_type
converting cue time series with time rate : 30
export NWB file : ./M99_8.nwb
~~~

  
