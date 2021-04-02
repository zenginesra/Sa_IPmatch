DESCRIPTION  
This folder contains MATLAB code files for simultaneous matching of near-fault records to 
a target response spectrum and to a target Instantaneous Power (IP) spectrum, as described in the manuscript:  

Esra Zengin, Norman Abrahamson
<esrazengin@gmail.com>
Zengin E., and Abrahamson, N.A. (2021) A Procedure for Matching the Near-Fault 
Ground Motions based on Spectral Accelerations and Instantaneous Power, Earthquake Spectra, in press. 

Main Steps:  

1) Change Matlab’s current directory to the example folder of interest and make sure that 
H1 and H2 components of the acceleration files and target response spectrum file are present 
in the folder. Target IPs are computed in the program by the IP_GMM.m function, so the user does not 
need to add the IP spectrum file. 

2) Open the main program runSa_IP_Match.m and input the required information:

accH1 = H1 component of the acceleration (in g) 
accH2 = H2 component of the acceleration (in g) 
dt_original = original time step of the record 
M = moment magnitude of the target scenario earthquake 
Rrup = rupture distance of the target scenario earthquake (km) 
targetSpec = target spectral accelerations (in g) 
Tall = target periods (s) 
T_lower = lower period for matching 
T_upper = upper period for matching 
targetPGA = target peak ground acceleration (in g) 
ksi = damping ratio 
matched_component = component used in the matching process

3) Run the main program. The program will create an output file that contains an adjusted acceleration 
time series. It can also produce two figures that include time histories and Sa- and IP- spectra of 
the record.  

*The example files and the scripts in the folder provide a demonstration of usage of the program. 
