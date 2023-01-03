# PhD
*Codes written for PhD work in seismology*
*Louisa Brotherson - 03/01/23.*
*This file lists the Matlab2021b codes used to process seismic and mechanical data for my PhD at the University of Liverpool.*
*These are codes that I solely, some of which have been used in published and in-prep papers.*
1. find_stickslip_v4.m: function which uses short-term average, long-term average (STA-LTA) algorithm to find large stress drops, which indicate the start of stick-slip (lab-generated earthquake). 
2. stick-slip_process.m: code to process mechanical Labview data (force, displacement, time, normal stress etc.) to calculate stress drops, coefficient of friction and slip distance. Reads from and writes to files, makes data visualisations of choice. Calculates stastical measures for data analytics
3. read_atf_acoustic.m: function to read in multiple Excel files from an inputted filepath into cell called "celldata" and return filepath, name of file andextension separately. Used to mine seismic data .atf files, which are large, complex data structures.
