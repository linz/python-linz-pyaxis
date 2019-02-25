pyaxis example
==============

This directory contains example data for the pyaxis program.  This data comes from
the 2015 local tie survey of 50423 Warkworth.  This data is provided for demonstration
purposes only - please refer to Land Information New Zealand for authoritative data
on this station.

To run this example use the command:

   pyaxis wark2015lt.adj

The output files will be created in the output subdirectory (these are already 
provided here for reference).  

This survey determines the relative locations of the antenna reference points of 
two antennae, a 12m dish 50243S001 WARKWORTH and a 30m dish 50243S002 WARKWORTH. 
It also includes the ITRF GNSS station 50243M001 WARKWORTH.

The observations to the antennae are defined in the files antenna12.csv and antenna30.csv.
Within these files the antenna targets have codes such as 30W4, which is the obsevation to
target 4 on arc W with the antenna rotated to position labelled 30 (300 degrees).

This naming scheme matches the axis_target_re in the configuration file.
