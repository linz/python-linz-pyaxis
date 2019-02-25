linz-pyaxis package
===================

Provides the pyaxis script for adjusting antenna local tie surveys.

This requires the `python-linz-adjustment <http://github.com/linz/python-linz-adjustment>`_
and `python-linz-geodetic <http://github.com/linz/pyhon-linz-geodetic>`_ packages.  

The adds functionality to the LINZ.Geodetic.Adjustment python module
for calculating the Invariant Point of a radio telescope, which is the 
intersection of the primary axis with the common perpendicular to the 
secondary axis.  

This implementation assumes that the primary telescope axis is vertical.  
It would be a minor modification to adapt this to telescopes with a horizontal 
primary axis.  

This software is an extension of the python-linz-adjustment module which adds the 
functionality to calculate the antennae reference points.  It is otherw

Survey requirements
===================

This adjustment method requires a 3 dimensional survey of one or more radio telescopes.  
The survey is assumed to include both fixed control marks on the ground (including GNSS 
survey marks) and marks attached to the radio telescope antennae.

The adjustment can use any observation types supported by the python-linz-adjustment software.

During the survey the antennae marks are observed for different positions of the antenna. 
The antennae positions are in a series of arcs, each of which involves rotation the antenna 
about either the primary or secondary axis, but not both axes.  The antennae is stopped in a 
set of positions around the arc (eg 30 degree steps), and in each position the targets are 
observed.  This is repeated for multiple arcs.

This software assumes that the observations to the targets are labelled by assigning each target
a different name for each position of each arc.  For example the LINZ surveys give each target 
position a name ``rrat`` where ``rr`` is a code for the rotation angle (eg 03, 06), ``a`` is a code
defining the arc (eg X, Y, Z), and ``t`` is a code for the specific antenna target.  The software
can include observations to more than one antennae - each arc code is associated with a specific 
antenna.

The software uses a `configuration file <https://raw.githubusercontent.com/linz/python-linz-pyaxis/master/LINZ/Geodetic/AxisAdjustment.adj>`_ (combined with the basic `adjustment configuration 
file <https://raw.githubusercontent.com/linz/python-linz-adjustment/master/LINZ/Geodetic/Adjustment.adj>`_) 
to how the adjustment is run.  In particular this defines the encoding of the target names.  The 
station naming must be done in such a way that the antenna target observations can be identified
using a `python regular expression <https://docs.python.org/2/library/re.html>`_ (defined by the axis_target_re parameter).
Any station names not matching this regular expression are assumed to belong to the survey control points.

Adjustment methodology
======================

The adjustment runs in three phases: 

* phase 1: unconstrained least squares adjustment in which each target and control station XYZ position is calculated. This ignores the fact that the target positions on each must lie on arcs - each target position is calculated independently.  This requires that the observations are sufficient to locate every target position in 3 dimensions.

* phase 2: antenna initial parameter estimation.  This phase uses the calculated target positions to calculate for each target its position on the antenna, and the primary and secondary axis rotations for each position on the arc.  This is not a rigorous least squares adjustment of the observations - it is just using uncorrelated calculated positions of the targets.

* phase 3: in this phase the least squares adjustment is reparameterised in terms of the target locations on the antennae and the antenna primary and secondary axis rotations, the antenna parameters (location of the ARP, direction of the primary rotation axis, offset and orthogonality), and the control station XYZ positions.  This allows the location of the ARP to be rigorously calculated within the control network.

The pyaxis adjustment can also set or calculate a prism constant for each target in the adjustment.  The prism constant is a correction added to distances measured to the target.  These are defined by target_calibration adjustment configuration parameters for each target.  Where these are not defined for a target the constant is assumed to be zero (ie already applied to the input observations)

The adjustment sinex_output_file option can be used to generate a SINEX file defining the relative positions of the ARP to specific marks such as GNSS reference station antennae.

Using the pyaxis program
========================

To use the pyaxis program the first step is to compile the observation data files.  The formats are described in the `documentation for python-linz-adjustment <https://github.com/linz/python-linz-adjustment>`_.  The antenna target station
names must follow a convention as described above.

If the data do not define the positions of all marks (for example if it doesn't include a SINEX file of XYZ coordinates for some stations), then a station coordinate file will be needed to define at least some of the station coordinates.  This also follows the format defined by the python-linz-adjustment software.  Not all the initial station positions need to be defined - the software will attempt to calculate positions from the observation data.

Finally the adjustment needs a configuration file.  A sample configuration file can be created using the command:
::
    pyaxis -c myconfig.cfg
    
This will create a configuration file naemd myconfig.cfg. This needs to be edited to define the parameters for the survey.  
In particular the parameters most likely to need editing are:

* Input files:
  ::
   coordinate_file
   data_file 
   
* Adjustment configuration:
  ::
   axis_target_re
   antenna_arcs
   target_calibration

* Output files:
  ::
   listing_file
   residual_csv_file
   output_coordinate_file
   sinex_output_file
   sinex_site_id
   sinex_header_info
   






 
