linz-pyaxis package
===================

Provides the pyaxis script for adjusting antenna local tie surveys.

This requires on the linz-adjustment and linz-geodetic packages.  

The adds functionality to the LINZ.Geodetic.Adjustment python module
for calculating the Invariant Point of a radio telescope, which is the 
intersection of the primary axis with the common perpendicular to the 
secondary axis.  

This implementation assumes that the primary axis is vertical.  It would
be a minor modification to adapt this to telescopes with a horizontal 
primary axis.  

