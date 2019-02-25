#!/usr/bin/python

# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
import re
import math
import numpy as np
from numpy import linalg
from LINZ.Geodetic.Adjustment import Adjustment, Options, Plugin
from LINZ.Geodetic.Station import Station
from LINZ.Geodetic.Ellipsoid import GRS80
from LINZ.Geodetic import Sinex

def veclen(x):
    return math.sqrt(np.vdot(x,x))

# Define parameters of antenna axes, antenna orientation (rotation about axis),
# an target position along the axis.  These are assigned to targets, arcs, etc, 
# but can be reassigned so that several objects are using the same parameter, as
# target mark locations are used to define axes, which are then merged to arc axes,
# and then again to antenna axes.
#
# Each parameter has three main functions:
#
#  initParams, which sets up the parameters it uses
#  paramCount, counts the parameters used
#  apply, which applies the parameters to a coordinate and optionally observation equations
#    (Note that TargetParams has getXYZ instead of apply)
#  update, which updates the parameters from the least squares result parameter vector


class AngleParam( object ):
    '''
    Defines a rotation of an antenna around the Z axis.  The angle is the 
    rotation applied relative to an arbitrary 0 position to get to the
    rotated position.  

    angle is the rotation in radians
    rmat is the rotation matrix, such that 

    uvw' = uvw.dot(rmat.T)

    where uvw is the unrotation position in terms of an axis based coord
    system, uvw' is the rotated position.  uvw and uvw' are (3,) arrays,
    rmat is a (3,3) array
    '''

    def __init__( self, name, angleName=None ):
        self.name=name
        self.angleName=angleName or name
        self.angleParam=None
        self.setAngle( 0.0 )

    def setAngle( self, angle ):
        self.angle=angle
        ca=math.cos(angle)
        sa=math.sin(angle)
        self.rmat=np.array(((ca,-sa,0.0),(sa,ca,0.0),(0.0,0.0,1.0)))

    def initParams( self, adjustment, calculate=False ):
        self.angleParam=None
        if calculate:
            self.angleParam=adjustment.addParameter( self.angleName )
            adjustment.addParameterUpdate( self.updateAngle )

    def updateAngle( self, params, paramcvr=None ):
        if self.angleParam is not None:
            self.setAngle(self.angle+params[self.angleParam])
        return True

    def apply( self, uvw, obseq=None ):
        '''
        Updates a position and optionally observation equations to account
        for the rotation.  The position is a (3,) array, and the 
        observation equations a (3,nprm) array.

        The parameters for the axis are the rotation about the u and v axes.
        (w is the rotation axis of the antenna axis)

        '''
        if obseq is not None:
            if self.angleParam is not None:
                obseq[0,self.angleParam]=-uvw[1]
                obseq[1,self.angleParam]=uvw[0]
            t1=self.rmat.dot(obseq)
            obseq[:]=t1

        uvw2=uvw.dot(self.rmat.T)
        return uvw2

class AxisParams( object ):

    paramNames=['u rotation','v rotation','centre X','centre Y','centre Z']

    def __init__( self, name, axes, centre ):
        '''
        axes is a mapping from axes u,v,w to x,y,z where w is the axis of rotation of 
        the antenna axis and u,v orthogonal axes.  Centre is the xyz position of the reference
        point of the antenna.

        uvw=(xyz-centre).dot(axes.T)
        xyx=uvw.dot(axes)+centre
        '''
        self.name=name
        self.axes=axes.copy()
        self.centre=centre.copy()
        self.params=[None]*5

    def initParams( self, adjustment, calcAlignment=False, calcCentre=False ):
        '''
        Initiallize adjustment parameters for the antenna axis.  This can include
        alignment (reorientation of axes around u,v axes), rotation (around w axis)
        and offset of the centre

        Parameters are for incremental changes to each of these values
        '''

        self.params=[None]*5
        if calcAlignment:
            self.params[0]=adjustment.addParameter(self.name+' U rotation')
            self.params[1]=adjustment.addParameter(self.name+' V rotation')
        if calcCentre:
            self.params[2]=adjustment.addParameter(self.name+' centre X')
            self.params[3]=adjustment.addParameter(self.name+' centre Y')
            self.params[4]=adjustment.addParameter(self.name+' centre Z')
        if calcAlignment or calcCentre:
            adjustment.addParameterUpdate( self.updateAxis )


    def xyzCoordinates( self ):
        return self.centre,self.params[2:5]

    def apply( self, uvw, obseq=None ):
        '''
        Updates a position and optionally observation equations relative to the 
        axis axes to the global system.  The position is a (3,) array, and the 
        observation equations a (3,nprm) array.

        The parameters for the axis are the rotation about the u and v axes.
        (w is the rotation axis of the antenna axis)
        '''
        if obseq is not None:
            if self.params[0] is not None:
                p0,p1=self.params[:2]
                obseq[1,p0]=-uvw[2]
                obseq[2,p0]=uvw[1]
                obseq[2,p1]=-uvw[0]
                obseq[0,p1]=uvw[2]
            t1=self.axes.T.dot(obseq)
            obseq[:]=t1

            if self.params[2] is not None:
                obseq[0,self.params[2]]=1.0
                obseq[1,self.params[3]]=1.0
                obseq[2,self.params[4]]=1.0
            
        xyz=self.axes.T.dot(uvw)+self.centre

        return xyz

    def updateAxis( self, params, paramcvr=None ):
        if self.params[0] is not None:
            angle=params[self.params[0]]
            ca=math.cos(angle)
            sa=math.sin(angle)
            self.axes=np.array(((1.0,0.0,0.0),(0.0,ca,sa),(0.0,-sa,ca))).dot(self.axes)
            angle=params[self.params[1]]
            ca=math.cos(angle)
            sa=math.sin(angle)
            self.axes=np.array(((ca,0.0,-sa),(0.0,1.0,0.0),(sa,0.0,ca))).dot(self.axes)
        if self.params[2] is not None:
            self.centre += params[[self.params[2:5]]]
        return True

    @staticmethod
    def nearestPoint( axis1, axis2 ):
        '''
        Find the nearest point of approach of two axes.  

        Returns the nearest points w1, w2

        Based on parameterised point c1+t1.w1 on axis1 and c2+t2.w2 on axis2.  The
        difference between them must be orthogonal to w1 and w2.  Solve for t1, t2.
        '''
        
        c1=axis1.centre
        c2=axis2.centre
        w1=axis1.axes[2]
        w2=axis2.axes[2]

        w1w2=w1.dot(w2)

        A=np.array(((1.0,-w1w2),(-w1w2,1.0)))
        b=np.array(((c2-c1).dot(w1),(c1-c2).dot(w2)))

        t1,t2=linalg.solve(A,b)

        p1=c1+t1*w1 
        p2=c2+t2*w2

        return p1,p2


class ElevationAxisParams( AngleParam ):
    '''
    Parameters defining the configuration of the elevation axis relative to the azimuth axis.  
    Converts the uvw vector in terms of the centre and orientation of the elevation axis to 
    uvw in terms of the azimuth axis, where the azimuth axis u vector is perpendicular to the 
    azimuth and elevation rotation axes (direction of cross product of elevation with azimuth
    rotation axis), and is the u vector for both coordinate systems.

    Configuration parameters are the offset of the elevation axis from the azimuth axis in 
    the u direction, and the rotation of the axis w vector relative to the v vector of the azimuth
    coordinate system. 
    '''

    rot1=np.array(((0.0,0.0,1.0),(0.0,-1.0,0.0),(1.0,0.0,0.0)))
    rot2=np.array(((0.0,0.0,1.0),(1.0,0.0,0.0),(0.0,1.0,0.0)))

    def __init__( self, name, offset, angle ):
        AngleParam.__init__( self, name, angleName=name+' misalignment' )
        self.setAngle(angle)
        self.offset=offset
        self.offsetParam=None

    def initParams( self, adjustment, calcOffset=False, calcAngle=False ):
        self.offsetParam=None
        if calcOffset:
            self.offsetParam=adjustment.addParameter(self.name+' offset')
            adjustment.addParameterUpdate( self.updateOffset )
        AngleParam.initParams( self, adjustment, calcAngle )

    def updateOffset( self, params, paramcvr=None ):
        '''
        Update values based on param vector from adjustment
        '''
        if self.offsetParam is not None:
            self.offset += params[self.offsetParam]
        return True

    def apply( self, uvw, obseq=None):
        '''
        Convert the position in terms of the elevation axis to a position in terms of 
        the azimuth axis.
        Optionally add the observation equations defining the position in terms of 
        incremental changes to offset, radius, and angle.

        Observation equations are (3,nprm) array.
        '''

        # First convert elevation u,v,w to w,-v,u so that misalignment is
        # rotation around w axis (which is how AngleParam is defined)

        uvw=self.rot1.dot(uvw)
        # Save obseq so that can replace contents with final version
        obseq0=obseq
        if obseq is not None:
            obseq=self.rot1.dot(obseq)

        # Apply the misalignment rotation
        uvw=AngleParam.apply( self,uvw,obseq)

        # Convert to azimuth u,v,w ...
        #  u,v,w -> v,w,u

        uvw=self.rot2.dot(uvw)
        if obseq is not None:
            obseq=self.rot2.dot(obseq)

        # Apply the axis offset (in the u direction)

        uvw[0] += self.offset

        if obseq is not None:
            obseq0[:,:]=obseq
            if self.offsetParam is not None:
                obseq0[0,self.offsetParam]=1.0

        return uvw

class ElevationArcParams( object ):
    '''
    Class combining the parameters relating to an elevation arc.

    Does not provide initiallization or update - assumes this is done for
    the individual components.  Used to provide an arc parameter set for 
    TargetMark on elevation arcs.
    '''

    def __init__( self, name, azimuthAxisParams, elevationAxisParams, arcOrientationParams ):
        self.azimuthAxisParams=azimuthAxisParams
        self.elevationAxisParams=elevationAxisParams
        self.arcOrientationParams=arcOrientationParams

    def apply( self, uvw, obseq=None ):
        uvw=self.elevationAxisParams.apply( uvw, obseq )
        uvw=self.arcOrientationParams.apply( uvw, obseq )
        uvw=self.azimuthAxisParams.apply( uvw, obseq )
        return uvw

class TargetParams( AngleParam ):
    '''
    Defines the location of a target on an antenna relative to the reference point of the antenna and
    the orientation of its axes (u,v,w, where w is along axis of rotation).  Defined in polar coordinates.

    offset:  distance along the rotation axis from the centre
    radius:  distance of target from rotation axis
    angle:   angle of line from rotation axis to target relative to axis.
    '''

    def __init__( self, name ):
        AngleParam.__init__(self,name,angleName='Target '+name+' angle')
        self.offset=0.0
        self.radius=0.0
        self.tgtParams=[None,None]

    def initParams( self, adjustment, calcOffset=False, calcRadius=False, calcAngle=False ):
        self.tgtParams=[None,None]
        if calcOffset:
            self.tgtParams[0]=adjustment.addParameter('Target '+self.name+' offset')
        if calcRadius:
            self.tgtParams[1]=adjustment.addParameter('Target '+self.name+' radius')
        if calcOffset or calcRadius:
            adjustment.addParameterUpdate( self.update )
        AngleParam.initParams( self, adjustment, calcAngle )

    def update( self, params, paramcvr=None ):
        '''
        Update values based on param vector from adjustment
        '''
        if self.tgtParams[0] is not None:
            self.offset += params[self.tgtParams[0]]
        if self.tgtParams[1] is not None:
            self.radius += params[self.tgtParams[1]]
        return True

    def getXYZ( self, obseq=None):
        '''
        Calculate the position of the mark in the unrotated uvw system of the antenna axis,
        and optionally add the observation equations defining the position in terms of 
        incremental changes to offset, radius, and angle.

        Observation equations are (3,nprm) array.
        '''
        uvw=np.array((self.radius,0.0,self.offset))

        if obseq is not None:
            if self.tgtParams[0] is not None:
                obseq[2,self.tgtParams[0]]=1.0
            if self.tgtParams[1] is not None:
                obseq[0,self.tgtParams[1]]=1.0
        uvw=AngleParam.apply(self,uvw,obseq)
        return uvw

class TargetCalibration( object ):

    def __init__( self, name, correction=0.0, calculate=False ):
        self.name=name
        self.correction=correction
        self.initcorrection=correction
        self.calculate=calculate
        self.prmno=None

    def initParams( self, adjustment ):
        self.prmno=None
        if self.calculate:
            self.prmno=adjustment.addParameter('Target '+self.name+' prism calibration')
            adjustment.addParameterUpdate( self.update )
            
    def update( self, params, paramcvr=None ):
        if self.prmno is not None:
            self.correction += params[self.prmno]
        return True

class TargetMark( object ):
    '''
    Representing an observed location of a specific physical target at a named
    orientation observed as part of an arc of observations as the antenna is 
    rotated around an axis.

    In the survey this location is represented by a mark (observed point)
    '''

    def __init__( self, mark, target, arc, orientation ):
        self.mark=mark
        self.target=target
        self.arc=arc
        self.orientation=orientation

    def calculateXyz( self, obseq=None ):
        xyz=self.target.params.getXYZ( obseq )
        xyz=self.orientation.apply( xyz, obseq )
        xyz=self.arc.axisParams.apply( xyz, obseq )
        return xyz

class TargetMarkAdjustment( object ):

    def __init__( self, adjustment, marks ):
        self.adjustment=adjustment
        self.write=adjustment.write
        self.options=adjustment.options
        self.marks=marks
        self.updates=[]
        self.parameters=[]

    def addParameter( self, paramName ):
        self.parameters.append(paramName)
        return len(self.parameters)-1

    def addParameterUpdate( self, func ):
        self.updates.append( func )

    def solve( self ):
        '''
        Solve equations for a list of marks.  Assumes that the target 
        parameters have already been initiallized.
        '''

        nprm=len(self.parameters)
        marks=self.marks
        nobs=len(marks)*3
        maxIterations=self.options.axisMaxIterations
        convergence=self.options.axisConvergence
        writeIterations=self.options.writeTargetAdjustmentIterations
        writeTotalChange=self.options.debugTargetAdjustmentTotalChange

        iteration=0
        ssr0=0
        obsval0=np.zeros((nobs,))
        converged=False
        totalChange=np.zeros((nprm,))
        while iteration < maxIterations:
            iteration += 1
            obseq=np.zeros((nobs,nprm))
            obsval=np.zeros((nobs,))
            for imrk,mark in enumerate(marks):
                iobs=imrk*3
                mrkobs=obseq[iobs:iobs+3,:]
                xyz=mark.calculateXyz( mrkobs )
                obsval[iobs:iobs+3]=mark.mark.xyz()-xyz

            residuals=linalg.norm(obsval.reshape((nobs // 3,3)),axis=1)
            imax=residuals.argmax()
            if writeIterations:
                self.write("     Iteration {2}: Maximum coordinate residual {0:.5f}m at {1}\n".
                       format(residuals[imax], marks[imax].mark.code(),iteration))
            ssr=veclen(obsval)
            zero=np.where((obseq == 0.0).all(axis=0))[0]
            if zero.size > 0:
                zerolist=' '.join((str(index) for index in zero))
                self.write("**** Parameters {0} not used in observations\n".format(zerolist))
                break
            paramValues,resid,rank,sval=linalg.lstsq(obseq,obsval)
            if rank < nprm:
                self.write("**** Axis solution singular (rank {0} < {1})\n".format(rank,nprm))
                break
            for update in self.updates:
                update(paramValues)
            totalChange += paramValues

            criteria= np.abs(obsval-obsval0).max() 
            if criteria < convergence:
                converged=True
                break
            obsval0=obsval

        if writeTotalChange:
            self.write("     Total parameter changes in adjustment\n")
            prmNames={}
            for p in params:
                prmNames.update(p.parameters())
            for i, v in enumerate(totalChange):
                self.write("       {0:3d} {1:<30s} {2:15.8f}\n".format(i,prmNames[i],v))

        if not converged:
            self.write("**** Failed to converge on axis solution after {0} iterations\n".format(iteration))
        else:
            self.write("     Target axis solution converged after {0} iterations\n".format(iteration))
        residuals=linalg.norm(obsval0.reshape((nobs // 3,3)),axis=1)
        imax=residuals.argmax()
        self.write("     Maximum coordinate residual {0:.5f}m at {1}\n".
                   format(residuals[imax], marks[imax].mark.code()))

        return converged

class Orientation( AngleParam ):

    def __init__( self, name, code=None ):
        AngleParam.__init__( self, name )
        self.code=code or name

class Target( object ):
    '''
    Represents the observations of a target on one arc.  Used to calculate the
    antenna axis based on this mark, before compiling all marks for the arc.
    '''

    def __init__( self, name, adjustment ):
        self.name=name
        self.write=adjustment.write
        self.adjustment=adjustment
        self.marks=[]
        self.nrow0=0
        self.params=TargetParams(self.name)
        self.axisParams=None

    def addAngle( self, mark, arc, orientationName ):
        # Create a new observation of the target.  Initially each target has
        # its own Orientation parameters, and axis parameters for the target
        targetMark=TargetMark(mark,self,self,Orientation('Orientation '+orientationName+' for '+self.name,orientationName)) 
        self.marks.append( targetMark )
        return targetMark

    def setup( self ):
        self.write(' Setting up target {0}\n'.format(self.name))
        self.write('   Target has {0} observations\n'.format(len(self.marks)))
        if len(self.marks) < 3:
            raise RuntimeError('Not enough observations to calculate target '+self.name)

        # Sort the targets by the angle code

        self.marks.sort(key=lambda x: x.orientation.name)

        # Initial calculation of axes.  Start with the most 
        # distant pair of marks...

        maxdist=0.0
        mark1=None
        mark2=None
        for i,a1 in enumerate(self.marks[:-1]):
            for a2 in self.marks[i+1:]:
                distance=a1.mark.distanceTo(a2.mark)
                if distance > maxdist:
                    mark1,mark2,maxdist=a1,a2,distance

        # Find a third mark with a maximum offset with the line between
        # the first two

        uvec=mark1.mark.vectorTo(mark2.mark)
        uvec /= veclen(uvec)
        mark3=None
        maxdist=0
        for a3 in self.marks:
            if a3 == mark1 or a3 == mark2:
                continue
            offset=mark1.mark.vectorTo( a3.mark )
            offset-=np.vdot(uvec,offset)*uvec
            distance=veclen(offset)
            if distance > maxdist:
                mark3,maxdist=a3,distance

        # Now use these three marks to estimate an arc axis, and centre
        # uvec is unit vector from a1 to a2
        # vvec is perp unit vector in plane of a3 
        # wvec is unit vector in direction of proposed axis
        # axis_axes transform from XYZ to axis coordinates

        wvec=np.cross(uvec,mark1.mark.vectorTo(mark3.mark))
        wvec /= veclen(wvec)
        vvec = np.cross(wvec,uvec)
        axis_axes=np.vstack((uvec,vvec,wvec))

        # Convert midpoints of vectors a1-a2 and a1-a3 to axis coordinate system

        v2=mark1.mark.vectorTo(mark2.mark).dot(axis_axes.T)/2
        v3=mark1.mark.vectorTo(mark3.mark).dot(axis_axes.T)/2

        # Now need to find intersection of bisectors of each vector
        # Bisector is defined by
        #    x = v2[0] - v2[1] * t
        #    y = v2[1] + v2[0] * t
        # Form equations to solve for t2, t3 at which x value are the same

        A=np.array(((v2[1],-v3[1]),(-v2[0],v3[0])))
        b=np.array(((v2[0]-v3[0],),(v2[1]-v3[1],)))
        t2=linalg.solve(A,b)[0,0]
        vc=v2+t2*np.array([-v2[1],v2[0],0])
        centre=mark1.mark.xyz()+axis_axes.T.dot(vc)

        # Rotate the axes to set the u vector to point from the centre to the
        # first mark.

        uvec = self.marks[0].mark.xyz()-centre
        uvec -= (uvec.dot(wvec))*uvec
        uvec /= veclen(uvec)
        vvec = np.cross(wvec,uvec)
        axis_axes=np.vstack((uvec,vvec,wvec))

        # Assign an angle to each target, and determine a mean radius

        radii=[]
        offsets=[]
        for a in self.marks:
            axyz=a.mark.xyz()-centre
            ua=uvec.dot(axyz)
            va=vvec.dot(axyz)
            radii.append(math.hypot(ua,va))
            a.orientation.setAngle(math.atan2(va,ua))
            offsets.append(wvec.dot(axyz))

        self.params.setAngle(0.0)
        self.params.radius=np.mean(radii)

        # Now refine to get best fit axes..
        #
        # Iterative least squares with parameters corresponding to
        #   X,Y,Z of centre
        #   rotation around u and v axes to change axis direction
        #   radius of target arc
        #   angle of each target on the arc

        axisParams=AxisParams( 'Target '+self.name+' arc axes', axis_axes, centre )
        self.axisParams=axisParams

        tmadjustment=TargetMarkAdjustment( self.adjustment, self.marks )
        self.params.initParams( tmadjustment, calcRadius=True )
        axisParams.initParams( tmadjustment, calcAlignment=True, calcCentre=True )
        for mark in self.marks:
            mark.axis=self.axisParams
            mark.orientation.initParams( tmadjustment, calculate=True )
        tmadjustment.solve()

    def alignAxes( self, axisParams, ownOrientations=False ):
        '''
        Reorient current target axes to best align with specified set.
        May require flipping axes, as well as corresponding change to
        angles for targets and observed angles.
        
        Updates the target axis parameters with the new parameters

        Optionally can also update the target angles - this should only
        be done when they belong to this target only...
        '''
        # Check if axes need to be flipped..

        axes=axisParams.axes
        centre=axisParams.centre
        uvw = self.axisParams.axes

        axisAlignment=math.degrees(math.acos(min(1.0,abs(axes[2].dot(uvw[2])))))

        if uvw[2].dot(axes[2]) < 0:
            self.params.setAngle(-self.params.angle)
            self.params.offset=-self.params.offset
            if ownOrientations:
                for mark in self.marks:
                    mark.orientation.setAngle(-mark.orientation.angle)

        offset=(self.axisParams.centre-centre).dot(axes[2])
        self.params.offset += offset

        # Determine alignment of current u axis in new system and
        # subtract from target angle

        cs=uvw[0].dot(axes[0])
        sn=uvw[0].dot(axes[1])
        angle=math.atan2(sn,cs)
        self.params.setAngle(self.params.angle + angle)

        self.axisParams=axisParams
        
    def countCommonOrientations( self, orientations ):
        '''
        Count how many of keys of angles dictionary are used for this
        target.
        '''
        return sum((1 for mark in self.marks if mark.orientation.code in orientations))

    def mergeOrientations( self, orientations, name ):
        '''
        Adjust observed orientations against target orientations to best match supplied
        list of observed orientations, and add any new orientations to the list.
        Replace targets orientations with those in the list to provide a common
        set of antenna orientation angles for all the targets.

        Angles is a dicitionary keyed on the orientation angle name.
        '''

        diff0=0.0
        if len(orientations) > 0:
            diffs=[]
            for mark in self.marks:
                code=mark.orientation.code
                if code in orientations:
                    diffs.append(mark.orientation.angle-orientations[code].angle)

            if not diffs:
                raise RuntimeError('Cannot merge orientations for target '+self.name+' - no common orientations')

            diff0=diffs[0]
            diffs=np.remainder(np.array(diffs)-diff0+math.pi,2*math.pi)+diff0-math.pi
            diff0=diffs.mean()
        
        self.params.setAngle( self.params.angle + diff0 )
        for mark in self.marks:
            code = mark.orientation.code
            if code not in orientations:
                orientation=Orientation('Orientation '+code+' for '+name,code)
                orientation.setAngle(mark.orientation.angle-diff0)
                orientations[code]=orientation
            mark.orientation=orientations[code]

        return orientations

class Arc( object ):

    def __init__( self, name, adjustment ):
        self.name=name
        self.adjustment=adjustment
        self.write=adjustment.write
        self.targets={}
        self.orientations=None
        self.marks=[]
        self.axisParams=None
        self.isAzimuth=False
        # azimuthAngle used for elevation arcs to define orientation of 
        # azimuth axis..
        self.azimuthAngle=None

    def getTarget( self, targetName ):
        if targetName not in self.targets:
            self.targets[targetName]=Target(targetName,self.adjustment)
        return self.targets[targetName]

    def addTargetMark( self, mark, targetName, orientationName ):
        target=self.getTarget(targetName)
        targetMark=target.addAngle(mark,self,orientationName)
        self.marks.append(targetMark)
        return targetMark

    def setup( self ):
        # Set up the individual target arcs
        for targetName in sorted(self.targets):
            target=self.targets[targetName]
            target.setup()

        # Combine target arcs to a single axis
        # First determine a good set of initial orthogonal axes for the arc.
        
        # For the axis w vector (along rotation axis), first combine the 
        # w vectors for each target

        self.write(" Combining targets axes to arc axis\n")
        wvec=None
        targets=self.targets.values()
        for t in targets:
            tw=t.axisParams.axes[2]*t.params.radius
            if wvec is None:
                wvec=tw
            else:
                if tw.dot(wvec) < 0:
                    tw=-tw
                wvec+=tw

        # Next find the nearest and farthest centre in the direction 
        # of the w vector, and add the vector between them
        #
        # (This has been disabled as want to look at offsets of centres
        # relative to preferred axis)

        targets.sort( key=lambda t: wvec.dot(t.axisParams.centre) )
        centre=targets[0].axisParams.centre
        if False:
            wvec += (targets[-1].axisParams.centre - targets[0].axisParams.centre)
            wvec /= veclen(wvec)

        # Check if this is an azimuth axis, and if so ensure that the vector is pointing up

        lon,lat,hgt=GRS80.geodetic(centre)
        local_enu=GRS80.enu_axes(lon,lat)
        vdot=wvec.dot(local_enu[2])
        if vdot < 0.0:
            wvec=-wvec
            vdot=-vdot
        self.isAzimuth=vdot > 0.7 # Approx cos(45)!

        # Form an orthonormal set of axis using the plane of the target with
        # the greatest radius
        wvec /= veclen(wvec)
        targets.sort( key=lambda t: t.params.radius )
        uvec=np.cross(targets[-1].axisParams.axes[1],wvec)
        uvec /= veclen(uvec)
        vvec = np.cross(wvec,uvec)
        uvw=np.vstack((uvec,vvec,wvec ))

        # Find a good mean centre, and check offsets of targets
        offsets=[]
        targets.sort(key=lambda t: t.name)
        for t in targets:
            offsets.append(uvw.dot(t.axisParams.centre-centre))
        offsets=np.array(offsets)
        offset=offsets.mean(axis=0)
        centre += offset[0]*uvw[0]
        centre += offset[1]*uvw[1]
        offsets -= offset
        self.write("     Centre offsets (X,Y,distance)\n")
        for t,offset in zip(targets,offsets):
            ds=math.hypot(offset[0],offset[1])
            self.write("       {0:4s} {1:7.4f} {2:7.4f} {3:7.4f}\n"
                       .format(t.name,offset[0],offset[1],ds))

        # Now set up axis parameters for the arc
        axisParams=AxisParams( 'Arc '+self.name,uvw, centre )

        # Finally align each target with the new axes

        targets.sort( key=lambda t: len(t.marks), reverse=True)
        orientations={}

        for i in range(len(targets)):
            if i > 0:
                targets[i:].sort( 
                    key=lambda t: t.countCommonOrientations(orientations),
                    reverse=True)
            targets[i].alignAxes(axisParams,ownOrientations=True)
            orientations=targets[i].mergeOrientations(orientations,'Arc '+self.name)

        self.orientations=orientations.values()
        self.axisParams=axisParams

        # Now look at refining arc axes by least squares adjustment...
        #
        # Parameters are:
        #   centre and orientation of axes
        #   offset, angle, and radius of each target, apart from 
        #      offset and angle of the first target
        #   each orientation angle apart from the first

        # Set up the marks to reference this arc for axis parameters
        for m in self.marks:
            m.arc=self

        tmadjustment=TargetMarkAdjustment( self.adjustment, self.marks )
        axisParams.initParams( tmadjustment, calcAlignment=True, calcCentre=True )
        targets[0].params.initParams( tmadjustment, calcRadius=True )
        for t in targets[1:]:
            t.params.initParams( tmadjustment, calcRadius=True, calcOffset=True, calcAngle=True )
        for o in self.orientations:
            o.initParams( tmadjustment, calculate=True )

        tmadjustment.solve()

        # Write summary information

        names=' '.join(sorted(self.targets.keys()))
        orientationCodes=' '.join(sorted((o.code for o in self.orientations)))

        self.write(" Arc is for {0} axis\n".format('azimuth' if self.isAzimuth else 'elevation'))
        self.write(" Arc uses targets: {0}\n".format(names))
        self.write(" Arc uses orientations: {0}\n".format(orientationCodes))

    def maxMarkOffset( self ):
        '''
        Check maximum mark offset with current parameters..
        '''
        maxmark=None
        maxdiff=0.0
        for mark in self.marks:
            xyzo=mark.mark.xyz()
            xyzc=mark.calculateXyz()
            diff=veclen(xyzo-xyzc)
            if maxmark is None or maxdiff < diff:
                maxdiff=diff
                maxmark=mark
        return maxdiff,maxmark.mark.code()

    def alignAxes( self, axisParams ):
        '''
        Realign arc to new axes
        '''
        for t in self.targets.values():
            t.alignAxes( axisParams )
            
        if axisParams.axes[2].dot(self.axisParams.axes[2]) < 0:
            for o in self.orientations:
                o.setAngle( -o.angle )

        self.axisParams=axisParams

    def countCommonTargets( self, targetParams ):
        return sum((1 for t in self.targets.values() if t.name in targetParams))


    def mergeTargetParams( self, targetParams, apply=True, reverse=False ):
        '''
        Combine targets for an arc to ensure a common set of parameters.
        Adjusts centre of arc axisParams to match offsets, adjusts targets
        to match angles.  Replaces target parameters with common parameters
        from the targets dictionary if apply is True.  If reverse=True then
        the value is tested with the offsets and angles reversed, 

        Returns an updated list of targets adding those from this arc not 
        already included, and measure of the greatest difference between
        this arc and the supplied list of targets.
        '''

        assert not (apply and reverse)
        if reverse:
            apply=False

        diff0=0
        offset0=0
        maxError=0.0
        if len(targetParams) > 0:
            factor=-1.0 if reverse else 1.0
            diffs=[]
            offsets=[]
            radii=[]
            for t in self.targets.values():
                if t.name in targetParams:
                    params=targetParams[t.name]
                    diffs.append(params.angle - t.params.angle*factor)
                    offsets.append(params.offset - t.params.offset*factor )
                    radii.append((params.radius,t.params.radius))

            if apply and len(diffs) == 0:
                raise RuntimeError('Cannot combine targets for arc '+self.name+
                                   ': no common targets')

            diff0=diffs[0]
            diffs=np.remainder(np.array(diffs)-diff0+math.pi,2*math.pi)+diff0-math.pi
            diff0=diffs.mean()
            diffs -= diff0

            offsets=np.array(offsets)
            offset0=offsets.mean()
            offsets -= offset0

            for a,o,r in zip(diffs,offsets,radii):
                r0,r1=r
                error=veclen([r0-r1*math.cos(a),r1*math.sin(a),o])
                if error > maxError:
                    maxError=error
                
        if apply:
            for t in self.targets.values():
                if t.name not in targetParams:
                    params=TargetParams(t.name)
                    params.radius=t.params.radius
                    params.offset=t.params.offset+offset0
                    params.setAngle(t.params.angle+diff0)
                    targetParams[t.name]=params
                t.params=targetParams[t.name]
            for o in self.orientations:
                o.setAngle(o.angle-diff0)

        return targetParams,maxError

class Antenna( object ):

    def __init__( self, name, adjustment ):
        self.name=name
        self.adjustment=adjustment
        self.write=adjustment.write
        self.arcs={}
        self.targetMarks=[]
        self.elevationArcs=[]
        self.azimuthArcs=[]
        self.azimuthAxisParams=None
        self.elevationTargetParams=[]
        self.azimuthTargetParams=[]
        self.paramList=[]

    def arc( self, arcName ):
        if arcName not in self.arcs:
            self.arcs[arcName]=Arc(arcName,self.adjustment)
        return self.arcs[arcName]

    def addTargetMark( self, mark, targetName, arcName, orientationName ):
        arc=self.arc(arcName)
        # Intially each mark gets its own orientation parameter and each target
        # its own arc parameter.  These get combined during the axis setup
        # phase.
        targetMark=arc.addTargetMark( mark, targetName, orientationName )

        self.targetMarks.append(targetMark)
        return targetMark

    def setup( self ):
        self.write("\nSetting up axes for antenna {0}\n".format(self.name))
        # Derive the rotation axes and targets for each arc
        for arcName in sorted(self.arcs.keys()):
            self.write("\nInitiallizing arc {0}\n".format(arcName))
            arc=self.arcs[arcName]
            arc.setup()

        # Separate out azimuth and elevation arcs
        azimuthArcs=[a for a in self.arcs.values() if a.isAzimuth]
        elevationArcs=[a for a in self.arcs.values() if not a.isAzimuth]
        if len(azimuthArcs) == 0:
            raise RuntimeError('No azimuth arcs observed for antenna '+self.name)
        if len(elevationArcs) == 0:
            raise RuntimeError('No elevation arcs observed for antenna '+self.name)

        self.elevationArcs=elevationArcs
        self.azimuthArcs=azimuthArcs

        # Align the elevation arcs to ensure that the axis is running in a consistent 
        # direction wrt the antenna.  Sorting is to have max common targets for each
        # test

        self.write("\n Merging elevation axis targets\n");
        elevationArcs.sort(key=lambda a: len(a.targets),reverse=True)
        targetParams={}
        for i in range(len(elevationArcs)):
            if i > 0:
                elevationArcs[i:].sort(
                    key=lambda a: a.countCommonTargets(targetParams),
                    reverse=True)
            arc=elevationArcs[i]

            if i > 0:
                targets,maxErrorF=arc.mergeTargetParams( targetParams, apply=False )
                targets,maxErrorR=arc.mergeTargetParams( targetParams, apply=False, reverse=True )
                if maxErrorR < maxErrorF:
                    self.write("     Reversing axes for arc {0}\n".format(arc.name))
                    axes=arc.axisParams.axes.copy()
                    axes[1]=-axes[1]
                    axes[2]=-axes[2]
                    arcAxisParams=AxisParams(arc.name, axes, arc.axisParams.centre )
                    arc.alignAxes(arcAxisParams)

            targets, maxError=arc.mergeTargetParams( targetParams )
            self.write("     Arc {0}: maximum offset merging targets {1:.4f}\n"
                       .format(arc.name,maxError))


        self.elevationTargetParams=targetParams.values()
        self.elevationTargetParams.sort(key=lambda t: t.name)

        # Combine azimuth arc axis and target parameters
        # Calculate a centre that is best aligned with elevation axes
        # (ie trial Antenna Reference Point)


        wvec=azimuthArcs[0].axisParams.axes[2]
        centre=azimuthArcs[0].axisParams.centre
        for arc in azimuthArcs[1:]:
            wvec += arc.axisParams.axes[2]
            centre+=arc.axisParams.centre
        wvec /= len(azimuthArcs)
        centre /= len(azimuthArcs)

        lon,lat,hgt=GRS80.geodetic(centre)
        local_enu=GRS80.enu_axes(lon,lat)
        uvec=np.cross(local_enu[1],wvec)
        uvec /= veclen(uvec)
        vvec=np.cross(wvec,uvec)
        axes=np.vstack((uvec,vvec,wvec))

        # Determine the ARP and reset the centre of the
        # elevation arcs to match it. Also find the orientation of the 
        # elevation axis for each arc.

        # offset is offset of elevation axis from azimuth axis
        # misalignment is offset of elevation arc from perpendicular to
        # normal axis
        # elChangeOffsets offsets to elevation axis target position to move centre
        # elAngleOffsets offsets to elevation axis target to rotate onto new axes

        axisParams=AxisParams(self.name,axes,centre)
        centre=[0,0,0]
        offsets=[]
        misalignments=[]
        elOffsets=[]
        for arc in elevationArcs:
            p1,p2=AxisParams.nearestPoint(axisParams,arc.axisParams)
            centre += p1
            arcAxes=arc.axisParams.axes

            # arcu, arcv, arcw are the arc axes rotated to align arcu with
            # with the normal to the two axes (ie in the direction of the offset 
            # between them)
            arcw=arcAxes[2]
            arcu=np.cross(arcw,wvec)
            arcu /= veclen(arcu)
            arcv = np.cross(arcw,arcu)
            offset=arcu.dot(p2-p1)
            orientation=math.atan2(vvec.dot(arcu),uvec.dot(arcu))
            misalignment=math.atan2(wvec.dot(arcw),-wvec.dot(arcv))
            offsets.append(offset)
            misalignments.append(misalignment)

            # Calculate the offset of the elevation axis normal point from
            # the current centre of the elevation axis parameters
            elOffset=arcAxes[2].dot(arc.axisParams.centre-p2)
            elOffsets.append(elOffset)
            # Determine the orientation, offset, and misalignment 
            # of the elevation axis
            elAngleOffset=-math.atan2( arcAxes[1].dot(arcu), arcAxes[0].dot(arcu) )
            for arcOrientation in arc.orientations:
                arcOrientation.setAngle(arcOrientation.angle+elAngleOffset)

            self.write("     Elevation arc {0}: orientation {1:.4f} offset {2:.4f} m  misalignment {3:.4f} degrees\n".format(
                arc.name,math.degrees(orientation),offset,math.degrees(misalignment)))

            arc.azimuthAngle=Orientation('Arc '+arc.name+' azimuth')
            arc.azimuthAngle.setAngle(orientation)

            # Not needed as replacing with elevation axis params, but just checking...
            arc.axisParams=AxisParams('Arc '+arc.name,np.vstack((arcu,arcv,arcw)),p2)

        elOffset=np.mean(elOffsets)
        for tgt in self.elevationTargetParams:
            tgt.offset += elOffset

        centre /= len(elevationArcs)
        offset=np.mean(offsets)
        misalignment=np.mean(misalignments)
        self.elevationAxisParams=ElevationAxisParams( self.name+' elevation axis',offset,misalignment)
        
        # Combine azimuth arcs onto common axes based on ARP

        axisParams.centre=centre
        self.write("\n Merging azimuth arcs\n")
        azimuthArcs.sort(key=lambda a: len(a.targets),reverse=True)
        targetParams={}
        for i in range(len(azimuthArcs)):
            if i > 0:
                azimuthArcs[i:].sort(
                    key=lambda a: a.countCommonTargets(targetParams),
                    reverse=True)
            arc=azimuthArcs[i]
            arc.alignAxes( axisParams )
            targets,maxError=arc.mergeTargetParams( targetParams )
            self.write("     Arc {0}: maximum offset merging targets {1:.4f}\n"
                       .format(arc.name,maxError))

        self.azimuthAxisParams=axisParams
        self.azimuthTargetParams=targetParams.values()
        self.azimuthTargetParams.sort(key=lambda t: t.name)

        # Now construct the elevation axis parameters..
        # These are 
        #    the azimuth orientation for each arc (defined as arc.azimuthAngle),
        #    the offset of the elevation axis from the azimuth axis
        #    the misalignment of the elevation axis with the normal to the azimuth axis

        for arc in elevationArcs:
            arc.axisParams=ElevationArcParams( arc.name+' elevation arc parameters',
                                              axisParams, self.elevationAxisParams, arc.azimuthAngle )

        # Now solve the equations
        self.write("\n Solving for all antenna parameters\n")

        tmadjustment=TargetMarkAdjustment( self.adjustment, self.targetMarks )
        self.initParams( tmadjustment )
        solved=tmadjustment.solve()

        self.axesDefined=True

        # Add the antenna IVP to the network as a station

        station=self.adjustment.stations().get( self.name )
        xyz=self.azimuthAxisParams.xyzCoordinates()[0]
        if station is None:
            self.adjustment.stations().addStation(
                Station( self.name, name=self.name+' IVP',xyz=xyz)
                )
        else:
            station.setXYZ(xyz)

        return solved

    def initParams( self, adjustment, calculate=True ):
        '''
        Initiallize all parameters for the antenna.  

        Returns the count of parameter rows in the observation equations, and the
        list of parameter objects used.
        '''
        # ARP and orientation of azimuth azis
        self.azimuthAxisParams.initParams(adjustment,calcAlignment=calculate,calcCentre=calculate)
        # Configuration of orientation axis
        self.elevationAxisParams.initParams(adjustment,calcAngle=calculate,calcOffset=calculate)
        # Orientation of azimuth axis for each elevation arc
        for arc in self.elevationArcs:
            arc.azimuthAngle.initParams(adjustment,calculate=calculate)

        # Target parameters for each target.  Note that for first target don't calculate
        # angle as otherwise ambiguity between angle and orientation parameters.
        for targetlist in (self.azimuthTargetParams,self.elevationTargetParams):
            targetlist[0].initParams(adjustment,calcOffset=calculate,calcRadius=calculate)
            for t in targetlist[1:]:
                t.initParams(adjustment,calcOffset=calculate,calcRadius=calculate,calcAngle=calculate)

        # Orientations around axis for each arc setting
        for arclist in (self.azimuthArcs,self.elevationArcs):
            for arc in arclist:
                for o in arc.orientations:
                    o.initParams( adjustment, calculate=calculate)

    def ivpCoordinates( self ):
        '''
        Return the coordinates of the centre and the 
        corresponding parameter ids
        '''
        return self.azimuthAxisParams.xyzCoordinates()

    def report( self, adjustment ):
        '''
        Report on antenna.  First report on antenna configuration.  This
        is defined by the centre and orientation of the azimuth axis, and the
        offset and misalignment of the elevation axis.
        '''

        azimuthParams=self.azimuthAxisParams
        elevationParams=self.elevationAxisParams
        
        centre=azimuthParams.centre
        uvw=azimuthParams.axes

        obseq=np.zeros((7,adjustment.nparam))
        for i in range(3):
            obseq[i,azimuthParams.params[i+2]]=1.0

        # Determine the deflection of azimuth axis from the vertical at centre
        lon,lat,hgt=GRS80.geodetic(centre)
        cenu=GRS80.enu_axes(lon,lat)
        radtosec=math.degrees(1.0)*3600
        de=uvw[2].dot(cenu[0])*radtosec
        dn=uvw[2].dot(cenu[1])*radtosec

        for i in range(2):
            obseq[i+3,azimuthParams.params[i]]=radtosec

        # Determine local deflection of the vertical.  Crudely done by using
        # target marks

        xieta=np.array([0.0,0.0])
        for m in self.targetMarks:
            xieta += m.mark.xieta()
        xieta *= 3600/len(self.targetMarks)

        offset=elevationParams.offset
        misalignment=elevationParams.angle*3600
        if offset < 0:
            offset=-offset
            misalignment=-misalignment

        obseq[5,elevationParams.offsetParam]=1.0
        obseq[6,elevationParams.angleParam]=3600.0

        prmcovar=adjustment.covariance()
        covar=obseq.dot(prmcovar.dot(obseq.T))
        seu=adjustment.seu
        covar *= seu*seu
        stderr=np.sqrt(covar.diagonal())

        self.write("\nAntenna invariant reference point\n")
        for i,prm in enumerate('XYZ'):
            self.write("     {0} coordinate: {1:13.4f} +/- {2:6.4f}\n".format(
                prm,centre[i],stderr[i]))

        self.write("\nAntenna azimuth axis deflection from geodetic vertical\n")
        self.write("   East deflection:  {0:6.2f} +/- {1:6.2f} seconds of arc\n"
                   .format(de,stderr[3]))
        self.write("   North deflection: {0:6.2f} +/- {1:6.2f} seconds of arc\n"
                   .format(dn,stderr[4]))
        self.write("   (Local geodetic deflection of vertical {0:.2f} {1:.2f})\n"
                   .format(xieta[1],xieta[0]))

        self.write("\nElevation axis configuration\n")
        self.write("    Offset from azimuth: {0:7.4f} +/- {1:6.4f} metres\n"
                  .format(offset,stderr[5]))
        self.write("    Non-orthogonality:   {0:5.2f}   +/- {1:4.2f}   seconds of arc\n"
                  .format(misalignment,stderr[6]))
        

        self.write("\nCorrelation of parameters\n")
        
        stderr=stderr.reshape((7,1))
        covar /= stderr
        covar /= stderr.T

        for i in range(7):
            self.write("    ")
            for j in range(7):
                self.write(" {0:6.3f}".format(covar[i,j]))
            self.write("\n")
        self.write("\n")

class AxisAdjustment( Plugin ):

    # Phases of adjustment
    NETWORK_ADJUSTMENT=1
    ANTENNA_ESTIMATION=2
    FINAL_ADJUSTMENT=3

    # Used for SINEX header
    softwareName='pyaxis software: Land Information New Zealand'

    def initPlugin( self ):
        self.axesDefined=False
        self.antennae={}
        self.targetMarks={}
        self.targetCalibrations={}
        self.antennaNprm=0
        self.phase=self.NETWORK_ADJUSTMENT

    def pluginOptions( self ):
        return dict(
            axisTargetRe=None,
            axisMaxIterations=10,
            axisConvergence=0.00001,
            skipPhase1Adjustment=False,
            antennaeArcs={},
            targetAdjustmentCsv=None,
            targetCalibrations={},
            writeTargetAdjustmentIterations=False,
            debugTargetAdjustmentTotalChange=False
        )

    def write( self, *params, **opts ):
        self.adjustment.write( *params, **opts )

    def stations( self ):
        return self.adjustment.stations
        
    def setConfigOption( self, item, value ):
        options=self.options
        if item == 'axis_target_re':
            options.axisTargetRe=re.compile(value)
        elif item == 'axis_max_iterations':
            options.axisMaxIterations=int(value)
        elif item == 'axis_convergence':
            options.axisConvergence=float(value)
        elif item == 'skip_phase1_adjustment':
            options.skipPhase1Adjustment=options.boolOption(value)
        elif item == 'antenna_arcs':
            parts=value.split()
            if len(parts) > 0:
                antname=parts.pop(0)
                options.antennaeArcs[antname]=parts
        elif item == 'target_adjustment_csv_file':
            options.targetAdjustmentCsv=value
        elif item == 'write_target_adjustment_iterations':
            options.writeTargetAdjustmentIterations=options.boolOption(value)
        elif item == 'debug_target_adjustment_change':
            options.debugTargetAdjustmentTotalChange=options.boolOption(value)

        elif item == 'target_calibration':
            parts=value.lower().split()
            if len(parts) < 2:
                raise RuntimeError("Invalid target_calibration: "+value)
            target=parts.pop(0)
            calculate=False
            correction=0.0
            if parts[0] == 'calculate':
                calculate=True
                parts.pop(0)
                if len(parts) == 1 and re.match(r'[+-]?\d+(\.\d+)$',parts[0]):
                    correction=float(parts[0])
                elif len(parts) > 0 or not calculate:
                    raise RuntimeError("Invalid target_calibration: "+value)
            options.targetCalibrations[target]=TargetCalibration(target,correction,calculate)
        else:
            return False
        return True
    
#        # Override adjustENU value
#        if self.adjustENU:
#            raise RuntimeError("adjust_enu is not a valid option in an axis adjustment")

    def defineAntennae( self ):
        if self.axesDefined:
            return
        tgtre=self.options.axisTargetRe
        # Create lookup from arc name to antenna name
        antennaeArcs=self.options.antennaeArcs
        arcAntenna={}
        for a in antennaeArcs:
            for arc in antennaeArcs[a]:
                arcAntenna[arc]=a
        
        # Assign marks to antennae targets

        targets=set()
        for code in self.adjustment.usedStations():
            m=tgtre.match(code)
            if m is not None:
                s=self.stations().get(code)
                arcName=m.group('arc')
                targetName=m.group('target')
                orientationName=m.group('angle')
                if arcName not in arcAntenna:
                    continue
                antennaName=arcAntenna[arcName]
                if antennaName not in self.antennae:
                    self.antennae[antennaName]=Antenna(antennaName,self)
                antenna=self.antennae[antennaName]
                antenna.addTargetMark( s, targetName, arcName, orientationName )
                self.targetMarks[code]=targetName
                targets.add(targetName)

        for target in targets:
            self.targetCalibrations[target]=(
                self.options.targetCalibrations.get(target,TargetCalibration(target)))
 
    def estimateAntennaeParameters( self ):
        solved=True
        for a in self.antennae.values():
            if not a.setup():
                solved=False
        if not solved:
            raise RuntimeError("Cannot estimate initial antennae parameters")

    def setupParameters( self ):
        '''
        Do nothing unless this is the final phase of the adjustment.
        In the final phase add the antenna parameters and target
        calibration parameters to the adjustment
        '''
        if self.phase==self.FINAL_ADJUSTMENT:
            for a in sorted(self.antennae.keys()):
                antenna=self.antennae[a]
                antenna.initParams(self.adjustment)
            
            for calibration in self.targetCalibrations.values():
                calibration.initParams(self.adjustment)

    def setupStationCoordMapping( self ):
        '''
        Do nothing unless this is the final phase of the adjustment.
        In the final phase set up the mapping from parameters to station coordinates. 
        '''
        adjustment=self.adjustment
        if self.phase==self.FINAL_ADJUSTMENT:
            mapping=adjustment.coordParamMapping
            for a in sorted(self.antennae.keys()):
                antenna=self.antennae[a]
                for m in antenna.targetMarks:
                    obseq=np.zeros((3,adjustment.nparam))
                    xyz=m.calculateXyz( obseq )
                    # Reset mark coordinate to value calculated from parameters
                    m.mark.setXYZ(xyz)
                    # Extract non-zero rows from obseq
                    prmnos=np.where((obseq != 0).any(axis=0))[0]
                    obsnonzero=obseq[:,prmnos].T
                    mapping[m.mark.code()]=(prmnos,obsnonzero,False)
                # Add coordinate mapping for antenna IVP so that
                # parameters can be output in SINEX file
                prmnos=antenna.ivpCoordinates()[1]
                mapping[a]=(prmnos,np.identity(3),False)

    def updateObservationEquation( self, obs, obseq ):
        if obs.obstype.code == 'SD':
            for i,o in enumerate(obs.obsvalues):
                target=self.targetMarks.get(o.trgtstn,None)
                if target is not None:
                    calibration=self.targetCalibrations[target]
                    obseq.obsres += calibration.correction
                    if self.phase==self.FINAL_ADJUSTMENT and calibration.prmno is not None:
                        obseq.obseq[i,calibration.prmno]=-1.0

    def writeHeader( self, header ):
        self.write("\n"+("="*70)+"\n")
        self.write(header)
        self.write("\n")

    def phase1( self ):
        self.defineAntennae()
        self.writeHeader("Phase 1: Unconstrained network adjustment\n")
        self.phase=self.NETWORK_ADJUSTMENT
        self.adjustment.calculateSolution()
        self.adjustment.writeSummaryStats()
        self.adjustment.writeResidualSummary()

    def phase2( self ):
        self.writeHeader("Phase 2: Antenna initial parameter estimation\n")
        self.phase=self.ANTENNA_ESTIMATION
        if len(self.antennae) > 0:
            self.estimateAntennaeParameters()
        else:
            self.write("*** No antennae defined in adjustment\n")
            self.write("    Check axis_target_re adjustment parameter\n")
            
    def phase3( self ):
        self.saveUnconstrainedCoords()
        self.writeHeader("Phase 3: Final constrained adjustment\n")
        self.phase=self.FINAL_ADJUSTMENT
        self.adjustment.calculateSolution()
        # Update IVP station coordinates
        for name,a in self.antennae.iteritems():
            xyz=a.ivpCoordinates()[0]
            self.stations().get(name).setXYZ(xyz)

    def saveUnconstrainedCoords( self ):
        self.unconstrainedCoords={}
        for a in self.antennae.values():
            for tm in a.targetMarks:
                self.unconstrainedCoords[tm.mark.code()]=np.array(tm.mark.xyz())

    def calculateSolution( self ):
        '''
        Replace adjustment calculate solution function with three phase calculation
        '''
        try:
            if not self.options.skipPhase1Adjustment:
                self.phase1()
            self.phase2()
            self.phase3()
        except Exception as ex:
            self.writeHeader("Axis adjustment failed")
            self.write(ex.message)
            raise
            sys.exit()
    
    def report( self ):
        listCalibrations=False
        for c in self.targetCalibrations.values():
            if c.calculate or c.correction != 0.0:
                listCalibrations=True
                break
        if listCalibrations:
            self.listTargetCalibrations()

        for a in sorted(self.antennae.keys()):
            self.writeHeader("Antenna "+a+" report")
            self.antennae[a].report(self.adjustment)

    def listTargetCalibrations( self ):
        self.writeHeader("Target calibrations")
        covariance=self.adjustment.covariance()
        seu=self.adjustment.seu
        self.write("\nTarget  calibration  error\n")
        for target in sorted(self.targetCalibrations):
            calibration=self.targetCalibrations[target]
            self.write("  {0:4s}  {1:8.4f}     "
                       .format(target,calibration.correction))
            prmno=calibration.prmno
            if prmno is None:
                self.write("  -")
            else:
                self.write("{0:5.4f}".format(math.sqrt(covariance[prmno,prmno])*seu))
            self.write("\n")

    def writeTargetAdjustments( self ):
        if self.options.targetAdjustmentCsv is None:
            return
        offsetsvalid=True
        for target in sorted(self.targetCalibrations):
            calibration=self.targetCalibrations[target]
            if abs(calibration.correction-calibration.initcorrection) > 0.0001:
                if offsetsvalid:
                    self.write("\n*** Cannot calculate target adjustments in {0}\n"
                               .format(self.options.targetAdjustmentCsv))
                    offsetsvalid=False
                self.write("   Target {0} value {1:.4f} differs from original value by more than 0.0001m\n"
                          .format(target,calibration.correction))
                self.write("   Set adjustment option: target_calibration {0} calculate {1:.4f}\n"
                          .format(target,calibration.correction))

        if not offsetsvalid:
            return

        with open(self.options.targetAdjustmentCsv,'w') as csvf:
            csvf.write("mark,antenna,axis,target,arc,de,dn,du\n")
            for a in sorted(self.antennae):
                antenna=self.antennae[a]
                centre,params=antenna.ivpCoordinates()
                lon,lat,hgt=GRS80.geodetic(centre)
                enu=GRS80.enu_axes(lon,lat)
                for m in antenna.targetMarks:
                    code=m.mark.code()
                    dxyz=m.mark.xyz()-self.unconstrainedCoords[code]
                    denu=enu.dot(dxyz)
                    axis='azimuth' if m.arc.isAzimuth else 'elevation'
                    csvf.write("{0},{1},{2},{3},{4},{5:.4f},{6:.4f},{7:.4f}\n".format(
                        code,a,axis,m.target.name,m.arc.name,denu[0],denu[1],denu[2]))

    def writeOutputFiles( self ):
        self.writeTargetAdjustments()

def main():
    from LINZ.Geodetic.LocalGeoidModelPlugin import LocalGeoidModelPlugin
    from LINZ.Geodetic.StationLocatorPlugin import StationLocatorPlugin
    from LINZ.Geodetic.SetupHeightPlugin import SetupHeightPlugin
    plugins=[LocalGeoidModelPlugin,StationLocatorPlugin,SetupHeightPlugin,AxisAdjustment]
    Adjustment.main(plugins=plugins,software=AxisAdjustment.softwareName)

if __name__=='__main__':
    main()

