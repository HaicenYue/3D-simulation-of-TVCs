from PySteppables import *
import CompuCell
import sys

from PlayerPython import *
import CompuCellSetup
from math import *
import numpy

theta1 = 0.0
theta2 = 0.0
omega1 = 0.005
omega2 = 0.1
sigma = 0.1
index = 0.0
class DirectionUpdate(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
    def start(self):
        print "DirectionUpdate: Tis function is called once before simulation"
        global theta1
        global theta2
        theta1 = 0.0
        theta2 = 0.0
    def step(self,mcs):
        randn = numpy.random.normal(0.0, 1.0, 2)
        global theta1
        global theta2
        for cell in self.cellListByType(self.LEADER):  
            xCOM_L = cell.xCOM
            yCOM_L = cell.yCOM
            zCOM_L = cell.zCOM
        for cell in self.cellListByType(self.TRAILER):  
            xCOM_T = cell.xCOM
            yCOM_T = cell.yCOM
            zCOM_T = cell.zCOM
        distbt = math.sqrt((xCOM_L-xCOM_T)**2+(yCOM_L-yCOM_T)**2)
        costhe0 = (xCOM_L-xCOM_T)/distbt
        sinthe0 = (yCOM_L-yCOM_T)/distbt
        if sinthe0<0:
            anglelink = -math.acos(costhe0)
        else:
            anglelink = math.acos(costhe0)
        if theta1>math.pi:
            theta1 = theta1-2*math.pi
        if theta2>math.pi:
            theta2 = theta2-2*math.pi
        if theta1<-math.pi:
            theta1 = theta1+2*math.pi
        if theta2<-math.pi:
            theta2 = theta2+2*math.pi
        theta1 = theta1 - omega1*theta1 + sigma*randn[0]
        theta2 = theta2 - omega2*(theta2-anglelink) + sigma*randn[1]
    def finish(self):
#         self.pW1.savePlotAsPNG("polarity angle.png", 1000,1000)  
        return
    

from PySteppables import *
import CompuCell
import sys

from PlayerPython import *
import CompuCellSetup
from math import *
import numpy


class UpdateForce(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
    def start(self):
        pass
    def step(self,mcs):
        for cell in self.cellListByType(self.LEADER):  
            xCOM_L = cell.xCOM
            yCOM_L = cell.yCOM
            zCOM_L = cell.zCOM
            Pro_L = self.chemotaxisPlugin.addChemotaxisData(cell, "PROTRUSION_L")
            Ret_L = self.chemotaxisPlugin.addChemotaxisData(cell, "RETRACTION_L")
            Pro_L.setLambda(180.0)
            Ret_L.setLambda(20.0) 
        for cell in self.cellListByType(self.TRAILER):  
            xCOM_T = cell.xCOM
            yCOM_T = cell.yCOM
            zCOM_T = cell.zCOM
            Pro_T = self.chemotaxisPlugin.addChemotaxisData(cell, "PROTRUSION_T")
            Ret_T = self.chemotaxisPlugin.addChemotaxisData(cell, "RETRACTION_T")
            Pro_T.setLambda(150.0)
            Ret_T.setLambda(70.0)   
        field1 = self.getConcentrationField('PROTRUSION_L')
        field2 = self.getConcentrationField('RETRACTION_L')
        field3 = self.getConcentrationField('PROTRUSION_T')
        field4 = self.getConcentrationField('RETRACTION_T')
        for x, y, z in self.everyPixel():
            d_L = math.sqrt((x-xCOM_L)**2+(y-yCOM_L)**2+(z-10.0)**2)
            r_L = math.sqrt((x-xCOM_L)**2+(y-yCOM_L)**2)
            r_LT = math.sqrt((xCOM_L-xCOM_T)**2+(yCOM_L-yCOM_T)**2)
            cos_betaL = (x-xCOM_L)/r_L
            sin_betaL = (y-yCOM_L)/r_L
            cos_L = cos_betaL*math.cos(theta1)+sin_betaL*math.sin(theta1)
            length_L = 1.0
            field1[x, y, z] = r_L*((1.0+math.tanh((cos_L)/0.01))/2.0)
            field2[x, y, z] = -d_L*((1.0+math.tanh(-(cos_L)/0.01))/2.0)
            d_T = math.sqrt((x-xCOM_T)**2+(y-yCOM_T)**2+(z-10.0)**2)
            r_T = math.sqrt((x-xCOM_T)**2+(y-yCOM_T)**2)
            cos_betaT = (x-xCOM_T)/r_T
            sin_betaT = (y-yCOM_T)/r_T
            cos_T = cos_betaT*math.cos(theta2)+sin_betaT*math.sin(theta2)
            length_T = 1.0
            field3[x, y, z] = r_T*((1.0+math.tanh((cos_T)/0.01))/2.0)*cos_T**2
            field4[x, y, z] = -d_T*((1.0+math.tanh(-(cos_T)/0.01))/2.0)
    def finish(self):
        # this function may be called at the end of simulation - used very infrequently though
        return
    

from PySteppables import *
import CompuCell
import sys
import math
import numpy
class TwoCellChemoSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        pass
    def step(self,mcs): 
        try:
            fileHandle1, fullFileName = self.openFileInSimulationOutputDirectory("output.txt", "a")
        except IOError:
            print "Could not open file ", "YOUR FILE NAME", " for writing. "
            return
        for cell in self.cellListByType(self.LEADER):  
            xCOM_L = cell.xCOM
            yCOM_L = cell.yCOM
            zCOM_L = cell.zCOM
            Volume_L = cell.volume
            Radius_L = (Volume_L*3./4/math.pi)**(1./3)
            Surface_L = cell.surface
            Sphericity_L = 4*math.pi*Radius_L**2/Surface_L
            EpidemisContact_L = 0
            Contact = 0
            for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                if neighbor:
                    if neighbor.type == 1:
                        EpidemisContact_L += commonSurfaceArea
                    if neighbor.type == 3:
                        Contact += commonSurfaceArea
            print >> fileHandle1,mcs, xCOM_L, yCOM_L, zCOM_L, Sphericity_L, EpidemisContact_L, Contact,Surface_L,Volume_L
        for cell in self.cellListByType(self.TRAILER):  
            xCOM_T = cell.xCOM
            yCOM_T = cell.yCOM
            zCOM_T = cell.zCOM
            Volume_T = cell.volume
            Radius_T = (Volume_T*3./4/math.pi)**(1./3)
            Surface_T = cell.surface
            Sphericity_T = 4*math.pi*Radius_T**2/Surface_T
            EpidemisContact_T = 0
            Contact = 0
            for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                if neighbor:
                    if neighbor.type == 1:
                        EpidemisContact_T += commonSurfaceArea
                    if neighbor.type == 3:
                        Contact += commonSurfaceArea
            print >> fileHandle1, mcs,xCOM_T, yCOM_T, zCOM_T, Sphericity_T, EpidemisContact_T, Contact,Surface_T,Volume_T
    def finish(self):
        pass
        


