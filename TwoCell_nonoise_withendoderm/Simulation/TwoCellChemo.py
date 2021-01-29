
import sys
from os import environ
from os import getcwd
import string

sys.path.append(environ["PYTHON_MODULE_PATH"])

import CompuCellSetup

sim,simthread = CompuCellSetup.getCoreSimulationObjects()
        
# add extra attributes here
        
CompuCellSetup.initializeSimulationObjects(sim,simthread)
CompuCellSetup.doNotOutputField("PROTRUSION_L")
CompuCellSetup.doNotOutputField("RETRACTION_L")
CompuCellSetup.doNotOutputField("PROTRUSION_T")
CompuCellSetup.doNotOutputField("RETRACTION_T")
# Definitions of additional Python-managed fields go here
        
#Add Python steppables here
steppableRegistry=CompuCellSetup.getSteppableRegistry()
        
from TwoCellChemoSteppables import TwoCellChemoSteppable
steppableInstance=TwoCellChemoSteppable(sim,_frequency=10)
steppableRegistry.registerSteppable(steppableInstance)
        

# from TwoCellChemoSteppables import DirectionUpdate
# instanceOfDirectionUpdate=DirectionUpdate(_simulator=sim,_frequency=1)
# steppableRegistry.registerSteppable(instanceOfDirectionUpdate)


from TwoCellChemoSteppables import UpdateForce
instanceOfUpdateForce=UpdateForce(_simulator=sim,_frequency=10)
steppableRegistry.registerSteppable(instanceOfUpdateForce)

CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)
        
        