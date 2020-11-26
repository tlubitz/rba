#!/usr/bin/env python
import sys
#from rbastructure.NewControler import RBA_newControler
#import rbastructure.NewControler
print('START2')
from rbastructure.NewControler import RBA_newControler

#model_name = sys.argv[1]
#path = 'models/' + model_name
path = sys.argv[1]
print("le path: ", path)
Simulation = RBA_newControler(path)

## Set medium to first condition, simulate and record results.##
Simulation.setMedium({'M_Glucose': 10})
Simulation.findMaxGrowthRate()
Simulation.recordResults('Glucose')


## Write recorded results to RBA_SimulationData object ##
Simulation.writeResults(session_name='LeTest')

## Export results in various formats ##
Simulation.SimulationData.exportCSV()
Simulation.SimulationData.exportEscherMap(type='investment')
print('result!')

sys.stdout.flush()
