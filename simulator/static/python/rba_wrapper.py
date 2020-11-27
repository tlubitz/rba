#!/usr/bin/env python
import sys
from .rbastructure.NewControler import RBA_newControler

class Wrapper:

    def __init__(self):
        self.path = False
        self.Simulation = False

    def set_path(self, rpath):
        self.path = rpath
        print('Set path: ', self.path)

    def create_simulation(self):
        self.Simulation = RBA_newControler(self.path)
        print('Created Sim: ', self.Simulation)
    
    def set_default_parameters(self):
        ## Set medium to first condition, simulate and record results. ##
        self.Simulation.setMedium({'M_Glucose': 10})
        self.Simulation.findMaxGrowthRate()
        self.Simulation.recordResults('Glucose')
        print('Set default params: ', self.Simulation)

    def write_results(self):
        ## Write recorded results to RBA_SimulationData object ##
        self.Simulation.writeResults(session_name='Test')

    def get_csv(self):
        ## Export results in CSV ##
        csv = self.Simulation.SimulationData.exportCSV()
        return csv
        
    def get_eschermap(self):
        ## Export results as Escher Map ##
        em = self.Simulation.SimulationData.exportEscherMap(type='investment')
        print(em)
        sys.stdout.flush()
        print(em)
        return em