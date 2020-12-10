#!/usr/bin/env python
import sys
from .rba_Session import RBA_Session

class Wrapper:

    def __init__(self):
        self.path = False
        self.Simulation = False

    def set_path(self, rpath):
        self.path = rpath

    def create_simulation(self):
        self.Simulation = RBA_Session(self.path)
    
    def set_default_parameters(self):
        ## Set medium to first condition, simulate and record results. ##
        self.Simulation.setMedium({'M_Glucose': 10})
        self.Simulation.findMaxGrowthRate()
        self.Simulation.recordResults('Glucose')

    def write_results(self):
        ## Write recorded results to RBA_SimulationData object ##
        self.Simulation.writeResults(session_name='Test')

    def get_csv(self):
        ## Export results in CSV ##
        csv = self.Simulation.SimulationData.exportCSV()
        return csv
        
    def get_eschermap(self):
        ## Export results as Escher Map ##
        self.Simulation.SimulationData.exportEscherMap(type='fluxes')
        em = self.Simulation.SimulationData.getEscherMap()

        return em

    def get_proteomap(self):
        ## Export results in CSV ##
        self.Simulation.SimulationData.exportProteoMap()
        proteomap = self.Simulation.SimulationData.getProteoMap()
        
        return proteomap

    def get_sbtab(self):
        ## Export results in CSV ##
        self.Simulation.SimulationData.exportSBtab(filename_SBtab='Sbtab_Results_Glucose_Screen')
        sbtab = self.Simulation.SimulationData.getSBtabDoc()
        return sbtab

