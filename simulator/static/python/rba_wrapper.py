#!/usr/bin/env python
import sys
from .rba_session import SessionRBA

class Wrapper:

    def __init__(self):
        self.path = False
        self.Simulation = False

    def set_path(self, rpath):
        self.path = rpath

    def create_simulation(self):
        self.Simulation = SessionRBA(xml_dir=self.path,lp_solver="glpk")

    def set_default_parameters(self):
        ## Set medium to first condition, simulate and record results. ##
        self.Simulation.set_medium({'M_Glucose': 10})
        self.Simulation.find_max_growth_rate()
        self.Simulation.record_results('Glucose')

    def write_results(self):
        ## Write recorded results to RBA_SimulationData object ##
        self.Simulation.write_results(session_name='Test')

    def get_csv(self):
        ## Export results in CSV ##
        csv = self.Simulation.SimulationData.export_csv(return_result=True)
        return csv

    def get_eschermap(self):
        ## Export results as Escher Map ##
        em = self.Simulation.SimulationData.export_escher_map(type='fluxes',return_result=True)

        return em

    def get_proteomap(self):
        ## Export results in CSV ##
        proteomap = self.Simulation.SimulationData.export_proteo_map(return_result=True)

        return proteomap

    def get_sbtab(self):
        ## Export results in CSV ##
        sbtab = self.Simulation.SimulationData.export_sbtab(filename='Sbtab_Results_Glucose_Screen',return_result=True)
        return sbtab
