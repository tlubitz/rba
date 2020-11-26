# python 2/3 compatibility
from __future__ import division, print_function
import sys
import os.path
import numpy
import pandas
import copy
import json
import jxmlease
import xml.etree.ElementTree as ET

# package imports
import rba
from rbastructure.model_Data import ModelData
from rbastructure.information_block import InformationBlock
from rbastructure.description_block import DescriptionBlock


class RBA_SimulationParameters(ModelData):
    """
    Class holding information on simulations with the model.
    Inherits from from rbastructure.model_Data ModelData

    Attributes
    ----------
    EnzymeCapacities : rbastructure.InformationBlock object holding enzyme capacities.

    ProcessCapacities : rbastructure.InformationBlock object holding process capacities.

    CompartmentCapacities : rbastructure.InformationBlock object holding compartment capacities.


    Methods
    ----------
    __init__:

    fromSimulationResults:
    Arguments: 'Controller'= rbastructure.Controler object

    exportSBtab_OneFile

    """

    def __init__(self, StaticData):

        self.EnzymeCapacities = InformationBlock()
        self.ProcessCapacities = InformationBlock()
        self.CompartmentCapacities = InformationBlock()
        self.EnzymeCapacities.fromDict({})
        self.ProcessCapacities.fromDict({})
        self.CompartmentCapacities.fromDict({})

    def fromSimulationResults(self, Controller):

        for effi in list(Controller.Parameters['EnzymeEfficiencies'].index):
            if effi not in self.EnzymeCapacities.Elements:
                self.EnzymeCapacities.Elements.update({effi: {}})
            self.EnzymeCapacities.Elements[effi].update({'ID': effi})
            for run in list(Controller.Parameters['EnzymeEfficiencies']):
                self.EnzymeCapacities.Elements[effi].update(
                    {run: json.dumps(Controller.Parameters['EnzymeEfficiencies'].loc[effi, run])})

        for procapa in list(Controller.Parameters['NetProcessEfficiencies'].index):
            if procapa not in self.ProcessCapacities.Elements:
                self.ProcessCapacities.Elements.update({procapa: {}})
            self.ProcessCapacities.Elements[procapa].update({'ID': procapa})
            for run in list(Controller.Parameters['NetProcessEfficiencies']):
                self.ProcessCapacities.Elements[procapa].update(
                    {run: json.dumps(Controller.Parameters['NetProcessEfficiencies'].loc[procapa, run])})

        for compcap in list(Controller.Parameters['CompartmentCapacities'].index):
            if compcap not in self.CompartmentCapacities.Elements:
                self.CompartmentCapacities.Elements.update({compcap: {}})
            self.CompartmentCapacities.Elements[compcap].update({'ID': compcap})
            for run in list(Controller.Parameters['CompartmentCapacities']):
                self.CompartmentCapacities.Elements[compcap].update(
                    {run: json.dumps(Controller.Parameters['CompartmentCapacities'].loc[compcap, run])})

    def exportSBtab_OneFile(self):
        from sbtab import SBtab

        EnzymeCapacitiesTable = self.EnzymeCapacities.toSBtab_forDoc(
            'RBA_parameters', 'Quantity', 'Enzyme Capacities', 'RBAparameters', 'Stored RBAparameters')
        ProcessCapacitiesTable = self.ProcessCapacities.toSBtab_forDoc(
            'RBA_parameters', 'Quantity', 'Process Capacities', 'RBAparameters', 'Stored RBAparameters')
        CompartmentCapacitiesTable = self.CompartmentCapacities.toSBtab_forDoc(
            'RBA_parameters', 'Quantity', 'Compartment Capacities', 'RBAparameters', 'Stored RBAparameters')

        EnzymeCapacitiesTable.filename = 'EnzymeCapacities.tsv'
        EnzymeCapacitiesTable.change_attribute('Text', 'Information on RBA enzyme efficiencies')
        EnzymeCapacitiesTable.unset_attribute('Unit')
        EnzymeCapacitiesTable.change_attribute('TableID', 'EnzymeCapacity')
        EnzymeCapacitiesTable.unset_attribute('DocumentName')
        EnzymeCapacitiesTable.unset_attribute('Document')
        EnzymeCapacitiesTable.unset_attribute('Date')
        EnzymeCapacitiesTable.unset_attribute('SBtabVersion')

        ProcessCapacitiesTable.filename = 'ProcessCapacities.tsv'
        ProcessCapacitiesTable.change_attribute('Text', 'Information on RBA process efficiencies')
        ProcessCapacitiesTable.unset_attribute('Unit')
        ProcessCapacitiesTable.change_attribute('TableID', 'ProcessCapacity')
        ProcessCapacitiesTable.unset_attribute('DocumentName')
        ProcessCapacitiesTable.unset_attribute('Document')
        ProcessCapacitiesTable.unset_attribute('Date')
        ProcessCapacitiesTable.unset_attribute('SBtabVersion')

        CompartmentCapacitiesTable.filename = 'CompartmentCapacities.tsv'
        CompartmentCapacitiesTable.change_attribute(
            'Text', 'Information on RBA compartment capacities')
        CompartmentCapacitiesTable.unset_attribute('Unit')
        CompartmentCapacitiesTable.change_attribute('TableID', 'CompartmentCapacity')
        CompartmentCapacitiesTable.unset_attribute('DocumentName')
        CompartmentCapacitiesTable.unset_attribute('Document')
        CompartmentCapacitiesTable.unset_attribute('Date')
        CompartmentCapacitiesTable.unset_attribute('SBtabVersion')

        Out = SBtab.SBtabDocument('RBAparameters', EnzymeCapacitiesTable, 'RBA_parameters.tsv')
        Out.add_sbtab(ProcessCapacitiesTable)
        Out.add_sbtab(CompartmentCapacitiesTable)

        Out.change_attribute('DocumentName', 'RBA parameters')
        Out.name = 'RBA_parameters'
        Out.write()
