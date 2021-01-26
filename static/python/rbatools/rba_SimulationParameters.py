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
from sbtab import SBtab

# package imports
import rba
from .parameter_block import ParameterBlock


class RBA_SimulationParameters(object):
    """
    Class holding information on simulations with the model.

    Attributes
    ----------
    EnzymeCapacities : rbatools.ParameterBlock object holding enzyme capacities.

    ProcessCapacities : rbatools.ParameterBlock object holding process capacities.

    CompartmentCapacities : rbatools.ParameterBlock object holding compartment capacities.


    Methods
    ----------
    __init__:

    fromSimulationResults:
    Arguments: 'Controller'= rbatools.Controler object

    exportSBtab_OneFile

    """

    def __init__(self, StaticData):

        self.EnzymeCapacities_FW = ParameterBlock()
        self.EnzymeCapacities_BW = ParameterBlock()
        self.ProcessCapacities = ParameterBlock()
        self.CompartmentCapacities = ParameterBlock()
        self.TargetValues = ParameterBlock()      
        self.Media = ParameterBlock()
        self.EnzymeCapacities_FW.fromDict({})
        self.EnzymeCapacities_BW.fromDict({})
        self.ProcessCapacities.fromDict({})
        self.CompartmentCapacities.fromDict({})
        self.Media.fromDict({})        
        self.TargetValues.fromDict({})

    def fromSimulationResults(self, Controller):

        for effi in list(Controller.Parameters['EnzymeEfficiencies_FW'].index):
            if effi not in self.EnzymeCapacities_FW.Elements:
                self.EnzymeCapacities_FW.Elements.update({effi: {}})
            self.EnzymeCapacities_FW.Elements[effi].update({'ID': effi})
            for run in list(Controller.Parameters['EnzymeEfficiencies_FW']):
                self.EnzymeCapacities_FW.Elements[effi].update(
                    {run: Controller.Parameters['EnzymeEfficiencies_FW'].loc[effi, run]})

        for effi in list(Controller.Parameters['EnzymeEfficiencies_BW'].index):
            if effi not in self.EnzymeCapacities_BW.Elements:
                self.EnzymeCapacities_BW.Elements.update({effi: {}})
            self.EnzymeCapacities_BW.Elements[effi].update({'ID': effi})
            for run in list(Controller.Parameters['EnzymeEfficiencies_BW']):
                self.EnzymeCapacities_BW.Elements[effi].update(
                    {run: Controller.Parameters['EnzymeEfficiencies_BW'].loc[effi, run]})

        for procapa in list(Controller.Parameters['ProcessEfficiencies'].index):
            if procapa not in self.ProcessCapacities.Elements:
                self.ProcessCapacities.Elements.update({procapa: {}})
            self.ProcessCapacities.Elements[procapa].update({'ID': procapa})
            for run in list(Controller.Parameters['ProcessEfficiencies']):
                self.ProcessCapacities.Elements[procapa].update(
                    {run: Controller.Parameters['ProcessEfficiencies'].loc[procapa, run]})

        for compcap in list(Controller.Parameters['CompartmentCapacities'].index):
            if compcap not in self.CompartmentCapacities.Elements:
                self.CompartmentCapacities.Elements.update({compcap: {}})
            self.CompartmentCapacities.Elements[compcap].update({'ID': compcap})
            for run in list(Controller.Parameters['CompartmentCapacities']):
                self.CompartmentCapacities.Elements[compcap].update(
                    {run: Controller.Parameters['CompartmentCapacities'].loc[compcap, run]})
        
        for species in list(Controller.Parameters['Medium'].index):
            if species not in self.Media.Elements:
                self.Media.Elements.update({species: {}})
            self.Media.Elements[species].update({'ID': species})
            for run in list(Controller.Parameters['Medium']):
                self.Media.Elements[species].update(
                    {run: Controller.Parameters['Medium'].loc[species, run]})
        
        for targetvalue in list(Controller.Parameters['TargetValues'].index):
            if targetvalue not in self.TargetValues.Elements:
                self.TargetValues.Elements.update({targetvalue: {}})
            self.TargetValues.Elements[targetvalue].update({'ID': targetvalue})
            for run in list(Controller.Parameters['TargetValues']):
                self.TargetValues.Elements[targetvalue].update(
                    {run: Controller.Parameters['TargetValues'].loc[targetvalue, run]})      

    def exportSBtab(self, filename_SBtab):

        EnzymeCapacitiesTable_FW = self.EnzymeCapacities_FW.toSBtab(
            table_id='enzyme_forward_capacity', table_type='QuantityMatrix', table_name='Enzyme forward-capacities')
        EnzymeCapacitiesTable_BW = self.EnzymeCapacities_BW.toSBtab(
            table_id='enzyme_backward_capacity', table_type='QuantityMatrix', table_name='Enzyme backward-capacities')
        ProcessCapacitiesTable = self.ProcessCapacities.toSBtab(
            table_id='Machine_capacity', table_type='QuantityMatrix', table_name='Machine capacities')
        CompartmentCapacitiesTable = self.CompartmentCapacities.toSBtab(
            table_id='compartment_capacity', table_type='QuantityMatrix', table_name='Compartment capacities')
        TargetValuesTable = self.TargetValues.toSBtab(
            table_id='target_values', table_type='QuantityMatrix', table_name='Target values')            
        MediaConcentrationTable = self.Media.toSBtab(
            table_id='medium_composition', table_type='QuantityMatrix', table_name='Medium composition')

        TargetValuesTable.filename = 'TargetValues.tsv'
        TargetValuesTable.change_attribute('QuantityType', 'target_value')
        TargetValuesTable.change_attribute('Unit', '')
        TargetValuesTable.change_attribute(
            'Text', 'Values of cellular targets, at given growth-rate and medium.')
        TargetValuesTable.change_attribute('TableID', 'target_values')
        TargetValuesTable.unset_attribute('Date')
        TargetValuesTable.unset_attribute('SBtabVersion')

        EnzymeCapacitiesTable_FW.filename = 'EnzymeForwardCapacities.tsv'
        EnzymeCapacitiesTable_FW.change_attribute('QuantityType', 'enzyme_capacity')
        EnzymeCapacitiesTable_FW.change_attribute('Unit', '1/h')
        EnzymeCapacitiesTable_FW.change_attribute(
            'Text', 'Enzyme forward capacities (relating enzyme concentrations to reaction rates) may depend on model parameters such as growth rate and will therefore vary between simulation runs.')
        EnzymeCapacitiesTable_FW.change_attribute('TableID', 'enzyme_forward_capacity')
        EnzymeCapacitiesTable_FW.unset_attribute('Date')
        EnzymeCapacitiesTable_FW.unset_attribute('SBtabVersion')

        EnzymeCapacitiesTable_BW.filename = 'EnzymeBackwardCapacities.tsv'
        EnzymeCapacitiesTable_BW.change_attribute('QuantityType', 'enzyme_capacity')
        EnzymeCapacitiesTable_BW.change_attribute('Unit', '1/h')
        EnzymeCapacitiesTable_BW.change_attribute(
            'Text', 'Enzyme backward capacities (relating enzyme concentrations to reaction rates) may depend on model parameters such as growth rate and will therefore vary between simulation runs.')
        EnzymeCapacitiesTable_BW.change_attribute('TableID', 'enzyme_backward_capacity')
        EnzymeCapacitiesTable_BW.unset_attribute('Date')
        EnzymeCapacitiesTable_BW.unset_attribute('SBtabVersion')

        ProcessCapacitiesTable.filename = 'MachineCapacities.tsv'
        ProcessCapacitiesTable.change_attribute('QuantityType', 'machinery_capacity')
        ProcessCapacitiesTable.change_attribute('Unit', '1/h')
        ProcessCapacitiesTable.change_attribute(
            'Text', 'Machine capacities (relating machine concentrations to process rates)  may depend on model parameters such as growth rate and will therefore vary between simulation runs.')
        ProcessCapacitiesTable.change_attribute('TableID', 'machinery_capacity')
        ProcessCapacitiesTable.unset_attribute('Date')
        ProcessCapacitiesTable.unset_attribute('SBtabVersion')

        CompartmentCapacitiesTable.filename = 'CompartmentCapacities.tsv'
        CompartmentCapacitiesTable.change_attribute('QuantityType', 'compartment_capacity')
        CompartmentCapacitiesTable.change_attribute('Unit', 'mmol/gDW')
        CompartmentCapacitiesTable.change_attribute(
            'Text', 'Compartment capacities (defining the maximal macromolecule density in a compartment)  may depend on model parameters such as growth rate and will therefore vary between simulation runs.')
        CompartmentCapacitiesTable.change_attribute('TableID', 'compartment_capacity')
        CompartmentCapacitiesTable.unset_attribute('Date')
        CompartmentCapacitiesTable.unset_attribute('SBtabVersion')

        MediaConcentrationTable.filename = 'CompartmentCapacities.tsv'
        MediaConcentrationTable.change_attribute('QuantityType', 'concentration')
        MediaConcentrationTable.change_attribute('Unit', 'mmol/l')
        MediaConcentrationTable.change_attribute('Text', 'Composition of growth medium in mmol/l.')
        MediaConcentrationTable.change_attribute('TableID', 'medium_composition')
        MediaConcentrationTable.unset_attribute('Date')
        MediaConcentrationTable.unset_attribute('SBtabVersion')

        if filename is not None:
            filename_SBtab = filename
        else:
            filename_SBtab = 'RBA_parameters'

        if add_links:
            EnzymeCapacitiesTable_FW.add_column(column_list=['!ElementID']+[str('(!'+'Enzyme/'+entry+'!)')
                                                                            for entry in list(EnzymeCapacitiesTable_FW.to_data_frame()['ID'])], position=1)
            EnzymeCapacitiesTable_BW.add_column(column_list=['!ElementID']+[str('(!'+'Enzyme/'+entry+'!)')
                                                                            for entry in list(EnzymeCapacitiesTable_BW.to_data_frame()['ID'])], position=1)
            ProcessCapacitiesTable.add_column(column_list=['!ElementID']+[str('(!'+'Process/'+entry+'!)')
                                                                          for entry in list(ProcessCapacitiesTable.to_data_frame()['ID'])], position=1)
            CompartmentCapacitiesTable.add_column(column_list=['!ElementID']+[str('(!'+'Compartment/'+entry+'!)')
                                                                              for entry in list(CompartmentCapacitiesTable.to_data_frame()['ID'])], position=1)
            MediaConcentrationTable.add_column(column_list=['!ElementID']+[str('(!'+'Compound/'+entry+'!)')
                                                                           for entry in list(MediaConcentrationTable.to_data_frame()['ID'])], position=1)
            TargetValuesTable.add_column(column_list=['!ElementID']+[str('(!'+'CellTarget/'+entry+'!)')
                                                                     for entry in list(TargetValuesTable.to_data_frame()['ID'])], position=1)            
            filename_SBtab += '_HTML'
        else:
            EnzymeCapacitiesTable_FW.add_column(
                column_list=['!ElementID']+list(EnzymeCapacitiesTable_FW.to_data_frame()['ID']), position=1)
            EnzymeCapacitiesTable_BW.add_column(
                column_list=['!ElementID']+list(EnzymeCapacitiesTable_BW.to_data_frame()['ID']), position=1)
            ProcessCapacitiesTable.add_column(
                column_list=['!ElementID']+list(ProcessCapacitiesTable.to_data_frame()['ID']), position=1)
            CompartmentCapacitiesTable.add_column(
                column_list=['!ElementID']+list(CompartmentCapacitiesTable.to_data_frame()['ID']), position=1)
            MediaConcentrationTable.add_column(
                column_list=['!ElementID']+list(MediaConcentrationTable.to_data_frame()['ID']), position=1)
            TargetValuesTable.add_column(
                column_list=['!ElementID']+list(TargetValuesTable.to_data_frame()['ID']), position=1)

        EnzymeCapacitiesTable_FW.remove_column(position=2)
        EnzymeCapacitiesTable_BW.remove_column(position=2)
        ProcessCapacitiesTable.remove_column(position=2)
        CompartmentCapacitiesTable.remove_column(position=2)
        TargetValuesTable.remove_column(position=2)
        MediaConcentrationTable.remove_column(position=2)

        Out = SBtab.SBtabDocument(name='RBAparameters', sbtab_init=None,
                                  filename=str(filename_SBtab+'.tsv'))

        Out.add_sbtab(MediaConcentrationTable)
        Out.add_sbtab(EnzymeCapacitiesTable_FW)
        Out.add_sbtab(EnzymeCapacitiesTable_BW)
        Out.add_sbtab(ProcessCapacitiesTable)
        Out.add_sbtab(CompartmentCapacitiesTable)
        Out.add_sbtab(TargetValuesTable)

        Out.change_attribute('DocumentName', 'RBA parameters')
        Out.name = filename
        Out.change_attribute('DocumentType', 'rba-simulation-parameters')
        Out.write()
