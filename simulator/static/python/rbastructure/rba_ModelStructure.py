# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import sys
import os.path
import numpy
import pandas
import copy
import json
import jxmlease
#import xml.etree.ElementTree as ET
import libsbml
# package imports
import rba
from ..rba.core.constraint_blocks import ConstraintBlocks

from .model_Data import ModelData
from .infoMatrices import InfoMatrices
from .description_block import DescriptionBlock
from .metabolite_block import MetaboliteBlock
from .module_block import ModuleBlock
from .process_block import ProcessBlock
from .reaction_block import ReactionBlock
from .enzyme_block import EnzymeBlock
from .protein_block import ProteinBlock
from .compartment_block import CompartmentBlock
from .metabolite_constraints import MetaboliteConstraints
from .density_constraints import DensityConstraints
from .process_constraints import ProcessConstraints
from .enzyme_constraints import EnzymeConstraints
from .statistics_block import StatisticsBlock
from .target_block import TargetBlock


class RBA_ModelStructure(ModelData):
    """
    Class holding information on model-structure.

    Attributes
    ----------
    GeneralInfo : rbastructure.description_block.Description_block
         Model description
    MetaboliteInfo : rbastructure.metabolite_block.Metabolite_block
         Metabolite information
    ModuleInfo : rbastructure.module_block.Module_block
         Module information
    ProcessInfo : rbastructure.process_block.Process_block
         Process information
    ReactionInfo : rbastructure.reaction_block.Reaction_block
         Reaction information
    EnzymeInfo : rbastructure.enzyme_block.Enzyme_block
         Enzyme information
    ProteinInfo : rbastructure.protein_block.Protein_block
         Protein information
    CompartmentInfo : rbastructure.compartment_block.Compartment_block
         Compartment information
    MetaboliteConstraintsInfo : rbastructure.metabolite_constraints.Metabolite_constraints
         Metabolite-constraint information
    DensityConstraintsInfo : rbastructure.density_constraints.Density_constraints
         Density-constraint information
    ProcessConstraintsInfo : rbastructure.process_constraints.Process_constraints
         Process-constraint information
    EnzymeConstraintsInfo : rbastructure.enzyme_constraints.Enzyme_constraints
         Enzyme-constraint information
    ModelStatistics : rbastructure.statistics_block.Statistics_block
         Statistics on Model
    TargetInfo : rbastructure.statistics_block.Statistics_block
         Target information
    ProteinMatrix : numpy.array
        Matrix mapping proteins to consumers (enzymes and process-machineries)
    MediumDependencies : dict
        Dictionary with boundary metabolites as keys
        and a list of constraint_IDs, which are affected by its concentration, as values.
    MuDependencies : list
        List of constraint_IDs, which are affected by the growth-rate

    Methods
    ----------
    fromFiles:
    fromJSON:
    fromSBtab:
    generateMatrices:
    toXML:
    exportJSON:
    exportSBtab:
    exportSBtab_OneFile:
    """

    def fromFiles(self, xml_dir):

        UniprotFile = importUniprotFile(xml_dir)
        GeneMap = importGeneAnnotations(xml_dir)
        Info = importModelInfo(xml_dir)
        SBMLfile = str('Not There')
        if Info['Value']['SBML-file'] != 'Not Provided':
            SBMLfile = importSbmlFile(xml_dir, str(Info['Value']['SBML-file']))

        MetaboliteAnnotations = importMetaboliteAnnotations(xml_dir)
        ReactionAnnotations = importReactionAnnotations(xml_dir)

        print('')
        print('Generating model-structure')
        print('...')

        model = rba.RbaModel.from_xml(xml_dir)
        Zero_matrix = rba.ConstraintMatrix(model)
        Zero_matrix.build_matrices(0)
        constraints = sortConstraints(Zero_matrix, model)

        MetaboliteInfo = MetaboliteBlock()
        ModuleInfo = ModuleBlock()
        ProcessInfo = ProcessBlock()
        ReactionInfo = ReactionBlock()
        EnzymeInfo = EnzymeBlock()
        ProteinInfo = ProteinBlock()
        CompartmentInfo = CompartmentBlock()
        TargetInfo = TargetBlock()

        TargetInfo.fromFiles(model)
        MetaboliteInfo.fromFiles(model, Info, MetaboliteAnnotations, SBMLfile)
        ModuleInfo.fromFiles(model, SBMLfile)
        ProcessInfo.fromFiles(model, Info)
        ReactionInfo.fromFiles(model, Info, ReactionAnnotations, SBMLfile, MetaboliteInfo)
        EnzymeInfo.fromFiles(model, Info)
        ProteinInfo.fromFiles(model, GeneMap, Info, UniprotFile)
        CompartmentInfo.fromFiles(model, Info)

        self.GeneralInfo = DescriptionBlock()
        self.GeneralInfo.fromFiles(Info)
        self.MetaboliteInfo = MetaboliteInfo
        self.ModuleInfo = ModuleInfo
        self.ProcessInfo = ProcessInfo
        self.TargetInfo = TargetInfo

        for protein in ProteinInfo.Elements.keys():
            AssociatedEnzyme = findAssociatedEnzyme(EnzymeInfo.Elements, protein)
            ProteinInfo.Elements[protein]['associatedEnzymes'] = AssociatedEnzyme['Enz']
            ProteinInfo.Elements[protein]['associatedReactions'] = AssociatedEnzyme['Rx']

        CB = ConstraintBlocks(model)
        for enzyme in EnzymeInfo.Elements.keys():
            EnzymeInfo.Elements[enzyme]['Isozymes'] = findIsozymes(
                enzyme, CB, ReactionInfo.Elements, EnzymeInfo.Elements[enzyme]['Reaction'])

        for rx in ReactionInfo.Elements.keys():
            if ReactionInfo.Elements[rx]['Enzyme'] is not '':
                ReactionInfo.Elements[rx]['Compartment_Machinery'] = EnzymeInfo.Elements[ReactionInfo.Elements[rx]
                                                                                         ['Enzyme']]['EnzymeCompartment']

        for comp in CompartmentInfo.Elements.keys():
            ContainedProteins = findContainedProteins(comp, ProteinInfo.Elements)
            containedEnzymes = findContainedEnzymes(
                ContainedProteins, ProteinInfo.Elements, EnzymeInfo.Elements)
            CompartmentInfo.Elements[comp]['associatedProteins'] = ContainedProteins
            CompartmentInfo.Elements[comp]['associatedReactions'] = containedEnzymes['containedReactions']
            CompartmentInfo.Elements[comp]['associatedEnzymes'] = containedEnzymes['containedEnzymes']

        self.ReactionInfo = ReactionInfo
        self.EnzymeInfo = EnzymeInfo
        self.ProteinInfo = ProteinInfo
        self.CompartmentInfo = CompartmentInfo

        self.MetaboliteConstraintsInfo = MetaboliteConstraints()
        self.DensityConstraintsInfo = DensityConstraints()
        self.ProcessConstraintsInfo = ProcessConstraints()
        self.EnzymeConstraintsInfo = EnzymeConstraints()

        self.MetaboliteConstraintsInfo.fromFiles(constraints, Zero_matrix)
        self.DensityConstraintsInfo.fromFiles(model, constraints, Zero_matrix)
        self.ProcessConstraintsInfo.fromFiles(model, constraints, Zero_matrix)
        self.EnzymeConstraintsInfo.fromFiles(model, constraints, Zero_matrix)

        BioConstraintStats = StatsConstraintsBiological(self.MetaboliteConstraintsInfo.Elements,
                                                        self.EnzymeConstraintsInfo.Elements,
                                                        self.DensityConstraintsInfo.Elements,
                                                        self.ProcessConstraintsInfo.Elements)

        MathConstraintStats = StatsConstraintsMathematical(self.MetaboliteConstraintsInfo.Elements,
                                                           self.EnzymeConstraintsInfo.Elements,
                                                           self.DensityConstraintsInfo.Elements,
                                                           self.ProcessConstraintsInfo.Elements,
                                                           self.EnzymeInfo.Elements,
                                                           self.ReactionInfo.Elements,
                                                           self.ProcessInfo.Elements)

        FullOverview = generateOverview(self.ReactionInfo.overview(),
                                        self.MetaboliteInfo.overview(),
                                        self.ModuleInfo.overview(),
                                        self.EnzymeInfo.overview(),
                                        self.ProteinInfo.overview(),
                                        self.ProcessInfo.overview(),
                                        self.CompartmentInfo.overview(),
                                        self.TargetInfo.overview(),
                                        BioConstraintStats,
                                        MathConstraintStats)

        self.ModelStatistics = StatisticsBlock()
        self.ModelStatistics.derive(FullOverview)

        self.ProteinMatrix = generateProteinMatrix(self)
        self.ProteinGeneMatrix = generateProteinGeneMatrix(self)
        self.MediumDependencies, self.MuDependencies = findParameterDependencies(self)

    def fromJSON(self, inputString):
        Block = json.loads(inputString)
        self.ModelStatistics = StatisticsBlock()
        self.GeneralInfo = DescriptionBlock()
        self.ProcessInfo = ProcessBlock()
        self.CompartmentInfo = CompartmentBlock()
        self.MetaboliteInfo = MetaboliteBlock()
        self.TargetInfo = TargetBlock()
        self.ModuleInfo = ModuleBlock()
        self.EnzymeInfo = EnzymeBlock()
        self.ProteinInfo = ProteinBlock()
        self.ReactionInfo = ReactionBlock()
        self.DensityConstraintsInfo = DensityConstraints()
        self.ProcessConstraintsInfo = ProcessConstraints()
        self.MetaboliteConstraintsInfo = MetaboliteConstraints()
        self.EnzymeConstraintsInfo = EnzymeConstraints()

        self.ModelStatistics.fromDict(Block['ModelStatistics'])
        self.GeneralInfo.fromDict(Block['ModelInformation'])
        self.ProcessInfo.fromDict(Block['Process'])
        self.CompartmentInfo.fromDict(Block['Compartment'])
        self.MetaboliteInfo.fromDict(Block['Metabolite'])
        self.ModuleInfo.fromDict(Block['Module'])
        self.EnzymeInfo.fromDict(Block['Enzyme'])
        self.ProteinInfo.fromDict(Block['Protein'])
        self.ReactionInfo.fromDict(Block['Reaction'])
        self.DensityConstraintsInfo.fromDict(Block['DensityConstraint'])
        self.ProcessConstraintsInfo.fromDict(Block['ProcessConstraint'])
        self.MetaboliteConstraintsInfo.fromDict(Block['MetaboliteConstraint'])
        self.EnzymeConstraintsInfo.fromDict(Block['EnzymeConstraint'])
        self.TargetInfo.fromDict(Block['Target'])
        self.ProteinMatrix = Block['ProteinMatrix']
        self.ProteinMatrix['Matrix'] = numpy.array(self.ProteinMatrix['Matrix'])
        self.ProteinGeneMatrix = Block['ProteinGeneMatrix']
        self.ProteinGeneMatrix['Matrix'] = numpy.array(self.ProteinGeneMatrix['Matrix'])
        self.MediumDependencies = Block['MediumDependencies']
        self.MuDependencies = Block['MuDependencies']

    def fromSBtab(self, file_name):
        from sbtab import SBtab
        sbtab_file = open(file_name, 'r')
        file_content = sbtab_file.read()
        sbtab_file.close()
        SBtabDoc = SBtab.SBtabDocument('RBAstructure', file_content, file_name)
        for table in SBtabDoc.sbtabs:
            df = table.to_data_frame()
            if table.table_name == 'Model_Information':
                self.GeneralInfo = DescriptionBlock()
                self.GeneralInfo.fromDict(df.to_dict('index'))
            if table.table_name == 'Metabolites':
                self.MetaboliteInfo = MetaboliteBlock()
                self.MetaboliteInfo.fromDict(df.to_dict('index'))
            if table.table_name == 'Reactions':
                self.ReactionInfo = ReactionBlock()
                self.ReactionInfo.fromDict(df.to_dict('index'))
            if table.table_name == 'Enzymes':
                self.EnzymeInfo = EnzymeBlock()
                self.EnzymeInfo.fromDict(df.to_dict('index'))
            if table.table_name == 'Proteins':
                self.ProteinInfo = ProteinBlock()
                self.ProteinInfo.fromDict(df.to_dict('index'))
            if table.table_name == 'Cell_Compartments':
                self.CompartmentInfo = CompartmentBlock()
                self.CompartmentInfo.fromDict(df.to_dict('index'))
            if table.table_name == 'Modules':
                self.ModuleInfo = ModuleBlock()
                self.ModuleInfo.fromDict(df.to_dict('index'))
            if table.table_name == 'Processes':
                self.ProcessInfo = ProcessBlock()
                self.ProcessInfo.fromDict(df.to_dict('index'))
            if table.table_name == 'Density_Constraints':
                self.DensityConstraintsInfo = DensityConstraints()
                self.DensityConstraintsInfo.fromDict(df.to_dict('index'))
            if table.table_name == 'Process_Capacity_Constraints':
                self.ProcessConstraintsInfo = ProcessConstraints()
                self.ProcessConstraintsInfo.fromDict(df.to_dict('index'))
            if table.table_name == 'Enzyme_Capacity_Constraints':
                self.EnzymeConstraintsInfo = EnzymeConstraints()
                self.EnzymeConstraintsInfo.fromDict(df.to_dict('index'))
            if table.table_name == 'Metabolite_Constraints':
                self.MetaboliteConstraintsInfo = MetaboliteConstraints()
                self.MetaboliteConstraintsInfo.fromDict(df.to_dict('index'))

        BioConstraintStats = StatsConstraintsBiological(self.MetaboliteConstraintsInfo.Elements,
                                                        self.EnzymeConstraintsInfo.Elements,
                                                        self.DensityConstraintsInfo.Elements,
                                                        self.ProcessConstraintsInfo.Elements)

        MathConstraintStats = StatsConstraintsMathematical(self.MetaboliteConstraintsInfo.Elements,
                                                           self.EnzymeConstraintsInfo.Elements,
                                                           self.DensityConstraintsInfo.Elements,
                                                           self.ProcessConstraintsInfo.Elements,
                                                           self.EnzymeInfo.Elements,
                                                           self.ReactionInfo.Elements,
                                                           self.ProcessInfo.Elements)

        FullOverview = generateOverview(self.ReactionInfo.overview(),
                                        self.MetaboliteInfo.overview(),
                                        self.ModuleInfo.overview(),
                                        self.EnzymeInfo.overview(),
                                        self.ProteinInfo.overview(),
                                        self.ProcessInfo.overview(),
                                        self.CompartmentInfo.overview(),
                                        BioConstraintStats,
                                        MathConstraintStats)

        self.ModelStatistics = StatisticsBlock()
        self.ModelStatistics.derive(FullOverview)

    def generateMatrices(self):
        self.InfoMatrices = InfoMatrices(self)

#    def toXML(self):
#        x = htmlStyle(self)
#        root = ET.fromstring(jxmlease.emit_xml(x, encoding='utf-8'))
#        m = ET.tostring(root, 'utf-8')
#        return(m)

    def exportJSON(self, path):
        Block = {'ModelInformation': self.GeneralInfo.Elements,
                 'ModelStatistics': self.ModelStatistics.Elements,
                 'Process': self.ProcessInfo.Elements,
                 'Compartment': self.CompartmentInfo.Elements,
                 'Metabolite': self.MetaboliteInfo.Elements,
                 'Target': self.TargetInfo.Elements,
                 'Module': self.ModuleInfo.Elements,
                 'Enzyme': self.EnzymeInfo.Elements,
                 'Protein': self.ProteinInfo.Elements,
                 'Reaction': self.ReactionInfo.Elements,
                 'DensityConstraint': self.DensityConstraintsInfo.Elements,
                 'ProcessConstraint': self.ProcessConstraintsInfo.Elements,
                 'MetaboliteConstraint': self.MetaboliteConstraintsInfo.Elements,
                 'EnzymeConstraint': self.EnzymeConstraintsInfo.Elements,
                 'ProteinMatrix': self.ProteinMatrix,
                 'ProteinGeneMatrix': self.ProteinGeneMatrix,
                 'MediumDependencies': self.MediumDependencies,
                 'MuDependencies': self.MuDependencies}
        Block['ProteinMatrix']['Matrix'] = Block['ProteinMatrix']['Matrix'].tolist()
        Block['ProteinGeneMatrix']['Matrix'] = Block['ProteinGeneMatrix']['Matrix'].tolist()
        JSONstring = json.dumps(Block, default=JSON_Int64_compensation)
        filename = path + '/ModelStructure.json'
        f = open(filename, 'w')
        f.write(JSONstring)
        f.close()
        return(JSONstring)

    def exportSBtab(self):

        self.ModelStatistics.toSBtab('RBA model', 'ModelStatistics',
                                     'Key_Numbers', 'RBAstructure', 'Statistics on RBAmodel')
        self.GeneralInfo.toSBtab('RBA model', 'ModelInformation',
                                 'Model_Information', 'RBAstructure', 'Information on RBAmodel')
        self.MetaboliteInfo.toSBtab('RBA model', 'Metabolite', 'Metabolites',
                                    'RBAstructure', 'Information on Metabolites in RBAmodel')
        self.ReactionInfo.toSBtab('RBA model', 'Reaction', 'Reactions',
                                  'RBAstructure', 'Information on Reactions in RBAmodel')
        self.EnzymeInfo.toSBtab('RBA model', 'Enzyme', 'Enzymes',
                                'RBAstructure', 'Information on Enzymes in RBAmodel')
        self.ProteinInfo.toSBtab('RBA model', 'Protein', 'Proteins',
                                 'RBAstructure', 'Information on Proteins in RBAmodel')
        self.CompartmentInfo.toSBtab('RBA model', 'Compartment', 'Cell_Compartments',
                                     'RBAstructure', 'Information on Compartments in RBAmodel')
        self.ModuleInfo.toSBtab('RBA model', 'Module', 'Modules',
                                'RBAstructure', 'Information on Modules in RBAmodel')
        self.ProcessInfo.toSBtab('RBA model', 'Process', 'Processes',
                                 'RBAstructure', 'Information on Processes in RBAmodel')
        self.MetaboliteConstraintsInfo.toSBtab(
            'RBA model', 'MetaboliteConstraint', 'Metabolite Constraints', 'RBAstructure', 'Information on MetaboliteConstraints in RBAmodel')
        self.DensityConstraintsInfo.toSBtab('RBA model', 'DensityConstraint', 'Density Constraints',
                                            'RBAstructure', 'Information on DensityConstraints in RBAmodel')
        self.ProcessConstraintsInfo.toSBtab('RBA model', 'ProcessConstraint', 'Process Capacity Constraints',
                                            'RBAstructure', 'Information on ProcessConstraints in RBAmodel')
        self.EnzymeConstraintsInfo.toSBtab('RBA model', 'EnzymeConstraint', 'Enzyme Capacity Constraints',
                                           'RBAstructure', 'Information on EnzymeConstraints in RBAmodel')

    def exportSBtab_OneFile(self):
        from sbtab import SBtab
        GeneralInfoTable = self.GeneralInfo.toSBtab_forDoc(
            'RBA model', 'Quantity', 'Model Information', 'RBAstructure', 'Information on RBAmodel', ['Property', 'Value'])
        GeneralInfoTable.filename = 'ModelInformation.tsv'
        GeneralInfoTable.change_attribute('Text', 'Information on RBAmodel')
        GeneralInfoTable.unset_attribute('Unit')
        GeneralInfoTable.change_attribute('TableID', 'ModelInformation')
        GeneralInfoTable.unset_attribute('DocumentName')
        GeneralInfoTable.unset_attribute('Document')
        GeneralInfoTable.unset_attribute('Date')
        GeneralInfoTable.unset_attribute('SBtabVersion')
        StatsTable = self.ModelStatistics.toSBtab_forDoc(
            'RBA model', 'Quantity', 'Key Numbers', 'RBAstructure', 'Statistics on RBAmodel', ['Item', 'Count_number'])
        StatsTable.filename = 'ModelStatistics.tsv'
        StatsTable.change_attribute('Text', 'Statistics on RBAmodel')
        StatsTable.unset_attribute('Unit')
        StatsTable.change_attribute('TableID', 'KeyNumber')
        StatsTable.unset_attribute('DocumentName')
        StatsTable.unset_attribute('Document')
        StatsTable.unset_attribute('Date')
        StatsTable.unset_attribute('SBtabVersion')
        if len(list(self.MetaboliteInfo.Elements.keys())) > 0:
            MetaboliteTable = self.MetaboliteInfo.toSBtab_forDoc('RBA_model', 'Quantity', 'Metabolites', 'RBAstructure', 'Information on Metabolites in RBAmodel', [
                                                                 'ID', 'Name', 'Annotations', 'Reactions', 'Boundary', 'Role', 'Compartment'])
            MetaboliteTable.filename = 'Metabolite.tsv'
            MetaboliteTable.change_attribute('Text', 'Information on Metabolites in RBAmodel')
            MetaboliteTable.unset_attribute('Unit')
            MetaboliteTable.change_attribute('TableID', 'Metabolite')
            MetaboliteTable.unset_attribute('DocumentName')
            MetaboliteTable.unset_attribute('Document')
            MetaboliteTable.unset_attribute('Date')
            MetaboliteTable.unset_attribute('SBtabVersion')

        if len(list(self.ReactionInfo.Elements.keys())) > 0:
            ReactionTable = self.ReactionInfo.toSBtab_forDoc('RBA_model', 'Quantity', 'Reactions', 'RBAstructure', 'Information on Reactions in RBAmodel', [
                                                             'ID', 'Enzyme_compartment', 'Name', 'Annotations', 'Formula', 'Reactants', 'Products', 'Reversible', 'Type', 'Metabolite_compartment', 'Enzyme', 'Isoenzyme_reactions'])
            ReactionTable.filename = 'Reaction.tsv'
            ReactionTable.change_attribute('Text', 'Information on Reactions in RBAmodel')
            ReactionTable.unset_attribute('Unit')
            ReactionTable.change_attribute('TableID', 'Reaction')
            ReactionTable.unset_attribute('DocumentName')
            ReactionTable.unset_attribute('Document')
            ReactionTable.unset_attribute('Date')
            ReactionTable.unset_attribute('SBtabVersion')

        if len(list(self.EnzymeInfo.Elements.keys())) > 0:
            EnzymeTable = self.EnzymeInfo.toSBtab_forDoc('RBA_model', 'Quantity', 'Enzymes', 'RBAstructure', 'Information on Enzymes in RBAmodel', [
                                                         'ID', 'Annotations', 'Reaction', 'Isoenzymes', 'Other_enzymatic_activities', 'Subunit_composition', 'Compartment'])
            EnzymeTable.filename = 'Enzyme.tsv'
            EnzymeTable.change_attribute('Text', 'Information on Enzymes in RBAmodel')
            EnzymeTable.unset_attribute('Unit')
            EnzymeTable.change_attribute('TableID', 'Enzyme')
            EnzymeTable.unset_attribute('DocumentName')
            EnzymeTable.unset_attribute('Document')
            EnzymeTable.unset_attribute('Date')
            EnzymeTable.unset_attribute('SBtabVersion')
        if len(list(self.ProteinInfo.Elements.keys())) > 0:
            ProteinTable = self.ProteinInfo.toSBtab_forDoc('RBA_model', 'Quantity', 'Proteins', 'RBAstructure', 'Information on Proteins in RBAmodel', [
                                                           'ID', 'Annotations', 'Function', 'Name', 'Contributes_to_reaction', 'Contained_in_enzyme', 'Compartment', 'Composition', 'Length', 'Weight', 'Production_and_processing', 'Contributes_to_process'])
            ProteinTable.filename = 'Protein.tsv'
            ProteinTable.change_attribute('Text', 'Information on Proteins in RBAmodel')
            ProteinTable.unset_attribute('Unit')
            ProteinTable.change_attribute('TableID', 'Protein')
            ProteinTable.unset_attribute('DocumentName')
            ProteinTable.unset_attribute('Document')
            ProteinTable.unset_attribute('Date')
            ProteinTable.unset_attribute('SBtabVersion')

        if len(list(self.CompartmentInfo.Elements.keys())) > 0:
            CompartmentTable = self.CompartmentInfo.toSBtab_forDoc(
                'RBA_model', 'Quantity', 'Cell Compartments', 'RBAstructure', 'Information on Compartments in RBAmodel', ['ID', 'Proteins', 'Reactions', 'Enzymes'])
            CompartmentTable.filename = 'Compartment.tsv'
            CompartmentTable.change_attribute('Text', 'Information on Compartments in RBAmodel')
            CompartmentTable.unset_attribute('Unit')
            CompartmentTable.change_attribute('TableID', 'CellCompartment')
            CompartmentTable.unset_attribute('DocumentName')
            CompartmentTable.unset_attribute('Document')
            CompartmentTable.unset_attribute('Date')
            CompartmentTable.unset_attribute('SBtabVersion')

        if len(list(self.ModuleInfo.Elements.keys())) > 0:
            ModuleTable = self.ModuleInfo.toSBtab_forDoc(
                'RBA_model', 'Quantity', 'Cell Modules', 'RBAstructure', 'Information on Modules in RBAmodel', ['ID', 'Name'])
            ModuleTable.filename = 'Module.tsv'
            ModuleTable.change_attribute('Text', 'Information on Modules in RBAmodel')
            ModuleTable.unset_attribute('Unit')
            ModuleTable.change_attribute('TableID', 'CellModule')
            ModuleTable.unset_attribute('DocumentName')
            ModuleTable.unset_attribute('Document')
            ModuleTable.unset_attribute('Date')
            ModuleTable.unset_attribute('SBtabVersion')

        if len(list(self.ProcessInfo.Elements.keys())) > 0:
            ProcessTable = self.ProcessInfo.toSBtab_forDoc('RBA_model', 'Quantity', 'Processes', 'RBAstructure', 'Information on Processes in RBAmodel', [
                                                           'ID', 'Initiation', 'Machinery_subunit_composition'])
            ProcessTable.filename = 'Process.tsv'
            ProcessTable.change_attribute('Text', 'Information on Processes in RBAmodel')
            ProcessTable.unset_attribute('Unit')
            ProcessTable.change_attribute('TableID', 'Process')
            ProcessTable.unset_attribute('DocumentName')
            ProcessTable.unset_attribute('Document')
            ProcessTable.unset_attribute('Date')
            ProcessTable.unset_attribute('SBtabVersion')

        if len(list(self.MetaboliteConstraintsInfo.Elements.keys())) > 0:
            MetaboliteConstraintTable = self.MetaboliteConstraintsInfo.toSBtab_forDoc(
                'RBA_model', 'Quantity', 'Metabolite Constraints', 'RBAstructure', 'Information on MetaboliteConstraints in RBAmodel', ['ID', 'Metabolite', 'Mathematical_type'])
            MetaboliteConstraintTable.filename = 'MetaboliteConstraint.tsv'
            MetaboliteConstraintTable.change_attribute(
                'Text', 'Information on MetaboliteConstraints in RBAmodel')
            MetaboliteConstraintTable.unset_attribute('Unit')
            MetaboliteConstraintTable.change_attribute('TableID', 'MetaboliteConstraint')
            MetaboliteConstraintTable.unset_attribute('DocumentName')
            MetaboliteConstraintTable.unset_attribute('Document')
            MetaboliteConstraintTable.unset_attribute('Date')
            MetaboliteConstraintTable.unset_attribute('SBtabVersion')

        if len(list(self.DensityConstraintsInfo.Elements.keys())) > 0:
            DensityConstraintTable = self.DensityConstraintsInfo.toSBtab_forDoc(
                'RBA_model', 'Quantity', 'Density Constraints', 'RBAstructure', 'Information on DensityConstraints in RBAmodel', ['ID', 'Compartment', 'Mathematical_type'])
            DensityConstraintTable.filename = 'DensityConstraint.tsv'
            DensityConstraintTable.change_attribute(
                'Text', 'Information on DensityConstraints in RBAmodel')
            DensityConstraintTable.unset_attribute('Unit')
            DensityConstraintTable.change_attribute('TableID', 'DensityConstraint')
            DensityConstraintTable.unset_attribute('DocumentName')
            DensityConstraintTable.unset_attribute('Document')
            DensityConstraintTable.unset_attribute('Date')
            DensityConstraintTable.unset_attribute('SBtabVersion')

        if len(list(self.ProcessConstraintsInfo.Elements.keys())) > 0:
            ProcessConstraintTable = self.ProcessConstraintsInfo.toSBtab_forDoc(
                'RBA_model', 'Quantity', 'Process Capacity Constraints', 'RBAstructure', 'Information on ProcessConstraints in RBAmodel', ['ID', 'Process', 'Mathematical_type'])
            ProcessConstraintTable.filename = 'ProcessConstraint.tsv'
            ProcessConstraintTable.change_attribute(
                'Text', 'Information on ProcessCapacityConstraints in RBAmodel')
            ProcessConstraintTable.unset_attribute('Unit')
            ProcessConstraintTable.change_attribute('TableID', 'ProcessCapacityConstraint')
            ProcessConstraintTable.unset_attribute('DocumentName')
            ProcessConstraintTable.unset_attribute('Document')
            ProcessConstraintTable.unset_attribute('Date')
            ProcessConstraintTable.unset_attribute('SBtabVersion')

        if len(list(self.EnzymeConstraintsInfo.Elements.keys())) > 0:
            EnzymeConstraintTable = self.EnzymeConstraintsInfo.toSBtab_forDoc('RBA_model', 'Quantity', 'Enzyme Capacity Constraints', 'RBAstructure', 'Information on EnzymeConstraints in RBAmodel', [
                                                                              'ID', 'Enzyme', 'Reaction', 'Direction', 'Mathematical_type'])
            EnzymeConstraintTable.filename = 'EnzymeConstraint.tsv'
            EnzymeConstraintTable.change_attribute(
                'Text', 'Information on EnzymeCapacityConstraints in RBAmodel')
            EnzymeConstraintTable.unset_attribute('Unit')
            EnzymeConstraintTable.change_attribute('TableID', 'EnzymeCapacityConstraint')
            EnzymeConstraintTable.unset_attribute('DocumentName')
            EnzymeConstraintTable.unset_attribute('Document')
            EnzymeConstraintTable.unset_attribute('Date')
            EnzymeConstraintTable.unset_attribute('SBtabVersion')

        if len(list(self.TargetInfo.Elements.keys())) > 0:
            TargetTable = self.TargetInfo.toSBtab_forDoc('RBA_model', 'Quantity', 'Cell Targets', 'RBAstructure', 'Information on Targets in RBAmodel', [
                                                         'ID', 'Group', 'Type', 'TargetSpecies', 'TargetValue'])
            TargetTable.filename = 'Target.tsv'
            TargetTable.change_attribute('Text', 'Information on CellTargets in RBAmodel')
            TargetTable.unset_attribute('Unit')
            TargetTable.change_attribute('TableID', 'CellTarget')
            TargetTable.unset_attribute('DocumentName')
            TargetTable.unset_attribute('Document')
            TargetTable.unset_attribute('Date')
            TargetTable.unset_attribute('SBtabVersion')

        Out = SBtab.SBtabDocument('RBAstructure', GeneralInfoTable, 'RBA_model.tsv')
        Out.add_sbtab(StatsTable)
        if len(list(self.MetaboliteInfo.Elements.keys())) > 0:
            Out.add_sbtab(MetaboliteTable)
        if len(list(self.ReactionInfo.Elements.keys())) > 0:
            Out.add_sbtab(ReactionTable)
        if len(list(self.EnzymeInfo.Elements.keys())) > 0:
            Out.add_sbtab(EnzymeTable)
        if len(list(self.ProteinInfo.Elements.keys())) > 0:
            Out.add_sbtab(ProteinTable)
        if len(list(self.CompartmentInfo.Elements.keys())) > 0:
            Out.add_sbtab(CompartmentTable)
        if len(list(self.ModuleInfo.Elements.keys())) > 0:
            Out.add_sbtab(ModuleTable)
        if len(list(self.ProcessInfo.Elements.keys())) > 0:
            Out.add_sbtab(ProcessTable)
        if len(list(self.TargetInfo.Elements.keys())):
            Out.add_sbtab(TargetTable)
        if len(list(self.MetaboliteConstraintsInfo.Elements.keys())):
            Out.add_sbtab(MetaboliteConstraintTable)
        if len(list(self.DensityConstraintsInfo.Elements.keys())):
            Out.add_sbtab(DensityConstraintTable)
        if len(list(self.ProcessConstraintsInfo.Elements.keys())):
            Out.add_sbtab(ProcessConstraintTable)
        if len(list(self.EnzymeConstraintsInfo.Elements.keys())):
            Out.add_sbtab(EnzymeConstraintTable)
        Out.name = 'RBA_model'
        Out.change_attribute('DocumentName', 'RBA model')
        Out.write()

    def exportSBtab_withLinks(self):
        from sbtab import SBtab
        GeneralInfoTable = self.GeneralInfo.toSBtab_forDoc(
            'RBA model', 'Quantity', 'Model Information', 'RBAstructure', 'Information on RBAmodel', ['Property', 'Value'])
        GeneralInfoTable.filename = 'ModelInformation.tsv'
        GeneralInfoTable.change_attribute('Text', 'Information on RBAmodel')
        GeneralInfoTable.unset_attribute('Unit')
        GeneralInfoTable.change_attribute('TableID', 'ModelInformation')
        GeneralInfoTable.unset_attribute('DocumentName')
        GeneralInfoTable.unset_attribute('Document')
        GeneralInfoTable.unset_attribute('Date')
        GeneralInfoTable.unset_attribute('SBtabVersion')
        StatsTable = self.ModelStatistics.toSBtab_forDoc(
            'RBA model', 'Quantity', 'Key Numbers', 'RBAstructure', 'Statistics on RBAmodel', ['Item', 'Count_number'])
        StatsTable.filename = 'ModelStatistics.tsv'
        StatsTable.change_attribute('Text', 'Statistics on RBAmodel')
        StatsTable.unset_attribute('Unit')
        StatsTable.change_attribute('TableID', 'KeyNumber')
        StatsTable.unset_attribute('DocumentName')
        StatsTable.unset_attribute('Document')
        StatsTable.unset_attribute('Date')
        StatsTable.unset_attribute('SBtabVersion')
        if len(list(self.MetaboliteInfo.Elements.keys())) > 0:
            MetaboliteBlock_forChanges = copy.deepcopy(self.MetaboliteInfo)
            for k in list(MetaboliteBlock_forChanges.Elements.keys()):
                for index, id in enumerate(MetaboliteBlock_forChanges.Elements[k]['ReactionsInvolvedWith']):
                    MetaboliteBlock_forChanges.Elements[k]['ReactionsInvolvedWith'][index] = '(!' + \
                        'Reaction'+'/'+id+'!)'
                for id in list(MetaboliteBlock_forChanges.Elements[k]['OtherIDs'].keys()):
                    MetaboliteBlock_forChanges.Elements[k]['OtherIDs'][id] = str('(!identifiers:'+id+'/'+str(
                        MetaboliteBlock_forChanges.Elements[k]['OtherIDs'][id])+'|'+str(MetaboliteBlock_forChanges.Elements[k]['OtherIDs'][id])+'!)')

            MetaboliteTable = MetaboliteBlock_forChanges.toSBtab_forDoc('RBA_model', 'Quantity', 'Metabolites', 'RBAstructure', 'Information on Metabolites in RBAmodel', [
                'ID', 'Name', 'Annotations', 'Reactions', 'Boundary', 'Role', 'Compartment'])
            MetaboliteTable.filename = 'Metabolite.tsv'
            MetaboliteTable.change_attribute('Text', 'Information on Metabolites in RBAmodel')
            MetaboliteTable.unset_attribute('Unit')
            MetaboliteTable.change_attribute('TableID', 'Metabolite')
            MetaboliteTable.unset_attribute('DocumentName')
            MetaboliteTable.unset_attribute('Document')
            MetaboliteTable.unset_attribute('Date')
            MetaboliteTable.unset_attribute('SBtabVersion')

        if len(list(self.ReactionInfo.Elements.keys())) > 0:
            ReactionBlock_forChanges = copy.deepcopy(self.ReactionInfo)
            for k in list(ReactionBlock_forChanges.Elements.keys()):
                oldEnz = ReactionBlock_forChanges.Elements[k]['Enzyme']
                ReactionBlock_forChanges.Elements[k]['Enzyme'] = '(!'+'Enzyme'+'/'+oldEnz+'!)'
                for id in list(ReactionBlock_forChanges.Elements[k]['Reactants'].keys()):
                    ReactionBlock_forChanges.Elements[k]['Reactants']['(!'+'Metabolite' +
                                                                      '/'+id+'!)'] = ReactionBlock_forChanges.Elements[k]['Reactants'].pop(id)
                for id in list(ReactionBlock_forChanges.Elements[k]['Products'].keys()):
                    ReactionBlock_forChanges.Elements[k]['Products']['(!'+'Metabolite' +
                                                                     '/'+id+'!)'] = ReactionBlock_forChanges.Elements[k]['Products'].pop(id)
                for index, id in enumerate(ReactionBlock_forChanges.Elements[k]['Twins']):
                    ReactionBlock_forChanges.Elements[k]['Twins'][index] = '(!' + \
                        'Reaction'+'/'+id+'!)'
                if 'ProtoID' in list(ReactionBlock_forChanges.Elements[k]['OtherIDs'].keys()):
                    ReactionBlock_forChanges.Elements[k]['OtherIDs'].pop('ProtoID')
                for id in list(ReactionBlock_forChanges.Elements[k]['OtherIDs'].keys()):
                    ReactionBlock_forChanges.Elements[k]['OtherIDs'][id] = str('(!identifiers:'+id+'/'+str(
                        ReactionBlock_forChanges.Elements[k]['OtherIDs'][id])+'|'+str(ReactionBlock_forChanges.Elements[k]['OtherIDs'][id])+'!)')

            ReactionTable = ReactionBlock_forChanges.toSBtab_forDoc('RBA_model', 'Quantity', 'Reactions', 'RBAstructure', 'Information on Reactions in RBAmodel', [
                'ID', 'Enzyme_compartment', 'Name', 'Annotations', 'Formula', 'Reactants', 'Products', 'Reversible', 'Type', 'Metabolite_compartment', 'Enzyme', 'Isoenzyme_reactions'])
            ReactionTable.filename = 'Reaction.tsv'
            ReactionTable.change_attribute('Text', 'Information on Reactions in RBAmodel')
            ReactionTable.unset_attribute('Unit')
            ReactionTable.change_attribute('TableID', 'Reaction')
            ReactionTable.unset_attribute('DocumentName')
            ReactionTable.unset_attribute('Document')
            ReactionTable.unset_attribute('Date')
            ReactionTable.unset_attribute('SBtabVersion')

        if len(list(self.EnzymeInfo.Elements.keys())) > 0:
            EnzymeBlock_forChanges = copy.deepcopy(self.EnzymeInfo)
            for k in list(EnzymeBlock_forChanges.Elements.keys()):
                oldRx = EnzymeBlock_forChanges.Elements[k]['Reaction']
                EnzymeBlock_forChanges.Elements[k]['Reaction'] = '(!'+'Reaction'+'/'+oldRx+'!)'
                for index, id in enumerate(EnzymeBlock_forChanges.Elements[k]['Isozymes']):
                    EnzymeBlock_forChanges.Elements[k]['Isozymes'][index] = '(!' + \
                        'Enzyme'+'/'+id+'!)'
                for index, id in enumerate(EnzymeBlock_forChanges.Elements[k]['IdenticalEnzymes']):
                    EnzymeBlock_forChanges.Elements[k]['IdenticalEnzymes'][index] = '(!' + \
                        'Enzyme'+'/'+id+'!)'
                for id in list(EnzymeBlock_forChanges.Elements[k]['Subunits'].keys()):
                    EnzymeBlock_forChanges.Elements[k]['Subunits']['(!'+'Protein'+'/' +
                                                                   id+'!)'] = EnzymeBlock_forChanges.Elements[k]['Subunits'].pop(id)

            EnzymeTable = EnzymeBlock_forChanges.toSBtab_forDoc('RBA_model', 'Quantity', 'Enzymes', 'RBAstructure', 'Information on Enzymes in RBAmodel', [
                'ID', 'Annotations', 'Reaction', 'Isoenzymes', 'Other_enzymatic_activities', 'Subunit_composition', 'Compartment'])
            EnzymeTable.filename = 'Enzyme.tsv'
            EnzymeTable.change_attribute('Text', 'Information on Enzymes in RBAmodel')
            EnzymeTable.unset_attribute('Unit')
            EnzymeTable.change_attribute('TableID', 'Enzyme')
            EnzymeTable.unset_attribute('DocumentName')
            EnzymeTable.unset_attribute('Document')
            EnzymeTable.unset_attribute('Date')
            EnzymeTable.unset_attribute('SBtabVersion')

        if len(list(self.ProteinInfo.Elements.keys())) > 0:
            ProteinBlock_forChanges = copy.deepcopy(self.ProteinInfo)
            for k in list(ProteinBlock_forChanges.Elements.keys()):
                oldComp = ProteinBlock_forChanges.Elements[k]['Compartment']
                ProteinBlock_forChanges.Elements[k]['Compartment'] = '(!' + \
                    'Compartment'+'/'+oldComp+'!)'
                for index, id in enumerate(ProteinBlock_forChanges.Elements[k]['associatedReactions']):
                    ProteinBlock_forChanges.Elements[k]['associatedReactions'][index] = '(!' + \
                        'Reaction'+'/'+id+'!)'
                for index, id in enumerate(ProteinBlock_forChanges.Elements[k]['associatedEnzymes']):
                    ProteinBlock_forChanges.Elements[k]['associatedEnzymes'][index] = '(!' + \
                        'Enzyme'+'/'+id+'!)'
                for index, id in enumerate(ProteinBlock_forChanges.Elements[k]['SupportsProcess']):
                    ProteinBlock_forChanges.Elements[k]['SupportsProcess'][index] = '(!' + \
                        'Process'+'/'+id+'!)'
                for id in list(ProteinBlock_forChanges.Elements[k]['ProcessRequirements'].keys()):
                    ProteinBlock_forChanges.Elements[k]['ProcessRequirements']['(!'+'Process'+'/' +
                                                                               id+'!)'] = ProteinBlock_forChanges.Elements[k]['ProcessRequirements'].pop(id)

                for id in list(ProteinBlock_forChanges.Elements[k]['ExternalIDs'].keys()):
                    if id == 'UniprotID':
                        ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id] = str('(!identifiers:uniprot/'+str(
                            ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id])+'|'+str(ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id])+'!)')
                    if id == 'ECnumber':
                        if ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id] != 'nan':
                            ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id] = str('(!identifiers:ec-code/'+str(
                                ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id]).split('EC ')[1]+'|'+str(ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id])+'!)')

            ProteinTable = ProteinBlock_forChanges.toSBtab_forDoc('RBA_model', 'Quantity', 'Proteins', 'RBAstructure', 'Information on Proteins in RBAmodel', [
                'ID', 'Annotations', 'Function', 'Name', 'Contributes_to_reaction', 'Contained_in_enzyme', 'Compartment', 'Composition', 'Length', 'Weight', 'Production_and_processing', 'Contributes_to_process'])
            ProteinTable.filename = 'Protein.tsv'
            ProteinTable.change_attribute('Text', 'Information on Proteins in RBAmodel')
            ProteinTable.unset_attribute('Unit')
            ProteinTable.change_attribute('TableID', 'Protein')
            ProteinTable.unset_attribute('DocumentName')
            ProteinTable.unset_attribute('Document')
            ProteinTable.unset_attribute('Date')
            ProteinTable.unset_attribute('SBtabVersion')

        if len(list(self.CompartmentInfo.Elements.keys())) > 0:
            CompartmentBlock_forChanges = copy.deepcopy(self.CompartmentInfo)
            for k in list(CompartmentBlock_forChanges.Elements.keys()):
                for index, id in enumerate(CompartmentBlock_forChanges.Elements[k]['associatedReactions']):
                    CompartmentBlock_forChanges.Elements[k]['associatedReactions'][index] = '(!' + \
                        'Reaction'+'/'+id+'!)'
                for index, id in enumerate(CompartmentBlock_forChanges.Elements[k]['associatedEnzymes']):
                    CompartmentBlock_forChanges.Elements[k]['associatedEnzymes'][index] = '(!' + \
                        'Enzyme'+'/'+id+'!)'
                for index, id in enumerate(CompartmentBlock_forChanges.Elements[k]['associatedProteins']):
                    CompartmentBlock_forChanges.Elements[k]['associatedProteins'][index] = '(!' + \
                        'Protein'+'/'+id+'!)'

            CompartmentTable = CompartmentBlock_forChanges.toSBtab_forDoc(
                'RBA_model', 'Quantity', 'Cell Compartments', 'RBAstructure', 'Information on Compartments in RBAmodel', ['ID', 'Proteins', 'Reactions', 'Enzymes'])
            CompartmentTable.filename = 'Compartment.tsv'
            CompartmentTable.change_attribute('Text', 'Information on Compartments in RBAmodel')
            CompartmentTable.unset_attribute('Unit')
            CompartmentTable.change_attribute('TableID', 'CellCompartment')
            CompartmentTable.unset_attribute('DocumentName')
            CompartmentTable.unset_attribute('Document')
            CompartmentTable.unset_attribute('Date')
            CompartmentTable.unset_attribute('SBtabVersion')

        if len(list(self.ModuleInfo.Elements.keys())) > 0:
            ModuleTable = self.ModuleInfo.toSBtab_forDoc(
                'RBA_model', 'Quantity', 'Cell Modules', 'RBAstructure', 'Information on Modules in RBAmodel', ['ID', 'Name'])
            ModuleTable.filename = 'Module.tsv'
            ModuleTable.change_attribute('Text', 'Information on Modules in RBAmodel')
            ModuleTable.unset_attribute('Unit')
            ModuleTable.change_attribute('TableID', 'CellModule')
            ModuleTable.unset_attribute('DocumentName')
            ModuleTable.unset_attribute('Document')
            ModuleTable.unset_attribute('Date')
            ModuleTable.unset_attribute('SBtabVersion')

        if len(list(self.ProcessInfo.Elements.keys())) > 0:
            ProcessBlock_forChanges = copy.deepcopy(self.ProcessInfo)
            for k in list(ProcessBlock_forChanges.Elements.keys()):
                for compo in ProcessBlock_forChanges.Elements[k]['Components']:
                    for id in ProcessBlock_forChanges.Elements[k]['Components'][compo]['Reactants']:
                        ProcessBlock_forChanges.Elements[k]['Components'][compo]['Reactants']['(!'+'Metabolite' +
                                                                                              '/'+id+'!)'] = ProcessBlock_forChanges.Elements[k]['Components'][compo]['Reactants'].pop(id)
                    for id in ProcessBlock_forChanges.Elements[k]['Components'][compo]['Products']:
                        ProcessBlock_forChanges.Elements[k]['Components'][compo]['Products']['(!'+'Metabolite' +
                                                                                             '/'+id+'!)'] = ProcessBlock_forChanges.Elements[k]['Components'][compo]['Products'].pop(id)
                for id in list(ProcessBlock_forChanges.Elements[k]['Composition'].keys()):
                    ProcessBlock_forChanges.Elements[k]['Composition']['(!'+'Process'+'/' +
                                                                       id+'!)'] = ProcessBlock_forChanges.Elements[k]['Composition'].pop(id)

            ProcessTable = ProcessBlock_forChanges.toSBtab_forDoc('RBA_model', 'Quantity', 'Processes', 'RBAstructure', 'Information on Processes in RBAmodel', [
                'ID', 'Initiation', 'Machinery_subunit_composition'])
            ProcessTable.filename = 'Process.tsv'
            ProcessTable.change_attribute('Text', 'Information on Processes in RBAmodel')
            ProcessTable.unset_attribute('Unit')
            ProcessTable.change_attribute('TableID', 'Process')
            ProcessTable.unset_attribute('DocumentName')
            ProcessTable.unset_attribute('Document')
            ProcessTable.unset_attribute('Date')
            ProcessTable.unset_attribute('SBtabVersion')

        if len(list(self.MetaboliteConstraintsInfo.Elements.keys())) > 0:
            MetaboliteConstraintsBlock_forChanges = copy.deepcopy(self.MetaboliteConstraintsInfo)
            for k in list(MetaboliteConstraintsBlock_forChanges.Elements.keys()):
                oldMet = MetaboliteConstraintsBlock_forChanges.Elements[k]['AssociatedMetabolite']
                MetaboliteConstraintsBlock_forChanges.Elements[k]['AssociatedMetabolite'] = '(!' + \
                    'Metabolite'+'/'+oldMet+'!)'

            MetaboliteConstraintTable = MetaboliteConstraintsBlock_forChanges.toSBtab_forDoc(
                'RBA_model', 'Quantity', 'Metabolite Constraints', 'RBAstructure', 'Information on MetaboliteConstraints in RBAmodel', ['ID', 'Metabolite', 'Mathematical_type'])
            MetaboliteConstraintTable.filename = 'MetaboliteConstraint.tsv'
            MetaboliteConstraintTable.change_attribute(
                'Text', 'Information on MetaboliteConstraints in RBAmodel')
            MetaboliteConstraintTable.unset_attribute('Unit')
            MetaboliteConstraintTable.change_attribute('TableID', 'MetaboliteConstraint')
            MetaboliteConstraintTable.unset_attribute('DocumentName')
            MetaboliteConstraintTable.unset_attribute('Document')
            MetaboliteConstraintTable.unset_attribute('Date')
            MetaboliteConstraintTable.unset_attribute('SBtabVersion')

        if len(list(self.DensityConstraintsInfo.Elements.keys())) > 0:
            CompartmentConstraintsBlock_forChanges = copy.deepcopy(self.DensityConstraintsInfo)
            for k in list(CompartmentConstraintsBlock_forChanges.Elements.keys()):
                oldComp = CompartmentConstraintsBlock_forChanges.Elements[k]['AssociatedCompartment']
                CompartmentConstraintsBlock_forChanges.Elements[k]['AssociatedCompartment'] = '(!' + \
                    'Compartment'+'/'+oldComp+'!)'

            DensityConstraintTable = CompartmentConstraintsBlock_forChanges.toSBtab_forDoc(
                'RBA_model', 'Quantity', 'Density Constraints', 'RBAstructure', 'Information on DensityConstraints in RBAmodel', ['ID', 'Compartment', 'Mathematical_type'])
            DensityConstraintTable.filename = 'DensityConstraint.tsv'
            DensityConstraintTable.change_attribute(
                'Text', 'Information on DensityConstraints in RBAmodel')
            DensityConstraintTable.unset_attribute('Unit')
            DensityConstraintTable.change_attribute('TableID', 'DensityConstraint')
            DensityConstraintTable.unset_attribute('DocumentName')
            DensityConstraintTable.unset_attribute('Document')
            DensityConstraintTable.unset_attribute('Date')
            DensityConstraintTable.unset_attribute('SBtabVersion')

        if len(list(self.ProcessConstraintsInfo.Elements.keys())) > 0:
            ProcessConstraintsBlock_forChanges = copy.deepcopy(self.ProcessConstraintsInfo)
            for k in list(ProcessConstraintsBlock_forChanges.Elements.keys()):
                oldComp = ProcessConstraintsBlock_forChanges.Elements[k]['AssociatedProcess']
                ProcessConstraintsBlock_forChanges.Elements[k]['AssociatedProcess'] = '(!' + \
                    'Process'+'/'+oldComp+'!)'

            ProcessConstraintTable = ProcessConstraintsBlock_forChanges.toSBtab_forDoc(
                'RBA_model', 'Quantity', 'Process Capacity Constraints', 'RBAstructure', 'Information on ProcessConstraints in RBAmodel', ['ID', 'Process', 'Mathematical_type'])
            ProcessConstraintTable.filename = 'ProcessConstraint.tsv'
            ProcessConstraintTable.change_attribute(
                'Text', 'Information on ProcessCapacityConstraints in RBAmodel')
            ProcessConstraintTable.unset_attribute('Unit')
            ProcessConstraintTable.change_attribute('TableID', 'ProcessCapacityConstraint')
            ProcessConstraintTable.unset_attribute('DocumentName')
            ProcessConstraintTable.unset_attribute('Document')
            ProcessConstraintTable.unset_attribute('Date')
            ProcessConstraintTable.unset_attribute('SBtabVersion')

        if len(list(self.EnzymeConstraintsInfo.Elements.keys())) > 0:
            EnzymeConstraintsBlock_forChanges = copy.deepcopy(self.EnzymeConstraintsInfo)
            for k in list(EnzymeConstraintsBlock_forChanges.Elements.keys()):
                oldComp = EnzymeConstraintsBlock_forChanges.Elements[k]['AssociatedEnzyme']
                EnzymeConstraintsBlock_forChanges.Elements[k]['AssociatedEnzyme'] = '(!' + \
                    'Enzyme'+'/'+oldComp+'!)'

            EnzymeConstraintTable = EnzymeConstraintsBlock_forChanges.toSBtab_forDoc('RBA_model', 'Quantity', 'Enzyme Capacity Constraints', 'RBAstructure', 'Information on EnzymeConstraints in RBAmodel', [
                'ID', 'Enzyme', 'Reaction', 'Direction', 'Mathematical_type'])
            EnzymeConstraintTable.filename = 'EnzymeConstraint.tsv'
            EnzymeConstraintTable.change_attribute(
                'Text', 'Information on EnzymeCapacityConstraints in RBAmodel')
            EnzymeConstraintTable.unset_attribute('Unit')
            EnzymeConstraintTable.change_attribute('TableID', 'EnzymeCapacityConstraint')
            EnzymeConstraintTable.unset_attribute('DocumentName')
            EnzymeConstraintTable.unset_attribute('Document')
            EnzymeConstraintTable.unset_attribute('Date')
            EnzymeConstraintTable.unset_attribute('SBtabVersion')

        if len(list(self.TargetInfo.Elements.keys())) > 0:
            TargetBlock_forChanges = copy.deepcopy(self.TargetInfo)
            for k in list(TargetBlock_forChanges.Elements.keys()):
                oldTargSpec = TargetBlock_forChanges.Elements[k]['TargetSpecies']
                if oldTargSpec in list(self.MetaboliteInfo.Elements.keys()):
                    TargetBlock_forChanges.Elements[k]['TargetSpecies'] = '(!' + \
                        'Metabolite'+'/'+oldTargSpec+'!)'
                if oldTargSpec in list(self.ReactionInfo.Elements.keys()):
                    TargetBlock_forChanges.Elements[k]['TargetSpecies'] = '(!' + \
                        'Reaction'+'/'+oldTargSpec+'!)'
                if oldTargSpec in list(self.ProteinInfo.Elements.keys()):
                    TargetBlock_forChanges.Elements[k]['TargetSpecies'] = '(!' + \
                        'Protein'+'/'+oldTargSpec+'!)'

            TargetTable = TargetBlock_forChanges.toSBtab_forDoc('RBA_model', 'Quantity', 'Cell Targets', 'RBAstructure', 'Information on Targets in RBAmodel', [
                'ID', 'Group', 'Type', 'TargetSpecies', 'TargetValue'])
            TargetTable.filename = 'Target.tsv'
            TargetTable.change_attribute('Text', 'Information on CellTargets in RBAmodel')
            TargetTable.unset_attribute('Unit')
            TargetTable.change_attribute('TableID', 'CellTarget')
            TargetTable.unset_attribute('DocumentName')
            TargetTable.unset_attribute('Document')
            TargetTable.unset_attribute('Date')
            TargetTable.unset_attribute('SBtabVersion')

        Out = SBtab.SBtabDocument('RBAstructure_withLinks',
                                  GeneralInfoTable, 'RBA_model_withLinks.tsv')
        Out.add_sbtab(StatsTable)
        if len(list(self.MetaboliteInfo.Elements.keys())) > 0:
            Out.add_sbtab(MetaboliteTable)
        if len(list(self.ReactionInfo.Elements.keys())) > 0:
            Out.add_sbtab(ReactionTable)
        if len(list(self.EnzymeInfo.Elements.keys())) > 0:
            Out.add_sbtab(EnzymeTable)
        if len(list(self.ProteinInfo.Elements.keys())) > 0:
            Out.add_sbtab(ProteinTable)
        if len(list(self.CompartmentInfo.Elements.keys())) > 0:
            Out.add_sbtab(CompartmentTable)
        if len(list(self.ModuleInfo.Elements.keys())) > 0:
            Out.add_sbtab(ModuleTable)
        if len(list(self.ProcessInfo.Elements.keys())) > 0:
            Out.add_sbtab(ProcessTable)
        if len(list(self.TargetInfo.Elements.keys())):
            Out.add_sbtab(TargetTable)
        if len(list(self.MetaboliteConstraintsInfo.Elements.keys())):
            Out.add_sbtab(MetaboliteConstraintTable)
        if len(list(self.DensityConstraintsInfo.Elements.keys())):
            Out.add_sbtab(DensityConstraintTable)
        if len(list(self.ProcessConstraintsInfo.Elements.keys())):
            Out.add_sbtab(ProcessConstraintTable)
        if len(list(self.EnzymeConstraintsInfo.Elements.keys())):
            Out.add_sbtab(EnzymeConstraintTable)
        Out.name = 'RBA_model_withLinks'
        Out.change_attribute('DocumentName', 'RBA model with links')
        Out.write()


def findParameterDependencies(ModelStructure):
    MedDepts = {}
    MuDepts = []
    for eK in ModelStructure.EnzymeConstraintsInfo.Elements.keys():
        e = ModelStructure.EnzymeConstraintsInfo.Elements[eK]
        for pf in e['CapacityParameter']:
            iV = list(pf.values())[0]['IndependentVariable']
            if list(pf.values())[0]['FunctionType'] != 'constant':
                if iV.startswith('M_'):
                    if iV.rsplit('_', 1)[0] in MedDepts.keys():
                        MedDepts[iV.rsplit('_', 1)[0]].append(e['ID'])
                    elif iV.rsplit('_', 1)[0] not in MedDepts.keys():
                        MedDepts.update({iV.rsplit('_', 1)[0]: [e['ID']]})
                if iV == 'growth_rate':
                    if e['ID'] not in MuDepts:
                        MuDepts.append(e['ID'])
    for dK in ModelStructure.DensityConstraintsInfo.Elements.keys():
        d = ModelStructure.DensityConstraintsInfo.Elements[dK]
        for pf in d['CapacityParameter']:
            iV = list(pf.values())[0]['IndependentVariable']
            if list(pf.values())[0]['FunctionType'] != 'constant':
                if iV.startswith('M_'):
                    if iV.rsplit('_', 1)[0] in MedDepts.keys():
                        MedDepts[iV.rsplit('_', 1)[0]].append(d['ID'])
                    elif iV.rsplit('_', 1)[0] not in MedDepts.keys():
                        MedDepts.update({iV.rsplit('_', 1)[0]: [d['ID']]})
                if iV == 'growth_rate':
                    if d['ID'] not in MuDepts:
                        MuDepts.append(d['ID'])
    for pK in ModelStructure.ProcessConstraintsInfo.Elements.keys():
        p = ModelStructure.ProcessConstraintsInfo.Elements[pK]
        for pf in p['CapacityParameter']:
            iV = list(pf.values())[0]['IndependentVariable']
            if list(pf.values())[0]['FunctionType'] != 'constant':
                if iV.startswith('M_'):
                    if iV.rsplit('_', 1)[0] in MedDepts.keys():
                        MedDepts[iV.rsplit('_', 1)[0]].append(p['ID'])
                    elif iV.rsplit('_', 1)[0] not in MedDepts.keys():
                        MedDepts.update({iV.rsplit('_', 1)[0]: [p['ID']]})
                if iV == 'growth_rate':
                    if p['ID'] not in MuDepts:
                        MuDepts.append(p['ID'])
    return(MedDepts, MuDepts)


def generateProteinGeneMatrix(ModelStructure):
    uniqueProteins = []
    uniqueProteinMap = {}
    Proteins = list(ModelStructure.ProteinInfo.Elements.keys())
    for i in ModelStructure.ProteinInfo.Elements.keys():
        if ModelStructure.ProteinInfo.Elements[i]['ProtoID'] not in list(uniqueProteinMap.keys()):
            uniqueProteinMap.update({ModelStructure.ProteinInfo.Elements[i]['ProtoID']: []})
            uniqueProteins.append(ModelStructure.ProteinInfo.Elements[i]['ProtoID'])
        uniqueProteinMap[ModelStructure.ProteinInfo.Elements[i]['ProtoID']].append(i)
    ProteinProteinMatrix = numpy.zeros(
        (len(list(uniqueProteinMap.keys())), len(list(ModelStructure.ProteinInfo.Elements.keys()))))
    for u in list(uniqueProteinMap.keys()):
        row_ind = uniqueProteins.index(u)
        for i in uniqueProteinMap[u]:
            col_ind = Proteins.index(i)
            ProteinProteinMatrix[row_ind, col_ind] = 1
    return({'Matrix': numpy.array(ProteinProteinMatrix), 'Proteins': Proteins, 'ProtoProteins': uniqueProteins})


def generateProteinMatrix(ModelStructure):
    Proteins = list(ModelStructure.ProteinInfo.Elements.keys())
    # print(list(ModelStructure.ProcessInfo.Elements.keys()))
    Processes = [ModelStructure.ProcessInfo.Elements[i]['ID'] +
                 '_machinery' for i in list(ModelStructure.ProcessInfo.Elements.keys())]
    Enzymes = list(ModelStructure.EnzymeInfo.Elements.keys())
    Consumers = list(set(list(Enzymes+Processes)))
    ProteinMatrix = numpy.zeros((len(Proteins), len(Consumers)))
    for p in Proteins:
        if len(ModelStructure.ProteinInfo.Elements[p]['SupportsProcess']) > 0:
            # print(list(ModelStructure.ProteinInfo.Elements[p]['SupportsProcess']))
            for pc in list(ModelStructure.ProteinInfo.Elements[p]['SupportsProcess']):
                coeff = 0
                row_ind = Proteins.index(p)
                col_ind = Consumers.index(
                    ModelStructure.ProcessInfo.Elements[pc]['ID']+'_machinery')
                coeff = ModelStructure.ProcessInfo.Elements[pc]['Composition'][p]
                ProteinMatrix[row_ind, col_ind] += coeff
        if len(ModelStructure.ProteinInfo.Elements[p]['associatedEnzymes']) > 0:
            for ez in list(ModelStructure.ProteinInfo.Elements[p]['associatedEnzymes']):
                coeff = 0
                row_ind = Proteins.index(p)
                col_ind = Consumers.index(ez)
                coeff = ModelStructure.EnzymeInfo.Elements[ez]['Subunits'][p]['StochFac']
                ProteinMatrix[row_ind, col_ind] += coeff
    return({'Matrix': numpy.array(ProteinMatrix), 'Consumers': Consumers, 'Proteins': Proteins})


def importBiggMetabolites(xml_dir):
    if os.path.isfile(str(xml_dir+'/bigg_models_metabolites.txt')):
        return(pandas.read_csv(str(xml_dir+'/bigg_models_metabolites.txt'), sep='\t', index_col=0))
    else:
        sys.exit('\n Required BiGG Metabolite File "bigg_models_metabolites.txt" not found.\n' +
                 ' Please provide in input-directory\n' +
                 ' To be found under: http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt\n')


def importBiggReactions(xml_dir):
    if os.path.isfile(str(xml_dir+'/bigg_models_reactions.txt')):
        return(pandas.read_csv(str(xml_dir+'/bigg_models_reactions.txt'), sep='\t', index_col=0))
    else:
        sys.exit('\n Required BiGG Reaction File "bigg_models_reactions.txt "not found.\n' +
                 ' Please provide in input-directory\n' +
                 ' To be found under: http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt\n')


def importUniprotFile(xml_dir):
    if os.path.isfile(str(xml_dir+'/uniprot.csv')):
        return(pandas.read_csv(str(xml_dir+'/uniprot.csv'), sep='\t'))
    else:
        print('\n Uniprot-file "uniprot.csv" not found.\n' +
              ' Continuing without additional information...\n')
        return(str('Not There'))


def importSbmlFile(xml_dir, filename):
    if os.path.isfile(str(xml_dir+'/'+filename)):
        SBfile = libsbml.readSBML(str(xml_dir+'/'+filename))
        if SBfile.getNumErrors() > 0:
            SBfile.printErrors()
            print('Invalid SBML')
            return(str('Not There'))
        else:
            sbml = SBfile
            return(sbml)
    else:
        print('\n SBML-file {} not found.\n' +
              ' Continuing without additional information...\n'.format(filename))
        return(str('Not There'))


def importGeneAnnotations(xml_dir):
    if os.path.isfile(str(xml_dir+'/GeneAnnotations.csv')):
        out = pandas.read_csv(str(xml_dir+'/GeneAnnotations.csv'), sep=',', index_col=0)
        if len(list(out)) == 0:
            print('WARNING: Your file "GeneAnnotations.csv" seems to be empty or has the wrong delimiter (comma required).')
        return(out)
    else:
        print('\n No Gene-annotation file "GeneAnnotations.csv" provided.\n' +
              ' Continuing without additional information...\n')
        return(str('Not There'))


def importReactionAnnotations(xml_dir):
    if os.path.isfile(str(xml_dir+'/ReactionAnnotations.csv')):
        out = pandas.read_csv(str(xml_dir+'/ReactionAnnotations.csv'), sep=',', index_col=0)
        if len(list(out)) == 0:
            print('WARNING: Your file "ReactionAnnotations.csv" seems to be empty or has the wrong delimiter (comma required).')
        return(out)
    else:
        print('\n No Reaction-annotation file "ReactionAnnotations.csv" provided.\n' +
              ' Continuing without additional information...\n')
        return(str('Not There'))


def importMetaboliteAnnotations(xml_dir):
    if os.path.isfile(str(xml_dir+'/MetaboliteAnnotations.csv')):
        out = pandas.read_csv(str(xml_dir+'/MetaboliteAnnotations.csv'), sep=',', index_col=0)
        if len(list(out)) == 0:
            print('WARNING: Your file "MetaboliteAnnotations.csv" seems to be empty or has the wrong delimiter (comma required).')
        return(out)
    else:
        print('\n No Reaction-annotation file "MetaboliteAnnotations.csv" provided.\n' +
              ' Continuing without additional information...\n')
        return(str('Not There'))


def importModelInfo(xml_dir):
    if os.path.isfile(str(xml_dir+'/ModelInformation.csv')):
        out = pandas.read_csv(str(xml_dir+'/ModelInformation.csv'),
                              sep=',', header=0)
        out.index = list(out['Key'])
        if len(list(out)) == 0:
            print('WARNING: Your file "ModelInformation.csv" seems to be empty or has the wrong delimiter (comma required).')
        return(out)
    else:
        print('\n No model-info file "ModelInformation.csv" provided.\n' +
              ' Using dummy-information\n')
        return(pandas.DataFrame([['Name', 'ModelName'], ['Author', 'John Doe'], ['Organism', 'Life'], ['Reconstruction', 'GSMM'], ['SBML-file', 'Not Provided']], index=['Name', 'Author', 'Organism', 'Reconstruction', 'SBML-file'], columns=['Key', 'Value']))


def htmlStyle(structOriginal):
    struct = copy.deepcopy(structOriginal)

    for j in struct.ProcessInfo.Elements.keys():
        for i in struct.ProcessInfo.Elements[j]['Composition'].keys():
            struct.ProcessInfo.Elements[j]['Composition']['Protein##' +
                                                          i] = struct.ProcessInfo.Elements[j]['Composition'].pop(i)

    for j in struct.EnzymeInfo.Elements.keys():
        for i in struct.EnzymeInfo.Elements[j]['Subunits'].keys():
            struct.EnzymeInfo.Elements[j]['Subunits']['Protein##' +
                                                      i] = struct.EnzymeInfo.Elements[j]['Subunits'].pop(i)
        for i in range(len(struct.EnzymeInfo.Elements[j]['Isozymes'])):
            struct.EnzymeInfo.Elements[j]['Isozymes'][i] = 'Enzyme##' + \
                struct.EnzymeInfo.Elements[j]['Isozymes'][i]
        for i in range(len(struct.EnzymeInfo.Elements[j]['IdenticalEnzymes'])):
            struct.EnzymeInfo.Elements[j]['IdenticalEnzymes'][i] = 'Enzyme##' + \
                struct.EnzymeInfo.Elements[j]['IdenticalEnzymes'][i]
        for i in range(len(struct.EnzymeInfo.Elements[j]['EnzymeCompartment'])):
            struct.EnzymeInfo.Elements[j]['EnzymeCompartment'][i] = 'Compartment##' + \
                struct.EnzymeInfo.Elements[j]['EnzymeCompartment'][i]
        struct.EnzymeInfo.Elements[j]['Reaction'] = 'Reaction##' + \
            struct.EnzymeInfo.Elements[j]['Reaction']

    for j in struct.ReactionInfo.Elements.keys():
        for i in struct.ReactionInfo.Elements[j]['Reactants'].keys():
            struct.ReactionInfo.Elements[j]['Reactants']['Metabolite##' +
                                                         i] = struct.ReactionInfo.Elements[j]['Reactants'].pop(i)
        for i in struct.ReactionInfo.Elements[j]['Products'].keys():
            struct.ReactionInfo.Elements[j]['Products']['Metabolite##' +
                                                        i] = struct.ReactionInfo.Elements[j]['Products'].pop(i)
        for i in range(len(struct.ReactionInfo.Elements[j]['Twins'])):
            struct.ReactionInfo.Elements[j]['Twins'][i] = 'Reaction##' + \
                struct.ReactionInfo.Elements[j]['Twins'][i]
        struct.ReactionInfo.Elements[j]['Enzyme'] = 'Enzyme##' + \
            struct.ReactionInfo.Elements[j]['Enzyme']

    for j in struct.ProteinInfo.Elements.keys():
        for i in struct.ProteinInfo.Elements[j]['ProcessRequirements'].keys():
            struct.ProteinInfo.Elements[j]['ProcessRequirements']['Process##' +
                                                                  i] = struct.ProteinInfo.Elements[j]['ProcessRequirements'].pop(i)
        for i in range(len(struct.ProteinInfo.Elements[j]['associatedReactions'])):
            struct.ProteinInfo.Elements[j]['associatedReactions'][i] = 'Reaction##' + \
                struct.ProteinInfo.Elements[j]['associatedReactions'][i]
        for i in range(len(struct.ProteinInfo.Elements[j]['associatedEnzymes'])):
            struct.ProteinInfo.Elements[j]['associatedEnzymes'][i] = 'Enzyme##' + \
                struct.ProteinInfo.Elements[j]['associatedEnzymes'][i]
        for i in range(len(struct.ProteinInfo.Elements[j]['SupportsProcess'])):
            struct.ProteinInfo.Elements[j]['SupportsProcess'][i] = 'Process##' + \
                struct.ProteinInfo.Elements[j]['SupportsProcess'][i]
        struct.ProteinInfo.Elements[j]['Compartment'] = 'Compartment##' + \
            struct.ProteinInfo.Elements[j]['Compartment']

    for j in struct.CompartmentInfo.Elements.keys():
        for i in range(len(struct.CompartmentInfo.Elements[j]['associatedProteins'])):
            struct.CompartmentInfo.Elements[j]['associatedProteins'][i] = 'Protein##' + \
                struct.CompartmentInfo.Elements[j]['associatedProteins'][i]
        for i in range(len(struct.CompartmentInfo.Elements[j]['associatedReactions'])):
            struct.CompartmentInfo.Elements[j]['associatedReactions'][i] = 'Reaction##' + \
                struct.CompartmentInfo.Elements[j]['associatedReactions'][i]
        for i in range(len(struct.CompartmentInfo.Elements[j]['associatedEnzymes'])):
            struct.CompartmentInfo.Elements[j]['associatedEnzymes'][i] = 'Enzyme##' + \
                struct.CompartmentInfo.Elements[j]['associatedEnzymes'][i]

    for j in struct.MetaboliteInfo.Elements.keys():
        for i in range(len(struct.MetaboliteInfo.Elements[j]['ReactionsInvolvedWith'])):
            struct.MetaboliteInfo.Elements[j]['ReactionsInvolvedWith'][i] = 'Reaction##' + \
                struct.MetaboliteInfo.Elements[j]['ReactionsInvolvedWith'][i]

    for j in struct.MetaboliteConstraintsInfo.Elements.keys():
        struct.MetaboliteConstraintsInfo.Elements[j]['AssociatedMetabolite'] = 'Metabolite##' + \
            struct.MetaboliteConstraintsInfo.Elements[j]['AssociatedMetabolite']
    for j in struct.EnzymeConstraintsInfo.Elements.keys():
        struct.EnzymeConstraintsInfo.Elements[j]['AssociatedEnzyme'] = 'Enzyme##' + \
            struct.EnzymeConstraintsInfo.Elements[j]['AssociatedEnzyme']
    for j in struct.EnzymeConstraintsInfo.Elements.keys():
        struct.EnzymeConstraintsInfo.Elements[j]['AssociatedReaction'] = 'Reaction##' + \
            struct.EnzymeConstraintsInfo.Elements[j]['AssociatedReaction']
    for j in struct.ProcessConstraintsInfo.Elements.keys():
        struct.ProcessConstraintsInfo.Elements[j]['AssociatedProcess'] = 'Process##' + \
            struct.ProcessConstraintsInfo.Elements[j]['AssociatedProcess']
    for j in struct.DensityConstraintsInfo.Elements.keys():
        struct.DensityConstraintsInfo.Elements[j]['AssociatedCompartment'] = 'Compartment##' + \
            struct.DensityConstraintsInfo.Elements[j]['AssociatedCompartment']

    for i in struct.CompartmentInfo.Elements.keys():
        struct.CompartmentInfo.Elements['ID_' + i] = struct.CompartmentInfo.Elements.pop(i)
    for i in struct.ProcessInfo.Elements.keys():
        struct.ProcessInfo.Elements['ID_' + i] = struct.ProcessInfo.Elements.pop(i)
    for i in struct.MetaboliteInfo.Elements.keys():
        struct.MetaboliteInfo.Elements['ID_' + i] = struct.MetaboliteInfo.Elements.pop(i)
    for i in struct.ModuleInfo.Elements.keys():
        struct.ModuleInfo.Elements['ID_' + i] = struct.ModuleInfo.Elements.pop(i)
    for i in struct.EnzymeInfo.Elements.keys():
        struct.EnzymeInfo.Elements['ID_' + i] = struct.EnzymeInfo.Elements.pop(i)
    for i in struct.ProteinInfo.Elements.keys():
        struct.ProteinInfo.Elements['ID_' + i] = struct.ProteinInfo.Elements.pop(i)
    for i in struct.ReactionInfo.Elements.keys():
        struct.ReactionInfo.Elements['ID_' + i] = struct.ReactionInfo.Elements.pop(i)
    for i in struct.MetaboliteConstraintsInfo.Elements.keys():
        struct.MetaboliteConstraintsInfo.Elements['ID_' +
                                                  i] = struct.MetaboliteConstraintsInfo.Elements.pop(i)
    for i in struct.EnzymeConstraintsInfo.Elements.keys():
        struct.EnzymeConstraintsInfo.Elements['ID_' +
                                              i] = struct.EnzymeConstraintsInfo.Elements.pop(i)
    for i in struct.ProcessConstraintsInfo.Elements.keys():
        struct.ProcessConstraintsInfo.Elements['ID_' +
                                               i] = struct.ProcessConstraintsInfo.Elements.pop(i)
    for i in struct.DensityConstraintsInfo.Elements.keys():
        struct.DensityConstraintsInfo.Elements['ID_' +
                                               i] = struct.DensityConstraintsInfo.Elements.pop(i)

    struct.ProcessInfo.Elements.update(
        {'Description': 'abcdefghijklmnopqrstuvwxyz', 'Pictures': []})
    struct.CompartmentInfo.Elements.update(
        {'Description': 'abcdefghijklmnopqrstuvwxyz', 'Pictures': []})
    struct.MetaboliteInfo.Elements.update(
        {'Description': 'abcdefghijklmnopqrstuvwxyz', 'Pictures': []})
    struct.ModuleInfo.Elements.update({'Description': 'abcdefghijklmnopqrstuvwxyz', 'Pictures': []})
    struct.EnzymeInfo.Elements.update({'Description': 'abcdefghijklmnopqrstuvwxyz', 'Pictures': []})
    struct.ProteinInfo.Elements.update(
        {'Description': 'abcdefghijklmnopqrstuvwxyz', 'Pictures': []})
    struct.ReactionInfo.Elements.update(
        {'Description': 'abcdefghijklmnopqrstuvwxyz', 'Pictures': []})
    struct.DensityConstraintsInfo.Elements.update(
        {'Description': 'abcdefghijklmnopqrstuvwxyz', 'Pictures': []})
    struct.ProcessConstraintsInfo.Elements.update(
        {'Description': 'abcdefghijklmnopqrstuvwxyz', 'Pictures': []})
    struct.MetaboliteConstraintsInfo.Elements.update(
        {'Description': 'abcdefghijklmnopqrstuvwxyz', 'Pictures': []})
    struct.EnzymeConstraintsInfo.Elements.update(
        {'Description': 'abcdefghijklmnopqrstuvwxyz', 'Pictures': []})

    Block = {'ModelInformation': struct.GeneralInfo.JSONize(),
             'ModelStatistics': struct.ModelStatistics.JSONize(),
             'Process': struct.ProcessInfo.JSONize(),
             'Compartment': struct.CompartmentInfo.JSONize(),
             'Metabolite': struct.MetaboliteInfo.JSONize(),
             'Module': struct.ModuleInfo.JSONize(),
             'Enzyme': struct.EnzymeInfo.JSONize(),
             'Protein': struct.ProteinInfo.JSONize(),
             'Reaction': struct.ReactionInfo.JSONize(),
             'DensityConstraint': struct.DensityConstraintsInfo.JSONize(),
             'ProcessConstraint': struct.ProcessConstraintsInfo.JSONize(),
             'MetaboliteConstraint': struct.MetaboliteConstraintsInfo.JSONize(),
             'EnzymeConstraint': struct.EnzymeConstraintsInfo.JSONize()
             }

    return({'RBA_ModelData': {'StructuralInformation': Block}})


def generateOverview(StatsReactions, StatsMetabolites, StatsModules, StatsEnzymes, StatsProteins, StatsProcesses, StatsCompartments, StatsTargets, StatsConstraintsBiological, StatsConstraintsMathematical):
    out = {'nReactionsTotal': StatsReactions['nReactionsTotal'],
           'nReactionsUnique': StatsReactions['nReactionsUnique'],
           'nReactionsSpontaneous': StatsReactions['nReactionsSpontaneous'],
           'nReactionsEnzymatic': StatsReactions['nReactionsEnzymatic'],
           'nReactionsInternal': StatsReactions['nReactionsInternal'],
           'nReactionsExchange': StatsReactions['nReactionsExchange'],
           'nReactionsCompartmentTransport': StatsReactions['nReactionsCompartmentTransport'],
           'nReactionsReversible': StatsReactions['nReactionsReversible'],
           'nReactionsIrreversible': StatsReactions['nReactionsIrreversible'],
           'nMetabolitesTotal': StatsMetabolites['nMetabolitesTotal'],
           'nMetabolitesInternal': StatsMetabolites['nMetabolitesInternal'],
           'nMetabolitesExternal': StatsMetabolites['nMetabolitesExternal'],
           'nMetabolitesGrowthRelevant': StatsMetabolites['nMetabolitesGrowthRelevant'],
           'nBoundaryMetabolites': StatsMetabolites['nBoundaryMetabolites'],
           'nEnzymesTotal': StatsEnzymes['nEnzymesTotal'],
           'nEnzymesUnique': StatsEnzymes['nEnzymesUnique'],
           'nBioConstraintsMetabolite': StatsConstraintsBiological['nBioConstraintsMetabolite'],
           'nBioConstraintsCapacity': StatsConstraintsBiological['nBioConstraintsCapacity'],
           'nBioConstraintsProcess': StatsConstraintsBiological['nBioConstraintsProcess'],
           'nBioConstraintsDensity': StatsConstraintsBiological['nBioConstraintsDensity'],
           'nMathConstraintsVariables': StatsConstraintsMathematical['nMathConstraintsVariables'],
           'nMathConstraintsConstraints': StatsConstraintsMathematical['nMathConstraintsConstraints'],
           'nMathConstraintsEqualities': StatsConstraintsMathematical['nMathConstraintsEqualities'],
           'nMathConstraintsInequalities': StatsConstraintsMathematical['nMathConstraintsInequalities'],
           'nProteinsTotal': StatsProteins['nProteinsTotal'],
           'nProcessesTotal': StatsProcesses['nProcessesTotal'],
           'nModulesTotal': StatsModules['nModulesTotal'],
           'nCompartmentsTotal': StatsCompartments['nCompartmentsTotal']}
    out.update(StatsTargets)
    return(out)


def StatsConstraintsBiological(MetCs, CapCs, DenCs, EffCs):
    out = {'nBioConstraintsMetabolite': len(MetCs.keys()),
           'nBioConstraintsCapacity': len(CapCs.keys()),
           'nBioConstraintsProcess': len(EffCs.keys()),
           'nBioConstraintsDensity': len(DenCs.keys())}
    return(out)


def StatsConstraintsMathematical(MetCs, CapCs, DenCs, EffCs, Enzymes, Reactions, Processes):
    nVars = len(Enzymes.keys())+len(Reactions.keys())+len(Processes.keys())
    nConsts = len(MetCs.keys())+len(CapCs.keys())+len(DenCs.keys())+len(EffCs.keys())
    nEqC = 0
    nInC = 0
    for i in MetCs.keys():
        if MetCs[i]['Type'] == '<=':
            nInC += 1
        if MetCs[i]['Type'] == '=':
            nEqC += 1
    for i in CapCs.keys():
        if CapCs[i]['Type'] == '<=':
            nInC += 1
        if CapCs[i]['Type'] == '=':
            nEqC += 1
    for i in DenCs.keys():
        if DenCs[i]['Type'] == '<=':
            nInC += 1
        if DenCs[i]['Type'] == '=':
            nEqC += 1
    for i in EffCs.keys():
        if EffCs[i]['Type'] == '<=':
            nInC += 1
        if EffCs[i]['Type'] == '=':
            nEqC += 1
    out = {'nMathConstraintsVariables': nVars,
           'nMathConstraintsConstraints': nConsts,
           'nMathConstraintsEqualities': nEqC,
           'nMathConstraintsInequalities': nInC}
    return(out)


def findContainedProteins(comp, Prots):
    out = []
    for k in Prots.keys():
        if Prots[k]['Compartment'] == comp:
            out.append(k)
    return(out)


def findContainedEnzymes(ContainedProteins, Prots, Enzymes):
    rs = []
    enzs = []
    for k in ContainedProteins:
        rs = rs+Prots[k]['associatedReactions']
        enzs = enzs+Prots[k]['associatedEnzymes']
    out = {'containedEnzymes': list(numpy.unique(enzs)),
           'containedReactions': list(numpy.unique(rs))}
    return(out)


def findAssociatedEnzyme(Enzymes, protein):
    out1 = []
    out2 = []
    for e in Enzymes.keys():
        if protein in Enzymes[e]['Subunits'].keys():
            out1.append(e)
            out2.append(Enzymes[e]['Reaction'])
    out = {'Enz': out1,
           'Rx': out2}
    return(out)


def findIsozymes(ez, blocks, Reactions, rx):
    out = []
    twins = Reactions[rx]['Twins']
    if len(twins) > 0:
        for r in twins:
            if not type(r) == list:
                if not type(Reactions[r]['Enzyme']) == list:
                    out.append(Reactions[r]['Enzyme'])
    return(out)


def sortConstraints(matrix, model):
    metaboliteList = [model.metabolism.species._elements[m].id for m in range(
        len(model.metabolism.species._elements))]
    mets = {}
    capacity = {}
    efficiency = {}
    density = {}
    for j in range(len(matrix.row_names)):
        i = matrix.row_names[j]
        if i in metaboliteList:
            mets[i] = j
        if 'enzyme' in i:
            capacity[i] = j
        if i.startswith('P_'):
            efficiency[i] = j
        if '_density' in i:
            density[i] = j
    out = {'MetaboliteConsts': mets,
           'ProcessConsts': efficiency,
           'EnzymeConsts': capacity,
           'DensityConsts': density}
    return(out)


def JSON_Int64_compensation(o):
    if isinstance(o, numpy.int64):
        return int(o)
#    raise TypeError
