# python 2/3 compatibility
from __future__ import division, print_function
import sys
import os.path
import numpy
import pandas
import copy
import difflib
import scipy
import collections
import json

# package imports
#import rba
from . import rba
from .rba_SimulationData import RBA_SimulationData
from .rba_SimulationParameters import RBA_SimulationParameters
from .rba_ModelStructure import RBA_ModelStructure
from .rba_Problem import RBA_Problem
from .rba_Matrix import RBA_Matrix
from .rba_LP import RBA_LP
from .rba_FBA import RBA_FBA
from .rba_LogBook import RBA_LogBook
from .rba.model import RbaModel
from .rba.core.constraint_matrix import ConstraintMatrix
from .rba.core.solver import Solver



class RBA_Session(object):
    """
    Top level of the RBA API.

    Attributes
    ----------
    xml_dir : str
        Current Growth rate as numeric value
    model : rba.RbaModel
        Current Growth rate as numeric value
    matrices : rba.ConstraintMatrix
        Current Growth rate as numeric value
    solver : rba.Solver
        Current Growth rate as numeric value
    Problem : rbatools.RBA_Problem
        Current Growth rate as numeric value
    Medium : dict
        Current Growth rate as numeric value
    ModelStructure : rbatools.RBA_ModelStructure
        Current Growth rate as numeric value
    Results : dict
        Current Growth rate as numeric value
    Parameters : dict
        Current Growth rate as numeric value
    SimulationData : rbatools.RBA_SimulationData
        Current Growth rate as numeric value
    SimulationParameters : rbatools.RBA_SimulationParameters
        Current Growth rate as numeric value

    Methods
    ----------
    __init__(xml_dir)
        Creates controler object from files

    reloadModel()
        Reloads model from files

    recordResults(runName)
        Records Simulation output for further use.

    recordParameters(runName)
        Records Simulation parameters for further use.

    clearResults()
        Removes all previosly recorded results.

    clearParameters()
        Removes all previosly recorded parameters.

    writeResults(session_name='', digits=10)
        Creates SimulationData and SimulationParameters objects from recordings.

    setMu(Mu)
        Sets growth-rate to desired value.

    doSolve(runName='DontSave')
        Solves problem to find solution.

    findMaxGrowthRate(precision=0.00001, max=4, recording=False)
        Applies dichotomy-search to find the maximal feasible growth-rate.

    findMinMediumConcentration(metabolite, precision=0.00001, max=100, recording=False)
        Applies dichotomy-search to find the minimal feasible concentration of
        growth-substrate in medium.

    setMedium(changes)
        Sets the concentration of growth-substrate in medium.

    knockOut(gene)
        Emulates a gene knock out.

    FeasibleRange(*variables)
        Determines the feasible range of model variables.

    ConstraintSaturation(*constraints)
        Determines the saturation of model constraints at current solution.

    ParetoFront(variables, N)
        Determine Pareto front of two model variables.

    addProtein(input)
        Adds representation of individual proteins to problem.

    returnExchangeFluxes()

    """

    def __init__(self, xml_dir):
        """
        Creates RBA_Session object from files

        Parameters
        ----------
        xml_dir : str
            Path to the directory where rba-model files are located.
        """
        self.xml_dir = xml_dir
        self.LogBook = RBA_LogBook('Controler')

        if not hasattr(self, 'ModelStructure'):
            if os.path.isfile(str(self.xml_dir+'/ModelStructure.json')):
                self.ModelStructure = RBA_ModelStructure()
                with open(str(self.xml_dir+'/ModelStructure.json'), 'r') as myfile:
                    data = myfile.read()
                self.ModelStructure.fromJSON(inputString=data)
            else:
                self.build_ModelStructure()

        self.model = RbaModel.from_xml(input_dir=xml_dir)
        self.matrices = ConstraintMatrix(model=self.model)
        self.solver = Solver(matrix=self.matrices)

        self.LogBook.addEntry('Model loaded from {}.'.format(self.xml_dir))
        self.Problem = RBA_Problem(solver=self.solver)

        medium = pandas.read_csv(xml_dir+'/medium.tsv', sep='\t')
        self.Medium = dict(zip(list(medium.iloc[:, 0]), [float(i)
                                                         for i in list(medium.iloc[:, 1])]))

        self.Mu = self.Problem.Mu
        self.ExchangeMap = buildExchangeMap(self)
        self.ExchangeMap2 = buildExchangeMap2(self)

    def build_ModelStructure(self):
        self.ModelStructure = RBA_ModelStructure()
        self.ModelStructure.fromFiles(xml_dir=self.xml_dir)
        self.ModelStructure.exportJSON(path=self.xml_dir)

    def rebuild_from_model(self):
        """
        Rebuilds computational model-representation from own attribute "model" (rba.RbaModel-object).
        """
        self.LogBook.addEntry('Model rebuilt.')
        self.matrices = rba.ConstraintMatrix(model=self.model)
        self.solver = rba.Solver(matrix=self.matrices)
        self.Problem = RBA_Problem(solver=self.solver)
        self.setMedium(changes=self.Medium)

    def reloadModel(self):
        """
        Reloads model from xml-files and then rebuild computational model-representation.
        """
        self.LogBook.addEntry('Model reloaded from {}.'.format(self.xml_dir))
        self.model = rba.RbaModel.from_xml(input_dir=self.xml_dir)
        self.rebuild_from_model()

    def recordResults(self, runName):
        """
        Records Simulation output for further use.
        and strores them in own 'Results'-attribute as pandas.DataFrames in a dictionary with the respective run-name being a column in all DataFrames.

        Parameters
        ----------
        runName : str
            Name of observation.
            Serves as ID for all Data, originating from these.
        """
        self.LogBook.addEntry('Solution recorded under {}.'.format(runName))
        if not hasattr(self, 'Results'):
            self.Results = {'Reactions': pandas.DataFrame(index=list(self.ModelStructure.ReactionInfo.Elements.keys())),
                            'Enzymes': pandas.DataFrame(index=list(self.ModelStructure.EnzymeInfo.Elements.keys())),
                            'Processes': pandas.DataFrame(index=[self.ModelStructure.ProcessInfo.Elements[i]['ID']+'_machinery' for i in self.ModelStructure.ProcessInfo.Elements.keys()]),
                            'Proteins': pandas.DataFrame(index=list(self.ModelStructure.ProteinMatrix['Proteins'])),
                            'ProtoProteins': pandas.DataFrame(index=list(self.ModelStructure.ProteinGeneMatrix['ProtoProteins'])),
                            'Constraints': pandas.DataFrame(index=self.Problem.LP.row_names),
                            'SolutionType': pandas.DataFrame(index=['SolutionType']),
                            'Mu': pandas.DataFrame(index=['Mu']),
                            'ObjectiveFunction': pandas.DataFrame(index=self.Problem.LP.col_names),             
                            'ObjectiveValue': pandas.DataFrame(index=['ObjectiveValue']),
                            'ExchangeFluxes': pandas.DataFrame(index=list(self.ExchangeMap.keys()))}

        Exchanges = self.returnExchangeFluxes()
        for i in Exchanges.keys():
            self.Results['ExchangeFluxes'].loc[i, runName] = Exchanges[i]

        self.Results['Reactions'][runName] = self.Results['Reactions'].index.map(
            {i: self.Problem.SolutionValues[i] for i in list(self.Results['Reactions'].index)})
        self.Results['Enzymes'][runName] = self.Results['Enzymes'].index.map(
            {i: self.Problem.SolutionValues[i] for i in list(self.Results['Enzymes'].index)})
        self.Results['Processes'][runName] = self.Results['Processes'].index.map(
            {i: self.Problem.SolutionValues[i] for i in list(self.Results['Processes'].index)})
        self.Results['Constraints'][runName] = self.Results['Constraints'].index.map(
            {i: self.Problem.DualValues[i] for i in self.Problem.LP.row_names})
        self.Results['Proteins'][runName] = self.Results['Proteins'].index.map(
            ProteomeRecording(self, runName))
        self.Results['ProtoProteins'][runName] = self.Results['ProtoProteins'].index.map(
            ProtoProteomeRecording(self, runName, self.Results['Proteins']))
        self.Results['SolutionType'][runName] = self.Problem.SolutionType
        self.Results['Mu'][runName] = self.Problem.Mu
        self.Results['ObjectiveFunction'][runName] = list(self.Problem.getObjective().values())
        self.Results['ObjectiveValue'][runName] = self.Problem.ObjectiveValue

    def recordParameters(self, runName):
        """
        Records Simulation parameters (LP-coefficients etc.) for further use.
        and strores them in own 'Parameters'-attribute as pandas.DataFrames in a dictionary with the respective run-name being a column in all DataFrames.

        Parameters
        ----------
        runName : str
            Name of observation.
            Serves as ID for all Data, originating from these.
        """
        self.LogBook.addEntry('Coefficients recorded under {}.'.format(runName))
        EnzymeCapacities = self.get_parameter_values(
            parameter_type='enzyme_efficiencies', species=None, output_format='dict')
        ProcessCapacities = self.get_parameter_values(
            parameter_type='machine_efficiencies', species=None, output_format='dict')
        CompartmentCapacities = self.get_parameter_values(
            parameter_type='maximal_densities', species=None, output_format='dict')
        TargetValues = self.get_parameter_values(
            parameter_type='target_values', species=None, output_format='dict')
        if not hasattr(self, 'Parameters'):
            self.Parameters = {'EnzymeEfficiencies_FW': pandas.DataFrame(index=list(EnzymeCapacities.keys())),
                               'EnzymeEfficiencies_BW': pandas.DataFrame(index=list(EnzymeCapacities.keys())),
                               'ProcessEfficiencies': pandas.DataFrame(index=list(ProcessCapacities.keys())),
                               'CompartmentCapacities': pandas.DataFrame(index=list(CompartmentCapacities.keys())),
                               'Medium': pandas.DataFrame(index=self.Medium.keys()),
                               'TargetValues': pandas.DataFrame(index=[TargetValues[i]['Target_id'] for i in list(TargetValues.keys())])}
        
        self.Parameters['EnzymeEfficiencies_FW'][runName] = self.Parameters['EnzymeEfficiencies_FW'].index.map({i: list(
            EnzymeCapacities[i]['Forward'].values())[0] for i in list(EnzymeCapacities.keys()) if len(list(EnzymeCapacities[i]['Forward'].values())) > 0})
        self.Parameters['EnzymeEfficiencies_BW'][runName] = self.Parameters['EnzymeEfficiencies_BW'].index.map({i: list(
            EnzymeCapacities[i]['Backward'].values())[0] for i in list(EnzymeCapacities.keys()) if len(list(EnzymeCapacities[i]['Forward'].values())) > 0})
        self.Parameters['ProcessEfficiencies'][runName] = self.Parameters['ProcessEfficiencies'].index.map(
            {i: list(ProcessCapacities[i].values())[0] for i in list(ProcessCapacities.keys()) if len(list(ProcessCapacities[i].values())) > 0})
        self.Parameters['CompartmentCapacities'][runName] = self.Parameters['CompartmentCapacities'].index.map(
            {i: list(CompartmentCapacities[i].values())[0] for i in list(CompartmentCapacities.keys()) if len(list(CompartmentCapacities[i].values())) > 0})
        self.Parameters['Medium'][runName] = self.Parameters['Medium'].index.map(self.Medium)
        self.Parameters['TargetValues'][runName] = self.Parameters['TargetValues'].index.map(
            {TargetValues[i]['Target_id']: list(TargetValues[i]['Target_value'].values())[0] for i in list(TargetValues.keys()) if len(list(TargetValues[i]['Target_value'].values())) > 0})


    def clearResults(self):
        """
        Removes all previosly recorded results and deletes own 'Results'-attribute.
        """
        self.LogBook.addEntry('Results cleared.')

        delattr(self, 'Results')

    def clearParameters(self):
        """
        Removes all previosly recorded parameters and deletes own 'Parameters'-attribute.
        """

        self.LogBook.addEntry('Parameters cleared.')
        delattr(self, 'Parameters')

    def writeResults(self, session_name='', digits=5, loggingIntermediateSteps=False):
        """
        Creates SimulationData and SimulationParameters objects from recordings ('Results'.'Parameters').

        Stores them as rbatools.RBA_SimulationData
        and rbatools.RBA_SimulationParameters objects as attributes.
        Access via attributes .SimulationData and SimulationParameters respectively.

        Parameters
        ----------
        digits : int
            Number of decimal places in the numeric results
            Default: 10
        session_name : str
            Name of Simulation session.
            Default: ''
        """
        self.LogBook.addEntry('Data written under {}.'.format(session_name))
        if hasattr(self, 'Results'):
            self.Results['uniqueReactions'] = mapIsoReactions(Controller=self)
            self.Results['SolutionType'] = self.Results['SolutionType']
            self.Results['Mu'] = self.Results['Mu'].round(digits)
            self.Results['ObjectiveFunction'] = self.Results['ObjectiveFunction'].loc[(
                self.Results['ObjectiveFunction'] != 0).any(axis=1)].round(digits)
            self.Results['ObjectiveValue'] = self.Results['ObjectiveValue'].round(digits)
            self.Results['Proteins'] = self.Results['Proteins'].round(digits)
            self.Results['uniqueReactions'] = self.Results['uniqueReactions'].round(digits)
            self.Results['Reactions'] = self.Results['Reactions'].round(digits)
            self.Results['Enzymes'] = self.Results['Enzymes'].round(digits)
            self.Results['Processes'] = self.Results['Processes'].round(digits)
            self.Results['Constraints'] = self.Results['Constraints'].round(digits)
            self.Results['ExchangeFluxes'] = self.Results['ExchangeFluxes'].round(digits)

            self.SimulationData = RBA_SimulationData(StaticData=self.ModelStructure)
            self.SimulationData.fromSimulationResults(Controller=self, session_name=session_name)

        if hasattr(self, 'Parameters'):
            self.Parameters['EnzymeEfficiencies_FW'] = self.Parameters['EnzymeEfficiencies_FW'].round(
                digits)
            self.Parameters['EnzymeEfficiencies_BW'] = self.Parameters['EnzymeEfficiencies_BW'].round(
                digits)
            self.Parameters['ProcessEfficiencies'] = self.Parameters['ProcessEfficiencies'].round(
                digits)
            self.Parameters['CompartmentCapacities'] = self.Parameters['CompartmentCapacities'].round(
                digits)
            self.Parameters['TargetValues'] = self.Parameters['TargetValues'].round(digits)                
            self.Parameters['Medium'] = self.Parameters['Medium'].loc[(
                self.Parameters['Medium'] != 0).any(axis=1)].round(digits)                
            self.SimulationParameters = RBA_SimulationParameters(StaticData=self.ModelStructure)
            self.SimulationParameters.fromSimulationResults(Controller=self)

    def returnExchangeFluxes(self):
        """
        Generates a dictonary with the exchang-rates of boundary-metabolites.

        Returns
        -------
        Dictonary with exchange-keys and respective -rates.
        """
        out = {}
        for j in self.ExchangeMap.keys():
            netflux = 0
            for k in self.ExchangeMap[j].keys():
                netflux += self.ExchangeMap[j][k]*self.Problem.SolutionValues[k]
            if netflux != 0:
                out[j] = netflux
        return(out)

    def setMu(self, Mu, loggingIntermediateSteps=False):
        """
        Sets growth-rate to desired value.

        Parameters
        ----------
        Mu : float
            Growth rate
        """
        self.LogBook.addEntry('Growth-rate changed:{} --> {}'.format(self.Mu, float(Mu)))
        self.Problem.setMu(Mu=float(Mu), ModelStructure=self.ModelStructure,
                           logging=loggingIntermediateSteps)
        self.Mu = float(Mu)

    def doSolve(self, runName='DontSave', loggingIntermediateSteps=False):
        """
        Solves problem to find solution.

        Does the same as rbatools.RBA_Problem.solveLP().
        Just has some automatic option for results-recording.

        Parameters
        ----------
        runName : str
            Name of observation.
            Serves as ID for all data, originating from this run.
            Special values :
                'DontSave' : Results are not recorded
                'Auto' : Results are automatically recorded
                         and appended to existing ones.
                    Named with number.
                Any other string: Results are recorded under this name.
            Default: 'DontSave'
        """

        self.Problem.solveLP(logging=loggingIntermediateSteps)
        if self.Problem.Solved:
            if runName is not 'DontSave':
                if runName is 'Auto':
                    if hasattr(self, 'Results'):
                        name = str(self.Results['Reactions'].shape[1]+1)
                    if not hasattr(self, 'Results'):
                        name = '1'
                if runName is not 'Auto':
                    name = runName
                self.recordResults(runName=name)

    def findMaxGrowthRate(self, precision=0.00001, max=4, start_value=None, recording=False, loggingIntermediateSteps=False):
        """
        Applies dichotomy-search to find the maximal feasible growth-rate.

        Parameters
        ----------
        precision : float
            Numberic precision with which maximum is approximated.
            Default : 0.00001
        max : float
            Defines the highest growth rate to be screened for.
            Default=4
        recording : bool
            Records intermediate feasible solutions
            while approaching the maximum growth-rate.
            Default : False

        Returns
        -------
        maximum feasible growth rate as float.
        """

        minMu = 0
        maxMu = max
        if start_value is None:
            testMu = maxMu
        else:
            testMu = start_value
        iteration = 0
        while (maxMu - minMu) > precision:
            self.setMu(Mu=testMu)
            self.Problem.solveLP(logging=loggingIntermediateSteps)
            if self.Problem.Solved:
                iteration += 1
                if recording:
                    self.recordResults('DichotomyMu_iteration_'+str(iteration))
                minMu = testMu
            else:
                maxMu = testMu
            testMu = numpy.mean([maxMu, minMu])
        self.LogBook.addEntry('Maximal growth-rate found to be: {}.'.format(minMu))
        if minMu == max:
            print('Warning: Maximum growth rate might exceed specified range. Try rerunning this method with larger max-argument.')
        self.setMu(Mu=minMu)
        self.Problem.solveLP(logging=False)
        self.Problem.SolutionType = 'GrowthRate_maximization'
        return(minMu)

    def findMinMediumConcentration(self, metabolite, precision=0.00001, max=100, recording=False, loggingIntermediateSteps=False):
        """
        Applies dichotomy-search to find the minimal feasible concentration of
        growth-substrate in medium, at a previously set growth-rate.

        Parameters
        ----------
        metabolite : str
            ID of metabolite in medium.
        precision : float
            Numberic precision with which minimum is approximated.
            Default : 0.00001
        max : float
            Defines the highest concentration rate to be screened for.
            Default=100
        recording : bool
            Records intermediate feasible solutions
            while approaching the minimum concentration.
            Default : False

        Returns
        -------
        minimum feasible growth-substrate concentration as float.
        """

        minConc = 0.0
        maxConc = max
        testConc = minConc
        iteration = 0
        oldConc = self.Medium[metabolite]
        while (maxConc - minConc) > precision:
            self.setMedium(changes={metabolite: testConc})
            self.Problem.solveLP(logging=loggingIntermediateSteps)
            if self.Problem.Solved:
                iteration += 1
                if recording:
                    run_name = 'Dichotomy_'+metabolite+'_' + \
                        str(testConc)+'_iteration_'+str(iteration)
                    self.recordResults(run_name)
                maxConc = testConc
            else:
                minConc = testConc
            testConc = numpy.mean([maxConc, minConc])
        self.LogBook.addEntry(
            'Minimal required {} concentration found to be: {}.'.format(metabolite, maxConc))
        self.setMedium(changes={metabolite: oldConc})
        return(maxConc)

    def setMedium(self, changes, loggingIntermediateSteps=False):
        """
        Sets the concentration of specified growth-substrate(s) in medium.

        Parameters
        ----------
        changes : dict
            Keys : ID of metabolite(s) in medium.
            Values : New concention(s)
        """

        for species in (changes.keys()):
            self.Medium[species] = float(changes[species])

        self.Problem.ClassicRBAmatrix.set_medium(self.Medium)
        self.Problem.ClassicRBAmatrix.build_matrices(self.Mu)

        inputMatrix = RBA_Matrix()
        inputMatrix.loadMatrix(matrix=self.Problem.ClassicRBAmatrix)
        self.Problem.LP.updateMatrix(matrix=inputMatrix, Ainds=MediumDependentCoefficients_A(
            self), Binds=[], CTinds=[], LBinds=None, UBinds=None)

    def knockOut(self, gene, loggingIntermediateSteps=False):
        """
        Emulates a gene knock out.
        Constrains all variables in the LP-problem (enzymes, other machineries), which require this gene(s), to zero.

        Parameters
        ----------
        gene : str
            ID of model-protein to be knocked out.
            Can either be gene-identifier, represented as ID or ProtoID of proteins in rbatools.protein_bloc.ProteinBlock.Elements class (depends on whether protein-isoforms are considered).
        """

        if type(gene) is str:
            genes = [gene]
        if type(gene) is list:
            genes = gene
        isoform_genes = [g for g in genes if g in list(self.ModelStructure.ProteinInfo.Elements.keys(
        ))]+[i for g in genes for i in self.ModelStructure.ProteinInfo.Elements.keys() if self.ModelStructure.ProteinInfo.Elements[i]['ProtoID'] == g]
        for g in isoform_genes:
            self.LogBook.addEntry('Gene {} knocked out.'.format(g))
            ConsumersEnzymes = self.ModelStructure.ProteinInfo.Elements[g]['associatedEnzymes']
            for i in ConsumersEnzymes:
                LikeliestVarName = difflib.get_close_matches(i, self.Problem.LP.col_names, 1)[0]
                self.Problem.setLB(inputDict={LikeliestVarName: 0},
                                   logging=loggingIntermediateSteps)
                self.Problem.setUB(inputDict={LikeliestVarName: 0},
                                   logging=loggingIntermediateSteps)
            ConsumersProcess = self.ModelStructure.ProteinInfo.Elements[g]['SupportsProcess']
            for i in ConsumersProcess:
                LikeliestVarName = difflib.get_close_matches(
                    str(self.ModelStructure.ProcessInfo.Elements[i]['ID']+'_machinery'), self.Problem.LP.col_names, 1)[0]
                self.Problem.setLB(inputDict={LikeliestVarName: 0},
                                   logging=loggingIntermediateSteps)
                self.Problem.setUB(inputDict={LikeliestVarName: 0},
                                   logging=loggingIntermediateSteps)

    def FeasibleRange(self, *variables, loggingIntermediateSteps=False):
        """
        Determines the feasible range of model variables.

        Parameters
        ----------
        variables : str or list of str
            Specifies variable(s) for which the feasible range is to be determined.
            Optional input:
                If not provided all model-variables are taken

        Returns
        -------
        Dictionary with variable-names as keys and other dictionaries as values.
        The 'inner' dictionaries hold keys 'Min' and 'Max'
        with values representing lower and upper bound of feasible range respectively.
        E.g. : {'variableA':{'Min':42 , 'Max':9000},
                'variableB':{'Min':-9000 , 'Max':-42}}
        """

        if len(list(variables)) > 0:
            if isinstance(variables[0], list):
                VariablesInQuestion = variables[0]
            if isinstance(variables[0], str):
                VariablesInQuestion = [variables[0]]
        if len(list(variables)) == 0:
            VariablesInQuestion = self.Problem.LP.col_names
        out = {}
        for i in VariablesInQuestion:
            min = numpy.nan
            max = numpy.nan
            self.Problem.clearObjective(logging=loggingIntermediateSteps)
#            self.Problem.setObjectiveCoefficients(inputDict=dict(
#                zip(self.Problem.LP.col_names, [0.0]*len(self.Problem.LP.col_names))))
            self.Problem.setObjectiveCoefficients(
                inputDict={i: 1.0}, logging=loggingIntermediateSteps)
            self.Problem.solveLP(logging=loggingIntermediateSteps)
            if self.Problem.Solved:
                min = self.Problem.SolutionValues[i]
            self.Problem.setObjectiveCoefficients(
                inputDict={i: -1.0}, logging=loggingIntermediateSteps)
            self.Problem.solveLP(logging=loggingIntermediateSteps)
            if self.Problem.Solved:
                max = self.Problem.SolutionValues[i]
            out.update({i: {'Min': min, 'Max': max}})
            self.LogBook.addEntry(
                'Feasible-range of {} determined to be between {} and {}.'.format(i, min, max))

        return(out)

    def ConstraintSaturation(self, *constraints):
        """
        Determines the saturation of model constraints at current solution.

        Parameters
        ----------
        constraints : str or list of str
            Specifies constraints(s) for which the saturation is to be determined.
            Optional input:
                If not provided all model-constraints are taken

        Returns
        -------
        Pandas DataFrame with constraint-names as indices and the columns 'LHS', 'RHS', and 'Saturation'.
        'LHS': The sum over the respoctive constraint-row multiplied elementwise with the solution vector.
        'RHS': The value of the problem's righthand side, correesponding to the respective constraint.
        'Saturation': The saturation of the respective constraint ('LHS'/'RHS').
        (Equality constraints are always saturated)
        """
        if len(list(constraints)) > 0:
            if isinstance(constraints[0], list):
                ConstraintsInQuestion = constraints[0]
            if isinstance(constraints[0], str):
                ConstraintsInQuestion = [constraints[0]]
        if len(list(constraints)) == 0:
            ConstraintsInQuestion = self.Problem.LP.row_names

        rhs = self.Problem.getRighthandSideValue(ConstraintsInQuestion)
        lhs = self.Problem.calculateLefthandSideValue(ConstraintsInQuestion)
        RHS = list(rhs.values())
        LHS = list(lhs.values())
        Out = pandas.DataFrame(columns=['LHS', 'RHS', 'Saturation'], index=ConstraintsInQuestion)
        for i in ConstraintsInQuestion:
            lhval = LHS[self.Problem.LP.rowIndicesMap[i]]
            rhval = RHS[self.Problem.LP.rowIndicesMap[i]]
            sat = numpy.nan
            if rhval != 0:
                sat = lhval/rhval
            Out.loc[i, 'LHS'] = lhval
            Out.loc[i, 'RHS'] = rhval
            Out.loc[i, 'Saturation'] = sat
            self.LogBook.addEntry(
                'Saturation of constraint {} determined to be {}.'.format(i, sat))
        return(Out)

    def ParetoFront(self, variables, N, signV2='max', CropFeasibleRange=False, loggingIntermediateSteps=False):
        """
        Determine Pareto front of two model variables.

        Parameters
        ----------
        variables : list of str
            List holding two variable IDs between which the tradeoffs
            are determined as Pareto front.

        N : int
            Number of Increments within the feasible range of the first entry to
            argument 'variables'.
        signV2 : str
            'max': Variable 2 is maximised
            'min': Variable 2 is minimised
        Returns
        -------
        Pandas DataFrame with columns named after the two input variables
        and 'N' rows. Each row represents an interval on the Pareto front.
        Entries on each row are the X and Y coordinate on the Pareto front,
        representing the values of the two variables.
        """

        if variables[0] not in self.Problem.LP.col_names:
            print('Chosen Element not among problem variables')
            return
        if variables[1] not in self.Problem.LP.col_names:
            print('Chosen Element not among problem variables')
            return
        FR = self.FeasibleRange(variables[0])
        cMin = FR[variables[0]]['Min']
        cMax = FR[variables[0]]['Max']
        if CropFeasibleRange == float:
            if cMin < 0:
                cMin -= cMin*CropFeasibleRange
            elif cMin > 0:
                cMin += cMin*CropFeasibleRange
            if cMax < 0:
                cMax += cMax*CropFeasibleRange
            elif cMax > 0:
                cMax -= cMax*FeasibleRange
        elif CropFeasibleRange == int:
            if cMin < 0:
                cMin -= cMin*CropFeasibleRange
            elif cMin > 0:
                cMin += cMin*CropFeasibleRange
            if cMax < 0:
                cMax += cMax*CropFeasibleRange
            elif cMax > 0:
                cMax -= cMax*FeasibleRange

        concentrations = [float(cMin+(cMax-cMin)*i/N) for i in range(N+1)]

        Out = pandas.DataFrame(columns=[variables[0], variables[1]])
        oldLB = self.Problem.getLB(variables[0])
        oldUB = self.Problem.getUB(variables[0])
        iteration = -1
        for conc in concentrations:
            iteration += 1
            self.Problem.setLB(inputDict={variables[0]: conc}, logging=loggingIntermediateSteps)
            self.Problem.setUB(inputDict={variables[0]: conc}, logging=loggingIntermediateSteps)
            self.Problem.clearObjective(logging=loggingIntermediateSteps)
            if signV2 == 'max':
                self.Problem.setObjectiveCoefficients(
                    inputDict={variables[1]: -1}, logging=loggingIntermediateSteps)
            if signV2 == 'min':
                self.Problem.setObjectiveCoefficients(
                    inputDict={variables[1]: 1}, logging=loggingIntermediateSteps)
            self.Problem.solveLP(logging=loggingIntermediateSteps)
            if self.Problem.Solved:
                max = abs(self.Problem.ObjectiveValue)
            else:
                max = numpy.nan
            self.Problem.setLB(inputDict=oldLB, logging=loggingIntermediateSteps)
            self.Problem.setUB(inputDict=oldUB, logging=loggingIntermediateSteps)
            Out.loc[iteration, variables[0]] = conc
            Out.loc[iteration, variables[1]] = max
        self.LogBook.addEntry(
            'Pareto-front between {} and {} determined.'.format(variables[0], variables[1]))
        return(Out)

    def addExchangeReactions(self):
        """
        Adds explicit exchange-reactions of boundary-metabolites to RBA-problem, named R_EX_ followed by metabolite name (without M_ prefix).
        """
        Mets_external = [m.id for m in self.model.metabolism.species if m.boundary_condition]
        Mets_internal = [m.id for m in self.model.metabolism.species if not m.boundary_condition]
        Reactions = [r.id for r in self.model.metabolism.reactions]
        full_S = rba.core.metabolism.build_S(
            Mets_external+Mets_internal, self.model.metabolism.reactions)
        S_M_ext = full_S[:len(Mets_external), ].toarray()
        col_indices_toremove = []
        for i in range(S_M_ext.shape[1]):
            s_col_uniques = list(set(list(S_M_ext[:, i])))
            if len(s_col_uniques) == 1:
                if s_col_uniques[0] == 0:
                    col_indices_toremove.append(i)
        RemainingReactions = [i for i in Reactions if Reactions.index(
            i) not in col_indices_toremove]
        S_ext = numpy.delete(S_M_ext, col_indices_toremove, axis=1)
        A = numpy.concatenate((S_ext, numpy.eye(len(Mets_external))), axis=1, out=None)
        ColNames = RemainingReactions+[str('R_EX_'+i.split('M_')[-1]) for i in Mets_external]
        # print(str('R_EX_'+i.split('M_')[-1]))
        LBs = list([self.Problem.LP.LB[self.Problem.LP.col_names.index(i)]
                    for i in RemainingReactions]+[-10000]*len(Mets_external))
        UBs = list([self.Problem.LP.UB[self.Problem.LP.col_names.index(i)]
                    for i in RemainingReactions]+[10000]*len(Mets_external))
        b = [0]*len(Mets_external)
        f = list([self.Problem.LP.f[self.Problem.LP.col_names.index(i)]
                  for i in RemainingReactions]+[0]*len(Mets_external))

        ExchangeMatrix = RBA_Matrix()
        ExchangeMatrix.A = scipy.sparse.coo_matrix(A)
        ExchangeMatrix.b = numpy.array([0]*len(Mets_external))
        ExchangeMatrix.f = numpy.array(f)
        ExchangeMatrix.LB = numpy.array(LBs)
        ExchangeMatrix.UB = numpy.array(UBs)
        ExchangeMatrix.row_signs = ['E']*len(Mets_external)
        ExchangeMatrix.row_names = Mets_external
        ExchangeMatrix.col_names = ColNames
        ExchangeMatrix.mapIndices()
        self.Problem.LP.addMatrix(matrix=ExchangeMatrix)

        self.ExchangeReactionMap = dict(
            zip(Mets_external, [str('R_EX_'+i.split('M_')[-1]) for i in Mets_external]))

    def buildFBA(self, type='classic', objective='classic', maintenanceToBM=False):
        """
        Derives and constructs FBA-problem from RBA-problem and stores it under attribute 'FBA'.

        Parameters
        ----------
        type : str
        objective : str
        maintenanceToBM : boolean
        """
        RBAproblem = self.Problem.LP
        A = RBAproblem.A.toarray()
        if type == 'classic':
            Cols2remove = list(set([RBAproblem.col_names.index(i) for i in RBAproblem.col_names if not i.startswith('R_') and not i.startswith('M_') and not i.endswith('_synthesis')]
                                   + [RBAproblem.col_names.index(i) for i in RBAproblem.col_names if '_duplicate_' in i]
                                   + [RBAproblem.col_names.index(i) for i in RBAproblem.col_names if 'enzyme' in i]))
            Rows2remove = [RBAproblem.row_names.index(
                i) for i in RBAproblem.row_names if not i.startswith('M_')]
        elif type == 'parsi':
            Cols2remove = list(set([RBAproblem.col_names.index(i) for i in RBAproblem.col_names if not i.startswith(
                'R_') and not i.startswith('M_') and not i.endswith('_synthesis')]+[RBAproblem.col_names.index(i) for i in RBAproblem.col_names if '_duplicate_' in i]))
            Rows2remove = [RBAproblem.row_names.index(
                i) for i in RBAproblem.row_names if not i.startswith('R_') and not i.startswith('M_')]
        if objective == 'classic':
            if 'R_maintenance_atp' in RBAproblem.col_names:
                Cols2remove.append(RBAproblem.col_names.index('R_maintenance_atp'))
        Anew = numpy.delete(A, Cols2remove, axis=1)
        col_namesNew = list(numpy.delete(RBAproblem.col_names, Cols2remove))
        LBnew = numpy.delete(RBAproblem.LB, Cols2remove)
        UBnew = numpy.delete(RBAproblem.UB, Cols2remove)
        fNew = numpy.delete(RBAproblem.f, Cols2remove)
        Anew2 = numpy.delete(Anew, Rows2remove, axis=0)
        row_namesNew = list(numpy.delete(RBAproblem.row_names, Rows2remove))
        row_signsNew = list(numpy.delete(RBAproblem.row_signs, Rows2remove))
        bNew = numpy.delete(RBAproblem.b, Rows2remove)
        trnaInds = [i for i in range(len(row_namesNew)) if row_namesNew[i].startswith(
            'M_') and 'trna' in row_namesNew[i]]
        # bNew[trnaInds] = 0

        if objective == 'targets':
            col_namesNew.append('R_BIOMASS_targetsRBA')
            LBnew = numpy.append(LBnew, 0)
            UBnew = numpy.append(UBnew, 10000)
            fNew = numpy.append(fNew, 0)
            BMrxnCol = numpy.ones((len(row_namesNew), 1))
            BMrxnCol[:, 0] = bNew
            if maintenanceToBM:
                MaintenanceTarget = LBnew[col_namesNew.index('R_maintenance_atp')]
                BMrxnCol[row_namesNew.index('M_atp_c')] += MaintenanceTarget
                BMrxnCol[row_namesNew.index('M_h2o_c')] += MaintenanceTarget
                BMrxnCol[row_namesNew.index('M_adp_c')] -= MaintenanceTarget
                BMrxnCol[row_namesNew.index('M_pi_c')] -= MaintenanceTarget
                BMrxnCol[row_namesNew.index('M_h_c')] -= MaintenanceTarget
                LBnew[col_namesNew.index('R_maintenance_atp')] = 0
            Anew2 = numpy.append(Anew2, -BMrxnCol, axis=1)
            bNew = numpy.array([0]*Anew2.shape[0])

        Matrix1 = RBA_Matrix()
        Matrix1.A = scipy.sparse.coo_matrix(Anew2)
        Matrix1.b = bNew
        Matrix1.LB = LBnew
        Matrix1.UB = UBnew
        Matrix1.row_signs = row_signsNew
        Matrix1.row_names = row_namesNew
        Matrix1.col_names = col_namesNew
        Matrix1.f = fNew
        if type == 'classic':
            Matrix1.b = numpy.array([0]*len(row_signsNew))
            LP1 = RBA_LP()
            LP1.loadMatrix(Matrix1)
        elif type == 'parsi':
            MetaboliteRows = {i: Matrix1.row_names.index(
                i) for i in Matrix1.row_names if i.startswith('M_')}
            EnzymeCols = {i: Matrix1.col_names.index(
                i) for i in Matrix1.col_names if i.startswith('R_') and '_enzyme' in i}
            Matrix2 = RBA_Matrix()
            Matrix2.A = scipy.sparse.coo_matrix(numpy.zeros((len(MetaboliteRows), len(EnzymeCols))))
            Matrix2.b = numpy.array(Matrix1.b[list(MetaboliteRows.values())])
            Matrix2.LB = numpy.array(Matrix1.LB[list(EnzymeCols.values())])
            Matrix2.UB = numpy.array(Matrix1.UB[list(EnzymeCols.values())])
            Matrix2.f = numpy.array(Matrix1.f[list(EnzymeCols.values())])
            Matrix2.row_signs = [Matrix1.row_signs[i] for i in list(MetaboliteRows.values())]
            Matrix2.row_names = list(MetaboliteRows.keys())
            Matrix2.col_names = list(EnzymeCols.keys())
            Matrix2.mapIndices()
            Matrix1.b = numpy.array([0]*len(bNew))
            LP1 = RBA_LP()
            LP1.loadMatrix(Matrix1)
            LP1.updateMatrix(Matrix2)
        self.FBA = RBA_FBA(LP1)



    def addProtein(self, input):
        """
        Adds representation of individual proteins to problem.

        Parameters
        ----------
        input : dict or str
            If input is str it has to be the ID of a protein in the model.
            Then this protein is added to the problem an creates:
                One constraint named Protein_'ID' (equality).
                One variable named TotalLevel_'ID' representing the total amount.
                One variable named Free_'ID'_'respectiveCompartment', this
                represents the fraction of the protein not assuming any function.
                It however consumes resources for synthesis (precursors and processes),
                which are the same as defined in the model files.
                And takes up space i the compartment as specified in the model-files
                for the protein.

            If input is dict it has to have two keys; 'ID' and 'UnusedProteinFraction'.
            By specifying this input one can define that the unused franction of the protein
            can also reside in other compartments and which processes it requires.

            The value to 'ID' is the ID of a protein in the model.
            The value to 'UnusedProteinFraction' is another dictionary.
            This can have several keys which must be model-compartments.
            For each of the keys the value is a dict holding IDs of model-processes as Keys
            and process requirements as Values (numerical).
            This specifies which processes each of the compartment-species of the protein
            requires.

            This generates the same constraint and TotalLevel-variable as with the simple input,
            however a variable representing each of the compartment-species for the unused fraction
            is added and incorporates the specific process requirements.

            E.g: input = {'ID': 'proteinA',
                          'UnusedProteinFraction':{'Cytoplasm':{'Translation':100}, {'Folding':10}],
                                                   'Membrane':{'Translation':100}, {'Folding':20}, {'Secretion':100}
                                                   }
                          }
                This adds 'proteinA' to the model, where the unused fraction can reside either in
                the Cytoplasm or the Membrane. However while the cytosolic-species only requires the
                processes 'Translation' and 'Folding'; the membrane-bound species also requires 'Secretion'
                and occupies more folding capacity.
                Then the constraint 'Protein_proteinA' is added and the 3 variables
                'TotalLevel_proteinA', 'Free_proteinA_Cytoplasm' and 'Free_proteinA_Membrane'.
        """

        if type(input) is str:
            input = {'ID': input}

        if 'ID' not in list(input.keys()):
            print('Error, no protein ID provided')
            return

        if input['ID'] not in list(self.ModelStructure.ProteinInfo.Elements.keys()):
            print('Error, protein not in model')
            return

        if 'UnusedProteinFraction' not in list(input.keys()):
            input.update({'UnusedProteinFraction':
                          {self.ModelStructure.ProteinInfo.Elements[input['ID']]['Compartment']:
                           self.ModelStructure.ProteinInfo.Elements[input['ID']]['ProcessRequirements']}})

        self.LogBook.addEntry('Protein {} added with specifications {}.'.format(
            input['ID'], str(json.dumps(input))))

        Muindexlist = []

        ## Building RBA_Matrix-object for new constraint-row, representing protein ##
        UsedProtein = RBA_Matrix()
        UsedProtein.A = scipy.sparse.coo_matrix(
            buildUsedProteinConstraint(Controler=self, protein=input['ID']))
        UsedProtein.b = numpy.array([float(0)])
        UsedProtein.f = numpy.array(self.Problem.LP.f)
        UsedProtein.LB = numpy.array(self.Problem.LP.LB)
        UsedProtein.UB = numpy.array(self.Problem.LP.UB)
        UsedProtein.row_signs = ['E']
        UsedProtein.row_names = ['Protein_'+input['ID']]
        UsedProtein.col_names = self.Problem.LP.col_names
        ## Add used protein row to problem ##
        self.Problem.LP.addMatrix(matrix=UsedProtein)
        ## Add used protein row to reference Matrix (Mu == 1) ##
        self.Problem.MuOneMatrix.addMatrix(matrix=UsedProtein)

        ## Building RBA_Matrix-object for new variable-col, representing total level of protein ##
        TotProtein = RBA_Matrix()
        TotProtein.A = scipy.sparse.coo_matrix(numpy.array(numpy.matrix(
            numpy.array([float(0)]*self.Problem.LP.A.shape[0]+[float(-1)])).transpose()))
        TotProtein.f = numpy.array([float(0)])
        TotProtein.LB = numpy.array([float(0)])
        TotProtein.UB = numpy.array([float(100000.0)])
        TotProtein.b = numpy.array(list(self.Problem.LP.b)+list(UsedProtein.b))
        TotProtein.row_signs = self.Problem.LP.row_signs+UsedProtein.row_signs
        TotProtein.row_names = self.Problem.LP.row_names+UsedProtein.row_names
        TotProtein.col_names = ['TotalLevel_'+input['ID']]
        ## Add total protein col to problem ##
        self.Problem.LP.addMatrix(matrix=TotProtein)
        ## Add total protein col to reference Matrix (Mu == 1) ##
        self.Problem.MuOneMatrix.addMatrix(matrix=TotProtein)

        ## Building RBA_Matrix-object for new variable-col,##
        ## representing each compartment-species of the protein ##
        for comp_species in list(input['UnusedProteinFraction'].keys()):
            ## Initiate RBA_Matrix object##
            UnusedProtein = RBA_Matrix()
            UnusedProtein.col_names = ['Free_'+input['ID']+'_'+comp_species]

            ## Extract required processes for protein and the respective demand ##
            ProcIDs = list(input['UnusedProteinFraction'][comp_species].keys())
            Preq = list(input['UnusedProteinFraction'][comp_species].values())
            ProcessCost = dict(
                zip([self.ModelStructure.ProcessInfo.Elements[k]['ID'] for k in ProcIDs], Preq))

            ## Get required charged trna buildingblocks and their stoichiometry in protein ##
            composition = self.ModelStructure.ProteinInfo.Elements[input['ID']]['AAcomposition']
            ## Extract the composition of charged trnas in terms of metabolic species ##
            species = self.ModelStructure.ProcessInfo.Elements['Translation']['Components']
            ## Determine required metabolites and their stoichiometry in protein ##
            MetaboliteCost = buildCompositionofUnusedProtein(
                species=species, composition=composition)

            ## Assemble process and metabolite requirements into stoichiometric coloumn vector ##
            ## And add to RBA_Matrix object ##
            colToAdd = numpy.array(numpy.matrix(numpy.array(list(MetaboliteCost.values())+list(ProcessCost.values()) +
                                                            [float(1)]+[self.ModelStructure.ProteinInfo.Elements[input['ID']]['AAnumber']])).transpose())
            UnusedProtein.A = scipy.sparse.coo_matrix(colToAdd)
            ## Add other information to RBA_Matrix object ##
            UnusedProtein.row_names = list(MetaboliteCost.keys())+[str(pc+'_capacity') for pc in list(
                ProcessCost.keys())]+['Protein_'+input['ID']]+[str(comp_species + '_density')]
            UnusedProtein.b = numpy.zeros(len(UnusedProtein.row_names))
            UnusedProtein.row_signs = ['E']*len(UnusedProtein.row_names)
            UnusedProtein.LB = numpy.array([float(0)])
            UnusedProtein.UB = numpy.array([float(100000.0)])
            UnusedProtein.f = numpy.array([float(0)])
            self.ProteinDilutionIndices = list(
                zip(list(MetaboliteCost.keys()), UnusedProtein.col_names*len(list(MetaboliteCost.keys()))))
            ## Add free protein col to problem ##
            self.Problem.LP.addMatrix(matrix=UnusedProtein)
            ## Add free protein col to reference Matrix (Mu == 1) ##
            self.Problem.MuOneMatrix.addMatrix(matrix=UnusedProtein)

            ## Find coefficients of unused protein column, subject to dilution (Metabolite and Process cost) ##
            ## And add them to MuDepIndices_A ##
            nonZeroEntries = numpy.where(UnusedProtein.A != 0)[0]
            self.Problem.MuDepIndices_A += [(UnusedProtein.row_names[i], UnusedProtein.col_names[0]) for i in nonZeroEntries if UnusedProtein.row_names[i]
                                            != 'Protein_'+input['ID'] and UnusedProtein.row_names[i] not in self.Problem.CompartmentDensities]

            self.setMu(self.Problem.Mu)

    ## !!! ##
    def eukaryoticDensities(self, totalAA=3.1, CompartmentRelationships=True, CompartmentComponents=False):
        Compartments = ['n', 'mIM', 'vM', 'mIMS', 'm', 'erM', 'mOM', 'x', 'c', 'cM', 'gM']
        Signs = ['L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L']
        totalAA = 3.1*0.71
        m_mIM = 0.66
        m_mIMS = 2
        m_mOM = 8
        DensityIndices = [self.Problem.LP.row_names.index(
            i) for i in self.Problem.CompartmentDensities]

        A = self.Problem.LP.A.toarray()
        A[numpy.min(DensityIndices):numpy.max(DensityIndices)+1, :] /= totalAA
        self.Problem.LP.A = scipy.sparse.coo_matrix(A)

        A0 = self.Problem.MuOneMatrix.A.toarray()
        A0[numpy.min(DensityIndices):numpy.max(DensityIndices)+1, :] /= totalAA
        self.Problem.MuOneMatrix.A = scipy.sparse.coo_matrix(A0)

        CompartmentMatrix = RBA_Matrix()
        A = numpy.ones((len(Compartments)+1, len(Compartments)))
        Eye = -numpy.eye(len(Compartments))
        A[0:len(Compartments), :] = Eye
        CompartmentMatrix.A = scipy.sparse.coo_matrix(A)
        CompartmentMatrix.b = numpy.array([float(0)]*len(Compartments)+[float(1)])
        CompartmentMatrix.f = numpy.array([float(0)]*len(Compartments))
        CompartmentMatrix.LB = numpy.array([float(0)]*len(Compartments))
        CompartmentMatrix.UB = numpy.array([float(1)]*len(Compartments))
        CompartmentMatrix.row_signs = ['L']*len(Compartments)+['E']
        # CompartmentMatrix.row_signs = ['E']*(len(Compartments)+1)
        CompartmentMatrix.row_names = ['n_density', 'mIM_density', 'vM_density', 'mIMS_density', 'm_density',
                                       'erM_density', 'mOM_density', 'x_density', 'cM_density', 'gM_density', 'c_density', 'TotalCapacity']
        CompartmentMatrix.col_names = ['F_n', 'F_mIM', 'F_vM', 'F_mIMS',
                                       'F_m', 'F_erM', 'F_mOM', 'F_x', 'F_cM', 'F_gM', 'F_c']
        # CompartmentMatrix.row_signs[CompartmentMatrix.col_names.index('F_m')]='E'

        if CompartmentRelationships:
            Anew = numpy.zeros((A.shape[0]+3, A.shape[1]))
            Anew[0:A.shape[0], :] = A
            CompartmentMatrix.row_names += ['m_mIM', 'm_mIMS', 'm_mOM']
            CompartmentMatrix.row_signs += ['E', 'E', 'E']
            CompartmentMatrix.b = numpy.array(list(CompartmentMatrix.b)+[float(0)]*3)
            Anew[CompartmentMatrix.row_names.index(
                'm_mIM'), CompartmentMatrix.col_names.index('F_m')] = float(1)
            Anew[CompartmentMatrix.row_names.index(
                'm_mIMS'), CompartmentMatrix.col_names.index('F_m')] = float(1)
            Anew[CompartmentMatrix.row_names.index(
                'm_mOM'), CompartmentMatrix.col_names.index('F_m')] = float(1)
            Anew[CompartmentMatrix.row_names.index(
                'm_mIM'), CompartmentMatrix.col_names.index('F_mIM')] = -m_mIM
            Anew[CompartmentMatrix.row_names.index(
                'm_mIMS'), CompartmentMatrix.col_names.index('F_mIMS')] = -m_mIMS
            Anew[CompartmentMatrix.row_names.index(
                'm_mOM'), CompartmentMatrix.col_names.index('F_mOM')] = -m_mOM
            CompartmentMatrix.A = scipy.sparse.coo_matrix(Anew)

        self.Problem.LP.addMatrix(matrix=CompartmentMatrix)
        self.Problem.MuOneMatrix.addMatrix(matrix=CompartmentMatrix)

        if CompartmentComponents:
            AlipidsA = numpy.zeros((7, len(Compartments)))
            Alipids = RBA_Matrix()
            Alipids.col_names = ['F_n', 'F_mIM', 'F_vM', 'F_mIMS',
                                 'F_m', 'F_erM', 'F_mOM', 'F_x', 'F_cM', 'F_gM', 'F_c']
            Alipids.row_names = ['M_pc_SC_c', 'M_pe_SC_c', 'M_ptd1ino_SC_c',
                                 'M_ps_SC_c', 'M_clpn_SC_m', 'M_pa_SC_c', 'M_ergst_c']
            Alipids.row_signs += ['E', 'E', 'E', 'E', 'E', 'E', 'E']
            Alipids.b = numpy.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
            Alipids.LB = numpy.array([float(0)]*len(Compartments))
            Alipids.UB = numpy.array([float(1)]*len(Compartments))
            Alipids.f = numpy.array([float(0)]*len(Compartments))
            AlipidsA[Alipids.row_names.index('M_pc_SC_c'), Alipids.col_names.index(
                'F_mIM')] = -0.0000883*totalAA
            AlipidsA[Alipids.row_names.index('M_pe_SC_c'), Alipids.col_names.index(
                'F_mIM')] = -0.00005852*totalAA
            AlipidsA[Alipids.row_names.index('M_ptd1ino_SC_c'),
                     Alipids.col_names.index('F_mIM')] = -0.00003377*totalAA
            AlipidsA[Alipids.row_names.index('M_ps_SC_c'), Alipids.col_names.index(
                'F_mIM')] = -0.00000873*totalAA
            AlipidsA[Alipids.row_names.index('M_clpn_SC_m'),
                     Alipids.col_names.index('F_mIM')] = -0.00002*totalAA
            AlipidsA[Alipids.row_names.index('M_pa_SC_c'), Alipids.col_names.index(
                'F_mIM')] = -0.0000039*totalAA
            AlipidsA[Alipids.row_names.index(
                'M_ergst_c'), Alipids.col_names.index('F_mIM')] = -0.008547*totalAA
            self.Problem.MuDepIndices_A += [('M_pc_SC_c', 'F_mIM'), ('M_pe_SC_c', 'F_mIM'), ('M_ptd1ino_SC_c', 'F_mIM'),
                                            ('M_ps_SC_c', 'F_mIM'), ('M_clpn_SC_m', 'F_mIM'), ('M_pa_SC_c', 'F_mIM'), ('M_ergst_c', 'F_mIM')]

            AlipidsA[Alipids.row_names.index(
                'M_pc_SC_c'), Alipids.col_names.index('F_mOM')] = -0.000636*totalAA
            AlipidsA[Alipids.row_names.index('M_pe_SC_c'), Alipids.col_names.index(
                'F_mOM')] = -0.0004822*totalAA
            AlipidsA[Alipids.row_names.index('M_ptd1ino_SC_c'),
                     Alipids.col_names.index('F_mOM')] = -0.0001289*totalAA
            AlipidsA[Alipids.row_names.index('M_ps_SC_c'), Alipids.col_names.index(
                'F_mOM')] = -0.0000167*totalAA
            AlipidsA[Alipids.row_names.index('M_clpn_SC_m'), Alipids.col_names.index(
                'F_mOM')] = -0.00004467*totalAA
            AlipidsA[Alipids.row_names.index('M_pa_SC_c'), Alipids.col_names.index(
                'F_mOM')] = -0.0000696*totalAA

            self.Problem.MuDepIndices_A += [('M_pc_SC_c', 'F_mOM'), ('M_pe_SC_c', 'F_mOM'), ('M_ptd1ino_SC_c',
                                                                                             'F_mOM'), ('M_ps_SC_c', 'F_mOM'), ('M_clpn_SC_m', 'F_mOM'), ('M_pa_SC_c', 'F_mOM')]

            Alipids.A = scipy.sparse.coo_matrix(AlipidsA)
            Alipids.mapIndices()
            self.Problem.LP.updateMatrix(Alipids, Ainds=[('M_pc_SC_c', 'F_mIM'), ('M_pe_SC_c', 'F_mIM'), ('M_ptd1ino_SC_c', 'F_mIM'), ('M_ps_SC_c', 'F_mIM'), ('M_clpn_SC_m', 'F_mIM'), ('M_pa_SC_c', 'F_mIM'), (
                'M_ergst_c', 'F_mIM'), ('M_pc_SC_c', 'F_mOM'), ('M_pe_SC_c', 'F_mOM'), ('M_ptd1ino_SC_c', 'F_mOM'), ('M_ps_SC_c', 'F_mOM'), ('M_clpn_SC_m', 'F_mOM'), ('M_pa_SC_c', 'F_mOM')])
            self.Problem.MuOneMatrix.updateMatrix(Alipids, Ainds=[('M_pc_SC_c', 'F_mIM'), ('M_pe_SC_c', 'F_mIM'), ('M_ptd1ino_SC_c', 'F_mIM'), ('M_ps_SC_c', 'F_mIM'), ('M_clpn_SC_m', 'F_mIM'), (
                'M_pa_SC_c', 'F_mIM'), ('M_ergst_c', 'F_mIM'), ('M_pc_SC_c', 'F_mOM'), ('M_pe_SC_c', 'F_mOM'), ('M_ptd1ino_SC_c', 'F_mOM'), ('M_ps_SC_c', 'F_mOM'), ('M_clpn_SC_m', 'F_mOM'), ('M_pa_SC_c', 'F_mOM')])

    ## !!! ##
    def eukaryoticDensities2(self, totalAA=3.1, CompartmentRelationships=True, CompartmentComponents=False):
        Compartments = ['n', 'mIM', 'vM', 'mIMS', 'm', 'erM', 'mOM', 'x', 'c', 'cM', 'gM']
        totalAA = 3.1*0.69
        m_mIM = 1.11
        m_mIMS = 0.7
        m_mOM = 7.2
        DensityIndices = [self.Problem.LP.row_names.index(
            i) for i in self.Problem.CompartmentDensities]

        A = self.Problem.LP.A.toarray()
        A[numpy.min(DensityIndices):numpy.max(DensityIndices)+1, :] /= totalAA
        self.Problem.LP.A = scipy.sparse.coo_matrix(A)

        A0 = self.Problem.MuOneMatrix.A.toarray()
        A0[numpy.min(DensityIndices):numpy.max(DensityIndices)+1, :] /= totalAA
        self.Problem.MuOneMatrix.A = scipy.sparse.coo_matrix(A0)

        CompartmentMatrix = RBA_Matrix()
        A = numpy.ones((len(Compartments)+1, len(Compartments)))
        Eye = -numpy.eye(len(Compartments))
        A[0:len(Compartments), :] = Eye
        CompartmentMatrix.A = scipy.sparse.coo_matrix(A)
        CompartmentMatrix.b = numpy.array([float(0)]*len(Compartments)+[float(1)])
        CompartmentMatrix.f = numpy.array([float(0)]*len(Compartments))
        CompartmentMatrix.LB = numpy.array([float(0)]*len(Compartments))
        CompartmentMatrix.UB = numpy.array([float(1)]*len(Compartments))
        CompartmentMatrix.row_signs = ['L']*(len(Compartments)+1)
        # CompartmentMatrix.row_signs = ['E']*(len(Compartments)+1)
        CompartmentMatrix.row_names = ['n_density', 'mIM_density', 'vM_density', 'mIMS_density', 'm_density',
                                       'erM_density', 'mOM_density', 'x_density', 'cM_density', 'gM_density', 'c_density', 'TotalCapacity']
        CompartmentMatrix.col_names = ['F_n', 'F_mIM', 'F_vM', 'F_mIMS',
                                       'F_m', 'F_erM', 'F_mOM', 'F_x', 'F_cM', 'F_gM', 'F_c']
        # CompartmentMatrix.row_signs[CompartmentMatrix.col_names.index('F_m')]='E'

        if CompartmentRelationships:
            Anew = numpy.zeros((A.shape[0]+3, A.shape[1]))
            Anew[0:A.shape[0], :] = A
            CompartmentMatrix.row_names += ['m_mIM', 'm_mIMS', 'm_mOM']
            CompartmentMatrix.row_signs += ['E', 'E', 'E']
            CompartmentMatrix.b = numpy.array(list(CompartmentMatrix.b)+[float(0)]*3)
            Anew[CompartmentMatrix.row_names.index(
                'm_mIM'), CompartmentMatrix.col_names.index('F_m')] = float(1)
            Anew[CompartmentMatrix.row_names.index(
                'm_mIMS'), CompartmentMatrix.col_names.index('F_m')] = float(1)
            Anew[CompartmentMatrix.row_names.index(
                'm_mOM'), CompartmentMatrix.col_names.index('F_m')] = float(1)
            Anew[CompartmentMatrix.row_names.index(
                'm_mIM'), CompartmentMatrix.col_names.index('F_mIM')] = -m_mIM
            Anew[CompartmentMatrix.row_names.index(
                'm_mIMS'), CompartmentMatrix.col_names.index('F_mIMS')] = -m_mIMS
            Anew[CompartmentMatrix.row_names.index(
                'm_mOM'), CompartmentMatrix.col_names.index('F_mOM')] = -m_mOM
            CompartmentMatrix.A = scipy.sparse.coo_matrix(Anew)

        self.Problem.LP.addMatrix(matrix=CompartmentMatrix)
        self.Problem.MuOneMatrix.addMatrix(matrix=CompartmentMatrix)

        if CompartmentComponents:

            PC_mIM = 0.0000883
            PE_mIM = 0.00005852
            PI_mIM = 0.00003377
            PS_mIM = 0.00000873
            CL_mIM = 0.00002
            PA_mIM = 0.0000039
            ES_mIM = 0.008547

            PC_mOM = 0.000636
            PE_mOM = 0.0004822
            PI_mOM = 0.0001289
            PS_mOM = 0.0000167
            CL_mOM = 0.00004467
            PA_mOM = 0.0000696
            ES_mOM = 0.0

            ConstraintMatrix = numpy.zeros((7, 0))
            Alipids = RBA_Matrix()
            Alipids.col_names = []
            Alipids.row_names = ['M_pc_SC_c', 'M_pe_SC_c', 'M_ptd1ino_SC_c',
                                 'M_ps_SC_c', 'M_clpn_SC_m', 'M_pa_SC_c', 'M_ergst_c']
            Alipids.row_signs = [
                self.Problem.LP.row_signs[self.Problem.LP.row_names.index(i)] for i in Alipids.row_names]
            Alipids.b = numpy.array(
                [self.Problem.LP.b[self.Problem.LP.row_names.index(i)] for i in Alipids.row_names])
            Alipids.LB = numpy.array([])
            Alipids.UB = numpy.array([])
            Alipids.f = numpy.array([])
            MudepIndices = []
            for pc in self.ModelStructure.ProcessInfo.Elements.keys():
                if self.ModelStructure.ProcessInfo.Elements[pc]['ID'] not in self.Problem.LP.col_names:
                    continue
                ConstraintMatrixNew = numpy.zeros(
                    (ConstraintMatrix.shape[0], ConstraintMatrix.shape[1]+1))
                ConstraintMatrixNew[:, 0:ConstraintMatrix.shape[1]] = ConstraintMatrix
                Alipids.col_names.append(self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
#                Alipids.LB = numpy.array(list(Alipids.LB).append(list(self.Problem.LP.LB)[
#                                         self.Problem.LP.col_names.index(self.ModelStructure.ProcessInfo.Elements[pc]['ID'])]))
#                Alipids.UB = numpy.array(list(Alipids.UB).append(list(self.Problem.LP.UB)[
#                                         self.Problem.LP.col_names.index(self.ModelStructure.ProcessInfo.Elements[pc]['ID'])]))
#                Alipids.f = numpy.array(list(Alipids.f).append(list(self.Problem.LP.f)[
#                                        self.Problem.LP.col_names.index(self.ModelStructure.ProcessInfo.Elements[pc]['ID'])]))
                Alipids.LB = numpy.concatenate([Alipids.LB, numpy.array(
                    list(self.Problem.LP.LB)[self.Problem.LP.col_names.index(self.ModelStructure.ProcessInfo.Elements[pc]['ID'])])])
                Alipids.UB = numpy.concatenate([Alipids.UB, numpy.array(
                    list(self.Problem.LP.UB)[self.Problem.LP.col_names.index(self.ModelStructure.ProcessInfo.Elements[pc]['ID'])])])
                Alipids.f = numpy.concatenate([Alipids.f, numpy.array(
                    list(self.Problem.LP.f)[self.Problem.LP.col_names.index(self.ModelStructure.ProcessInfo.Elements[pc]['ID'])])])
                for p in self.ModelStructure.ProcessInfo.Elements[pc]['Composition'].keys():
                    lE = sum(list(self.ModelStructure.ProteinInfo.Elements[p]['AAcomposition'].values(
                    )))*self.ModelStructure.ProcessInfo.Elements[pc]['Composition'][p]
                    if self.ModelStructure.ProteinInfo.Elements[p]['Compartment'] == 'mOM':
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pc_SC_c'), ConstraintMatrix.shape[1]] -= PC_mOM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pe_SC_c'), ConstraintMatrix.shape[1]] -= PE_mOM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_ptd1ino_SC_c'), ConstraintMatrix.shape[1]] -= PI_mOM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_ps_SC_c'), ConstraintMatrix.shape[1]] -= PS_mOM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_clpn_SC_m'), ConstraintMatrix.shape[1]] -= CL_mOM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pa_SC_c'), ConstraintMatrix.shape[1]] -= PA_mOM*lE/totalAA
                        MudepIndices += ('M_pc_SC_c',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                        MudepIndices += ('M_pe_SC_c',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                        MudepIndices += ('M_ptd1ino_SC_c',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                        MudepIndices += ('M_ps_SC_c',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                        MudepIndices += ('M_clpn_SC_m',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                        MudepIndices += ('M_pa_SC_c',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                    elif self.ModelStructure.ProteinInfo.Elements[p]['Compartment'] == 'mIM':
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pc_SC_c'), ConstraintMatrix.shape[1]] -= PC_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pe_SC_c'), ConstraintMatrix.shape[1]] -= PE_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_ptd1ino_SC_c'), ConstraintMatrix.shape[1]] -= PI_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_ps_SC_c'), ConstraintMatrix.shape[1]] -= PS_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_clpn_SC_m'), ConstraintMatrix.shape[1]] -= CL_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pa_SC_c'), ConstraintMatrix.shape[1]] -= PA_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_ergst_c'), ConstraintMatrix.shape[1]] -= ES_mIM*lE/totalAA
                        MudepIndices += ('M_pc_SC_c',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                        MudepIndices += ('M_pe_SC_c',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                        MudepIndices += ('M_ptd1ino_SC_c',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                        MudepIndices += ('M_ps_SC_c',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                        MudepIndices += ('M_clpn_SC_m',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                        MudepIndices += ('M_pa_SC_c',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                        MudepIndices += ('M_ergst_c',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                ConstraintMatrix = ConstraintMatrixNew

            for e in self.ModelStructure.EnzymeInfo.Elements.keys():
                if e not in self.Problem.LP.col_names:
                    continue
                ConstraintMatrixNew = numpy.zeros(
                    (ConstraintMatrix.shape[0], ConstraintMatrix.shape[1]+1))
                ConstraintMatrixNew[:, 0:ConstraintMatrix.shape[1]] = ConstraintMatrix
                Alipids.col_names.append(e)
#                xnew = list(self.Problem.LP.LB)[self.Problem.LP.col_names.index(e)]
                Alipids.LB = numpy.concatenate([Alipids.LB, numpy.array(
                    list(self.Problem.LP.LB)[self.Problem.LP.col_names.index(e)])])
                Alipids.UB = numpy.concatenate([Alipids.UB, numpy.array(
                    list(self.Problem.LP.UB)[self.Problem.LP.col_names.index(e)])])
                Alipids.f = numpy.concatenate([Alipids.f, numpy.array(
                    list(self.Problem.LP.f)[self.Problem.LP.col_names.index(e)])])
                # Alipids.LB = numpy.array(list(Alipids.LB).append(xnew))
#                Alipids.UB = numpy.array(list(Alipids.UB).append(
#                    list(self.Problem.LP.UB)[self.Problem.LP.col_names.index(e)]))
                # Alipids.f = numpy.array(list(Alipids.f).append(
                #    list(self.Problem.LP.f)[self.Problem.LP.col_names.index(e)]))
                for p in self.ModelStructure.EnzymeInfo.Elements[e]['Subunits'].keys():
                    lE = sum(
                        list(self.ModelStructure.ProteinInfo.Elements[p]['AAcomposition'].values()))
                    lE *= self.ModelStructure.EnzymeInfo.Elements[e]['Subunits'][p]['StochFac']
                    if self.ModelStructure.ProteinInfo.Elements[p]['Compartment'] == 'mOM':
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pc_SC_c'), ConstraintMatrix.shape[1]] -= PC_mOM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pe_SC_c'), ConstraintMatrix.shape[1]] -= PE_mOM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_ptd1ino_SC_c'), ConstraintMatrix.shape[1]] -= PI_mOM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_ps_SC_c'), ConstraintMatrix.shape[1]] -= PS_mOM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_clpn_SC_m'), ConstraintMatrix.shape[1]] -= CL_mOM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pa_SC_c'), ConstraintMatrix.shape[1]] -= PA_mOM*lE/totalAA
                        MudepIndices += ('M_pc_SC_c', e)
                        MudepIndices += ('M_pe_SC_c', e)
                        MudepIndices += ('M_ptd1ino_SC_c', e)
                        MudepIndices += ('M_ps_SC_c', e)
                        MudepIndices += ('M_clpn_SC_m', e)
                        MudepIndices += ('M_pa_SC_c', e)
                    elif self.ModelStructure.ProteinInfo.Elements[p]['Compartment'] == 'mIM':
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pc_SC_c'), ConstraintMatrix.shape[1]] -= PC_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pe_SC_c'), ConstraintMatrix.shape[1]] -= PE_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_ptd1ino_SC_c'), ConstraintMatrix.shape[1]] -= PI_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_ps_SC_c'), ConstraintMatrix.shape[1]] -= PS_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_clpn_SC_m'), ConstraintMatrix.shape[1]] -= CL_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pa_SC_c'), ConstraintMatrix.shape[1]] -= PA_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_ergst_c'), ConstraintMatrix.shape[1]] -= ES_mIM*lE/totalAA
                        MudepIndices += ('M_pc_SC_c', e)
                        MudepIndices += ('M_pe_SC_c', e)
                        MudepIndices += ('M_ptd1ino_SC_c', e)
                        MudepIndices += ('M_ps_SC_c', e)
                        MudepIndices += ('M_clpn_SC_m', e)
                        MudepIndices += ('M_pa_SC_c', e)
                        MudepIndices += ('M_ergst_c', e)
                ConstraintMatrix = ConstraintMatrixNew

            self.Problem.MuDepIndices_A += MudepIndices
            Alipids.A = scipy.sparse.coo_matrix(ConstraintMatrix)
            Alipids.mapIndices()
            self.Problem.LP.updateMatrix(Alipids, Ainds=MudepIndices)
            self.Problem.LP.updateMatrix(MuOneMatrix, Ainds=MudepIndices)

    ## !!! ##
    def eukaryoticDensities3(self, totalAA=3.1, VolumeFraction=False, CompartmentRelationships=True, CompartmentComponents=False):
        Compartments = ['n', 'mIM', 'vM', 'mIMS', 'm', 'erM', 'mOM', 'x', 'c', 'cM', 'gM']
        totalAA = 3.1*0.91
        m_mIM = 0.66
        m_mIMS = 2
        m_mOM = 8
        DensityIndices = [self.Problem.LP.row_names.index(
            i) for i in self.Problem.CompartmentDensities]

        A = self.Problem.LP.A.toarray()
        # A[numpy.min(DensityIndices):numpy.max(DensityIndices)+1, :] /= totalAA
        # self.Problem.LP.A = scipy.sparse.coo_matrix(A)

        A0 = self.Problem.MuOneMatrix.A.toarray()
        # A0[numpy.min(DensityIndices):numpy.max(DensityIndices)+1, :] /= totalAA
        # self.Problem.MuOneMatrix.A = scipy.sparse.coo_matrix(A0)

        OccupationMatrix = RBA_Matrix()
        # A = numpy.ones((len(Compartments)+1, len(Compartments)))
        A = -numpy.eye(len(Compartments))
        # Eye = -numpy.eye(len(Compartments))
        # A[0:len(Compartments), :] = Eye
        OccupationMatrix.A = scipy.sparse.coo_matrix(A)
#        OccupationMatrix.b = numpy.array([-0.209*totalAA]+[float(0)]*(len(Compartments)-1)+[float(totalAA)])
        OccupationMatrix.b = numpy.array([-0.209*totalAA]+[float(0)]*(len(Compartments)-1))
        OccupationMatrix.f = numpy.array([float(0)]*len(Compartments))
        OccupationMatrix.LB = numpy.array([float(0)]*len(Compartments))
        OccupationMatrix.UB = numpy.array([float(totalAA)]*len(Compartments))
        # OccupationMatrix.row_signs = ['E']*(len(Compartments))+['L']
        OccupationMatrix.row_signs = ['E']*(len(Compartments))
        # OccupationMatrix.row_names = ['n_density', 'mIM_density', 'vM_density', 'mIMS_density', 'm_density',
        #                              'erM_density', 'mOM_density', 'x_density', 'cM_density', 'gM_density', 'c_density', 'TotalProtein']
        OccupationMatrix.row_names = ['n_density', 'mIM_density', 'vM_density', 'mIMS_density',
                                      'm_density', 'erM_density', 'mOM_density', 'x_density', 'cM_density', 'gM_density', 'c_density']
        OccupationMatrix.col_names = ['O_n', 'O_mIM', 'O_vM', 'O_mIMS',
                                      'O_m', 'O_erM', 'O_mOM', 'O_x', 'O_cM', 'O_gM', 'O_c']
        # CompartmentMatrix.row_signs[CompartmentMatrix.col_names.index('F_m')]='E'
        OccupationMatrix.mapIndices()
        self.Problem.LP.addMatrix(matrix=OccupationMatrix)
        self.Problem.MuOneMatrix.addMatrix(matrix=OccupationMatrix)

        CompartmentMatrix = RBA_Matrix()
        if VolumeFraction:
            A = numpy.eye(len(Compartments))*5/float(totalAA)
        else:
            A = numpy.eye(len(Compartments))/float(totalAA)
        CompartmentMatrix.A = scipy.sparse.coo_matrix(A)
        CompartmentMatrix.b = numpy.array([float(0)]*len(Compartments))
        CompartmentMatrix.f = numpy.array([float(0)]*len(Compartments))
        CompartmentMatrix.LB = numpy.array([float(0)]*len(Compartments))
        CompartmentMatrix.UB = numpy.array([float(totalAA)]*len(Compartments))
        CompartmentMatrix.row_signs = ['L']*(len(Compartments))
#        CompartmentMatrix.row_signs = ['E']*(len(Compartments))
        CompartmentMatrix.row_names = ['n_volume', 'mIM_volume', 'vM_volume', 'mIMS_volume',
                                       'm_volume', 'erM_volume', 'mOM_volume', 'x_volume', 'cM_volume', 'gM_volume', 'c_volume']
        CompartmentMatrix.col_names = ['O_n', 'O_mIM', 'O_vM', 'O_mIMS',
                                       'O_m', 'O_erM', 'O_mOM', 'O_x', 'O_cM', 'O_gM', 'O_c']
        CompartmentMatrix.mapIndices()
        self.Problem.LP.addMatrix(matrix=CompartmentMatrix)
        self.Problem.MuOneMatrix.addMatrix(matrix=CompartmentMatrix)

        VolumeMatrix = RBA_Matrix()
        A = numpy.ones((len(Compartments)+1, len(Compartments)))
        Eye = -numpy.eye(len(Compartments))
        A[0:len(Compartments), :] = Eye
        # A[len(Compartments), [1, 5, 6, 8, 9]] = 0
        # A[len(Compartments), 8] = 0
        VolumeMatrix.A = scipy.sparse.coo_matrix(A)
        VolumeMatrix.b = numpy.array([float(0)]*len(Compartments)+[float(1)])
        VolumeMatrix.f = numpy.array([float(0)]*len(Compartments))
        VolumeMatrix.LB = numpy.array([float(0)]*len(Compartments))
        VolumeMatrix.UB = numpy.array([float(1)]*len(Compartments))
        VolumeMatrix.row_signs = ['L']*(len(Compartments))+['E']
        # VolumeMatrix.row_signs = ['E']*(len(Compartments))+['E']
        VolumeMatrix.row_names = ['n_volume', 'mIM_volume', 'vM_volume', 'mIMS_volume', 'm_volume',
                                  'erM_volume', 'mOM_volume', 'x_volume', 'cM_volume', 'gM_volume', 'c_volume', 'TotalVolume']
        VolumeMatrix.col_names = ['F_n', 'F_mIM', 'F_vM', 'F_mIMS',
                                  'F_m', 'F_erM', 'F_mOM', 'F_x', 'F_cM', 'F_gM', 'F_c']

        if not CompartmentRelationships:
            VolumeMatrix.mapIndices()
            self.Problem.LP.addMatrix(matrix=VolumeMatrix)
            self.Problem.MuOneMatrix.addMatrix(matrix=VolumeMatrix)

        if CompartmentRelationships:
            Anew = numpy.zeros((A.shape[0]+3, A.shape[1]))
            Anew[0:A.shape[0], :] = A
            VolumeMatrix.row_names += ['m_mIM', 'm_mIMS', 'm_mOM']
            VolumeMatrix.row_signs += ['E', 'E', 'E']
            VolumeMatrix.b = numpy.array(list(VolumeMatrix.b)+[float(0)]*3)
            Anew[VolumeMatrix.row_names.index(
                'm_mIM'), VolumeMatrix.col_names.index('F_m')] = float(1)
            Anew[VolumeMatrix.row_names.index(
                'm_mIMS'), VolumeMatrix.col_names.index('F_m')] = float(1)
            Anew[VolumeMatrix.row_names.index(
                'm_mOM'), VolumeMatrix.col_names.index('F_m')] = float(1)
            Anew[VolumeMatrix.row_names.index(
                'm_mIM'), VolumeMatrix.col_names.index('F_mIM')] = -m_mIM
            Anew[VolumeMatrix.row_names.index(
                'm_mIMS'), VolumeMatrix.col_names.index('F_mIMS')] = -m_mIMS
            Anew[VolumeMatrix.row_names.index(
                'm_mOM'), VolumeMatrix.col_names.index('F_mOM')] = -m_mOM
            VolumeMatrix.A = scipy.sparse.coo_matrix(Anew)

            VolumeMatrix.mapIndices()
            self.Problem.LP.addMatrix(matrix=VolumeMatrix)
            self.Problem.MuOneMatrix.addMatrix(matrix=VolumeMatrix)

        if CompartmentComponents:

            PC_mIM = 0.0000883
            PE_mIM = 0.00005852
            PI_mIM = 0.00003377
            PS_mIM = 0.00000873
            CL_mIM = 0.00002
            PA_mIM = 0.0000039
            ES_mIM = 0.008547

            PC_mOM = 0.000636
            PE_mOM = 0.0004822
            PI_mOM = 0.0001289
            PS_mOM = 0.0000167
            CL_mOM = 0.00004467
            PA_mOM = 0.0000696
            ES_mOM = 0.0

            PC_vM = 0.0003635
            PE_vM = 0.4156
            PI_vM = 0.0001297
            PS_vM = 0.00003435
            CL_vM = 0.0000068
            PA_vM = 0.0000186
            ES_vM = 0.0142

            PC_n = 0.000055
            PE_n = 0.000035
            PI_n = 0.000017
            PS_n = 0.0000072
            CL_n = 0.0
            PA_n = 0.0000031
            ES_n = 0.0086

            PC_gM = 0.00043
            PE_gM = 0.00044
            PI_gM = 0.00041
            PS_gM = 0.0
            CL_gM = 0.00022
            PA_gM = 0.0
            ES_gM = 0.0

            PC_n = 0.0
            PE_n = 0.0
            PI_n = 0.0
            PS_n = 0.0
            CL_n = 0.0
            PA_n = 0.0
            ES_n = 0.0

            PC_gM = 0.0
            PE_gM = 0.0
            PI_gM = 0.0
            PS_gM = 0.0
            CL_gM = 0.0
            PA_gM = 0.0
            ES_gM = 0.0

            PC_vM = 0.0
            PE_vM = 0.0
            PI_vM = 0.0
            PS_vM = 0.0
            CL_vM = 0.0
            PA_vM = 0.0
            ES_vM = 0.0

            PC_mIM = 0.0
            PE_mIM = 0.0
            PI_mIM = 0.0
            PS_mIM = 0.0
            CL_mIM = 0.0
            PA_mIM = 0.0
            ES_mIM = 0.0

            PC_mOM = 0.0
            PE_mOM = 0.0
            PI_mOM = 0.0
            PS_mOM = 0.0
            CL_mOM = 0.0
            PA_mOM = 0.0
            ES_mOM = 0.0

            Alipids = RBA_Matrix()
            Alipids.col_names = ['F_mIM', 'F_mOM', 'F_vM', 'F_n', 'F_gM']
            Alipids.row_names = ['M_pc_SC_c', 'M_pe_SC_c', 'M_ptd1ino_SC_c',
                                 'M_ps_SC_c', 'M_clpn_SC_m', 'M_pa_SC_c', 'M_ergst_c']
            Alipids.row_signs = [
                self.Problem.LP.row_signs[self.Problem.LP.row_names.index(i)] for i in Alipids.row_names]
            Alipids.b = numpy.array(
                [self.Problem.LP.b[self.Problem.LP.row_names.index(i)] for i in Alipids.row_names])
            Alipids.LB = numpy.array([0, 0, 0, 0, 0])
            Alipids.UB = numpy.array([1, 1, 1, 1, 1])
            Alipids.f = numpy.array([0, 0, 0, 0, 0])
            LipidMatrix = numpy.zeros((7, 5))

            LipidMatrix[Alipids.row_names.index(
                'M_pc_SC_c'), Alipids.col_names.index('F_mIM')] = PC_mIM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_pe_SC_c'), Alipids.col_names.index('F_mIM')] = PE_mIM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ptd1ino_SC_c'), Alipids.col_names.index('F_mIM')] = PI_mIM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ps_SC_c'), Alipids.col_names.index('F_mIM')] = PS_mIM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_clpn_SC_m'), Alipids.col_names.index('F_mIM')] = CL_mIM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_pa_SC_c'), Alipids.col_names.index('F_mIM')] = PA_mIM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ergst_c'), Alipids.col_names.index('F_mIM')] = ES_mIM/totalAA

            LipidMatrix[Alipids.row_names.index(
                'M_pc_SC_c'), Alipids.col_names.index('F_mOM')] = PC_mOM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_pe_SC_c'), Alipids.col_names.index('F_mOM')] = PE_mOM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ptd1ino_SC_c'), Alipids.col_names.index('F_mOM')] = PI_mOM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ps_SC_c'), Alipids.col_names.index('F_mOM')] = PS_mOM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_clpn_SC_m'), Alipids.col_names.index('F_mOM')] = CL_mOM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_pa_SC_c'), Alipids.col_names.index('F_mOM')] = PA_mOM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ergst_c'), Alipids.col_names.index('F_mOM')] = ES_mOM/totalAA

            LipidMatrix[Alipids.row_names.index(
                'M_pc_SC_c'), Alipids.col_names.index('F_vM')] = PC_vM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_pe_SC_c'), Alipids.col_names.index('F_vM')] = PE_vM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ptd1ino_SC_c'), Alipids.col_names.index('F_vM')] = PI_vM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ps_SC_c'), Alipids.col_names.index('F_vM')] = PS_vM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_clpn_SC_m'), Alipids.col_names.index('F_vM')] = CL_vM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_pa_SC_c'), Alipids.col_names.index('F_vM')] = PA_vM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ergst_c'), Alipids.col_names.index('F_vM')] = ES_vM/totalAA

            LipidMatrix[Alipids.row_names.index(
                'M_pc_SC_c'), Alipids.col_names.index('F_n')] = PC_n/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_pe_SC_c'), Alipids.col_names.index('F_n')] = PE_n/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ptd1ino_SC_c'), Alipids.col_names.index('F_n')] = PI_n/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ps_SC_c'), Alipids.col_names.index('F_n')] = PS_n/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_clpn_SC_m'), Alipids.col_names.index('F_n')] = CL_n/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_pa_SC_c'), Alipids.col_names.index('F_n')] = PA_n/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ergst_c'), Alipids.col_names.index('F_n')] = ES_n/totalAA

            LipidMatrix[Alipids.row_names.index(
                'M_pc_SC_c'), Alipids.col_names.index('F_gM')] = PC_gM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_pe_SC_c'), Alipids.col_names.index('F_gM')] = PE_gM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ptd1ino_SC_c'), Alipids.col_names.index('F_gM')] = PI_gM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ps_SC_c'), Alipids.col_names.index('F_gM')] = PS_gM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_clpn_SC_m'), Alipids.col_names.index('F_gM')] = CL_gM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_pa_SC_c'), Alipids.col_names.index('F_gM')] = PA_gM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ergst_c'), Alipids.col_names.index('F_gM')] = ES_gM/totalAA

            MudepIndices = [('M_pc_SC_c', i) for i in Alipids.col_names]+[('M_pe_SC_c', i) for i in Alipids.col_names]+[('M_ptd1ino_SC_c', i) for i in Alipids.col_names]+[('M_ps_SC_c', i)
                                                                                                                                                                           for i in Alipids.col_names]+[('M_clpn_SC_m', i) for i in Alipids.col_names]+[('M_pa_SC_c', i) for i in Alipids.col_names]+[('M_ergst_c', i) for i in Alipids.col_names]
            self.Problem.MuDepIndices_A += MudepIndices
            Alipids.A = scipy.sparse.coo_matrix(LipidMatrix)
            Alipids.mapIndices()
            self.Problem.LP.updateMatrix(Alipids, Ainds=MudepIndices)
            self.Problem.MuOneMatrix.updateMatrix(Alipids, Ainds=MudepIndices)

    ## !!! ##
    def eukaryoticDensities4(self, CompartmentRelationships=True):
        Compartments = ['n', 'mIM', 'vM', 'mIMS', 'm', 'erM', 'mOM', 'x', 'c', 'cM', 'gM']
        totalAA = 3.1*0.91
        DensityIndices = [self.Problem.LP.row_names.index(
            i) for i in self.Problem.CompartmentDensities]

        A = self.Problem.LP.A.toarray()

        A0 = self.Problem.MuOneMatrix.A.toarray()

        OccupationMatrix = RBA_Matrix()
        A = numpy.ones((len(Compartments)+1, len(Compartments)))
        Eye = -numpy.eye(len(Compartments))
        A[0:len(Compartments), :] = Eye

        OccupationMatrix.b = numpy.array(list([float(0)]*len(Compartments))+[totalAA])
        OccupationMatrix.f = numpy.array([float(0)]*len(Compartments))
        OccupationMatrix.LB = numpy.array([float(0)]*len(Compartments))
        OccupationMatrix.UB = numpy.array([float(totalAA)]*len(Compartments))
        OccupationMatrix.row_signs = ['E']*(len(Compartments)+1)
        OccupationMatrix.row_names = ['n_density', 'mIM_density', 'vM_density', 'mIMS_density',
                                      'm_density', 'erM_density', 'mOM_density', 'x_density', 'cM_density', 'gM_density', 'c_density', 'O_total']
        OccupationMatrix.col_names = ['O_n', 'O_mIM', 'O_vM', 'O_mIMS',
                                      'O_m', 'O_erM', 'O_mOM', 'O_x', 'O_cM', 'O_gM', 'O_c']

        if CompartmentRelationships:
            m_mIM = 0.5
            m_mIMS = 1
            m_mOM = 5
            Anew = numpy.zeros((A.shape[0]+3, A.shape[1]))
            Anew[0:A.shape[0], :] = A
            OccupationMatrix.row_names += ['m_mIM', 'm_mIMS', 'm_mOM']
            OccupationMatrix.row_signs += ['E', 'E', 'E']
            OccupationMatrix.b = numpy.array(list(OccupationMatrix.b)+[float(0)]*3)
            Anew[OccupationMatrix.row_names.index(
                'm_mIM'), OccupationMatrix.col_names.index('O_m')] = float(1)
            Anew[OccupationMatrix.row_names.index(
                'm_mIMS'), OccupationMatrix.col_names.index('O_m')] = float(1)
            Anew[OccupationMatrix.row_names.index(
                'm_mOM'), OccupationMatrix.col_names.index('O_m')] = float(1)
            Anew[OccupationMatrix.row_names.index(
                'm_mIM'), OccupationMatrix.col_names.index('O_mIM')] = -m_mIM
            Anew[OccupationMatrix.row_names.index(
                'm_mIMS'), OccupationMatrix.col_names.index('O_mIMS')] = -m_mIMS
            Anew[OccupationMatrix.row_names.index(
                'm_mOM'), OccupationMatrix.col_names.index('O_mOM')] = -m_mOM

            OccupationMatrix.A = scipy.sparse.coo_matrix(Anew)
        else:
            OccupationMatrix.A = scipy.sparse.coo_matrix(A)

        OccupationMatrix.mapIndices()
        self.Problem.LP.addMatrix(matrix=OccupationMatrix)
        self.Problem.MuOneMatrix.addMatrix(matrix=OccupationMatrix)
        # {'Index':{'Param1':'+','Param2':'+','Param2':'-'}}

        #Type: 'Sum'#
        # {'Index':'Param1'}
        self.Problem.MuDependencies['FromParameters']['b'].update(
            {'n_density': 'AAres_PG_nucleus_DNA'})
        self.Problem.MuDependencies['FromParameters']['b'].update(
            {'O_total': {'Equation': 'amino_acid_concentration_total - AAres_PG_secreted_Euk', 'Variables': ['amino_acid_concentration_total', 'AAres_PG_secreted_Euk']}})
        self.Problem.MuDependencies['FromMatrix']['b'].remove('n_density')
        self.Problem.MuDependencies['FromMatrix']['b'].remove('mIM_density')
        self.Problem.MuDependencies['FromMatrix']['b'].remove('vM_density')
        self.Problem.MuDependencies['FromMatrix']['b'].remove('mIMS_density')
        self.Problem.MuDependencies['FromMatrix']['b'].remove('m_density')
        self.Problem.MuDependencies['FromMatrix']['b'].remove('erM_density')
        self.Problem.MuDependencies['FromMatrix']['b'].remove('mOM_density')
        self.Problem.MuDependencies['FromMatrix']['b'].remove('x_density')
        self.Problem.MuDependencies['FromMatrix']['b'].remove('cM_density')
        self.Problem.MuDependencies['FromMatrix']['b'].remove('gM_density')
        self.Problem.MuDependencies['FromMatrix']['b'].remove('c_density')

    ## !!! ##
    def eukaryoticDensities_calibration(self, CompartmentRelationships=False, mitoProportions={}, amino_acid_concentration_total='amino_acid_concentration_total'):
        Compartments = ['n', 'mIM', 'vM', 'mIMS', 'm', 'erM', 'mOM', 'x', 'c', 'cM', 'gM']
        totalAA_parameter = amino_acid_concentration_total
        totalAA = 3.1*0.91
        DensityIndices = [self.Problem.LP.row_names.index(
            i) for i in self.Problem.CompartmentDensities]

        A = self.Problem.LP.A.toarray()

        A0 = self.Problem.MuOneMatrix.A.toarray()

        OccupationMatrix = RBA_Matrix()
        A = numpy.ones((len(Compartments)+1, len(Compartments)))
        Eye = -numpy.eye(len(Compartments))
        A[0:len(Compartments), :] = Eye

        OccupationMatrix.b = numpy.array(list([float(0)]*len(Compartments))+[totalAA])
        OccupationMatrix.f = numpy.array([float(0)]*len(Compartments))
        OccupationMatrix.LB = numpy.array([float(0)]*len(Compartments))
        OccupationMatrix.UB = numpy.array([float(totalAA)]*len(Compartments))
        OccupationMatrix.row_signs = ['E']*(len(Compartments)+1)
        OccupationMatrix.row_names = ['n_density', 'mIM_density', 'vM_density', 'mIMS_density',
                                      'm_density', 'erM_density', 'mOM_density', 'x_density', 'cM_density', 'gM_density', 'c_density', 'O_total']
        OccupationMatrix.col_names = ['O_n', 'O_mIM', 'O_vM', 'O_mIMS',
                                      'O_m', 'O_erM', 'O_mOM', 'O_x', 'O_cM', 'O_gM', 'O_c']

        if CompartmentRelationships:
            if len(list(mitoProportions.keys())) == 3:
                m_mIM = mitoProportions['m_mIM']
                m_mIMS = mitoProportions['m_mIMS']
                m_mOM = mitoProportions['m_mOM']
                Anew = numpy.zeros((A.shape[0]+3, A.shape[1]))
                Anew[0:A.shape[0], :] = A
                OccupationMatrix.row_names += ['m_mIM', 'm_mIMS', 'm_mOM']
                OccupationMatrix.row_signs += ['E', 'E', 'E']
                OccupationMatrix.b = numpy.array(list(OccupationMatrix.b)+[float(0)]*3)
                Anew[OccupationMatrix.row_names.index(
                    'm_mIM'), OccupationMatrix.col_names.index('O_m')] = float(1)
                Anew[OccupationMatrix.row_names.index(
                    'm_mIMS'), OccupationMatrix.col_names.index('O_m')] = float(1)
                Anew[OccupationMatrix.row_names.index(
                    'm_mOM'), OccupationMatrix.col_names.index('O_m')] = float(1)
                Anew[OccupationMatrix.row_names.index(
                    'm_mIM'), OccupationMatrix.col_names.index('O_mIM')] = -m_mIM
                Anew[OccupationMatrix.row_names.index(
                    'm_mIMS'), OccupationMatrix.col_names.index('O_mIMS')] = -m_mIMS
                Anew[OccupationMatrix.row_names.index(
                    'm_mOM'), OccupationMatrix.col_names.index('O_mOM')] = -m_mOM

                OccupationMatrix.A = scipy.sparse.coo_matrix(Anew)
        else:
            OccupationMatrix.A = scipy.sparse.coo_matrix(A)

        OccupationMatrix.mapIndices()
        self.Problem.LP.addMatrix(matrix=OccupationMatrix)
        self.Problem.MuOneMatrix.addMatrix(matrix=OccupationMatrix)
        # {'Index':{'Param1':'+','Param2':'+','Param2':'-'}}

        #Type: 'Sum'#
        # {'Index':'Param1'}
        self.Problem.MuDependencies['FromParameters']['b'].update(
            {'n_density': {'Equation': '-nonenzymatic_proteins_n/inverse_average_protein_length', 'Variables': ['nonenzymatic_proteins_n', 'inverse_average_protein_length']}})
        self.Problem.MuDependencies['FromParameters']['b'].update({'mIM_density': {
                                                                  'Equation': '-nonenzymatic_proteins_mIM/inverse_average_protein_length', 'Variables': ['nonenzymatic_proteins_mIM', 'inverse_average_protein_length']}})
        self.Problem.MuDependencies['FromParameters']['b'].update({'vM_density': {
                                                                  'Equation': '-nonenzymatic_proteins_vM/inverse_average_protein_length', 'Variables': ['nonenzymatic_proteins_vM', 'inverse_average_protein_length']}})
        self.Problem.MuDependencies['FromParameters']['b'].update({'mIMS_density': {
                                                                  'Equation': '-nonenzymatic_proteins_mIMS/inverse_average_protein_length', 'Variables': ['nonenzymatic_proteins_mIMS', 'inverse_average_protein_length']}})
        self.Problem.MuDependencies['FromParameters']['b'].update(
            {'m_density': {'Equation': '-nonenzymatic_proteins_m/inverse_average_protein_length', 'Variables': ['nonenzymatic_proteins_m', 'inverse_average_protein_length']}})
        self.Problem.MuDependencies['FromParameters']['b'].update({'erM_density': {
                                                                  'Equation': '-nonenzymatic_proteins_erM/inverse_average_protein_length', 'Variables': ['nonenzymatic_proteins_erM', 'inverse_average_protein_length']}})
        self.Problem.MuDependencies['FromParameters']['b'].update({'mOM_density': {
                                                                  'Equation': '-nonenzymatic_proteins_mOM/inverse_average_protein_length', 'Variables': ['nonenzymatic_proteins_mOM', 'inverse_average_protein_length']}})
        self.Problem.MuDependencies['FromParameters']['b'].update(
            {'x_density': {'Equation': '-nonenzymatic_proteins_x/inverse_average_protein_length', 'Variables': ['nonenzymatic_proteins_x', 'inverse_average_protein_length']}})
        self.Problem.MuDependencies['FromParameters']['b'].update({'cM_density': {
                                                                  'Equation': '-nonenzymatic_proteins_cM/inverse_average_protein_length', 'Variables': ['nonenzymatic_proteins_cM', 'inverse_average_protein_length']}})
        self.Problem.MuDependencies['FromParameters']['b'].update({'gM_density': {
                                                                  'Equation': '-nonenzymatic_proteins_gM/inverse_average_protein_length', 'Variables': ['nonenzymatic_proteins_gM', 'inverse_average_protein_length']}})
        self.Problem.MuDependencies['FromParameters']['b'].update(
            {'c_density': {'Equation': '-nonenzymatic_proteins_c/inverse_average_protein_length', 'Variables': ['nonenzymatic_proteins_c', 'inverse_average_protein_length']}})
        self.Problem.MuDependencies['FromParameters']['b'].update({'O_total': {'Equation': '{} - nonenzymatic_proteins_Secreted/inverse_average_protein_length'.format(totalAA_parameter), 'Variables': [
                                                                  totalAA_parameter, 'nonenzymatic_proteins_Secreted', 'inverse_average_protein_length']}})

    # !!! deal with hardcoded parameter_names... !!!
    def estimate_specific_Kapps(self, proteomicsData, flux_bounds, mu, biomass_function=None, target_biomass_function=True, parsimonious_fba=True):
        """
        Parameters
        ----------
        proteomicsData : pandas.DataFrame (in mmol/gDW)
        flux_bounds : pandas.DataFrame  (in mmol/(gDW*h))
        mu : float (in 1/h)
        biomass_function : str
        target_biomass_function : bool
        atp_maintenance_to_biomassfunction : bool
        eukaryotic : bool
        """
        from scipy.stats.mstats import gmean

        old_model = copy.deepcopy(self.model)
        for i in self.model.targets.target_groups._elements_by_id['translation_targets'].concentrations._elements:
            if i.species == 'average_protein_c':
                new_agg = rba.xml.parameters.Aggregate(id_='total_protein', type_='multiplication')
                new_agg.function_references.append(rba.xml.parameters.FunctionReference(
                    function='amino_acid_concentration_total'))
                new_agg.function_references.append(rba.xml.parameters.FunctionReference(
                    function='inverse_average_protein_length'))
                self.model.parameters.aggregates._elements.append(new_agg)
                i.value = 'total_protein'
            else:
                self.model.targets.target_groups._elements_by_id['translation_targets'].concentrations._elements.remove(
                    i)
        for i in self.model.targets.target_groups._elements_by_id['transcription_targets'].concentrations._elements:
            if i.species == 'mrna':
                new_agg = rba.xml.parameters.Aggregate(id_='total_rna', type_='multiplication')
                new_agg.function_references.append(rba.xml.parameters.FunctionReference(
                    function='RNA_massfraction_CarbonLimitation'))
                new_agg.function_references.append(
                    rba.xml.parameters.FunctionReference(function='RNA_inversemillimolarweight'))
                self.model.parameters.aggregates._elements.append(new_agg)
                i.value = 'total_rna'
            else:
                self.model.targets.target_groups._elements_by_id['transcription_targets'].concentrations._elements.remove(
                    i)

        self.rebuild_from_model()
        self.setMedium(self.Medium)

        self.addExchangeReactions()
        self.setMu(mu)

        if target_biomass_function:
            self.buildFBA(objective='targets', maintenanceToBM=True)
            BMfunction = 'R_BIOMASS_targetsRBA'
        else:
            self.buildFBA(objective='classic', maintenanceToBM=False)
            BMfunction = biomass_function

        for j in [i for i in self.Medium.keys() if self.Medium[i] == 0]:
            Exrxn = 'R_EX_'+j.split('M_')[-1]+'_e'
            self.FBA.setUB({Exrxn: 0})

        rxn_LBs = {}
        rxn_UBs = {}
        for rx in flux_bounds['Reaction_ID']:
            lb = flux_bounds.loc[flux_bounds['Reaction_ID'] == rx, 'LB'].values[0]
            ub = flux_bounds.loc[flux_bounds['Reaction_ID'] == rx, 'UB'].values[0]
            if not pandas.isna(lb):
                rxn_LBs.update({rx: lb})
            if not pandas.isna(ub):
                rxn_UBs.update({rx: ub})
        self.FBA.setLB(rxn_LBs)
        self.FBA.setUB(rxn_UBs)

        self.FBA.clearObjective()
        self.FBA.setObjectiveCoefficients({BMfunction: -1})
        self.FBA.solveLP(feasibleStatuses=[1, 2, 3, 5, 6])
        BMfluxOld = self.FBA.SolutionValues[BMfunction]

        if parsimonious_fba:
            self.FBA.parsimonise()
            self.FBA.setLB(rxn_LBs)
            self.FBA.setUB(rxn_UBs)
            self.FBA.setLB({BMfunction: BMfluxOld})
            self.FBA.setUB({BMfunction: BMfluxOld})
            self.FBA.solveLP(feasibleStatuses=[1, 2, 3, 5, 6])

        FluxDistribution = pandas.DataFrame(index=list(
            self.FBA.SolutionValues.keys()), columns=['FluxValues'])
        FluxDistribution['FluxValues'] = list(self.FBA.SolutionValues.values())
        BMfluxNew = self.FBA.SolutionValues[BMfunction]

        ProtoIDmap = {}
        for i in self.ModelStructure.ProteinInfo.Elements.keys():
            ProtoID = self.ModelStructure.ProteinInfo.Elements[i]['ProtoID']
            if ProtoID in list(proteomicsData['ID']):
                if not pandas.isna(proteomicsData.loc[proteomicsData['ID'] == ProtoID, 'copy_number'].values[0]):
                    if proteomicsData.loc[proteomicsData['ID'] == ProtoID, 'copy_number'].values[0] != 0:
                        if ProtoID in ProtoIDmap.keys():
                            ProtoIDmap[ProtoID]['ModelProteins'].append(i)
                        else:
                            ProtoIDmap.update(
                                {ProtoID: {'ModelProteins': [i], 'CopyNumber': proteomicsData.loc[proteomicsData['ID'] == ProtoID, 'copy_number'].values[0]}})

        ReactionMap = {}
        for i in self.ModelStructure.ReactionInfo.Elements.keys():
            if '_duplicate_' in i:
                continue
            else:
                if i in list(FluxDistribution.index):
                    if FluxDistribution.loc[i, 'FluxValues'] != 0:
                        ReactionMap.update({i: {'ModelReactions': list(
                            [i]+self.ModelStructure.ReactionInfo.Elements[i]['Twins']), 'Flux': FluxDistribution.loc[i, 'FluxValues']}})

        IsoReaction2ProtoReaction = {}
        for i in ReactionMap.keys():
            for j in ReactionMap[i]['ModelReactions']:
                IsoReaction2ProtoReaction[j] = i

        EnzymeMap = {}
        for i in self.ModelStructure.EnzymeInfo.Elements.keys():
            if self.ModelStructure.EnzymeInfo.Elements[i]['Reaction'] in IsoReaction2ProtoReaction:
                CompositionDict = {self.ModelStructure.ProteinInfo.Elements[j]['ProtoID']: self.ModelStructure.EnzymeInfo.Elements[
                    i]['Subunits'][j] for j in self.ModelStructure.EnzymeInfo.Elements[i]['Subunits'].keys()}
                ProtoReaction = IsoReaction2ProtoReaction[self.ModelStructure.EnzymeInfo.Elements[i]['Reaction']]
                CopyNumbers = []
                Stoichiometries = []
                EnzymeNumbers = []
                for j in CompositionDict.keys():
                    if j in ProtoIDmap.keys():
                        CopyNumbers.append(ProtoIDmap[j]['CopyNumber'])
                        Stoichiometries.append(CompositionDict[j])
                        EnzymeNumbers.append(ProtoIDmap[j]['CopyNumber']/CompositionDict[j])
                GM_enzymenumber = 0
                if len(EnzymeNumbers) > 0:
                    GM_enzymenumber = gmean(numpy.array(EnzymeNumbers))
                EnzymeMap.update(
                    {i: {'ProtoReaction': ProtoReaction, 'EnzymeNumber': GM_enzymenumber}})

        EnzymeMap2 = {}
        for i in ReactionMap.keys():
            totalIsoEnzymeNumber = 0
            for j in ReactionMap[i]['ModelReactions']:
                respectiveEnzyme = self.ModelStructure.ReactionInfo.Elements[j]['Enzyme']
                if respectiveEnzyme in EnzymeMap.keys():
                    totalIsoEnzymeNumber += EnzymeMap[respectiveEnzyme]['EnzymeNumber']
            for j in ReactionMap[i]['ModelReactions']:
                respectiveEnzyme = self.ModelStructure.ReactionInfo.Elements[j]['Enzyme']
                if respectiveEnzyme in EnzymeMap.keys():
                    concentration = EnzymeMap[respectiveEnzyme]['EnzymeNumber']
                    if concentration != 0:
                        if numpy.isfinite(concentration):
                            specificFlux = ReactionMap[i]['Flux'] * \
                                EnzymeMap[respectiveEnzyme]['EnzymeNumber']/totalIsoEnzymeNumber
                            EnzymeMap2.update({respectiveEnzyme: {'CopyNumber': EnzymeMap[respectiveEnzyme]['EnzymeNumber'],
                                                                  'Concentration': concentration, 'Flux': specificFlux, 'Kapp': abs(specificFlux/concentration)}})

        self.model = old_model
        self.rebuild_from_model()
        self.setMedium(self.Medium)

        out = pandas.DataFrame()
        for i in EnzymeMap2.keys():
            # if EnzymeMap2[i]['CopyNumber'] == 0:
            #    continue
            out.loc[i, 'Enzyme_ID'] = i
            out.loc[i, 'CopyNumber'] = EnzymeMap2[i]['CopyNumber']
            out.loc[i, 'Concentration'] = EnzymeMap2[i]['Concentration']
            out.loc[i, 'Flux'] = EnzymeMap2[i]['Flux']
            out.loc[i, 'Kapp'] = EnzymeMap2[i]['Kapp']

        return(out)

    def estimate_default_Kapps(self, target_mu, compartment_densities_and_PGs=None, flux_bounds=None, plateau_limit=4, mu_approximation_precision=0.0001, transporter_to_lumen_coefficient=10, default_kapp_LB=0, default_kapp_UB=1000000, start_val=None, densities_to_fix=None, eukaryotic=False):
        """
        Parameters
        ----------
        target_mu : float
        compartment_densities_and_PGs : pandas.DataFrame
        flux_bounds : pandas.DataFrame
        """
        old_model = copy.deepcopy(self.model)

        orig_enz = self.model.parameters.functions._elements_by_id[
            'default_efficiency'].parameters._elements_by_id['CONSTANT'].value

        out = pandas.DataFrame()
        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
            self.model.parameters.functions._elements_by_id[str(
                'fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density']
            self.model.parameters.functions._elements_by_id[str(
                'fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction']
        self.rebuild_from_model()
        self.addExchangeReactions()
        self.setMedium(self.Medium)
        if eukaryotic:
            self.eukaryoticDensities_calibration(CompartmentRelationships=False)
        else:
            if densities_to_fix is None:
                comp_density_rows = self.Problem.CompartmentDensities
            else:
                comp_density_rows = densities_to_fix
            self.Problem.setConstraintType(
                dict(zip(comp_density_rows, ['E']*len(comp_density_rows))))
        # self.Problem.setConstraintType(
        #    dict(zip(['mOM_density', 'mIMS_density', 'mIM_density', 'm_density'], ['L'])))

        rxn_LBs = {}
        rxn_UBs = {}
        for rx in flux_bounds['Reaction_ID']:
            lb = flux_bounds.loc[flux_bounds['Reaction_ID'] == rx, 'LB'].values[0]
            ub = flux_bounds.loc[flux_bounds['Reaction_ID'] == rx, 'UB'].values[0]
            if not pandas.isna(lb):
                rxn_LBs.update({rx: lb})
            if not pandas.isna(ub):
                rxn_UBs.update({rx: ub})
        self.Problem.setLB(rxn_LBs)
        self.Problem.setUB(rxn_UBs)

        kapp_LB = default_kapp_LB
        if default_kapp_UB is not None:
            kapp_UB = default_kapp_UB
        else:
            kapp_UB = orig_enz*1000
        # new_kapp = (kapp_UB+kapp_LB)/2

        if start_val is not None:
            new_kapp = start_val
        else:
            new_kapp = orig_enz

        Mu_pred = self.findMaxGrowthRate(precision=0.0001, max=1)

        Mus = []
        Mus_Error = []
        Kapps = []
        last_Mu = numpy.nan
        plateau_count = 0
        if abs(target_mu - Mu_pred) > mu_approximation_precision:
            while abs(target_mu - Mu_pred) > mu_approximation_precision:
                if plateau_count >= plateau_limit:
                    break
                self.model.parameters.functions._elements_by_id[
                    'default_efficiency'].parameters._elements_by_id['CONSTANT'].value = new_kapp
                self.model.parameters.functions._elements_by_id['default_transporter_efficiency'].parameters._elements_by_id[
                    'CONSTANT'].value = transporter_to_lumen_coefficient*new_kapp
                self.rebuild_from_model()
                self.addExchangeReactions()
                self.setMedium(self.Medium)
                self.Problem.setLB(rxn_LBs)
                self.Problem.setUB(rxn_UBs)
                if eukaryotic:
                    self.eukaryoticDensities_calibration(CompartmentRelationships=False)
                else:
                    self.Problem.setConstraintType(
                        dict(zip(comp_density_rows, ['E']*len(comp_density_rows))))

                Mu_pred = self.findMaxGrowthRate()
                Mus_Error.append(abs(target_mu - Mu_pred))
                Mus.append(Mu_pred)
                Kapps.append(new_kapp)
                if Mu_pred > target_mu:
                    new_kapp_prelim = kapp_LB+(0.5*abs(kapp_LB-new_kapp))
                    kapp_UB = new_kapp
                elif Mu_pred < target_mu:
                    new_kapp_prelim = kapp_UB-(0.5*abs(new_kapp-kapp_UB))
                    kapp_LB = new_kapp
                new_kapp = new_kapp_prelim
                if len(Mus) > 2:
                    if Mus[-2] == Mu_pred:
                        plateau_count += 1
                    else:
                        plateau_count = 0
        else:
            Mus.append(Mu_pred)
            Mus_Error.append(abs(target_mu - Mu_pred))
            Kapps.append(
                self.model.parameters.functions._elements_by_id['default_efficiency'].parameters._elements_by_id['CONSTANT'].value)
        # self.model = old_model
        self.rebuild_from_model()
        self.setMedium(self.Medium)
        out = pandas.DataFrame()
        out['Mu'] = Mus
        out['delta_Mu'] = Mus_Error
        out['default_efficiency'] = Kapps
        out['default_transporter_efficiency'] = [transporter_to_lumen_coefficient*i for i in Kapps]
        return(out)

    def inject_default_kapps(self, default_kapp, default_transporter_kapp):
        if numpy.isfinite(default_kapp):
            self.model.parameters.functions._elements_by_id[
                'default_efficiency'].parameters._elements_by_id['CONSTANT'].value = default_kapp
        if numpy.isfinite(default_transporter_kapp):
            self.model.parameters.functions._elements_by_id[
                'default_transporter_efficiency'].parameters._elements_by_id['CONSTANT'].value = default_transporter_kapp
        self.rebuild_from_model()

    def inject_process_capacities(self, process_efficiencies):
        """
        Parameters
        ----------
        process_efficiencies : pandas.DataFrame(columns=['Process','Parameter','Value'])
        """
        for i in process_efficiencies.index:
            if numpy.isfinite(process_efficiencies.loc[i, 'Value']):
                if process_efficiencies.loc[i, 'Process'] in self.model.processes.processes._elements_by_id.keys():
                    if not pandas.isna(process_efficiencies.loc[i, 'Value']):
                        self.model.processes.processes._elements_by_id[process_efficiencies.loc[i,
                                                                                                'Process']].machinery.capacity.value = process_efficiencies.loc[i, 'Parameter']
                        const = rba.xml.parameters.Function(process_efficiencies.loc[i, 'Parameter'], 'constant', parameters={
                                                            'CONSTANT': process_efficiencies.loc[i, 'Value']}, variable=None)
                        if process_efficiencies.loc[i, 'Parameter'] not in self.model.parameters.functions._elements_by_id.keys():
                            self.model.parameters.functions.append(const)
                        else:
                            self.model.parameters.functions._elements_by_id[const.id].parameters._elements_by_id[
                                'CONSTANT'].value = process_efficiencies.loc[i, 'Value']
        self.rebuild_from_model()

    def inject_specific_kapps(self, specific_kapps, round_to_digits=0):
        """
        Parameters
        ----------
        specific_kapps : pandas.DataFrame
        """
        parameterized = []
        if 'Enzyme_ID' in list(specific_kapps.columns):
            for enz in list(specific_kapps['Enzyme_ID']):
                if not pandas.isna(specific_kapps.loc[specific_kapps['Enzyme_ID'] == enz, 'Kapp'].values[0]):
                    if numpy.isfinite(specific_kapps.loc[specific_kapps['Enzyme_ID'] == enz, 'Kapp'].values[0]):
                        if enz not in parameterized:
                            all_enzs = self.ModelStructure.EnzymeInfo.Elements[enz]['Isozymes']
                            all_enzs.append(enz)
                            parameterized += all_enzs
                            if len(all_enzs) == 1:
                                proto_enz = all_enzs[0]
                            else:
                                proto_enz = [i for i in all_enzs if not '_duplicate_' in i][0]
                            val = round(specific_kapps.loc[specific_kapps['Enzyme_ID']
                                                           == enz, 'Kapp'].values[0], round_to_digits)
                            const = rba.xml.parameters.Function(
                                str(proto_enz + '_kapp__constant'), 'constant', parameters={'CONSTANT': val}, variable=None)
                            if str(proto_enz + '_kapp__constant') not in self.model.parameters.functions._elements_by_id.keys():
                                self.model.parameters.functions.append(const)
                            else:
                                # self.model.parameters.functions._elements_by_id[const.id] = const
                                self.model.parameters.functions._elements_by_id[
                                    const.id].parameters._elements_by_id['CONSTANT'].value = val
                            count = 0
                            # self.model.parameters.functions._elements_by_id['default_efficiency'].parameters._elements_by_id['CONSTANT'].value = default_kapp
                            for e in self.model.enzymes.enzymes:
                                if e.id in all_enzs:
                                    count += 1
                                    e.forward_efficiency = str(proto_enz + '_kapp__constant')
                                    e.backward_efficiency = str(proto_enz + '_kapp__constant')
                                    if count == len(all_enzs):
                                        break
            self.rebuild_from_model()

    def findMinMediumConcentration(self, metabolite, precision=0.00001, max=100, recording=False, loggingIntermediateSteps=False):
        """
        Applies dichotomy-search to find the minimal feasible concentration of
        growth-substrate in medium, at a previously set growth-rate.
        Parameters
        ----------
        metabolite : str
            ID of metabolite in medium.
        precision : float
            Numberic precision with which minimum is approximated.
            Default : 0.00001
        max : float
            Defines the highest concentration rate to be screened for.
            Default=100
        recording : bool
            Records intermediate feasible solutions
            while approaching the minimum concentration.
            Default : False
        Returns
        -------
        minimum feasible growth-substrate concentration as float.
        """

        minConc = 0.0
        maxConc = max
        testConc = minConc
        iteration = 0
        oldConc = self.Medium[metabolite]
        while (maxConc - minConc) > precision:
            self.setMedium(changes={metabolite: testConc})
            self.Problem.solveLP(logging=loggingIntermediateSteps)
            if self.Problem.Solved:
                iteration += 1
                if recording:
                    run_name = 'Dichotomy_'+metabolite+'_' + \
                        str(testConc)+'_iteration_'+str(iteration)
                    self.recordResults(run_name)
                maxConc = testConc
            else:
                minConc = testConc
            testConc = numpy.mean([maxConc, minConc])
        self.LogBook.addEntry(
            'Minimal required {} concentration found to be: {}.'.format(metabolite, maxConc))
        self.setMedium(changes={metabolite: oldConc})
        return(maxConc)

    def addProtein(self, input):
        """
        Adds representation of individual proteins to problem.
        Parameters
        ----------
        input : dict or str
            If input is str it has to be the ID of a protein in the model.
            Then this protein is added to the problem an creates:
                One constraint named Protein_'ID' (equality).
                One variable named TotalLevel_'ID' representing the total amount.
                One variable named Free_'ID'_'respectiveCompartment', this
                represents the fraction of the protein not assuming any function.
                It however consumes resources for synthesis (precursors and processes),
                which are the same as defined in the model files.
                And takes up space i the compartment as specified in the model-files
                for the protein.
            If input is dict it has to have two keys; 'ID' and 'UnusedProteinFraction'.
            By specifying this input one can define that the unused franction of the protein
            can also reside in other compartments and which processes it requires.
            The value to 'ID' is the ID of a protein in the model.
            The value to 'UnusedProteinFraction' is another dictionary.
            This can have several keys which must be model-compartments.
            For each of the keys the value is a dict holding IDs of model-processes as Keys
            and process requirements as Values (numerical).
            This specifies which processes each of the compartment-species of the protein
            requires.
            This generates the same constraint and TotalLevel-variable as with the simple input,
            however a variable representing each of the compartment-species for the unused fraction
            is added and incorporates the specific process requirements.
            E.g: input = {'ID': 'proteinA',
                          'UnusedProteinFraction':{'Cytoplasm':{'Translation':100}, {'Folding':10}],
                                                   'Membrane':{'Translation':100}, {'Folding':20}, {'Secretion':100}
                                                   }
                          }
                This adds 'proteinA' to the model, where the unused fraction can reside either in
                the Cytoplasm or the Membrane. However while the cytosolic-species only requires the
                processes 'Translation' and 'Folding'; the membrane-bound species also requires 'Secretion'
                and occupies more folding capacity.
                Then the constraint 'Protein_proteinA' is added and the 3 variables
                'TotalLevel_proteinA', 'Free_proteinA_Cytoplasm' and 'Free_proteinA_Membrane'.
        """

        if type(input) is str:
            input = {'ID': input}

        if 'ID' not in list(input.keys()):
            print('Error, no protein ID provided')
            return

        if input['ID'] not in list(self.ModelStructure.ProteinInfo.Elements.keys()):
            print('Error, protein not in model')
            return

        if 'UnusedProteinFraction' not in list(input.keys()):
            input.update({'UnusedProteinFraction':
                          {self.ModelStructure.ProteinInfo.Elements[input['ID']]['Compartment']:
                           self.ModelStructure.ProteinInfo.Elements[input['ID']]['ProcessRequirements']}})

        self.LogBook.addEntry('Protein {} added with specifications {}.'.format(
            input['ID'], str(json.dumps(input))))

        Muindexlist = []

        ## Building RBA_Matrix-object for new constraint-row, representing protein ##
        UsedProtein = RBA_Matrix()
        UsedProtein.A = scipy.sparse.coo_matrix(
            buildUsedProteinConstraint(Controler=self, protein=input['ID']))
        UsedProtein.b = numpy.array([float(0)])
        UsedProtein.f = numpy.array(self.Problem.LP.f)
        UsedProtein.LB = numpy.array(self.Problem.LP.LB)
        UsedProtein.UB = numpy.array(self.Problem.LP.UB)
        UsedProtein.row_signs = ['E']
        UsedProtein.row_names = ['Protein_'+input['ID']]
        UsedProtein.col_names = self.Problem.LP.col_names
        ## Add used protein row to problem ##
        self.Problem.LP.addMatrix(matrix=UsedProtein)
        ## Add used protein row to reference Matrix (Mu == 1) ##
        self.Problem.MuOneMatrix.addMatrix(matrix=UsedProtein)

        ## Building RBA_Matrix-object for new variable-col, representing total level of protein ##
        TotProtein = RBA_Matrix()
        TotProtein.A = scipy.sparse.coo_matrix(numpy.array(numpy.matrix(
            numpy.array([float(0)]*self.Problem.LP.A.shape[0]+[float(-1)])).transpose()))
        TotProtein.f = numpy.array([float(0)])
        TotProtein.LB = numpy.array([float(0)])
        TotProtein.UB = numpy.array([float(100000.0)])
        TotProtein.b = numpy.array(list(self.Problem.LP.b)+list(UsedProtein.b))
        TotProtein.row_signs = self.Problem.LP.row_signs+UsedProtein.row_signs
        TotProtein.row_names = self.Problem.LP.row_names+UsedProtein.row_names
        TotProtein.col_names = ['TotalLevel_'+input['ID']]
        ## Add total protein col to problem ##
        self.Problem.LP.addMatrix(matrix=TotProtein)
        ## Add total protein col to reference Matrix (Mu == 1) ##
        self.Problem.MuOneMatrix.addMatrix(matrix=TotProtein)

        ## Building RBA_Matrix-object for new variable-col,##
        ## representing each compartment-species of the protein ##
        for comp_species in list(input['UnusedProteinFraction'].keys()):
            ## Initiate RBA_Matrix object##
            UnusedProtein = RBA_Matrix()
            UnusedProtein.col_names = ['Free_'+input['ID']+'_'+comp_species]

            ## Extract required processes for protein and the respective demand ##
            ProcIDs = list(input['UnusedProteinFraction'][comp_species].keys())
            Preq = list(input['UnusedProteinFraction'][comp_species].values())
            ProcessCost = dict(
                zip([self.ModelStructure.ProcessInfo.Elements[k]['ID'] for k in ProcIDs], Preq))

            ## Get required charged trna buildingblocks and their stoichiometry in protein ##
            composition = self.ModelStructure.ProteinInfo.Elements[input['ID']]['AAcomposition']
            ## Extract the composition of charged trnas in terms of metabolic species ##
            species = self.ModelStructure.ProcessInfo.Elements['Translation']['Components']
            ## Determine required metabolites and their stoichiometry in protein ##
            MetaboliteCost = buildCompositionofUnusedProtein(
                species=species, composition=composition)

            ## Assemble process and metabolite requirements into stoichiometric coloumn vector ##
            ## And add to RBA_Matrix object ##
            colToAdd = numpy.array(numpy.matrix(numpy.array(list(MetaboliteCost.values())+list(ProcessCost.values()) +
                                                            [float(1)]+[self.ModelStructure.ProteinInfo.Elements[input['ID']]['AAnumber']])).transpose())
            UnusedProtein.A = scipy.sparse.coo_matrix(colToAdd)
            ## Add other information to RBA_Matrix object ##
            UnusedProtein.row_names = list(MetaboliteCost.keys())+[str(pc+'_capacity') for pc in list(
                ProcessCost.keys())]+['Protein_'+input['ID']]+[str(comp_species + '_density')]
            UnusedProtein.b = numpy.zeros(len(UnusedProtein.row_names))
            UnusedProtein.row_signs = ['E']*len(UnusedProtein.row_names)
            UnusedProtein.LB = numpy.array([float(0)])
            UnusedProtein.UB = numpy.array([float(100000.0)])
            UnusedProtein.f = numpy.array([float(0)])
            self.ProteinDilutionIndices = list(
                zip(list(MetaboliteCost.keys()), UnusedProtein.col_names*len(list(MetaboliteCost.keys()))))
            ## Add free protein col to problem ##
            self.Problem.LP.addMatrix(matrix=UnusedProtein)
            ## Add free protein col to reference Matrix (Mu == 1) ##
            self.Problem.MuOneMatrix.addMatrix(matrix=UnusedProtein)

            ## Find coefficients of unused protein column, subject to dilution (Metabolite and Process cost) ##
            ## And add them to MuDepIndices_A ##
            nonZeroEntries = numpy.where(UnusedProtein.A != 0)[0]
            self.Problem.MuDepIndices_A += [(UnusedProtein.row_names[i], UnusedProtein.col_names[0]) for i in nonZeroEntries if UnusedProtein.row_names[i]
                                            != 'Protein_'+input['ID'] and UnusedProtein.row_names[i] not in self.Problem.CompartmentDensities]

            self.setMu(self.Problem.Mu)

    ## !!! ##
    def eukaryoticDensities(self, totalAA=3.1, CompartmentRelationships=True, CompartmentComponents=False):
        Compartments = ['n', 'mIM', 'vM', 'mIMS', 'm', 'erM', 'mOM', 'x', 'c', 'cM', 'gM']
        Signs = ['L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L']
        totalAA = 3.1*0.71
        m_mIM = 0.66
        m_mIMS = 2
        m_mOM = 8
        DensityIndices = [self.Problem.LP.row_names.index(
            i) for i in self.Problem.CompartmentDensities]

        A = self.Problem.LP.A.toarray()
        A[numpy.min(DensityIndices):numpy.max(DensityIndices)+1, :] /= totalAA
        self.Problem.LP.A = scipy.sparse.coo_matrix(A)

        A0 = self.Problem.MuOneMatrix.A.toarray()
        A0[numpy.min(DensityIndices):numpy.max(DensityIndices)+1, :] /= totalAA
        self.Problem.MuOneMatrix.A = scipy.sparse.coo_matrix(A0)

        CompartmentMatrix = RBA_Matrix()
        A = numpy.ones((len(Compartments)+1, len(Compartments)))
        Eye = -numpy.eye(len(Compartments))
        A[0:len(Compartments), :] = Eye
        CompartmentMatrix.A = scipy.sparse.coo_matrix(A)
        CompartmentMatrix.b = numpy.array([float(0)]*len(Compartments)+[float(1)])
        CompartmentMatrix.f = numpy.array([float(0)]*len(Compartments))
        CompartmentMatrix.LB = numpy.array([float(0)]*len(Compartments))
        CompartmentMatrix.UB = numpy.array([float(1)]*len(Compartments))
        CompartmentMatrix.row_signs = ['L']*len(Compartments)+['E']
        # CompartmentMatrix.row_signs = ['E']*(len(Compartments)+1)
        CompartmentMatrix.row_names = ['n_density', 'mIM_density', 'vM_density', 'mIMS_density', 'm_density',
                                       'erM_density', 'mOM_density', 'x_density', 'cM_density', 'gM_density', 'c_density', 'TotalCapacity']
        CompartmentMatrix.col_names = ['F_n', 'F_mIM', 'F_vM', 'F_mIMS',
                                       'F_m', 'F_erM', 'F_mOM', 'F_x', 'F_cM', 'F_gM', 'F_c']
        # CompartmentMatrix.row_signs[CompartmentMatrix.col_names.index('F_m')]='E'

        if CompartmentRelationships:
            Anew = numpy.zeros((A.shape[0]+3, A.shape[1]))
            Anew[0:A.shape[0], :] = A
            CompartmentMatrix.row_names += ['m_mIM', 'm_mIMS', 'm_mOM']
            CompartmentMatrix.row_signs += ['E', 'E', 'E']
            CompartmentMatrix.b = numpy.array(list(CompartmentMatrix.b)+[float(0)]*3)
            Anew[CompartmentMatrix.row_names.index(
                'm_mIM'), CompartmentMatrix.col_names.index('F_m')] = float(1)
            Anew[CompartmentMatrix.row_names.index(
                'm_mIMS'), CompartmentMatrix.col_names.index('F_m')] = float(1)
            Anew[CompartmentMatrix.row_names.index(
                'm_mOM'), CompartmentMatrix.col_names.index('F_m')] = float(1)
            Anew[CompartmentMatrix.row_names.index(
                'm_mIM'), CompartmentMatrix.col_names.index('F_mIM')] = -m_mIM
            Anew[CompartmentMatrix.row_names.index(
                'm_mIMS'), CompartmentMatrix.col_names.index('F_mIMS')] = -m_mIMS
            Anew[CompartmentMatrix.row_names.index(
                'm_mOM'), CompartmentMatrix.col_names.index('F_mOM')] = -m_mOM
            CompartmentMatrix.A = scipy.sparse.coo_matrix(Anew)

        self.Problem.LP.addMatrix(matrix=CompartmentMatrix)
        self.Problem.MuOneMatrix.addMatrix(matrix=CompartmentMatrix)

        if CompartmentComponents:
            AlipidsA = numpy.zeros((7, len(Compartments)))
            Alipids = RBA_Matrix()
            Alipids.col_names = ['F_n', 'F_mIM', 'F_vM', 'F_mIMS',
                                 'F_m', 'F_erM', 'F_mOM', 'F_x', 'F_cM', 'F_gM', 'F_c']
            Alipids.row_names = ['M_pc_SC_c', 'M_pe_SC_c', 'M_ptd1ino_SC_c',
                                 'M_ps_SC_c', 'M_clpn_SC_m', 'M_pa_SC_c', 'M_ergst_c']
            Alipids.row_signs += ['E', 'E', 'E', 'E', 'E', 'E', 'E']
            Alipids.b = numpy.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
            Alipids.LB = numpy.array([float(0)]*len(Compartments))
            Alipids.UB = numpy.array([float(1)]*len(Compartments))
            Alipids.f = numpy.array([float(0)]*len(Compartments))
            AlipidsA[Alipids.row_names.index('M_pc_SC_c'), Alipids.col_names.index(
                'F_mIM')] = -0.0000883*totalAA
            AlipidsA[Alipids.row_names.index('M_pe_SC_c'), Alipids.col_names.index(
                'F_mIM')] = -0.00005852*totalAA
            AlipidsA[Alipids.row_names.index('M_ptd1ino_SC_c'),
                     Alipids.col_names.index('F_mIM')] = -0.00003377*totalAA
            AlipidsA[Alipids.row_names.index('M_ps_SC_c'), Alipids.col_names.index(
                'F_mIM')] = -0.00000873*totalAA
            AlipidsA[Alipids.row_names.index('M_clpn_SC_m'),
                     Alipids.col_names.index('F_mIM')] = -0.00002*totalAA
            AlipidsA[Alipids.row_names.index('M_pa_SC_c'), Alipids.col_names.index(
                'F_mIM')] = -0.0000039*totalAA
            AlipidsA[Alipids.row_names.index(
                'M_ergst_c'), Alipids.col_names.index('F_mIM')] = -0.008547*totalAA
            self.Problem.MuDepIndices_A += [('M_pc_SC_c', 'F_mIM'), ('M_pe_SC_c', 'F_mIM'), ('M_ptd1ino_SC_c', 'F_mIM'),
                                            ('M_ps_SC_c', 'F_mIM'), ('M_clpn_SC_m', 'F_mIM'), ('M_pa_SC_c', 'F_mIM'), ('M_ergst_c', 'F_mIM')]

            AlipidsA[Alipids.row_names.index(
                'M_pc_SC_c'), Alipids.col_names.index('F_mOM')] = -0.000636*totalAA
            AlipidsA[Alipids.row_names.index('M_pe_SC_c'), Alipids.col_names.index(
                'F_mOM')] = -0.0004822*totalAA
            AlipidsA[Alipids.row_names.index('M_ptd1ino_SC_c'),
                     Alipids.col_names.index('F_mOM')] = -0.0001289*totalAA
            AlipidsA[Alipids.row_names.index('M_ps_SC_c'), Alipids.col_names.index(
                'F_mOM')] = -0.0000167*totalAA
            AlipidsA[Alipids.row_names.index('M_clpn_SC_m'), Alipids.col_names.index(
                'F_mOM')] = -0.00004467*totalAA
            AlipidsA[Alipids.row_names.index('M_pa_SC_c'), Alipids.col_names.index(
                'F_mOM')] = -0.0000696*totalAA

            self.Problem.MuDepIndices_A += [('M_pc_SC_c', 'F_mOM'), ('M_pe_SC_c', 'F_mOM'), ('M_ptd1ino_SC_c',
                                                                                             'F_mOM'), ('M_ps_SC_c', 'F_mOM'), ('M_clpn_SC_m', 'F_mOM'), ('M_pa_SC_c', 'F_mOM')]

            Alipids.A = scipy.sparse.coo_matrix(AlipidsA)
            Alipids.mapIndices()
            self.Problem.LP.updateMatrix(Alipids, Ainds=[('M_pc_SC_c', 'F_mIM'), ('M_pe_SC_c', 'F_mIM'), ('M_ptd1ino_SC_c', 'F_mIM'), ('M_ps_SC_c', 'F_mIM'), ('M_clpn_SC_m', 'F_mIM'), ('M_pa_SC_c', 'F_mIM'), (
                'M_ergst_c', 'F_mIM'), ('M_pc_SC_c', 'F_mOM'), ('M_pe_SC_c', 'F_mOM'), ('M_ptd1ino_SC_c', 'F_mOM'), ('M_ps_SC_c', 'F_mOM'), ('M_clpn_SC_m', 'F_mOM'), ('M_pa_SC_c', 'F_mOM')])
            self.Problem.MuOneMatrix.updateMatrix(Alipids, Ainds=[('M_pc_SC_c', 'F_mIM'), ('M_pe_SC_c', 'F_mIM'), ('M_ptd1ino_SC_c', 'F_mIM'), ('M_ps_SC_c', 'F_mIM'), ('M_clpn_SC_m', 'F_mIM'), (
                'M_pa_SC_c', 'F_mIM'), ('M_ergst_c', 'F_mIM'), ('M_pc_SC_c', 'F_mOM'), ('M_pe_SC_c', 'F_mOM'), ('M_ptd1ino_SC_c', 'F_mOM'), ('M_ps_SC_c', 'F_mOM'), ('M_clpn_SC_m', 'F_mOM'), ('M_pa_SC_c', 'F_mOM')])

    ## !!! ##
    def eukaryoticDensities2(self, totalAA=3.1, CompartmentRelationships=True, CompartmentComponents=False):
        Compartments = ['n', 'mIM', 'vM', 'mIMS', 'm', 'erM', 'mOM', 'x', 'c', 'cM', 'gM']
        totalAA = 3.1*0.69
        m_mIM = 1.11
        m_mIMS = 0.7
        m_mOM = 7.2
        DensityIndices = [self.Problem.LP.row_names.index(
            i) for i in self.Problem.CompartmentDensities]

        A = self.Problem.LP.A.toarray()
        A[numpy.min(DensityIndices):numpy.max(DensityIndices)+1, :] /= totalAA
        self.Problem.LP.A = scipy.sparse.coo_matrix(A)

        A0 = self.Problem.MuOneMatrix.A.toarray()
        A0[numpy.min(DensityIndices):numpy.max(DensityIndices)+1, :] /= totalAA
        self.Problem.MuOneMatrix.A = scipy.sparse.coo_matrix(A0)

        CompartmentMatrix = RBA_Matrix()
        A = numpy.ones((len(Compartments)+1, len(Compartments)))
        Eye = -numpy.eye(len(Compartments))
        A[0:len(Compartments), :] = Eye
        CompartmentMatrix.A = scipy.sparse.coo_matrix(A)
        CompartmentMatrix.b = numpy.array([float(0)]*len(Compartments)+[float(1)])
        CompartmentMatrix.f = numpy.array([float(0)]*len(Compartments))
        CompartmentMatrix.LB = numpy.array([float(0)]*len(Compartments))
        CompartmentMatrix.UB = numpy.array([float(1)]*len(Compartments))
        CompartmentMatrix.row_signs = ['L']*(len(Compartments)+1)
        # CompartmentMatrix.row_signs = ['E']*(len(Compartments)+1)
        CompartmentMatrix.row_names = ['n_density', 'mIM_density', 'vM_density', 'mIMS_density', 'm_density',
                                       'erM_density', 'mOM_density', 'x_density', 'cM_density', 'gM_density', 'c_density', 'TotalCapacity']
        CompartmentMatrix.col_names = ['F_n', 'F_mIM', 'F_vM', 'F_mIMS',
                                       'F_m', 'F_erM', 'F_mOM', 'F_x', 'F_cM', 'F_gM', 'F_c']
        # CompartmentMatrix.row_signs[CompartmentMatrix.col_names.index('F_m')]='E'

        if CompartmentRelationships:
            Anew = numpy.zeros((A.shape[0]+3, A.shape[1]))
            Anew[0:A.shape[0], :] = A
            CompartmentMatrix.row_names += ['m_mIM', 'm_mIMS', 'm_mOM']
            CompartmentMatrix.row_signs += ['E', 'E', 'E']
            CompartmentMatrix.b = numpy.array(list(CompartmentMatrix.b)+[float(0)]*3)
            Anew[CompartmentMatrix.row_names.index(
                'm_mIM'), CompartmentMatrix.col_names.index('F_m')] = float(1)
            Anew[CompartmentMatrix.row_names.index(
                'm_mIMS'), CompartmentMatrix.col_names.index('F_m')] = float(1)
            Anew[CompartmentMatrix.row_names.index(
                'm_mOM'), CompartmentMatrix.col_names.index('F_m')] = float(1)
            Anew[CompartmentMatrix.row_names.index(
                'm_mIM'), CompartmentMatrix.col_names.index('F_mIM')] = -m_mIM
            Anew[CompartmentMatrix.row_names.index(
                'm_mIMS'), CompartmentMatrix.col_names.index('F_mIMS')] = -m_mIMS
            Anew[CompartmentMatrix.row_names.index(
                'm_mOM'), CompartmentMatrix.col_names.index('F_mOM')] = -m_mOM
            CompartmentMatrix.A = scipy.sparse.coo_matrix(Anew)

        self.Problem.LP.addMatrix(matrix=CompartmentMatrix)
        self.Problem.MuOneMatrix.addMatrix(matrix=CompartmentMatrix)

        if CompartmentComponents:

            PC_mIM = 0.0000883
            PE_mIM = 0.00005852
            PI_mIM = 0.00003377
            PS_mIM = 0.00000873
            CL_mIM = 0.00002
            PA_mIM = 0.0000039
            ES_mIM = 0.008547

            PC_mOM = 0.000636
            PE_mOM = 0.0004822
            PI_mOM = 0.0001289
            PS_mOM = 0.0000167
            CL_mOM = 0.00004467
            PA_mOM = 0.0000696
            ES_mOM = 0.0

            ConstraintMatrix = numpy.zeros((7, 0))
            Alipids = RBA_Matrix()
            Alipids.col_names = []
            Alipids.row_names = ['M_pc_SC_c', 'M_pe_SC_c', 'M_ptd1ino_SC_c',
                                 'M_ps_SC_c', 'M_clpn_SC_m', 'M_pa_SC_c', 'M_ergst_c']
            Alipids.row_signs = [
                self.Problem.LP.row_signs[self.Problem.LP.row_names.index(i)] for i in Alipids.row_names]
            Alipids.b = numpy.array(
                [self.Problem.LP.b[self.Problem.LP.row_names.index(i)] for i in Alipids.row_names])
            Alipids.LB = numpy.array([])
            Alipids.UB = numpy.array([])
            Alipids.f = numpy.array([])
            MudepIndices = []
            for pc in self.ModelStructure.ProcessInfo.Elements.keys():
                if self.ModelStructure.ProcessInfo.Elements[pc]['ID'] not in self.Problem.LP.col_names:
                    continue
                ConstraintMatrixNew = numpy.zeros(
                    (ConstraintMatrix.shape[0], ConstraintMatrix.shape[1]+1))
                ConstraintMatrixNew[:, 0:ConstraintMatrix.shape[1]] = ConstraintMatrix
                Alipids.col_names.append(self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
#                Alipids.LB = numpy.array(list(Alipids.LB).append(list(self.Problem.LP.LB)[
#                                         self.Problem.LP.col_names.index(self.ModelStructure.ProcessInfo.Elements[pc]['ID'])]))
#                Alipids.UB = numpy.array(list(Alipids.UB).append(list(self.Problem.LP.UB)[
#                                         self.Problem.LP.col_names.index(self.ModelStructure.ProcessInfo.Elements[pc]['ID'])]))
#                Alipids.f = numpy.array(list(Alipids.f).append(list(self.Problem.LP.f)[
#                                        self.Problem.LP.col_names.index(self.ModelStructure.ProcessInfo.Elements[pc]['ID'])]))
                Alipids.LB = numpy.concatenate([Alipids.LB, numpy.array(
                    list(self.Problem.LP.LB)[self.Problem.LP.col_names.index(self.ModelStructure.ProcessInfo.Elements[pc]['ID'])])])
                Alipids.UB = numpy.concatenate([Alipids.UB, numpy.array(
                    list(self.Problem.LP.UB)[self.Problem.LP.col_names.index(self.ModelStructure.ProcessInfo.Elements[pc]['ID'])])])
                Alipids.f = numpy.concatenate([Alipids.f, numpy.array(
                    list(self.Problem.LP.f)[self.Problem.LP.col_names.index(self.ModelStructure.ProcessInfo.Elements[pc]['ID'])])])
                for p in self.ModelStructure.ProcessInfo.Elements[pc]['Composition'].keys():
                    lE = sum(list(self.ModelStructure.ProteinInfo.Elements[p]['AAcomposition'].values(
                    )))*self.ModelStructure.ProcessInfo.Elements[pc]['Composition'][p]
                    if self.ModelStructure.ProteinInfo.Elements[p]['Compartment'] == 'mOM':
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pc_SC_c'), ConstraintMatrix.shape[1]] -= PC_mOM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pe_SC_c'), ConstraintMatrix.shape[1]] -= PE_mOM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_ptd1ino_SC_c'), ConstraintMatrix.shape[1]] -= PI_mOM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_ps_SC_c'), ConstraintMatrix.shape[1]] -= PS_mOM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_clpn_SC_m'), ConstraintMatrix.shape[1]] -= CL_mOM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pa_SC_c'), ConstraintMatrix.shape[1]] -= PA_mOM*lE/totalAA
                        MudepIndices += ('M_pc_SC_c',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                        MudepIndices += ('M_pe_SC_c',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                        MudepIndices += ('M_ptd1ino_SC_c',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                        MudepIndices += ('M_ps_SC_c',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                        MudepIndices += ('M_clpn_SC_m',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                        MudepIndices += ('M_pa_SC_c',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                    elif self.ModelStructure.ProteinInfo.Elements[p]['Compartment'] == 'mIM':
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pc_SC_c'), ConstraintMatrix.shape[1]] -= PC_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pe_SC_c'), ConstraintMatrix.shape[1]] -= PE_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_ptd1ino_SC_c'), ConstraintMatrix.shape[1]] -= PI_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_ps_SC_c'), ConstraintMatrix.shape[1]] -= PS_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_clpn_SC_m'), ConstraintMatrix.shape[1]] -= CL_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pa_SC_c'), ConstraintMatrix.shape[1]] -= PA_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_ergst_c'), ConstraintMatrix.shape[1]] -= ES_mIM*lE/totalAA
                        MudepIndices += ('M_pc_SC_c',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                        MudepIndices += ('M_pe_SC_c',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                        MudepIndices += ('M_ptd1ino_SC_c',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                        MudepIndices += ('M_ps_SC_c',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                        MudepIndices += ('M_clpn_SC_m',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                        MudepIndices += ('M_pa_SC_c',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                        MudepIndices += ('M_ergst_c',
                                         self.ModelStructure.ProcessInfo.Elements[pc]['ID'])
                ConstraintMatrix = ConstraintMatrixNew

            for e in self.ModelStructure.EnzymeInfo.Elements.keys():
                if e not in self.Problem.LP.col_names:
                    continue
                ConstraintMatrixNew = numpy.zeros(
                    (ConstraintMatrix.shape[0], ConstraintMatrix.shape[1]+1))
                ConstraintMatrixNew[:, 0:ConstraintMatrix.shape[1]] = ConstraintMatrix
                Alipids.col_names.append(e)
#                xnew = list(self.Problem.LP.LB)[self.Problem.LP.col_names.index(e)]
                Alipids.LB = numpy.concatenate([Alipids.LB, numpy.array(
                    list(self.Problem.LP.LB)[self.Problem.LP.col_names.index(e)])])
                Alipids.UB = numpy.concatenate([Alipids.UB, numpy.array(
                    list(self.Problem.LP.UB)[self.Problem.LP.col_names.index(e)])])
                Alipids.f = numpy.concatenate([Alipids.f, numpy.array(
                    list(self.Problem.LP.f)[self.Problem.LP.col_names.index(e)])])
                # Alipids.LB = numpy.array(list(Alipids.LB).append(xnew))
#                Alipids.UB = numpy.array(list(Alipids.UB).append(
#                    list(self.Problem.LP.UB)[self.Problem.LP.col_names.index(e)]))
                # Alipids.f = numpy.array(list(Alipids.f).append(
                #    list(self.Problem.LP.f)[self.Problem.LP.col_names.index(e)]))
                for p in self.ModelStructure.EnzymeInfo.Elements[e]['Subunits'].keys():
                    lE = sum(
                        list(self.ModelStructure.ProteinInfo.Elements[p]['AAcomposition'].values()))
                    lE *= self.ModelStructure.EnzymeInfo.Elements[e]['Subunits'][p]['StochFac']
                    if self.ModelStructure.ProteinInfo.Elements[p]['Compartment'] == 'mOM':
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pc_SC_c'), ConstraintMatrix.shape[1]] -= PC_mOM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pe_SC_c'), ConstraintMatrix.shape[1]] -= PE_mOM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_ptd1ino_SC_c'), ConstraintMatrix.shape[1]] -= PI_mOM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_ps_SC_c'), ConstraintMatrix.shape[1]] -= PS_mOM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_clpn_SC_m'), ConstraintMatrix.shape[1]] -= CL_mOM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pa_SC_c'), ConstraintMatrix.shape[1]] -= PA_mOM*lE/totalAA
                        MudepIndices += ('M_pc_SC_c', e)
                        MudepIndices += ('M_pe_SC_c', e)
                        MudepIndices += ('M_ptd1ino_SC_c', e)
                        MudepIndices += ('M_ps_SC_c', e)
                        MudepIndices += ('M_clpn_SC_m', e)
                        MudepIndices += ('M_pa_SC_c', e)
                    elif self.ModelStructure.ProteinInfo.Elements[p]['Compartment'] == 'mIM':
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pc_SC_c'), ConstraintMatrix.shape[1]] -= PC_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pe_SC_c'), ConstraintMatrix.shape[1]] -= PE_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_ptd1ino_SC_c'), ConstraintMatrix.shape[1]] -= PI_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_ps_SC_c'), ConstraintMatrix.shape[1]] -= PS_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_clpn_SC_m'), ConstraintMatrix.shape[1]] -= CL_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_pa_SC_c'), ConstraintMatrix.shape[1]] -= PA_mIM*lE/totalAA
                        ConstraintMatrixNew[Alipids.col_names.index(
                            'M_ergst_c'), ConstraintMatrix.shape[1]] -= ES_mIM*lE/totalAA
                        MudepIndices += ('M_pc_SC_c', e)
                        MudepIndices += ('M_pe_SC_c', e)
                        MudepIndices += ('M_ptd1ino_SC_c', e)
                        MudepIndices += ('M_ps_SC_c', e)
                        MudepIndices += ('M_clpn_SC_m', e)
                        MudepIndices += ('M_pa_SC_c', e)
                        MudepIndices += ('M_ergst_c', e)
                ConstraintMatrix = ConstraintMatrixNew

            self.Problem.MuDepIndices_A += MudepIndices
            Alipids.A = scipy.sparse.coo_matrix(ConstraintMatrix)
            Alipids.mapIndices()
            self.Problem.LP.updateMatrix(Alipids, Ainds=MudepIndices)
            self.Problem.LP.updateMatrix(MuOneMatrix, Ainds=MudepIndices)

    ## !!! ##
    def eukaryoticDensities3(self, totalAA=3.1, VolumeFraction=False, CompartmentRelationships=True, CompartmentComponents=False):
        Compartments = ['n', 'mIM', 'vM', 'mIMS', 'm', 'erM', 'mOM', 'x', 'c', 'cM', 'gM']
        totalAA = 3.1*0.91
        m_mIM = 0.66
        m_mIMS = 2
        m_mOM = 8
        DensityIndices = [self.Problem.LP.row_names.index(
            i) for i in self.Problem.CompartmentDensities]

        A = self.Problem.LP.A.toarray()
        # A[numpy.min(DensityIndices):numpy.max(DensityIndices)+1, :] /= totalAA
        # self.Problem.LP.A = scipy.sparse.coo_matrix(A)

        A0 = self.Problem.MuOneMatrix.A.toarray()
        # A0[numpy.min(DensityIndices):numpy.max(DensityIndices)+1, :] /= totalAA
        # self.Problem.MuOneMatrix.A = scipy.sparse.coo_matrix(A0)

        OccupationMatrix = RBA_Matrix()
        # A = numpy.ones((len(Compartments)+1, len(Compartments)))
        A = -numpy.eye(len(Compartments))
        # Eye = -numpy.eye(len(Compartments))
        # A[0:len(Compartments), :] = Eye
        OccupationMatrix.A = scipy.sparse.coo_matrix(A)
#        OccupationMatrix.b = numpy.array([-0.209*totalAA]+[float(0)]*(len(Compartments)-1)+[float(totalAA)])
        OccupationMatrix.b = numpy.array([-0.209*totalAA]+[float(0)]*(len(Compartments)-1))
        OccupationMatrix.f = numpy.array([float(0)]*len(Compartments))
        OccupationMatrix.LB = numpy.array([float(0)]*len(Compartments))
        OccupationMatrix.UB = numpy.array([float(totalAA)]*len(Compartments))
        # OccupationMatrix.row_signs = ['E']*(len(Compartments))+['L']
        OccupationMatrix.row_signs = ['E']*(len(Compartments))
        # OccupationMatrix.row_names = ['n_density', 'mIM_density', 'vM_density', 'mIMS_density', 'm_density',
        #                              'erM_density', 'mOM_density', 'x_density', 'cM_density', 'gM_density', 'c_density', 'TotalProtein']
        OccupationMatrix.row_names = ['n_density', 'mIM_density', 'vM_density', 'mIMS_density',
                                      'm_density', 'erM_density', 'mOM_density', 'x_density', 'cM_density', 'gM_density', 'c_density']
        OccupationMatrix.col_names = ['O_n', 'O_mIM', 'O_vM', 'O_mIMS',
                                      'O_m', 'O_erM', 'O_mOM', 'O_x', 'O_cM', 'O_gM', 'O_c']
        # CompartmentMatrix.row_signs[CompartmentMatrix.col_names.index('F_m')]='E'
        OccupationMatrix.mapIndices()
        self.Problem.LP.addMatrix(matrix=OccupationMatrix)
        self.Problem.MuOneMatrix.addMatrix(matrix=OccupationMatrix)

        CompartmentMatrix = RBA_Matrix()
        if VolumeFraction:
            A = numpy.eye(len(Compartments))*5/float(totalAA)
        else:
            A = numpy.eye(len(Compartments))/float(totalAA)
        CompartmentMatrix.A = scipy.sparse.coo_matrix(A)
        CompartmentMatrix.b = numpy.array([float(0)]*len(Compartments))
        CompartmentMatrix.f = numpy.array([float(0)]*len(Compartments))
        CompartmentMatrix.LB = numpy.array([float(0)]*len(Compartments))
        CompartmentMatrix.UB = numpy.array([float(totalAA)]*len(Compartments))
        CompartmentMatrix.row_signs = ['L']*(len(Compartments))
#        CompartmentMatrix.row_signs = ['E']*(len(Compartments))
        CompartmentMatrix.row_names = ['n_volume', 'mIM_volume', 'vM_volume', 'mIMS_volume',
                                       'm_volume', 'erM_volume', 'mOM_volume', 'x_volume', 'cM_volume', 'gM_volume', 'c_volume']
        CompartmentMatrix.col_names = ['O_n', 'O_mIM', 'O_vM', 'O_mIMS',
                                       'O_m', 'O_erM', 'O_mOM', 'O_x', 'O_cM', 'O_gM', 'O_c']
        CompartmentMatrix.mapIndices()
        self.Problem.LP.addMatrix(matrix=CompartmentMatrix)
        self.Problem.MuOneMatrix.addMatrix(matrix=CompartmentMatrix)

        VolumeMatrix = RBA_Matrix()
        A = numpy.ones((len(Compartments)+1, len(Compartments)))
        Eye = -numpy.eye(len(Compartments))
        A[0:len(Compartments), :] = Eye
        # A[len(Compartments), [1, 5, 6, 8, 9]] = 0
        # A[len(Compartments), 8] = 0
        VolumeMatrix.A = scipy.sparse.coo_matrix(A)
        VolumeMatrix.b = numpy.array([float(0)]*len(Compartments)+[float(1)])
        VolumeMatrix.f = numpy.array([float(0)]*len(Compartments))
        VolumeMatrix.LB = numpy.array([float(0)]*len(Compartments))
        VolumeMatrix.UB = numpy.array([float(1)]*len(Compartments))
        VolumeMatrix.row_signs = ['L']*(len(Compartments))+['E']
        # VolumeMatrix.row_signs = ['E']*(len(Compartments))+['E']
        VolumeMatrix.row_names = ['n_volume', 'mIM_volume', 'vM_volume', 'mIMS_volume', 'm_volume',
                                  'erM_volume', 'mOM_volume', 'x_volume', 'cM_volume', 'gM_volume', 'c_volume', 'TotalVolume']
        VolumeMatrix.col_names = ['F_n', 'F_mIM', 'F_vM', 'F_mIMS',
                                  'F_m', 'F_erM', 'F_mOM', 'F_x', 'F_cM', 'F_gM', 'F_c']

        if not CompartmentRelationships:
            VolumeMatrix.mapIndices()
            self.Problem.LP.addMatrix(matrix=VolumeMatrix)
            self.Problem.MuOneMatrix.addMatrix(matrix=VolumeMatrix)

        if CompartmentRelationships:
            Anew = numpy.zeros((A.shape[0]+3, A.shape[1]))
            Anew[0:A.shape[0], :] = A
            VolumeMatrix.row_names += ['m_mIM', 'm_mIMS', 'm_mOM']
            VolumeMatrix.row_signs += ['E', 'E', 'E']
            VolumeMatrix.b = numpy.array(list(VolumeMatrix.b)+[float(0)]*3)
            Anew[VolumeMatrix.row_names.index(
                'm_mIM'), VolumeMatrix.col_names.index('F_m')] = float(1)
            Anew[VolumeMatrix.row_names.index(
                'm_mIMS'), VolumeMatrix.col_names.index('F_m')] = float(1)
            Anew[VolumeMatrix.row_names.index(
                'm_mOM'), VolumeMatrix.col_names.index('F_m')] = float(1)
            Anew[VolumeMatrix.row_names.index(
                'm_mIM'), VolumeMatrix.col_names.index('F_mIM')] = -m_mIM
            Anew[VolumeMatrix.row_names.index(
                'm_mIMS'), VolumeMatrix.col_names.index('F_mIMS')] = -m_mIMS
            Anew[VolumeMatrix.row_names.index(
                'm_mOM'), VolumeMatrix.col_names.index('F_mOM')] = -m_mOM
            VolumeMatrix.A = scipy.sparse.coo_matrix(Anew)

            VolumeMatrix.mapIndices()
            self.Problem.LP.addMatrix(matrix=VolumeMatrix)
            self.Problem.MuOneMatrix.addMatrix(matrix=VolumeMatrix)

        if CompartmentComponents:

            PC_mIM = 0.0000883
            PE_mIM = 0.00005852
            PI_mIM = 0.00003377
            PS_mIM = 0.00000873
            CL_mIM = 0.00002
            PA_mIM = 0.0000039
            ES_mIM = 0.008547

            PC_mOM = 0.000636
            PE_mOM = 0.0004822
            PI_mOM = 0.0001289
            PS_mOM = 0.0000167
            CL_mOM = 0.00004467
            PA_mOM = 0.0000696
            ES_mOM = 0.0

            PC_vM = 0.0003635
            PE_vM = 0.4156
            PI_vM = 0.0001297
            PS_vM = 0.00003435
            CL_vM = 0.0000068
            PA_vM = 0.0000186
            ES_vM = 0.0142

            PC_n = 0.000055
            PE_n = 0.000035
            PI_n = 0.000017
            PS_n = 0.0000072
            CL_n = 0.0
            PA_n = 0.0000031
            ES_n = 0.0086

            PC_gM = 0.00043
            PE_gM = 0.00044
            PI_gM = 0.00041
            PS_gM = 0.0
            CL_gM = 0.00022
            PA_gM = 0.0
            ES_gM = 0.0

            PC_n = 0.0
            PE_n = 0.0
            PI_n = 0.0
            PS_n = 0.0
            CL_n = 0.0
            PA_n = 0.0
            ES_n = 0.0

            PC_gM = 0.0
            PE_gM = 0.0
            PI_gM = 0.0
            PS_gM = 0.0
            CL_gM = 0.0
            PA_gM = 0.0
            ES_gM = 0.0

            PC_vM = 0.0
            PE_vM = 0.0
            PI_vM = 0.0
            PS_vM = 0.0
            CL_vM = 0.0
            PA_vM = 0.0
            ES_vM = 0.0

            PC_mIM = 0.0
            PE_mIM = 0.0
            PI_mIM = 0.0
            PS_mIM = 0.0
            CL_mIM = 0.0
            PA_mIM = 0.0
            ES_mIM = 0.0

            PC_mOM = 0.0
            PE_mOM = 0.0
            PI_mOM = 0.0
            PS_mOM = 0.0
            CL_mOM = 0.0
            PA_mOM = 0.0
            ES_mOM = 0.0

            Alipids = RBA_Matrix()
            Alipids.col_names = ['F_mIM', 'F_mOM', 'F_vM', 'F_n', 'F_gM']
            Alipids.row_names = ['M_pc_SC_c', 'M_pe_SC_c', 'M_ptd1ino_SC_c',
                                 'M_ps_SC_c', 'M_clpn_SC_m', 'M_pa_SC_c', 'M_ergst_c']
            Alipids.row_signs = [
                self.Problem.LP.row_signs[self.Problem.LP.row_names.index(i)] for i in Alipids.row_names]
            Alipids.b = numpy.array(
                [self.Problem.LP.b[self.Problem.LP.row_names.index(i)] for i in Alipids.row_names])
            Alipids.LB = numpy.array([0, 0, 0, 0, 0])
            Alipids.UB = numpy.array([1, 1, 1, 1, 1])
            Alipids.f = numpy.array([0, 0, 0, 0, 0])
            LipidMatrix = numpy.zeros((7, 5))

            LipidMatrix[Alipids.row_names.index(
                'M_pc_SC_c'), Alipids.col_names.index('F_mIM')] = PC_mIM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_pe_SC_c'), Alipids.col_names.index('F_mIM')] = PE_mIM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ptd1ino_SC_c'), Alipids.col_names.index('F_mIM')] = PI_mIM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ps_SC_c'), Alipids.col_names.index('F_mIM')] = PS_mIM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_clpn_SC_m'), Alipids.col_names.index('F_mIM')] = CL_mIM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_pa_SC_c'), Alipids.col_names.index('F_mIM')] = PA_mIM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ergst_c'), Alipids.col_names.index('F_mIM')] = ES_mIM/totalAA

            LipidMatrix[Alipids.row_names.index(
                'M_pc_SC_c'), Alipids.col_names.index('F_mOM')] = PC_mOM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_pe_SC_c'), Alipids.col_names.index('F_mOM')] = PE_mOM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ptd1ino_SC_c'), Alipids.col_names.index('F_mOM')] = PI_mOM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ps_SC_c'), Alipids.col_names.index('F_mOM')] = PS_mOM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_clpn_SC_m'), Alipids.col_names.index('F_mOM')] = CL_mOM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_pa_SC_c'), Alipids.col_names.index('F_mOM')] = PA_mOM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ergst_c'), Alipids.col_names.index('F_mOM')] = ES_mOM/totalAA

            LipidMatrix[Alipids.row_names.index(
                'M_pc_SC_c'), Alipids.col_names.index('F_vM')] = PC_vM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_pe_SC_c'), Alipids.col_names.index('F_vM')] = PE_vM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ptd1ino_SC_c'), Alipids.col_names.index('F_vM')] = PI_vM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ps_SC_c'), Alipids.col_names.index('F_vM')] = PS_vM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_clpn_SC_m'), Alipids.col_names.index('F_vM')] = CL_vM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_pa_SC_c'), Alipids.col_names.index('F_vM')] = PA_vM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ergst_c'), Alipids.col_names.index('F_vM')] = ES_vM/totalAA

            LipidMatrix[Alipids.row_names.index(
                'M_pc_SC_c'), Alipids.col_names.index('F_n')] = PC_n/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_pe_SC_c'), Alipids.col_names.index('F_n')] = PE_n/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ptd1ino_SC_c'), Alipids.col_names.index('F_n')] = PI_n/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ps_SC_c'), Alipids.col_names.index('F_n')] = PS_n/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_clpn_SC_m'), Alipids.col_names.index('F_n')] = CL_n/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_pa_SC_c'), Alipids.col_names.index('F_n')] = PA_n/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ergst_c'), Alipids.col_names.index('F_n')] = ES_n/totalAA

            LipidMatrix[Alipids.row_names.index(
                'M_pc_SC_c'), Alipids.col_names.index('F_gM')] = PC_gM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_pe_SC_c'), Alipids.col_names.index('F_gM')] = PE_gM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ptd1ino_SC_c'), Alipids.col_names.index('F_gM')] = PI_gM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ps_SC_c'), Alipids.col_names.index('F_gM')] = PS_gM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_clpn_SC_m'), Alipids.col_names.index('F_gM')] = CL_gM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_pa_SC_c'), Alipids.col_names.index('F_gM')] = PA_gM/totalAA
            LipidMatrix[Alipids.row_names.index(
                'M_ergst_c'), Alipids.col_names.index('F_gM')] = ES_gM/totalAA

            MudepIndices = [('M_pc_SC_c', i) for i in Alipids.col_names]+[('M_pe_SC_c', i) for i in Alipids.col_names]+[('M_ptd1ino_SC_c', i) for i in Alipids.col_names]+[('M_ps_SC_c', i)
                                                                                                                                                                           for i in Alipids.col_names]+[('M_clpn_SC_m', i) for i in Alipids.col_names]+[('M_pa_SC_c', i) for i in Alipids.col_names]+[('M_ergst_c', i) for i in Alipids.col_names]
            self.Problem.MuDepIndices_A += MudepIndices
            Alipids.A = scipy.sparse.coo_matrix(LipidMatrix)
            Alipids.mapIndices()
            self.Problem.LP.updateMatrix(Alipids, Ainds=MudepIndices)
            self.Problem.MuOneMatrix.updateMatrix(Alipids, Ainds=MudepIndices)

    ## !!! ##
    def eukaryoticDensities4(self, CompartmentRelationships=True):
        Compartments = ['n', 'mIM', 'vM', 'mIMS', 'm', 'erM', 'mOM', 'x', 'c', 'cM', 'gM']
        totalAA = 3.1*0.91
        DensityIndices = [self.Problem.LP.row_names.index(
            i) for i in self.Problem.CompartmentDensities]

        A = self.Problem.LP.A.toarray()

        A0 = self.Problem.MuOneMatrix.A.toarray()

        OccupationMatrix = RBA_Matrix()
        A = numpy.ones((len(Compartments)+1, len(Compartments)))
        Eye = -numpy.eye(len(Compartments))
        A[0:len(Compartments), :] = Eye

        OccupationMatrix.b = numpy.array(list([float(0)]*len(Compartments))+[totalAA])
        OccupationMatrix.f = numpy.array([float(0)]*len(Compartments))
        OccupationMatrix.LB = numpy.array([float(0)]*len(Compartments))
        OccupationMatrix.UB = numpy.array([float(totalAA)]*len(Compartments))
        OccupationMatrix.row_signs = ['E']*(len(Compartments)+1)
        OccupationMatrix.row_names = ['n_density', 'mIM_density', 'vM_density', 'mIMS_density',
                                      'm_density', 'erM_density', 'mOM_density', 'x_density', 'cM_density', 'gM_density', 'c_density', 'O_total']
        OccupationMatrix.col_names = ['O_n', 'O_mIM', 'O_vM', 'O_mIMS',
                                      'O_m', 'O_erM', 'O_mOM', 'O_x', 'O_cM', 'O_gM', 'O_c']

        if CompartmentRelationships:
            m_mIM = 0.5
            m_mIMS = 1
            m_mOM = 5
            Anew = numpy.zeros((A.shape[0]+3, A.shape[1]))
            Anew[0:A.shape[0], :] = A
            OccupationMatrix.row_names += ['m_mIM', 'm_mIMS', 'm_mOM']
            OccupationMatrix.row_signs += ['E', 'E', 'E']
            OccupationMatrix.b = numpy.array(list(OccupationMatrix.b)+[float(0)]*3)
            Anew[OccupationMatrix.row_names.index(
                'm_mIM'), OccupationMatrix.col_names.index('O_m')] = float(1)
            Anew[OccupationMatrix.row_names.index(
                'm_mIMS'), OccupationMatrix.col_names.index('O_m')] = float(1)
            Anew[OccupationMatrix.row_names.index(
                'm_mOM'), OccupationMatrix.col_names.index('O_m')] = float(1)
            Anew[OccupationMatrix.row_names.index(
                'm_mIM'), OccupationMatrix.col_names.index('O_mIM')] = -m_mIM
            Anew[OccupationMatrix.row_names.index(
                'm_mIMS'), OccupationMatrix.col_names.index('O_mIMS')] = -m_mIMS
            Anew[OccupationMatrix.row_names.index(
                'm_mOM'), OccupationMatrix.col_names.index('O_mOM')] = -m_mOM

            OccupationMatrix.A = scipy.sparse.coo_matrix(Anew)
        else:
            OccupationMatrix.A = scipy.sparse.coo_matrix(A)

        OccupationMatrix.mapIndices()
        self.Problem.LP.addMatrix(matrix=OccupationMatrix)
        self.Problem.MuOneMatrix.addMatrix(matrix=OccupationMatrix)
        # {'Index':{'Param1':'+','Param2':'+','Param2':'-'}}

        #Type: 'Sum'#
        # {'Index':'Param1'}
        self.Problem.MuDependencies['FromParameters']['b'].update(
            {'n_density': 'AAres_PG_nucleus_DNA'})
        self.Problem.MuDependencies['FromParameters']['b'].update(
            {'O_total': {'Equation': 'amino_acid_concentration_total - AAres_PG_secreted_Euk', 'Variables': ['amino_acid_concentration_total', 'AAres_PG_secreted_Euk']}})
        self.Problem.MuDependencies['FromMatrix']['b'].remove('n_density')
        self.Problem.MuDependencies['FromMatrix']['b'].remove('mIM_density')
        self.Problem.MuDependencies['FromMatrix']['b'].remove('vM_density')
        self.Problem.MuDependencies['FromMatrix']['b'].remove('mIMS_density')
        self.Problem.MuDependencies['FromMatrix']['b'].remove('m_density')
        self.Problem.MuDependencies['FromMatrix']['b'].remove('erM_density')
        self.Problem.MuDependencies['FromMatrix']['b'].remove('mOM_density')
        self.Problem.MuDependencies['FromMatrix']['b'].remove('x_density')
        self.Problem.MuDependencies['FromMatrix']['b'].remove('cM_density')
        self.Problem.MuDependencies['FromMatrix']['b'].remove('gM_density')
        self.Problem.MuDependencies['FromMatrix']['b'].remove('c_density')

    ## !!! ##
    def eukaryoticDensities_calibration(self, CompartmentRelationships=False, mitoProportions={}, amino_acid_concentration_total='amino_acid_concentration_total'):
        Compartments = ['n', 'mIM', 'vM', 'mIMS', 'm', 'erM', 'mOM', 'x', 'c', 'cM', 'gM']
        totalAA_parameter = amino_acid_concentration_total
        totalAA = 3.1*0.91
        DensityIndices = [self.Problem.LP.row_names.index(
            i) for i in self.Problem.CompartmentDensities]

        A = self.Problem.LP.A.toarray()

        A0 = self.Problem.MuOneMatrix.A.toarray()

        OccupationMatrix = RBA_Matrix()
        A = numpy.ones((len(Compartments)+1, len(Compartments)))
        Eye = -numpy.eye(len(Compartments))
        A[0:len(Compartments), :] = Eye

        OccupationMatrix.b = numpy.array(list([float(0)]*len(Compartments))+[totalAA])
        OccupationMatrix.f = numpy.array([float(0)]*len(Compartments))
        OccupationMatrix.LB = numpy.array([float(0)]*len(Compartments))
        OccupationMatrix.UB = numpy.array([float(totalAA)]*len(Compartments))
        OccupationMatrix.row_signs = ['E']*(len(Compartments)+1)
        OccupationMatrix.row_names = ['n_density', 'mIM_density', 'vM_density', 'mIMS_density',
                                      'm_density', 'erM_density', 'mOM_density', 'x_density', 'cM_density', 'gM_density', 'c_density', 'O_total']
        OccupationMatrix.col_names = ['O_n', 'O_mIM', 'O_vM', 'O_mIMS',
                                      'O_m', 'O_erM', 'O_mOM', 'O_x', 'O_cM', 'O_gM', 'O_c']

        if CompartmentRelationships:
            if len(list(mitoProportions.keys())) == 3:
                m_mIM = mitoProportions['m_mIM']
                m_mIMS = mitoProportions['m_mIMS']
                m_mOM = mitoProportions['m_mOM']
                Anew = numpy.zeros((A.shape[0]+3, A.shape[1]))
                Anew[0:A.shape[0], :] = A
                OccupationMatrix.row_names += ['m_mIM', 'm_mIMS', 'm_mOM']
                OccupationMatrix.row_signs += ['E', 'E', 'E']
                OccupationMatrix.b = numpy.array(list(OccupationMatrix.b)+[float(0)]*3)
                Anew[OccupationMatrix.row_names.index(
                    'm_mIM'), OccupationMatrix.col_names.index('O_m')] = float(1)
                Anew[OccupationMatrix.row_names.index(
                    'm_mIMS'), OccupationMatrix.col_names.index('O_m')] = float(1)
                Anew[OccupationMatrix.row_names.index(
                    'm_mOM'), OccupationMatrix.col_names.index('O_m')] = float(1)
                Anew[OccupationMatrix.row_names.index(
                    'm_mIM'), OccupationMatrix.col_names.index('O_mIM')] = -m_mIM
                Anew[OccupationMatrix.row_names.index(
                    'm_mIMS'), OccupationMatrix.col_names.index('O_mIMS')] = -m_mIMS
                Anew[OccupationMatrix.row_names.index(
                    'm_mOM'), OccupationMatrix.col_names.index('O_mOM')] = -m_mOM

                OccupationMatrix.A = scipy.sparse.coo_matrix(Anew)
        else:
            OccupationMatrix.A = scipy.sparse.coo_matrix(A)

        OccupationMatrix.mapIndices()
        self.Problem.LP.addMatrix(matrix=OccupationMatrix)
        self.Problem.MuOneMatrix.addMatrix(matrix=OccupationMatrix)
        # {'Index':{'Param1':'+','Param2':'+','Param2':'-'}}

        #Type: 'Sum'#
        # {'Index':'Param1'}
        self.Problem.MuDependencies['FromParameters']['b'].update(
            {'n_density': {'Equation': '-nonenzymatic_proteins_n/inverse_average_protein_length', 'Variables': ['nonenzymatic_proteins_n', 'inverse_average_protein_length']}})
        self.Problem.MuDependencies['FromParameters']['b'].update({'mIM_density': {
                                                                  'Equation': '-nonenzymatic_proteins_mIM/inverse_average_protein_length', 'Variables': ['nonenzymatic_proteins_mIM', 'inverse_average_protein_length']}})
        self.Problem.MuDependencies['FromParameters']['b'].update({'vM_density': {
                                                                  'Equation': '-nonenzymatic_proteins_vM/inverse_average_protein_length', 'Variables': ['nonenzymatic_proteins_vM', 'inverse_average_protein_length']}})
        self.Problem.MuDependencies['FromParameters']['b'].update({'mIMS_density': {
                                                                  'Equation': '-nonenzymatic_proteins_mIMS/inverse_average_protein_length', 'Variables': ['nonenzymatic_proteins_mIMS', 'inverse_average_protein_length']}})
        self.Problem.MuDependencies['FromParameters']['b'].update(
            {'m_density': {'Equation': '-nonenzymatic_proteins_m/inverse_average_protein_length', 'Variables': ['nonenzymatic_proteins_m', 'inverse_average_protein_length']}})
        self.Problem.MuDependencies['FromParameters']['b'].update({'erM_density': {
                                                                  'Equation': '-nonenzymatic_proteins_erM/inverse_average_protein_length', 'Variables': ['nonenzymatic_proteins_erM', 'inverse_average_protein_length']}})
        self.Problem.MuDependencies['FromParameters']['b'].update({'mOM_density': {
                                                                  'Equation': '-nonenzymatic_proteins_mOM/inverse_average_protein_length', 'Variables': ['nonenzymatic_proteins_mOM', 'inverse_average_protein_length']}})
        self.Problem.MuDependencies['FromParameters']['b'].update(
            {'x_density': {'Equation': '-nonenzymatic_proteins_x/inverse_average_protein_length', 'Variables': ['nonenzymatic_proteins_x', 'inverse_average_protein_length']}})
        self.Problem.MuDependencies['FromParameters']['b'].update({'cM_density': {
                                                                  'Equation': '-nonenzymatic_proteins_cM/inverse_average_protein_length', 'Variables': ['nonenzymatic_proteins_cM', 'inverse_average_protein_length']}})
        self.Problem.MuDependencies['FromParameters']['b'].update({'gM_density': {
                                                                  'Equation': '-nonenzymatic_proteins_gM/inverse_average_protein_length', 'Variables': ['nonenzymatic_proteins_gM', 'inverse_average_protein_length']}})
        self.Problem.MuDependencies['FromParameters']['b'].update(
            {'c_density': {'Equation': '-nonenzymatic_proteins_c/inverse_average_protein_length', 'Variables': ['nonenzymatic_proteins_c', 'inverse_average_protein_length']}})
        self.Problem.MuDependencies['FromParameters']['b'].update({'O_total': {'Equation': '{} - nonenzymatic_proteins_Secreted/inverse_average_protein_length'.format(totalAA_parameter), 'Variables': [
                                                                  totalAA_parameter, 'nonenzymatic_proteins_Secreted', 'inverse_average_protein_length']}})

    # !!! deal with hardcoded parameter_names... !!!
    def estimate_specific_Kapps(self, proteomicsData, flux_bounds, mu, biomass_function=None, target_biomass_function=True, parsimonious_fba=True):
        """
        Parameters
        ----------
        proteomicsData : pandas.DataFrame (in mmol/gDW)
        flux_bounds : pandas.DataFrame  (in mmol/(gDW*h))
        mu : float (in 1/h)
        biomass_function : str
        target_biomass_function : bool
        atp_maintenance_to_biomassfunction : bool
        eukaryotic : bool
        """
        from scipy.stats.mstats import gmean

        old_model = copy.deepcopy(self.model)
        for i in self.model.targets.target_groups._elements_by_id['translation_targets'].concentrations._elements:
            if i.species == 'average_protein_c':
                new_agg = rba.xml.parameters.Aggregate(id_='total_protein', type_='multiplication')
                new_agg.function_references.append(rba.xml.parameters.FunctionReference(
                    function='amino_acid_concentration_total'))
                new_agg.function_references.append(rba.xml.parameters.FunctionReference(
                    function='inverse_average_protein_length'))
                self.model.parameters.aggregates._elements.append(new_agg)
                i.value = 'total_protein'
            else:
                self.model.targets.target_groups._elements_by_id['translation_targets'].concentrations._elements.remove(
                    i)
        for i in self.model.targets.target_groups._elements_by_id['transcription_targets'].concentrations._elements:
            if i.species == 'mrna':
                new_agg = rba.xml.parameters.Aggregate(id_='total_rna', type_='multiplication')
                new_agg.function_references.append(rba.xml.parameters.FunctionReference(
                    function='RNA_massfraction_CarbonLimitation'))
                new_agg.function_references.append(
                    rba.xml.parameters.FunctionReference(function='RNA_inversemillimolarweight'))
                self.model.parameters.aggregates._elements.append(new_agg)
                i.value = 'total_rna'
            else:
                self.model.targets.target_groups._elements_by_id['transcription_targets'].concentrations._elements.remove(
                    i)

        self.rebuild_from_model()
        self.setMedium(self.Medium)

        self.addExchangeReactions()
        self.setMu(mu)

        if target_biomass_function:
            self.buildFBA(objective='targets', maintenanceToBM=True)
            BMfunction = 'R_BIOMASS_targetsRBA'
        else:
            self.buildFBA(objective='classic', maintenanceToBM=False)
            BMfunction = biomass_function

        for j in [i for i in self.Medium.keys() if self.Medium[i] == 0]:
            Exrxn = 'R_EX_'+j.split('M_')[-1]+'_e'
            self.FBA.setUB({Exrxn: 0})

        rxn_LBs = {}
        rxn_UBs = {}
        for rx in flux_bounds['Reaction_ID']:
            lb = flux_bounds.loc[flux_bounds['Reaction_ID'] == rx, 'LB'].values[0]
            ub = flux_bounds.loc[flux_bounds['Reaction_ID'] == rx, 'UB'].values[0]
            if not pandas.isna(lb):
                rxn_LBs.update({rx: lb})
            if not pandas.isna(ub):
                rxn_UBs.update({rx: ub})
        self.FBA.setLB(rxn_LBs)
        self.FBA.setUB(rxn_UBs)

        self.FBA.clearObjective()
        self.FBA.setObjectiveCoefficients({BMfunction: -1})
        self.FBA.solveLP(feasibleStatuses=[1, 2, 3, 5, 6])
        BMfluxOld = self.FBA.SolutionValues[BMfunction]

        if parsimonious_fba:
            self.FBA.parsimonise()
            self.FBA.setLB(rxn_LBs)
            self.FBA.setUB(rxn_UBs)
            self.FBA.setLB({BMfunction: BMfluxOld})
            self.FBA.setUB({BMfunction: BMfluxOld})
            self.FBA.solveLP(feasibleStatuses=[1, 2, 3, 5, 6])

        FluxDistribution = pandas.DataFrame(index=list(
            self.FBA.SolutionValues.keys()), columns=['FluxValues'])
        FluxDistribution['FluxValues'] = list(self.FBA.SolutionValues.values())
        BMfluxNew = self.FBA.SolutionValues[BMfunction]

        ProtoIDmap = {}
        for i in self.ModelStructure.ProteinInfo.Elements.keys():
            ProtoID = self.ModelStructure.ProteinInfo.Elements[i]['ProtoID']
            if ProtoID in list(proteomicsData['ID']):
                if not pandas.isna(proteomicsData.loc[proteomicsData['ID'] == ProtoID, 'copy_number'].values[0]):
                    if proteomicsData.loc[proteomicsData['ID'] == ProtoID, 'copy_number'].values[0] != 0:
                        if ProtoID in ProtoIDmap.keys():
                            ProtoIDmap[ProtoID]['ModelProteins'].append(i)
                        else:
                            ProtoIDmap.update(
                                {ProtoID: {'ModelProteins': [i], 'CopyNumber': proteomicsData.loc[proteomicsData['ID'] == ProtoID, 'copy_number'].values[0]}})

        ReactionMap = {}
        for i in self.ModelStructure.ReactionInfo.Elements.keys():
            if '_duplicate_' in i:
                continue
            else:
                if i in list(FluxDistribution.index):
                    if FluxDistribution.loc[i, 'FluxValues'] != 0:
                        ReactionMap.update({i: {'ModelReactions': list(
                            [i]+self.ModelStructure.ReactionInfo.Elements[i]['Twins']), 'Flux': FluxDistribution.loc[i, 'FluxValues']}})

        IsoReaction2ProtoReaction = {}
        for i in ReactionMap.keys():
            for j in ReactionMap[i]['ModelReactions']:
                IsoReaction2ProtoReaction[j] = i

        EnzymeMap = {}
        for i in self.ModelStructure.EnzymeInfo.Elements.keys():
            if self.ModelStructure.EnzymeInfo.Elements[i]['Reaction'] in IsoReaction2ProtoReaction:
                CompositionDict = {self.ModelStructure.ProteinInfo.Elements[j]['ProtoID']: self.ModelStructure.EnzymeInfo.Elements[
                    i]['Subunits'][j] for j in self.ModelStructure.EnzymeInfo.Elements[i]['Subunits'].keys()}
                ProtoReaction = IsoReaction2ProtoReaction[self.ModelStructure.EnzymeInfo.Elements[i]['Reaction']]
                CopyNumbers = []
                Stoichiometries = []
                EnzymeNumbers = []
                for j in CompositionDict.keys():
                    if j in ProtoIDmap.keys():
                        CopyNumbers.append(ProtoIDmap[j]['CopyNumber'])
                        Stoichiometries.append(CompositionDict[j])
                        EnzymeNumbers.append(ProtoIDmap[j]['CopyNumber']/CompositionDict[j])
                GM_enzymenumber = 0
                if len(EnzymeNumbers) > 0:
                    GM_enzymenumber = gmean(numpy.array(EnzymeNumbers))
                EnzymeMap.update(
                    {i: {'ProtoReaction': ProtoReaction, 'EnzymeNumber': GM_enzymenumber}})

        EnzymeMap2 = {}
        for i in ReactionMap.keys():
            totalIsoEnzymeNumber = 0
            for j in ReactionMap[i]['ModelReactions']:
                respectiveEnzyme = self.ModelStructure.ReactionInfo.Elements[j]['Enzyme']
                if respectiveEnzyme in EnzymeMap.keys():
                    totalIsoEnzymeNumber += EnzymeMap[respectiveEnzyme]['EnzymeNumber']
            for j in ReactionMap[i]['ModelReactions']:
                respectiveEnzyme = self.ModelStructure.ReactionInfo.Elements[j]['Enzyme']
                if respectiveEnzyme in EnzymeMap.keys():
                    concentration = EnzymeMap[respectiveEnzyme]['EnzymeNumber']
                    if concentration != 0:
                        if numpy.isfinite(concentration):
                            specificFlux = ReactionMap[i]['Flux'] * \
                                EnzymeMap[respectiveEnzyme]['EnzymeNumber']/totalIsoEnzymeNumber
                            EnzymeMap2.update({respectiveEnzyme: {'CopyNumber': EnzymeMap[respectiveEnzyme]['EnzymeNumber'],
                                                                  'Concentration': concentration, 'Flux': specificFlux, 'Kapp': abs(specificFlux/concentration)}})

        self.model = old_model
        self.rebuild_from_model()
        self.setMedium(self.Medium)

        out = pandas.DataFrame()
        for i in EnzymeMap2.keys():
            # if EnzymeMap2[i]['CopyNumber'] == 0:
            #    continue
            out.loc[i, 'Enzyme_ID'] = i
            out.loc[i, 'CopyNumber'] = EnzymeMap2[i]['CopyNumber']
            out.loc[i, 'Concentration'] = EnzymeMap2[i]['Concentration']
            out.loc[i, 'Flux'] = EnzymeMap2[i]['Flux']
            out.loc[i, 'Kapp'] = EnzymeMap2[i]['Kapp']

        return(out)

    def estimate_default_Kapps(self, target_mu, compartment_densities_and_PGs=None, flux_bounds=None, plateau_limit=4, mu_approximation_precision=0.0001, transporter_to_lumen_coefficient=10, default_kapp_LB=0, default_kapp_UB=1000000, start_val=None, densities_to_fix=None, eukaryotic=False):
        """
        Parameters
        ----------
        target_mu : float
        compartment_densities_and_PGs : pandas.DataFrame
        flux_bounds : pandas.DataFrame
        """
        old_model = copy.deepcopy(self.model)

        orig_enz = self.model.parameters.functions._elements_by_id[
            'default_efficiency'].parameters._elements_by_id['CONSTANT'].value

        out = pandas.DataFrame()
        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
            self.model.parameters.functions._elements_by_id[str(
                'fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density']
            self.model.parameters.functions._elements_by_id[str(
                'fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction']
        self.rebuild_from_model()
        self.addExchangeReactions()
        self.setMedium(self.Medium)
        if eukaryotic:
            self.eukaryoticDensities_calibration(CompartmentRelationships=False)
        else:
            if densities_to_fix is None:
                comp_density_rows = self.Problem.CompartmentDensities
            else:
                comp_density_rows = densities_to_fix
            self.Problem.setConstraintType(
                dict(zip(comp_density_rows, ['E']*len(comp_density_rows))))
        # self.Problem.setConstraintType(
        #    dict(zip(['mOM_density', 'mIMS_density', 'mIM_density', 'm_density'], ['L'])))

        rxn_LBs = {}
        rxn_UBs = {}
        for rx in flux_bounds['Reaction_ID']:
            lb = flux_bounds.loc[flux_bounds['Reaction_ID'] == rx, 'LB'].values[0]
            ub = flux_bounds.loc[flux_bounds['Reaction_ID'] == rx, 'UB'].values[0]
            if not pandas.isna(lb):
                rxn_LBs.update({rx: lb})
            if not pandas.isna(ub):
                rxn_UBs.update({rx: ub})
        self.Problem.setLB(rxn_LBs)
        self.Problem.setUB(rxn_UBs)

        kapp_LB = default_kapp_LB
        if default_kapp_UB is not None:
            kapp_UB = default_kapp_UB
        else:
            kapp_UB = orig_enz*1000
        #new_kapp = (kapp_UB+kapp_LB)/2

        if start_val is not None:
            new_kapp = start_val
        else:
            new_kapp = orig_enz

        Mu_pred = self.findMaxGrowthRate(precision=0.0001, max=1)

        Mus = []
        Mus_Error = []
        Kapps = []
        last_Mu = numpy.nan
        plateau_count = 0
        if abs(target_mu - Mu_pred) > mu_approximation_precision:
            while abs(target_mu - Mu_pred) > mu_approximation_precision:
                if plateau_count >= plateau_limit:
                    break
                self.model.parameters.functions._elements_by_id[
                    'default_efficiency'].parameters._elements_by_id['CONSTANT'].value = new_kapp
                self.model.parameters.functions._elements_by_id['default_transporter_efficiency'].parameters._elements_by_id[
                    'CONSTANT'].value = transporter_to_lumen_coefficient*new_kapp
                self.rebuild_from_model()
                self.addExchangeReactions()
                self.setMedium(self.Medium)
                self.Problem.setLB(rxn_LBs)
                self.Problem.setUB(rxn_UBs)
                if eukaryotic:
                    self.eukaryoticDensities_calibration(CompartmentRelationships=False)
                else:
                    self.Problem.setConstraintType(
                        dict(zip(comp_density_rows, ['E']*len(comp_density_rows))))

                Mu_pred = self.findMaxGrowthRate()
                Mus_Error.append(abs(target_mu - Mu_pred))
                Mus.append(Mu_pred)
                Kapps.append(new_kapp)
                if Mu_pred > target_mu:
                    new_kapp_prelim = kapp_LB+(0.5*abs(kapp_LB-new_kapp))
                    kapp_UB = new_kapp
                elif Mu_pred < target_mu:
                    new_kapp_prelim = kapp_UB-(0.5*abs(new_kapp-kapp_UB))
                    kapp_LB = new_kapp
                new_kapp = new_kapp_prelim
                if len(Mus) > 2:
                    if Mus[-2] == Mu_pred:
                        plateau_count += 1
                    else:
                        plateau_count = 0
        else:
            Mus.append(Mu_pred)
            Mus_Error.append(abs(target_mu - Mu_pred))
            Kapps.append(
                self.model.parameters.functions._elements_by_id['default_efficiency'].parameters._elements_by_id['CONSTANT'].value)
        #self.model = old_model
        self.rebuild_from_model()
        self.setMedium(self.Medium)
        out = pandas.DataFrame()
        out['Mu'] = Mus
        out['delta_Mu'] = Mus_Error
        out['default_efficiency'] = Kapps
        out['default_transporter_efficiency'] = [transporter_to_lumen_coefficient*i for i in Kapps]
        return(out)

    def inject_default_kapps(self, default_kapp, default_transporter_kapp):
        if numpy.isfinite(default_kapp):
            self.model.parameters.functions._elements_by_id[
                'default_efficiency'].parameters._elements_by_id['CONSTANT'].value = default_kapp
        if numpy.isfinite(default_transporter_kapp):
            self.model.parameters.functions._elements_by_id[
                'default_transporter_efficiency'].parameters._elements_by_id['CONSTANT'].value = default_transporter_kapp
        self.rebuild_from_model()

    def inject_process_capacities(self, process_efficiencies):
        """
        Parameters
        ----------
        process_efficiencies : pandas.DataFrame(columns=['Process','Parameter','Value'])
        """
        for i in process_efficiencies.index:
            if numpy.isfinite(process_efficiencies.loc[i, 'Value']):
                if process_efficiencies.loc[i, 'Process'] in self.model.processes.processes._elements_by_id.keys():
                    if not pandas.isna(process_efficiencies.loc[i, 'Value']):
                        self.model.processes.processes._elements_by_id[process_efficiencies.loc[i,
                                                                                                'Process']].machinery.capacity.value = process_efficiencies.loc[i, 'Parameter']
                        const = rba.xml.parameters.Function(process_efficiencies.loc[i, 'Parameter'], 'constant', parameters={
                                                            'CONSTANT': process_efficiencies.loc[i, 'Value']}, variable=None)
                        if process_efficiencies.loc[i, 'Parameter'] not in self.model.parameters.functions._elements_by_id.keys():
                            self.model.parameters.functions.append(const)
                        else:
                            self.model.parameters.functions._elements_by_id[const.id].parameters._elements_by_id[
                                'CONSTANT'].value = process_efficiencies.loc[i, 'Value']
        self.rebuild_from_model()

    def inject_specific_kapps(self, specific_kapps, round_to_digits=0):
        """
        Parameters
        ----------
        specific_kapps : pandas.DataFrame
        """
        parameterized = []
        if 'Enzyme_ID' in list(specific_kapps.columns):
            for enz in list(specific_kapps['Enzyme_ID']):
                if not pandas.isna(specific_kapps.loc[specific_kapps['Enzyme_ID'] == enz, 'Kapp'].values[0]):
                    if numpy.isfinite(specific_kapps.loc[specific_kapps['Enzyme_ID'] == enz, 'Kapp'].values[0]):
                        if enz not in parameterized:
                            all_enzs = self.ModelStructure.EnzymeInfo.Elements[enz]['Isozymes']
                            all_enzs.append(enz)
                            parameterized += all_enzs
                            if len(all_enzs) == 1:
                                proto_enz = all_enzs[0]
                            else:
                                proto_enz = [i for i in all_enzs if not '_duplicate_' in i][0]
                            val = round(specific_kapps.loc[specific_kapps['Enzyme_ID']
                                                           == enz, 'Kapp'].values[0], round_to_digits)
                            const = rba.xml.parameters.Function(
                                str(proto_enz + '_kapp__constant'), 'constant', parameters={'CONSTANT': val}, variable=None)
                            if str(proto_enz + '_kapp__constant') not in self.model.parameters.functions._elements_by_id.keys():
                                self.model.parameters.functions.append(const)
                            else:
                                #self.model.parameters.functions._elements_by_id[const.id] = const
                                self.model.parameters.functions._elements_by_id[
                                    const.id].parameters._elements_by_id['CONSTANT'].value = val
                            count = 0
                            #self.model.parameters.functions._elements_by_id['default_efficiency'].parameters._elements_by_id['CONSTANT'].value = default_kapp
                            for e in self.model.enzymes.enzymes:
                                if e.id in all_enzs:
                                    count += 1
                                    e.forward_efficiency = str(proto_enz + '_kapp__constant')
                                    e.backward_efficiency = str(proto_enz + '_kapp__constant')
                                    if count == len(all_enzs):
                                        break
            self.rebuild_from_model()

    def get_parameter_values(self, parameter_type, species=None, output_format='dict'):
        if parameter_type == 'medium_composition':
            if species is None:
                results = self.Medium
            elif type(species) is str:
                results = {species: self.Medium[species]}
            elif type(species) is list:
                results = {sp: self.Medium[sp] for sp in species}

        elif parameter_type == 'machine_efficiencies':
            if species is None:
                parameter_names = {process_name: self.model.processes.processes._elements_by_id[self.ModelStructure.ProcessInfo.Elements[
                    process_name]['ID']].machinery.capacity.value for process_name in self.ModelStructure.ProcessInfo.Elements.keys()}
            elif type(species) is str:
                parameter_names = {
                    species: self.model.processes.processes._elements_by_id[self.ModelStructure.ProcessInfo.Elements[species]['ID']].machinery.capacity.value}
            elif type(species) is list:
                parameter_names = {
                    sp: self.model.processes.processes._elements_by_id[self.ModelStructure.ProcessInfo.Elements[sp]['ID']].machinery.capacity.value for sp in species}
            results = {pn: self.get_parameter_value(
                parameter=parameter_names[pn]) for pn in parameter_names}

        elif parameter_type == 'enzyme_efficiencies' or parameter_type == 'enzyme_efficiencies_forward' or parameter_type == 'enzyme_efficiencies_backward':
            if species is None:
                parameter_names = {enzyme_name: {'Forward': self.model.enzymes.enzymes._elements_by_id[enzyme_name].forward_efficiency, 'Backward': self.model.enzymes.enzymes._elements_by_id[
                    enzyme_name].backward_efficiency} for enzyme_name in self.ModelStructure.EnzymeInfo.Elements.keys()}
            elif type(species) is str:
                parameter_names = {species: {'Forward': self.model.enzymes.enzymes._elements_by_id[
                    species].forward_efficiency, 'Backward': self.model.enzymes.enzymes._elements_by_id[species].backward_efficiency}}
            elif type(species) is list:
                parameter_names = {enzyme_name: {'Forward': self.model.enzymes.enzymes._elements_by_id[enzyme_name].forward_efficiency,
                                                 'Backward': self.model.enzymes.enzymes._elements_by_id[enzyme_name].backward_efficiency} for enzyme_name in species}
            if parameter_type == 'enzyme_efficiencies':
                results = {pn: {'Forward': self.get_parameter_value(parameter=parameter_names[pn]['Forward']), 'Backward': self.get_parameter_value(
                    parameter=parameter_names[pn]['Backward'])} for pn in parameter_names.keys()}
            elif parameter_type == 'enzyme_efficiencies_forward':
                results = {pn: self.get_parameter_value(
                    parameter=parameter_names[pn]['Forward']) for pn in parameter_names.keys()}
            elif parameter_type == 'enzyme_efficiencies_backward':
                results = {pn: self.get_parameter_value(
                    parameter=parameter_names[pn]['Backward']) for pn in parameter_names.keys()}

        elif parameter_type == 'maximal_densities':
            density_dict = {i.compartment: self.get_parameter_value(
                parameter=i.upper_bound) for i in self.model.density.target_densities}
            if species is None:
                results = density_dict
            elif type(species) is str:
                results = {species: density_dict[species]
                           for sp in [species] if sp in density_dict.keys()}
            elif type(species) is list:
                results = {sp: density_dict[sp] for sp in species if sp in density_dict.keys()}

        elif parameter_type == 'target_values':
            target_dict = {self.ModelStructure.TargetInfo.Elements[target_ID]['TargetEntity']: {'Target_id': target_ID, 'Target_value': self.get_parameter_value(
                parameter=self.ModelStructure.TargetInfo.Elements[target_ID]['TargetValue'])} for target_ID in self.ModelStructure.TargetInfo.Elements.keys()}
            if species is None:
                results = target_dict
            elif type(species) is str:
                results = {species: target_dict[species]
                           for sp in [species] if sp in target_dict.keys()}
            elif type(species) is list:
                results = {sp: target_dict[sp] for sp in species if sp in target_dict.keys()}

        if output_format == 'dict':
            return(results)
        if output_format == 'json':
            return(json.dumps(results))


    def get_parameter_definition(self, parameter):
        if parameter in self.model.parameters.functions._elements_by_id.keys():
            function = self.model.parameters.functions._elements_by_id[parameter]
            expression = parse_function(function)
        elif parameter in self.model.parameters.aggregates._elements_by_id.keys():
            function_id_list = get_function_list_of_aggregate(
                aggregate=self.model.parameters.aggregates._elements_by_id[parameter])
            expression = parse_aggregate(aggregate=self.model.parameters.aggregates._elements_by_id[parameter], function_list=[
                                         self.model.parameters.functions._elements_by_id[f_id] for f_id in function_id_list])
        else:
            return({})
        return(expression)

    def get_parameter_value(self, parameter):
        if parameter in self.model.parameters.functions._elements_by_id.keys():
            function = self.model.parameters.functions._elements_by_id[parameter]
            expression = parse_function_with_parameter_values(function)
        elif parameter in self.model.parameters.aggregates._elements_by_id.keys():
            function_id_list = get_function_list_of_aggregate(
                aggregate=self.model.parameters.aggregates._elements_by_id[parameter])
            expression = parse_aggregate_with_parameter_values(aggregate=self.model.parameters.aggregates._elements_by_id[parameter], function_list=[
                self.model.parameters.functions._elements_by_id[f_id] for f_id in function_id_list])
        else:
            return({parameter: numpy.nan})
        variable_values = {}
        for v in expression[parameter]['Variables']:
            if v == 'growth_rate':
                variable_values[v] = self.Mu
            elif v in self.Medium.keys():
                variable_values[v] = self.Medium[v]
            elif v.endswith('_e'):
                if v[:-2] in self.Medium.keys():
                    variable_values[v] = self.Medium[v[:-2]]
            else:
                variable_values = {}
                return({parameter: numpy.nan})
        result = evaluate_expression(expression_dictionary=expression,
                                     variable_values=variable_values)
        return(result)


def get_parameter_value_from_model(function, parameter_ID):
    return(function.parameters._elements_by_id[parameter_ID].value)


def make_paramter_function_specific(function_ID, parameter, return_normal=False):
    if return_normal:
        return(str(parameter))
    else:
        return(str('{}__parameter__{}'.format(function_ID, parameter)))


def parse_function(function):
    independent_variable = function.variable
    function_ID = function.id
    if function.type == 'constant':
        eq = make_paramter_function_specific(
            function_ID=function_ID, parameter='CONSTANT', return_normal=True)
        latex_string = str(make_paramter_function_specific(
            function_ID=function_ID, parameter='CONSTANT', return_normal=True))
        function_parameter_values = {'CONSTANT': get_parameter_value_from_model(
            function=function, parameter_ID='CONSTANT')}   

    elif function.type == 'exponential':
        eq = 'e**({}*{})'.format(make_paramter_function_specific(function_ID=function_ID,
                                                                    parameter='RATE', return_normal=True), str(independent_variable))
        latex_string = str('e^{'+str(make_paramter_function_specific(function_ID=function_ID,
                                                                    parameter='RATE', return_normal=True)) + ' '+str(independent_variable)+'}')
        function_parameter_values = {'RATE': get_parameter_value_from_model(
            function=function, parameter_ID='RATE')}

    elif function.type == 'linear':
        eq = str('{}+{}*{}'.format(make_paramter_function_specific(function_ID=function_ID, parameter='LINEAR_CONSTANT', return_normal=True),
                                   make_paramter_function_specific(function_ID=function_ID, parameter='LINEAR_COEF', return_normal=True), str(independent_variable)))
        latex_string = str(make_paramter_function_specific(function_ID=function_ID, parameter='LINEAR_CONSTANT', return_normal=True) +
                           make_paramter_function_specific(function_ID=function_ID, parameter='LINEAR_COEF', return_normal=True)+' '+str(independent_variable))
        function_parameter_values = {'LINEAR_CONSTANT': get_parameter_value_from_model(function=function, parameter_ID='LINEAR_CONSTANT'),
                                     'LINEAR_COEF': get_parameter_value_from_model(function=function, parameter_ID='LINEAR_COEF'),
                                     'X_MIN': get_parameter_value_from_model(function=function, parameter_ID='X_MIN'),
                                     'X_MAX': get_parameter_value_from_model(function=function, parameter_ID='X_MAX'),
                                     'Y_MIN': get_parameter_value_from_model(function=function, parameter_ID='Y_MIN'),
                                     'Y_MAX': get_parameter_value_from_model(function=function, parameter_ID='Y_MAX'), }

    elif function.type == 'michaelisMenten':
        eq = str('{}*{}/({}+{})'.format(make_paramter_function_specific(function_ID=function_ID, parameter='kmax', return_normal=True),
                                        str(independent_variable), str(independent_variable), make_paramter_function_specific(function_ID=function_ID, parameter='Km', return_normal=True)))
        function_parameter_values = {'kmax': get_parameter_value_from_model(function=function, parameter_ID='kmax'),
                                     'Km': get_parameter_value_from_model(function=function, parameter_ID='Km'),
                                     'Y_MIN': get_parameter_value_from_model(function=function, parameter_ID='Y_MIN')}

    return({function_ID: {'Type': function.type, 'Equation': eq, 'Variables': [str(independent_variable)], 'Function_parameters': function_parameter_values}})



def parse_function_with_parameter_values(function):
    independent_variable = function.variable
    function_ID = function.id
    if function.type == 'constant':
        return({function_ID: {'Equation': '{}'.format(str(get_parameter_value_from_model(function=function, parameter_ID='CONSTANT'))), 'Variables': []}})

    elif function.type == 'exponential':
        return({function_ID: {'Equation': '{}**({}*{})'.format(str(numpy.e), str(get_parameter_value_from_model(function=function, parameter_ID='RATE')), str(independent_variable)), 'Variables': [str(independent_variable)]}})

    elif function.type == 'linear':
        return({function_ID: {'Equation': str('{}+{}*{}'.format(str(get_parameter_value_from_model(function=function, parameter_ID='LINEAR_CONSTANT')), str(get_parameter_value_from_model(function=function, parameter_ID='LINEAR_COEF')), str(independent_variable))), 'Variables': [str(independent_variable)]}})

    elif function.type == 'michaelisMenten':
        return({function_ID: {'Equation': str('{}*{}/({}+{})'.format(str(get_parameter_value_from_model(function=function, parameter_ID='kmax')), str(independent_variable), str(get_parameter_value_from_model(function=function, parameter_ID='Km')), str(independent_variable))), 'Variables': [str(independent_variable)]}})        


def get_parameter_of_function(function, parameter):
    return(function.parameters._elements_by_id[parameter])


def join_functions_multiplicatively(parsed_function_list):
    term_list = []
    variable_list = []
    for function in parsed_function_list:
        function_ID = list(function.keys())[0]
        term_list.append(str('('+function[function_ID]['Equation']+')'))
        variable_list += function[function_ID]['Variables']
    return({'Type': 'Aggregate', 'Equation': '*'.join(term_list), 'Variables': list(set(variable_list))})

def get_function_list_of_aggregate(aggregate):
    return([agg.function for agg in aggregate.function_references._elements])


def parse_aggregate_with_parameter_values(aggregate, function_list):
    aggregate_ID = aggregate.id
    if aggregate.type == 'multiplication':
        parsed_function_list = [parse_function_with_parameter_values(
            function) for function in function_list]
        return({aggregate_ID: join_functions_multiplicatively(parsed_function_list=parsed_function_list)})
    else:
        return({aggregate_ID: {'Equation': '', 'Variables': []}})


def parse_aggregate(aggregate, function_list):
    aggregate_ID = aggregate.id
    if aggregate.type == 'multiplication':
        parsed_function_list = [parse_function(function) for function in function_list]
        result = {aggregate_ID: join_functions_multiplicatively(
            parsed_function_list=parsed_function_list)}
        result[aggregate_ID]['Multiplicative Terms'] = [f.id for f in function_list]
        return(result)
    else:
        return({aggregate_ID: {'Type': 'Aggregate', 'Equation': '', 'Variables': [], 'Multiplicative Terms': []}})

# def transform_to_latex(equation):
#

def MediumDependentCoefficients_A(Controler):
    out = {}
    MedDepRxns = [list(i.keys()) for i in list(Controler.ExchangeMap.values())]
    MedDepRxnsFlatted = list(set([item for sublist in MedDepRxns for item in sublist]))
    for i in Controler.ModelStructure.EnzymeConstraintsInfo.Elements.keys():
        if Controler.ModelStructure.EnzymeConstraintsInfo.Elements[i]['AssociatedReaction'] in MedDepRxnsFlatted:
            nonConst = False
            for j in Controler.ModelStructure.EnzymeConstraintsInfo.Elements[i]['CapacityParameter']:
                if list(j.values())[0]['FunctionType'] != 'constant':
                    nonConst = True
            if nonConst:
                if Controler.ModelStructure.EnzymeConstraintsInfo.Elements[i]['AssociatedReaction'] in list(out.keys()):
                    out[Controler.ModelStructure.EnzymeConstraintsInfo.Elements[i]
                        ['AssociatedReaction']].append(i)
                else:
                    out.update(
                        {Controler.ModelStructure.EnzymeConstraintsInfo.Elements[i]['AssociatedReaction']: [i]})
    return([(out[i][0], Controler.ModelStructure.ReactionInfo.Elements[i]['Enzyme'])for i in out.keys()])


def buildCompositionofUnusedProtein(species, composition):
    products = {}
    reactants = {}
    for j in composition.keys():
        stoch = composition[j]
        if len(list(species[j]['Products'].keys())) > 0:
            for k in species[j]['Products'].keys():
                if k in products:
                    products[k] += species[j]['Products'][k]*stoch
                else:
                    products.update({k: species[j]['Products'][k]*stoch})
        if len(list(species[j]['Reactants'].keys())) > 0:
            for k in species[j]['Reactants'].keys():
                if k in reactants:
                    reactants[k] += species[j]['Reactants'][k]*stoch
                else:
                    reactants.update({k: species[j]['Reactants'][k]*stoch})
    uniqueMets = list(set(list(products.keys())+list(reactants.keys())))
    NetMets = {}
    for j in uniqueMets:
        if j in list(products.keys()):
            produced = products[j]
        elif j not in list(products.keys()):
            produced = 0
        if j in list(reactants.keys()):
            consumed = reactants[j]
        elif j not in list(reactants.keys()):
            consumed = 0
        NetStoch = produced-consumed
        NetMets.update({j: float(NetStoch)})
    return(NetMets)


def buildUsedProteinConstraint(Controler, protein):
    ProteinConstraint = numpy.zeros(Controler.Problem.LP.A.shape[1])
    ConsumersEnzymes = Controler.ModelStructure.ProteinInfo.Elements[protein]['associatedEnzymes']
    ConsumersProcess = Controler.ModelStructure.ProteinInfo.Elements[protein]['SupportsProcess']
    for j in ConsumersEnzymes:
        StochFactor = Controler.ModelStructure.EnzymeInfo.Elements[j]['Subunits'][protein]['StochFac']
        LikeliestVarName = difflib.get_close_matches(j, Controler.Problem.LP.col_names, 1)[0]
        ColIndex = Controler.Problem.LP.colIndicesMap[LikeliestVarName]
        ProteinConstraint[ColIndex] = float(StochFactor)
    for j in ConsumersProcess:
        LikeliestVarName = difflib.get_close_matches(str(
            Controler.ModelStructure.ProcessInfo.Elements[j]['ID']+'_machinery'), Controler.Problem.LP.col_names, 1)[0]
        StochFactor = Controler.ModelStructure.ProcessInfo.Elements[j]['Composition'][protein]
        ColIndex = Controler.Problem.LP.colIndicesMap[LikeliestVarName]
        ProteinConstraint[ColIndex] = float(StochFactor)
    return(ProteinConstraint)


def QualitativeMediumChange(Controller, changes, species):
    QualitativeMediumChange = False
    if float(Controller.Medium[species]) == float(0):
        if float(changes[species]) != float(0):
            boundValue = 1000.0
            QualitativeMediumChange = True
        else:
            return([QualitativeMediumChange])
    if float(Controller.Medium[species]) != float(0):
        if float(changes[species]) == float(0):
            boundValue = 0.0
            QualitativeMediumChange = True
        else:
            return([QualitativeMediumChange])
    return([QualitativeMediumChange, float(boundValue)])


def findExchangeReactions(Controller, species):
    Reactions = list(Controller.ExchangeMap[species].keys())
    exchanges = {}
    for i in Reactions:
        if Controller.ExchangeMap[species][i] > 0:
            exchanges.update({i: 'Product'})
        elif Controller.ExchangeMap[species][i] < 0:
            exchanges.update({i: 'Reactant'})
    return(exchanges)


def determineCoefficient(x, changes, species):
    multiplicativeFactors = []
    for k in x:
        result = 1
        type = list(k.values())[0]['FunctionType']
        pars = list(k.values())[0]['FunctionParameters']
        if type == 'constant':
            result = numpy.float64(pars['C'])
        if type == 'exponential':
            L = 1
            if 'Lambda' in list(pars.keys()):
                L = numpy.float64(pars['Lambda'])
            result = numpy.exp(float(changes[species])*L)
        if type == 'indicator':
            maxi = numpy.inf
            mini = -numpy.inf
            if 'xMax' in list(pars.keys()):
                maxi = numpy.float64(pars['xMax'])
            if 'xMin' in list(pars.keys()):
                mini = numpy.float64(pars['xMin'])
            result = (float(changes[species]) > mini) and (float(changes[species]) < maxi)
        if type == 'linear':
            X_maxi = numpy.inf
            X_mini = -numpy.inf
            Y_maxi = numpy.inf
            Y_mini = -numpy.inf
            A = 1
            C = 0
            if 'A' in list(pars.keys()):
                A = numpy.float64(pars['A'])
            if 'C' in list(pars.keys()):
                C = numpy.float64(pars['C'])
            if 'xMin' in list(pars.keys()):
                X_mini = numpy.float64(pars['xMin'])
            if 'xMax' in list(pars.keys()):
                X_maxi = numpy.float64(pars['xMax'])
            if 'yMin' in list(pars.keys()):
                Y_mini = numpy.float64(pars['yMin'])
            if 'yMax' in list(pars.keys()):
                Y_maxi = numpy.float64(pars['yMax'])
            X = float(changes[species])
            if float(changes[species]) < X_mini:
                X = X_mini
            if float(changes[species]) > X_maxi:
                X = X_maxi
            Y = A*X + C
            result = Y
            if Y < Y_mini:
                result = Y_mini
            if Y > Y_maxi:
                result = Y_maxi
        if type == 'michaelisMenten':
            Y_mini = -numpy.inf
            KM = 0
            VM = 1
            if 'Km' in list(pars.keys()):
                KM = numpy.float64(pars['Km'])
            if 'Vmax' in list(pars.keys()):
                VM = numpy.float64(pars['Vmax'])
            if 'yMin' in list(pars.keys()):
                Y_mini = numpy.float64(pars['yMin'])
            Y = VM*float(changes[species])/(float(changes[species])+KM)
            result = Y
            if Y < Y_mini:
                result = Y_mini
        if type == 'competitiveInhibition':
            Y_mini = -numpy.inf
            KM = 0
            VM = 1
            KI = 0
            I = 0
            if 'Ki' in list(pars.keys()):
                KI = numpy.float64(pars['Ki'])
            if 'I' in list(pars.keys()):
                I = numpy.float64(pars['I'])
            if 'Km' in list(pars.keys()):
                KM = numpy.float64(pars['Km'])
            if 'Vmax' in list(pars.keys()):
                VM = numpy.float64(pars['Vmax'])
            if 'yMin' in list(pars.keys()):
                Y_mini = numpy.float64(pars['yMin'])
            Y = VM*float(changes[species])/(float(changes[species])+KM*(1+I/KI))
            result = Y
            if Y < Y_mini:
                result = Y_mini
        if type == 'inverse':
            C = 1
            if 'C' in list(pars.keys()):
                C = numpy.float64(pars['C'])
            result = 1
            if float(changes[species]) is not 0:
                result = C/float(changes[i])
        multiplicativeFactors.append(result)
        value = numpy.prod(numpy.array(multiplicativeFactors))
    return(float(value))


def ProtoProteomeRecording(Controller, run, Proteinlevels):
    out = {}
    for i in list(Controller.ModelStructure.ProteinGeneMatrix['ProtoProteins']):
        row_ind = list(Controller.ModelStructure.ProteinGeneMatrix['ProtoProteins']).index(i)
        respective_row = Controller.ModelStructure.ProteinGeneMatrix['Matrix'][row_ind]
        nonZero = list(numpy.where(respective_row != 0)[0])
        level = 0
        for j in nonZero:
            id = Controller.ModelStructure.ProteinGeneMatrix['Proteins'][j]
            level += Proteinlevels.loc[id, run]
        out.update({i: level})
    return(out)


def ProteomeRecording(Controller, run):

    EnzDF = pandas.DataFrame(index=Controller.Problem.Enzymes)
    PrcDF = pandas.DataFrame(index=Controller.Problem.Processes)
    EnzDF[run] = [Controller.Problem.SolutionValues[i]for i in Controller.Problem.Enzymes]
    PrcDF[run] = [Controller.Problem.SolutionValues[i]for i in Controller.Problem.Processes]

    ProteinProteinMatrix = numpy.array(
        Controller.ModelStructure.ProteinMatrix['Matrix']).astype(numpy.float64)
    C = Controller.ModelStructure.ProteinMatrix['Consumers']
    Consumers = []
    for i in C:
        if i.startswith('P_'):
            # Consumers.append(str(i+'_machinery'))
            Consumers.append(str(i))
        if not i.startswith('P_'):
            Consumers.append(i)
    Proteins = Controller.ModelStructure.ProteinMatrix['Proteins']
    DF = pandas.concat([EnzDF, PrcDF], axis=0)
    ProteinLevels = pandas.DataFrame(index=Proteins)
    vector = numpy.nan_to_num(DF[run].reindex(Consumers))
    Level = ProteinProteinMatrix.dot(vector)
    ProteinLevels[run] = Level
    addedProts = [col for col in Controller.Problem.LP.col_names if col.startswith('TotalLevel_')]
    if len(addedProts) > 0:
        for p in addedProts:
            protID = p.split('TotalLevel_')[1]
            ProteinLevels[run].loc[protID] = Controller.Problem.SolutionValues[p]
    return(dict(zip(list(ProteinLevels.index), list(ProteinLevels[run]))))

def mapIsoReactions(Controller):
    if hasattr(Controller, 'Results'):
        out = pandas.DataFrame()
        for run in list(Controller.Results['Reactions'].columns):
            rf = dict(zip(list(Controller.Results['Reactions'].index), list(
                Controller.Results['Reactions'][run])))
            rf = {k: v for k, v in rf.items() if v != 0.}
            rf_merged = collections.defaultdict(float)
            for reac_id, flux_val in rf.items():
                if "duplicate" in reac_id:
                    last_idx = reac_id.index('duplicate') - 1
                    rf_merged[reac_id[:last_idx]] += flux_val
                else:
                    rf_merged[reac_id] += flux_val
            if len(list(out)) == 0:
                out[run] = list(rf_merged.values())
                out.index = list(rf_merged.keys())
            else:
                runDF = pandas.DataFrame(list(rf_merged.values()),
                                         index=list(rf_merged.keys()), columns=[run])
                runDF = runDF.reindex(list(set(list(out.index)).union(
                    set(list(rf_merged.keys())))), fill_value=0)
                out = out.reindex(list(set(list(out.index)).union(
                    set(list(rf_merged.keys())))), fill_value=0)
                out = out.join(runDF, how='outer')
        return(out)


def buildExchangeMap(Controller):
    """
    Returns a map of all metabolites, the corresponding transport-reactions and stoichiometires;
    exchanged with the medium.
    {Metabolite1 : {ExchangeReaction1 : stoch-coefficient1 , ExchangeReaction2 : stoch-coefficient2},
    {Metabolite2 : {ExchangeReaction1 : stoch-coefficient1 , ExchangeReaction2 : stoch-coefficient2}}

    Metabolite1 - ... MetaboliteN : All metabolite-species in the medium (see medium.tsv file)
    ExchangeReaction1 - ... ExchangeReactionN : All metabolic reactions, which exchange the respective metabolite with the medium.
    stoch-coefficient : Stochiometric coefficient with which the respective metabolite is exchanged by the corresponding reaction.
    (Negative when reaction transports metabolite out of the cell; and positive when inside the cell.)

    Parameters
    ----------
    Controller : rbatools.NewControler.RBA_newControler

    Returns
    -------
    Dict.
    """
    BoundaryMetabolites = [i for i in list(Controller.ModelStructure.MetaboliteInfo.Elements.keys(
    )) if Controller.ModelStructure.MetaboliteInfo.Elements[i]['boundary']]
    ExchangeMap = {}
    for bM in BoundaryMetabolites:
        for rxn in Controller.ModelStructure.MetaboliteInfo.Elements[bM]['ReactionsInvolvedWith']:
            Reactants = list(
                Controller.ModelStructure.ReactionInfo.Elements[rxn]['Reactants'].keys())
            Products = list(Controller.ModelStructure.ReactionInfo.Elements[rxn]['Products'].keys())
            if len(list(set(list(Reactants+Products)))) > 1:
                for met in list(set(list(Reactants+Products))):
                    # if met != bM:
                    if met == bM:
                        MediumSpecies = findExchangeMetInMedium(met, Controller.Medium)
                        if met in Reactants:
                            stochCoeff = - \
                                Controller.ModelStructure.ReactionInfo.Elements[rxn]['Reactants'][met]
                        elif met in Products:
                            stochCoeff = Controller.ModelStructure.ReactionInfo.Elements[rxn]['Products'][met]
                        if MediumSpecies in list(ExchangeMap.keys()):
                            ExchangeMap[MediumSpecies].update({rxn: stochCoeff})
                        else:
                            ExchangeMap[MediumSpecies] = {rxn: stochCoeff}
    return(ExchangeMap)


def buildExchangeMap2(Controller):
    """
    Returns a map of all metabolites, the corresponding transport-reactions and stoichiometires;
    exchanged with the medium.
    {Metabolite1 : {ExchangeReaction1 : stoch-coefficient1 , ExchangeReaction2 : stoch-coefficient2},
    {Metabolite2 : {ExchangeReaction1 : stoch-coefficient1 , ExchangeReaction2 : stoch-coefficient2}}

    Metabolite1 - ... MetaboliteN : All metabolite-species in the medium (see medium.tsv file)
    ExchangeReaction1 - ... ExchangeReactionN : All metabolic reactions, which exchange the respective metabolite with the medium.
    stoch-coefficient : Stochiometric coefficient with which the respective metabolite is exchanged by the corresponding reaction.
    (Negative when reaction transports metabolite out of the cell; and positive when inside the cell.)

    Parameters
    ----------
    Controller : rbatools.NewControler.RBA_newControler

    Returns
    -------
    Dict.
    """
    BoundaryMetabolites = [i for i in list(Controller.ModelStructure.MetaboliteInfo.Elements.keys(
    )) if Controller.ModelStructure.MetaboliteInfo.Elements[i]['boundary']]
    ExchangeMap = {}
    for bM in BoundaryMetabolites:
        for rxn in Controller.ModelStructure.MetaboliteInfo.Elements[bM]['ReactionsInvolvedWith']:
            Reactants = list(
                Controller.ModelStructure.ReactionInfo.Elements[rxn]['Reactants'].keys())
            Products = list(Controller.ModelStructure.ReactionInfo.Elements[rxn]['Products'].keys())
            if len(list(set(list(Reactants+Products)))) == 1:
                for met in list(set(list(Reactants+Products))):
                    # if met != bM:
                    if met == bM:
                        MediumSpecies = findExchangeMetInMedium(met, Controller.Medium)
                        if met in Reactants:
                            stochCoeff = - \
                                Controller.ModelStructure.ReactionInfo.Elements[rxn]['Reactants'][met]
                        elif met in Products:
                            stochCoeff = Controller.ModelStructure.ReactionInfo.Elements[rxn]['Products'][met]
                        if MediumSpecies in list(ExchangeMap.keys()):
                            ExchangeMap[MediumSpecies].update({rxn: stochCoeff})
                        else:
                            ExchangeMap[MediumSpecies] = {rxn: stochCoeff}
    return(ExchangeMap)


def findExchangeMetInMedium(metabolite, Medium):
    """
    Returns the most likely species in the Medium, for any Metabolic species.
    Parameters
    ----------
    metabolite : str
    Medium : dict
    -------
    Most likely ID as str
    """
    if metabolite.endswith('_e'):
        out = difflib.get_close_matches('_e'.join(metabolite.split('_e')[: -1]), Medium, 1)
    else:
        out = difflib.get_close_matches(metabolite, Medium, 1)
    if len(out) > 0:
        return(out[0])
    else:
        return('')


def buildCompositionofUnusedProtein(species, composition):
    products = {}
    reactants = {}
    for j in composition.keys():
        stoch = composition[j]
        if len(list(species[j]['Products'].keys())) > 0:
            for k in species[j]['Products'].keys():
                if k in products:
                    products[k] += species[j]['Products'][k]*stoch
                else:
                    products.update({k: species[j]['Products'][k]*stoch})
        if len(list(species[j]['Reactants'].keys())) > 0:
            for k in species[j]['Reactants'].keys():
                if k in reactants:
                    reactants[k] += species[j]['Reactants'][k]*stoch
                else:
                    reactants.update({k: species[j]['Reactants'][k]*stoch})
    uniqueMets = list(set(list(products.keys())+list(reactants.keys())))
    NetMets = {}
    for j in uniqueMets:
        if j in list(products.keys()):
            produced = products[j]
        elif j not in list(products.keys()):
            produced = 0
        if j in list(reactants.keys()):
            consumed = reactants[j]
        elif j not in list(reactants.keys()):
            consumed = 0
        NetStoch = produced-consumed
        NetMets.update({j: float(NetStoch)})
    return(NetMets)


def buildUsedProteinConstraint(Controler, protein):
    ProteinConstraint = numpy.zeros(Controler.Problem.LP.A.shape[1])
    ConsumersEnzymes = Controler.ModelStructure.ProteinInfo.Elements[protein]['associatedEnzymes']
    ConsumersProcess = Controler.ModelStructure.ProteinInfo.Elements[protein]['SupportsProcess']
    for j in ConsumersEnzymes:
        StochFactor = Controler.ModelStructure.EnzymeInfo.Elements[j]['Subunits'][protein]['StochFac']
        LikeliestVarName = difflib.get_close_matches(j, Controler.Problem.LP.col_names, 1)[0]
        ColIndex = Controler.Problem.LP.colIndicesMap[LikeliestVarName]
        ProteinConstraint[ColIndex] = float(StochFactor)
    for j in ConsumersProcess:
        LikeliestVarName = difflib.get_close_matches(str(
            Controler.ModelStructure.ProcessInfo.Elements[j]['ID']+'_machinery'), Controler.Problem.LP.col_names, 1)[0]
        StochFactor = Controler.ModelStructure.ProcessInfo.Elements[j]['Composition'][protein]
        ColIndex = Controler.Problem.LP.colIndicesMap[LikeliestVarName]
        ProteinConstraint[ColIndex] = float(StochFactor)
    return(ProteinConstraint)


def determineCoefficient(x, changes, species):
    multiplicativeFactors = []
    for k in x:
        result = 1
        type = list(k.values())[0]['FunctionType']
        pars = list(k.values())[0]['FunctionParameters']
        if type == 'constant':
            result = numpy.float64(pars['C'])
        if type == 'exponential':
            L = 1
            if 'Lambda' in list(pars.keys()):
                L = numpy.float64(pars['Lambda'])
            result = numpy.exp(float(changes[species])*L)
        if type == 'indicator':
            maxi = numpy.inf
            mini = -numpy.inf
            if 'xMax' in list(pars.keys()):
                maxi = numpy.float64(pars['xMax'])
            if 'xMin' in list(pars.keys()):
                mini = numpy.float64(pars['xMin'])
            result = (float(changes[species]) > mini) and (float(changes[species]) < maxi)
        if type == 'linear':
            X_maxi = numpy.inf
            X_mini = -numpy.inf
            Y_maxi = numpy.inf
            Y_mini = -numpy.inf
            A = 1
            C = 0
            if 'A' in list(pars.keys()):
                A = numpy.float64(pars['A'])
            if 'C' in list(pars.keys()):
                C = numpy.float64(pars['C'])
            if 'xMin' in list(pars.keys()):
                X_mini = numpy.float64(pars['xMin'])
            if 'xMax' in list(pars.keys()):
                X_maxi = numpy.float64(pars['xMax'])
            if 'yMin' in list(pars.keys()):
                Y_mini = numpy.float64(pars['yMin'])
            if 'yMax' in list(pars.keys()):
                Y_maxi = numpy.float64(pars['yMax'])
            X = float(changes[species])
            if float(changes[species]) < X_mini:
                X = X_mini
            if float(changes[species]) > X_maxi:
                X = X_maxi
            Y = A*X + C
            result = Y
            if Y < Y_mini:
                result = Y_mini
            if Y > Y_maxi:
                result = Y_maxi
        if type == 'michaelisMenten':
            Y_mini = -numpy.inf
            KM = 0
            VM = 1
            if 'Km' in list(pars.keys()):
                KM = numpy.float64(pars['Km'])
            if 'Vmax' in list(pars.keys()):
                VM = numpy.float64(pars['Vmax'])
            if 'yMin' in list(pars.keys()):
                Y_mini = numpy.float64(pars['yMin'])
            Y = VM*float(changes[species])/(float(changes[species])+KM)
            result = Y
            if Y < Y_mini:
                result = Y_mini
        if type == 'competitiveInhibition':
            Y_mini = -numpy.inf
            KM = 0
            VM = 1
            KI = 0
            I = 0
            if 'Ki' in list(pars.keys()):
                KI = numpy.float64(pars['Ki'])
            if 'I' in list(pars.keys()):
                I = numpy.float64(pars['I'])
            if 'Km' in list(pars.keys()):
                KM = numpy.float64(pars['Km'])
            if 'Vmax' in list(pars.keys()):
                VM = numpy.float64(pars['Vmax'])
            if 'yMin' in list(pars.keys()):
                Y_mini = numpy.float64(pars['yMin'])
            Y = VM*float(changes[species])/(float(changes[species])+KM*(1+I/KI))
            result = Y
            if Y < Y_mini:
                result = Y_mini
        if type == 'inverse':
            C = 1
            if 'C' in list(pars.keys()):
                C = numpy.float64(pars['C'])
            result = 1
            if float(changes[species]) is not 0:
                result = C/float(changes[i])
        multiplicativeFactors.append(result)
        value = numpy.prod(numpy.array(multiplicativeFactors))

    return(float(value))