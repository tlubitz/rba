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
import rba
from .rba_SimulationData import RBA_SimulationData
from .rba_SimulationParameters import RBA_SimulationParameters
from .rba_ModelStructure import RBA_ModelStructure
from .rba_Problem import RBA_Problem
from .rba_Matrix import RBA_Matrix
from .rba_LogBook import RBA_LogBook
from ..rba import model
from ..rba.core import ConstraintMatrix
from ..rba.core import Solver


class RBA_newControler(object):
    """
    Top level of the RBA API.

    Attributes
    ----------
    xml_dir : str
        Current Growth rate as numeric value
    model : rba.RbaModel
        Current Growth rate as numeric value
    matrices : ConstraintMatrix
        Current Growth rate as numeric value
    solver : rba.Solver
        Current Growth rate as numeric value
    Problem : rbastructure.RBA_Problem
        Current Growth rate as numeric value
    Medium : dict
        Current Growth rate as numeric value
    ModelStructure : rbastructure.RBA_ModelStructure
        Current Growth rate as numeric value
    Results : dict
        Current Growth rate as numeric value
    Parameters : dict
        Current Growth rate as numeric value
    SimulationData : rbastructure.RBA_SimulationData
        Current Growth rate as numeric value
    SimulationParameters : rbastructure.RBA_SimulationParameters
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
    """

    def __init__(self, xml_dir):
        """
        Creates controler object from files

        Parameters
        ----------
        xml_dir : str
            Pstartedath to the directory where rba-model files are located.
        """
        self.xml_dir = xml_dir
        self.LogBook = RBA_LogBook('Controler')

        if not hasattr(self, 'ModelStructure'):
            self.ModelStructure = RBA_ModelStructure()
            if os.path.isfile(str(self.xml_dir+'/ModelStructure.json')):
                with open(str(self.xml_dir+'/ModelStructure.json'), 'r') as myfile:
                    data = myfile.read()
                self.ModelStructure.fromJSON(inputString=data)
            else:
                self.ModelStructure.fromFiles(xml_dir=self.xml_dir)
                self.ModelStructure.exportJSON(path=self.xml_dir)

        self.model = model.RbaModel.from_xml(input_dir=xml_dir)
        self.matrices = ConstraintMatrix(model=self.model)
        self.solver = Solver(matrix=self.matrices)

        self.LogBook.addEntry('Model loaded from {}.'.format(self.xml_dir))
        self.Problem = RBA_Problem(solver=self.solver)

        medium = pandas.read_csv(xml_dir+'/medium.tsv', sep='\t')
        self.Medium = dict(zip(list(medium.ix[:, 0]), [float(i) for i in list(medium.ix[:, 1])]))

        self.Mu = self.Problem.Mu
        self.ExchangeMap = buildExchangeMap(self)

    def reloadModel(self):
        """
        Reloads model from files
        """
        self.LogBook.addEntry('Model reloaded from {}.'.format(self.xml_dir))
        self.model = rba.RbaModel.from_xml(input_dir=self.xml_dir)
        self.matrices = ConstraintMatrix(model=self.model)
        self.solver = Solver(matrix=self.matrices)
        self.Problem = RBA_Problem(solver=self.solver)

    def recordResults(self, runName):
        """
        Records Simulation output for further use.

        Can be written with 'writeResults' method.

        Parameters
        ----------
        runName : str
            Name of observation.
            Serves as ID for all Data, originating from these.
        """
        self.LogBook.addEntry('Solution recorded under {}.'.format(runName))
        if not hasattr(self, 'Results'):
            self.Results = {'Reactions': pandas.DataFrame(index=self.Problem.Reactions),
                            'Enzymes': pandas.DataFrame(index=self.Problem.Enzymes),
                            'Processes': pandas.DataFrame(index=self.Problem.Processes),
                            'Proteins': pandas.DataFrame(index=list(self.ModelStructure.ProteinMatrix['Proteins'])),
                            'ProtoProteins': pandas.DataFrame(index=list(self.ModelStructure.ProteinGeneMatrix['ProtoProteins'])),
                            'Constraints': pandas.DataFrame(index=self.Problem.LP.row_names),
                            'ObjectiveFunction': pandas.DataFrame(index=self.Problem.LP.col_names),
                            'Mu': pandas.DataFrame(index=['Mu']),
                            'ObjectiveValue': pandas.DataFrame(index=['ObjectiveValue'])}
        self.Results['Reactions'][runName] = [self.Problem.SolutionValues[i]
                                              for i in self.Problem.Reactions]
        self.Results['Enzymes'][runName] = [self.Problem.SolutionValues[i]
                                            for i in self.Problem.Enzymes]
        self.Results['Processes'][runName] = [self.Problem.SolutionValues[i]
                                              for i in self.Problem.Processes]
        self.Results['Constraints'][runName] = [self.Problem.DualValues[i]
                                                for i in self.Problem.LP.row_names]
        self.Results['Proteins'][runName] = ProteomeRecording(self, runName)
        self.Results['ProtoProteins'][runName] = ProtoProteomeRecording(
            self, runName, self.Results['Proteins'])
        self.Results['Mu'][runName] = self.Problem.Mu
        self.Results['ObjectiveValue'][runName] = self.Problem.ObjectiveValue
        self.Results['ObjectiveFunction'][runName] = list(self.Problem.getObjective().values())

    def returnExchangeFluxes(self):
        out = {}
        for j in self.ExchangeMap.keys():
            netflux = 0
            for k in self.ExchangeMap[j].keys():
                netflux += self.ExchangeMap[j][k]*self.Problem.SolutionValues[k]
            if netflux != 0:
                out[j] = netflux
        return(out)

    def recordParameters(self, runName):
        """
        Records Simulation parameters for further use.

        Can be written with 'writeResults' method.

        Parameters
        ----------
        runName : str
            Name of observation.
            Serves as ID for all Data, originating from these.
        """
        self.LogBook.addEntry('Coefficients recorded under {}.'.format(runName))
        if not hasattr(self, 'Parameters'):
            EnzymeCapacities = self.Problem.getEnzymeCapacities()
            ProcessCapacities = self.Problem.getProcessCapacities()
            CompartmentCapacities = self.Problem.getCompartmentCapacities()

            self.Parameters = {'EnzymeEfficiencies': pandas.DataFrame(index=list(EnzymeCapacities.keys())),
                               'NetProcessEfficiencies': pandas.DataFrame(index=list(ProcessCapacities.keys())),
                               'CompartmentCapacities': pandas.DataFrame(index=list(CompartmentCapacities.keys()))}

        self.Parameters['EnzymeEfficiencies'][runName+'_FW'] = [EnzymeCapacities[i]
                                                                ['Forward'] for i in list(EnzymeCapacities.keys())]
        self.Parameters['EnzymeEfficiencies'][runName+'_BW'] = [EnzymeCapacities[i]
                                                                ['Backward'] for i in list(EnzymeCapacities.keys())]
        self.Parameters['NetProcessEfficiencies'][runName] = [ProcessCapacities[i]
                                                              for i in list(ProcessCapacities.keys())]
        self.Parameters['CompartmentCapacities'][runName] = [CompartmentCapacities[i]
                                                             for i in list(CompartmentCapacities.keys())]

    def clearResults(self):
        """
        Removes all previosly recorded results.
        """
        self.LogBook.addEntry('Results cleared.')

        delattr(self, 'Results')

    def clearParameters(self):
        """
        Removes all previosly recorded parameters.
        """

        self.LogBook.addEntry('Parameters cleared.')
        delattr(self, 'Parameters')

    def writeResults(self, session_name='', digits=25, loggingIntermediateSteps=False):
        """
        Creates SimulationData and SimulationParameters objects from recordings.

        Stores them as rbastructure.RBA_SimulationData
        and rbastructure.RBA_SimulationParameters objects as attributes.
        Access via "NameOfInstance".SimulationData and
        "NameOfInstance".RBA_SimulationParameters respectively.

        Parameters
        ----------
        digits : int
            Number of decimal places in the results
            Default: 10
        session_name : str
            Name of Simulation session.
            Default: ''
        """
        self.LogBook.addEntry('Data written under {}.'.format(session_name))
        if hasattr(self, 'Results'):
            #            self.Results['Proteins'] = predictProteome(Controller=self)
            self.Results['uniqueReactions'] = mapIsoReactions(Controller=self)
            self.Results['Mu'] = self.Results['Mu'].round(digits)
            self.Results['ObjectiveValue'] = self.Results['ObjectiveValue'].round(digits)
            self.Results['Proteins'] = self.Results['Proteins'].round(digits)
            self.Results['uniqueReactions'] = self.Results['uniqueReactions'].round(digits)
            self.Results['Reactions'] = self.Results['Reactions'].round(digits)
            self.Results['Enzymes'] = self.Results['Enzymes'].round(digits)
            self.Results['Processes'] = self.Results['Processes'].round(digits)
            self.Results['Constraints'] = self.Results['Constraints'].round(digits)

            self.SimulationData = RBA_SimulationData(StaticData=self.ModelStructure)
            self.SimulationData.fromSimulationResults(Controller=self, session_name=session_name)

        if hasattr(self, 'Parameters'):
            self.Parameters['EnzymeEfficiencies'] = self.Parameters['EnzymeEfficiencies'].round(
                digits)
            self.Parameters['NetProcessEfficiencies'] = self.Parameters['NetProcessEfficiencies'].round(
                digits)
            self.Parameters['CompartmentCapacities'] = self.Parameters['CompartmentCapacities'].round(
                digits)
            if deleteZeros:
                self.Parameters['EnzymeEfficiencies'] = self.Parameters['EnzymeEfficiencies'][(
                    self.Parameters['EnzymeEfficiencies'].T != 0).any()]
                self.Parameters['NetProcessEfficiencies'] = self.Parameters['NetProcessEfficiencies'][(
                    self.Parameters['NetProcessEfficiencies'].T != 0).any()]
                self.Parameters['CompartmentCapacities'] = self.Parameters['CompartmentCapacities'][(
                    self.Parameters['CompartmentCapacities'].T != 0).any()]
            self.SimulationParameters = RBA_SimulationParameters(StaticData=self.ModelStructure)
            self.SimulationParameters.fromSimulationResults(Controller=self)

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

        Does the same as rbastructure.RBA_Problem.solveLP().
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

    def findMaxGrowthRate(self, precision=0.00001, max=4, recording=False, loggingIntermediateSteps=False):
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
        testMu = maxMu
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
        return(minMu)

    def findMinMediumConcentration(self, metabolite, precision=0.00001, max=100, recording=False, loggingIntermediateSteps=False):
        """
        Applies dichotomy-search to find the minimal feasible concentration of
        growth-substrate in medium.

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
        Sets the concentration of growth-substrate in medium.

        Parameters
        ----------
        changes : dict
            Keys : ID of metabolite in medium.
            Vaules : New concention
        """

        for species in (changes.keys()):
            self.LogBook.addEntry(
                'Concentration of {} changed: {} --> {}.'.format(species, self.Medium[species], float(changes[species])))
            ChangeCharacteristics = QualitativeMediumChange(
                Controller=self, changes=changes, species=species)
            if ChangeCharacteristics[0]:
                Exchanges = findExchangeReactions(Controller=self, species=species)
                for exRx in list(Exchanges.keys()):
                    if Exchanges[exRx] is 'Product':
                        self.Problem.setUB(
                            inputDict={exRx: ChangeCharacteristics[1]}, logging=loggingIntermediateSteps)
                    elif Exchanges[exRx] is 'Reactant':
                        if self.ModelStructure.ReactionInfo.Elements[exRx]['Reversible']:
                            self.Problem.setLB(
                                inputDict={exRx: -ChangeCharacteristics[1]}, logging=loggingIntermediateSteps)

            # List of coefficients (named after Matrix-constraints) depend on species#
            MediumDependentCoefficients = []
            if species in list(self.ModelStructure.MediumDependencies.keys()):
                MediumDependentCoefficients = self.ModelStructure.MediumDependencies[species]

            for coeff in MediumDependentCoefficients:
                if coeff in list(self.ModelStructure.EnzymeConstraintsInfo.Elements.keys()):
                    parameters = self.ModelStructure.EnzymeConstraintsInfo.Elements[coeff]['CapacityParameter']
                    rowName = coeff
                    colName = self.ModelStructure.EnzymeConstraintsInfo.Elements[coeff]['AssociatedEnzyme']
                elif coeff in list(self.ModelStructure.ProcessConstraintsInfo.Elements.keys()):
                    parameters = self.ModelStructure.ProcessConstraintsInfo.Elements[coeff]['CapacityParameter']
                    rowName = coeff
                    colName = self.ModelStructure.ProcessConstraintsInfo.Elements[coeff]['AssociatedProcess']
                elif coeff in list(self.ModelStructure.DensityConstraintsInfo.Elements.keys()):
                    parameters = self.ModelStructure.DensityConstraintsInfo.Elements[coeff]['CapacityParameter']
                    rowName = coeff

                # Evaluate Parameter Formula#
                value = determineCoefficient(x=parameters, changes=changes, species=species)

                if coeff in list(self.ModelStructure.DensityConstraintsInfo.Elements.keys()):
                    self.Problem.setRighthandSideValue(
                        inputDict={rowName: value}, logging=loggingIntermediateSteps)
                elif coeff in list(self.ModelStructure.EnzymeConstraintsInfo.Elements.keys()):
                    self.Problem.setProblemCoefficients(
                        inputDict={rowName: {colName: -value}}, logging=loggingIntermediateSteps)

            self.Medium[species] = float(changes[species])

    def knockOut(self, gene, loggingIntermediateSteps=False):
        """
        Emulates a gene knock out.

        Parameters
        ----------
        gene : str
            ID of model-protein to be knocked out.
        """

        if type(gene) is str:
            genes = [gene]
        if type(gene) is list:
            genes = gene
        for g in genes:
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

    def FeasibleRange(self, loggingIntermediateSteps=False, *variables):
        """
        ***
        ATTENTION! XXX! Changed order of parameters given to function because
        as of new *args must be declared LAST when also having a parameter without
        specification!
        ***

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

    def eukaryoticDensities(self, totalAA=3.1, CompartmentRelationships=True, CompartmentComponents=False):
        Compartments = ['n', 'mIM', 'vM', 'mIMS', 'm', 'erM', 'mOM', 'x', 'c', 'cM', 'gM']
        totalAA = 3.1*0.71
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
        #CompartmentMatrix.row_signs = ['L']*(len(Compartments)+1)
        CompartmentMatrix.row_signs = ['E']*(len(Compartments)+1)
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
        #CompartmentMatrix.row_signs = ['E']*(len(Compartments)+1)
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
                #Alipids.LB = numpy.array(list(Alipids.LB).append(xnew))
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
    associatedReactions = {k: Controller.ModelStructure.MetaboliteInfo.Elements[k]['ReactionsInvolvedWith'] for k in list(
        Controller.ModelStructure.MetaboliteInfo.Elements.keys()) if species in k}
    exchanges = {}
    for met in list(associatedReactions.keys()):
        for rxn in associatedReactions[met]:
            Reaction = Controller.ModelStructure.ReactionInfo.Elements[rxn]
            if 'Transport' in Reaction['Type']:
                isProduct = False
                isReactant = False
                for r in list(Reaction['Reactants'].keys()):
                    if str(species + '_') in r:
                        isReactant = True
                for p in list(Reaction['Products'].keys()):
                    if str(species + '_') in p:
                        isProduct = True
                if isProduct and not isReactant:
                    if len(Reaction['Formula'].split(' <=> ')[0])*len(Reaction['Formula'].split(' <=> ')[1]) == 0:
                        exchanges.update({rxn: 'Product'})
                elif isReactant and not isProduct:
                    if len(Reaction['Formula'].split(' <=> ')[0])*len(Reaction['Formula'].split(' <=> ')[1]) == 0:
                        exchanges.update({rxn: 'Reactant'})
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


def predictProteome(Controller):
    ProteinProteinMatrix = Controller.ModelStructure.ProteinMatrix['Matrix'].astype(numpy.float64)
    C = Controller.ModelStructure.ProteinMatrix['Consumers']
    Consumers = []
    for i in C:
        if i.startswith('P_'):
            Consumers.append(str(i+'_machinery'))
        if not i.startswith('P_'):
            Consumers.append(i)
    Proteins = Controller.ModelStructure.ProteinMatrix['Proteins']
    DF = pandas.concat([Controller.Results['Enzymes'], Controller.Results['Processes']], axis=0)
    ProteinLevels = pandas.DataFrame(index=Proteins)
    for run in list(DF.columns):
        vector = numpy.nan_to_num(DF[run].reindex(Consumers))
        Level = ProteinProteinMatrix.dot(vector)
        ProteinLevels[run] = Level
    return(ProteinLevels)


def ProtoProteomeRecording(Controller, run, Proteinlevels):
    out = []
    for i in list(Controller.ModelStructure.ProteinGeneMatrix['ProtoProteins']):
        row_ind = list(Controller.ModelStructure.ProteinGeneMatrix['ProtoProteins']).index(i)
        # print(row_ind)
        nonZero = list(numpy.where(
            Controller.ModelStructure.ProteinGeneMatrix['Matrix'][row_ind, ] != 0)[0])
        level = 0
        for j in nonZero:
            id = Controller.ModelStructure.ProteinGeneMatrix['Proteins'][j]
            level += Proteinlevels.loc[id, run]
        out.append(level)
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
            ProteinLevels[run].ix[protID] = Controller.Problem.SolutionValues[p]
    return(list(ProteinLevels[run]))


def mapIsoReactions(Controller):
    if hasattr(Controller, 'Results'):
        out = pandas.DataFrame()
        for run in list(Controller.Results['Reactions'].columns):
            rf = dict(zip(list(Controller.Results['Reactions'].index),
                          list(Controller.Results['Reactions'][run])))
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
    Exchanges = [i for i in Controller.ModelStructure.ReactionInfo.Elements.keys() if Controller.ModelStructure.ReactionInfo.Elements[i]
                 ['Type'] == 'Transport (-->)' and len(list(set(Controller.ModelStructure.ReactionInfo.Elements[i]['Compartment_Species']))) == 1]
    ExchangeMap = {}
    for i in Exchanges:
        for j in Controller.ModelStructure.ReactionInfo.Elements[i]['Reactants'].keys():
            d = -Controller.ModelStructure.ReactionInfo.Elements[i]['Reactants'][j]
            if findExchangeMetInMedium(j, Controller.Medium) not in ExchangeMap.keys():
                ExchangeMap.update({findExchangeMetInMedium(j, Controller.Medium): {i: d}})
            else:
                ExchangeMap[findExchangeMetInMedium(j, Controller.Medium)].update({i: d})
        for j in Controller.ModelStructure.ReactionInfo.Elements[i]['Products'].keys():
            d = Controller.ModelStructure.ReactionInfo.Elements[i]['Products'][j]
            if findExchangeMetInMedium(j, Controller.Medium) not in ExchangeMap.keys():
                ExchangeMap.update({findExchangeMetInMedium(j, Controller.Medium): {i: d}})
            else:
                ExchangeMap[findExchangeMetInMedium(j, Controller.Medium)].update({i: d})
    return(ExchangeMap)


def findExchangeMetInMedium(metabolite, Medium):
    return(metabolite)
    # for i in Medium.keys():
    #    if i in metabolite:
    #        return(i)