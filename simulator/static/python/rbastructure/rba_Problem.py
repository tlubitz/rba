from __future__ import division, print_function
import numpy
import scipy
import copy
import cplex
from .rba_Matrix import RBA_Matrix
from .rba_LP import RBA_LP
from .rba_LogBook import RBA_LogBook


class RBA_Problem(object):
    """
    Class holding RBA-problem as mathematical manifestation of the RBA-model.

    Attributes
    ----------
    ClassicRBAmatrix : rba.solver.matrix object
    Mu : float
        Current Growth rate as numeric value
    LP : rbastructure.LP object
    Enzyme_FWcapacities : list
        List of constraint IDs, which represent forward efficiencies of enzymes
        Created by method 'extractConstraintTypes'
    Enzyme_BWcapacities : list
        List of constraint IDs, which represent backward efficiencies of enzymes
        Created by method 'extractConstraintTypes'
    ProcessCapacities : list
        List of constraint IDs, which represent net efficiencies of processes
        Created by method 'extractConstraintTypes'
    Metabolites : list
        List of constraint IDs, which represent mass balances of metabolites
        Created by method 'extractConstraintTypes'
    CompartmentDensities : list
        List of constraint IDs, which represent compartment capacities
        Created by method 'extractConstraintTypes'
    Reactions : list
        List of constraint IDs, which represent metabolic reactions
        Created by method 'extractVariableTypes'
    Enzymes : list
        List of constraint IDs, which represent enzymes, associated with metabolic reactions
        Created by method 'extractVariableTypes'
    Processes : list
        List of constraint IDs, which represent process machineries
        Created by method 'extractVariableTypes'
    MuDepIndices_A : list
        List of tuples holding rows (constraint IDs) and columns (variable IDs) of
        constraint-matrix coefficients (LHS), which depend on the growth rate.
        Created by method 'findMudependencies'
    MuDepIndices_b: list
        List of constraint IDs, whos RHS depend on the growth rate
        Created by method 'findMudependencies'
    MuDepIndices_LB: list
        List of variable IDs, whos lower-bounds depend on the growth rate
        Created by method 'findMudependencies'
    MuDepIndices_UB: list
        List of variable IDs, whos upper-bounds depend on the growth rate
        Created by method 'findMudependencies'
    Solved: bool
        Booelean indicating wheter Problem has been successfully solved,
        according to the acceptable solution statuses provided to the solveLP-method.
        Created by method 'solveLP'
    SolutionStatus: int
        Numerical value indicating the solution-status of the problem.
        Consult CPLEX documentation for their meaning.
        Created by method 'solveLP'
    ObjectiveValue: float
        Numeric value of the objective function after optimisation by CPLEX.
        Created by method 'solveLP'
    SolutionValues: dict
        Solution vector after optimisation by CPLEX.
        (Dictionary with variable IDs as keys and numeric values as values)
        Created by method 'solveLP'
    DualValues: dict
        Vector of dual-values after optimisation by CPLEX.
        (Dictionary with constraint IDs as keys and numeric values as values)
        Created by method 'solveLP'

    Methods
    ----------
    __init__(solver)
        Initiates rba_Problem from supplied rba.solver object.

    extractConstraintTypes()
        Extracts information on the different constraint types in the standard RBA-matrix.

    extractVariableTypes()
        Extracts information on the different variable types in the standard RBA-matrix.

    findMudependencies()
        Extracts information on the growth-rate dependent LP-coefficients.

    BuildClassicMatrix(Mu)
        Builds standard RBA-matrix according to growth-rate

    solveLP(feasibleStatuses=[1, 2, 5, 6])
        Solves Linear RBA problem.

    setMu(Mu, keepParameters=None)
        Changes growth-rate of problem and sets all associated coefficients.

    resetCPLEXparams()
        Sets cplex.Cplex() parameters to predefined values.

    setConstraintType(inputDict)
        Sets objective function coefficients.

    getConstraintType(*constraints)
        Returns type of constraints.

    getObjective(*variables)
        Returns objective coefficient of variables in problem.

    setObjectiveCoefficients(inputDict)
        Sets objective function coefficients.

    clearObjective()
        Sets all coefficients of the objective function to zero.

    invertObjective()
        Changes sign (optimisation-sense) of objective function.

    getRighthandSideValue(*constraints)
        Extracts coefficients of problem's righthand side (B-vector).

    setRighthandSideValue(inputDict)
        Set coefficients of the problems' RHS (b-vector).

    calculateLefthandSideValue(*constraints)
        Calculates value of problem's lefthand side
        (after multiplying with solution-vector).

    getProblemCoefficients(*inputTuples)
        Returns coefficients of LHS of problem.

    setProblemCoefficients(inputDict)
        Set coefficients of the problems' LHS (constraint matrix).

    getUB(*variables)
        Returns upper bounds of problem variables.

    setUB(inputDict)
        Set upper-bounds of the problem variables.

    getLB(*variables)
        Returns lower bounds of problem variables.

    setLB(inputDict)
        Set lower-bounds of the problem variables.

    getEnzymeCapacities(*Enzyme)
        Returns of capacity coefficients oof enzymes in problem.

    getProcessCapacities(*Process)
        Returns of capacity coefficients of process machineries in problem.

    getCompartmentCapacities(*Compartment)
        Returns of capacity coefficients of compartments in problem.

    """

    def __init__(self, solver):
        """
        Initiates rba_Problem from supplied rba.solver object.

        Makes sure the problem is consistent and optimisable.

        Parameters
        ----------
        solver : rba.solver
        """

        ## Import solver-object information ##
        self.ClassicRBAmatrix = solver.matrix
        self.Mu = 0
        ## Set Mu of solver to 0 ##
        self.LogBook = RBA_LogBook('Problem')
        self.BuildClassicMatrix(Mu=self.Mu)
        ## Initiate LP-object ##
        self.LP = RBA_LP()
        ## Transfer solver information to LP-object ##
        self.LP.loadMatrix(matrix=self.ClassicRBAmatrix)
        self.LP.buildCPLEX_LP()
        ## Solve for Mu=0 to check if problem is consistent##
        self.LogBook.addEntry('Problem initiated.')
        self.solveLP()
        if not self.Solved:
            ## Problem is wrongly defined and can not be solved ##
            print(' INCONSISTENT RBA-PROBLEM; NOT SOLVABLE --> ABORTED!\n')
            self.LogBook.addEntry('Problem inconsistent at growth-rate 0.')
            return
        ## If problem is consistent the growth-rate dependent indices ##
        ## and constraint/variable types are extractet ##
        self.findMudependencies()
        self.extractConstraintTypes()
        self.extractVariableTypes()

    def extractConstraintTypes(self):
        """
        Extracts information on the different constraint types in the standard RBA-matrix.

        AddsAttributes: Enzyme_FWcapacities, Enzyme_BWcapacities, ProcessCapacities,
                        Metabolites and CompartmentDensities
        """

        self.Enzyme_FWcapacities = [i for i in self.LP.row_names if i.endswith(
            '_forward_capacity') and i.startswith('R_')]
        self.Enzyme_BWcapacities = [i for i in self.LP.row_names if i.endswith(
            '_backward_capacity') and i.startswith('R_')]
        self.ProcessCapacities = [i for i in self.LP.row_names if i.endswith(
            '_capacity') and i.startswith('P_')]
        self.Metabolites = [i for i in self.LP.row_names if i.startswith('M_')]
        self.CompartmentDensities = [i for i in self.LP.row_names if i.endswith('_density')]

    def extractVariableTypes(self):
        """
        Extracts information on the different variable types in the standard RBA-matrix.

        AddsAttributes: Enzymes, Reactions and Processes
        """

        self.Reactions = [i for i in self.LP.col_names if i.startswith(
            'R_') and not i.endswith('_enzyme')]
        self.Enzymes = [i for i in self.LP.col_names if i.startswith(
            'R_') and i.endswith('_enzyme')]
        self.Processes = [i for i in self.LP.col_names if i.startswith(
            'P_') and i.endswith('_machinery')]

    def findMudependencies(self):
        """
        Extracts information on the growth-rate dependent LP-coefficients.

        AddsAttributes: MuDepIndices_A, MuDepIndices_b, MuDepIndices_LB
                        and MuDepIndices_UB
        """

        ## Construct 2 rba.solver.matrix objects for growth-rate 0 and 1.1 ##
        MuOneMatrix = copy.deepcopy(self.ClassicRBAmatrix)
        MuOneMatrix.build_matrices(1.0)
        self.MuOneMatrix = RBA_LP()
        self.MuOneMatrix.loadMatrix(MuOneMatrix)

        M0 = copy.deepcopy(self.ClassicRBAmatrix)
        M11 = copy.deepcopy(self.ClassicRBAmatrix)
        M0.build_matrices(0)
        M11.build_matrices(1.1)
        ## Build arrays from matrices ##
        A_0 = M0.A.toarray()
        A_11 = M11.A.toarray()
        ## Find index pairs at which the two constraint-matrices differ ##
        MuDeps = numpy.where(A_11 != A_0)
        ## Transform list of numeric indices in list of tuples with row- and col names ##
        self.MuDepIndices_A = list(zip([self.ClassicRBAmatrix.row_names[i] for i in MuDeps[0]], [
                                   self.ClassicRBAmatrix.col_names[j] for j in MuDeps[1]]))
#        self.MuDepIndices_Adict = {i: (self.LP.rowIndicesMap[i[0]], self.LP.colIndicesMap[i[1]]) for i in self.MuDepIndices_A}

        ## Find rows at which the two righthandsides differ ##
        self.MuDepIndices_b = [n for n in self.ClassicRBAmatrix.row_names if M0.b[self.ClassicRBAmatrix.row_names.index(
            n)] != M11.b[self.ClassicRBAmatrix.row_names.index(n)]]
        ## Find columns at which the variable bounds differ ##
        self.MuDepIndices_LB = [n for n in self.ClassicRBAmatrix.col_names if M0.LB[self.ClassicRBAmatrix.col_names.index(
            n)] != M11.LB[self.ClassicRBAmatrix.col_names.index(n)]]
        self.MuDepIndices_UB = [n for n in self.ClassicRBAmatrix.col_names if M0.UB[self.ClassicRBAmatrix.col_names.index(
            n)] != M11.UB[self.ClassicRBAmatrix.col_names.index(n)]]

    def BuildClassicMatrix(self, Mu):
        """
        Builds standard RBA-matrix according to growth-rate

        Parameters
        ----------
        Mu : float
            Growth rate
        """

        self.Mu = Mu
        self.ClassicRBAmatrix.build_matrices(Mu)

    def solveLP(self, feasibleStatuses=[1, 2, 5, 6], logging=True):
        """
        Solves Linear RBA problem.

        Optimises RBA-LP with CPLEX.
        When cplex-solution status is amongst the user-defined feasible statuses;
        boolean 'Solved' is set to True and 'ObjectiveValue', 'SolutionValues' and 'DualValues'
        are stored as attributes.

        Parameters
        ----------
        feasibleStatuses : list of int
            List with identifiers of acceptable solution statuses.
            (consult ILOG-CPLEX documentation for information on them).
            Default: feasibleStatuses=[1, 2, 5, 6]
        logging : bool
            Wheter to write change to log-file or not
        """

        self.Solved = False
        ## Solve cplex LP ##
        self.LP.cplexLP.solve()
        ## Determin solution-status ##
        self.SolutionStatus = self.LP.cplexLP.solution.get_status()
        ## Check if solution status is amongst acceptable ones ##
        if logging:
            self.LogBook.addEntry(
                'Problem solved with solution status: {}'.format(self.SolutionStatus))
        if self.SolutionStatus in feasibleStatuses:
            ## Extract solution-data ##
            self.Solved = True
            self.ObjectiveValue = self.LP.cplexLP.solution.get_objective_value()
            self.SolutionValues = dict(
                zip(self.LP.col_names, self.LP.cplexLP.solution.get_values()))
            self.DualValues = dict(
                zip(self.LP.row_names, self.LP.cplexLP.solution.get_dual_values()))

    def setMu_classicProblem(self, Mu, keepParameters=None, logging=True):
        """
        Changes growth-rate of problem and sets all associated coefficients.

        Can be provided with 'keepParameters' argument
        to define problem coefficients which should remain unchanged.

        Parameters
        ----------
        Mu : float
            Growth rate
        keepParameters : dict
            Dictionary indicating which elements of the linear problem should
            not be affected when setting growth-rate. Possible keys of dictionary:
                'LHS': List of index-tuples (constraint-ID,variable-ID),
                       indicating elements of the lefthandside (constraint-matrix).
                'RHS': List of constraint IDs, indicating elements of the righthandside (b-vector).
                'LB': List of variable-IDs, indicating elements of the lower-bound vector.
                'UB': List of variable-IDs, indicating elements of the upper-bound vector.
            Default: keepParameters=None
        logging : bool
            Wheter to write change to log-file or not
                """
        ## Set Mu of solver-object ##
        self.BuildClassicMatrix(Mu=float(Mu))
        ## Construct new LP from solver object ##
        inputMatrix = RBA_Matrix()
        inputMatrix.loadMatrix(matrix=self.ClassicRBAmatrix)
        ## Define Indices which define elements of the new LP to be passed ##
        ##  on to the new LP; as all indices which depend on Mu ##

        A_idxs = self.MuDepIndices_A
        B_idxs = self.MuDepIndices_b
        LB_idxs = self.MuDepIndices_LB
        UB_idxs = self.MuDepIndices_UB
        if keepParameters is not None:
            ## If indices are passed, which define elements to be not changed when setting Mu ##
            ## these indices are removed from the update-from-new-to-old indices ##
            if 'LHS' in list(keepParameters.keys()):
                A_idxs = set(self.MuDepIndices_A)-set(keepParameters['LHS'])
            if 'RHS' in list(keepParameters.keys()):
                B_idxs = set(self.MuDepIndices_b)-set(keepParameters['RHS'])
            if 'LB' in list(keepParameters.keys()):
                LB_idxs = set(self.MuDepIndices_LB)-set(keepParameters['LB'])
            if 'UB' in list(keepParameters.keys()):
                UB_idxs = set(self.MuDepIndices_UB)-set(keepParameters['UB'])
        ## Pass new matrix and indices of elements to update to LP.updateMatrix-method ##
        self.LP.updateMatrix(matrix=inputMatrix, Ainds=A_idxs, Binds=B_idxs,
                             CTinds=[], LBinds=LB_idxs, UBinds=UB_idxs)

    def setMu_modifiedProblem(self, Mu, ModelStructure, keepParameters=None, logging=True):
        """
        Changes growth-rate of problem and sets all associated coefficients.

        Can be provided with 'keepParameters' argument
        to define problem coefficients which should remain unchanged.

        Parameters
        ----------
        Mu : float
            Growth rate
        ModelStructure : RBA_ModellStructure object.
        keepParameters : dict
            Dictionary indicating which elements of the linear problem should
            not be affected when setting growth-rate. Possible keys of dictionary:
                'LHS': List of index-tuples (constraint-ID,variable-ID),
                       indicating elements of the lefthandside (constraint-matrix).
                'RHS': List of constraint IDs, indicating elements of the righthandside (b-vector).
                'LB': List of variable-IDs, indicating elements of the lower-bound vector.
                'UB': List of variable-IDs, indicating elements of the upper-bound vector.
            Default: keepParameters=None
        logging : bool
            Wheter to write change to log-file or not
                """

        ## Set Mu of solver-object ##
        # self.BuildClassicMatrix(Mu=float(Mu))
        ## Construct new LP from solver object ##
        # inputMatrix = RBA_Matrix()
        # inputMatrix.loadMatrix(matrix=self.ClassicRBAmatrix)
        ## Define Indices which define elements of the new LP to be passed ##
        ##  on to the new LP; as all indices which depend on Mu ##

        inputMatrix = RBA_Matrix()
        inputMatrix.loadMatrix(matrix=self.MuOneMatrix)
        inputMatrix.scaleLHS(float(Mu))
        inputMatrix.tofloat64()

        A_idxs = self.MuDepIndices_A
        B_idxs = self.MuDepIndices_b
        LB_idxs = self.MuDepIndices_LB
        UB_idxs = self.MuDepIndices_UB
        if keepParameters is not None:
            ## If indices are passed, which define elements to be not changed when setting Mu ##
            ## these indices are removed from the update-from-new-to-old indices ##
            if 'LHS' in list(keepParameters.keys()):
                A_idxs = set(self.MuDepIndices_A)-set(keepParameters['LHS'])
            if 'RHS' in list(keepParameters.keys()):
                B_idxs = set(self.MuDepIndices_b)-set(keepParameters['RHS'])
            if 'LB' in list(keepParameters.keys()):
                LB_idxs = set(self.MuDepIndices_LB)-set(keepParameters['LB'])
            if 'UB' in list(keepParameters.keys()):
                UB_idxs = set(self.MuDepIndices_UB)-set(keepParameters['UB'])
        ## Pass new matrix and indices of elements to update to LP.updateMatrix-method ##
#        self.LP.updateMatrix(matrix=inputMatrix, Ainds=A_idxs, Binds=B_idxs,
#                             CTinds=[], LBinds=LB_idxs, UBinds=UB_idxs)
        self.LP.updateMatrix(matrix=inputMatrix, Ainds=A_idxs, Binds=[],
                             CTinds=[], LBinds=[], UBinds=[])

#        for coeff in ModelStructure.MuDependencies:
#            if coeff in list(ModelStructure.EnzymeConstraintsInfo.Elements.keys()):
#                parameters = ModelStructure.EnzymeConstraintsInfo.Elements[coeff]['CapacityParameter']
#                rowName = coeff
#                colName = ModelStructure.EnzymeConstraintsInfo.Elements[coeff]['AssociatedEnzyme']
#                value = determineCoefficient(x=parameters, changes={'Mu': float(Mu)}, species='Mu')
#                self.setProblemCoefficients(inputDict={rowName: {colName: -value}}, logging=logging)
#            elif coeff in list(ModelStructure.ProcessConstraintsInfo.Elements.keys()):
#                parameters = ModelStructure.ProcessConstraintsInfo.Elements[coeff]['CapacityParameter']
#                rowName = coeff
#                colName = str(
#                    ModelStructure.ProcessConstraintsInfo.Elements[coeff]['AssociatedProcess'])+str('_machinery')
#                capa = determineCoefficient(x=parameters, changes={'Mu': float(Mu)}, species='Mu')
#                ProcessInQuestion = ModelStructure.ProcessConstraintsInfo.Elements[coeff]['AssociatedProcess']
#                ProcessDict = {ModelStructure.ProcessInfo.Elements[i]['ID']: i for i in list(
#                    ModelStructure.ProcessInfo.Elements.keys())}
#                InfoOnProcess = ModelStructure.ProcessInfo.Elements[ProcessDict[ProcessInQuestion]]
#                eigencost = 0
#                for subunit in list(InfoOnProcess['Composition'].keys()):
#                    if subunit in list(ModelStructure.ProteinInfo.Elements.keys()):
#                        req = ModelStructure.ProteinInfo.Elements[subunit]['ProcessRequirements']
#                        if ProcessDict[ProcessInQuestion] in list(req.keys()):
#                            eigencost += (req[ProcessDict[ProcessInQuestion]] *
#                                          InfoOnProcess['Composition'][subunit])
#                netcapacity = capa-float(Mu)*eigencost
#                self.setProblemCoefficients(
#                    inputDict={rowName: {colName: -netcapacity}}, logging=logging)
#            elif coeff in list(ModelStructure.DensityConstraintsInfo.Elements.keys()):
#                parameters = ModelStructure.DensityConstraintsInfo.Elements[coeff]['CapacityParameter']
#                rowName = coeff
#                value = determineCoefficient(x=parameters, changes={'Mu': float(Mu)}, species='Mu')
#                self.setRighthandSideValue(inputDict={rowName: value}, logging=logging)

        self.BuildClassicMatrix(Mu=float(Mu))
        ## Construct new LP from solver object ##
        inputMatrix = RBA_Matrix()
        inputMatrix.loadMatrix(matrix=self.ClassicRBAmatrix)
        process_Inds = [(i.split('_machinery')[0]+'_capacity', i) for i in self.Processes]
#        process_Inds = [('P_TA_capacity', 'P_TA_machinery'), ('P_CHP_capacity',
#                                                              'P_CHP_machinery'), ('P_SEC_capacity', 'P_SEC_machinery')]
        self.LP.updateMatrix(matrix=inputMatrix, Ainds=process_Inds, Binds=B_idxs,
                             CTinds=[], LBinds=[], UBinds=[])

    def setMu(self, Mu, ModelStructure, keepParameters=None, logging=True):
        """
        Changes growth-rate of problem and sets all associated coefficients.

        Can be provided with 'keepParameters' argument
        to define problem coefficients which should remain unchanged.

        Parameters
        ----------
        Mu : float
            Growth rate
        ModelStructure : RBA_ModellStructure object.
        keepParameters : dict
            Dictionary indicating which elements of the linear problem should
            not be affected when setting growth-rate. Possible keys of dictionary:
                'LHS': List of index-tuples (constraint-ID,variable-ID),
                       indicating elements of the lefthandside (constraint-matrix).
                'RHS': List of constraint IDs, indicating elements of the righthandside (b-vector).
                'LB': List of variable-IDs, indicating elements of the lower-bound vector.
                'UB': List of variable-IDs, indicating elements of the upper-bound vector.
            Default: keepParameters=None
        logging : bool
            Wheter to write change to log-file or not
                """

        self.Mu = float(Mu)

        if self.LP.row_names == self.ClassicRBAmatrix.row_names and self.LP.col_names == self.ClassicRBAmatrix.col_names:
            self.setMu_classicProblem(Mu=float(Mu), keepParameters=keepParameters, logging=logging)
        else:
            #self.setMu_classicProblem(Mu=float(Mu), keepParameters=keepParameters, logging=logging)
            self.setMu_modifiedProblem(Mu=float(Mu), ModelStructure=ModelStructure,
                                       keepParameters=keepParameters, logging=logging)

        self.resetCPLEXparams()

    def resetCPLEXparams(self):
        """
        Sets cplex.Cplex() parameters to predefined values.

        Settings:
            cplex.Cplex.parameters.feasopt.tolerance.set(1e-9)
            cplex.Cplex.parameters.simplex.tolerances.feasibility.set(1e-9)
            cplex.Cplex.parameters.simplex.tolerances.optimality.set(1e-9)
            cplex.Cplex.parameters.simplex.tolerances.markowitz.set(0.1)
            cplex.Cplex.parameters.barrier.convergetol.set(1e-9)
            cplex.Cplex.parameters.read.scale.set(1)
            cplex.Cplex.set_results_stream(None)
            cplex.Cplex.set_log_stream(None)
            cplex.Cplex.set_warning_stream(None)
        """
        self.LP.cplexLP.parameters.feasopt.tolerance.set(1e-9)
        self.LP.cplexLP.parameters.simplex.tolerances.feasibility.set(1e-9)
        self.LP.cplexLP.parameters.simplex.tolerances.optimality.set(1e-9)
        self.LP.cplexLP.parameters.simplex.tolerances.markowitz.set(0.1)
        self.LP.cplexLP.parameters.barrier.convergetol.set(1e-9)
        self.LP.cplexLP.parameters.read.scale.set(1)
        self.LP.cplexLP.set_results_stream(None)
        self.LP.cplexLP.set_log_stream(None)
        self.LP.cplexLP.set_warning_stream(None)

    def getConstraintType(*constraints):
        """
        Extracts type of constraints.

        Parameters
        ----------
        constraints : str or list of str or None
            Constraints to retreive the type for.
            Either constraint ID or list of constraint IDs to specify the type
            of which constraint to look up.
            This is an optional input; if not provided all constraint are looked up.
        Returns
        ----------
        dict
            Dictionary with constraint-IDs as keys and
            type identification-characters as values.
        """
        if len(list(constraints)) > 0:
            if isinstance(constraints[0], list):
                names = constraints[0]
            if isinstance(constraints[0], str):
                names = [constraints[0]]
        if len(list(constraints)) == 0:
            names = self.LP.row_names
        Types = [self.LP.cplexLP.linear_constraints.get_senses(
            self.LP.rowIndicesMap[c]) for c in names]
        return(dict(zip(names, Types)))

    def setConstraintType(inputDict, logging=True):
        """
        Sets type of constraints.

        E: = or L: <=

        Parameters
        ----------
        inputDict : dict
            Dictionary with constraint-IDs as keys and type identification-character as values.
            ({'col1':'E','col2':'L', ...}).
        logging : bool
            Wheter to write change to log-file or not
        """
        if logging:
            for c in list(inputDict.keys()):
                self.LogBook.addEntry('Constraint-type {} changed:{} --> {}'.format(c,
                                                                                    self.LP.row_signs[self.LP.rowIndicesMap[c]], inputDict[c]))

        ##Update in cplex.Cplex LP##
        self.LP.cplexLP.linear_constraints.set_senses(
            list(zip(inputDict.keys(), inputDict.values())))
        ##Transfer changes to rbastructure.RBA_LP object##
        self.LP.row_signs = self.LP.cplexLP.linear_constraints.get_senses()

    def getObjective(self, *variables):
        """
        Returns objective coefficient of variables in problem.

        Parameters
        ----------
        variables : str or list of str or None
            Variables to retreive the objective coefficients for.
            Either variable ID or list of variable IDs to specify
            the coefficients of which variables to look up.
            This is an optional input; if not provided all variables are looked up.

        Returns
        ----------
        dict
            Dictionary with variable-IDs as keys and
            objective coefficients as values.
        """

        if len(list(variables)) > 0:
            if isinstance(variables[0], list):
                vrs = variables[0]
                vls = []
            if isinstance(variables[0], str):
                vrs = [variables[0]]
                vls = []
            for v in vrs:
                vls.append(self.LP.cplexLP.objective.get_linear()[
                           numpy.where(numpy.array(self.LP.col_names) == v)[0][0]])
        if len(list(variables)) == 0:
            vrs = self.LP.col_names
            vls = self.LP.cplexLP.objective.get_linear()
        OF = dict(zip(vrs, vls))
        return(OF)

    def setObjectiveCoefficients(self, inputDict, logging=True):
        """
        Sets objective function coefficients.

        Parameters
        ----------
        inputDict : dict
            Dictionary with variable-IDs as keys and new numeric values as values.
            ({'col1':42,'col2':9000, ...}).
        logging : bool
            Wheter to write change to log-file or not
        """
        if logging:
            for v in list(inputDict.keys()):
                self.LogBook.addEntry('Objective coefficient {} changed:{} --> {}'.format(v,
                                                                                          self.LP.f[self.LP.colIndicesMap[v]], inputDict[v]))

        ##Update in cplex.Cplex LP##
        self.LP.cplexLP.objective.set_linear(
            zip(list(inputDict.keys()), [float(i) for i in list(inputDict.values())]))
        ##Transfer changes to rbastructure.RBA_LP object##
        self.LP.f = self.LP.cplexLP.objective.get_linear()

    def clearObjective(self, logging=True):
        """
        Sets all coefficients of the objective function to zero.

        Parameters
        ----------
        logging : bool
            Wheter to write change to log-file or not
        """

        if logging:
            self.LogBook.addEntry('Objective cleared.')
        ##Update in cplex.Cplex LP##
        self.LP.cplexLP.objective.set_linear(
            zip(self.LP.col_names, [float(0)]*len(self.LP.col_names)))
        ##Transfer changes to rbastructure.RBA_LP object##
        self.LP.f = self.LP.cplexLP.objective.get_linear()

    def invertObjective(self, logging=True):
        """
        Changes sign (optimisation-sense) of objective function.

        Parameters
        ----------
        logging : bool
            Wheter to write change to log-file or not
        """

        if logging:
            self.LogBook.addEntry('Objective inverted.')
        current_objective = self.getObjective()
        negative_objective = zip(list(current_objective.keys()),
                                 [-float(x) for x in list(current_objective.values())])
        ##Update in cplex.Cplex LP##
        self.LP.cplexLP.objective.set_linear(negative_objective)
        ##Transfer changes to rbastructure.RBA_LP object##
        self.LP.f = self.LP.cplexLP.objective.get_linear()

    def getRighthandSideValue(self, *constraints):
        """
        Extracts coefficients of problem's righthand side (B-vector).

        Parameters
        ----------
        constraints : str or list of str or None
            Constraints to retreive the objective coefficients for.
            Either constraint ID or list of constraint IDs to specify the RHS
            of which constraint to look up.
            This is an optional input; if not provided all constraint are looked up.
        Returns
        ----------
        dict
            Dictionary with constraint-IDs as keys and
            RHS-values as values.
        """

        if len(list(constraints)) > 0:
            if isinstance(constraints[0], list):
                names = constraints[0]
            if isinstance(constraints[0], str):
                names = [constraints[0]]
        if len(list(constraints)) == 0:
            names = self.LP.row_names
        Bs = [self.LP.cplexLP.linear_constraints.get_rhs(self.LP.rowIndicesMap[c]) for c in names]
        return(dict(zip(names, Bs)))

    def setRighthandSideValue(self, inputDict, logging=True):
        """
        Set coefficients of the problems' RHS (b-vector).
        Parameters
        ----------
        inputDict : dict
            Dictionary with constraint-IDs as keys and new numeric values as values.
            ({'row1':42,'row2':9000, ...}).
        logging : bool
            Wheter to write change to log-file or not
        """
        if logging:
            for c in list(inputDict.keys()):
                self.LogBook.addEntry('RHS coefficient {} changed:{} --> {}'.format(v,
                                                                                    self.LP.b[self.LP.rowIndicesMap[c]], inputDict[c]))

        ##Update in cplex.Cplex LP##
        self.LP.cplexLP.linear_constraints.set_rhs(
            list(zip(list(inputDict.keys()), [float(i) for i in list(inputDict.values())])))
        ##Transfer changes to rbastructure.RBA_LP object##
        self.LP.b = self.LP.cplexLP.linear_constraints.get_rhs()

    def calculateLefthandSideValue(self, *constraints):
        """
        Calculates value of problem's lefthand side
        (after multiplying with solution-vector).

        Parameters
        ----------
        constraints : str or list of str or None
            Constraints to retreive the LHS-value for.
            Either constraint ID or list of constraint IDs to specify the LHS-value
            of which constraint to look up.
            This is an optional input; if not provided all constraint are looked up.
        Returns
        ----------
        dict
            Dictionary with constraint-IDs as keys and
            LHS-value as values.
        """
        if len(list(constraints)) > 0:
            if isinstance(constraints[0], list):
                names = constraints[0]
            if isinstance(constraints[0], str):
                names = [constraints[0]]
        if len(list(constraints)) == 0:
            names = self.LP.row_names
        Sol = numpy.array(list(self.SolutionValues.values()))
        Amat = self.LP.A.toarray()
        multRes = Amat.dot(Sol)
        out = {c: multRes[self.LP.row_names.index(c)] for c in names}
        return(out)

    def getProblemCoefficients(self, *inputTuples):
        """
        Returns coefficients of LHS of problem.

        Parameters
        ----------
        inputTuples : tuple or list of tuples.
            Tuples hold row and column indices.
            [('row1','col1'),('row2','col2'),...] or ('row1','col1').

        Returns
        ----------
        dict
            Dictionary with index tuples as keys and
            matrix coefficients as values.
        """

        if len(list(inputTuples)) > 0:
            if isinstance(inputTuples[0], list):
                tuples = inputTuples[0]
            elif isinstance(inputTuples[0], tuple):
                tuples = [inputTuples[0]]
            else:
                print('Error: Please provide tuple or list of tuples as input')
                return
        if len(list(inputTuples)) == 0:
            print('Error: Please provide tuple or list of tuples as input')
            return
        return(dict(zip(tuples, [self.LP.cplexLP.linear_constraints.get_coefficients(self.LP.rowIndicesMap[tup[0]], self.LP.colIndicesMap[tup[1]]) for tup in tuples])))

    def setProblemCoefficients(self, inputDict, logging=True):
        """
        Set coefficients of the problems' LHS (constraint matrix).

        Parameters
        ----------
        inputDict : dict
            Dict with index-tuples ('row1','col1') as keys
            and new numeric values as values.
            ({('row1','col1'):42,('row2','col2'):9000, ...}).
        logging : bool
            Wheter to write change to log-file or not
        """

        variables = []
        constraints = []
        coefficients = []
        Changes = []
        for const in list(inputDict.keys()):
            for var in list(inputDict[const].keys()):
                if logging:
                    self.LogBook.addEntry('LHS coefficient {} changed:{} --> {}'.format(
                        (const, var), self.getProblemCoefficients((const, var))[(const, var)], inputDict[const][var]))
                constraints.append(const)
                variables.append(var)
                coefficients.append(numpy.float64(inputDict[const][var]))
        Changes = list(zip(constraints, variables, coefficients))
        ##Update in cplex.Cplex LP##
        self.LP.cplexLP.linear_constraints.set_coefficients(Changes)
        ##Transfer changes to rbastructure.RBA_LP object##
        self.LP.A = convertCPLEXmatrix_to_Sparse(self)

    def getUB(self, *variables):
        """
        Returns upper bounds of problem variables.

        Parameters
        ----------
        variables : str list of str or None
            Variables to retreive the objective coefficients for.
            Either variable ID or list of variable IDs to specify
            the coefficients of which variables to look up.
            This is an optional input; if not provided all variables are looked up.

        Returns
        ----------
        dict
            Dictionary with variable-IDs as keys and
            upper bounds as values.
        """

        if len(list(variables)) > 0:
            if isinstance(variables[0], list):
                names = variables[0]
            if isinstance(variables[0], str):
                names = [variables[0]]
        if len(list(variables)) == 0:
            names = self.LP.col_names
        bound = [self.LP.cplexLP.variables.get_upper_bounds(
            self.LP.colIndicesMap[v]) for v in names]
        return(dict(zip(names, bound)))

    def getLB(self, *variables):
        """
        Returns lower bounds of problem variables.

        Parameters
        ----------
        variables : str list of str or None
            Variables to retreive the objective coefficients for.
            Either variable ID or list of variable IDs to specify
            the coefficients of which variables to look up.
            This is an optional input; if not provided all variables are looked up.

        Returns
        ----------
        dict
            Dictionary with variable-IDs as keys and
            lower bounds as values.
        """

        if len(list(variables)) > 0:
            if isinstance(variables[0], list):
                names = variables[0]
            if isinstance(variables[0], str):
                names = [variables[0]]
        if len(list(variables)) == 0:
            names = self.LP.col_names
        bound = [self.LP.cplexLP.variables.get_lower_bounds(
            self.LP.colIndicesMap[v]) for v in names]
        return(dict(zip(names, bound)))

    def setLB(self, inputDict, logging=True):
        """
        Set lower-bounds of the problem variables.

        Parameters
        ----------
        inputDict : dict
            Dictionary with variable-IDs as keys and new numeric values as values.
            ({'col1':42,'col2':9000, ...}).
        logging : bool
            Wheter to write change to log-file or not
        """

        if logging:
            for v in list(inputDict.keys()):
                self.LogBook.addEntry('LowerBound {} changed:{} --> {}'.format(v,
                                                                               self.LP.LB[self.LP.colIndicesMap[v]], inputDict[v]))
        ##Update in cplex.Cplex LP##
        self.LP.cplexLP.variables.set_lower_bounds(
            zip(list(inputDict.keys()), [float(i) for i in list(inputDict.values())]))
        ##Transfer changes to rbastructure.RBA_LP object##
        self.LP.LB = self.LP.cplexLP.variables.get_lower_bounds()

    def setUB(self, inputDict, logging=True):
        """
        Set upper-bounds of the problem variables.

        Parameters
        ----------
        inputDict : dict
            Dictionary with variable-IDs as keys and new numeric values as values.
            ({'col1':42,'col2':9000, ...}).
        logging : bool
            Wheter to write change to log-file or not
        """

        if logging:
            for v in list(inputDict.keys()):
                self.LogBook.addEntry('UpperBound {} changed:{} --> {}'.format(v,
                                                                               self.LP.UB[self.LP.colIndicesMap[v]], inputDict[v]))
        ##Update in cplex.Cplex LP##
        # self.LP.cplexLP.variables.set_upper_bounds(zip(list(inputDict.keys()), [self.LP.colIndicesMap[v] for v in list(inputDict.keys())]))
        self.LP.cplexLP.variables.set_upper_bounds(
            zip(list(inputDict.keys()), [float(i) for i in list(inputDict.values())]))
        ##Transfer changes to rbastructure.RBA_LP object##
        self.LP.UB = self.LP.cplexLP.variables.get_upper_bounds()

    def getEnzymeCapacities(self, *Enzyme):
        """
        Returns of capacity coefficients oof enzymes in problem.

        Parameters
        ---------
        Enzyme : str list of str or None
            Enzymes to retreive the efficiency coefficients for.
            Either enzyme ID or list of enzyme IDs to specify
            the efficiencies of which enzymes to look up.
            This is an optional input; if not provided all enzymes are looked up.

        Returns
        ----------
        dict
            Dictionary with enzyme ID as keys and dictionaries as values.
            The "inner" dictionaries hold keys 'Forward' and 'Backward' with forward and backward efficiencies as numeric values respectively.
            {'R_abcdef_enzyme':{'Froward':3000000,'Backward':3000000}} (negative values).
            If enzyme is irreversible 'Backward' is numpy.nan.
        """

        import difflib
        if len(list(Enzyme)) > 0:
            if isinstance(Enzyme[0], list):
                EnzymesInQuestion = Enzyme[0]
            if isinstance(Enzyme[0], str):
                EnzymesInQuestion = [Enzyme[0]]
        if len(list(Enzyme)) == 0:
            EnzymesInQuestion = self.Enzymes
        A = scipy.sparse.lil_matrix(self.LP.A)
        out = {}
        for i in EnzymesInQuestion:
            likelyFWkapps = difflib.get_close_matches(i, self.Enzyme_FWcapacities, 1)
            if len(likelyFWkapps) > 0:
                FWkapp = likelyFWkapps[0]
                BWkapp = FWkapp.replace('_forward_capacity', '_backward_capacity')
                FW = A[self.LP.rowIndicesMap[FWkapp], self.LP.colIndicesMap[i]]
                if BWkapp in self.Enzyme_BWcapacities:
                    BW = A[self.LP.rowIndicesMap[BWkapp], self.LP.colIndicesMap[i]]
                else:
                    BW = numpy.nan
            else:
                BW = numpy.nan
                FW = numpy.nan
            out.update({i: {'Forward': FW, 'Backward': BW}})
        return(out)

    # def setEnzymeCapacities(self, inputDict):

    def getProcessCapacities(self, *Process):
        """
        Returns of capacity coefficients of process machineries in problem.

        Parameters
        ---------
        Process : str list of str or None
            Processes to retreive the efficiency coefficients for.
            Either process ID or list of process IDs to specify
            the efficiencies of which processes to look up.
            This is an optional input; if not provided all processes are looked up.

        Returns
        ----------
        dict
            Dictionary with process ID as keys and their efficiencies as values.
        """

        import difflib
        if len(list(Process)) > 0:
            if isinstance(Process[0], list):
                ProcessesInQuestion = Process[0]
            if isinstance(Process[0], str):
                ProcessesInQuestion = [Process[0]]
        if len(list(Process)) == 0:
            ProcessesInQuestion = self.Processes
        IndexPairs = [(self.LP.rowIndicesMap[difflib.get_close_matches(i, self.ProcessCapacities, 1)[
            0]], self.LP.colIndicesMap[i]) for i in ProcessesInQuestion]
        A = scipy.sparse.lil_matrix(self.LP.A)
        out = {i: A[IndexPairs[ProcessesInQuestion.index(i)]] for i in ProcessesInQuestion}
        return(out)

    # def setProcessCapacities(self, inputDict):

    def getCompartmentCapacities(self, *Compartment):
        """
        Returns of capacity coefficients of compartments in problem.

        Parameters
        ---------
        Compartment : str list of str or None
            Compartments to retreive the capacity coefficients for.
            Either compartment ID or list of compartment IDs to specify
            the capacities of which compartments to look up.
            This is an optional input; if not provided all compartments are looked up.

        Returns
        ----------
        dict
            Dictionary with compartments ID as keys and their capacities as values.
        """

        if len(list(Compartment)) > 0:
            if isinstance(Compartment[0], list):
                CompartmentsInQuestion = Compartment[0]
            if isinstance(Compartment[0], str):
                CompartmentsInQuestion = [Compartment[0]]
        if len(list(Compartment)) == 0:
            CompartmentsInQuestion = self.CompartmentDensities
        out = {i: self.LP.b[self.LP.rowIndicesMap[i]] for i in CompartmentsInQuestion}
        return(out)

    # def setCompartmentCapacities(self, inputDict):


def convertCPLEXmatrix_to_Sparse(inputStructure):
    import scipy
    Ma = inputStructure.LP.cplexLP.linear_constraints.get_rows()
    Anew = numpy.zeros((inputStructure.LP.cplexLP.linear_constraints.get_num(),
                        inputStructure.LP.cplexLP.variables.get_num()))
    rowIndex = 0
    for m in Ma:
        Anew[rowIndex, m.ind] = m.val
        rowIndex += 1
    return(scipy.sparse.coo_matrix(Anew))


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
            X = float(min(max(float(changes[species]), X_mini), X_maxi))
            Y = float(A*X + C)
            result = float(min(max(Y, Y_mini), Y_maxi))

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
            Y = float(VM * float(changes[species]) / (float(changes[species]) + KM))
            result = float(max(Y, Y_mini))

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
            result = max(Y, Y_mini)

        if type == 'inverse':
            C = 1
            if 'C' in list(pars.keys()):
                C = numpy.float64(pars['C'])
            try:
                result = C/float(changes[i])
            except KeyError:
                print('variable is 0, impossible to do inversion')

        multiplicativeFactors.append(result)
    value = numpy.prod(numpy.array(multiplicativeFactors))
    return(float(value))
