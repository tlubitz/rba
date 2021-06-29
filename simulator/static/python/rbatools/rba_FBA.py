from __future__ import division, print_function
import numpy
import scipy
import copy
import cplex
from .rba_Matrix import RBA_Matrix
from .rba_LP import RBA_LP
from .rba_LogBook import RBA_LogBook


class RBA_FBA(object):
    def __init__(self, FBA):
        self.LP = FBA
        self.LP.buildCPLEX_LP()
        self.parsimonious = False

    def parsimonise(self):
        A = self.LP.A.toarray()
        reversibleRxns = [i for i in self.LP.col_names if not i.startswith(
            'R_EX_') and self.LP.LB[self.LP.col_names.index(i)] < 0]
        Anew = numpy.zeros((len(self.LP.row_names), len(reversibleRxns)))
        fNew = numpy.zeros(len(reversibleRxns))
        LBnew = numpy.zeros(len(reversibleRxns))
        UBnew = numpy.zeros(len(reversibleRxns))
        addIndexCount = 0
        newCols = []
        for i in reversibleRxns:
            colIndex = self.LP.col_names.index(i)
            LBold = self.LP.LB[colIndex]
            rxnCol = A[:, colIndex]
            Anew[:, addIndexCount] = -rxnCol
            fNew[addIndexCount] = 1
            LBnew[addIndexCount] = 0
            UBnew[addIndexCount] = -LBold
            newCols.append(str('Backward_'+i))
            self.LP.LB[colIndex] = 0
            addIndexCount += 1

        ReverseMatrix = RBA_Matrix()
        ReverseMatrix.A = scipy.sparse.coo_matrix(Anew)
        ReverseMatrix.b = self.LP.b
        ReverseMatrix.f = fNew
        ReverseMatrix.LB = LBnew
        ReverseMatrix.UB = UBnew
        ReverseMatrix.row_signs = self.LP.row_signs
        ReverseMatrix.row_names = self.LP.row_names
        ReverseMatrix.col_names = newCols
        ReverseMatrix.mapIndices()
        self.LP.addMatrix(matrix=ReverseMatrix)
        internalRxns = [i for i in self.LP.col_names if not i.startswith(
            'R_EX_') and not i.startswith('R_maintenance_') and not i.startswith('R_BIOMASS_')]
        for i in range(len(self.LP.col_names)):
            if self.LP.col_names[i] in internalRxns:
                self.LP.f[i] = 1
            else:
                self.LP.f[i] = 0
        self.LP.buildCPLEX_LP()
        self.parsimonious = True

    def solveLP(self, feasibleStatuses=[1]):
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
            Default: feasibleStatuses=[1]
        logging : bool
            Wheter to write change to log-file or not
        """

        self.Solved = False
        ## Solve cplex LP ##
        self.LP.cplexLP.solve()
        ## Determin solution-status ##
        self.SolutionStatus = self.LP.cplexLP.solution.get_status()
        ## Check if solution status is amongst acceptable ones ##
        if self.SolutionStatus in feasibleStatuses:
            ## Extract solution-data ##
            self.Solved = True
            self.ObjectiveValue = self.LP.cplexLP.solution.get_objective_value()
            SolVals = dict(zip(self.LP.col_names, self.LP.cplexLP.solution.get_values()))
            self.SolutionValues = {}
            for i in SolVals.keys():
                if not i.startswith('Backward_'):
                    if str('Backward_'+i) in SolVals.keys():
                        val = SolVals[i]-SolVals[str('Backward_'+i)]
                    else:
                        val = SolVals[i]
                    self.SolutionValues.update({i: val})
            self.DualValues = dict(
                zip(self.LP.row_names, self.LP.cplexLP.solution.get_dual_values()))

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

    def setConstraintType(inputDict):
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
        ##Update in cplex.Cplex LP##
        self.LP.cplexLP.linear_constraints.set_senses(
            list(zip(inputDict.keys(), inputDict.values())))
        ##Transfer changes to rbatools.RBA_LP object##
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

    def setObjectiveCoefficients(self, inputDict):
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
        ##Update in cplex.Cplex LP##
        self.LP.cplexLP.objective.set_linear(
            zip(list(inputDict.keys()), [float(i) for i in list(inputDict.values())]))
        ##Transfer changes to rbatools.RBA_LP object##
        self.LP.f = self.LP.cplexLP.objective.get_linear()

    def clearObjective(self):
        """
        Sets all coefficients of the objective function to zero.

        Parameters
        ----------
        logging : bool
            Wheter to write change to log-file or not
        """

        ##Update in cplex.Cplex LP##
        self.LP.cplexLP.objective.set_linear(
            zip(self.LP.col_names, [float(0)]*len(self.LP.col_names)))
        ##Transfer changes to rbatools.RBA_LP object##
        self.LP.f = self.LP.cplexLP.objective.get_linear()

    def invertObjective(self):
        """
        Changes sign (optimisation-sense) of objective function.

        Parameters
        ----------
        logging : bool
            Wheter to write change to log-file or not
        """

        current_objective = self.getObjective()
        negative_objective = zip(list(current_objective.keys()),
                                 [-float(x) for x in list(current_objective.values())])
        ##Update in cplex.Cplex LP##
        self.LP.cplexLP.objective.set_linear(negative_objective)
        ##Transfer changes to rbatools.RBA_LP object##
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

    def setRighthandSideValue(self, inputDict):
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
        ##Update in cplex.Cplex LP##
        self.LP.cplexLP.linear_constraints.set_rhs(
            list(zip(list(inputDict.keys()), [float(i) for i in list(inputDict.values())])))
        ##Transfer changes to rbatools.RBA_LP object##
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

    def setProblemCoefficients(self, inputDict):
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
                constraints.append(const)
                variables.append(var)
                coefficients.append(numpy.float64(inputDict[const][var]))
        Changes = list(zip(constraints, variables, coefficients))
        ##Update in cplex.Cplex LP##
        self.LP.cplexLP.linear_constraints.set_coefficients(Changes)
        ##Transfer changes to rbatools.RBA_LP object##
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

    def setLB(self, inputDict):
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

        ##Update in cplex.Cplex LP##
        self.LP.cplexLP.variables.set_lower_bounds(
            zip(list(inputDict.keys()), [float(i) for i in list(inputDict.values())]))
        ##Transfer changes to rbatools.RBA_LP object##
        self.LP.LB = self.LP.cplexLP.variables.get_lower_bounds()

    def setUB(self, inputDict):
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

        ##Update in cplex.Cplex LP##
        # self.LP.cplexLP.variables.set_upper_bounds(zip(list(inputDict.keys()), [self.LP.colIndicesMap[v] for v in list(inputDict.keys())]))
        self.LP.cplexLP.variables.set_upper_bounds(
            zip(list(inputDict.keys()), [float(i) for i in list(inputDict.values())]))
        ##Transfer changes to rbatools.RBA_LP object##
        self.LP.UB = self.LP.cplexLP.variables.get_upper_bounds()

    def exportEscherMap(self, Filename):
        import json
        with open(Filename, 'w') as fout:
            fout.write(json.dumps({i[2:]: self.SolutionValues[i]
                                   for i in self.SolutionValues.keys()}, indent=4))


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
