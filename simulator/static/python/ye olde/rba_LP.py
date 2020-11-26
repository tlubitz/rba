from __future__ import division, print_function
import numpy
import scipy
import cplex
import copy
import itertools
from rbastructure.rba_Matrix import RBA_Matrix


class RBA_LP(RBA_Matrix):
    """

    Attributes
    ----------
    A : scipy.sparse.coo_matrix
        Lefthandside of constraint Matrix (aka Constraint Matrix)
    b : numpy.array
        Righthandside of Constraint Matrix
    row_signs : list
        Type of constraints ('E' for equality 'L' for lower-or-equal inequality) --> Ax=b or Ax<=b
    f : numyp.array
        Objective function linear get_coefficients
    LB : numpy.array
        Lower bounds of decision-variables
    UB : numpy.array
        Upper bounds of decision-variables
    row_names : list
        Names of constraints
    col_names : list
        Names of decision-variables
    rowIndicesMap : dict
        Dictionary mapping constraint names to their numeric index (generated automatically)
    colIndicesMap : dict
        Dictionary mapping variable names to their numeric index (generated automatically)
    cplex: bool
        Boolean wheter a CPLEX problem has been created
    cplexLP: cplex.CPLEX
        CPLEx object to optimise with cplex

    Methods
    ----------
    __init__(*Matrix)

    updateMatrix(matrix, Ainds=None, Binds=None, CTinds=None, LBinds=None, UBinds=None)
        Overwrites coefficients with new values from argument 'matrix'.

    addMatrix(matrix)
        Merges the Problem with the input one.

    buildCPLEX_LP()
        Constructs a CPLEX-compliant LP-object.
    """

    def __init__(self, *Matrix):
        """
        If provided with an RBA_Matrix object as argument, assumes its attributes.
        If not provided with an RBA_Matrix object, initiates as empty RBA_LP object.
        """
        if isinstance(Matrix, RBA_Matrix):
            self.loadMatrix(matrix=Matrix)
            self.mapIndices()
        else:
            RBA_Matrix.__init__(self)
            self.mapIndices()
        self.cplex = False

    def updateMatrix(self, matrix, Ainds=None, Binds=None, CTinds=None, LBinds=None, UBinds=None):
        """
        Overwrites coefficients with new values from argument 'matrix'.

        Parameters
        ----------
        matrix : rbastructure.RBA_Matrix
            The matrix with elements to be added
        Ainds : list of tuples
            List of index-pair tuples [(row_coeff1,col_coeff1),(row_coeff2,col_coeff2),...],
            specifying which elements are updated.
            Default is None
            (then all elements of argument matrix (present in both matrices) are taken to update)
        Binds : list
            List of constraint-IDs for indices of the RHS, specifying which RHS Values are updated.
            Default is None (then all constraints (present in both matrices) are taken to update)
        CTinds : list
            List of constraint-IDs, specifying which row_signs are updated.
            Default is None (then all constraint types (present in both matrices) are taken to update)
        LBinds : list
            List of variable-IDs, specifying which lower-bounds values are updated.
            Default is None (then all variables (present in both matrices) are taken to update)
        UBinds : list
            List of variable-IDs, specifying which upper-bounds values are updated.
            Default is None (then all variables (present in both matrices) are taken to update)
        """

        ## Check if indices to update are provided, if not use all indices ##
        if Ainds is None:
            Ainds = list(zip(matrix.row_names, matrix.col_names))
        if Binds is None:
            Binds = matrix.row_names
        if CTinds is None:
            CTinds = matrix.row_names
        if LBinds is None:
            LBinds = matrix.col_names
        if UBinds is None:
            UBinds = matrix.col_names

        ####### Update constraint-matrix LHS (A) #######
        if matrix.row_names == self.row_names and matrix.col_names == self.col_names:
            ## If old and new matrix have same elements and order of indices ##
            x1, x2 = zip(*Ainds)
            #Find numeric indices of elements to update#
            rowsOld = [self.rowIndicesMap[i] for i in x1]
            colsOld = [self.colIndicesMap[i] for i in x2]
            newA = scipy.sparse.lil_matrix(matrix.A)
            oldA = scipy.sparse.lil_matrix(self.A)
            #Overwrite old elements at indices with corresponding elements from new matrix#
            oldA[tuple(rowsOld), tuple(colsOld)] = newA[tuple(rowsOld), tuple(colsOld)]
        else:  # If old and new matrix do not have same elements and order of indices ##
            ## Find elements (index pairs) which are in the old, as well in th new matrix. ##
            #            intersectionIndices = [(i[0], i[1]) for i in Ainds if i[0] in matrix.row_names and i[0]in self.row_names and i[1] in matrix.col_names and i[1] in self.col_names]
            #            oldIndPairs = set(itertools.product(self.row_names, self.col_names))
            #            newIndPairs = set(itertools.product(matrix.row_names, matrix.col_names))
            #            intersectionIndices = list(set(Ainds).intersection(oldIndPairs, newIndPairs))
            #            x1, x2 = zip(*intersectionIndices)
            x1, x2 = zip(*Ainds)
            ## Find the numeric indices of the intersecting elements to update in both matrices ##
            rowsOld = [self.rowIndicesMap[i] for i in x1]
            colsOld = [self.colIndicesMap[i] for i in x2]
            rowsNew = [matrix.rowIndicesMap[i] for i in x1]
            colsNew = [matrix.colIndicesMap[i] for i in x2]
            newA = scipy.sparse.lil_matrix(matrix.A)
            oldA = scipy.sparse.lil_matrix(self.A)
            #Overwrite old elements at indices with corresponding elements from new matrix#
            oldA[tuple(rowsOld), tuple(colsOld)] = newA[tuple(rowsNew), tuple(colsNew)]
        self.A = scipy.sparse.coo_matrix(oldA)

        ## Update RHS (b)##
        if matrix.row_names == self.row_names:
            ## If old and new LPs have same rows and row-order ##
            #Find numeric indices of rows to update (same for old and new matrix)#
            rowsNew = [matrix.rowIndicesMap[i] for i in Binds]
            #Overwrite old elements at row-indices with corresponding new elements#
            self.b[rowsNew] = matrix.b[rowsNew]
        else:
            x = [i for i in Binds if i in matrix.row_names and i in self.row_names]
            #Find numeric indices of rows to update (for old and new matrix individually)#
            rowsNew = [matrix.rowIndicesMap[i] for i in x]
            rowsOld = [self.rowIndicesMap[i] for i in x]
            #Overwrite old elements at row-indices with corresponding new elements#
            self.b[rowsOld] = matrix.b[rowsNew]

        ## Update Constraint type ##
        if matrix.row_names == self.row_names:
            rowsNew = [matrix.rowIndicesMap[i] for i in CTinds]
            self.b[rowsNew] = matrix.b[rowsNew]
        else:
            x = [i for i in CTinds if i in matrix.row_names and i in self.row_names]
            rowsNew = [matrix.rowIndicesMap[i] for i in x]
            rowsOld = [self.rowIndicesMap[i] for i in x]
            RSign = numpy.array(self.row_signs)
            RSign[rowsOld] = numpy.array(matrix.row_signs)[rowsNew]
            self.row_signs = list(RSign)

        ## Update LB##
        oLB = numpy.array(self.LB)
        nLB = numpy.array(matrix.LB)
        if matrix.col_names == self.col_names:
            colsNew = [matrix.colIndicesMap[i] for i in LBinds]
            oLB[colsNew] = nLB[colsNew]
        else:
            x = [i for i in LBinds if i in matrix.col_names and i in self.col_names]
            colsOld = [self.colIndicesMap[i] for i in x]
            colsNew = [matrix.colIndicesMap[i] for i in x]
            oLB[colsOld] = nLB[colsNew]
        self.LB = oLB

        ## Update UB##
        oUB = numpy.array(self.UB)
        nUB = numpy.array(matrix.UB)
        if matrix.col_names == self.col_names:
            colsNew = [matrix.colIndicesMap[i] for i in UBinds]
            oUB[colsNew] = nUB[colsNew]
        else:
            x = [i for i in UBinds if i in matrix.col_names and i in self.col_names]
            colsOld = [self.colIndicesMap[i] for i in x]
            colsNew = [matrix.colIndicesMap[i] for i in x]
            oUB[colsOld] = nUB[colsNew]
        self.UB = oUB

        self.buildCPLEX_LP()

    def addMatrix(self, matrix):
        """
        Merges the Problem with the one provided as input-argument to this method.

        Matrix elements unique to input-matrix are added.
        Elements occuring in both matrices are overwritten with the value from new matrix.
        Generates CPLEX problem from merged problem.

        Parameters
        ----------
        matrix : rbastructure.RBA_Matrix
            The matrix with elements to be added
        """
        matrix.mapIndices()
        oldA = copy.deepcopy(self.A.toarray())
        matrix.A = matrix.A.toarray().astype('float64')
        ## Determine union of old- and new matrix's row-names.##
        ## Make sure the new names, not present in the old matrix are added to the end.##
        ## Same thing also with column-names ##
        ## Initiate compound RBA matrix and adapt elements to the compound LP with the new dimensions ##
        compoundProblem = RBA_Matrix()
        compoundProblem.row_names = list(
            self.row_names + list(set(matrix.row_names)-set(self.row_names)))
        compoundProblem.col_names = list(
            self.col_names + list(set(matrix.col_names)-set(self.col_names)))
        compoundProblem.A = numpy.zeros(
            (len(compoundProblem.row_names), len(compoundProblem.col_names)))
        compoundProblem.AtoLiL()
        compoundProblem.b = numpy.zeros(len(compoundProblem.row_names))
        compoundProblem.row_signs = ['E']*len(compoundProblem.row_names)
        compoundProblem.f = numpy.zeros(len(compoundProblem.col_names))
        compoundProblem.LB = numpy.zeros(len(compoundProblem.col_names))
        compoundProblem.UB = numpy.zeros(len(compoundProblem.col_names))
        compoundProblem.mapIndices()
        ## Since it has been made sure that the indices present in the original problem, ##
        ## are present in the compound problem at the beginning in the exact same order; ##
        ## One can now just put the information of the original problem at the beginning ##
        compoundProblem.A[0:oldA.shape[0], 0:oldA.shape[1]] = oldA
        compoundProblem.b[0:oldA.shape[0]] = copy.deepcopy(self.b)
        compoundProblem.f[0:oldA.shape[1]] = copy.deepcopy(self.f)
        compoundProblem.LB[0:oldA.shape[1]] = copy.deepcopy(self.LB)
        compoundProblem.UB[0:oldA.shape[1]] = copy.deepcopy(self.UB)
        compoundProblem.row_signs[0:oldA.shape[0]] = copy.deepcopy(self.row_signs)

        ## Find numeric indices of new-matrix elements in the new matrix ##
        NewMatrixRowIndices = tuple([matrix.rowIndicesMap[i] for i in matrix.row_names])
        NewMatrixColIndices = tuple([matrix.colIndicesMap[i] for i in matrix.col_names])

        ## Find numeric indices of new-matrix elements in the compound matrix ##
        CompoundMatrixRowIndices = tuple([compoundProblem.rowIndicesMap[i]
                                          for i in matrix.row_names])
        CompoundMatrixColIndices = tuple([compoundProblem.colIndicesMap[i]
                                          for i in matrix.col_names])

        ## Transfer new-matrix elements to compound problem ##
        compoundProblem.A[CompoundMatrixRowIndices,
                          CompoundMatrixColIndices] = matrix.A[NewMatrixRowIndices, NewMatrixColIndices]
        compoundProblem.f[list(CompoundMatrixColIndices)] = matrix.f[list(NewMatrixColIndices)]
        compoundProblem.LB[list(CompoundMatrixColIndices)] = matrix.LB[list(NewMatrixColIndices)]
        compoundProblem.UB[list(CompoundMatrixColIndices)] = matrix.UB[list(NewMatrixColIndices)]
        compoundProblem.b[list(CompoundMatrixRowIndices)] = matrix.b[list(NewMatrixRowIndices)]

        for i in range(len(NewMatrixRowIndices)):
            compoundProblem.row_signs[CompoundMatrixRowIndices[i]
                                      ] = matrix.row_signs[NewMatrixRowIndices[i]]

        ## Overwrite old matrix with compound problem. And rebuild CPLEX-LP if there was one before##
        self.loadMatrix(compoundProblem)
        if self.cplex:
            self.buildCPLEX_LP()

    def buildCPLEX_LP(self):
        """
        Constructs a CPLEX-compliant LP-object.
        """
        lhs = self.A.tolil()
        rows = []
        for nz_ind, data in zip(lhs.rows, lhs.data):
            rows.append(cplex.SparsePair(nz_ind, data))
        # define problem
        cpxLP = cplex.Cplex()
        cpxLP.variables.add(obj=list(self.f), ub=list(self.UB),
                            lb=list(self.LB), names=list(self.col_names))
        cpxLP.linear_constraints.add(lin_expr=rows,
                                     rhs=self.b,
                                     senses=self.row_signs,
                                     names=self.row_names)
        cpxLP.objective.set_sense(cpxLP.objective.sense.minimize)
        cpxLP.parameters.feasopt.tolerance.set(1e-9)
        cpxLP.parameters.simplex.tolerances.feasibility.set(1e-9)
        cpxLP.parameters.simplex.tolerances.optimality.set(1e-9)
        cpxLP.parameters.simplex.tolerances.markowitz.set(0.01)
        cpxLP.parameters.barrier.convergetol.set(1e-9)
        cpxLP.parameters.read.scale.set(1)
        cpxLP.set_results_stream(None)
        cpxLP.set_log_stream(None)
        cpxLP.set_warning_stream(None)
        self.cplexLP = cpxLP
        self.cplex = True


def convertCPLEXmatrix_to_Sparse(inputStructure):
    import scipy
    Ma = inputStructure.cplexLP.linear_constraints.get_rows()
    Anew = numpy.zeros((inputStructure.cplexLP.linear_constraints.get_num(),
                        inputStructure.cplexLP.variables.get_num()))
    rowIndex = 0
    for m in Ma:
        Anew[rowIndex, m.ind] = m.val
        rowIndex += 1
    return(scipy.sparse.coo_matrix(Anew))
