# python 2/3 compatibility
from __future__ import division, print_function

import copy
import json
import pandas

# package imports
from .information_block import InformationBlock


class DescriptionBlock(InformationBlock):
    """
    Class holding general information on model.

    Author, metabolic reconstruction etc...

    Attributes
    ----------
    Elements : Dictionary information on model.


    """

    def __init__(self):
        self.Elements = {}

    def fromFiles(self, File):
        for i in File.index.tolist():
            self.Elements[i] = File.ix[i][1]

    def addEntries(self, Dict):
        for i in Dict.keys():
            self.Elements[i] = Dict[i]

    def JSONize(self):
        Block = self.Elements
        block2 = copy.deepcopy(Block)
        for i in Block.keys():
            block2[i] = json.dumps(Block[i])
        return(block2)

    def toDataFrame(self):
        Block = self.Elements
        if len(list(Block.keys())) > 0:
            fields = ['Measure', 'Value']
            TableOut = pandas.DataFrame(columns=fields, index=list(Block.keys()))
            for i in list(Block.keys()):
                Var = json.dumps(i, default=JSON_Int64_compensation)
                Val = json.dumps(Block[i], default=JSON_Int64_compensation)
                TableOut.ix[i, 'Measure'] = Var
                TableOut.ix[i, 'Value'] = Val
            return TableOut
        if len(list(Block.keys())) == 0:
            return pandas.DataFrame()

    def toDataFrame_SBtabCompatibility(self, NameList):
        Block = self.Elements
        if len(list(Block.keys())) > 0:
            fields = ['Measure', 'Value']
            TableOut = pandas.DataFrame(columns=fields, index=list(Block.keys()))
            for i in list(Block.keys()):
                Var = json.dumps(i, default=JSON_Int64_compensation)
                Val = json.dumps(Block[i], default=JSON_Int64_compensation)
                TableOut.ix[i, 'Measure'] = '"'+Var.replace('"', '')+'"'
                TableOut.ix[i, 'Value'] = '"'+Val.replace('"', '')+'"'
            if len(list(NameList)) == len(list(TableOut)):
                TableOut.columns = list(NameList)
            return TableOut
        if len(list(Block.keys())) == 0:
            return pandas.DataFrame()

    def toDataFrame_RunInfo(self):
        Block = self.Elements
        if len(list(Block.keys())) > 0:
            runs = list(Block[list(Block.keys())[0]].keys())
            fields = ['Property']+runs
            TableOut = pandas.DataFrame(columns=fields, index=list(Block.keys()))
            for i in list(Block.keys()):
                Var = json.dumps(i, default=JSON_Int64_compensation)
#                    TableOut.ix[i,'Measure']=Var.replace('"','')
                if "'" in Var:
                    Var = Var.replace("'", "")
                TableOut.ix[i, 'Property'] = Var
                for j in runs:
                    Val = json.dumps(Block[i][j], default=JSON_Int64_compensation)
                    if "'" in Val:
                        Val = Val.replace("'", "")
                    TableOut.ix[i, j] = Val.replace('"', '')
#                        TableOut.ix[i,j]=Val
            return TableOut
        if len(list(Block.keys())) == 0:
            return pandas.DataFrame()

    def toSBtab_RunInfo_forDoc(self, document_name, table_type, table_name, document, unit, *NameList):
        SBtab_Colnames = []
        if len(list(NameList)) > 0:
            if len(list(NameList[0])) > 0:
                SBtab_Colnames = list(NameList[0])
        from sbtab import SBtab
#          print(SBtab_Colnames)
#          DF=self.toDataFrame_RunInfo(SBtab_Colnames)
        DF = self.toDataFrame_RunInfo()
        return(SBtab.SBtabTable.from_data_frame(DF, table_type, document_name, table_name, document, unit, sbtab_version='1.0'))

    def toSBtab_RunInfo(self, document_name, table_type, table_name, document, unit):
        from sbtab import SBtab
        DF = self.toDataFrame_RunInfo()
        SB = SBtab.SBtabTable.from_data_frame(
            DF, document_name, table_type, table_name, document, unit, sbtab_version='1.0')
        SB.write(table_type+'.tsv')

#     def toSBtab(self,document_name, table_type, table_name,document, unit):
#          from sbtab import SBtab
#          DF=self.toDataFrame()
#          SB=SBtab.SBtabTable.from_data_frame(DF, document_name, table_type, table_name, document, unit, sbtab_version='1.0')
#          SB.write(table_type)
#          f=open(document_name+'.tsv','w')
#          f.write(SB.table_string)
#          f.close()


def JSON_Int64_compensation(o):
    import numpy
    if isinstance(o, numpy.int64):
        return int(o)
    raise TypeError
