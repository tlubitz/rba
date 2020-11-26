# python 2/3 compatibility
from __future__ import division, print_function

import copy
import json
import pandas

from rbastructure.information_block import InformationBlock


class StatisticsBlock(InformationBlock):
    """
    Class holding model-statistics.

    Brief summary of key-numbers of model.

    Attributes
    ----------
    Elements : Dictionary different numbers on model.

    """

    def derive(self, Dict):
        self.Elements = Dict

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
#                    TableOut.ix[i,'Measure']=Var.replace('"','')
#                    TableOut.ix[i,'Value']=Val.replace('"','')
                if "'" in Var:
                    Var = Var.replace("'", "")
                if "'" in Val:
                    Val = Val.replace("'", "")
                TableOut.ix[i, 'Measure'] = Var
                TableOut.ix[i, 'Value'] = Val
# WORKS                    TableOut.ix[i,'Measure']='"'+Var+'"'
# WORKS                    TableOut.ix[i,'Value']='"'+Val+'"'
            if len(list(NameList)) == len(list(TableOut)):
                TableOut.columns = list(NameList)
            return TableOut
        if len(list(Block.keys())) == 0:
            return pandas.DataFrame()

#     def toSBtab(self,document_name, table_type, table_name,document, unit):
#          from sbtab import SBtab
#          DF=self.toDataFrame()
#          SB=SBtab.SBtabTable.from_data_frame(DF, document_name, table_type, table_name, document, unit, sbtab_version='1.0')
#          SB.write(table_type)
#          f=open(document_name+'.tsv','w')
#          f.write(SB.table_string)
#          f.close()


def JSON_Int64_compensation(o):
    if isinstance(o, numpy.int64):
        return int(o)
    raise TypeError
