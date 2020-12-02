# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import json
import copy
import pandas
import numpy


class InformationBlock(object):
    """
    Class holding information on the model.

    Parent to block-classes in ModelStructure.

    """

    def fromDict(self, Dict):
        self.Elements = Dict

    def JSONize(self):
        Block = self.Elements
        block2 = copy.deepcopy(Block)
        for i in list(Block.keys()):
            if type(Block[i]) is dict:
                for j in list(Block[i].keys()):
                    block2[i][j] = json.dumps(Block[i][j], default=JSON_Int64_compensation)
            else:
                block2[i] = json.dumps(Block[i], default=JSON_Int64_compensation)
        return(block2)

    def toDataFrame(self):
        Block = self.Elements
        if len(list(Block.keys())) > 0:
            fields = list(Block[list(Block.keys())[0]].keys())
            if 'index' in fields:
                fields.remove('index')
            TableOut = pandas.DataFrame(index=list(Block.keys()), columns=fields)
            for i in list(Block.keys()):
                for j in fields:
                    #                    intString = json.dumps(Block[i][j], default=JSON_Int64_compensation)
                    TableOut.ix[i, j] = Block[i][j]
            return TableOut
        else:
            return pandas.DataFrame()

    def toDataFrame_SBtabCompatibility(self, NameList):
        Block = self.Elements
        if len(list(Block.keys())) > 0:
            fields = list(Block[list(Block.keys())[0]].keys())
            if 'index' in fields:
                fields.remove('index')
            TableOut = pandas.DataFrame(columns=fields)
            for i in list(Block.keys()):
                for j in fields:
                    entry = Block[i][j]
                    if isinstance(entry, str):
                        if "'" in entry:
                            entry = entry.replace("'", " ")
                        if "," in entry:
                            entry = entry.replace(",", " ")
                    intString = json.dumps(entry, default=JSON_Int64_compensation)
#                         TableOut.ix[i,j]=intString.replace('"','')
#                         TableOut.ix[i,j]=intString
# OLD WORKS                    TableOut.ix[i, j] = '"'+intString+'"'
                    if "'" in intString:
                        intString = intString.replace("'", "")
                    TableOut.ix[i, j] = intString
            if len(list(NameList)) == len(list(TableOut)):
                TableOut.columns = list(NameList)
            return TableOut
        if len(list(Block.keys())) == 0:
            return pandas.DataFrame()

    def toSBtab_forDoc(self, document_name, table_type, table_name, document, unit, *NameList):
        SBtab_Colnames = []
        if len(list(NameList)) > 0:
            if len(list(NameList[0])) > 0:
                SBtab_Colnames = list(NameList[0])
        from sbtab import SBtab
        DF = self.toDataFrame_SBtabCompatibility(SBtab_Colnames)
        return(SBtab.SBtabTable.from_data_frame(DF, table_type, document_name, table_name, document, unit, sbtab_version='1.0'))

    def toSBtab(self, document_name, table_type, table_name, document, unit):
        from sbtab import SBtab
        DF = self.toDataFrame_SBtabCompatibility()
        SB = SBtab.SBtabTable.from_data_frame(
            DF, document_name, table_type, table_name, document, unit, sbtab_version='1.0')
        SB.write(table_type+'.tsv')


def JSON_Int64_compensation(o):
    if isinstance(o, numpy.int64):
        return int(o)
    raise TypeError
