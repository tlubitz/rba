# python 2/3 compatibility
from __future__ import division, print_function
import numpy
import xml.etree.ElementTree as ET
import copy

# package imports
import rba

class Data_Controler(object):

     def __init__(self,StaticData,SimulationData):
          self.Structural=StaticData
          self.SimulationData=SimulationData


     def makeXML(self):
          if hasattr(self,'Structural'):
               if hasattr(self,'SimulationData'):
                    new=ReferencesBetweenInformation(self.Structural,self.SimulationData)
                    S=ET.fromstring(new['Structural'].toXML())
                    D=ET.fromstring(new['SimulationData'].toXML())[0]
                    S.insert(1,D)
               else:
                    S=ET.fromstring(new['Structural'].toXML())
          else:
               if hasattr(self,'SimulationData'):
                    S=ET.fromstring(new['SimulationData'].toXML())
          return(ET.tostring(S,'utf-8'))


#     def makeSBtab(self,*options)


#     def exportProteome(self)


#     def exportFluxome(self)


def ReferencesBetweenInformation(Structural,SimulationData):
      static=copy.deepcopy(Structural)
      dynamic=copy.deepcopy(SimulationData)
      for i in list(static.MetaboliteConstraintsInfo.Elements.keys()):
            if i in list(dynamic.MetaboliteConstraintData.Elements.keys()):
                  static.MetaboliteConstraintsInfo.Elements[i].update({'DataElement':'MetaboliteConstraintData##'+i+'_Data'})
                  dynamic.MetaboliteConstraintData.Elements[i].update({'StructureElement':'MetaboliteConstraint##'+i})
      for i in list(static.ProcessConstraintsInfo.Elements.keys()):
            if i in list(dynamic.ProcessConstraintData.Elements.keys()):
                  static.ProcessConstraintsInfo.Elements[i].update({'DataElement':'ProcessConstraintData##'+i+'_Data'})
                  dynamic.ProcessConstraintData.Elements[i].update({'StructureElement':'ProcessConstraint##'+i})
      for i in list(static.DensityConstraintsInfo.Elements.keys()):
            if i in list(dynamic.DensityConstraintData.Elements.keys()):
                  static.DensityConstraintsInfo.Elements[i].update({'DataElement':'DensityConstraintData##'+i+'_Data'})
                  dynamic.DensityConstraintData.Elements[i].update({'StructureElement':'DensityConstraint##'+i})
      for i in list(static.EnzymeConstraintsInfo.Elements.keys()):
            if i in list(dynamic.EnzymeConstraintData.Elements.keys()):
                  static.EnzymeConstraintsInfo.Elements[i].update({'DataElement':'EnzymeConstraintData##'+i+'_Data'})
                  dynamic.EnzymeConstraintData.Elements[i].update({'StructureElement':'EnzymeConstraint##'+i})
      for i in list(static.ReactionInfo.Elements.keys()):
            if i in list(dynamic.ReactionData.Elements.keys()):
                 static.ReactionInfo.Elements[i].update({'DataElement':'ReactionData##'+i+'_Data'})
                 dynamic.ReactionData.Elements[i].update({'StructureElement':'Reaction##'+i})
      for i in list(static.EnzymeInfo.Elements.keys()):
            if i in list(dynamic.EnzymeData.Elements.keys()):
                 static.EnzymeInfo.Elements[i].update({'DataElement':'EnzymeData##'+i+'_Data'})
                 dynamic.EnzymeData.Elements[i].update({'StructureElement':'Enzyme##'+i})
      for i in list(static.ProcessInfo.Elements.keys()):
            if i in list(dynamic.ProcessData.Elements.keys()):
                 static.ProcessInfo.Elements[i].update({'DataElement':'ProcessData##'+i+'_Data'})
                 dynamic.ProcessData.Elements[i].update({'StructureElement':'Process##'+i})
      for i in list(static.ProteinInfo.Elements.keys()):
            if i in list(dynamic.ProteinData.Elements.keys()):
                 static.ProteinInfo.Elements[i].update({'DataElement':'ProteinData##'+i+'_Data'})
                 dynamic.ProteinData.Elements[i].update({'StructureElement':'Protein##'+i})
      return({'Structural': static , 'SimulationData': dynamic})
