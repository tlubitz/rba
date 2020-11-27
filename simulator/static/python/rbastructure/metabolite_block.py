# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import numpy
import json
import urllib3
import pandas
import libsbml

# package imports
from ..rba.core.constraint_blocks import ConstraintBlocks
from .element_block import ElementBlock


class MetaboliteBlock(ElementBlock):
    """
    Class holding information on the metabolites in the model.

  Attributes
  ----------
  Elements : dict
      Each model-meatbolite is represented by a key.
      The values, holding information on each metabolite, are dicts with predefined keys:
          'ID' : meatbolite ID in model (type str)
          'OtherIDs' : identifiers of this metabolite in other namespaces (BiGG, KEGG ...) (type dict)
          'Name' : Name according to BiGG (type str)
          'ReactionsInvolvedWith' : Reactions which produce or consume this metabolite (type list)
          'boundary' : Boundary metabolite (type boolean)
          'Type' :  Type of metabolite (internal exernal or biomass-precursor) (type str)
          'Compartment' : Location of meatbolite (type str)
    """

    def fromFiles(self, model, Info, MetaboliteAnnotations, sbml):
        """
        Derive reaction-info from RBA-model.

        Input
        -----
        RBA-model
        Dataframe, holding BiGG-metabolite information.

        Returns
        -------
        Dictionary with metabolite-info.

        """
        blocks = ConstraintBlocks(model)
        targetMetabolites = findTargetMetabolites(model)
        self.Elements = {}
        index = 0
        http = urllib3.PoolManager()
        reconstruction = Info.ix['Reconstruction', 1]
        if type(sbml) is not str:
            if type(sbml.model) is libsbml.Model:
                sbmlIDMap = [species.id for species in sbml.model.species]
        for m in range(len(model.metabolism.species._elements)):
            i = model.metabolism.species._elements[m].id
            index += 1
            IDdict = {}
            speciesname = ' '
            speciescompartment = i[-1]
            if type(sbml) is not str:
                if type(sbml.model) is libsbml.Model:
                    if i in sbmlIDMap:
                        IDdict.update(getMetaboliteAnnotationsFromSBML(sbmlIDMap.index(i), sbml))
                        speciesname = sbml.model.species[sbmlIDMap.index(i)].name
                        speciescompartment = sbml.model.species[sbmlIDMap.index(i)].compartment
            if type(MetaboliteAnnotations) is pandas.core.frame.DataFrame:
                IDdict.update(readMetaboliteAnnotations(i, MetaboliteAnnotations))
            if i in blocks.metabolism.external:
                self.Elements[i] = {'ID': i,
                                    #                                      'Name': findMetInfo(i[2:],metabolitesBiGG) ,
                                    'Name': speciesname,
                                    'OtherIDs': IDdict,
                                    'index': index,
                                    'ReactionsInvolvedWith': associatedReactions(i, blocks),
                                    'boundary': model.metabolism.species._elements[m].boundary_condition,
                                    'Type': 'external',
                                    'Compartment': speciescompartment}
            elif i in blocks.metabolism.internal:
                typ = 'internal'
                if checkForTarget(i, targetMetabolites):
                    typ = 'precursor'
                self.Elements[i] = {'ID': i,
                                    'Name': speciesname,
                                    'OtherIDs': IDdict,
                                    'ReactionsInvolvedWith': associatedReactions(i, blocks),
                                    'index': index,
                                    'boundary': model.metabolism.species._elements[m].boundary_condition,
                                    'Type': typ,
                                    'Compartment': speciescompartment}

    def overview(self):
        """
        Derive statistics on metabolites.

        Returns
        -------
        Dictionary with general numbers on metabolites.

        """
        nT = len(self.Elements.keys())
        nI = 0
        nE = 0
        nBio = 0
        nBound = 0
        for i in self.Elements.keys():
            if self.Elements[i]['Type'] == 'internal':
                nI += 1
            if self.Elements[i]['Type'] == 'external':
                nE += 1
            if self.Elements[i]['Type'] == 'precursor':
                nBio += 1
            if self.Elements[i]['boundary']:
                nBound += 1
        out = {'nMetabolitesTotal': nT,
               'nMetabolitesInternal': nI,
               'nMetabolitesExternal': nE,
               'nMetabolitesGrowthRelevant': nBio,
               'nBoundaryMetabolites': nBound}
        return(out)


def associatedReactions(metabolite, blocks):
    out = []
    if metabolite in blocks.metabolism.internal:
        Srow = blocks.metabolism.S.toarray()[blocks.metabolism.internal.index(metabolite), :]
        out = list(numpy.asarray(blocks.metabolism.reactions)[numpy.nonzero(Srow)])
    return(out)


def getMetaboliteAnnotationsFromSBML(index, sbml):
    out = {}
    for a in sbml.model.species[index].getAnnotationString().split('\n'):
        if 'rdf:resource="http://identifiers.org/' in a:
            annotationType = a.split(
                'rdf:resource="http://identifiers.org/')[1].split('"/>')[0].split('/')
            out.update({annotationType[0]: annotationType[1]})
    return(out)


def readMetaboliteAnnotations(m, MetaboliteAnnotations):
    AnnotationKeys = list(MetaboliteAnnotations)
    AnnotationIDs = [numpy.nan]*len(AnnotationKeys)
    for i in AnnotationKeys:
        if m in list(MetaboliteAnnotations[i]):
            metabolite = list(MetaboliteAnnotations[i]).index(m)
            AnnotationIDs = list(MetaboliteAnnotations.ix[metabolite])
    return(dict(zip(AnnotationKeys, AnnotationIDs)))


def findOtherIDs(met_name, http, reconstruction):
    #     if met_name[-2] is '_': met_name=met_name.rsplit('_',1)[0]
    if not reconstruction is 'GSMM':
        response = http.request('GET', 'http://bigg.ucsd.edu/api/v2/models/' +
                                reconstruction + '/metabolites/' + met_name)
    else:
        response = http.request('GET', 'http://bigg.ucsd.edu/api/v2/models/' +
                                'universal' + '/metabolites/' + met_name)
    try:
        x = json.loads(response.data.decode('utf-8'))
        out = {'BiGG': 'NA', 'KEGG': 'NA', 'CHEBI': 'NA', 'BioCyc': 'NA', 'SEED': 'NA', 'Name': ''}
        out['BiGG'] = met_name
        if 'KEGG Compound' in list(x['database_links'].keys()):
            out['KEGG'] = [d['id'] for d in x['database_links']['KEGG Compound']]
        if 'CHEBI' in list(x['database_links'].keys()):
            out['CHEBI'] = [d['id'] for d in x['database_links']['CHEBI']]
        if 'BioCyc' in list(x['database_links'].keys()):
            out['BioCyc'] = [d['id'] for d in x['database_links']['BioCyc']]
        if 'BioCyc' in list(x['database_links'].keys()):
            out['SEED'] = [d['id'] for d in x['database_links']['SEED Compound']]
        out['Name'] = x['name']
    except:
        out = {'BiGG': 'NA', 'KEGG': 'NA', 'CHEBI': 'NA', 'BioCyc': 'NA', 'SEED': 'NA', 'Name': ''}
#     print(out)
    return(out)


def findMetInfo(metaboliteID, metabolitesBiGG):
    if metaboliteID in list(metabolitesBiGG.index):
        out = str(metabolitesBiGG.ix[metaboliteID]['name'])
    else:
        out = str('')
    return(out)


def findTargetMetabolites(model):
    out = []
    for j in range(len(model.targets.target_groups._elements)):
        if model.targets.target_groups._elements[j].id == 'metabolite_production':
            for k in range(len(model.targets.target_groups._elements[j].concentrations._elements)):
                out.append(
                    model.targets.target_groups._elements[j].concentrations._elements[k].species)
    return out


def checkForTarget(i, targetMets):
    if i in targetMets:
        return(True)
    else:
        return(False)
