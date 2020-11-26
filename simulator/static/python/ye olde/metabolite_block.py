# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import numpy
import json
import urllib3
# package imports
from rba.core.constraint_blocks import ConstraintBlocks
from rbastructure.element_block import ElementBlock


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

    def fromFiles(self, model, Info):
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
        for m in range(len(model.metabolism.species._elements)):
            i = model.metabolism.species._elements[m].id
            index += 1
            IDdict = {}
            KEGGid = findOtherIDs(i[2:], http, reconstruction)
            if KEGGid['BiGG'] is not 'NA':
                IDdict.update({'BiGG.Metabolite': KEGGid['BiGG']})
            if KEGGid['KEGG'] is not 'NA':
                IDdict.update({'KEGG.Metabolite': KEGGid['KEGG']})
            if KEGGid['CHEBI'] is not 'NA':
                IDdict.update({'CHEBI.Metabolite': KEGGid['CHEBI']})
            if KEGGid['BioCyc'] is not 'NA':
                IDdict.update({'BioCyc.Metabolite': KEGGid['BioCyc']})
            if KEGGid['SEED'] is not 'NA':
                IDdict.update({'SEED.Metabolite': KEGGid['SEED']})
            if i in blocks.metabolism.external:
                self.Elements[i] = {'ID': i,
                                    #                                      'Name': findMetInfo(i[2:],metabolitesBiGG) ,
                                    'Name': KEGGid['Name'],
                                    'OtherIDs': IDdict,
                                    'index': index,
                                    'ReactionsInvolvedWith': associatedReactions(i, blocks),
                                    'boundary': model.metabolism.species._elements[m].boundary_condition,
                                    'Type': 'external',
                                    'Compartment': i[-1]}
            elif i in blocks.metabolism.internal:
                typ = 'internal'
                if checkForTarget(i, targetMetabolites):
                    typ = 'precursor'
                self.Elements[i] = {'ID': i,
                                    'Name': KEGGid['Name'],
                                    'OtherIDs': IDdict,
                                    'ReactionsInvolvedWith': associatedReactions(i, blocks),
                                    'index': index,
                                    'boundary': model.metabolism.species._elements[m].boundary_condition,
                                    'Type': typ,
                                    'Compartment': i[-1]}

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
