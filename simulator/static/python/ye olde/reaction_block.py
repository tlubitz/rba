# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import numpy
import json
import urllib3
# package imports
from rba.core.constraint_blocks import ConstraintBlocks
from rbastructure.element_block import ElementBlock


class ReactionBlock(ElementBlock):
    """
    Class holding information on the reactions in the model.

   Attributes
   ----------
   Elements : dict
       Each model-enzyme is represented by a key.
       The values, holding information on each enzyme, are dicts with predefined keys:
           'ID' : reaction ID in model (type str)
           'Compartment_Machinery' : Localisation of enzyme-subunit, catalysing this reaction. (type list)
           'Name' : name according to BiGG (type str)
           'OtherIDs' : Other reaction names (eg. BiGG, KEGG) (type list)
           'Formula' : Reaction formula as string (type str)
           'Reactants' : Which metabolites does this reaction consume of and many (type dict)
           'Products' : Which metabolites does this reaction produce of and many (type dict)
           'Type' : Type of reaction ('normal' or 'transport') (type str)
           'Compartment_Species' : Location of the metabolites involved with this reaction (type list)
           'Enzyme' : Enzyme catalysing this reaction (type str)
           'Twins' : Isoreactions of this reactions (catalysed by iso-enzymes) (type list)
    """

    def fromFiles(self, model, Info):
        """
        Derive reaction-info from RBA-model.

        Input
        -----
        RBA-model
        Dataframe, holding BiGG-reaction information.

        Returns
        -------
        Dictionary with reaction-info.

        """
        blocks = ConstraintBlocks(model)
        self.Elements = {}
        self.BiGGids = []
        index = 0
        http = urllib3.PoolManager()
        reconstruction = Info.ix['Reconstruction', 1]
        for i in range(len(blocks.metabolism.reactions)):
            Reactants = associatedReactants(i, blocks)
            Products = associatedProducts(i, blocks)
            Reversibility = checkReversibility(i, blocks)
            CompartmentInfo = findCompartment(i, blocks, Reactants, Products, Reversibility)
            Twins = findTwinRxns(i, blocks)
            biggID = deriveProto_ID(i, blocks, Twins)
            IDdict = {'ProtoID': biggID}
            KEGGid = findKEGGID(biggID, http, reconstruction)
            if KEGGid['BiGG'] is not 'NA':
                IDdict.update({'BiGG.Reaction': KEGGid['BiGG']})
            if KEGGid['KEGG'] is not 'NA':
                IDdict.update({'KEGG.Reaction': KEGGid['KEGG']})
            if KEGGid['BioCyc'] is not 'NA':
                IDdict.update({'BioCyc.Reaction': KEGGid['BioCyc']})
            self.BiGGids.append(biggID)
            index += 1
            self.Elements[blocks.metabolism.reactions[i]] = {'ID': blocks.metabolism.reactions[i],
                                                             'Compartment_Machinery': [],
                                                             'Name': KEGGid['Name'],
                                                             'OtherIDs': IDdict,
                                                             'index': index,
                                                             'Formula': Reactants['rSide'] + ' <=> ' + Products['pSide'],
                                                             'Reactants': Reactants['reactants'],
                                                             'Products': Products['products'],
                                                             'Reversible': Reversibility['Reversible'],
                                                             'Type': CompartmentInfo['type'],
                                                             'Compartment_Species': CompartmentInfo['comp'],
                                                             'Enzyme': findAssociatedEnzyme(i, blocks),
                                                             'Twins': Twins}

    def overview(self):
        """
        Derive statistics on reactions.

        Returns
        -------
        Dictionary with general numbers on reactions.

        """
        nTot = len(list(self.Elements.keys()))
        nUnique = numpy.nan
        nSpont = 0
        nEnzyme = 0
        nInternalTransport = 0
        nExchange = 0
        nInternal = 0
        nRev = 0
        nIrrev = 0
        BiGGIDs = []
        for i in list(self.Elements.keys()):
            if len(self.Elements[i]['Enzyme']) > 0:
                nEnzyme += 1
            else:
                nSpont += 1
            if self.Elements[i]['Type'] == 'Normal':
                nInternal += 1
            if self.Elements[i]['Type'] == 'Transport':
                if len(self.Elements[i]['Reactants']) == 0 or len(self.Elements[i]['Products']) == 0:
                    nExchange += 1
                else:
                    nInternalTransport += 1
            if self.Elements[i]['Reversible']:
                nRev += 1
            else:
                nIrrev += 1
#               BiGGIDs.append(self.Elements[i]['OtherIDs']['BiGG.Reaction'])
        nUnique = len(list(numpy.unique(self.BiGGids)))
        out = {'nReactionsTotal': nTot,
               'nReactionsUnique': nUnique,
               'nReactionsSpontaneous': nSpont,
               'nReactionsEnzymatic': nEnzyme,
               'nReactionsInternal': nInternal,
               'nReactionsExchange': nExchange,
               'nReactionsCompartmentTransport': nInternalTransport,
               'nReactionsReversible': nRev,
               'nReactionsIrreversible': nIrrev}
        return(out)


def deriveProto_ID(reaction, blocks, Twins):
    """
    Derive BiGG-ID from reactionID.
    Relies on the assumption that reactionID in RBA-model equals an 'R_', followed by the BiGG-ID.
    When the reaction has twins (duplicated reactions due to isozymes),
    the one without an appended '_2', '_3' ... needs to be found first.

    Returns
    -------
    String with derived BiGG-ID.
    """
    if len(Twins) > 0:  # Check if isozymic reactions exist
        # Find "original" reaction (shortest name)
        b = min(Twins+[blocks.metabolism.reactions[reaction]], key=len)
        out = b  # Remove 'M_'-prefix
        if b.startswith('R_'):
            out = b[2:]
    else:
        out = blocks.metabolism.reactions[reaction]
        if out.startswith('R_'):
            out = out[2:]
    return(out)


def findKEGGID(rx_name, http, reconstruction):
    #     if not reconstruction is 'GSMM':
    #        response = http.request('GET', 'http://bigg.ucsd.edu/api/v2/models/' + reconstruction +'/reactions/' + rx_name)
    if not reconstruction is 'GSMM':
        response = http.request('GET', 'http://bigg.ucsd.edu/api/v2/models/' +
                                reconstruction + '/reactions/' + rx_name)
    else:
        response = http.request('GET', 'http://bigg.ucsd.edu/api/v2/models/' +
                                'universal' + '/reactions/' + rx_name)
    try:
        x = json.loads(response.data.decode('utf-8'))
        out = {'BiGG': 'NA', 'KEGG': 'NA', 'BioCyc': 'NA', 'Name': ''}
        out['BiGG'] = rx_name
        out['Name'] = x['name']
        if 'KEGG Reaction' in list(x['database_links'].keys()):
            out['KEGG'] = x['database_links']['KEGG Reaction'][0]['id']
            out['BioCyc'] = x['database_links']['BioCyc'][0]['id']
        return(out)
    except:
        return({'BiGG': 'NA', 'KEGG': 'NA', 'BioCyc': 'NA', 'Name': ''})
#     else:
#        return('NA')


def findRxnName(rxnid, reactionsBiGG):
    """
    Retreive (descriptive) name of reaction from BiGG-file.

    Returns
    -------
    String with name (when found), otherwise empty.
    """
    if rxnid in list(reactionsBiGG.index):
        return(str(reactionsBiGG.ix[rxnid]['name']))
    else:
        return(str(''))


def associatedReactants(i, blocks):
    """
    Derive information of reactant-side of reaction.

    Returns
    -------
    'reactants': Dictionary with reactants and stoichiometric factors.
    'rSide': String, representing reactant-side of reaction-formula.
    """
    Scol = blocks.metabolism.S.toarray()[:, i]
    rS = list(numpy.asarray(blocks.metabolism.internal)[Scol < 0])
    sF = [int(i) for i in list(abs(Scol[Scol < 0]))]
    eq = ""
    for k in range(len(rS)):
        if k > 0:
            eq = eq+' + '
            eq = eq+'{} '.format(sF[k])
        eq = eq+rS[k]
    out = {'reactants': dict(zip(rS, sF)),
           'rSide': eq}
    return(out)


def associatedProducts(i, blocks):
    """
    Derive information of product-side of reaction.

    Returns
    -------
    'products': Dictionary with products and stoichiometric factors.
    'pSide': String, representing product-side of reaction-formula.
    """
    Scol = blocks.metabolism.S.toarray()[:, i]
    pS = list(numpy.asarray(blocks.metabolism.internal)[Scol > 0])
    sF = [int(i) for i in list(abs(Scol[Scol > 0]))]
    eq = ""
    for k in range(len(pS)):
        if k > 0:
            eq = eq+' + '
            eq = eq+'{} '.format(sF[k])
        eq = eq+pS[k]
    out = {'products': dict(zip(pS, sF)),
           'pSide': eq}
    return(out)


def checkReversibility(rx, blocks):
    """
    Information on default reaction flux-bounds and reversibility.

    Returns
    -------
    Dictionary with numerical values on flux-bounds and boolean for reversibility.
    """
    LB = blocks.metabolism._lb[rx].__dict__['value']
    UB = blocks.metabolism._ub[rx].__dict__['value']
    out = {'Reversible': True,
           'UB': UB,
           'LB': LB}
    if LB == 0:
        out['Reversible'] = False
    return(out)


def findCompartment(rx, blocks, aR, aP, rR):
    """
    Derive information on compartment aspects of the reaction.

    Returns
    -------
    'type': String 'Transport' or 'Normal' (kind of reaction).
    'comp': compartment of SBML-model. arrow between compartments when reaction is 'Transport'.
    """
    r = list(aR['reactants'].keys())
    if len(r) > 0:
        rComp = r[0]
    elif len(r) == 0:
        rComp = ' '
    p = list(aP['products'].keys())
    if len(p) > 0:
        pComp = p[0]
    else:
        pComp = ' '
    if rComp[-1] == pComp[-1]:
        typ = 'Normal'
        comp = rComp[-1]
    else:
        typ = 'Transport'
        comp = rComp[-1] + '-->' + pComp[-1]
        if rR['Reversible'] == 'True':
            comp = rComp[-1] + '<==>' + pComp[-1]
    out = {'type': typ,
           'comp': comp}
    return(out)


def findAssociatedEnzyme(rx, blocks):
    """
    Return enzyme species, associated with reaction.

    Returns
    -------
    String with enzymeID, empty string if reaction is spontaneous.
    """
    RN = blocks.metabolism.reactions[rx]
    if RN in blocks.enzymes.reaction_catalyzed:
        return(blocks.enzymes.ids[blocks.enzymes.reaction_catalyzed.index(RN)])
    else:
        return(str(''))


def findTwinRxns(rx, blocks):
    """
    Find Twin reactions (identical reactions, catalyzed by different (iso-)enzymes)

    Returns
    -------
    List of iso-reactions.
    """
    out = []
    if 'duplicate' in blocks.metabolism.reactions[rx]:
        for x in blocks.metabolism.reactions:
            if blocks.metabolism.reactions[rx].rsplit('_duplicate')[0] in x:
                if blocks.metabolism.reactions[rx] is not x:
                    out.append(x)
    else:
        for x in blocks.metabolism.reactions:
            if blocks.metabolism.reactions[rx]+'_duplicate' in x:
                out.append(x)
    return(out)
