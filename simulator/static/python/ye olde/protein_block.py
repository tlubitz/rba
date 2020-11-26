# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import numpy
# package imports
from rbastructure.element_block import ElementBlock


class ProteinBlock(ElementBlock):
    """
    Class holding information on the proteins in the model.

   Attributes
   ----------
   Elements : dict
       Each model-protein is represented by a key.
       The values, holding information on each protein, are dicts with predefined keys:
           'ID' : protein ID in model (type str)
           'ExternalIDs' : identifiers of this protein in other namespaces (Locus-tag, Gene-symbol, ...) (type dict)
           'Function' : associated function, according to Uniprot (type str)
           'Name' : name, according to Uniprot (type list)
           'associatedEnzymes' : Enzymes, this protein is a subunit of (type list)
           'Compartment' : Location of the protein (type str)
           'AAcomposition' : Which amino-acids the protein is composed of and what is the stoichiometry(type dict)
           'AAnumber' : Length of protein (type int)
           'Weight' : Weight of protein (type float)
           'ProcessRequirements' : Which processes the protein requires for synthesis and maintenance and how much (type dict)
           'SupportsProcess' : Process-machineries, this protein is a subunit of (type list)
    """

    def fromFiles(self, model, UniprotData=None, IDmap, Info):
        """
        Derive reaction-info from RBA-model.

        Input
        -----
        RBA-model
        Dataframe, holding uniprot-file information.
        Dataframe, holding ID-map information.

        Returns
        -------
        Dictionary with protein-info.

        """
        Prots = getProteinList(model)
        if IDmap is 'Not There':
            IDmap = pandas.DataFrame(index=Prots)
        self.Elements = {}
        index = 0
        for i in range(len(Prots)):
            RequiredProcess = ProcessRequirements(model, Prots[i])
            Compostion = getProteinComposition(model, i)
            genes = {'ids': {},
                     'desc': '',
                     'name': '',
                     'weight': '',
                     'length': ''}
            if UniprotData is not None:
                genes = mineUniprotFile(UniprotData, IDmap, Prots[i])
            index += 1
            self.Elements[Prots[i]] = {'ID': Prots[i],
                                       'ExternalIDs': genes['ids'],
                                       'Function': genes['desc'],
                                       'Name': genes['name'],
                                       'associatedReactions': [],
                                       'index': index,
                                       'associatedEnzymes': [],
                                       'Compartment': getCompartment(model, i),
                                       'AAcomposition': Compostion['AAcomp'],
                                       'AAnumber': genes['length'],
                                       'Weight': genes['weight'],
                                       'ProcessRequirements': ProcessCost(model, Compostion['AAcomp'], RequiredProcess),
                                       'SupportsProcess':  RequiredProcess['Part']}

    def overview(self):
        """
        Derive statistics on proteins.

        Returns
        -------
        Dictionary with general numbers on proteins.

        """
        nT = len(self.Elements.keys())
        out = {'nProteinsTotal': nT}
        return(out)


def getProteinList(model):
    out = []
    for e in model.proteins.macromolecules._elements:
        out.append(e.id)
    return(out)


def mineUniprotFile(UniprotData, IDmap, protein):
    differentIDs = {}
    function = ''
    name = ''
    length = numpy.nan
    mass = numpy.nan
    IDlist = [' ' + m + ' ' + str(n) + ' ' for m,
              n in zip(list(UniprotData['Entry']), list(UniprotData['Gene names']))]
    ProteinRowList = [i for i, s in enumerate(IDlist) if str(' ' + protein + ' ') in s]
    if len(ProteinRowList) > 0:
        ProteinRow = ProteinRowList[0]
        uniprotID = UniprotData['Entry'][ProteinRow]
        differentIDs.update({'UniprotID': uniprotID})
        # and len(str(UniprotData['Length'][ProteinRow])) >= 1:
        if str(UniprotData['Length'][ProteinRow]) is not 'nan':
            length = UniprotData['Length'][ProteinRow]
        # and len(str(UniprotData['Mass'][ProteinRow])) >= 1:
        if str(UniprotData['Mass'][ProteinRow]) is not 'nan':
            mass = UniprotData['Mass'][ProteinRow]
        # and len(str(UniprotData['EC number'][ProteinRow])) >= 1:
        if str(UniprotData['EC number'][ProteinRow]) is not 'nan':
            differentIDs.update({'ECnumber': 'EC '+str(UniprotData['EC number'][ProteinRow])})
        # and type(UniprotData['Function [CC]'][ProteinRow]) is 'str':
        if str(UniprotData['Function [CC]'][ProteinRow]) is not 'nan':
            function = UniprotData['Function [CC]'][ProteinRow]
        # and len(str(UniprotData['Protein names'][ProteinRow])) >= 1:
        if str(UniprotData['Protein names'][ProteinRow]) is not 'nan':
            name = UniprotData['Protein names'][ProteinRow]
        if type(IDmap) is not str:
            for j in range(IDmap.shape[1]):
                if '##()##' not in list(IDmap)[j]:
                    differentIDs[list(IDmap)[j]] = IDmap.ix[uniprotID, j]
                if '##()##' in list(IDmap)[j]:
                    differentIDs[list(IDmap)[j].split('##()##')[1]] = list(
                        IDmap)[j].split('##()##')[0] + '###' + IDmap.ix[uniprotID, j]

    out = {'ids': differentIDs,
           'desc': function,
           'name': name,
           'weight': mass,
           'length': length}

    return(out)


def getCompartment(model, protein):
    return(model.proteins.macromolecules._elements[protein].__dict__['compartment'])


def getProteinComposition(model, protein):
    out = {}
    numberAA = 0
    MacroMolecules = model.proteins.macromolecules._elements
    for a in range(len(MacroMolecules[protein].composition._elements)):
        out[MacroMolecules[protein].composition._elements[a].component] = int(
            round(MacroMolecules[protein].composition._elements[a].stoichiometry, 3))  # round(...,3) added#
        numberAA += MacroMolecules[protein].composition._elements[a].stoichiometry
    out = {'AAcomp': out,
           'AAnum': int(numberAA)}
    return(out)


def ProcessRequirements(model, protein):
    out1 = []
    out2 = []
    Processes = model.processes.processes._elements
    for p in range(len(Processes)):
        for p2 in range(len(Processes[p].processings.productions._elements)):
            if Processes[p].processings.productions._elements[p2].set == 'protein':
                for inp in range(len(Processes[p].processings.productions._elements[p2].inputs)):
                    if Processes[p].processings.productions._elements[p2].inputs._elements[inp].__dict__['species'] == protein:
                        out1.append(Processes[p].name)
        for p3 in range(len(Processes[p].machinery.machinery_composition.reactants._elements)):
            if Processes[p].machinery.machinery_composition.reactants._elements[p3].species == protein:
                out2.append(Processes[p].name)
    out = {'Req': out1,
           'Part': out2}
    return(out)


def ProcessCost(model, AminoAcid, req):
    out = {}
    Processes = model.processes.processes._elements
    ProcessingMaps = model.processes.processing_maps._elements
    for Pref in req['Req']:
        cost = 0
        for p in range(len(Processes)):
            if Pref == Processes[p].name:
                Pmap = Processes[p].processings.productions._elements[0].processing_map
        for pm in range(len(ProcessingMaps)):
            AAcost = 0
            if ProcessingMaps[pm].id == Pmap:
                for a in AminoAcid.keys():
                    for ap in range(len(ProcessingMaps[pm].component_processings._elements)):
                        if ProcessingMaps[pm].component_processings._elements[ap].component == a:
                            AAcost = AminoAcid[a] * \
                                ProcessingMaps[pm].component_processings._elements[ap].machinery_cost
                            cost += AAcost
            out[Pref] = round(cost, 3)  # round(...,3) added#
    return(out)


def ProcessCost_Old(model, AminoAcid, req):
    out = {}
    Processes = model.processes.processes._elements
    ProcessingMaps = model.processes.processing_maps._elements
    for Pref in req['Req']:
        cost = 0
        for p in range(len(Processes)):
            if Pref == Processes[p].name:
                Pmap = Processes[p].processings.productions._elements[0].processing_map
        for pm in range(len(ProcessingMaps)):
            AAcost = 0
            if ProcessingMaps[pm].id == Pmap:
                for a in AminoAcid.keys():
                    for ap in range(len(ProcessingMaps[pm].component_processings._elements)):
                        if ProcessingMaps[pm].component_processings._elements[ap].component == a:
                            AAcost = AminoAcid[a] * \
                                ProcessingMaps[pm].component_processings._elements[ap].machinery_cost
            cost += AAcost
            out[Pref] = round(cost, 3)  # round(...,3) added#
    return(out)
