import rba
import numpy
import pandas
from rbatools.rba_Session import RBA_Session
from scipy.stats.mstats import gmean


class RBA_Calibrator(object):
    def __init__(self, xml_dir):
        self.rbaSession = RBA_Session(xml_dir)

    def estimate_specific_Kapps(self, proteomicsData, flux_bounds, mu, biomass_function=None, target_biomass_function=True):
        """
        Parameters
        ----------
        proteomicsData : pandas.DataFrame (in mmol/gDW)
        flux_bounds : pandas.DataFrame  (in mmol/(gDW*h))
        mu : float (in 1/h)
        biomass_function : str
        target_biomass_function : bool
        atp_maintenance_to_biomassfunction : bool
        eukaryotic : bool
        """
        Avogadro_constant = 6.022e23

        self.rbaSession.addExchangeReactions()
        self.rbaSession.setMu(mu)

        if target_biomass_function:
            self.rbaSession.buildFBA(objective='targets', maintenanceToBM=True)
            BMfunction = 'R_BIOMASS_targetsRBA'
        else:
            self.rbaSession.buildFBA(objective='classic', maintenanceToBM=False)
            BMfunction = biomass_function

        for j in [i for i in self.rbaSession.Medium.keys() if self.rbaSession.Medium[i] == 0]:
            Exrxn = 'R_EX_'+j.split('M_')[-1]+'_e'
            self.rbaSession.FBA.setUB({Exrxn: 0})

        rxn_LBs = {}
        rxn_UBs = {}
        for rx in flux_bounds['Reaction_ID']:
            lb = flux_bounds.loc[flux_bounds['Reaction_ID'] == rx, 'LB'].values[0]
            ub = flux_bounds.loc[flux_bounds['Reaction_ID'] == rx, 'UB'].values[0]
            if not pandas.isna(lb):
                rxn_LBs.update({rx: lb})
            if not pandas.isna(ub):
                rxn_UBs.update({rx: ub})
        self.rbaSession.FBA.setLB(rxn_LBs)
        self.rbaSession.FBA.setUB(rxn_UBs)

        self.rbaSession.FBA.clearObjective()
        self.rbaSession.FBA.setObjectiveCoefficients({BMfunction: -1})
        self.rbaSession.FBA.solveLP(feasibleStatuses=[1, 2, 3, 5, 6])
        BMfluxOld = self.rbaSession.FBA.SolutionValues[BMfunction]

        self.rbaSession.FBA.parsimonise()
        self.rbaSession.FBA.setLB(rxn_LBs)
        self.rbaSession.FBA.setUB(rxn_UBs)
        self.rbaSession.FBA.setLB({BMfunction: BMfluxOld})
        self.rbaSession.FBA.setUB({BMfunction: BMfluxOld})
        self.rbaSession.FBA.solveLP(feasibleStatuses=[1, 2, 3, 5, 6])

        FluxDistribution = pandas.DataFrame(index=list(
            self.rbaSession.FBA.SolutionValues.keys()), columns=['FluxValues'])
        FluxDistribution['FluxValues'] = list(self.rbaSession.FBA.SolutionValues.values())
        BMfluxNew = self.rbaSession.FBA.SolutionValues[BMfunction]

        ProtoIDmap = {}
        for i in self.rbaSession.ModelStructure.ProteinInfo.Elements.keys():
            ProtoID = self.rbaSession.ModelStructure.ProteinInfo.Elements[i]['ProtoID']
            if ProtoID in list(proteomicsData.index):
                if not pandas.isna(proteomicsData.loc[ProtoID, 'copy_number']):
                    if proteomicsData.loc[ProtoID, 'copy_number'] != 0:
                        if ProtoID in ProtoIDmap.keys():
                            ProtoIDmap[ProtoID]['ModelProteins'].append(i)
                        else:
                            ProtoIDmap.update(
                                {ProtoID: {'ModelProteins': [i], 'CopyNumber': proteomicsData.loc[ProtoID, 'copy_number']}})

        ReactionMap = {}
        for i in self.rbaSession.ModelStructure.ReactionInfo.Elements.keys():
            if '_duplicate_' in i:
                continue
            else:
                if i in list(FluxDistribution.index):
                    if FluxDistribution.loc[i, 'FluxValues'] != 0:
                        ReactionMap.update({i: {'ModelReactions': list(
                            [i]+self.rbaSession.ModelStructure.ReactionInfo.Elements[i]['Twins']), 'Flux': FluxDistribution.loc[i, 'FluxValues']}})

        IsoReaction2ProtoReaction = {}
        for i in ReactionMap.keys():
            for j in ReactionMap[i]['ModelReactions']:
                IsoReaction2ProtoReaction[j] = i

        EnzymeMap = {}
        for i in self.rbaSession.ModelStructure.EnzymeInfo.Elements.keys():
            if self.rbaSession.ModelStructure.EnzymeInfo.Elements[i]['Reaction'] in IsoReaction2ProtoReaction:
                CompositionDict = {self.rbaSession.ModelStructure.ProteinInfo.Elements[j]['ProtoID']: self.rbaSession.ModelStructure.EnzymeInfo.Elements[
                    i]['Subunits'][j] for j in self.rbaSession.ModelStructure.EnzymeInfo.Elements[i]['Subunits'].keys()}
                ProtoReaction = IsoReaction2ProtoReaction[self.rbaSession.ModelStructure.EnzymeInfo.Elements[i]['Reaction']]
                CopyNumbers = []
                Stoichiometries = []
                EnzymeNumbers = []
                for j in CompositionDict.keys():
                    if j in ProtoIDmap.keys():
                        CopyNumbers.append(ProtoIDmap[j]['CopyNumber'])
                        Stoichiometries.append(CompositionDict[j])
                        EnzymeNumbers.append(ProtoIDmap[j]['CopyNumber']/CompositionDict[j])
                GM_enzymenumber = 0
                if len(EnzymeNumbers) > 0:
                    GM_enzymenumber = gmean(numpy.array(EnzymeNumbers))
                EnzymeMap.update(
                    {i: {'ProtoReaction': ProtoReaction, 'EnzymeNumber': GM_enzymenumber}})

        EnzymeMap2 = {}
        for i in ReactionMap.keys():
            totalIsoEnzymeNumber = 0
            for j in ReactionMap[i]['ModelReactions']:
                respectiveEnzyme = self.rbaSession.ModelStructure.ReactionInfo.Elements[j]['Enzyme']
                if respectiveEnzyme in EnzymeMap.keys():
                    totalIsoEnzymeNumber += EnzymeMap[respectiveEnzyme]['EnzymeNumber']
            for j in ReactionMap[i]['ModelReactions']:
                respectiveEnzyme = self.rbaSession.ModelStructure.ReactionInfo.Elements[j]['Enzyme']
                if respectiveEnzyme in EnzymeMap.keys():
                    if EnzymeMap[respectiveEnzyme]['EnzymeNumber'] != 0:
                        specificFlux = ReactionMap[i]['Flux'] * \
                            EnzymeMap[respectiveEnzyme]['EnzymeNumber']/totalIsoEnzymeNumber
                        concentration = EnzymeMap[respectiveEnzyme]['EnzymeNumber'] / \
                            Avogadro_constant
                        EnzymeMap2.update({respectiveEnzyme: {'CopyNumber': EnzymeMap[respectiveEnzyme]['EnzymeNumber'],
                                                              'Concentration': concentration, 'Flux': specificFlux, 'Kapp': abs(specificFlux/concentration)}})

        self.specific_Kapps = pandas.DataFrame()
        for i in EnzymeMap2.keys():
            if EnzymeMap2[i]['CopyNumber'] == 0:
                continue
            self.specific_Kapps.loc[i, 'Enzyme_ID'] = i
            self.specific_Kapps.loc[i, 'CopyNumber'] = EnzymeMap2[i]['CopyNumber']
            self.specific_Kapps.loc[i, 'Concentration'] = EnzymeMap2[i]['Concentration']
            self.specific_Kapps.loc[i, 'Flux'] = EnzymeMap2[i]['Flux']
            self.specific_Kapps.loc[i, 'Kapp'] = EnzymeMap2[i]['Kapp']

    def estimate_default_Kapps(self, target_mu, compartment_densities_and_PGs=None, flux_bounds=None, eukaryotic=False, plateau_limit=4, mu_approximation_precision=0.0001, transporter_to_lumen_coefficient=10, default_kapp_LB=0, default_kapp_UB=None):
        """
        Parameters
        ----------
        target_mu : float
        compartment_densities : pandas.DataFrame
        compartment_PGs : pandas.DataFrame
        flux_bounds : pandas.DataFrame
        """
        orig_enz = self.rbaSession.model.parameters.functions._elements_by_id[
            'default_efficiency'].parameters._elements_by_id['CONSTANT'].value

        out = pandas.DataFrame()
        for comp in list(compartment_densities_and_PGs['Compartment_ID']):
            self.rbaSession.model.parameters.functions._elements_by_id[str(
                'fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = 0.01*compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density']
            self.rbaSession.model.parameters.functions._elements_by_id[str('fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = 0.01 * \
                compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID']
                                                  == comp, 'PG_fraction']
        self.rbaSession.rebuild_from_model()
        self.rbaSession.addExchangeReactions()

        rxn_LBs = {}
        rxn_UBs = {}
        for rx in flux_bounds['Reaction_ID']:
            lb = flux_bounds.loc[flux_bounds['Reaction_ID'] == rx, 'LB'].values[0]
            ub = flux_bounds.loc[flux_bounds['Reaction_ID'] == rx, 'UB'].values[0]
            if not pandas.isna(lb):
                rxn_LBs.update({rx: lb})
            if not pandas.isna(ub):
                rxn_UBs.update({rx: ub})
        self.rbaSession.Problem.setLB(rxn_LBs)
        self.rbaSession.Problem.setUB(rxn_UBs)

        # self.rbaSession.setMedium(medium)
        if eukaryotic:
            self.rbaSession.eukaryoticDensities_calibration(CompartmentRelationships=False)

        kapp_LB = default_kapp_LB
        if default_kapp_UB is not None:
            kapp_UB = default_kapp_UB
        else:
            kapp_UB = orig_enz*1000
        new_kapp = (kapp_UB+kapp_LB)/2

        self.rbaSession.model.parameters.functions._elements_by_id[
            'default_efficiency'].parameters._elements_by_id['CONSTANT'].value = new_kapp
        self.rbaSession.model.parameters.functions._elements_by_id['default_transporter_efficiency'].parameters._elements_by_id[
            'CONSTANT'].value = transporter_to_lumen_coefficient*new_kapp
        Mu_pred = self.rbaSession.findMaxGrowthRate()

        Mus = []
        Mus_Error = []
        Kapps = []
        last_Mu = numpy.nan
        plateau_count = 0

        while abs(target_mu - Mu_pred) > mu_approximation_precision:
            if plateau_count >= plateau_limit:
                break
            self.rbaSession.model.parameters.functions._elements_by_id[
                'default_efficiency'].parameters._elements_by_id['CONSTANT'].value = new_kapp
            self.rbaSession.model.parameters.functions._elements_by_id['default_transporter_efficiency'].parameters._elements_by_id[
                'CONSTANT'].value = transporter_to_lumen_coefficient*new_kapp
            self.rbaSession.rebuild_from_model()
            self.rbaSession.addExchangeReactions()
            self.rbaSession.Problem.setLB(rxn_LBs)
            self.rbaSession.Problem.setUB(rxn_UBs)
            Mu_pred = self.rbaSession.findMaxGrowthRate()
            Mus_Error.append(abs(target_mu - Mu_pred))
            Mus.append(Mu_pred)
            Kapps.append(new_kapp)
            if Mu_pred > target_mu:
                new_kapp_prelim = kapp_LB+(0.5*abs(kapp_LB-new_kapp))
                kapp_UB = new_kapp
            elif Mu_pred < target_mu:
                new_kapp_prelim = kapp_UB-(0.5*abs(new_kapp-kapp_UB))
                kapp_LB = new_kapp
            new_kapp = new_kapp_prelim
            if len(Mus) > 2:
                if Mus[-2] == Mu_pred:
                    plateau_count += 1
                else:
                    plateau_count = 0

        self.default_kapp_estimation = pandas.DataFrame()
        self.default_kapp_estimation['Mu'] = Mus
        self.default_kapp_estimation['delta_Mu'] = Mus_Error
        self.default_kapp_estimation['default_efficiency'] = Kapps
        self.default_kapp_estimation['default_transporter_efficiency'] = [
            transporter_to_lumen_coefficient*i for i in Kapps]

    def inject_specific_kapps(self, specific_kapps, round_to_digits=0):
        """
        Parameters
        ----------
        specific_kapps : pandas.DataFrame
        """
        parameterized = []
        for enz in list(specific_kapps['Enzyme_ID']):
            if not pandas.isna(specific_kapps.loc[specific_kapps['Enzyme_ID'] == enz, 'Kapp'].values[0]):
                if enz not in parameterized:
                    all_enzs = self.rbaSession.ModelStructure.EnzymeInfo.Elements[enz]['Isozymes']
                    all_enzs.append(enz)
                    parameterized += all_enzs
                    if len(all_enzs) == 1:
                        proto_enz = all_enzs[0]
                    else:
                        proto_enz = [i for i in all_enzs if not '_duplicate_' in i][0]
                    val = round(specific_kapps.loc[specific_kapps['Enzyme_ID']
                                                   == enz, 'Kapp'].values[0], round_to_digits)
                    const = rba.xml.parameters.Function(
                        str(proto_enz + '_kapp__constant'), 'constant', parameters={'CONSTANT': val}, variable=None)
                    self.rbaSession.model.parameters.functions.append(const)
                    count = 0
                    for e in self.rbaSession.model.enzymes.enzymes:
                        if e.id in all_enzs:
                            count += 1
                            e.forward_efficiency = str(proto_enz + '_kapp__constant')
                            e.backward_efficiency = str(proto_enz + '_kapp__constant')
                            if count == len(all_enzs):
                                break
