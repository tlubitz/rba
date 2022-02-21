
import rba
import copy
import pandas
import time
import numpy
import seaborn
import matplotlib.pyplot as plt
from rbatools.rba_Session import RBA_Session
from sklearn.linear_model import LinearRegression
# import matplotlib.pyplot as plt


def find_ribosomal_proteins(rba_session, model_processes=['TranslationC', 'TranslationM'], external_annotations=None):
    out = []
    for i in model_processes:
        out += [rba_session.ModelStructure.ProteinInfo.Elements[j]['ProtoID']
                for j in list(rba_session.ModelStructure.ProcessInfo.Elements[i]['Composition'].keys()) if j in rba_session.ModelStructure.ProteinInfo.Elements.keys()]
    if external_annotations is not None:
        out += list(external_annotations['ID'])
    return(list(set(out)))


def build_model_compartment_map(rba_session):
    out = {rba_session.ModelStructure.ProteinInfo.Elements[i]['ProtoID']: rba_session.ModelStructure.ProteinInfo.Elements[i]['Compartment'] for i in list(
        rba_session.ModelStructure.ProteinInfo.Elements.keys())}
    return(out)


def build_compartment_annotations(Compartment_Annotations_external, model_protein_compartment_map):
    for i in Compartment_Annotations_external.index:
        if Compartment_Annotations_external.loc[i, 'ID'] in list(model_protein_compartment_map.keys()):
            Compartment_Annotations_external.loc[i, 'modelproteinannotation'] = 1
        else:
            Compartment_Annotations_external.loc[i, 'modelproteinannotation'] = 0
    Compartment_Annotations_internal = pandas.DataFrame()
    Compartment_Annotations_internal['ID'] = list(model_protein_compartment_map.keys())
    Compartment_Annotations_internal['ModelComp'] = list(model_protein_compartment_map.values())
    Compartment_Annotations = pandas.concat(
        [Compartment_Annotations_internal, Compartment_Annotations_external.loc[Compartment_Annotations_external['modelproteinannotation'] == 0, ['ID', 'ModelComp']]], axis=0)
    return(Compartment_Annotations)


def build_dataset_annotations(input, ID_column, Uniprot, Compartment_Annotations, model_protein_compartment_map, ribosomal_proteins):
    print('riboprots-----------------')
    print(ribosomal_proteins)
    out = pandas.DataFrame()
    for g in list(input[ID_column]):
        out.loc[g, 'ID'] = g
        matches = [i for i in list(Uniprot.loc[pandas.isna(
            Uniprot['Gene names']) == False, 'Gene names']) if g in i]
        mass_prot = numpy.nan
        if len(matches) > 0:
            mass_prot = len(Uniprot.loc[Uniprot['Gene names'] == matches[0], 'Sequence'].values[0])
        out.loc[g, 'AA_residues'] = mass_prot
        if g in list(Compartment_Annotations['ID']):
            out.loc[g, 'Location'] = Compartment_Annotations.loc[Compartment_Annotations['ID']
                                                                 == g, 'ModelComp'].values[0]
        in_model = 0
        if g in model_protein_compartment_map.keys():
            in_model = 1
        is_ribosomal = 0
        if g in ribosomal_proteins:
            is_ribosomal = 1
        out.loc[g, 'InModel'] = in_model
        out.loc[g, 'IsRibosomal'] = is_ribosomal
    return(out)


def build_full_annotations_from_dataset_annotations(annotations_list):
    out = pandas.concat(annotations_list, axis=0)
    index = out.index
    is_duplicate = index.duplicated(keep="first")
    not_duplicate = ~is_duplicate
    out = out[not_duplicate]
    return(out)


def infer_copy_numbers_from_reference_copy_numbers(fold_changes, absolute_data, matching_column_in_fold_change_data, matching_column_in_absolute_data, conditions_in_fold_change_data_to_restore):
    out = pandas.DataFrame()
    for i in list(absolute_data['Gene']):
        if i in list(fold_changes['Gene']):
            FoldChange_match = fold_changes.loc[fold_changes['Gene']
                                                == i, matching_column_in_fold_change_data].values[0]
            CopyNumber_match = absolute_data.loc[absolute_data['Gene']
                                                 == i, matching_column_in_absolute_data].values[0]
            if not pandas.isna(FoldChange_match):
                if not pandas.isna(CopyNumber_match):
                    out.loc[i, 'ID'] = i
                    out.loc[i, 'Absolute_Reference'] = CopyNumber_match/(2**FoldChange_match)
    for gene in list(out['ID']):
        Abs_Ref = out.loc[gene, 'Absolute_Reference']
        for condition in conditions_in_fold_change_data_to_restore:
            out.loc[gene, condition] = Abs_Ref * \
                (2**fold_changes.loc[fold_changes['Gene'] == gene, condition].values[0])
    return(out)


def add_annotations_to_proteome(input, ID_column, annotations):
    for i in input.index:
        if input.loc[i, ID_column] in annotations.index:
            input.loc[i, 'AA_residues'] = annotations.loc[input.loc[i, ID_column], 'AA_residues']
            input.loc[i, 'Location'] = annotations.loc[input.loc[i, ID_column], 'Location']
            input.loc[i, 'InModel'] = annotations.loc[input.loc[i, ID_column], 'InModel']
            input.loc[i, 'IsRibosomal'] = annotations.loc[input.loc[i, ID_column], 'IsRibosomal']
    return(input)


def determine_compartment_occupation(Data, Condition, mass_col='AA_residues', only_in_model=False, compartments_to_ignore=['DEF'], compartments_no_original_PG=[], ribosomal_proteins_as_extra_compartment=True):
    for i in compartments_to_ignore:
        Data = Data.loc[Data['Location'] != i]
    for i in compartments_no_original_PG:
        Data = Data.loc[(Data['Location'] != i) | (Data['InModel'] == 1)]
    if only_in_model:
        Data = Data.loc[Data['InModel'] >= 1]
    if ribosomal_proteins_as_extra_compartment:
        Data_R = Data.loc[Data['IsRibosomal'] == 1].copy()
        Data = Data.loc[Data['IsRibosomal'] == 0]
        Data_R_df = Data_R.loc[:, [Condition, mass_col, 'Location']]
        Data_R_df[Condition] = Data_R_df[Condition]*Data_R_df[mass_col]
        Ribosomal_sum = Data_R_df[Condition].sum()
    df = Data.loc[:, [Condition, mass_col, 'Location']]
    df[Condition] = df[Condition]*df[mass_col]
    out = pandas.DataFrame(df.groupby('Location').sum())
    if ribosomal_proteins_as_extra_compartment:
        out.loc['Ribosomes', Condition] = Ribosomal_sum
    out.loc['Total', Condition] = out[Condition].sum()
    out.loc[:, 'original_protein_fraction'] = out[Condition]/out.loc['Total', Condition]
    out.rename(columns={Condition: 'original_amino_acid_occupation'}, inplace=True)
    out.drop(columns=['AA_residues'], inplace=True)
    return(out)


def build_proteome_overview(input, condition, compartments_to_ignore=['DEF', 'DEFA', 'Def'], compartments_no_original_PG=['n', 'Secreted'], ribosomal_proteins_as_extra_compartment=True):
    out = determine_compartment_occupation(Data=input, Condition=condition, compartments_to_ignore=compartments_to_ignore,
                                           compartments_no_original_PG=compartments_no_original_PG, ribosomal_proteins_as_extra_compartment=ribosomal_proteins_as_extra_compartment, only_in_model=False)
    out_in_model = determine_compartment_occupation(Data=input, Condition=condition, compartments_to_ignore=compartments_to_ignore,
                                                    compartments_no_original_PG=compartments_no_original_PG, ribosomal_proteins_as_extra_compartment=ribosomal_proteins_as_extra_compartment, only_in_model=True)
    out['original_PG_fraction'] = 1-out_in_model['original_amino_acid_occupation'] / \
        out['original_amino_acid_occupation']
    return(out)


def determine_correction_factor_A(fractions_entirely_replaced_with_expected_value):
    expected_fraction_sum = 0
    for i in fractions_entirely_replaced_with_expected_value.keys():
        expected_fraction_sum += fractions_entirely_replaced_with_expected_value[i]
    factor = 1/(1-expected_fraction_sum)
    return(factor)


def determine_correction_factor_B(imposed_compartment_fractions):
    expected_fractions = 0
    for i in imposed_compartment_fractions.keys():
        expected_fractions += imposed_compartment_fractions[i]
    factor = 1-expected_fractions
    return(factor)


def determine_correction_factor_C(input, condition, reference_condition):
    return(input.loc[input['ID'] == 'Total_protein', condition].values[0]/input.loc[input['ID'] == 'Total_protein', reference_condition].values[0])


def correct_protein_fractions(input, factors, directly_corrected_compartments, imposed_compartment_fractions):
    out = input.copy()
    for c in out.index:
        if c in directly_corrected_compartments:
            out.loc[c, 'new_protein_fraction'] = out.loc[c,
                                                         'original_protein_fraction']*factors['A']*factors['B']
        elif c in imposed_compartment_fractions.keys():
            out.loc[c, 'new_protein_fraction'] = imposed_compartment_fractions[c]
    return(out)


def correct_PG_fraction(input, factors, compartments_no_original_PG, merged_compartments):
    out = input.copy()
    for c in out.index:
        if c == 'Total':
            continue
        else:
            if c in compartments_no_original_PG:
                original_fraction = out.loc[c, 'original_protein_fraction']
                out.loc[c, 'new_PG_fraction'] = 1 - ((factors['A']*factors['B']*original_fraction) /
                                                     out.loc[c, 'new_protein_fraction'])
            elif c in merged_compartments.keys():
                out.loc[c, 'new_PG_fraction'] = out.loc[c, 'original_PG_fraction']*out.loc[c, 'original_protein_fraction']/(
                    out.loc[c, 'original_protein_fraction']+out.loc[merged_compartments[c], 'original_protein_fraction'])
            else:
                out.loc[c, 'new_PG_fraction'] = out.loc[c, 'original_PG_fraction']
    return(out)


def merge_compartments(input, merged_compartments):
    out = input.copy()
    for c in merged_compartments.keys():
        out.loc[c, 'new_protein_fraction'] = out.loc[c, 'new_protein_fraction'] + \
            out.loc[merged_compartments[c], 'new_protein_fraction']
    return(out)


def calculate_new_total_PG_fraction(input):
    out = input.copy()
    fraction = 0
    for c in out.index:
        if c not in ['Total', 'Ribosomes']:
            fraction += out.loc[c, 'new_protein_fraction']*out.loc[c, 'new_PG_fraction']
    out.loc['Total', 'new_PG_fraction'] = fraction
    out.loc['Total', 'new_protein_fraction'] = 1
    return(out)


def determine_apparent_process_efficiencies(growth_rate, input, rba_session, proteome_summary, protein_data, condition, gene_id_col):
    process_efficiencies = pandas.DataFrame()
    for i in input.index:
        process_ID = input.loc[i, 'Process_ID']
        process_name = input.loc[i, 'Process_Name']
        process_client_compartments = input.loc[i, 'Client_Compartments'].split(' , ')
        constituting_proteins = {rba_session.ModelStructure.ProteinInfo.Elements[i]['ProtoID']: rba_session.ModelStructure.ProteinInfo.Elements[
            i]['AAnumber'] for i in rba_session.ModelStructure.ProcessInfo.Elements[process_name]['Composition'].keys()}
        Total_client_fraction = sum([proteome_summary.loc[i, 'new_protein_fraction']
                                     for i in process_client_compartments])
        n_AAs_in_machinery = 0
        machinery_size = 0
        for i in constituting_proteins.keys():
            if i in protein_data['ID']:
                protein_data.loc[protein_data['ID'] == i, ]
                n_AAs_in_machinery += protein_data.loc[protein_data['ID'] == i, condition].values[0] * \
                    protein_data.loc[protein_data['ID'] == i, 'AA_residues'].values[0]
                machinery_size += constituting_proteins[i]
        # right reference amounth?
        if n_AAs_in_machinery > 0:
            relative_Protein_fraction_of_machinery = n_AAs_in_machinery / \
                proteome_summary.loc['Total', 'original_amino_acid_occupation']
            specific_capacity = growth_rate*Total_client_fraction/relative_Protein_fraction_of_machinery
            apparent_capacity = specific_capacity*machinery_size
            # process_ID[process_name] = apparent_capacity
            process_efficiencies.loc[process_name, 'Process'] = process_ID
            process_efficiencies.loc[process_name, 'Parameter'] = str(
                process_ID+'_apparent_efficiency')
            process_efficiencies.loc[process_name, 'Value'] = apparent_capacity
    return(process_efficiencies)


def correction_pipeline(input, condition, compartments_to_ignore, compartments_no_original_PG, fractions_entirely_replaced_with_expected_value, imposed_compartment_fractions, directly_corrected_compartments, merged_compartments):
    out = build_proteome_overview(input=input, condition=condition, compartments_to_ignore=compartments_to_ignore,
                                  compartments_no_original_PG=compartments_no_original_PG, ribosomal_proteins_as_extra_compartment=True)
    factor_A = determine_correction_factor_A(fractions_entirely_replaced_with_expected_value={
                                             i: imposed_compartment_fractions[i] for i in fractions_entirely_replaced_with_expected_value})
    factor_B = determine_correction_factor_B(
        imposed_compartment_fractions=imposed_compartment_fractions)
    out = correct_protein_fractions(input=out, factors={
                                    'A': factor_A, 'B': factor_B}, directly_corrected_compartments=directly_corrected_compartments, imposed_compartment_fractions=imposed_compartment_fractions)
    out = correct_PG_fraction(input=out, factors={
                              'A': factor_A, 'B': factor_B}, compartments_no_original_PG=compartments_no_original_PG, merged_compartments=merged_compartments)
    out = merge_compartments(input=out, merged_compartments=merged_compartments)
    out = calculate_new_total_PG_fraction(input=out)
    out.to_csv(str('Correction_overview_'+condition+'.csv'))
    return({'Summary': out, 'Correction_factors': {'A': factor_A, 'B': factor_B}})


def build_input_for_default_kapp_estimation(input):
    out = pandas.DataFrame(columns=['Compartment_ID', 'Density', 'PG_fraction'])
    for i in input['Summary'].index:
        if i not in ['Total', 'Ribosomes']:
            out.loc[i, 'Compartment_ID'] = i
            out.loc[i, 'Density'] = input['Summary'].loc[i, 'new_protein_fraction']
            out.loc[i, 'PG_fraction'] = input['Summary'].loc[i, 'new_PG_fraction']
    return(out)


def flux_bounds_from_input(input, condition, specific_exchanges=None):
    flux_mean_df = input.loc[input['Type'] == 'ExchangeFlux_Mean', :]
    flux_mean_SE = input.loc[input['Type'] == 'ExchangeFlux_StandardError', :]
    out = pandas.DataFrame(columns=['Reaction_ID', 'LB', 'UB'])
    if specific_exchanges is None:
        exchanges_to_set = list(flux_mean_df['ID'])
    else:
        exchanges_to_set = specific_exchanges
    for rx in exchanges_to_set:
        mean_val = flux_mean_df.loc[flux_mean_df['ID'] == rx, condition].values[0]
        if not pandas.isna(mean_val):
            SE_val = flux_mean_SE.loc[flux_mean_SE['ID'] == str(rx+'_SE'), condition].values[0]
            out.loc[rx, 'Reaction_ID'] = rx
            if not pandas.isna(SE_val):
                lb = mean_val-SE_val
                ub = mean_val+SE_val
                if mean_val < 0:
                    out.loc[rx, 'LB'] = lb
                    if ub > 0:
                        out.loc[rx, 'UB'] = 0
                    else:
                        out.loc[rx, 'UB'] = ub
                elif mean_val > 0:
                    out.loc[rx, 'UB'] = ub
                    if lb < 0:
                        out.loc[rx, 'LB'] = 0
                    else:
                        out.loc[rx, 'LB'] = lb
                else:
                    out.loc[rx, 'LB'] = lb
                    out.loc[rx, 'UB'] = ub
            else:
                out.loc[rx, 'LB'] = mean_val
                out.loc[rx, 'UB'] = mean_val
    flux_dir_df = input.loc[input['Type'] == 'Flux_Direction', :]
    if specific_exchanges is None:
        exchanges_to_set = list(flux_dir_df['ID'])
    else:
        exchanges_to_set = specific_exchanges
    for rx in exchanges_to_set:
        out.loc[rx, 'Reaction_ID'] = rx
        if flux_dir_df.loc[flux_dir_df['ID'] == rx, condition].values[0] == 1:
            out.loc[rx, 'LB'] = 0
        elif flux_dir_df.loc[flux_dir_df['ID'] == rx, condition].values[0] == -1:
            out.loc[rx, 'UB'] = 0
        elif flux_dir_df.loc[flux_dir_df['ID'] == rx, condition].values[0] == 0:
            out.loc[rx, 'LB'] = 0
            out.loc[rx, 'UB'] = 0
    flux_upper_df = input.loc[input['Type'] == 'Flux_Upper_Bound', :]
    for rx in list(flux_upper_df['ID']):
        out.loc[rx, 'Reaction_ID'] = rx
        out.loc[rx, 'UB'] = flux_upper_df.loc[flux_upper_df['ID'] == rx, condition].values[0]
    flux_lower_df = input.loc[input['Type'] == 'Flux_Lower_Bound', :]
    for rx in list(flux_lower_df['ID']):
        out.loc[rx, 'Reaction_ID'] = rx
        out.loc[rx, 'LB'] = flux_lower_df.loc[flux_lower_df['ID'] == rx, condition].values[0]
    return(out)


def growth_Rate_from_input(input, condition):
    return(input.loc[input['Type'] == 'Growth_Rate', condition].values[0])


def proteome_fractions_from_input(input, condition):
    df = input.loc[input['Type'] == 'Expected_ProteomeFraction', :]
    return(dict(zip(list(df['ID']), list(df[condition]))))


def medium_concentrations_from_input(input, condition):
    df = input.loc[input['Type'] == 'Medium_Concentration', :]
    return(dict(zip(list(df['ID']), list(df[condition]))))


def build_input_proteome_for_specific_kapp_estimation(proteomics_data, condition):
    out = pandas.DataFrame()
    out['ID'] = proteomics_data['ID']
    out['copy_number'] = proteomics_data[condition]
    return(out)


def inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps=None, process_efficiencies=None, round_to_digits=0):
    """
    Parameters
    ----------
    specific_kapps : pandas.DataFrame(columns=['Enzyme_ID','Kapp'])
    default_kapps : {'default_kapp':value,'default_transporter_kapp':value}
    process_efficiencies : pandas.DataFrame(columns=['Process','Parameter','Value'])
    """
    if specific_kapps is not None:
        parameterized = []
        for enz in list(specific_kapps['Enzyme_ID']):
            if not pandas.isna(specific_kapps.loc[specific_kapps['Enzyme_ID'] == enz, 'Kapp'].values[0]):
                if enz not in parameterized:
                    all_enzs = rba_session.ModelStructure.EnzymeInfo.Elements[enz]['Isozymes']
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
                    if str(proto_enz + '_kapp__constant') not in rba_session.model.parameters.functions._elements_by_id.keys():
                        rba_session.model.parameters.functions.append(const)
                    else:
                        rba_session.model.parameters.functions._elements_by_id[const.id] = const
                    count = 0
                    for e in rba_session.model.enzymes.enzymes:
                        if e.id in all_enzs:
                            count += 1
                            e.forward_efficiency = str(proto_enz + '_kapp__constant')
                            e.backward_efficiency = str(proto_enz + '_kapp__constant')
                            if count == len(all_enzs):
                                break

    if default_kapps is not None:
        if type(default_kapps) is dict:
            rba_session.model.parameters.functions._elements_by_id[
                'default_efficiency'].parameters._elements_by_id['CONSTANT'].value = default_kapps['default_kapp']
            rba_session.model.parameters.functions._elements_by_id['default_transporter_efficiency'].parameters._elements_by_id[
                'CONSTANT'].value = default_kapps['default_transporter_kapp']

    if process_efficiencies is not None:
        for i in process_efficiencies.index:
            if process_efficiencies.loc[i, 'Process'] in rba_session.model.processes.processes._elements_by_id.keys():
                if not pandas.isna(process_efficiencies.loc[i, 'Value']):
                    rba_session.model.processes.processes._elements_by_id[process_efficiencies.loc[i,
                                                                                                   'Process']].machinery.capacity.value = process_efficiencies.loc[i, 'Parameter']
                    const = rba.xml.parameters.Function(process_efficiencies.loc[i, 'Parameter'], 'constant', parameters={
                        'CONSTANT': process_efficiencies.loc[i, 'Value']}, variable=None)
                    if process_efficiencies.loc[i, 'Parameter'] not in rba_session.model.parameters.functions._elements_by_id.keys():
                        rba_session.model.parameters.functions.append(const)
                    else:
                        rba_session.model.parameters.functions._elements_by_id[const.id] = const

    rba_session.rebuild_from_model()


def calibration_workflow(proteome,
                         condition,
                         reference_condition,
                         gene_ID_column,
                         definition_file,
                         rba_session,
                         process_efficiency_estimation_input=None,
                         default_kapps_provided=None):
    t0 = time.time()
    correction_results = correction_pipeline(input=proteome,
                                             condition=condition,
                                             compartments_to_ignore=['DEF', 'DEFA', 'Def'],
                                             compartments_no_original_PG=['n', 'Secreted'],
                                             fractions_entirely_replaced_with_expected_value=[
                                                 'Ribosomes'],
                                             imposed_compartment_fractions=proteome_fractions_from_input(
                                                 input=definition_file, condition=condition),
                                             directly_corrected_compartments=[
                                                 'c', 'cM', 'erM', 'gM', 'm', 'mIM', 'mIMS', 'mOM', 'vM', 'x'],
                                             merged_compartments={'c': 'Ribosomes'})

    # mumax0 = rba_session.findMaxGrowthRate()
    rba_session.setMedium(medium_concentrations_from_input(
        input=definition_file, condition=condition))

    # mumax1 = rba_session.findMaxGrowthRate()
    if process_efficiency_estimation_input is not None:
        process_efficiencies = determine_apparent_process_efficiencies(growth_rate=growth_Rate_from_input(
            input=definition_file, condition=condition), input=process_efficiency_estimation_input, rba_session=rba_session, protein_data=proteome, proteome_summary=correction_results['Summary'], condition=condition, gene_id_col=gene_ID_column)
        inject_estimated_efficiencies_into_model(
            rba_session, specific_kapps=None, default_kapps=None, process_efficiencies=process_efficiencies)
    else:
        process_efficiencies = None
    protein_scaling_coefficient = 1000 * determine_correction_factor_C(input=definition_file, condition=condition, reference_condition=reference_condition) * \
        correction_results['Correction_factors']['A'] * \
        correction_results['Correction_factors']['B']/6.022e23
#    protein_scaling_coefficient = 1000 * correction_results['Correction_factors']['A'] * correction_results['Correction_factors']['B']/6.022e23
    proteome[condition] *= protein_scaling_coefficient
    Specific_Kapps = rba_session.estimate_specific_Kapps(proteomicsData=build_input_proteome_for_specific_kapp_estimation(proteome, condition),
                                                         flux_bounds=flux_bounds_from_input(
                                                             input=definition_file, condition=condition, specific_exchanges=None),
                                                         mu=growth_Rate_from_input(
        input=definition_file, condition=condition),
        biomass_function=None,
        target_biomass_function=True)
    # Specific_Kapps.loc[(Specific_Kapps['Kapp'] <= 1000000) &
    #                   (Specific_Kapps['Kapp'] >= 1), 'Kapp'].hist()
    # plt.show()

    # mumax2 = rba_session.findMaxGrowthRate()
    if default_kapps_provided is None:
        Default_Kapps = rba_session.estimate_default_Kapps(target_mu=growth_Rate_from_input(input=definition_file, condition=condition), compartment_densities_and_PGs=build_input_for_default_kapp_estimation(
            correction_results), flux_bounds=flux_bounds_from_input(input=definition_file, condition=condition, specific_exchanges=None), mu_approximation_precision=0.01)
        inject_estimated_efficiencies_into_model(rba_session, specific_kapps=None, default_kapps={
                                                 'default_kapp': Default_Kapps.iloc[-1, 2], 'default_transporter_kapp': Default_Kapps.iloc[-1, 3]}, process_efficiencies=None)
    else:
        inject_estimated_efficiencies_into_model(
            rba_session, specific_kapps=None, default_kapps=default_kapps_provided, process_efficiencies=None)
        Default_Kapps = default_kapps_provided
    inject_estimated_efficiencies_into_model(
        rba_session, specific_kapps=Specific_Kapps, default_kapps=None, process_efficiencies=None)

    # mumax3 = rba_session.findMaxGrowthRate()

    compartment_densities_and_PGs = build_input_for_default_kapp_estimation(correction_results)
    for comp in list(compartment_densities_and_PGs['Compartment_ID']):
        rba_session.model.parameters.functions._elements_by_id[str(
            'fraction_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'Density']
        rba_session.model.parameters.functions._elements_by_id[str(
            'fraction_non_enzymatic_protein_'+comp)].parameters._elements_by_id['CONSTANT'].value = compartment_densities_and_PGs.loc[compartment_densities_and_PGs['Compartment_ID'] == comp, 'PG_fraction']
    rba_session.rebuild_from_model()
    rba_session.addExchangeReactions()
    rba_session.setMedium(medium_concentrations_from_input(
        input=definition_file, condition=condition))
    # FBs = flux_bounds_from_input(
    #    input=definition_file, condition=condition, specific_exchanges=None)
    #rba_session.Problem.setLB(dict(zip(list(FBs['Reaction_ID']), list(FBs['LB']))))
    # rba_session.Problem.setLB({FBs.loc[i, 'Reaction_ID']: FBs.loc[i, 'LB']
    #                           for i in FBs.index if not pandas.isna(FBs.loc[i, 'LB'])})
    # rba_session.Problem.setLB({FBs.loc[i, 'Reaction_ID']: FBs.loc[i, 'UB']
    #                           for i in FBs.index if not pandas.isna(FBs.loc[i, 'UB'])})
    #rba_session.Problem.setUB(dict(zip(list(FBs['Reaction_ID']), list(FBs['UB']))))
    rba_session.Problem.setLB({'R_EX_cys__L_e': 0, 'R_EX_met__L_e': 0})
    rba_session.Problem.setUB({'R_EX_cys__L_e': 0, 'R_EX_met__L_e': 0})
    mumax4 = rba_session.findMaxGrowthRate()
    rba_session.recordResults('Prokaryotic')
    prok_results = copy.deepcopy(rba_session.Results)
    rba_session2 = copy.copy(rba_session)
    rba_session2.eukaryoticDensities4(CompartmentRelationships=False)
    mumax5 = rba_session2.findMaxGrowthRate()
    rba_session2.recordResults('Eukaryotic')
    # print([Default_Kapps.iloc[-1, 2], Default_Kapps.iloc[-1, 3]])
    # print([growth_Rate_from_input(input=definition_file,
    #                              condition=condition), mumax0, mumax1, mumax2, mumax3, mumax4, mumax5])
    print(time.time() - t0)
    return({'Simulation_Results': prok_results, 'Simulation_Results_Euk': copy.deepcopy(rba_session2.Results), 'Proteome': build_input_proteome_for_specific_kapp_estimation(proteome, condition), 'Correction_Results': correction_results, 'Default_Kapps': Default_Kapps, 'Specific_Kapps': Specific_Kapps, 'Process_Efficiencies': process_efficiencies})
#    seaborn.violinplot(x=Specific_Kapps.loc[Specific_Kapps['Kapp'] <= 400000, 'Kapp'])
    # Specific_Kapps.loc[(Specific_Kapps['Kapp'] <= 1000000) &
    #                               (Specific_Kapps['Kapp'] >= 1), 'Kapp']).hist()
    # plt.show()
    # Test predictions
    # Given medium predict Mu, Exchanges and Proteome
    # Prokaryotic
    # Eukaryotic


# 1. import model and uniprot-file and compartment-annotation

## external_annotations for ribosomal-proteins!!! ##
## process-efficiency estimation input ##
## parse input-data properly and add Lahtvee information ##


print('---------------------START----------------------')
Input_Data = pandas.read_csv(
    'DataSetsYeastRBACalibration/Calibration_InputDefinition.csv', sep=';', decimal=',', index_col=0)
Process_Efficiency_Estimation_Input = pandas.read_csv(
    'DataSetsYeastRBACalibration/Process_Efficiency_Estimation_Input.csv', sep=';', decimal=',')
Simulation = RBA_Session('Yeast_iMM904_RBA_model')

Uniprot = pandas.read_csv('Yeast_iMM904_RBA_model/uniprot.csv', sep='\t')

Compartment_Annotations_external = pandas.read_csv(
    'DataSetsYeastRBACalibration/Manually_curated_Protein_Locations_for_Calibration.csv', index_col=None, sep=';')
Ribosomal_Proteins_Uniprot = pandas.read_csv(
    'DataSetsYeastRBACalibration/uniprot_ribosomal_proteins.csv', index_col=None, sep=';')

Hackett_Clim_FCs = pandas.read_csv('DataSetsYeastRBACalibration/Hacket_Clim_ProteinFCs.csv')

Lahtvee_REF = pandas.read_csv('DataSetsYeastRBACalibration/LahtveeRefProteomicsData.csv')
picogram_togram_coefficient = 1e12
Lahtvee_REF['Lahtvee_REF'] *= picogram_togram_coefficient
Lahtvee_REF = Lahtvee_REF.loc[pandas.isna(Lahtvee_REF['Lahtvee_REF']) == False]

ribosomal_proteins = find_ribosomal_proteins(rba_session=Simulation, model_processes=[
    'TranslationC', 'TranslationM'], external_annotations=Ribosomal_Proteins_Uniprot)

model_protein_compartment_map = build_model_compartment_map(rba_session=Simulation)

Compartment_Annotations = build_compartment_annotations(
    Compartment_Annotations_external=Compartment_Annotations_external, model_protein_compartment_map=model_protein_compartment_map)
print('Annotations to data')
annotations_Lahtvee = build_dataset_annotations(input=Lahtvee_REF, ID_column='Gene', Uniprot=Uniprot,
                                                Compartment_Annotations=Compartment_Annotations, model_protein_compartment_map=model_protein_compartment_map, ribosomal_proteins=ribosomal_proteins)
annotations_Hackett = build_dataset_annotations(input=Hackett_Clim_FCs, ID_column='Gene', Uniprot=Uniprot,
                                                Compartment_Annotations=Compartment_Annotations, model_protein_compartment_map=model_protein_compartment_map, ribosomal_proteins=ribosomal_proteins)
full_annotations = build_full_annotations_from_dataset_annotations(
    annotations_list=[annotations_Lahtvee, annotations_Hackett])


####### Bootstrapping-loop starts here #######
restored_Hackett_Data = infer_copy_numbers_from_reference_copy_numbers(fold_changes=Hackett_Clim_FCs, absolute_data=Lahtvee_REF, matching_column_in_fold_change_data='Hackett_C01',
                                                                       matching_column_in_absolute_data='Lahtvee_REF', conditions_in_fold_change_data_to_restore=['Hackett_C005', 'Hackett_C01', 'Hackett_C016', 'Hackett_C022', 'Hackett_C03'])

restored_Hackett_Data = add_annotations_to_proteome(
    input=restored_Hackett_Data, ID_column='ID', annotations=full_annotations)
Lahtvee_REF = add_annotations_to_proteome(
    input=Lahtvee_REF, ID_column='Gene', annotations=full_annotations)


# default_kapps_provided={'default_kapp':39673 , 'default_transporter_kapp':396730 }
# default_kapps_provided={'default_kapp':85449 , 'default_transporter_kapp':854490 }
# default_kapps_provided={'default_kapp':128174 , 'default_transporter_kapp':1281740 }
# default_kapps_provided={'default_kapp':280762 , 'default_transporter_kapp':2807620 }
# default_kapps_provided = {'default_kapp': 268555, 'default_transporter_kapp': 2685550}

Simulation = RBA_Session('Yeast_iMM904_RBA_model')
Calibration_Hackett_C005 = calibration_workflow(proteome=restored_Hackett_Data, condition='Hackett_C005', reference_condition='Lahtvee_REF', gene_ID_column='Gene',
                                                definition_file=Input_Data, rba_session=Simulation, process_efficiency_estimation_input=Process_Efficiency_Estimation_Input, default_kapps_provided={'default_kapp': 39673, 'default_transporter_kapp': 396730})

print('0.05')
print('')
print('')
print('')
print('')
Simulation = RBA_Session('Yeast_iMM904_RBA_model')
Calibration_Hackett_C01 = calibration_workflow(proteome=restored_Hackett_Data, condition='Hackett_C01', reference_condition='Lahtvee_REF', gene_ID_column='Gene',
                                               definition_file=Input_Data, rba_session=Simulation, process_efficiency_estimation_input=Process_Efficiency_Estimation_Input, default_kapps_provided={'default_kapp': 85449, 'default_transporter_kapp': 854490})
print('0.1')
print('')
print('')
print('')
print('')
Simulation = RBA_Session('Yeast_iMM904_RBA_model')
Calibration_Hackett_C016 = calibration_workflow(proteome=restored_Hackett_Data, condition='Hackett_C016', reference_condition='Lahtvee_REF', gene_ID_column='Gene',
                                                definition_file=Input_Data, rba_session=Simulation, process_efficiency_estimation_input=Process_Efficiency_Estimation_Input, default_kapps_provided={'default_kapp': 128174, 'default_transporter_kapp': 1281740})
print('0.16')
print('')
print('')
print('')
print('')
Simulation = RBA_Session('Yeast_iMM904_RBA_model')
Calibration_Hackett_C022 = calibration_workflow(proteome=restored_Hackett_Data, condition='Hackett_C022', reference_condition='Lahtvee_REF', gene_ID_column='Gene',
                                                definition_file=Input_Data, rba_session=Simulation, process_efficiency_estimation_input=Process_Efficiency_Estimation_Input, default_kapps_provided={'default_kapp': 280762, 'default_transporter_kapp': 2807620})
print('0.22')
print('')
print('')
print('')
print('')
Simulation = RBA_Session('Yeast_iMM904_RBA_model')
Calibration_Hackett_C03 = calibration_workflow(proteome=restored_Hackett_Data, condition='Hackett_C03', reference_condition='Lahtvee_REF', gene_ID_column='Gene',
                                               definition_file=Input_Data, rba_session=Simulation, process_efficiency_estimation_input=Process_Efficiency_Estimation_Input, default_kapps_provided={'default_kapp': 280762, 'default_transporter_kapp': 2807620})
print('0.3')

specKapps_005 = pandas.DataFrame(index=list(
    Calibration_Hackett_C005['Specific_Kapps']['Enzyme_ID']))
specKapps_005['Hackett_C005'] = list(Calibration_Hackett_C005['Specific_Kapps']['Kapp'])

specKapps_01 = pandas.DataFrame(index=list(Calibration_Hackett_C01['Specific_Kapps']['Enzyme_ID']))
specKapps_01['Hackett_C01'] = list(Calibration_Hackett_C01['Specific_Kapps']['Kapp'])

specKapps_016 = pandas.DataFrame(index=list(
    Calibration_Hackett_C016['Specific_Kapps']['Enzyme_ID']))
specKapps_016['Hackett_C016'] = list(Calibration_Hackett_C016['Specific_Kapps']['Kapp'])

specKapps_022 = pandas.DataFrame(index=list(
    Calibration_Hackett_C022['Specific_Kapps']['Enzyme_ID']))
specKapps_022['Hackett_C022'] = list(Calibration_Hackett_C022['Specific_Kapps']['Kapp'])

specKapps_03 = pandas.DataFrame(index=list(Calibration_Hackett_C03['Specific_Kapps']['Enzyme_ID']))
specKapps_03['Hackett_C03'] = list(Calibration_Hackett_C03['Specific_Kapps']['Kapp'])

all_spec_Kapps = pandas.concat(
    [specKapps_005, specKapps_01, specKapps_016, specKapps_022, specKapps_03], axis=1)
all_spec_Kapps['ID'] = all_spec_Kapps.index
all_spec_Kapps.to_csv('Specific_Kapps_out.csv', sep=';', decimal=',')

process_efficiencies_005 = pandas.DataFrame(index=list(
    Calibration_Hackett_C005['Process_Efficiencies']['Process']))
process_efficiencies_005['Hackett_C005'] = list(
    Calibration_Hackett_C005['Process_Efficiencies']['Value'])

process_efficiencies_01 = pandas.DataFrame(index=list(
    Calibration_Hackett_C01['Process_Efficiencies']['Process']))
process_efficiencies_01['Hackett_C01'] = list(
    Calibration_Hackett_C01['Process_Efficiencies']['Value'])

process_efficiencies_016 = pandas.DataFrame(index=list(
    Calibration_Hackett_C016['Process_Efficiencies']['Process']))
process_efficiencies_016['Hackett_C016'] = list(
    Calibration_Hackett_C016['Process_Efficiencies']['Value'])

process_efficiencies_022 = pandas.DataFrame(index=list(
    Calibration_Hackett_C022['Process_Efficiencies']['Process']))
process_efficiencies_022['Hackett_C022'] = list(
    Calibration_Hackett_C022['Process_Efficiencies']['Value'])

process_efficiencies_03 = pandas.DataFrame(index=list(
    Calibration_Hackett_C03['Process_Efficiencies']['Process']))
process_efficiencies_03['Hackett_C03'] = list(
    Calibration_Hackett_C03['Process_Efficiencies']['Value'])

all_process_efficiencies = pandas.concat(
    [process_efficiencies_005, process_efficiencies_01, process_efficiencies_016, process_efficiencies_022, process_efficiencies_03], axis=1)
all_process_efficiencies['ID'] = all_process_efficiencies.index
all_process_efficiencies.to_csv('Process_efficiencies_out.csv', sep=';', decimal=',')

########
########
Mus_o2 = [0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.28, 0.3, 0.35, 0.4]
O2_J = [0.8,  1.3, 2.5, 3.9, 5.3, 7, 7.4, 6.1, 5.1, 3.7]
Glc_J = [0.3, 0.6, 1.1, 1.7, 2.3, 2.8, 3.4, 4.5, 8.6, 11.1]
CO2_J = [0.8,  1.4, 2.7, 4.2, 5.7, 7.5, 8, 8.8, 14.9, 18.9]
EtOH_J = [0, 0, 0, 0, 0, 0, 0.11, 2.3, 9.5, 13.9]
Ac_J = [0, 0, 0, 0, 0, 0, 0.08, 0.41, 0.62, 0.6]
Glyc_J = [0, 0, 0, 0, 0, 0, 0, 0, 0.05, 0.15]
## Hackett#
Mu_Hackett = [0.0498630244, 0.1054314572, 0.154377453333333, 0.2126503108, 0.293841410333333]
Glc_Hackett = [0.7367, 1.5462, 2.1722, 5.1571, 9.5962]
EtOH_Hackett = [0.0127, 0.0529, 0.1084, 4.6066, 14.0672]
Ac_Hackett = [0.0017, 0.0031, 0.0052, 0.4433, 0.8851]
Glyc_Hackett = [0.0035, 0.0077, 0.0065, 0.0579, 0.1699]
conditions = ['Hackett_C005', 'Hackett_C01', 'Hackett_C016', 'Hackett_C022', 'Hackett_C03']


Mus_predicted = [Calibration_Hackett_C005['Simulation_Results']['Mu'].loc['Mu', 'Prokaryotic'],
                 Calibration_Hackett_C01['Simulation_Results']['Mu'].loc['Mu', 'Prokaryotic'],
                 Calibration_Hackett_C016['Simulation_Results']['Mu'].loc['Mu', 'Prokaryotic'],
                 Calibration_Hackett_C022['Simulation_Results']['Mu'].loc['Mu', 'Prokaryotic'],
                 Calibration_Hackett_C03['Simulation_Results']['Mu'].loc['Mu', 'Prokaryotic']]

Mus_predicted_euk = [Calibration_Hackett_C005['Simulation_Results_Euk']['Mu'].loc['Mu', 'Eukaryotic'],
                     Calibration_Hackett_C01['Simulation_Results_Euk']['Mu'].loc['Mu', 'Eukaryotic'],
                     Calibration_Hackett_C016['Simulation_Results_Euk']['Mu'].loc['Mu', 'Eukaryotic'],
                     Calibration_Hackett_C022['Simulation_Results_Euk']['Mu'].loc['Mu', 'Eukaryotic'],
                     Calibration_Hackett_C03['Simulation_Results_Euk']['Mu'].loc['Mu', 'Eukaryotic']]

Glc_Exchange_predicted = [abs(Calibration_Hackett_C005['Simulation_Results']['ExchangeFluxes'].loc['M_glc__D', 'Prokaryotic']),
                          abs(Calibration_Hackett_C01['Simulation_Results']
                              ['ExchangeFluxes'].loc['M_glc__D', 'Prokaryotic']),
                          abs(Calibration_Hackett_C016['Simulation_Results']
                              ['ExchangeFluxes'].loc['M_glc__D', 'Prokaryotic']),
                          abs(Calibration_Hackett_C022['Simulation_Results']
                              ['ExchangeFluxes'].loc['M_glc__D', 'Prokaryotic']),
                          abs(Calibration_Hackett_C03['Simulation_Results']['ExchangeFluxes'].loc['M_glc__D', 'Prokaryotic'])]

EtOH_Exchange_predicted = [abs(Calibration_Hackett_C005['Simulation_Results']['ExchangeFluxes'].loc['M_etoh', 'Prokaryotic']),
                           abs(Calibration_Hackett_C01['Simulation_Results']
                               ['ExchangeFluxes'].loc['M_etoh', 'Prokaryotic']),
                           abs(Calibration_Hackett_C016['Simulation_Results']
                               ['ExchangeFluxes'].loc['M_etoh', 'Prokaryotic']),
                           abs(Calibration_Hackett_C022['Simulation_Results']
                               ['ExchangeFluxes'].loc['M_etoh', 'Prokaryotic']),
                           abs(Calibration_Hackett_C03['Simulation_Results']['ExchangeFluxes'].loc['M_etoh', 'Prokaryotic'])]

Ac_Exchange_predicted = [abs(Calibration_Hackett_C005['Simulation_Results']['ExchangeFluxes'].loc['M_ac', 'Prokaryotic']),
                         abs(Calibration_Hackett_C01['Simulation_Results']
                             ['ExchangeFluxes'].loc['M_ac', 'Prokaryotic']),
                         abs(Calibration_Hackett_C016['Simulation_Results']
                             ['ExchangeFluxes'].loc['M_ac', 'Prokaryotic']),
                         abs(Calibration_Hackett_C022['Simulation_Results']
                             ['ExchangeFluxes'].loc['M_ac', 'Prokaryotic']),
                         abs(Calibration_Hackett_C03['Simulation_Results']['ExchangeFluxes'].loc['M_ac', 'Prokaryotic'])]

O2_Exchange_predicted = [abs(Calibration_Hackett_C005['Simulation_Results']['ExchangeFluxes'].loc['M_o2', 'Prokaryotic']),
                         abs(Calibration_Hackett_C01['Simulation_Results']
                             ['ExchangeFluxes'].loc['M_o2', 'Prokaryotic']),
                         abs(Calibration_Hackett_C016['Simulation_Results']
                             ['ExchangeFluxes'].loc['M_o2', 'Prokaryotic']),
                         abs(Calibration_Hackett_C022['Simulation_Results']
                             ['ExchangeFluxes'].loc['M_o2', 'Prokaryotic']),
                         abs(Calibration_Hackett_C03['Simulation_Results']['ExchangeFluxes'].loc['M_o2', 'Prokaryotic'])]

Glycerol_Exchange_predicted = [abs(Calibration_Hackett_C005['Simulation_Results']['ExchangeFluxes'].loc['M_glyc', 'Prokaryotic']),
                               abs(Calibration_Hackett_C01['Simulation_Results']
                                   ['ExchangeFluxes'].loc['M_glyc', 'Prokaryotic']),
                               abs(Calibration_Hackett_C016['Simulation_Results']
                                   ['ExchangeFluxes'].loc['M_glyc', 'Prokaryotic']),
                               abs(Calibration_Hackett_C022['Simulation_Results']
                                   ['ExchangeFluxes'].loc['M_glyc', 'Prokaryotic']),
                               abs(Calibration_Hackett_C03['Simulation_Results']['ExchangeFluxes'].loc['M_glyc', 'Prokaryotic'])]

###
Glc_Exchange_predicted_euk = [abs(Calibration_Hackett_C005['Simulation_Results_Euk']['ExchangeFluxes'].loc['M_glc__D', 'Eukaryotic']),
                              abs(Calibration_Hackett_C01['Simulation_Results_Euk']
                                  ['ExchangeFluxes'].loc['M_glc__D', 'Eukaryotic']),
                              abs(Calibration_Hackett_C016['Simulation_Results_Euk']
                                  ['ExchangeFluxes'].loc['M_glc__D', 'Eukaryotic']),
                              abs(Calibration_Hackett_C022['Simulation_Results_Euk']
                                  ['ExchangeFluxes'].loc['M_glc__D', 'Eukaryotic']),
                              abs(Calibration_Hackett_C03['Simulation_Results_Euk']['ExchangeFluxes'].loc['M_glc__D', 'Eukaryotic'])]

EtOH_Exchange_predicted_euk = [abs(Calibration_Hackett_C005['Simulation_Results_Euk']['ExchangeFluxes'].loc['M_etoh', 'Eukaryotic']),
                               abs(Calibration_Hackett_C01['Simulation_Results_Euk']
                                   ['ExchangeFluxes'].loc['M_etoh', 'Eukaryotic']),
                               abs(Calibration_Hackett_C016['Simulation_Results_Euk']
                                   ['ExchangeFluxes'].loc['M_etoh', 'Eukaryotic']),
                               abs(Calibration_Hackett_C022['Simulation_Results_Euk']
                                   ['ExchangeFluxes'].loc['M_etoh', 'Eukaryotic']),
                               abs(Calibration_Hackett_C03['Simulation_Results_Euk']['ExchangeFluxes'].loc['M_etoh', 'Eukaryotic'])]

Ac_Exchange_predicted_euk = [abs(Calibration_Hackett_C005['Simulation_Results_Euk']['ExchangeFluxes'].loc['M_ac', 'Eukaryotic']),
                             abs(Calibration_Hackett_C01['Simulation_Results_Euk']
                                 ['ExchangeFluxes'].loc['M_ac', 'Eukaryotic']),
                             abs(Calibration_Hackett_C016['Simulation_Results_Euk']
                                 ['ExchangeFluxes'].loc['M_ac', 'Eukaryotic']),
                             abs(Calibration_Hackett_C022['Simulation_Results_Euk']
                                 ['ExchangeFluxes'].loc['M_ac', 'Eukaryotic']),
                             abs(Calibration_Hackett_C03['Simulation_Results_Euk']['ExchangeFluxes'].loc['M_ac', 'Eukaryotic'])]

O2_Exchange_predicted_euk = [abs(Calibration_Hackett_C005['Simulation_Results_Euk']['ExchangeFluxes'].loc['M_o2', 'Eukaryotic']),
                             abs(Calibration_Hackett_C01['Simulation_Results_Euk']
                                 ['ExchangeFluxes'].loc['M_o2', 'Eukaryotic']),
                             abs(Calibration_Hackett_C016['Simulation_Results_Euk']
                                 ['ExchangeFluxes'].loc['M_o2', 'Eukaryotic']),
                             abs(Calibration_Hackett_C022['Simulation_Results_Euk']
                                 ['ExchangeFluxes'].loc['M_o2', 'Eukaryotic']),
                             abs(Calibration_Hackett_C03['Simulation_Results_Euk']['ExchangeFluxes'].loc['M_o2', 'Eukaryotic'])]

Glycerol_Exchange_predicted_euk = [abs(Calibration_Hackett_C005['Simulation_Results_Euk']['ExchangeFluxes'].loc['M_glyc', 'Eukaryotic']),
                                   abs(Calibration_Hackett_C01['Simulation_Results_Euk']
                                       ['ExchangeFluxes'].loc['M_glyc', 'Eukaryotic']),
                                   abs(Calibration_Hackett_C016['Simulation_Results_Euk']
                                       ['ExchangeFluxes'].loc['M_glyc', 'Eukaryotic']),
                                   abs(Calibration_Hackett_C022['Simulation_Results_Euk']
                                       ['ExchangeFluxes'].loc['M_glyc', 'Eukaryotic']),
                                   abs(Calibration_Hackett_C03['Simulation_Results_Euk']['ExchangeFluxes'].loc['M_glyc', 'Eukaryotic'])]

###
fig, axs = plt.subplots(2, 3, figsize=(28, 7), sharex=True)

# plt.figure()
axs[0, 0].plot(Mu_Hackett, Mu_Hackett, color='lightgreen')
axs[0, 0].scatter(Mu_Hackett, Mus_predicted, color='black')
axs[0, 0].scatter(Mu_Hackett, Mus_predicted_euk, color='red')
axs[0, 0].legend(['Hackett', 'Prok.', 'Euk.'])
axs[0, 0].set_title('Predicted vs measured growth-rate')
axs[0, 0].set_ylabel('$\mu$ [$h^{-1}$]')
axs[0, 0].set_xlabel('$\mu$ [$h^{-1}$]')
# plt.show()
# plt.savefig(pp, format='pdf')

# plt.figure()
axs[0, 1].plot(Mus_o2, Glc_J, color='lightblue')
axs[0, 1].plot(Mu_Hackett, Glc_Hackett, color='lightgreen')
axs[0, 1].scatter(Mus_predicted, Glc_Exchange_predicted, color='black', alpha=0.8)
axs[0, 1].scatter(Mus_predicted, Glc_Exchange_predicted_euk, color='red', alpha=0.8)
# plt.scatter(Mus_Euk,[abs(i) for i in FluxFile_Euk.loc['M_glc__D',conditions].values.tolist()],color='red',alpha=0.8)
# plt.scatter(Mus_testcase,[abs(i) for i in FluxFile_testcase.loc['M_glc__D',conditions].values.tolist()],color='orange',alpha=0.8)
# Effect of Specific Growth Rate on Fermentative Capacity of Baker’s Yeast#
axs[0, 1].legend(['van Hoek', 'Hackett', 'Prok.', 'Euk.'])
axs[0, 1].set_title('Glucose-uptake rate')
axs[0, 1].set_xlabel('$\mu$ [$h^{-1}$]')
axs[0, 1].set_ylabel('$J^{Ex}$ [$mmol * g^{-1}_{DW} * h^{-1}$]')
# plt.show()
# plt.savefig(pp, format='pdf')

# plt.figure()
axs[0, 2].plot(Mus_o2, O2_J, color='lightblue')
# plt.plot(Mu_Hackett,Glc_Hackett,color='lightgreen')
axs[0, 2].scatter(Mus_predicted, O2_Exchange_predicted, color='black', alpha=0.8)
axs[0, 2].scatter(Mus_predicted, O2_Exchange_predicted_euk, color='red', alpha=0.8)
# plt.scatter(Mus_Euk,[abs(i) for i in FluxFile_Euk.loc['M_glc__D',conditions].values.tolist()],color='red',alpha=0.8)
# plt.scatter(Mus_testcase,[abs(i) for i in FluxFile_testcase.loc['M_glc__D',conditions].values.tolist()],color='orange',alpha=0.8)
# Effect of Specific Growth Rate on Fermentative Capacity of Baker’s Yeast#
axs[0, 2].legend(['van Hoek', 'Prok.', 'Euk.'])
axs[0, 2].set_title('Oxygen-uptake rate')
axs[0, 2].set_xlabel('$\mu$ [$h^{-1}$]')
axs[0, 2].set_ylabel('$J^{Ex}$ [$mmol * g^{-1}_{DW} * h^{-1}$]')
# plt.show()
# plt.savefig(pp, format='pdf')

# plt.figure()
axs[1, 0].plot(Mus_o2, EtOH_J, color='lightblue')
axs[1, 0].plot(Mu_Hackett, EtOH_Hackett, color='lightgreen')
axs[1, 0].scatter(Mus_predicted, EtOH_Exchange_predicted, color='black', alpha=0.8)
axs[1, 0].scatter(Mus_predicted, EtOH_Exchange_predicted_euk, color='red', alpha=0.8)
# plt.scatter(Mus_Euk,[abs(i) for i in FluxFile_Euk.loc['M_glc__D',conditions].values.tolist()],color='red',alpha=0.8)
# plt.scatter(Mus_testcase,[abs(i) for i in FluxFile_testcase.loc['M_glc__D',conditions].values.tolist()],color='orange',alpha=0.8)
# Effect of Specific Growth Rate on Fermentative Capacity of Baker’s Yeast#
axs[1, 0].legend(['van Hoek', 'Hackett', 'Prok.', 'Euk.'])
axs[1, 0].set_title('Ethanol-excretion rate')
axs[1, 0].set_xlabel('$\mu$ [$h^{-1}$]')
axs[1, 0].set_ylabel('$J^{Ex}$ [$mmol * g^{-1}_{DW} * h^{-1}$]')
# plt.show()
# plt.savefig(pp, format='pdf')

# plt.figure()
axs[1, 1].plot(Mus_o2, Ac_J, color='lightblue')
axs[1, 1].plot(Mu_Hackett, Ac_Hackett, color='lightgreen')
axs[1, 1].scatter(Mus_predicted, Ac_Exchange_predicted, color='black', alpha=0.8)
axs[1, 1].scatter(Mus_predicted, Ac_Exchange_predicted_euk, color='red', alpha=0.8)
# plt.scatter(Mus_Euk,[abs(i) for i in FluxFile_Euk.loc['M_ac',conditions].values.tolist()],color='red',alpha=0.8)
# plt.scatter(Mus_testcase,[abs(i) for i in FluxFile_testcase.loc['M_ac',conditions].values.tolist()],color='orange',alpha=0.8)
# Effect of Specific Growth Rate on Fermentative Capacity of Baker’s Yeast#
axs[1, 1].legend(['van Hoek', 'Hackett', 'Prok.', 'Euk.'])
axs[1, 1].set_title('Acetate-excretion rate')
axs[1, 1].set_xlabel('$\mu$ [$h^{-1}$]')
axs[1, 1].set_ylabel('$J^{Ex}$ [$mmol * g^{-1}_{DW} * h^{-1}$]')
# plt.show()
# plt.savefig(pp, format='pdf')

axs[1, 2].plot(Mus_o2, Glyc_J, color='lightblue')
axs[1, 2].plot(Mu_Hackett, Glyc_Hackett, color='lightgreen')
axs[1, 2].scatter(Mus_predicted, Glycerol_Exchange_predicted, color='black', alpha=0.8)
axs[1, 2].scatter(Mus_predicted, Glycerol_Exchange_predicted_euk, color='red', alpha=0.8)
# plt.scatter(Mus_Euk,[abs(i) for i in FluxFile_Euk.loc['M_ac',conditions].values.tolist()],color='red',alpha=0.8)
# plt.scatter(Mus_testcase,[abs(i) for i in FluxFile_testcase.loc['M_ac',conditions].values.tolist()],color='orange',alpha=0.8)
# Effect of Specific Growth Rate on Fermentative Capacity of Baker’s Yeast#
axs[1, 2].legend(['van Hoek', 'Hackett', 'Prok.', 'Euk.'])
axs[1, 2].set_title('Glycerol-excretion rate')
axs[1, 2].set_xlabel('$\mu$ [$h^{-1}$]')
axs[1, 2].set_ylabel('$J^{Ex}$ [$mmol * g^{-1}_{DW} * h^{-1}$]')

plt.show()

protein_comparison_005 = pandas.DataFrame()
for i in list(set(list(Calibration_Hackett_C005['Simulation_Results']['ProtoProteins'].index)+list(Calibration_Hackett_C005['Proteome']['ID']))):
    protein_comparison_005.loc[i, 'ID'] = i
    if i in list(Calibration_Hackett_C005['Simulation_Results']['ProtoProteins'].index):
        protein_comparison_005.loc[i, 'Predicted'] = 6.023e20 * \
            Calibration_Hackett_C005['Simulation_Results']['ProtoProteins'].loc[i].values[0]
    if i in list(Calibration_Hackett_C005['Proteome']['ID']):
        protein_comparison_005.loc[i, 'Measured'] = 6.023e20 * \
            Calibration_Hackett_C005['Proteome'].loc[Calibration_Hackett_C005['Proteome']
                                                     ['ID'] == i, 'copy_number'].values[0]

protein_comparison_01 = pandas.DataFrame()
for i in list(set(list(Calibration_Hackett_C01['Simulation_Results']['ProtoProteins'].index)+list(Calibration_Hackett_C01['Proteome']['ID']))):
    protein_comparison_01.loc[i, 'ID'] = i
    if i in list(Calibration_Hackett_C01['Simulation_Results']['ProtoProteins'].index):
        protein_comparison_01.loc[i, 'Predicted'] = 6.023e20 * \
            Calibration_Hackett_C01['Simulation_Results']['ProtoProteins'].loc[i].values[0]
    if i in list(Calibration_Hackett_C01['Proteome']['ID']):
        protein_comparison_01.loc[i, 'Measured'] = 6.023e20 * \
            Calibration_Hackett_C01['Proteome'].loc[Calibration_Hackett_C01['Proteome']
                                                    ['ID'] == i, 'copy_number'].values[0]

protein_comparison_016 = pandas.DataFrame()
for i in list(set(list(Calibration_Hackett_C016['Simulation_Results']['ProtoProteins'].index)+list(Calibration_Hackett_C016['Proteome']['ID']))):
    protein_comparison_016.loc[i, 'ID'] = i
    if i in list(Calibration_Hackett_C016['Simulation_Results']['ProtoProteins'].index):
        protein_comparison_016.loc[i, 'Predicted'] = 6.023e20 * \
            Calibration_Hackett_C016['Simulation_Results']['ProtoProteins'].loc[i].values[0]
    if i in list(Calibration_Hackett_C016['Proteome']['ID']):
        protein_comparison_016.loc[i, 'Measured'] = 6.023e20 * \
            Calibration_Hackett_C016['Proteome'].loc[Calibration_Hackett_C016['Proteome']
                                                     ['ID'] == i, 'copy_number'].values[0]

protein_comparison_022 = pandas.DataFrame()
for i in list(set(list(Calibration_Hackett_C022['Simulation_Results']['ProtoProteins'].index)+list(Calibration_Hackett_C022['Proteome']['ID']))):
    protein_comparison_022.loc[i, 'ID'] = i
    if i in list(Calibration_Hackett_C022['Simulation_Results']['ProtoProteins'].index):
        protein_comparison_022.loc[i, 'Predicted'] = 6.023e20 * \
            Calibration_Hackett_C022['Simulation_Results']['ProtoProteins'].loc[i].values[0]
    if i in list(Calibration_Hackett_C022['Proteome']['ID']):
        protein_comparison_022.loc[i, 'Measured'] = 6.023e20 * \
            Calibration_Hackett_C022['Proteome'].loc[Calibration_Hackett_C022['Proteome']
                                                     ['ID'] == i, 'copy_number'].values[0]

protein_comparison_03 = pandas.DataFrame()
for i in list(set(list(Calibration_Hackett_C03['Simulation_Results']['ProtoProteins'].index)+list(Calibration_Hackett_C03['Proteome']['ID']))):
    protein_comparison_03.loc[i, 'ID'] = i
    if i in list(Calibration_Hackett_C03['Simulation_Results']['ProtoProteins'].index):
        protein_comparison_03.loc[i, 'Predicted'] = 6.023e20 * \
            Calibration_Hackett_C03['Simulation_Results']['ProtoProteins'].loc[i].values[0]
    if i in list(Calibration_Hackett_C03['Proteome']['ID']):
        protein_comparison_03.loc[i, 'Measured'] = 6.023e20 * \
            Calibration_Hackett_C03['Proteome'].loc[Calibration_Hackett_C03['Proteome']
                                                    ['ID'] == i, 'copy_number'].values[0]

fig, axs = plt.subplots(2, 3, figsize=(28, 7), sharex=True)

predidcted_proteins = protein_comparison_005.loc[(pandas.isna(protein_comparison_005['Predicted']) == False) & (
    pandas.isna(protein_comparison_005['Measured']) == False), 'Predicted']
quantified_proteins = protein_comparison_005.loc[(pandas.isna(protein_comparison_005['Predicted']) == False) & (
    pandas.isna(protein_comparison_005['Measured']) == False), 'Measured']
x_reg = numpy.reshape(numpy.array(list(predidcted_proteins)), (len(list(predidcted_proteins)), 1))
y_reg = numpy.reshape(numpy.array(list(quantified_proteins)), (len(list(quantified_proteins)), 1))
regressor = LinearRegression(fit_intercept=False)
regressor.fit(x_reg, y_reg)
predictions = regressor.predict(x_reg)
axs[0, 0].scatter(numpy.log10(protein_comparison_005['Predicted']),
                  numpy.log10(protein_comparison_005['Measured']))
axs[0, 0].plot([11, 17], [11, 17], color='green')
axs[0, 0].plot(numpy.log10(x_reg), numpy.log10(predictions), color='red', linewidth=2)
axs[0, 0].legend(['Identity', str(regressor.coef_), 'Data'])
axs[0, 0].set_title('Protein-Protein (Log10) 0.05')
axs[0, 0].set_xlabel('Predicted')
axs[0, 0].set_ylabel('Measured')

predidcted_proteins = protein_comparison_01.loc[(pandas.isna(protein_comparison_01['Predicted']) == False) & (
    pandas.isna(protein_comparison_01['Measured']) == False), 'Predicted']
quantified_proteins = protein_comparison_01.loc[(pandas.isna(protein_comparison_01['Predicted']) == False) & (
    pandas.isna(protein_comparison_01['Measured']) == False), 'Measured']
x_reg = numpy.reshape(numpy.array(list(predidcted_proteins)), (len(list(predidcted_proteins)), 1))
y_reg = numpy.reshape(numpy.array(list(quantified_proteins)), (len(list(quantified_proteins)), 1))
regressor = LinearRegression(fit_intercept=False)
regressor.fit(x_reg, y_reg)
predictions = regressor.predict(x_reg)
axs[0, 1].scatter(numpy.log10(protein_comparison_01['Predicted']),
                  numpy.log10(protein_comparison_01['Measured']))
axs[0, 1].plot([11, 17], [11, 17], color='green')
axs[0, 1].plot(numpy.log10(x_reg), numpy.log10(predictions), color='red', linewidth=2)
axs[0, 1].legend(['Identity', str(regressor.coef_), 'Data'])
axs[0, 1].set_title('Protein-Protein (Log10) 0.1')
axs[0, 1].set_xlabel('Predicted')
axs[0, 1].set_ylabel('Measured')

predidcted_proteins = protein_comparison_016.loc[(pandas.isna(protein_comparison_016['Predicted']) == False) & (
    pandas.isna(protein_comparison_016['Measured']) == False), 'Predicted']
quantified_proteins = protein_comparison_016.loc[(pandas.isna(protein_comparison_016['Predicted']) == False) & (
    pandas.isna(protein_comparison_016['Measured']) == False), 'Measured']
x_reg = numpy.reshape(numpy.array(list(predidcted_proteins)), (len(list(predidcted_proteins)), 1))
y_reg = numpy.reshape(numpy.array(list(quantified_proteins)), (len(list(quantified_proteins)), 1))
regressor = LinearRegression(fit_intercept=False)
regressor.fit(x_reg, y_reg)
predictions = regressor.predict(x_reg)
axs[0, 2].scatter(numpy.log10(protein_comparison_016['Predicted']),
                  numpy.log10(protein_comparison_016['Measured']))
axs[0, 2].plot([11, 17], [11, 17], color='green')
axs[0, 2].plot(numpy.log10(x_reg), numpy.log10(predictions), color='red', linewidth=2)
axs[0, 2].legend(['Identity', str(regressor.coef_), 'Data'])
axs[0, 2].set_title('Protein-Protein (Log10) 0.16')
axs[0, 2].set_xlabel('Predicted')
axs[0, 2].set_ylabel('Measured')

predidcted_proteins = protein_comparison_022.loc[(pandas.isna(protein_comparison_022['Predicted']) == False) & (
    pandas.isna(protein_comparison_022['Measured']) == False), 'Predicted']
quantified_proteins = protein_comparison_022.loc[(pandas.isna(protein_comparison_022['Predicted']) == False) & (
    pandas.isna(protein_comparison_022['Measured']) == False), 'Measured']
x_reg = numpy.reshape(numpy.array(list(predidcted_proteins)), (len(list(predidcted_proteins)), 1))
y_reg = numpy.reshape(numpy.array(list(quantified_proteins)), (len(list(quantified_proteins)), 1))
regressor = LinearRegression(fit_intercept=False)
regressor.fit(x_reg, y_reg)
predictions = regressor.predict(x_reg)
axs[1, 0].scatter(numpy.log10(protein_comparison_022['Predicted']),
                  numpy.log10(protein_comparison_022['Measured']))
axs[1, 0].plot([11, 17], [11, 17], color='green')
axs[1, 0].plot(numpy.log10(x_reg), numpy.log10(predictions), color='red', linewidth=2)
axs[1, 0].legend(['Identity', str(regressor.coef_), 'Data'])
axs[1, 0].set_title('Protein-Protein (Log10) 0.22')
axs[1, 0].set_xlabel('Predicted')
axs[1, 0].set_ylabel('Measured')

predidcted_proteins = protein_comparison_03.loc[(pandas.isna(protein_comparison_03['Predicted']) == False) & (
    pandas.isna(protein_comparison_03['Measured']) == False), 'Predicted']
quantified_proteins = protein_comparison_03.loc[(pandas.isna(protein_comparison_03['Predicted']) == False) & (
    pandas.isna(protein_comparison_03['Measured']) == False), 'Measured']
x_reg = numpy.reshape(numpy.array(list(predidcted_proteins)), (len(list(predidcted_proteins)), 1))
y_reg = numpy.reshape(numpy.array(list(quantified_proteins)), (len(list(quantified_proteins)), 1))
regressor = LinearRegression(fit_intercept=False)
regressor.fit(x_reg, y_reg)
predictions = regressor.predict(x_reg)
axs[1, 1].scatter(numpy.log10(protein_comparison_03['Predicted']),
                  numpy.log10(protein_comparison_03['Measured']))
axs[1, 1].plot([11, 17], [11, 17], color='green')
axs[1, 1].plot(numpy.log10(x_reg), numpy.log10(predictions), color='red', linewidth=2)
axs[1, 1].legend(['Identity', str(regressor.coef_), 'Data'])
axs[1, 1].set_title('Protein-Protein (Log10) 0.3')
axs[1, 1].set_xlabel('Predicted')
axs[1, 1].set_ylabel('Measured')

plt.show()

defKapps = pandas.DataFrame()
defKapps.loc['Hackett_C005', 'ID'] = 'Hackett_C005'
defKapps.loc['Hackett_C005', 'Default Kapp'] = Calibration_Hackett_C005['Default_Kapps'].iloc[-1, 2]
defKapps.loc['Hackett_C01', 'ID'] = 'Hackett_C01'
defKapps.loc['Hackett_C01', 'Default Kapp'] = Calibration_Hackett_C01['Default_Kapps'].iloc[-1, 2]
defKapps.loc['Hackett_C016', 'ID'] = 'Hackett_C016'
defKapps.loc['Hackett_C016', 'Default Kapp'] = Calibration_Hackett_C016['Default_Kapps'].iloc[-1, 2]
defKapps.loc['Hackett_C022', 'ID'] = 'Hackett_C022'
defKapps.loc['Hackett_C022', 'Default Kapp'] = Calibration_Hackett_C022['Default_Kapps'].iloc[-1, 2]
defKapps.loc['Hackett_C03', 'ID'] = 'Hackett_C03'
defKapps.loc['Hackett_C03', 'Default Kapp'] = Calibration_Hackett_C03['Default_Kapps'].iloc[-1, 2]
defKapps.to_csv('Default_Kapps_out.csv', sep=';', decimal=',')
