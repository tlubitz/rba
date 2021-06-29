import numpy
from .rba_Session import RBA_Session
import pandas


def find_ribosomal_proteins(RBA_Session, model_processes=['TranslationC', 'TranslationM'], external_annotations=None):
    out = []
    for i in model_processes:
        out += [Simulation.ModelStructure.ProteinInfo.Elements[j]['ProtoID']
            for j in list(Simulation.ModelStructure.ProcessInfo.Elements[i]['Composition'].keys())]
    if external_annotations is not None:
        out += list(external_annotations['Ribosomal Proteins'])
    return(out)


def build_model_compartment_map(RBA_Session):
    out = {RBA_Session.ModelStructure.ProteinInfo.Elements[i]['ProtoID']: RBA_Session.ModelStructure.ProteinInfo.Elements[i]['Compartment'] for i in list(
        RBA_Session.ModelStructure.ProteinInfo.Elements.keys())}
    return(out)


def build_compartment_annotations(Compartment_Annotations_external, RBA_Session, model_protein_compartment_map):
    Compartment_Annotations_external.drop(columns=['Unnamed: 0'], inplace=True)
    Compartment_Annotations_internal = pandas.DataFrame()
    Compartment_Annotations_internal['ID'] = list(model_protein_compartment_map.keys())
    Compartment_Annotations_internal['ModelComp'] = list(model_protein_compartment_map.values())
    Compartment_Annotations = pandas.concat(
        [Compartment_Annotations_internal, Compartment_Annotations_external], axis=0)
    return(Compartment_Annotations)


def build_dataset_annotations(input, ID_column, Uniprot, Compartment_Annotations, model_protein_compartment_map):
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


def determine_apparent_process_efficiencies(growth_rate, input, rba_session, protein_data, condition, gene_id_col):
    process_efficiencies = pandas.DataFrame()
    for i in input.index:
        process_ID = input.loc[i, 'Process_ID']
        process_name = input.loc[i, 'Process_Name']
        process_client_compartments = input.loc[i, 'Process_Name'].split(' , ')
        constituting_proteins = {rba_session.ModelStructure.ProteinInfo.Elements[i]['ProtoID']: rba_session.ModelStructure.ProteinInfo.Elements[
            i]['AAnumber'] for i in rba_session.ModelStructure.ProcessInfo.Elements[process_name]['Composition'].keys()}
        Total_client_fraction = sum(
            [original_occupation_Lahtvee.loc[i, 'new_protein_fraction'] for i in process_client_compartments])
        n_AAs_in_machinery = 0
        machinery_size = 0
        for i in constituting_proteins.keys():
            if i in protein_data[gene_id_col]:
                protein_data.loc[protein_data[gene_id_col] == i, ]
                n_AAs_in_machinery += factor_A_Lahtvee*factor_B_Lahtvee * \
                    protein_data.loc[protein_data[gene_id_col] == i, condition] * \
                    protein_data.loc[protein_data[gene_id_col] == i, 'AA_residues']
                machinery_size += constituting_proteins[i]
        # right reference amounth?
        relative_Protein_fraction_of_machinery = n_AAs_in_machinery / \
            original_occupation_Lahtvee.loc['Total', 'original_amino_acid_occupation']
        specific_capacity = growth_rate*Total_client_fraction/relative_Protein_fraction_of_machinery
        apparent_capacity = specific_capacity*machinery_size
        process_ID[process_name] = apparent_capacity
        process_efficiencies.loc[process_name, 'Process'] = process_ID
        process_efficiencies.loc[process_name, 'Parameter'] = str(process_ID+'_apparent_efficiency')
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
    return({'Summary': out, 'Correction_factors': {'A': factor_A, 'B': factor_B}})


def build_input_for_default_kapp_estimation(input):
    out = pandas.DataFrame(columns=['Compartment_ID', 'Density', 'PG_fraction'])
    for i in input['Summary'].index:
        if i not in ['Total', 'Ribosomes']:
            out.loc[i, 'Compartment_ID'] = i
            out.loc[i, 'Density'] = input['Summary'].loc[i, 'new_protein_fraction']
            out.loc[i, 'PG_fraction'] = input['Summary'].loc[i, 'new_PG_fraction']
    return(out)


def calibration_workflow(input,
                         condition,
                         gene_ID_column,
                         imposed_compartment_fractions={
                             'Ribosomes': 0.117, 'Secreted': 0.047, 'n': 0.108},
                         growth_rate,
                         process_efficiency_estimation_input,
                         flux_bounds,
                         media_composition,
                         rba_session):
    correction_results = correction_pipeline(input=input,
                                             condition=condition,
                                             compartments_to_ignore=['DEF', 'DEFA', 'Def'],
                                             compartments_no_original_PG=['n', 'Secreted'],
                                             fractions_entirely_replaced_with_expected_value=[
                                                 'Ribosomes'],
                                             imposed_compartment_fractions=imposed_compartment_fractions,
                                             directly_corrected_compartments=[
                                                 'c', 'cM', 'erM', 'gM', 'm', 'mIM', 'mIMS', 'mOM', 'vM', 'x'],
                                             merged_compartments={'c': 'Ribosomes'})
    process_efficiencies = determine_apparent_process_efficiencies(growth_rate=growth_rate,
                                                                   input=process_efficiency_estimation_input,
                                                                   rba_session=rba_session,
                                                                   protein_data=input,
                                                                   condition=condition,
                                                                   gene_id_col=gene_ID_column)
    rba_Session_copy = copy.deepcopy(rba_session)
    rba_Session_copy.inject_process_capacities(process_efficiencies=process_efficiencies)
    # rba_Session_copy.setMedium(media_composition...)
    Default_Kapps = rba_Session_copy.estimate_default_Kapps(target_mu=growth_rate,
                                                            compartment_densities_and_PGs=build_input_for_default_kapp_estimation(
                                                                correction_results),
                                                            flux_bounds=flux_bounds)
    # Specific Kapps
    ## Correct Proteome beforehand !!!! ##
    Specific_Kapps = rba_Session_copy.estimate_specific_Kapps(proteomicsData=proteome,
                                                              flux_bounds=flux_bounds,
                                                              mu=growth_rate],
                                                              biomass_function = None,
                                                              target_biomass_function = True)
    rba_Session_copy.inject_specific_kapps(specific_kapps = Specific_Kapps, round_to_digits = 0)
    # rba_Session_copy.inject_default_kapps(default_kapp, default_transporter_kapp)
    # Test predictions
        # Given medium predict Mu, Exchanges and Proteome
            # Prokaryotic
            # Eukaryotic

# Think about elegant input of Parameters


# 1. import model and uniprot-file and compartment-annotation
dirInput='Yeast_iMM904_RBA_model'
Simulation=RBA_Session(dirInput)
Uniprot=pandas.read_csv('Yeast_iMM904_RBA_model/uniprot.csv', sep = '\t')

Compartment_Annotations_external=pandas.read_csv(
    'DataSetsYeastRBACalibration/Manually_curated_not_in_model_Protein_Locations.csv', index_col = None)


Compartment_Annotations=build_compartment_annotations(
    Compartment_Annotations_external = Compartment_Annotations_external, RBA_Session = Simulation)


ribosomal_proteins=find_ribosomal_proteins(RBA_Session = Simulation, model_processes = [
                                           'TranslationC', 'TranslationM'], external_annotations = None)

# 2. import Hackett fold-changes
Hackett_Clim_FCs=pandas.read_csv('DataSetsYeastRBACalibration/Hacket_Clim_ProteinFCs.csv')

########################################
#### Bootstrapping-loop starts here ####

# 3. import lahtvee copy-Numbers
#     Correct to bring reference to gram dry-weight (originally per picogram dry-weight)
#     Add noise to Lahtvee copy-Numbers

Lahtvee_REF=pandas.read_csv('DataSetsYeastRBACalibration/LahtveeRefProteomicsData.csv')
picogram_togram_coefficient=1e12
Lahtvee_REF['Lahtvee_REF'] *= picogram_togram_coefficient
Lahtvee_REF=Lahtvee_REF.loc[pandas.isna(Lahtvee_REF['Lahtvee_REF']) == False]

# 4.Add location-annotation and model-function annotiation to proteomes (should not be done every iteration...).

annotations_Lahtvee=build_dataset_annotations(input = Lahtvee_REF, ID_column = 'Gene ID', Uniprot = Uniprot,
                                                Compartment_Annotations = Compartment_Annotations, model_protein_compartment_map = model_protein_compartment_map)
annotations_Hackett=build_dataset_annotations(input = Hackett_Clim_FCs, ID_column = 'Gene', Uniprot = Uniprot,
                                                Compartment_Annotations = Compartment_Annotations, model_protein_compartment_map = model_protein_compartment_map)


full_annotations=pandas.concat([annotations_Lahtvee, annotations_Hackett], axis = 0)
index=full_annotations.index
is_duplicate=index.duplicated(keep = "first")
not_duplicate=~is_duplicate
full_annotations=full_annotations[not_duplicate]


# 5. restore absolute reference by mapping Lahtvee-REF condition onto Clim-0.1 condition and using the respective fold-changes.
# 6. restore all condition-proteomes with restore absolute reference.
Hackett_Clim_AbsoluteCopies_restored=infer_copy_numbers_from_reference_copy_numbers(
    fold_changes = Hackett_Clim_FCs, absolute_data = Lahtvee_REF, matching_column_in_fold_change_data = 'C0.11', matching_column_in_absolute_data = 'REF', conditions_in_fold_change_data_to_restore = ['C0.05', 'C0.11', 'C0.16', 'C0.22', 'C0.30'])

# 7. Add annotation-information to proteomes:


# Lahtvee_REF = add_annotations_to_proteome(input=Lahtvee_REF, ID_column='Gene ID', annotations=annotations)
# calibration_Lahtvee=calibration_workflow(input=Lahtvee_REF,condition='REF',gene_ID_column='Gene ID',growth_rate=0.1,process_efficiency_estimation_input,rba_session=Simulation)

Hackett_Clim_AbsoluteCopies_restored=add_annotations_to_proteome(
    input = Hackett_Clim_AbsoluteCopies_restored, ID_column = 'ID', annotations = full_annotations)
Lahtvee_REF=add_annotations_to_proteome(
    input = Lahtvee_REF, ID_column = 'Gene ID', annotations = full_annotations)

Calibration_Lahtvee=calibration_workflow(input, condition, gene_ID_column, growth_rate, process_efficiency_estimation_input, rba_session):


# 8. Determine compartment-occupation and respective PG-fraction among data-set
# 9. Determine compartment-specific copy-number correction-factor, protein- and PG-fraction, using expected protein-fractions as input.

correction_Lahtvee=correction_pipeline(input = Lahtvee_REF, condition = 'REF', compartments_to_ignore = ['DEF', 'DEFA', 'Def'], compartments_no_original_PG = ['n', 'Secreted'], fractions_entirely_replaced_with_expected_value = [
                                         'Ribosomes'], imposed_compartment_fractions = {'Ribosomes': 0.117, 'Secreted': 0.047, 'n': 0.108}, directly_corrected_compartments = ['c', 'cM', 'erM', 'gM', 'm', 'mIM', 'mIMS', 'mOM', 'vM', 'x'], merged_compartments = {'c': 'Ribosomes'})
correction_Hackett_005=correction_pipeline(input = Hackett_Clim_AbsoluteCopies_restored, condition = 'C0.05', compartments_to_ignore = ['DEF', 'DEFA', 'Def'], compartments_no_original_PG = ['n', 'Secreted'], fractions_entirely_replaced_with_expected_value = [
    'Ribosomes'], imposed_compartment_fractions = {'Ribosomes': 0.098, 'Secreted': 0.042, 'n': 0.1045}, directly_corrected_compartments = ['c', 'cM', 'erM', 'gM', 'm', 'mIM', 'mIMS', 'mOM', 'vM', 'x'], merged_compartments = {'c': 'Ribosomes'})
correction_Hackett_01=correction_pipeline(input = Hackett_Clim_AbsoluteCopies_restored, condition = 'C0.11', compartments_to_ignore = ['DEF', 'DEFA', 'Def'], compartments_no_original_PG = ['n', 'Secreted'], fractions_entirely_replaced_with_expected_value = [
    'Ribosomes'], imposed_compartment_fractions = {'Ribosomes': 0.117, 'Secreted': 0.047, 'n': 0.108}, directly_corrected_compartments = ['c', 'cM', 'erM', 'gM', 'm', 'mIM', 'mIMS', 'mOM', 'vM', 'x'], merged_compartments = {'c': 'Ribosomes'})
correction_Hackett_016=correction_pipeline(input = Hackett_Clim_AbsoluteCopies_restored, condition = 'C0.16', compartments_to_ignore = ['DEF', 'DEFA', 'Def'], compartments_no_original_PG = ['n', 'Secreted'], fractions_entirely_replaced_with_expected_value = [
    'Ribosomes'], imposed_compartment_fractions = {'Ribosomes': 0.134, 'Secreted': 0.05, 'n': 0.1115}, directly_corrected_compartments = ['c', 'cM', 'erM', 'gM', 'm', 'mIM', 'mIMS', 'mOM', 'vM', 'x'], merged_compartments = {'c': 'Ribosomes'})
correction_Hackett_022=correction_pipeline(input = Hackett_Clim_AbsoluteCopies_restored, condition = 'C0.22', compartments_to_ignore = ['DEF', 'DEFA', 'Def'], compartments_no_original_PG = ['n', 'Secreted'], fractions_entirely_replaced_with_expected_value = [
    'Ribosomes'], imposed_compartment_fractions = {'Ribosomes': 0.155, 'Secreted': 0.055, 'n': 0.115}, directly_corrected_compartments = ['c', 'cM', 'erM', 'gM', 'm', 'mIM', 'mIMS', 'mOM', 'vM', 'x'], merged_compartments = {'c': 'Ribosomes'})
correction_Hackett_03=correction_pipeline(input = Hackett_Clim_AbsoluteCopies_restored, condition = 'C0.30', compartments_to_ignore = ['DEF', 'DEFA', 'Def'], compartments_no_original_PG = ['n', 'Secreted'], fractions_entirely_replaced_with_expected_value = [
    'Ribosomes'], imposed_compartment_fractions = {'Ribosomes': 0.183, 'Secreted': 0.061, 'n': 0.1195}, directly_corrected_compartments = ['c', 'cM', 'erM', 'gM', 'm', 'mIM', 'mIMS', 'mOM', 'vM', 'x'], merged_compartments = {'c': 'Ribosomes'})

# 10 Determine correction-factor C:
# 11. Correct Copy-numbers respectively
# 12 Determine apparent process-capacities:
# Input: Process-ID , Process name and client-compartments and growth-rate.
process_efficiencies_Lahtvee=determine_apparent_process_efficiencies(
    growth_rate = 0.1, input = None, rba_session = Simulation, protein_data = Lahtvee_REF, condition = 'REF', gene_id_col = 'Gene ID')
# 13. Inject default process-capacities
Simulation.inject_process_capacities(process_efficiencies = process_efficiencies_Lahtvee)
# 14. Run default-kapp estimation
Default_Kapps=Simulation.estimate_default_Kapps(
    target_mu = growth_rate, compartment_densities_and_PGs = build_input_for_default_kapp_estimation(correction_Lahtvee['Summary']), flux_bounds = None)

# 15. Run specific-kapp estimation

# 16. Make RBA-prediction with injected process-efficiencies, default- and specific-kapps
#     Store Mu, Proteome and Fluxes

###############
# 1. Sort out model-version stuff and make calibration compatible with normal model.
#    1.1 Test

#    total AAs as protein target -->
#        Aggregate:
#        <functionReference function="amino_acid_concentrationHackett"/>
#        <functionReference function="inverse_average_protein_length"/>
#    total RNA as rna target -->
#        Aggregate:
#        <functionReference function="RNA_massfraction_CarbonLimitation"/>
#        <functionReference function="RNA_inversemillimolarweight"/>

# 2. Redo manual protein-locations annotation-file
