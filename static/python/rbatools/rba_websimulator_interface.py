import pandas
import rba
import math
import numpy
import os
import matplotlib.pyplot as plt
from .rba_Session import RBA_Session


class RBA_websimulator_interface(object):

    # ok
    def __init__(self, xml_dir, simplified_parameter_changes=True):
        self.rba_session = RBA_Session(xml_dir)
        if simplified_parameter_changes:
            self.original_parameter_values = {}
            self.current_parameter_values = {}
            self.changeable_model_parameters = pandas.read_csv(
                str(xml_dir+'/changeable_parameters.csv'))
            for pm in list(self.changeable_model_parameters.loc[self.changeable_model_parameters['Type'] != 'medium_composition', 'Parameter']):
                self.add_parameter_multiplier(model_parameter=pm)
                self.original_parameter_values[pm] = self.determine_parameter_values(
                    model_parameter=pm, x_min=0.0, x_max=2, intervals=100)
                self.current_parameter_values[pm] = self.determine_parameter_values(
                    model_parameter=pm, x_min=0.0, x_max=2, intervals=100)
        self.change_log = pandas.DataFrame(
            columns=['model_parameter', 'function_parameter', 'parameter_type', 'old_value', 'new_value', 'remark'])

    # ok
    def add_parameter_multiplier(self, model_parameter):
        if model_parameter in self.rba_session.model.parameters.functions._elements_by_id.keys():
            self.rba_session.model.parameters.functions._elements_by_id[model_parameter].id = str(
                model_parameter+'_original_definition')
            self.rba_session.model.parameters.functions._elements_by_id[str(
                model_parameter+'_original_definition')] = self.rba_session.model.parameters.functions._elements_by_id.pop(model_parameter)
            self.rba_session.model.parameters.functions.append(rba.xml.parameters.Function(
                str(model_parameter+'_multiplier'), 'constant', parameters={'CONSTANT': 1.0}, variable=None))
            self.rba_session.model.parameters.functions._elements = list(
                self.rba_session.model.parameters.functions._elements_by_id.values())
            self.rba_session.model.parameters.aggregates.append(
                rba.xml.parameters.Aggregate(model_parameter, 'multiplication'))
            self.rba_session.model.parameters.aggregates._elements_by_id[model_parameter].function_references.append(
                rba.xml.parameters.FunctionReference(str(model_parameter+'_original_definition')))
            self.rba_session.model.parameters.aggregates._elements_by_id[model_parameter].function_references.append(
                rba.xml.parameters.FunctionReference(str(model_parameter+'_multiplier')))
        elif model_parameter in self.rba_session.model.parameters.aggregates._elements_by_id.keys():
            self.rba_session.model.parameters.functions.append(rba.xml.parameters.Function(
                str(model_parameter+'_multiplier'), 'constant', parameters={'CONSTANT': 1.0}, variable=None))
            self.rba_session.model.parameters.aggregates._elements_by_id[model_parameter].function_references.append(
                rba.xml.parameters.FunctionReference(str(model_parameter+'_multiplier')))
        self.rba_session.rebuild_from_model()

    # ok
    def set_parameter_multiplier(self, model_parameter, parameter_type='', new_value=None, logging=True):
        self.set_parameter_value(model_parameter=str(model_parameter+'_multiplier'), function_parameter='CONSTANT',
                                 parameter_type=parameter_type, change_message='multiplicative change', new_value=new_value, logging=logging)
        self.current_parameter_values[model_parameter] = self.determine_parameter_values(
            model_parameter=model_parameter, x_min=0.0, x_max=2, intervals=100)

    # ok
    def set_parameter_value(self, model_parameter, function_parameter, parameter_type='', change_message='', new_value=None, logging=True):
        if new_value is not None:
            if model_parameter in self.rba_session.model.parameters.functions._elements_by_id.keys():
                old_value = self.rba_session.model.parameters.functions._elements_by_id[
                    model_parameter].parameters._elements_by_id[function_parameter].value
                self.rba_session.model.parameters.functions._elements_by_id[model_parameter]
                self.rba_session.model.parameters.functions._elements_by_id[
                    model_parameter].parameters._elements_by_id[function_parameter].value = new_value
                self.rba_session.rebuild_from_model()
                if logging:
                    change_df = pandas.DataFrame(
                        columns=['model_parameter', 'function_parameter', 'parameter_type', 'old_value', 'new_value', 'remark'])
                    change_df.loc[self.change_log.shape[0], 'model_parameter'] = model_parameter
                    change_df.loc[self.change_log.shape[0],
                                  'function_parameter'] = function_parameter
                    change_df.loc[self.change_log.shape[0], 'parameter_type'] = parameter_type
                    change_df.loc[self.change_log.shape[0], 'old_value'] = old_value
                    change_df.loc[self.change_log.shape[0], 'new_value'] = new_value
                    change_df.loc[self.change_log.shape[0], 'remark'] = change_message
                    self.change_log.loc[self.change_log.shape[0],
                                        :] = change_df.loc[self.change_log.shape[0], :]

    # ok
    def set_medium_component(self, species=None, new_value=None, logging=True):
        if species is not None:
            if species in self.rba_session.Medium:
                if new_value is not None:
                    old_value = self.rba_session.Medium[species]
                    self.rba_session.setMedium({species: new_value})
                    if logging:
                        change_df = pandas.DataFrame(
                            columns=['model_parameter', 'function_parameter', 'parameter_type', 'old_value', 'new_value', 'remark'])
                        change_df.loc[self.change_log.shape[0], 'model_parameter'] = species
                        change_df.loc[self.change_log.shape[0], 'function_parameter'] = ''
                        change_df.loc[self.change_log.shape[0],
                                      'parameter_type'] = 'medium_composition'
                        change_df.loc[self.change_log.shape[0], 'old_value'] = old_value
                        change_df.loc[self.change_log.shape[0], 'new_value'] = new_value
                        change_df.loc[self.change_log.shape[0], 'remark'] = ''
                        self.change_log.loc[self.change_log.shape[0],
                                            :] = change_df.loc[self.change_log.shape[0], :]

    # ok
    def undo_last_change(self):
        if self.change_log.shape[0] > 0:
            old_value = self.change_log.iloc[-1, :]['old_value']
            model_parameter = self.change_log.iloc[-1, :]['model_parameter']
            function_parameter = self.change_log.iloc[-1, :]['function_parameter']
            parameter_type = self.change_log.iloc[-1, :]['parameter_type']
            remark = self.change_log.iloc[-1, :]['remark']
            if parameter_type != 'medium_composition':
                if remark == 'multiplicative change':
                    self.set_parameter_multiplier(model_parameter.split('_multiplier')[
                                                  0], parameter_type='', new_value=old_value, logging=False)

                else:
                    self.set_parameter_value(model_parameter=model_parameter, parameter_type=parameter_type,
                                             function_parameter=function_parameter, new_value=old_value, logging=False)
            elif parameter_type == 'medium_composition':
                self.set_medium_component(species=model_parameter,
                                          new_value=old_value, logging=False)
            self.change_log = self.change_log.iloc[:-1, :]

    # ok
    def get_change_log(self):
        return(self.change_log)

    # ok
    def get_changeable_model_parameters(self):
        return(self.changeable_model_parameters[self.changeable_model_parameters['Type'] != 'medium_composition'])

    # ok
    def get_changeable_medium_species(self):
        return(self.changeable_model_parameters[self.changeable_model_parameters['Type'] == 'medium_composition'])

    # ok
    def determine_parameter_values(self, model_parameter, x_min=0, x_max=None, intervals=10):
        parameter_expression = self.rba_session.get_parameter_definition(model_parameter)
        if list(parameter_expression.values())[0]['Type'] is 'Aggregate':
            if x_max is None:
                df_list = {term: evaluate_expression(expression_dictionary=self.rba_session.get_parameter_definition(term), independent_variable=list(
                    numpy.linspace(x_min, 1, intervals))) for term in list(parameter_expression.values())[0]['Multiplicative Terms']}
            else:
                df_list = {term: evaluate_expression(expression_dictionary=self.rba_session.get_parameter_definition(term), independent_variable=list(
                    numpy.linspace(x_min, x_max, intervals))) for term in list(parameter_expression.values())[0]['Multiplicative Terms']}
            out = pandas.DataFrame(columns=[list(parameter_expression.values())[
                                   0]['Variables'][0], list(parameter_expression.keys())[0]])
            out[list(parameter_expression.values())[0]['Variables'][0]] = df_list[list(
                df_list.keys())[0]][list(parameter_expression.values())[0]['Variables'][0]]
            out[list(parameter_expression.keys())[0]] = [1]*out.shape[0]
            for term in list(df_list.keys()):
                out[list(parameter_expression.keys())[0]] *= df_list[term][term]
        else:
            if x_max is None:
                out = evaluate_expression(expression_dictionary=parameter_expression, independent_variable=list(
                    numpy.linspace(x_min, 1, intervals)))
            else:
                out = evaluate_expression(expression_dictionary=parameter_expression, independent_variable=list(
                    numpy.linspace(x_min, x_max, intervals)))
        return(out)

    # ok
    def replay_from_logfile(self, file_path=''):
        if os.path.isfile(file_path):
            previous_log = pandas.read_csv(file_path)
            for i in range(previous_log.shape[0]):
                new_value = previous_log.iloc[i, :]['new_value']
                model_parameter = previous_log.iloc[i, :]['model_parameter']
                function_parameter = previous_log.iloc[i, :]['function_parameter']
                parameter_type = previous_log.iloc[i, :]['parameter_type']
                remark = previous_log.iloc[i, :]['remark']
                if parameter_type != 'medium_composition':
                    if remark == 'multiplicative change':
                        self.set_parameter_multiplier(model_parameter.split('_multiplier')[
                                                      0], parameter_type='', new_value=new_value, logging=True)

                    else:
                        self.set_parameter_value(model_parameter=model_parameter, parameter_type=parameter_type,
                                                 function_parameter=function_parameter, new_value=new_value, logging=True)

                elif parameter_type == 'medium_composition':
                    self.set_medium_component(species=model_parameter,
                                              new_value=new_value, logging=True)

    def get_parameter_values_as_DataFrame(self, model_parameter):
        out = pandas.DataFrame(columns=[list(self.original_parameter_values[model_parameter].columns)[
                               0], 'Original values', 'Current values'])
        out[list(self.original_parameter_values[model_parameter].columns)[
            0]] = self.original_parameter_values[model_parameter][list(self.original_parameter_values[model_parameter].columns)[0]]
        out['Original values'] = self.original_parameter_values[model_parameter].iloc[:, 1]
        out['Current values'] = self.current_parameter_values[model_parameter].iloc[:, 1]
        return(out)

    def plot_parameter_values(self, model_parameter):
        df = self.get_parameter_values_as_DataFrame(model_parameter)
        plt.plot(df[list(df.columns)[0]], df['Original values'])
        plt.plot(df[list(df.columns)[0]], df['Current values'])
        plt.legend(['Original', 'Current'])
        plt.title(model_parameter)
        plt.xlabel(list(df.columns)[0])
        plt.ylabel(model_parameter)
        plt.show()


def evaluate_expression(expression_dictionary, independent_variable):
    if type(independent_variable) is list:
        x_vals = independent_variable
        x_vals.sort()
    else:
        x_vals = [independent_variable]
    if list(expression_dictionary.values())[0]['Type'] == 'constant':
        out = pandas.DataFrame(columns=[list(expression_dictionary.values())[
                               0]['Variables'][0], list(expression_dictionary.keys())[0]])
        out[list(expression_dictionary.values())[0]['Variables'][0]] = x_vals
        out[list(expression_dictionary.keys())[0]] = [list(expression_dictionary.values())
                                                      [0]['Function_parameters']['CONSTANT']]*len(x_vals)
    elif list(expression_dictionary.values())[0]['Type'] == 'exponential':
        parameter_dict = list(expression_dictionary.values())[0]['Function_parameters']
        parameter_dict['e'] = math.e
        out = pandas.DataFrame(columns=[list(expression_dictionary.values())[
                               0]['Variables'][0], list(expression_dictionary.keys())[0]])
        out[list(expression_dictionary.values())[0]['Variables'][0]] = x_vals
        out[list(expression_dictionary.keys())[0]] = [eval(list(expression_dictionary.values())[0]['Equation'], {
            list(expression_dictionary.values())[0]['Variables'][0]:x_val}, parameter_dict) for x_val in x_vals]
    elif list(expression_dictionary.values())[0]['Type'] == 'linear':
        parameter_dict = list(expression_dictionary.values())[0]['Function_parameters']
        out = pandas.DataFrame(columns=[list(expression_dictionary.values())[
                               0]['Variables'][0], list(expression_dictionary.keys())[0]])
        out[list(expression_dictionary.values())[0]['Variables'][0]] = x_vals
        effective_x_vals = [check_var_bounds(v=x, v_min=list(expression_dictionary.values())[0]['Function_parameters']['X_MIN'], v_max=list(
            expression_dictionary.values())[0]['Function_parameters']['X_MAX']) for x in x_vals]
        equation_values = [eval(list(expression_dictionary.values())[0]['Equation'], {list(
            expression_dictionary.values())[0]['Variables'][0]:x_val}, parameter_dict) for x_val in effective_x_vals]
        out[list(expression_dictionary.keys())[0]] = [check_var_bounds(v=ev, v_min=list(expression_dictionary.values())[
            0]['Function_parameters']['Y_MIN'], v_max=list(expression_dictionary.values())[0]['Function_parameters']['Y_MAX']) for ev in equation_values]

    elif list(expression_dictionary.values())[0]['Type'] == 'michaelisMenten':
        parameter_dict = list(expression_dictionary.values())[0]['Function_parameters']
        out = pandas.DataFrame(columns=[list(expression_dictionary.values())[
            0]['Variables'][0], list(expression_dictionary.keys())[0]])
        out[list(expression_dictionary.values())[0]['Variables'][0]] = x_vals
        equation_values = [eval(list(expression_dictionary.values())[0]['Equation'], {list(
            expression_dictionary.values())[0]['Variables'][0]:x_val}, parameter_dict) for x_val in x_vals]
        out[list(expression_dictionary.keys())[0]] = [check_var_bounds(v=ev, v_min=list(expression_dictionary.values())[
            0]['Function_parameters']['Y_MIN'], v_max=numpy.inf) for ev in equation_values]
    return(out)


def check_var_bounds(v, v_min, v_max):
    if v > v_max:
        return(v_max)
    if v < v_min:
        return(v_min)
    else:
        return(v)