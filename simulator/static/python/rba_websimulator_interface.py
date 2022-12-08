import pandas
import math
import numpy
import os
import matplotlib.pyplot as plt
from rbatools.rba_session import SessionRBA

class WebsimulatorInterfaceRBA(object):
    # OK
    def __init__(self, xml_dir, simplified_parameter_changes=True):
        self.rba_session = SessionRBA(xml_dir=xml_dir,lp_solver='swiglpk')
        if simplified_parameter_changes:
            self.original_parameter_values = {}
            self.current_parameter_values = {}
            self.changeable_model_parameters = pandas.read_csv(
                str(xml_dir+'/changeable_parameters.csv'))
            for pm in list(self.changeable_model_parameters.loc[self.changeable_model_parameters['Type'] != 'medium_composition', 'Parameter']):
                self.add_parameter_multiplier(model_parameter=pm,rebuild_model=False)
                self.original_parameter_values[pm] = self.determine_parameter_values(model_parameter=pm, x_min=0.0, x_max=2, intervals=100)
                self.current_parameter_values[pm] = self.determine_parameter_values(model_parameter=pm, x_min=0.0, x_max=2, intervals=100)
            self.rba_session.rebuild_from_model()
        self.change_log = pandas.DataFrame(columns=['model_parameter', 'function_parameter', 'parameter_type', 'old_value', 'new_value', 'remark'])

    # OK
    def add_parameter_multiplier(self, model_parameter,rebuild_model=True):
        self.rba_session.add_parameter_multiplier(model_parameter=model_parameter,rebuild_model=rebuild_model)

    # OK
    def set_parameter_multiplier(self, model_parameter, parameter_type='', new_value=None, logging=True,rebuild_model=True):
        self.set_parameter_value(model_parameter=str(model_parameter+'_multiplier'), function_parameter='CONSTANT',
                                 parameter_type=parameter_type, change_message='multiplicative change', new_value=new_value, logging=logging,rebuild_model=rebuild_model)
        self.current_parameter_values[model_parameter] = self.determine_parameter_values(model_parameter=model_parameter, x_min=0.0, x_max=2, intervals=100)

    #  OK
    def set_parameter_value(self, model_parameter, function_parameter='CONSTANT', parameter_type='', change_message='', new_value=None, logging=True,rebuild_model=True):
        if new_value is not None:
            self.rba_session.set_parameter_value(model_parameter=model_parameter, function_parameter=function_parameter, parameter_type=parameter_type, new_value=new_value,rebuild_model=rebuild_model)
            if logging:
                change_df = pandas.DataFrame(columns=['model_parameter', 'function_parameter', 'parameter_type', 'old_value', 'new_value', 'remark'])
                change_df.loc[self.change_log.shape[0], 'model_parameter'] = model_parameter
                change_df.loc[self.change_log.shape[0],'function_parameter'] = function_parameter
                change_df.loc[self.change_log.shape[0], 'parameter_type'] = parameter_type
                old_value = self.rba_session.model.parameters.functions._elements_by_id[model_parameter].parameters._elements_by_id[function_parameter].value
                change_df.loc[self.change_log.shape[0], 'old_value'] = old_value
                change_df.loc[self.change_log.shape[0], 'new_value'] = new_value
                change_df.loc[self.change_log.shape[0], 'remark'] = change_message
                self.change_log.loc[self.change_log.shape[0],:] = change_df.loc[self.change_log.shape[0], :]

    # OK
    def set_medium_component(self, species=None, new_value=None, logging=True):
        if species is not None:
            if species in self.rba_session.Medium:
                if new_value is not None:
                    old_value = self.rba_session.Medium[species]
                    self.rba_session.set_medium({species: new_value})
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

    # OK
    def undo_last_change(self):
        if self.change_log.shape[0] > 0:
            old_value = self.change_log.iloc[-1, :]['old_value']
            model_parameter = self.change_log.iloc[-1, :]['model_parameter']
            function_parameter = self.change_log.iloc[-1, :]['function_parameter']
            parameter_type = self.change_log.iloc[-1, :]['parameter_type']
            remark = self.change_log.iloc[-1, :]['remark']
            if parameter_type != 'medium_composition':
                if remark == 'multiplicative change':
                    self.set_parameter_multiplier(model_parameter.split('_multiplier')[0], parameter_type='', new_value=old_value, logging=False,rebuild_model=True)
                else:
                    self.set_parameter_value(model_parameter=model_parameter, parameter_type=parameter_type,
                                             function_parameter=function_parameter, new_value=old_value, logging=False,rebuild_model=True)
            elif parameter_type == 'medium_composition':
                self.set_medium_component(species=model_parameter,
                                          new_value=old_value, logging=False)
            self.change_log = self.change_log.iloc[:-1, :]

    # OK
    def get_change_log(self):
        return(self.change_log)

    # OK
    def get_changeable_model_parameters(self):
        return(self.changeable_model_parameters[self.changeable_model_parameters['Type'] != 'medium_composition'])

    # OK
    def get_changeable_medium_species(self):
        return(self.changeable_model_parameters[self.changeable_model_parameters['Type'] == 'medium_composition'])

    # OK
    def determine_parameter_values(self, model_parameter, x_min=0, x_max=1, intervals=10):
        parameter_expression = self.rba_session.get_parameter_definition(model_parameter)
        if len(parameter_expression[model_parameter]["Variables"])==1:
            independent_variable=parameter_expression[model_parameter]["Variables"][0]
        elif len(parameter_expression[model_parameter]["Variables"])>1:
            independent_variable=[i for i in parameter_expression[model_parameter]["Variables"] if i!= "growth_rate"][0]
        else:
            independent_variable="growth_rate"
        independent_variable_values=list(numpy.linspace(x_min, x_max, intervals))
        out=self.rba_session.get_parameter_evolution(model_parameter=model_parameter,x_values={independent_variable:independent_variable_values})
        return(out)

    # OK
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
                        self.set_parameter_multiplier(model_parameter.split('_multiplier')[0], parameter_type='', new_value=new_value,logging=True,rebuild_model=False)
                    else:
                        self.set_parameter_value(model_parameter=model_parameter, parameter_type=parameter_type,function_parameter=function_parameter, new_value=new_value, logging=True,rebuild_model=False)
                elif parameter_type == 'medium_composition':
                    self.set_medium_component(species=model_parameter,new_value=new_value, logging=True)
            self.rba_session.rebuild_from_model()
        return 'ok'

    # OK
    def get_parameter_values_as_DataFrame(self, model_parameter):
        out = pandas.DataFrame(columns=[list(self.original_parameter_values[model_parameter].columns)[
                               0], 'Original values', 'Current values'])
        out[list(self.original_parameter_values[model_parameter].columns)[
            0]] = self.original_parameter_values[model_parameter][list(self.original_parameter_values[model_parameter].columns)[0]]
        out['Original values'] = self.original_parameter_values[model_parameter].iloc[:, 1]
        out['Current values'] = self.current_parameter_values[model_parameter].iloc[:, 1]
        return(out)

    # OK
    def plot_parameter_values(self, model_parameter):
        df = self.get_parameter_values_as_DataFrame(model_parameter)
        plt.plot(df[list(df.columns)[0]], df['Original values'], '--', color='k', alpha=0.8)
        plt.plot(df[list(df.columns)[0]], df['Current values'], '-', color='k', alpha=1)
        plt.legend(['Original', 'Current'])
        plt.title(model_parameter)
        plt.xlabel(list(df.columns)[0])
        plt.ylabel(model_parameter)
        plt.show()

    # OK
    def get_plot_values(self, model_parameter):
        df = self.get_parameter_values_as_DataFrame(model_parameter)
        return(df)
