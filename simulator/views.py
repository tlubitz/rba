from django.http import HttpResponseRedirect
from django.shortcuts import render, render_to_response
from simulator.forms import UploadFileForm
from django.conf import settings
from django.http import HttpResponse
from django.views.decorators.http import require_POST
from django.shortcuts import redirect
from django.views.decorators.csrf import csrf_exempt

import os
import shutil
import zipfile
import json
import socket
#import cplex
import random
import matplotlib
#matplotlib.use('TkAgg')
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .static.python.RBA_tools_WIP.rbatools import rba_websimulator_interface


def index(request):
    '''
    main page handling rba file upload and fine tuning parameters
    '''
    # create essential variables
    for var in ['rbafilezip',
                'rbafilename',
                'emap_path',
                'proteomap_path',
                'sbtab_path',
                'dl_path',
                'log_path',
                'status',
                'errors',
                'plot_path']:
        if not request.session.get(var, None):
            request.session[var] = False

    if not request.session.get('first_sim', None):
        request.session['first_sim'] = True

    if not request.session.get('session_id', None):
        request.session['session_id'] = random.randint(0,1000)

    for var in ['error_code', 'csv_paths', 'model_species_list', 'model_parameters_list']:
        if not request.session.get(var, None):
            request.session[var] = []

    if request.method == 'POST':
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            if not (request.FILES['zipfile'].name).endswith('.zip'):
                request.session['error_code'].append('This is not a zip file. Please only use zip files.')
            else:
                # first create new directory
                media_directory = settings.MEDIA_ROOT
                current_rba_name = request.FILES['zipfile'].name[:-4]
                request.session['newdir'] = os.path.join(media_directory, current_rba_name)

                try: os.mkdir(request.session['newdir'])
                except: print('Could not create new directory for zip file contents.')

                # then unzip the zip file
                try:
                    newzip = zipfile.ZipFile(request.FILES['zipfile'], 'r')
                    newzip.extractall(request.session['newdir'])
                    request.session['rbafilename'] = request.FILES['zipfile'].name
                except: request.session['error_code'].append('Could not unzip file. Is it valid?')
                
                return HttpResponseRedirect('/simulator')
    else:
        form = UploadFileForm()

    if request.session['error_code'] != []: request.session['errors'] = True
    else:  request.session['errors'] = False
    request.session.modified = True

    return render(request, 'index.html', {'form': form,
                                          'status': request.session['status'],
                                          'errors': request.session['errors'],
                                          'rbafilename': request.session['rbafilename'],
                                          'error_code': request.session['error_code'],
                                          'emap_path': request.session['emap_path'],
                                          'proteomap_path': request.session['proteomap_path'],
                                          'csv_paths': request.session['csv_paths'],
                                          'log_path': request.session['log_path'],
                                          'sbtab_path': request.session['sbtab_path'],
                                          'dl_path': request.session['dl_path'],
                                          'plot_path': request.session['plot_path'],
                                          'model_parameters_list': request.session['model_parameters_list'],
                                          'model_species_list': request.session['model_species_list']})


def clearsession(request):
    '''
    clears all session variables
    '''
    # delete current project directory (only if it was an uploaded model)
    if not request.session['newdir'].startswith('simulator/static/python/models/') and not request.session['newdir'].startswith(os.getcwd()):
        try: shutil.rmtree(request.session['newdir'])
        except: print('Cannot delete %s' %request.session['newdir'])

    # delete session variables
    keys = ['rbafilezip', 'rbafilename', 'newdir', 'error_code', 'emap_path', 'proteomap_path', 'csv_paths', 'sbtab_path', 'dl_path', 'plot_path', 'status', 'errors', 'model_parameters_list', 'model_species_list']
    for key in keys:
        try: del request.session[key]
        except: print('Cannot delete %s' %(key))

    request.session['first_sim'] = 'Yes'
    request.session['session_id'] = None
    request.session.modified = True

    return HttpResponse('ok')


@csrf_exempt
def loadmodel(request):
    '''
    load an existing model from server
    '''
    if socket.gethostname() == 'timputer':
        pre_path = 'simulator/static/results/'
        mode = 'dev'
    else:
        pre_path = os.getcwd() + '/rba/static/'
        mode = 'prod'

    request.session['error_code'] = []    
    request.session['csv_mode'] = 'w'

    # identify model and set directory
    request.session['rbafilename'] = json.loads(list(request.POST.items())[0][0])['modelname']
    if socket.gethostname() == 'timputer': request.session['newdir'] = 'simulator/static/python/models/%s' %(request.session['rbafilename'][:-4])
    else: request.session['newdir'] = os.getcwd() + '/rba/static/python/models/%s' %(request.session['rbafilename'][:-4])

    # load model and prepare parameters for change
    wrapper = load_local(request.session['newdir'])
    parameter_values = wrapper.current_parameter_values

    try:
        # get logfile, save it, and create link to download
        try: os.mkdir(pre_path + '%s'%(request.session['rbafilename'][:-4]))
        except: pass
        log_path = pre_path + '%s/changelog_%s.csv' %(request.session['rbafilename'][:-4], request.session['session_id'])
        logfile_content = wrapper.get_change_log()
        if mode == 'dev': 
            logfile_content.to_csv(log_path, index=None, sep=',', mode=request.session['csv_mode'])
            request.session['log_path'] = '../static/results/%s/changelog_%s.csv'%(request.session['rbafilename'][:-4], request.session['session_id'])
        else:
            logfile_content.to_csv(log_path, index=None, sep=',', mode=request.session['csv_mode'])
            request.session['log_path'] = '../static/%s/changelog_%s.csv'%(request.session['rbafilename'][:-4], request.session['session_id'])
        request.session['csv_mode'] = 'w'
    except:
        request.session['error_code'].append('Could not create Logfile for this model.')

    # model download path
    if mode == 'dev': request.session['dl_path'] = '../static/python/models/%s/%s' %(request.session['rbafilename'][:-4], request.session['rbafilename'])
    else: request.session['dl_path'] = '/static/python/models/%s/%s' %(request.session['rbafilename'][:-4], request.session['rbafilename'])

    # parameters
    mpl = []
    for p in parameter_values:
        df = parameter_values[p]
        for i,j in df.iterrows():
            mpl.append([p, j[p]])
            break
    request.session['model_parameters_list'] = mpl

    # medium species
    mss = wrapper.get_changeable_medium_species()
    msl = []
    for i,j in mss.iterrows():
        msl.append([j['Parameter'], j['Type']])
    request.session['model_species_list'] = msl
    
    request.session.modified = True
    return HttpResponse('ok')


def load_local(path):
    '''
    loads wrapper model
    '''
    return rba_websimulator_interface.RBA_websimulator_interface(path)


@csrf_exempt
def simulate(request):
    '''
    Simulate model
    '''
    request.session['error_code'] = []    
    cplex_error = False
    if socket.gethostname() == 'timputer':
        pre_path = 'simulator/static/results/'
        mode = 'dev'
    else:
        pre_path = os.getcwd() + '/rba/static/'
        mode = 'prod'

    try:
        parameters = json.loads(list(request.POST.items())[0][0])['parameters']
        species = json.loads(list(request.POST.items())[0][0])['species']
    except:
        parameters = {}
        species = {}

    wrapper = load_local(request.session['newdir'])
    if request.session['first_sim'] == 'Nope':
        try:
            if mode == 'dev':
                success = wrapper.replay_from_logfile(file_path = 'simulator/static/results/%s/changelog_%s.csv'%(request.session['rbafilename'][:-4], request.session['session_id']))
            elif mode == 'prod':
                success = wrapper.replay_from_logfile(file_path = 'rba/static/%s/changelog_%s.csv'%(request.session['rbafilename'][:-4], request.session['session_id']))
            if success != 'ok': request.session['error_code'].append('Repeated simulation failed due to internal error.')                
        except:
            request.session['error_code'].append('Could not correctly replay in %s.'%os.getcwd())
    else: request.session['first_sim'] = 'Nope'


    if parameters == {} and species == {}:
        try: wrapper.set_default_parameters()
        except: request.session['error_code'].append('The default parameters could not be set. Is the model valid?')
    else:
        for s in species:
            wrapper.set_medium_component(s, new_value=float(species[s]))
        #valid = ['k_UpperGlycolysis_efficiency_foreward', 'k_CitratCycle_efficiency_foreward', 'k_AA_synthesis_efficiency_foreward', 'Membrane_density', 'Cytoplasm_density',
        #'nonenzymatic_proteins_Cytoplasm', 'k_fermentation_ATP_generation_efficiency_foreward', 'k_respiration_efficiency_foreward', 'k_fermentation_NAD_regeneration_efficiency_foreward', 'ribosome_capacity', 'maintenance_atp']
        for p in parameters:
            wrapper.set_parameter_multiplier(p, new_value=float(parameters[p]))
            wrapper.determine_parameter_values(model_parameter=p, x_min=0, x_max=1, intervals=100)
            

    try:
        wrapper.rba_session.findMaxGrowthRate()
        wrapper.rba_session.recordResults('M_Glucose')
        wrapper.rba_session.writeResults(session_name='Test')
    except: request.session['error_code'].append('Growth rate of model could not be determined.')
    
    try:
        # get logfile, save it, and create link to download
        try: os.mkdir(pre_path + '%s'%(request.session['rbafilename'][:-4]))
        except: pass
        log_path = pre_path + '%s/changelog_%s.csv' %(request.session['rbafilename'][:-4], request.session['session_id'])
        logfile_content = wrapper.get_change_log()
        if mode == 'dev': 
            logfile_content.to_csv(log_path, index=None, sep=',', mode=request.session['csv_mode'])
            request.session['log_path'] = '../static/results/%s/changelog_%s.csv'%(request.session['rbafilename'][:-4], request.session['session_id'])
        else:
            logfile_content.to_csv(log_path, index=None, sep=',', mode=request.session['csv_mode'])
            request.session['log_path'] = '../static/%s/changelog_%s.csv'%(request.session['rbafilename'][:-4], request.session['session_id'])
        request.session['csv_mode'] = 'w'
    except:
        request.session['error_code'].append('Could not create Logfile for this model.')

    try:
        # create CSV files, save em, and create links to download
        try: os.mkdir(pre_path + '%s'%(request.session['rbafilename'][:-4]))
        except: pass
        wrapper.rba_session.SimulationData.exportCSV()
        csv_files = wrapper.rba_session.SimulationData.getCSVFiles()
        for cf_key in csv_files:
            csv_file = csv_files[cf_key]
            csv_path = pre_path + '%s/%s_%s' %(request.session['rbafilename'][:-4], request.session['session_id'], cf_key)
            current_paths = request.session['csv_paths']
            if mode == 'dev': current_paths.append('../static/results/%s/%s_%s' %(request.session['rbafilename'][:-4], request.session['session_id'], cf_key))
            else: current_paths.append('../static/%s/%s_%s' %(request.session['rbafilename'][:-4], request.session['session_id'], cf_key))
            f = open(csv_path, 'w+')
            f.write(csv_file)
            f.close()
            request.session['csv_paths'] = current_paths
    except: request.session['error_code'].append('Could not create CSV output.')

    try:
        # create Escher map file, save it, and create link to download
        try: os.mkdir(pre_path + '%s'%(request.session['rbafilename'][:-4]))
        except: pass
        emap_path = pre_path + '%s/eschermap_%s.json' %(request.session['rbafilename'][:-4],request.session['session_id'])
        wrapper.rba_session.SimulationData.exportEscherMap()
        emap_content = wrapper.rba_session.SimulationData.getEscherMap()   
        f = open(emap_path, 'w+')
        f.write(emap_content)
        f.close()
        if mode == 'dev': request.session['emap_path'] = '../static/results/%s/eschermap_%s.json'%(request.session['rbafilename'][:-4], request.session['session_id'])
        else: request.session['emap_path'] = '../static/%s/eschermap_%s.json'%(request.session['rbafilename'][:-4], request.session['session_id'])
    except:
        request.session['error_code'].append('Could not create Escher Map for this model.')

    try:
        # create Proteomap file, save it, and create link to download
        try: os.mkdir(pre_path + '%s'%(request.session['rbafilename'][:-4]))
        except: pass
        proteomap_path = pre_path + '%s/proteomap_%s.tsv' %(request.session['rbafilename'][:-4], request.session['session_id'])
        wrapper.rba_session.SimulationData.exportProteoMap()
        proteomap_content = wrapper.rba_session.SimulationData.getProteoMap()
        f = open(proteomap_path, 'w+')
        f.write(proteomap_content)
        f.close()
        if mode == 'dev': request.session['proteomap_path'] = '../static/results/%s/proteomap_%s.tsv'%(request.session['rbafilename'][:-4], request.session['session_id'])
        else: request.session['proteomap_path'] = '../static/%s/proteomap_%s.tsv'%(request.session['rbafilename'][:-4], request.session['session_id'])
    except:
        request.session['error_code'].append('Could not create Proteomap for this model.')

    try:
        # create SBtab Document, save it, and create link to download
        sbtab_path = pre_path + '%s/sbtab_%s.tsv' %(request.session['rbafilename'][:-4], request.session['session_id'])
        wrapper.rba_session.SimulationData.exportSBtab(filename='Sbtab_Results_Glucose_Screen', rba=True)
        sbtab_content = wrapper.rba_session.SimulationData.getSBtabDoc()
        f = open(sbtab_path, 'w+')
        f.write(sbtab_content.to_str())
        f.close()
        if mode == 'dev': request.session['sbtab_path'] = '../static/results/%s/sbtab_%s.tsv'%(request.session['rbafilename'][:-4], request.session['session_id'])
        else: request.session['sbtab_path'] = '../static/%s/sbtab_%s.tsv'%(request.session['rbafilename'][:-4], request.session['session_id'])
    except:
        request.session['error_code'].append('Could not create SBtab for this model.')

    if request.session['first_sim'] == 'Yes':
        request.session['first_sim'] = 'Nope'

    request.session['plot_path'] = False
    request.session['status'] = True
    request.session.modified = True

    return HttpResponse('ok')


@csrf_exempt
def undolast(request):
    '''
    Undo the last simulation step
    '''
    request.session['error_code'] = []
    if socket.gethostname() == 'timputer':
        pre_path = 'simulator/static/results/'
        mode = 'dev'
    else:
        pre_path = os.getcwd() + '/rba/static/'
        mode = 'prod'

    try:    
        wrapper = load_local(request.session['newdir'])
        try:
            if mode == 'dev':
                success = wrapper.replay_from_logfile(file_path = 'simulator/static/results/%s/changelog_%s.csv'%(request.session['rbafilename'][:-4], request.session['session_id']))
            elif mode == 'prod':
                success = wrapper.replay_from_logfile(file_path = 'rba/static/%s/changelog_%s.csv'%(request.session['rbafilename'][:-4], request.session['session_id']))
            if not success: request.session['error_code'].append('Repeated simulation failed due to internal error.')                
        except:
            request.session['error_code'].append('Could not correctly replay in %s.'%os.getcwd())
    except:
        request.session['error_code'].append('Could not regenerate model wrapper for plotting.')

    try:
        wrapper.undo_last_change()
    except: request.session['error_code'].append('Could not undo last change')
    
    try:
        # get logfile, save it, and create link to download
        log_path = pre_path + '%s/changelog_%s.csv' %(request.session['rbafilename'][:-4], request.session['session_id'])
        logfile_content = wrapper.get_change_log()
        if mode == 'dev': 
            logfile_content.to_csv(log_path, index=None, sep=',', mode=request.session['csv_mode'])
            request.session['log_path'] = '../static/results/%s/changelog_%s.csv'%(request.session['rbafilename'][:-4], request.session['session_id'])
        else:
            logfile_content.to_csv(log_path, index=None, sep=',', mode=request.session['csv_mode'])
            request.session['log_path'] = '../static/%s/changelog_%s.csv'%(request.session['rbafilename'][:-4], request.session['session_id'])
    except:
        request.session['error_code'].append('Could not create Logfile for this model.')
   
    request.session['plot_path'] = False
    request.session['status'] = True
    request.session.modified = True

    return HttpResponse('ok')


@csrf_exempt
def plot(request):
    '''
    Plot a parameter
    '''
    request.session['error_code'] = []    
    if socket.gethostname() == 'timputer':
        pre_path = 'simulator/static/results/'
        mode = 'dev'
    else:
        pre_path = os.getcwd() + '/rba/static/'
        mode = 'prod'

    try: parameter = json.loads(list(request.POST.items())[0][0])['plot_parameter']
    except:
        request.session['error_code'].append('Parameter was not submitted succesfully.')
        parameter = None

    try:
        wrapper = load_local(request.session['newdir'])
        # Simulation does not seem to be necessary
        '''try:
            if mode == 'dev':
                success = wrapper.replay_from_logfile(file_path = 'simulator/static/results/%s/changelog_%s.csv'%(request.session['rbafilename'][:-4],request.session['session_id']))
            elif mode == 'prod':
                success = wrapper.replay_from_logfile(file_path = 'rba/static/%s/changelog_%s.csv'%(request.session['rbafilename'][:-4],request.session['session_id']))
            if not success: request.session['error_code'].append('Repeated simulation failed due to internal error.')                
        except:
            request.session['error_code'].append('Could not correctly replay in %s.'%os.getcwd())'''
    except:
        request.session['error_code'].append('Could not regenerate model wrapper for plotting.')

    try:
        try: os.mkdir(pre_path + '%s'%(request.session['rbafilename'][:-4]))
        except: pass
        df = wrapper.get_plot_values(model_parameter=parameter)
        fig,ax = plt.subplots()

        plt.plot(df[list(df.columns)[0]], df['Original values'])
        plt.plot(df[list(df.columns)[0]], df['Current values'])
        plt.legend(['Original', 'Current'])
        plt.title(parameter)
        plt.xlabel(list(df.columns)[0])
        plt.ylabel(parameter)
        if mode == 'dev':
            plot_path = 'simulator/static/results/%s'%request.session['rbafilename'][:-4]
        else:
            plot_path = os.getcwd() + '/rba/static/results/%s'%request.session['rbafilename'][:-4]
        
        plt.savefig(plot_path + '/plot_%s.png'%request.session['session_id'])
        if mode == 'dev':
            request.session['plot_path'] = '../static/results/%s/plot_%s.png'%(request.session['rbafilename'][:-4], request.session['session_id'])
        else:
            request.session['plot_path'] = '../static/results/%s/plot_%s.png'%(request.session['rbafilename'][:-4], request.session['session_id'])
        #plt.show()
        #plt.close()
    except:
        request.session['error_code'].append('Parameter could not be plotted.')

    request.session.modified = True

    return HttpResponse('ok')


@csrf_exempt
def prepare_log(request):
    '''
    Plot a parameter
    '''
    return 'ok'