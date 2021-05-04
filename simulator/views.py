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
import cplex
import random
import matplotlib
#matplotlib.use('TkAgg')
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .static.python.rbatools import rba_wrapper
from .static.python.rbatools import rba_websimulator_interface

global wrapper

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
                'model_parameters_list',
                'model_species_list',
                'plot_path',
                'wrapper']:
        if not request.session.get(var, None):
            request.session[var] = False

    for var in ['error_code', 'csv_paths']:
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
    global wrapper
    while request.session['wrapper'] == True:
        try:
            wrapper = False
            request.session['wrapper'] = False
        except: pass

    # delete current project directory (only if it was an uploaded model)
    if not request.session['newdir'].startswith('simulator/static/python/models/') and not request.session['newdir'].startswith('/home/TimoSan/'):
        try: shutil.rmtree(request.session['newdir'])
        except: print('Cannot delete %s' %request.session['newdir'])

    # delete session variables
    keys = ['rbafilezip', 'rbafilename', 'newdir', 'error_code', 'emap_path', 'proteomap_path', 'csv_paths', 'sbtab_path', 'dl_path', 'plot_path', 'status', 'errors', 'model_parameters_list', 'model_species_list']
    for key in keys:
        try: del request.session[key]
        except: print('Cannot delete %s' %(key))

    request.session.modified = True

    return HttpResponse('ok')


@csrf_exempt
def loadmodel(request):
    '''
    load an existing model from server
    '''
    request.session['csv_mode'] = 'w'
    # identify model and set directory
    request.session['rbafilename'] = json.loads(list(request.POST.items())[0][0])['modelname']
    if socket.gethostname() == 'timputer': request.session['newdir'] = 'simulator/static/python/models/%s' %(request.session['rbafilename'][:-4])
    else: request.session['newdir'] = '/home/TimoSan/rba/static/python/models/%s' %(request.session['rbafilename'][:-4])

    # load model and prepare parameters for change xxx
    global wrapper
    wrapper = rba_websimulator_interface.RBA_websimulator_interface(request.session['newdir'])
    request.session['wrapper'] = True
    parameter_values = wrapper.current_parameter_values

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


@csrf_exempt
def simulate(request):
    '''
    Simulate model
    '''
    global wrapper
    cplex_error = False
    if socket.gethostname() == 'timputer':
        pre_path = 'simulator/static/results/'
        mode = 'dev'
    else:
        pre_path = '/home/TimoSan/rba/static/'
        mode = 'prod'

    try:
        parameters = json.loads(list(request.POST.items())[0][0])['parameters']
        species = json.loads(list(request.POST.items())[0][0])['species']
    except:
        parameters = {}
        species = {}

    if parameters == {} and species == {}:
        while request.session['wrapper'] == True:
            try: wrapper.set_default_parameters()
            except: pass #request.session['error_code'].append('The default parameters could not be set. Is the model valid?')
    else:
        for s in species:
            wrapper.set_medium_component(s, new_value=float(species[s]))
        for p in parameters:
            wrapper.set_parameter_multiplier(p, new_value=float(parameters[p]))
            wrapper.determine_parameter_values(model_parameter=p, x_min=0, x_max=1, intervals=100)

    wrapper.rba_session.findMaxGrowthRate()
    try: wrapper.rba_session.recordResults('Glucose')
    except: pass
    wrapper.rba_session.writeResults(session_name='Test')

    try:
        # get logfile, save it, and create link to download
        try: os.mkdir(pre_path + '%s'%(request.session['rbafilename'][:-4]))
        except: pass
        log_path = pre_path + '%s/changelog.csv' %request.session['rbafilename'][:-4]
        logfile_content = wrapper.get_change_log()
        if mode == 'dev': 
            logfile_content.to_csv(log_path, header=None, index=None, sep=',', mode=request.session['csv_mode'])
            request.session['log_path'] = '../static/results/%s/changelog.csv'%request.session['rbafilename'][:-4]
        else:
            logfile_content.to_csv(log_path, header=None, index=None, sep=',', mode=request.session['csv_mode'])
            request.session['log_path'] = '../static/%s/changelog.csv'%request.session['rbafilename'][:-4]
        request.session['csv_mode'] = 'w'
    except:
        request.session['error_code'].append('Could not create Logfile for this model.')
   
    if mode == 'dev': request.session['dl_path'] = '../static/python/models/%s/%s' %(request.session['rbafilename'][:-4], request.session['rbafilename'])
    else: request.session['dl_path'] = '/static/python/models/%s/%s' %(request.session['rbafilename'][:-4], request.session['rbafilename'])

    request.session['plot_path'] = False
    request.session['status'] = True
    request.session.modified = True

    return HttpResponse('ok')


@csrf_exempt
def undolast(request):
    '''
    Undo the last simulation step
    '''
    if socket.gethostname() == 'timputer':
        pre_path = 'simulator/static/results/'
        mode = 'dev'
    else:
        pre_path = '/home/TimoSan/rba/static/'
        mode = 'prod'

    global wrapper
    try:
        wrapper.undo_last_change()
    except: request.session['error_code'].append('Could not undo last change')
    
    try:
        # get logfile, save it, and create link to download
        try: os.mkdir(pre_path + '%s'%(request.session['rbafilename'][:-4]))
        except: pass
        log_path = pre_path + '%s/changelog.csv' %request.session['rbafilename'][:-4]
        logfile_content = wrapper.get_change_log()
        if mode == 'dev': 
            logfile_content.to_csv(log_path, header=None, index=None, sep=',', mode=request.session['csv_mode'])
            request.session['log_path'] = '../static/results/%s/changelog.csv'%request.session['rbafilename'][:-4]
        else:
            logfile_content.to_csv(log_path, header=None, index=None, sep=',', mode=request.session['csv_mode'])
            request.session['log_path'] = '../static/%s/changelog.csv'%request.session['rbafilename'][:-4]
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
    if socket.gethostname() == 'timputer':
        pre_path = 'simulator/static/results/'
        mode = 'dev'
    else:
        pre_path = '/home/TimoSan/rba/static/'
        mode = 'prod'

    try: parameter = json.loads(list(request.POST.items())[0][0])['plot_parameter']
    except:
        request.session['error_code'].append('Parameter was not submitted succesfully.')
        parameter = None

    global wrapper
    #try:
    df = wrapper.get_plot_values(model_parameter=parameter)
    #fig = plt.figure(1)
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
        plot_path = '/home/TimoSan/rba/static/results/%s'%request.session['rbafilename'][:-4]
    
    #try: os.mkdir(plot_path)
    #except: print('Could not create new directory for plot results.')
    '''import pickle
    pickle.dump(fig, open(plot_path + '/FigureObject.fig.pickle', 'wb'))
    request.session['plot_path'] = '../static/results/%s/FigureObject.fig.pickle'%request.session['rbafilename'][:-4]'''
    plt.savefig(plot_path + '/plot.png')
    if mode == 'dev':
        request.session['plot_path'] = '../static/results/%s/plot.png'%request.session['rbafilename'][:-4]
    else:
        request.session['plot_path'] = '../static/results/%s/plot.png'%request.session['rbafilename'][:-4]
        #plt.show()
        #plt.close()
    #except:
    #    request.session['error_code'].append('Parameter could not be plotted.')

    request.session.modified = True

    return HttpResponse('ok')
