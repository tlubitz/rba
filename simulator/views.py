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

from .static.python.rbatools import rba_wrapper


def index(request):
    '''
    main page handling rba file upload and fine tuning parameters
    '''
    # create essential variables
    for var in ['rbafilezip', 'rbafilename', 'emap_path', 'proteomap_path', 'sbtab_path', 'status', 'errors']:
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
                                          'sbtab_path': request.session['sbtab_path']})


def clearsession(request):
    '''
    clears all session variables
    '''
    # delete current project directory (only if it was an uploaded model)
    if not request.session['newdir'].startswith('simulator/static/python/models/') and not request.session['newdir'].startswith('/home/TimoSan/'):
        try: shutil.rmtree(request.session['newdir'])
        except: print('Cannot delete %s' %request.session['newdir'])

    # delete session variables
    keys = ['rbafilezip', 'rbafilename', 'newdir', 'error_code', 'emap_path', 'proteomap_path', 'csv_paths', 'sbtab_path', 'status', 'errors']
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
    request.session['rbafilename'] = json.loads(list(request.POST.items())[0][0])['modelname']
    if socket.gethostname() == 'timputer': request.session['newdir'] = 'simulator/static/python/models/%s' %(request.session['rbafilename'][:-4])
    else: request.session['newdir'] = '/home/TimoSan/rba/static/python/models/%s' %(request.session['rbafilename'][:-4])

    return HttpResponse('ok')


def simulate(request):
    '''
    Simulate model
    '''
    if socket.gethostname() == 'timputer':
        pre_path = 'simulator/static/results/'
        mode = 'dev'
    else:
        pre_path = '/home/TimoSan/rba/static/'
        mode = 'prod'

    try:
        wrapper = rba_wrapper.Wrapper()
        wrapper.set_path(request.session['newdir'])
    except: print('Could not create RBA wrapper.\n')

    try: wrapper.create_simulation()
    except: request.session['error_code'].append('The simulation cannot be initialised. Is the model valid?')

    wrapper.set_default_parameters()
    try: wrapper.set_default_parameters()
    except: request.session['error_code'].append('The default parameters could not be set. Is the model valid?')
    
    try: wrapper.write_results()
    except: request.session['error_code'].append('Model could not be simulated.')

    try:
        # create CSV files, save em, and create links to download
        try: os.mkdir(prepath + '%s'%(request.session['rbafilename'][:-4]))
        except: pass
        csv_files = wrapper.get_csvs()
        for cf_key in csv_files:
            csv_file = csv_files[cf_key]
            csv_path = pre_path + '%s/%s' %(request.session['rbafilename'][:-4], cf_key)
            current_paths = request.session['csv_paths']
            if mode == 'dev': current_paths.append('../static/results/%s/%s' %(request.session['rbafilename'][:-4], cf_key))
            else: current_paths.append('../static/%s/%s' %(request.session['rbafilename'][:-4], cf_key))
            f = open(csv_path, 'w+')
            f.write(csv_file)
            f.close()
            request.session['csv_paths'] = current_paths
    except: request.session['error_code'].append('Could not create CSV output.')
   
    try:
        # create Escher map file, save it, and create link to download
        try: os.mkdir(pre_path + '%s'%(request.session['rbafilename'][:-4]))
        except: pass
        emap_path = pre_path + '%s/eschermap.json' %request.session['rbafilename'][:-4]
        emap_content = wrapper.get_eschermap()
        f = open(emap_path, 'w+')
        f.write(emap_content)
        f.close()
        if mode == 'dev': request.session['emap_path'] = '../static/results/%s/eschermap.json'%request.session['rbafilename'][:-4]
        else: request.session['emap_path'] = '../static/%s/eschermap.json'%request.session['rbafilename'][:-4]        
    except:
        request.session['error_code'].append('Could not create Escher Map for this model.')
   
    try:
        # create Proteomap file, save it, and create link to download
        try: os.mkdir(pre_path + '%s'%(request.session['rbafilename'][:-4]))
        except: pass
        proteomap_path = pre_path + '%s/proteomap.tsv' %request.session['rbafilename'][:-4]
        proteomap_content = wrapper.get_proteomap()
        f = open(proteomap_path, 'w+')
        f.write(proteomap_content)
        f.close()
        if mode == 'dev': request.session['proteomap_path'] = '../static/results/%s/proteomap.tsv'%request.session['rbafilename'][:-4]
        else: request.session['proteomap_path'] = '../static/%s/proteomap.tsv'%request.session['rbafilename'][:-4]
    except:
        request.session['error_code'].append('Could not create Proteomap for this model.')
   
    try:
        # create SBtab Document, save it, and create link to download
        try: os.mkdir(pre_path + '%s'%(request.session['rbafilename'][:-4]))
        except: pass
        sbtab_path = pre_path + '%s/sbtab.tsv' %request.session['rbafilename'][:-4]
        sbtab_content = wrapper.get_sbtab()
        f = open(sbtab_path, 'w+')
        f.write(sbtab_content.to_str())
        f.close()
        if mode == 'dev': request.session['sbtab_path'] = '../static/results/%s/sbtab.tsv'%request.session['rbafilename'][:-4]
        else: request.session['sbtab_path'] = '../static/%s/sbtab.tsv'%request.session['rbafilename'][:-4]
    except:
        request.session['error_code'].append('Could not create SBtab for this model.')
    
    request.session['status'] = True
    request.session.modified = True

    return HttpResponse('ok')
