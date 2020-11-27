from django.http import HttpResponseRedirect
from django.shortcuts import render, render_to_response
from simulator.forms import UploadFileForm
from django.conf import settings
from django.http import HttpResponse
from django.views.decorators.http import require_POST
from django.shortcuts import redirect
from .static.python import rba_wrapper

import os
import shutil
import zipfile


def index(request):
    '''
    main page handling rba file upload and fine tuning parameters
    '''
    # create essential variables
    for var in ['rbafilezip', 'rbafilename', 'emap_path']:
        if not request.session.get(var, None):
            request.session[var] = False
    if not request.session.get('error_code', None):
        request.session['error_code'] = ''

    if request.method == 'POST':
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            if not (request.FILES['zipfile'].name).endswith('.zip'):
                request.session['error_code'] += 'This is not a zip file. Please only use zip files.\n'
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
                except: request.session['error_code'] += 'Could not unzip file. Is it valid?\n'
                
                return HttpResponseRedirect('/simulator')

    else:
        form = UploadFileForm()

    request.session.modified = True

    return render(request, 'index.html', {'form': form,
                                          'rbafilename': request.session['rbafilename'],
                                          'error_code': request.session['error_code'],
                                          'emap_path': request.session['emap_path']})


def clearsession(request):
    '''
    clears all session variables
    '''
    # delete current project directory
    try: shutil.rmtree(request.session['newdir'])
    except: print('Cannot delete %s' %request.session['newdir'])

    # delete session variables
    keys = ['rbafilezip', 'rbafilename', 'newdir', 'error_code', 'emap_path']
    for key in keys:
        try: del request.session[key]
        except: print('Cannot delete %s' %(key))

    request.session.modified = True

    return HttpResponse('ok')


def simulate(request):
    '''
    Simulate model
    '''
    try:
        wrapper = rba_wrapper.Wrapper()
        wrapper.set_path(request.session['newdir'])
    except: print('Could not create RBA wrapper.\n')

    try: wrapper.create_simulation()
    except: request.session['error_code'] += 'The simulation cannot be initialised. Is the model valid?\n'

    try: wrapper.set_default_parameters()
    except: request.session['error_code'] += 'The default parameters could not be set. Is the model valid?\n'
    
    try: wrapper.write_results()
    except: request.session['error_code'] += 'Model could not be simulated.\n'

    try: csv_file = wrapper.get_csv()
    except: request.session['error_code'] += 'Could not create CSV output.\n'
   
    try:
        # create Escher map file, save it, and create link to download
        # print(settings.DEV) #A RE THESE PATHS ALSO OK IN PRODUCTION?
        os.mkdir('simulator/static/results/%s'%(request.session['rbafilename'][:-4]))
        emap_path = 'simulator/static/results/%s/eschermap.json' %request.session['rbafilename'][:-4]
        emap_content = wrapper.get_eschermap()
        f = open(emap_path, 'w+')
        f.write(emap_content)
        f.close()
        request.session['emap_path'] = '../static/results/%s/eschermap.json'%request.session['rbafilename'][:-4]
    except:
        request.session['error_code'] += 'Could not create Escher Map for this model.\n'

    request.session.modified = True

    return HttpResponse('ok')