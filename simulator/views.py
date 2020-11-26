from django.http import HttpResponseRedirect
from django.shortcuts import render, render_to_response
from simulator.forms import UploadFileForm
from django.conf import settings
from django.http import HttpResponse
from django.views.decorators.http import require_POST
from django.shortcuts import redirect

import os
import zipfile


def index(request):
    '''
    main page handling rba file upload and fine tuning parameters
    '''
    # create essential variables
    error_code = False
    for var in ['rbafilezip', 'rbafilename']:
        if not request.session.get(var, None):
            request.session[var] = False
    

    if request.method == 'POST':
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            if not (request.FILES['zipfile'].name).endswith('.zip'):
                error_code = 'This is not a zip file. Please only use zip files'
            else:
                # first create new directory
                media_directory = settings.MEDIA_ROOT
                current_rba_name = request.FILES['zipfile'].name[:-4]
                new_dir = os.path.join(media_directory, current_rba_name)
                try: os.makedir(new_dir)
                except: print('Could not create new directory for zip file contents.')

                # then unzip the zip file
                try:
                    newzip = zipfile.ZipFile(request.FILES['zipfile'], 'r')
                    newzip.extractall(new_dir)
                    request.session['rbafilename'] = request.FILES['zipfile'].name
                    return HttpResponseRedirect('/simulator')
                except:
                    error_code = 'Could not unzip file. Is it valid?'             

    else:
        form = UploadFileForm()

    request.session.modified = True

    return render(request, 'index.html', {'form': form,
                                          'rbafilename': request.session['rbafilename'],
                                          'error_code': error_code })


def clearsession(request):
    '''
    clears all session variables
    '''
    keys = ['rbafilezip', 'rbafilename']
    for key in keys:
        try: del request.session[key]
        except: print('Cannot delete %s' %(key))

    request.session.modified = True

    return HttpResponse('ok')
