from django import forms

class UploadFileForm(forms.Form):
    zipfile = forms.FileField()
