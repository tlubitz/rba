from django.urls import path
from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('clearsession/', views.clearsession, name='clearsession'),
    path('simulate/', views.simulate, name='simulate'),
    path('loadmodel/', views.loadmodel, name='loadmodel'),
    path('undolast/', views.undolast, name='undolast'),
    path('plot/', views.plot, name='plot'),
]