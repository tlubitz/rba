from django.contrib import admin
from django.urls import path
from django.views.generic import RedirectView
from django.conf import settings
from django.conf.urls.static import static
from django.urls import include

urlpatterns = [
    path('admin/', admin.site.urls),
    path('simulator/', include('simulator.urls')),
    path('', RedirectView.as_view(url='simulator/')),
] + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)