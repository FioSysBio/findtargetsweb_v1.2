from django.urls import path
from . import views

# app_name = "front"

urlpatterns = [
    path('findtargetsweb/', views.index, name='index'),
    path('findtargetsweb/passo1/', views.passo1, name='passo1'),
    path('findtargetsweb/passo2/', views.passo2, name='passo2'),
    path('findtargetsweb/download/', views.download, name='download')
]