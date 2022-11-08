from django.urls import path
from . import views

# app_name = "front"

urlpatterns = [
    path('', views.index, name='index'),
    path('passo1/', views.passo1, name='passo1'),
    path('passo2/', views.passo2, name='passo2'),
    path('download/', views.download, name='download')
]