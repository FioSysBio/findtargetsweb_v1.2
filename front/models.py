import datetime

from django.db import models

# Create your models here.

class Modelo (models.Model):
    nome = models.CharField(max_length=200)
    nome_arquivo = models.CharField(max_length=200)
    arquivo = models.BinaryField()
    data_criacao = models.DateTimeField(default=datetime.datetime.now(), blank=True)
