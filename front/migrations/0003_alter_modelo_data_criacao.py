# Generated by Django 4.1.1 on 2022-11-08 18:17

import datetime
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('front', '0002_alter_modelo_arquivo_alter_modelo_data_criacao'),
    ]

    operations = [
        migrations.AlterField(
            model_name='modelo',
            name='data_criacao',
            field=models.DateTimeField(blank=True, default=datetime.datetime(2022, 11, 8, 15, 17, 27, 382558)),
        ),
    ]
