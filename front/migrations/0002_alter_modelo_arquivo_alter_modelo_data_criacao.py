# Generated by Django 4.1.1 on 2022-11-10 18:55

import datetime
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('front', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='modelo',
            name='arquivo',
            field=models.BinaryField(),
        ),
        migrations.AlterField(
            model_name='modelo',
            name='data_criacao',
            field=models.DateTimeField(blank=True, default=datetime.datetime(2022, 11, 10, 15, 55, 16, 549957)),
        ),
    ]
