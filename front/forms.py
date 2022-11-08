from django import forms
from django.core.validators import FileExtensionValidator
from front.pipelineFindTargets import FindTargets

class SBMLFileForm(forms.Form):
    name = forms.CharField(max_length=100, widget=forms.TextInput(attrs={'placeholder': 'Your Name', 'tabindex': '1'}))
    email = forms.EmailField(widget=forms.EmailInput(attrs={'placeholder': 'Your E-mail address', 'tabindex': '2'}))
    organism = forms.ChoiceField(choices=FindTargets().list_organism(), widget=forms.Select(attrs={'tabindex': '3'}))
    file = forms.FileField(validators=[FileExtensionValidator(allowed_extensions=['xml', 'sbml'])])
    
class Passo1Form(forms.Form):
    method = forms.ChoiceField(choices=(("1", "FBA+FVA"), ("2", "Only FBA")), widget=forms.Select(attrs={'tabindex': '3'}))

class Passo2Form(forms.Form):
    method = forms.ChoiceField(choices=(("1", "FBA+FVA"), ("2", "Only FBA")), widget=forms.Select(attrs={'tabindex': '3'}))