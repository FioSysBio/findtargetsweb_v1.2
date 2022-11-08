import os

from django.core.cache import cache
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import render
from libsbml import SBMLReader

from front.forms import SBMLFileForm, Passo1Form
from front.pipelineFindTargets import MyThread, FindTargets


# Create your views here.
def index(request):
    if request.method == 'POST':
        form = SBMLFileForm(request.POST, request.FILES)

        organism_in = form.data["organism"]
        email_in = form.data["email"]
        name_in = form.data["name"]
        sbmlfile_in = request.FILES["file"]

        request.session['organism'] = organism_in
        request.session['email'] = email_in
        request.session['name'] = name_in
        # request.session['file'] = sbmlfile_in

        dict_return_full = {}
        dict_return_form = {
            'name': name_in,
            'email': email_in,
            'organism': organism_in,
            'file': sbmlfile_in,
        }

        if form.is_valid():

            # reader = SBMLReader()
            # document = reader.readSBMLFromString(sbmlfile_in.read().decode('utf-8'))

            #/home/thiago/projetos/fiocruz/findtargetsweb/front/sbmls
            # dir_sbmls = os.path.dirname(os.path.realpath(__file__)) + '/sbmls/'
            # fileSBML = dir_sbmls + document.model.name + '.xml'
            # if not os.path.exists(fileSBML):
            #     with open(fileSBML, 'wb+') as destination:
            #         for chuck in sbmlfile_in.chunks():
            #             destination.write(chuck)

            model = FindTargets().readModel(sbmlfile_in.read().decode('utf-8'))
            # model.solver = "glpk"
            dict_return_validate = FindTargets().validateModel(model)

            cache.set('model', model)

            if dict_return_validate['ehParaFazer']:
                request.session['ehParaFazer'] = True
            else:
                request.session['ehParaFazer'] = False

            dict_return_full.update(dict_return_form)
            dict_return_full.update(dict_return_validate)
            del dict_return_full['file']
            cache.set('dict_return_full', dict_return_full)

            return HttpResponseRedirect('/passo1/')
            # return HttpResponseRedirect('FindTargetsWEB/passo1/')  # servidor
    else:
        request.session.flush()
        cache.clear()
        # cache.close()
        form = SBMLFileForm()
    return render(request, 'index.html', {'form': form})
    # return render(request, 'index.html');


def passo1(request):
    if request.method == 'POST':  # Aqui está clicando no submit
        form = Passo1Form(request.POST, request.FILES)
        organism_in2 = request.session["organism"]
        email_in2 = request.session["email"]
        name_in2 = request.session["name"]
        method = form.data['method']
        # sbmlfile_in2 = request.session['file']
        modelPOG = cache.get('model')
        # modelPOG = FindTargets().readModel(sbmlfile_in2)

        dict_return_full = {}
        dict_return_form = {
            'name': name_in2,
            'email': email_in2,
            'organism': organism_in2,
            'method': method
        }

        if form.is_valid():

            if request.session['ehParaFazer']:
                t = MyThread(name_in2, organism_in2, email_in2, modelPOG, method)
                t.start()

            dict_return_full.update(dict_return_form)
            return HttpResponseRedirect('/passo2/')
            # return HttpResponseRedirect('FindTargetsWEB/passo2/')  # servidor

    else:  # Aqui está entrando na pagina passo1.html
        dict_return_full = cache.get('dict_return_full')
        organism_in2 = request.session["organism"]
        email_in2 = request.session["email"]
        name_in2 = request.session["name"]
        valInitialSolutionFmt = dict_return_full['valInitialSolutionFmt']
        message = dict_return_full['message']

        form = Passo1Form()

        dict_form = {
            'name': name_in2,
            'email': email_in2,
            'organism': organism_in2,
            'valInitialSolutionFmt': valInitialSolutionFmt,
            'message': message,
            'form': form
        }

    return render(request, 'passo1.html', dict_form)


def passo2(request):
    request.session.flush()
    cache.clear()
    cache.close()
    return render(request, 'passo2.html')


def download(request):
    dir_path = os.path.dirname(os.path.realpath(__file__))  # identifica o local real onde esse arquivo esta
    response = HttpResponse(open( "static/User_Manual_FindTargetsWeb_v1.1.pdf", 'rb').read())
    response['Content-Type'] = 'application/pdf'
    response['Content-Disposition'] = 'attachment; filename=User_Manual_FindTargetsWeb_v1.1.pdf'
    return response