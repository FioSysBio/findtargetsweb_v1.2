# coding=ISO-8859-1

from __future__ import print_function

import traceback
from datetime import datetime
from bioservices.kegg import KEGG
from bioservices.uniprot import UniProt
from bioservices.ncbiblast import NCBIblast
from bs4 import BeautifulSoup
from zipfile import ZipFile, ZIP_DEFLATED
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders

import xlwt
import requests
import os
import re
import sys
import cobra.manipulation
import pandas as pd
import time
import threading
import smtplib
from requests import Response

from front.models import Modelo


# CLASSE CRIADA COM O OBJETIVO DE MANTER A CHAMADA ASSINCRONA DO METODO DA MODELAGEM
class MyThread(threading.Thread):
    
    name = None
    organism = None
    email = None
    model = None
    method = None
    
    def __init__(self, name, organism, email, modelParam, method):
        threading.Thread.__init__(self)
        self.model = modelParam
        self.organism = organism
        self.name = name
        self.email = email
        self.method = method
        
    def run(self):
        print("A THREAD COMECOU!!!!!!!!")
        # FindTargets.sendMailTest(self);
        FindTargets().mainMethod(self.model, self.organism, self.name, self.email, self.method)

class FindTargets:

    # ADD ALL ORGANISMS MAPPED FROM KEGG TO ATTRIBUTE IN COMBOBOX FORM
    # USED TO FULL OUR COMBOBOX FROM INDEX APP
    def list_organism(self):
        modelos = Modelo.objects.all()
        dir_path = os.path.dirname(os.path.realpath(__file__)) # identifica o local real onde esse arquivo esta
        with open(dir_path+"/organismos.txt", "r") as infile:
            data = infile.read()
        my_list = data.splitlines()
        my_list_aux = []
        my_list_aux.append(("", " --- SELECT --- "))

        for itens in my_list:
            itens_splt = itens.split(",")
            tuple_add = (itens_splt[0], itens_splt[1])
            my_list_aux.append(tuple_add)

        if modelos:
            modelo = modelos.first()
            my_list_aux.append((modelo.nome_arquivo, modelo.nome))

        return tuple(my_list_aux)

    # METODO PARA MONTAR O OBJETO COBRA PARA TRAFEGO NA APLICACAO!
    def readModel(self, sbmlfile):
        # try:
        model = cobra.io.read_sbml_model(sbmlfile)
        # model = cobra.io.read_sbml_model('/home/thiago/projetos/fiocruz/modelos/e_coli_core.xml')
        return model
        # except:
        #     print('Arquivo SBML inválido')
        # finally:
        #     print('Processo finalizado')
        # return cobra.io.sbml3.read_sbml_model(sbmlfile)

    # MODEL VALIDATION METHOD TO CONTINUE OUR EXECUTION
    def validateModel(self, model):
        return_dict_val_model = {}

        numCasasDec = 6

        initialSolution = model.optimize()
        valInitialSolutionFmt = round(float(initialSolution.objective_value), numCasasDec)

        return_dict_val_model['valInitialSolutionFmt'] = valInitialSolutionFmt

        if initialSolution.status != 'optimal':
            return_dict_val_model['message'] = "ERROR! SELECTED MODEL HAS AN ERROR. PLEASE VERIFY AND TRY AGAIN!"
            return_dict_val_model['ehParaFazer'] = False
        else:
            if initialSolution.objective_value == 0.0 or initialSolution.objective_value == -0.0:
                return_dict_val_model['message'] = "SELECTED MODEL DOES NOT GENERATE BIOMASS. PLEASE VERIFY AND TRY AGAIN!"
                return_dict_val_model['ehParaFazer'] = False
            else:
                return_dict_val_model['message'] = "SELECTED MODEL GENERATES BIOMASS. PRESS 'SUBMIT' TO CONTINUE!"
                return_dict_val_model['ehParaFazer'] = True

        return return_dict_val_model

    #
    # MAIN METHOD THAT REALIZE ANALYSIS OF SELECTED NETWORK SBML FILE
    #
    def mainMethod(self, model, organismParam, name, email, method):

        numCasasDec = 6
        TIMEOUT_SECONDS = 500000
        model = model

        # BEFORE, WE NEED TO CREATE FOLDERS TO INCLUDE OUR FILES AND AFTER ZIP THEM
        dataHoraAtual = datetime.now()
        dataHoraAtualFmt = dataHoraAtual.strftime('%Y-%m-%d_%H.%M.%S.%f')
        dir_path = os.path.dirname(os.path.realpath(__file__))
        directory = dir_path+"/results/"+dataHoraAtualFmt
        static_dir = dir_path.replace('/front', '/static')

        # TESTS LINES!
        #directory = dir_path+"/results/"+"testesMeriguetiCCBH"
        #directory = dir_path+"/results/"+"testesMeriguetiPAO1"
        #directory = dir_path+"/results/"+"testesMeriguetiPAO1_2017"
        
        dir_blasts = directory+"/blasts"

        if not os.path.exists(directory):
            os.makedirs(directory)

        if not os.path.exists(dir_blasts):
            os.makedirs(dir_blasts)
        
        start = time.time()
        fileLog = open(directory+"//"+"LOG_EXECUCAO_"+str(dataHoraAtualFmt)+".txt", "w")
        fileLog.write('*** INICIO DA EXECUCAO *** \n')
        fileLog.write('Modelo/Organismo/Nome/Email/Metodo = {0}, {1}, {2}, {3}, {4}'.format(model, organismParam, name, email, method))
        fileLog.write('\n\n\n')
        try:
            ####################################################
            # ITEM 1 - ANALYSIS OF P. AERUGINOSA MODEL SELECTED
            ####################################################
            fileLog.write('INICIO DO ITEM 1\n')
            #self.reportModel(model, directory+"//model_data_"+str(model)+"_.xls")
            fileLog.write('GRAVOU ARQUIVO COM OS DADOS DE GENE/REACAO DA REDE\n')

            initialSolution = model.optimize()
            fileLog.write('GEROU FBA \n')
            fileLog.write('FBA = {0}\n'.format(initialSolution))
            valInitialSolutionFmt = round(float(initialSolution.objective_value), numCasasDec)
            
            ################################################
            # ITEM 2 - PRIORITIZATION OF EXISTING REACTIONS
            ################################################
            fileLog.write('INICIO DO ITEM 2\n')
            fva_result = cobra.flux_analysis.flux_variability_analysis(model, model.reactions[:len(model.reactions)])
            pd.DataFrame.from_dict(fva_result).round(numCasasDec).to_csv(directory+"/01-fva.csv")
            fileLog.write('GEROU FVA \n')
            listReactPerturb = []
            reacoesFVADifZero = open(directory+"//"+"02-reacoesFVADifZero.txt", "w")

            df_fva_result = pd.DataFrame(fva_result)
            df_fva_result = df_fva_result.sort_index()
            itemsFVAResult = df_fva_result.iterrows()
            itemsFBAResult = sorted(initialSolution.fluxes.items())
            # itemsFBAResult = sorted(initialSolution.fluxes.iteritems())
            fileLog.write('INICIO DA VERIFICACAO ENTRE FBA/FVA \n')
            
            for (keyFVA, valueFVA), (keyFBA, valueFBA) in list(zip(itemsFVAResult, itemsFBAResult)):
                minFVA = round(float(valueFVA["minimum"]), numCasasDec)
                maxFVA = round(float(valueFVA["maximum"]), numCasasDec)
                valueFBA = round(float(valueFBA), numCasasDec)

                # Filtrando reacoes que estao ativas no momento
                if maxFVA != 0:
                    # method --> ("1", "FBA+FVA"), ("2", "Only FBA")
                    if method == "1":
                        diferencaReacao = maxFVA - minFVA
                        if diferencaReacao == 0:
                            reacoesFVADifZero.write(str(keyFBA) + '\n')
                            listReactPerturb.append(keyFBA)
                    else:
                        reacoesFVADifZero.write(str(keyFBA) + '\n')
                        listReactPerturb.append(keyFBA)
                    
                    #diferencaReacao = maxFVA - minFVA
                    #reacoesFVADifZero.write("Resultado da diferenca entre max e min = " + str(diferencaReacao))
                    #reacoesFVADifZero.write('\n\n')

                if keyFBA != keyFVA:
                    raise Exception("Por favor, verifique, chave FBA/FVA nao batem na comparacao")

                else:
                    if valueFBA < minFVA:
                        #print(keyFBA)
                        #print(valueFBA)
                        #print(minFVA)
                        #print(maxFVA)
                        raise Exception("Por favor, verifique o modelo. Valor do FBA dessa reacao ser menor que o minimo descrito no FVA")

                    if valueFBA > maxFVA:
                        #print(keyFBA)
                        #print(valueFBA)
                        #print(minFVA)
                        #print(maxFVA)
                        raise Exception("Por favor, verifique o modelo. Valor do FBA dessa reacao ser maior que o maximo descrito no FVA")

            reacoesFVADifZero.close()
            fileLog.write("Lista de reacoes a analisar = " + str(len(listReactPerturb)) + '\n')
            fileLog.write("GENERATED FILE 01-fva.csv\n")
            fileLog.write("GENERATED FILE %s\n" % reacoesFVADifZero.name)
                    
            #####################################################
            # ITEM 3 - SIMULATION OF SINGLE KNOCKOUT OF REACTION
            #####################################################
            fileLog.write("INICIO DO ITEM 3\n")
            contador = 0
            lbAux = ubAux = 0

            reacoesZerandoFBA = []
            
            file_compound_reaction_ko = open(directory+"//react_biomass_zero_sbml_no_genes.txt", 'w', encoding="ISO-8859-1")
            fileLog.write('len(listReactPerturb) = {0}\n'.format(len(listReactPerturb)))
            for i in listReactPerturb:
                reacao = model.reactions.get_by_id(i)
                contador += 1
                
                #print(reacao)
                #print("Entrei no FOR for i in listReactPerturb:")
               

                if str(model) == "MODEL1507180020":
                    if i == "EX_fe2_e" or i == "FE2abc":
                        #print("reacao dando problema no arquivo da PAO1 2008")
                        #print("-----------------------------------")
                        continue
                
                FBABeforeDelete = model.optimize()
                
                valFBeforeDelete = round(float(FBABeforeDelete.objective_value), numCasasDec)
                if valFBeforeDelete != valInitialSolutionFmt:
                    raise Exception("1o IF - Erro! FBA nao bate com a primeira execucao")
                #print(valFBeforeDelete)
                
                #print(reacao.lower_bound)
                #print(reacao.upper_bound)
                
                lbAux = reacao.lower_bound
                ubAux = reacao.upper_bound
                reacao.lower_bound = 0
                reacao.upper_bound = 0

                #print(reacao.lower_bound)
                #print(reacao.upper_bound)
                
                FBAAfterDelete = model.optimize()
                #print(FBAAfterDelete)

                valFAfterDelete = round(float(FBAAfterDelete.objective_value), numCasasDec)
                #print(valFAfterDelete)
                
                if valFAfterDelete == 0.0 or valFAfterDelete == -0.0:
                    #print("vai gravar no arquivo")
                    #print("{0}\n".format(reacao.build_reaction_string(reacao.reaction)))
                    #file_compound_reaction_ko.write("{0}\n".format(reacao.build_reaction_string(reacao.reaction)))
                    file_compound_reaction_ko.write("{0}\n".format(reacao.build_reaction_string()))
                    reacoesZerandoFBA.append(reacao)

                reacao.upper_bound = ubAux
                reacao.lower_bound = lbAux
                lbAux = 0
                ubAux = 0

                FBAAfterRestoreModel = model.optimize()
                valFFormatAfterRestore = round(float(FBAAfterRestoreModel.objective_value), numCasasDec)
                if valFFormatAfterRestore != valInitialSolutionFmt:
                    raise Exception("2o IF - Erro! FBA nao bate com a primeira execucao")
                
                #print("-------------------------------------------------------------------")
            
            #print("Total de reacoes cujo FBA zerou = " + str(len(reacoesZerandoFBA)))
            fileLog.write('Total de reacoes cujo FBA zerou = {0}\n'.format(str(len(reacoesZerandoFBA))))
            file_compound_reaction_ko.close()
            
            #print("fechou o arquivo")
            #print("{0}".format(len(model.genes)))
            #print("{0}".format(len(model.genes) > 0))
            
            # AQUI ACONTECE UMA BIFURCACAO NA APLICACAO, ONDE ELE FAZ ESSES PASSOS CASO TENHAMOS GENES PARA
            # TRABALHAR, CASO CONTRARIO, ELE PULA ISSO TUDO, FAZ O CENARIO DESCRITO NO ELSE E SEGUE O FLUXO NO PASSO 7
            fileLog.write("INICIO DA BIFURCACAO\n")
            fileLog.write("if len(model.genes) = {0}\n".format(len(model.genes) > 0))
            if len(model.genes) > 0:
                
                #print("Extracao dos genes referentes as reacoes cujo FBA zerou - inicio")
                reacoesZerandoFBAArq = open(directory+"//"+"03-relacao_nocaute_reacao-gene.txt", "w")

                # Aqui pegamos as reacoes e extraimos os genes envolvidos com cada uma
                listaGenesAlvos = []
                for reacao in reacoesZerandoFBA:
                    reacoesZerandoFBAArq.write(str(reacao) + '\n')
                    for gene in reacao.genes:
                        #print(gene)
                        listaGenesAlvos.append(str(gene))
                        reacoesZerandoFBAArq.write(str(gene) + '\n')
                    #print("-----------------")
                    reacoesZerandoFBAArq.write('----------------\n\n')
                reacoesZerandoFBAArq.close()
                
                #print('fechou arquivo reacoesZerandoFBAArq')
                
                # Aqui com a lista de genes gerada a partir das reacoes, ordenamos e colocamos em um arquivo
                listaGenesAlvosSet = list(set(listaGenesAlvos))
                listaGenesAlvosSet.sort()
                arqGenesAlvos = open(directory+"//"+"04-genesAlvos.txt", "w")
                contador = 0
                for gene in listaGenesAlvosSet:
                    #print(gene)
                    contador += 1
                    arqGenesAlvos.write(str(contador) + " - " + str(gene) + "\n")
                arqGenesAlvos.close()
                #print("Total de genes encontrados por reacao = " + str(len(listaGenesAlvosSet)))
                #print("Extracao dos genes referentes as reacoes cujo FBA zerou - final")

                #print("GENERATED FILE %s" % reacoesZerandoFBAArq.name)
                #print("GENERATED FILE %s" % arqGenesAlvos.name)

                ##################################################
                # ITEM 4 - SIMULATION OF SINGLE KNOCKOUT OF GENES
                ##################################################
                fileLog.write("INICIO DO ITEM 4\n")
                fo = open(directory+"//"+"05-relacao_nocaute_gene-reacao.txt", "w")
                contador = 0
                listaGenesAlvosSet02 = []
                for i in model.genes:
                    with model:
                    #print(i)
                        if str(model) == "MODEL1507180020":
                            if str(i) == "PA5259" or str(i) == "PA5260":
                                #print('gene problematico na rede da PAO1 2008')
                                #print('---------------------------------------')
                                continue

                        FBABeforeDelete = model.optimize()
                        valFBeforeDelete = round(float(FBABeforeDelete.objective_value), numCasasDec)
                        if valFBeforeDelete != valInitialSolutionFmt:
                            valInitialSolutionFmt2 = round(float(initialSolution.objective_value), 4)
                            valFBeforeDelete2 = round(float(FBAAfterRestoreModel.objective_value), 4)
                            if valFBeforeDelete2 != valInitialSolutionFmt2:
                                raise Exception("1o IF - Erro! FBA nao bate com a primeira execucao")

                        genes_deleted = cobra.manipulation.knock_out_model_genes(model, [str(i)])
                        print(genes_deleted)

                        FBAAfterDelete = model.optimize()
                        valFAfterDelete = round(float(FBAAfterDelete.objective_value), numCasasDec)
                        if valFAfterDelete == 0.0:
                            contador += 1
                            # reacoes_assoc = (l.reaction for l in i.reactions)
                            fo.write(str(contador) + ". Gene " + str(i) + " inibido.\n")
                            for x in i.reactions:
                                fo.write("%s : %s" % (x.id, x.reaction))
                            # fo.write("Reacoes associadas:\n%s" % ("{\n" + ",\n".join(reacoes_assoc) + "\n}"))
                            fo.write('\n\n')
                            listaGenesAlvosSet02.append(str(i))

                    # cobra.manipulation.undelete_model_genes(model)
                    # cobra.manipulation.

                    # FBAAfterRestoreModel = model.optimize()
                    # valFFormatAfterRestore = round(float(FBAAfterRestoreModel.objective_value), numCasasDec)
                    # if valFFormatAfterRestore != valInitialSolutionFmt:
                    #     valInitialSolutionFmt2 = round(float(initialSolution.objective_value), 4)
                    #     valFFormatAfterRestore2 = round(float(FBAAfterRestoreModel.objective_value), 4)
                    #     if valFFormatAfterRestore2 != valInitialSolutionFmt2:
                    #         raise Exception("2o IF -Erro! FBA nao bate com a primeira execucao")
                    
                    #print('--------------------------------------')
                    
                fo.write("Total de genes inibidos que interrompem a geracao de biomassa = " + str(contador))
                fo.close()

                #print("Total de genes encontrados por nocaute direto = " + str(len(listaGenesAlvosSet02)))
                #print("GENERATED FILE %s" % fo.name)


                ############################################################
                # ITEM 5 - VERIFICATION AND UNIFICATION OF KNOCKOUT RESULTS
                ############################################################
                fileLog.write("INICIO DO ITEM 5\n")
                count = 0
                deParaGenes = open(directory+"//"+"06-deparaGenes.txt", "w")
                for geneAlvo in listaGenesAlvosSet:
                    #print(geneAlvo)
                    if geneAlvo in listaGenesAlvosSet02:
                        count += 1
                        #print("vai gravar")
                        deParaGenes.write(geneAlvo)
                        deParaGenes.write('\n')
                    #print("-------------------------------")
                deParaGenes.close()
                #print("Total de genes encontrados apos de/para = " + str(count))
                #print("GENERATED FILE %s" % deParaGenes.name)


                ##############################################################
                # ITEM 6 - SEARCH FOR THE CORRESPONDING EC NUMBER OF THE GENE
                ##############################################################
                fileLog.write("INICIO DO ITEM 6\n")
                k = KEGG(cache=True)
                k.TIMEOUT = TIMEOUT_SECONDS

                contador = 0
                with open(directory+"/"+"06-deparaGenes.txt", "r") as infile:
                    data = infile.read()
                my_list_genes = data.splitlines()

                fileGenesWithEC = open(directory+"/"+"07-genes_ECNumbers.txt", "w")
                fileAssocGenesEC = open(directory+"/"+"07-1-assoc_genes_ECNumbers.txt", "w")

                # AQUI ELE BUSCA TODOS OS ACRONIMOS DO KEGG RELACIONADOS COM O ORGANISMO SELECIONADO NA TELA
                df_list_organism = pd.read_excel(static_dir + "/list_organism_kegg_2022_11_07.xls")
                df_list_organism_filter = df_list_organism[df_list_organism['name_organism'].str.contains(organismParam)]
                list_acron_kegg = df_list_organism_filter['acron_organism'].tolist()
                fileLog.write('list_acron_kegg = {0}\n'.format(list_acron_kegg))
                fileLog.write('len my_list_genes = {0}\n'.format(len(my_list_genes)))
                
                for acron in list_acron_kegg:
                    for gene in my_list_genes:
                        source = acron+":"+gene
                        #fileLog.write('source = {0}\n'.format(source))
                        #source = "pae:"+gene # pae:PA0005
                        ec = k.link("enzyme", source)
                        #fileLog.write('ec = {0}\n'.format(ec))
                        
                        returnKeggEC = str(ec).split() # split por espaco para separar o que veio do unicode
                        if (len(returnKeggEC)) == 0 :
                            continue
                        contador += 1

                        for ecnumber in returnKeggEC:
                            if "ec:" in ecnumber: # caso tenha o inicio "ec:"
                                ecnumbersplitFinal = ecnumber.split(":") # novo split para prevalecer apenas o numero
                                fileGenesWithEC.write(ecnumbersplitFinal[-1])
                                fileGenesWithEC.write('\n')
                                fileAssocGenesEC.write('{0};{1}\n'.format(gene, ecnumbersplitFinal[-1]))

                fileGenesWithEC.close()
                fileAssocGenesEC.close()

                #print("GENERATED FILE %s" % fileGenesWithEC.name)
                #print("GENERATED FILE %s" % fileAssocGenesEC.name)

            else:
                fileLog.write("CAIU NO ELSE DE REDE SEM GENES\n")
                self.alternativeStepToGetECNumberWithoutGenes(directory)
            
            
            #######################################################
            # ITEM 7 - SEARCH FOR PROTEIN BY EC NUMBER NO DRUGBANK
            #######################################################
            fileLog.write("INICIO DO ITEM 7\n")
            with open(directory+"//"+"07-genes_ECNumbers.txt", "r") as infile:
                data = infile.read()
            my_list_ecnumbers = data.splitlines()
            my_list_ecnumbers = list(set(my_list_ecnumbers))
            my_list_ecnumbers.sort()

            # my_list_ecnumbers = ['2.5.1.7']

            # Aqui filtra os ECs encontrados no drugbank e simula a navegacao na pagina
            fileFilterEC = open(directory+"//"+"08-filter_ECNumbers_drugbank.txt", "w")

            for ec in my_list_ecnumbers:
                #print(ec)
                time.sleep(1) # esperar 1 segundo para acessar o link
                
                link01 = "https://www.drugbank.ca/unearth/q?utf8=%E2%9C%93&searcher=bio_entities&query="+str(ec)+"&approved=1&nutraceutical=1&illicit=1&investigational=1&withdrawn=1&experimental=1&us=0&ca=0&eu=0&commit=Apply+Filter"
                #link01 = "https://www.drugbank.ca/unearth/q?utf8=%E2%9C%93&searcher=targets&query="+str(ec)+"&approved=1&illicit=1&investigational=1&withdrawn=1&experimental=1&us=1&canada=1&eu=1&commit=Apply+Filter"
                r = requests.get(link01)

                # verificar se houve retorno da pagina acessada
                if r.status_code == 200:

                    # Instancia de bs4 do resultado da 1a tela do drugbank
                    soup = BeautifulSoup(r.text)

                    # filtro para a lista de todos os itens encontrados na busca do requests
                    list_found_ec = soup.findAll(attrs={"class": "search-result p-2 mb-4 p-sm-3 mb-sm-2"})

                    # caso tenha item, inicia as buscas para pegar o resto da informacao
                    if len(list_found_ec) > 0:

                        # itera na lista de itens encontrados no drugbank
                        for item in list_found_ec:
                            em = item.find('em')
                            ec_encontrado = em.text

                            #print(ec_encontrado, ec, ec_encontrado == ec)
                            if ec_encontrado != ec:
                                continue

                            # busca pelo link para acesso a 2a tela do drugbank
                            var = item.find('a')

                            # Montagem do link e acesso a 2a tela do drugbank para buscar as informacoes necessarias
                            link02 = "https://www.drugbank.ca"+str(var['href'])
                            r2 = requests.get(link02)
                            if r2.status_code == 200:

                                # Nova instancia de bs4 com as informacoes da 2a tela
                                # e pegando nome da proteina, organismo, uniprot id e drugbank id
                                soup2 = BeautifulSoup(r2.text)
                                organism_data_class = soup2.find(attrs={"class": "content-container"})
                                protein_name = organism_data_class.findAll('dd')[0].text
                                organism_name = organism_data_class.findAll('dd')[2].text
                                uniprot_id = organism_data_class.find('a').text
                                drugbank_id = str(var['href']).split("/")[-1]
                                fileFilterEC.write('{0};{1};{2};{3};{4}\n'.format(ec, protein_name, organism_name, uniprot_id, drugbank_id))

                            else:
                                fileFilterEC.write('{0};{1}\n'.format(ec, r.raise_for_status()))

                else:
                    fileFilterEC.write('{0};{1}\n'.format(ec, r.raise_for_status()))

            fileFilterEC.close()

            #print("GENERATED FILE %s" % fileFilterEC.name)
            
            ################################################
            # ITEM 8 - SEARCH FOR THE HOMOLOGUES AT UNIPROT
            ################################################
            fileLog.write("INICIO DO ITEM 8\n")
            with open(directory+"//"+"08-filter_ECNumbers_drugbank.txt", "r") as infile:
                data = infile.read()
            my_list_uniprotid = data.splitlines()
            my_list_uniprotid = list(set(my_list_uniprotid))
            my_list_uniprotid.sort()

            # Instancia de chamada ao UniProt
            u = UniProt()
            s = NCBIblast()

            u.TIMEOUT = TIMEOUT_SECONDS
            s.TIMEOUT = TIMEOUT_SECONDS

            #my_list_uniprotid = ["P42898"]
            
            # GET ALL JOBID FROM UNIPROT FOUND
            file_jobid = open(directory+"//"+"list_jobid.txt", "w", encoding="ISO-8859-1")
            start_time_first_req = time.time()
            for uniprotid in my_list_uniprotid:

                uniprotid = uniprotid.split(';')[3]
                #print(uniprotid)
                
                # buscando a sequencia da proteina por ID, vindo do drugbank
                findFasta = u.retrieve(uniprotid, "fasta")
                sequence = ""
                if isinstance(findFasta, str):
                    if (findFasta != ""):
                        sequence = findFasta
                    else:
                        continue
                
                elif isinstance(findFasta, Response):
                    if findFasta.ok:
                        sequence = findFasta.text
                    else:
                        continue
                    
                else:
                    continue
                
                sequence = sequence.split("\n", 1)[1]
                sequence = re.sub(r"\n", r"", sequence)
                
                # executando o job que faz o blast
                jobid = s.run(program="blastp", database="uniprotkb_swissprot", sequence=sequence, stype="protein", email="findtargetweb@gmail.com", alignments='1000')
                file_jobid.write("{0};{1}\n".format(uniprotid, jobid))
                
                time.sleep(1)
                
                #print("----------------------------------------------")

            elapsed_time = time.time() - start_time_first_req
            #print("TEMPO TOTAL PRIMEIRO PASSO = ", time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
            file_jobid.close()
            
            fileResultBlast = open(directory+"//"+"11-hitsEncontradosUniprot.txt", "w")
            
            # AFTER GENERATE ALL JOBID FROM BLAST, WE ITERATE ALL ITENS AND GET ALL XML RETURN
            with open(directory+"//"+"list_jobid.txt", "r") as infile:
                data = infile.read()
            my_list_jobid = data.splitlines()
            #my_list_jobid = list(set(my_list_jobid))

            arquivo = open(directory+"//"+"list_jobid.txt", 'r')

            for item in arquivo:
                #print(item)
                item = item.strip()
                uniprotid = item.split(';')[0]
                jobid = item.split(';')[1]

                url_status_blast = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/"
                url_result_blast = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/"
                last_url_result_blast = "/xml"

                start_time_first_req = time.time()

                condition = True
                count_error = 0
                url_1 = url_status_blast+jobid
                #print(url_1)
                count_refresh = 0
                while condition:
                    response = requests.get(url_1)
                    #print(response.text)
                    if response.text != "RUNNING":
                        if response.text == "FINISHED":
                            condition = False

                        if response.text == "ERROR" or response.text == "NOT_FOUND":
                            count_error += 1
                            #print("ERRO!", count_error)
                            if count_error == 3:
                                condition = False

                    time.sleep(3)
                    if count_refresh > 60:
                        condition = False

                    count_refresh += 1

                elapsed_time = time.time() - start_time_first_req
                #print("TEMPO TOTAL PRIMEIRO LINK = ", time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

                start_time_second_req = time.time()

                if response.text != 'FINISHED':
                    #fileResultBlast.write("{0}--{1}--{2}--{3}--{4}--{5}--{6}--{7}--{8}--{9}\n".format(
                    #    uniprotid, "ERRO", "ERRO", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                    #))
                    continue

                url_2 = url_result_blast+jobid+last_url_result_blast
                #print(url_2)
                response = requests.get(url_2)
                soup = BeautifulSoup(response.text)

                elapsed_time_2 = time.time() - start_time_second_req
                #print("TEMPO TOTAL SEGUNDO LINK = ", time.strftime("%H:%M:%S", time.gmtime(elapsed_time_2)))

                # ARMAZENO EM OUTRA PASTA OS ARQUIVOS DO BLAST NA INTEGRA
                #fileHitsBlastIntegra = open(dir_blasts+"//all_hits_"+uniprotid+".xml", "w")
                #fileHitsBlastIntegra.write(str(soup))
                #fileHitsBlastIntegra.close()

                # busca por todos os hits encontrados
                listAllHits = soup.findAll('hit')

                # itera os hits e, caso encontre pseudomonas, guarda as informacoes dele
                fileListAllHist = open(dir_blasts+"//hit_organism_found_"+uniprotid+".txt", 'w')

                percentSimilarHuman = 0.0
                for hit in listAllHits:
                    organism = hit['description'].split("=")[1].split("GN")[0].strip()
                    if "Homo sapiens" in organism:
                        percentSimilarHuman = float(hit.find('identity').string.strip())
                        break

                for hit in listAllHits:
                    organism = hit['description'].split("=")[1].split("GN")[0].strip()
                    fileListAllHist.write("{0}\n".format(organism))

                    #if organism == 'Pseudomonas aeruginosa (strain ATCC 15692 / DSM 22644 / CIP 104116 / JCM 14847 / LMG 12228 / 1C / PRS 101 / PAO1)' :
                    if organismParam in organism: # se o valor da combo for encontrado no hit, armazena
                        idPAO1 = hit['ac']
                        percentBlast = hit.find('identity').string.strip()
                        eValue = hit.find('expectation').string.strip()

                        #print(idPAO1)
                        # ###### TESTE ----- COMMENT ----- (SUBCELLULAR LOCATION) - deprecated - substituir por 'subcellular locations'
                        foundID = u.search(idPAO1)
                        # foundID = u.search(idPAO1, columns='entry name, id, genes, pathway, comment(FUNCTION), comment(CATALYTIC ACTIVITY), database(PDB), database(PSEUDOCAP), subcellular locations, ec')
                        foundIDResponse = ""
                        if isinstance(foundID, str):
                            if (foundID != ""):
                                foundIDResponse = foundID
                            else:
                                continue
                        else:
                            if foundID != 400:
                                foundIDResponse = foundID.text
                            else:
                                continue
                        
                        #print(foundIDResponse)
                        
                        #splitForListFoundID = list(str(foundID).split("\t"))
                        splitForListFoundID = list(str(foundIDResponse).split("\t"))
                        
                        #print(len(splitForListFoundID))
                        if len(splitForListFoundID) == 0:
                            continue 

                        # 0-uniprotid;;1-hit pseudomonas;;2-percentidentityblast;;3-evalue;;
                        # 4-Gene names;;5-Pathway;;6-Function;;7-Catalitic Activity;;
                        # 8-Localization;;9-ID PDB
                        #if float(percentSimilarHuman) < float(percentBlast): # so grava se o percentual de similaridade for menor do que com humanos
                        #    fileResultBlast.write("{0};;{1};;{2};;{3};;{4};;{5};;{6};;{7};;{8};;{9}\n".format(
                        #        uniprotid, idPAO1, percentBlast, eValue, 
                        #					splitForListFoundID[10], splitForListFoundID[11],
                        #        splitForListFoundID[12], splitForListFoundID[13], splitForListFoundID[16], splitForListFoundID[14]
                        #    ))
                        # c10 = splitForListFoundID[10]
                        # c11 = splitForListFoundID[11]
                        # c12 = splitForListFoundID[12]
                        # c13 = splitForListFoundID[13]
                        # c14 = splitForListFoundID[14]

                        string_resultante = ""
                        for valor in splitForListFoundID:
                            string_resultante += ";" + str(valor)

                        if float(percentSimilarHuman) < float(percentBlast): # so grava se o percentual de similaridade for menor do que com humanos
                            fileResultBlast.write("{0};{1};{2};{3}{4}\n".format(
                                uniprotid, idPAO1, percentBlast, eValue, string_resultante
                            ))

                        else: # Caso tenha hit com percent de humano maior, poe tudo zerado
                            fileResultBlast.write("{0};{1};{2};{3};{4};{5};{6};{7};{8};{9}\n".format(
                                uniprotid, idPAO1, percentBlast, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                            ))

                fileListAllHist.close()

            fileResultBlast.close()
            arquivo.close()

            ##print("GENERATED FILE %s" % fileListAllHist.name)
            #print("GENERATED FILE %s" % fileResultBlast.name)

            #####################################################
            # ITEM 9 - SEARCH FOR THE INHIBITORS FOR THE BACTERIA
            #####################################################
            fileLog.write("INICIO DO ITEM 9\n")
            # PEGANDO A ENTRADA DO ITEM POR MEIO DOS ARQUIVOS GERADOS E CONVERTENDO A DATAFRAME PARA FACILITAR
            my_list_uniprotid_drugbankid = []
            my_list_uniprotid_hit = []
            
            #data08 = pd.read_csv(directory + "//" + "08-filter_ECNumbers_drugbank.txt", sep=";", header=None)
            data08 = pd.DataFrame()
            try:
                data08 = pd.read_csv(directory + "//" + "08-filter_ECNumbers_drugbank.txt", sep=";", header=None)
            except:
                print('Note: file 08-filter... was empty. Skipping.')

            if not data08.empty:
                data08.columns = ["EC NUMBER", "PRODUCT", "ORGANISM NAME", "UNIPROTID", "DRUGBANKID"]
                data08 = data08.sort_values("EC NUMBER")
                my_list_uniprotid_drugbankid = data08["UNIPROTID"].tolist()
                data08.to_excel(directory + "//" + "08-filter_ECNumbers_drugbank.xls")

            data11 = pd.DataFrame()
            try:
                data11 = pd.read_csv(directory + "//" + "11-hitsEncontradosUniprot.txt", sep=";", header=None)
            except:
                print('Note: file 11-hits... was empty. Skipping.')
            
            #data11 = pd.read_csv(directory + "//" + "11-hitsEncontradosUniprot.txt", sep="--", header=None)
            if not data11.empty:
                data11.columns = ["UNIPROTID", "HIT_UNIPROTID", "PERCENT_IDENT_BLAST", "EVALUE", "GENE_NAME",
                                  "PATHWAY", "FUNCTION", "CATALYTIC ACTIVITY", "LOCALIZATION", "ID PDB"]
                data11 = data11.sort_values("UNIPROTID")
                my_list_uniprotid_hit = data11["UNIPROTID"].tolist()
                data11.to_excel(directory + "//" + "11-hits_Uniprot.xls")

            # GERACAO DOS ARQUIVOS DE SAIDA DESSE ITEM
            fileDataDrugs = open(directory + "//" + "13-list_inhibitors_per_target.txt", "w")
            fileInhibitorsDrugs = open(directory + "//" + "14-list_inhibitors_approved.txt", "w")

            for uniprot_drugbank in my_list_uniprotid_drugbankid:

                if uniprot_drugbank in my_list_uniprotid_hit:
                    data11Aux = data11[data11['UNIPROTID'] == uniprot_drugbank]

                    hit_toxicity = False
                    for item in data11Aux.itertuples():
                        validate = item.GENE_NAME
                        if validate == '0.0':
                            hit_toxicity = True
                            break

                    # CASO TENHA HIT MENOR QUE HUMANO, ESTARA CHEIO DE ZEROS NO RESULTADO
                    # SENDO ASSIM, ELE NAO DEVE SER CONSIDERADO UM BOM RESULTADO E NAO DEVE SER ARMAZENADO
                    if hit_toxicity:
                        continue

                    #print(uniprot_drugbank)
                    data08Aux = data08[data08["UNIPROTID"] == uniprot_drugbank]
                    drugbankid = data08Aux.iloc[0][4]

                    linkAccessForDrugsID = "https://www.drugbank.ca/biodb/bio_entities/" + str(drugbankid)
                    r = requests.get(linkAccessForDrugsID)

                    if r.status_code == 200:
                        soup = BeautifulSoup(r.text)
                        listDrugsFound = soup.find(attrs={"class": "table table-sm table-bordered datatable dt-responsive"})
                        listDrugsFound = listDrugsFound.find('tbody')
                        listDrugsFound = listDrugsFound.findAll('tr')  # lista com todos os farmacos para o drugbankid informado

                        for item in listDrugsFound:
                            drugbank_drug_id = item('td')[0].text
                            drug_name = item('td')[1].text
                            drug_group = item('td')[2].text
                            pharma_action = item('td')[3].text
                            actions = item('td')[4].text

                            # 0-uniprotid;;1-drugbankid;;2-drugbank_drug_id;;3-drug_name;;
                            # 4-drug_group;;5-pharma_action;;6-actions
                            fileDataDrugs.write("{0};{1};{2};{3};{4};{5};{6}\n".format(
                                uniprot_drugbank, drugbankid, drugbank_drug_id, drug_name,
                                drug_group, pharma_action, actions
                            ))

                            if pharma_action == 'yes' and actions == 'inhibitor':
                                # 0-uniprotid;;1-drugbankid;;2-drugbank_drug_id;;3-drug_name;;
                                # 4-drug_group;;5-pharma_action;;6-actions
                                fileInhibitorsDrugs.write("{0};{1};{2};{3};{4};{5};{6}\n".format(
                                    uniprot_drugbank, drugbankid, drugbank_drug_id, drug_name,
                                    drug_group, pharma_action, actions
                                ))

                    else:
                        continue
                        '''
                        fileDataDrugs.write("{0};{1};{2};{3};{4};{5};{6}\n".format(
                            uniprot_drugbank, drugbankid, r.raise_for_status(), r.raise_for_status(),
                            r.raise_for_status(), r.raise_for_status(), r.raise_for_status()
                        ))
                        '''

            fileDataDrugs.close()
            fileInhibitorsDrugs.close()

            data13 = pd.DataFrame()
            try:
                data13 = pd.read_csv(directory + "//" + "13-list_inhibitors_per_target.txt", sep=";", header=None)
            except:
                print('Note: file 13-list_inhi... was empty. Skipping.')

            #data13 = pd.read_csv(directory + "//" + "13-list_inhibitors_per_target.txt", sep=";", header=None)
            if not data13.empty:
                data13.columns = ["UNIPROTID", "DRUGBANKID", "DRUGBANKDRUGID", "DRUGNAME", "DRUGGROUP", "PHARMAACTION", "ACTIONS"]
                data13 = data13.sort_values("UNIPROTID")
                data13.to_excel(directory + "//" + "13-list_inhibitors_per_target.xls")

            data14 = pd.DataFrame()
            try:
                data14 = pd.read_csv(directory + "//" + "14-list_inhibitors_approved.txt", sep=";", header=None)
            except:
                print('Note: file 14-list_inhi... was empty. Skipping.')

            #data14 = pd.read_csv(directory + "//" + "14-list_inhibitors_approved.txt", sep=";", header=None)
            if not data14.empty:
                data14.columns = ["UNIPROTID", "DRUGBANKID", "DRUGBANKDRUGID", "DRUGNAME", "DRUGGROUP", "PHARMAACTION", "ACTIONS"]
                data14 = data14.sort_values("UNIPROTID")
                data14.to_excel(directory + "//" + "14-list_inhibitors_approved.xls")

            #print("GENERATED FILE %s" % fileDataDrugs.name)
            #print("GENERATED FILE %s" % fileInhibitorsDrugs.name)

            # CRIANDO A PLANILHA UNICA COM TODOS OS RESULTADOS JUNTOS
            # CONTENDO APENAS AS OCORRENCIAS EM TODAS AS 3 PLANILHAS GERADAS!
            try:
                df_08 = pd.read_excel(directory + "//" + "08-filter_ECNumbers_drugbank.xls")
                df_11 = pd.read_excel(directory + "//" + "11-hits_Uniprot.xls")
                df_13 = pd.read_excel(directory + "//" + "13-list_inhibitors_per_target.xls")
    
                df_merged = pd.merge(pd.merge(df_08, df_11, on='UNIPROTID', how='inner'), df_13, on=['UNIPROTID','DRUGBANKID'], how='inner').drop_duplicates().drop(columns=['LOCALIZATION']).reset_index(drop=True)
                df_merged.to_excel(directory + "//" + "summary_results.xls")
            except:
                print("UM DOS ARQUIVOS ESTA VAZIO. ARQUIVO NAO GERADO.")

            ##################################################################################
            # ITEM 10 - LAST ITENS FOR THE METHOD WHERE ZIP ALL DATAS GENERATED AND SEND MAIL
            ##################################################################################
            fileLog.write("INICIO DO ITEM 10\n")
            zip_file_report = self.zipReportsToSendMail(dir_path, dataHoraAtualFmt, directory)
            self.sendMailWithReportAttached(name, email, zip_file_report, model, method)
            fileLog.write("ACABOU!!! SEJA FELIZ!!!\n")

        except Exception as e:
            self.sendMailWithError(name, email, str(traceback.format_exception(None, e, e.__traceback__)))

        end = time.time()
        time_second = end - start
        time_minutes = time_second / 60
        time_hour = time_minutes / 60
        
        #print("TEMPO TOTAL DE EXECUCAO = {0} segundos".format(time_second))
        #print("TEMPO TOTAL DE EXECUCAO = {0} minutos".format(time_minutes))
        #print("TEMPO TOTAL DE EXECUCAO = {0} horas".format(time_hour))
        fileLog.close()

    # GENERATE ZIP FILE TO SEND MAIL TO USER
    def zipReportsToSendMail(self, dir_path, dataHoraAtualFmt, directory):
        zf = ZipFile(dir_path+"/results/"+dataHoraAtualFmt+'_results.zip', "w")
        for root, subdirs, files in os.walk(directory):
            for filename in files:
                file_extension = filename.split(".")[-1]
                if file_extension == 'xls':
                    zf.write(os.path.join(root, filename), filename, ZIP_DEFLATED)
        zf.close()

        return zf

    # SEND MAIL WITH RESULTS
    def sendMailWithReportAttached(self, name, email, zip_file_report, nameFile, method):
        methodSelected = ""
        if method == "1":
            methodSelected = "FBA+FVA"
        else:
            methodSelected = "Only FBA"
        
        if nameFile == "":
            nameFile = "NAME NOT FOUND"
        
        mensagem = "FINAL REPORT - FIND TARGETS WEB \n\n\n" \
                   "" \
                   "FILE: " + str(nameFile) + " | METHOD: " + str(methodSelected) +"\n" \
                   "" \
                   "YOUR ANALYSIS HAS FINISHED SUCESSFULLY. THE FOLLOWING FILES HAVE BEEN GENERATED:\n" \
                   "" \
                   "08-filter_ECNumbers_drugbank.xls\n" \
                   "FIELDS: EC NUMBER, PRODUCT, ORGANISM NAME, UNIPROTID, DRUGBANKID\n\n" \
                   "" \
                   "11-hits_Uniprot.xls\n" \
                   "FIELDS: UNIPROTID, HIT_UNIPROTID, PERCENT_IDENT_BLAST, EVALUE, GENE_NAME, PATHWAY, FUNCTION, CATALYTIC ACTIVITY, LOCALIZATION, ID PDB\n\n" \
                   "" \
                   "13-list_inhibitors_per_target.xls AND 14-list_inhibitors_approved.xls\n" \
                   "FIELDS: UNIPROTID	DRUGBANKID	DRUGBANKDRUGID	DRUGNAME	DRUGGROUP	PHARMAACTION	ACTIONS\n\n" \
                   "" \
                   "model_data.xls\n" \
                   "ALL GENES, REACTIONS AND METABOLITES IN THE SBML FILE USED IN ANALYSIS\n\n" \
                   "" \
                   "summary_results.xls\n" \
                   "FILE WITH A SUMMARY OF DATA FOUND IN 08, 11 AND 13 XLS FILES.\n\n" \

        # remetente    = 'findtargetsweb_fiocruz@hotmail.com'
        remetente = 'findtargetweb@gmail.com'
        # senha        = 'proccfiocruz1234'
        # senha = 'whzctyqjvqsxojcz'
        senha = 'ntucaeevfacubqyj'

        # Informacoes da mensagem
        destinatario = email
        assunto      = 'REPORT THERAPEUTICS TARGETS FROM YOUR NETWORK MODEL'
        
        msg = MIMEMultipart()
 
        msg['From'] = remetente
        msg['To'] = destinatario
        msg['Subject'] = assunto
         
        # Preparando a mensagem
        '''
        mensagem = '\r\n'.join([
          'From: %s' % remetente,
          'To: %s' % destinatario,
          'Subject: %s' % assunto,
          '',
          '%s' % mensagem
          ])
        '''
        mensagem = '\r\n'.join([
            '%s' % mensagem
        ])

        mensagem = mensagem.encode("UTF-8")
        
        msg.attach(MIMEText(mensagem.decode("UTF-8"), 'plain'))
        
        filename = zip_file_report.filename
        attachment = open(filename, "rb")
         
        part = MIMEBase('application', 'zip')
        part.set_payload((attachment).read())
        encoders.encode_base64(part)
        part.add_header('Content-Disposition', "attachment", filename=os.path.basename(filename))
        msg.attach(part)
         
        # Enviando o email (USANDO O SMTP DO HOTMAIL PARA ENVIAR)
        # server = smtplib.SMTP("smtp.live.com: 587")
        server = smtplib.SMTP("smtp.gmail.com: 587")
        server.starttls()
        server.login(remetente,senha)
        text = msg.as_string()
        server.sendmail(remetente, destinatario, text)
        server.quit()


    # SEND MAIL WITH ERRORS!
    def sendMailWithError(self, name, email, dsc_exception):
        mensagem = "ERROR! PLEASE CONTACT ADMINISTRATOR OR TRY AGAIN\n\n"
        mensagem = mensagem + dsc_exception

        # remetente    = 'findtargetsweb_fiocruz@hotmail.com'
        remetente    = 'findtargetweb@gmail.com'
        # senha        = 'proccfiocruz1234'
        # senha        = 'whzctyqjvqsxojcz'
        senha = 'ntucaeevfacubqyj'

        # Informacoes da mensagem
        destinatario = email
        assunto = 'REPORT THERAPEUTICS TARGETS FROM YOUR NETWORK MODEL - ERROR!'

        msg = MIMEMultipart()

        msg['From'] = remetente
        msg['To'] = destinatario
        msg['Subject'] = assunto

        # Preparando a mensagem
        mensagem = '\r\n'.join([
            '%s' % mensagem
        ])

        mensagem = mensagem.encode("UTF-8")

        msg.attach(MIMEText(mensagem.decode("UTF-8"), 'plain'))

        # Enviando o email (USANDO O SMTP DO GMAIL PARA ENVIAR)
        #server = smtplib.SMTP("smtp.live.com: 587")
        server = smtplib.SMTP("smtp.gmail.com: 587")
        server.starttls()
        server.login(remetente, senha)
        text = msg.as_string()
        server.sendmail(remetente, destinatario, text)
        server.quit()


    # METODO QUE ARMAZENA EM PLANILHA OS DADOS DE GENE E REACAO DE UM MODELO AVALIADO
    # @PARAMS => model = modelo | nomeArquivoModelo = nome do arquivo sbml
    def reportModel(self, model, nomeArquivoModelo):
        workbook_dados_modelo = xlwt.Workbook()
        sheet_genes = workbook_dados_modelo.add_sheet("GENES")
        sheet_genes.write(0, 0, "ID GENE")
        sheet_genes.write(0, 1, "GENE")
        sheet_genes.write(0, 2, "REACTIONS")
        contador_genes = 1
        for gene in model.genes:
            assoc = (l.id for l in gene.reactions)
            sheet_genes.write(contador_genes, 0, gene.id)
            sheet_genes.write(contador_genes, 1, gene.name)
            sheet_genes.write(contador_genes, 2, ", ".join(assoc))
            contador_genes += 1
        
        sheet_react = workbook_dados_modelo.add_sheet("REACTIONS")
        sheet_react.write(0, 0, "ID REACTION")
        sheet_react.write(0, 1, "NAME")
        sheet_react.write(0, 2, "COMPOSITION")
        sheet_react.write(0, 3, "LIMITS INF/SUP")
        sheet_react.write(0, 4, "SUBSYSTEM")
        sheet_react.write(0, 5, "GENES")
        sheet_react.write(0, 6, "METABOLITES")
        contador_react = 1
        
        for reacao in model.reactions:
            assoc_genes = (l.id for l in reacao.genes)
            assoc_metab = (l.id for l in reacao.metabolites)
            sheet_react.write(contador_react, 0, reacao.id)
            sheet_react.write(contador_react, 1, reacao.name)
            sheet_react.write(contador_react, 2, reacao.build_reaction_string(reacao.reaction))
            sheet_react.write(contador_react, 3, str(reacao.bounds))
            sheet_react.write(contador_react, 4, reacao.subsystem)
            sheet_react.write(contador_react, 5, ", ".join(assoc_genes))
            sheet_react.write(contador_react, 6, ", ".join(assoc_metab))
            contador_react += 1
            
        workbook_dados_modelo.save(str(nomeArquivoModelo))
        
    
    
    # METHOD TO GET EC NUMBER FROM SBML WITHOUT GENES MAPPED
    def alternativeStepToGetECNumberWithoutGenes(self, directory):
        with open(directory+"//react_biomass_zero_sbml_no_genes.txt", "r") as infile:
            data = infile.read()
        my_list_compound_sbml = data.splitlines()
        
        k = KEGG(False, True)
        k.TIMEOUT = 500000
        pd.options.display.max_colwidth = 1000
        pd.options.display.max_rows = 1000
        
        file_ecs = open(directory+"//07-genes_ECNumbers.txt", "w", encoding="ISO-8859-1")
        file_ecs_compound = open(directory + "//07-1-assoc_EC_compounds.txt", "w", encoding="ISO-8859-1")
        file_comp_not_found_from_reactant_sbml = open(directory+"//idcomp_notfound_kegg.txt", "w", encoding="ISO-8859-1")
        file_reaction_not_found_from_reactant_sbml = open(directory+"//idreaction_not_found_kegg.txt", "w", encoding="ISO-8859-1")
        
        #print("PARTE 1 METODO ALTERNATIVO")
        for compound in my_list_compound_sbml:
            #print("parte1", compound)
            compound_no_stoich = re.sub("\d+\.\d+", "", compound) # retirada de todos os valores estequiometricos
            
            if " - reduced " in compound:
                compound_no_stoich = compound_no_stoich.replace(" - reduced ", " ")
            
            param_splt = "" # preparo do split para separar a composicao quimica
            if "<=>" in compound_no_stoich:
                param_splt = "<=>"
            else:
                param_splt = "-->"
            
            compound_splt = compound_no_stoich.split(param_splt) # separacao da composicao quimica pelo reactante(0) e o produto(1) 
            
            # Caso seja uma reacao de insercao ou uma reacao de excrecao, ele ignora
            if compound_splt[0].strip() == "" or compound_splt[1].strip() == "":
                continue  
            
            # TRATAMENTO DOS DADOS DO REAGENTE DA COMPOSICAO ENCONTRADA
            reactant_sbml = compound_splt[0].strip()
        
            list_id_cpd_kegg = []
            len_reactant_sbml = 0
            
            # verifica os reagentes que possuem mais de um composto envolvido
            if " + " in reactant_sbml: #(ex.: acetyl-coa + atp + bicarbonate)
                reactant_sbml_splt = reactant_sbml.split(" + ")
                len_reactant_sbml = len(reactant_sbml_splt)
                
                # iteracao dentro dos compostos dos reactants (acetyl-coa + atp + bicarbonate)
                for cpd_reactant_sbml in reactant_sbml_splt:
                    #http://rest.kegg.jp/find/compound/acetyl-coa
                    result_id_cpd = k.find("compound", cpd_reactant_sbml)
                    result_id_cpd_splt = result_id_cpd.split('\n')
                    result_id_cpd_splt = filter(None, result_id_cpd_splt)
                    result_id_cpd_splt_filter = list(result_id_cpd_splt)
                    
                    # iteracao dentro do resultado encontrado para um dos compostos do reagente 
                    # (cpd:C00024 Acetyl-CoA; Acetyl coenzyme A)
                    if len(result_id_cpd_splt_filter) > 0:
                        for result_cpd in result_id_cpd_splt_filter:
                            local_result_cpd_splt = result_cpd.split('\t')
                            if ";" in local_result_cpd_splt[1]:
                                dsc_cpd_splt = local_result_cpd_splt[1].split(';')
                                for cpd_splt in dsc_cpd_splt:
                                    if cpd_splt.strip().lower() == cpd_reactant_sbml.strip().lower():
                                        list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                        break
                            else:
                                if cpd_reactant_sbml.lower() == local_result_cpd_splt[1].lower():
                                    list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                    break
            else:
                len_reactant_sbml = 1
                result_id_cpd = k.find("compound", reactant_sbml)
                result_id_cpd_splt = result_id_cpd.split('\n')
                result_id_cpd_splt = filter(None, result_id_cpd_splt)
                result_id_cpd_splt_filter = list(result_id_cpd_splt)
                
                if len(result_id_cpd_splt_filter) > 0:
                    for result_cpd in result_id_cpd_splt_filter:
                        local_result_cpd_splt = result_cpd.split('\t')
                        if ";" in local_result_cpd_splt[1]:
                            dsc_cpd_splt = local_result_cpd_splt[1].split(';')
                            for cpd_splt in dsc_cpd_splt:
                                if cpd_splt.strip().lower() == reactant_sbml.strip().lower():
                                    list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                    break
                        else:
                            if reactant_sbml.strip().lower() == local_result_cpd_splt[1].strip().lower():
                                list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                break
        
            # caso nao tenha encontrado todos os compostos, ele nao segue o fluxo
            if len(list_id_cpd_kegg) != len_reactant_sbml:
                file_comp_not_found_from_reactant_sbml.write("{0}\n".format(compound))
                continue
            
            # PREPARO DA APLICACAO PARA USAR OS IDS DE COMPOSTOS E BUSCAR PELAS REACOES DE CADA UM
            reactant_sbml_in_cpd = "+".join(list_id_cpd_kegg)
            result_link_reactions_cpd = k.link("reaction", str(reactant_sbml_in_cpd))
            result_link_reactions_cpd_splt = result_link_reactions_cpd.split('\n')
            result_link_reactions_cpd_splt = filter(None, result_link_reactions_cpd_splt)
            result_link_reactions_cpd_splt_filter = list(result_link_reactions_cpd_splt)
            df_link_reaction_cpd = pd.DataFrame(columns=['id_cpd', 'id_reaction'])
            
            index_link_react_cpd = 0
            if len(result_link_reactions_cpd_splt_filter) > 0:
                for item_result_link_reactions_cpd in result_link_reactions_cpd_splt_filter:
                    local_item_result_link = item_result_link_reactions_cpd.split('\t')
                    df_link_reaction_cpd.loc[index_link_react_cpd] = [local_item_result_link[0], local_item_result_link[1]]
                    index_link_react_cpd += 1
            else:
                file_reaction_not_found_from_reactant_sbml.write("{0}\n".format(compound))
                continue
            
            df_link_reaction_cpd = df_link_reaction_cpd.groupby("id_reaction").filter(lambda x: len(x) == len(reactant_sbml_splt))
            set_id_reaction_kegg = {x.id_reaction for x in df_link_reaction_cpd.itertuples()}
            
            if len(set_id_reaction_kegg) > 0:
                for item_id_react in set_id_reaction_kegg:
                    result_ec_number = k.link("enzyme", item_id_react)
                    result_ec_number_splt = result_ec_number.split('\n')
                    result_ec_number_splt = filter(None, result_ec_number_splt)
                    result_ec_number_filter = list(result_ec_number_splt)
                    
                    for result_ec in result_ec_number_filter:
                        local_result_ec = result_ec.split('\t')
                        local_result_ec[1] = local_result_ec[1].replace("ec:", "").strip()
                        file_ecs.write("{0}\n".format(local_result_ec[1]))
            
        file_comp_not_found_from_reactant_sbml.close()
        file_reaction_not_found_from_reactant_sbml.close()
        
        # AQUI ELE INICIA A BUSCA POR PRODUTO DO QUE NAO FOI ENCONTRADO POR REAGENTE
        with open(directory+"//idcomp_notfound_kegg.txt", "r") as infile:
            data = infile.read()
        my_list_compound_sbml_not_found_from_reactant = data.splitlines()
        
        file_comp_not_found_from_product_sbml = open(directory+"//idcomp_notfound_2_kegg.txt", "w", encoding="ISO-8859-1")
        file_reaction_not_found_from_product_sbml = open(directory+"//idreaction_notfound_2_kegg.txt", "w", encoding="ISO-8859-1")
        
        #print("PARTE 2 METODO ALTERNATIVO")
        for compound in my_list_compound_sbml_not_found_from_reactant:
            #print("parte2", compound)
            compound_no_stoich = re.sub("\d+\.\d+", "", compound) # retirada de todos os valores estequiometricos
            
            if "- reduced" in compound:
                compound_no_stoich = compound_no_stoich.replace("- reduced", " ")
            
            param_splt = "" # preparo do split para separar a composicao quimica
            if "<=>" in compound_no_stoich:
                param_splt = "<=>"
            else:
                param_splt = "-->"
            
            compound_splt = compound_no_stoich.split(param_splt) # separacao da composicao quimica pelo reactante(0) e o produto(1) 
            
            # Caso seja uma reacao de insercao ou uma reacao de excrecao, ele ignora
            if compound_splt[0].strip() == "" or compound_splt[1].strip() == "":
                continue  
            
            product_sbml = compound_splt[1].strip()
        
            list_id_cpd_kegg = []
            len_product_sbml = 0
            
            # verifica os reagentes que possuem mais de um composto envolvido
            if " + " in product_sbml: #(ex.: acetyl-coa + atp + bicarbonate)
                product_sbml_splt = product_sbml.split(" + ")
                len_product_sbml = len(product_sbml_splt)
                
                # iteracao dentro dos compostos dos reactants (acetyl-coa + atp + bicarbonate)
                for cpd_product_sbml in product_sbml_splt:
                    #http://rest.kegg.jp/find/compound/acetyl-coa
                    result_id_cpd = k.find("compound", cpd_product_sbml)
                    result_id_cpd_splt = result_id_cpd.split('\n')
                    result_id_cpd_splt = filter(None, result_id_cpd_splt)
                    result_id_cpd_splt_filter = list(result_id_cpd_splt)
                    
                    # iteracao dentro do resultado encontrado para um dos compostos do reagente 
                    # (cpd:C00024 Acetyl-CoA; Acetyl coenzyme A)
                    if len(result_id_cpd_splt_filter) > 0:
                        for result_cpd in result_id_cpd_splt_filter:
                            local_result_cpd_splt = result_cpd.split('\t')
                            if ";" in local_result_cpd_splt[1]:
                                dsc_cpd_splt = local_result_cpd_splt[1].split(';')
                                for cpd_splt in dsc_cpd_splt:
                                    if cpd_splt.strip().lower() == cpd_product_sbml.strip().lower():
                                        list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                        break
                            else:
                                if cpd_product_sbml.lower() == local_result_cpd_splt[1].lower():
                                    list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                    break
            else:
                len_product_sbml = 1
                result_id_cpd = k.find("compound", product_sbml)
                result_id_cpd_splt = result_id_cpd.split('\n')
                result_id_cpd_splt = filter(None, result_id_cpd_splt)
                result_id_cpd_splt_filter = list(result_id_cpd_splt)
                
                if len(result_id_cpd_splt_filter) > 0:
                    for result_cpd in result_id_cpd_splt_filter:
                        local_result_cpd_splt = result_cpd.split('\t')
                        if ";" in local_result_cpd_splt[1]:
                            dsc_cpd_splt = local_result_cpd_splt[1].split(';')
                            for cpd_splt in dsc_cpd_splt:
                                if cpd_splt.strip().lower() == product_sbml.strip().lower():
                                    list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                    break
                        else:
                            if product_sbml.strip().lower() == local_result_cpd_splt[1].strip().lower():
                                list_id_cpd_kegg.append(local_result_cpd_splt[0])
                                break
        
            # caso nao tenha encontrado todos os compostos, ele nao segue o fluxo
            if len(list_id_cpd_kegg) != len_product_sbml:
                file_comp_not_found_from_product_sbml.write("{0}\n".format(compound))
                continue
            
            # PREPARO DA APLICACAO PARA USAR OS IDS DE COMPOSTOS E BUSCAR PELAS REACOES DE CADA UM
            product_sbml_in_cpd = "+".join(list_id_cpd_kegg)
            result_link_reactions_cpd = k.link("reaction", str(product_sbml_in_cpd))
            result_link_reactions_cpd_splt = result_link_reactions_cpd.split('\n')
            result_link_reactions_cpd_splt = filter(None, result_link_reactions_cpd_splt)
            result_link_reactions_cpd_splt_filter = list(result_link_reactions_cpd_splt)
            df_link_reaction_cpd = pd.DataFrame(columns=['id_cpd', 'id_reaction'])
            
            index_link_react_cpd = 0
            if len(result_link_reactions_cpd_splt_filter) > 0:
                for item_result_link_reactions_cpd in result_link_reactions_cpd_splt_filter:
                    local_item_result_link = item_result_link_reactions_cpd.split('\t')
                    df_link_reaction_cpd.loc[index_link_react_cpd] = [local_item_result_link[0], local_item_result_link[1]]
                    index_link_react_cpd += 1
            else:
                file_reaction_not_found_from_product_sbml.write("{0}\n".format(compound))
                continue
            
            df_link_reaction_cpd = df_link_reaction_cpd.groupby("id_reaction").filter(lambda x: len(x) == len(product_sbml_splt))
            set_id_reaction_kegg = {x.id_reaction for x in df_link_reaction_cpd.itertuples()}
            
            if len(set_id_reaction_kegg) > 0:
                for item_id_prod in set_id_reaction_kegg:
                    result_ec_number = k.link("enzyme", item_id_prod)
                    result_ec_number_splt = result_ec_number.split('\n')
                    result_ec_number_splt = filter(None, result_ec_number_splt)
                    result_ec_number_filter = list(result_ec_number_splt)
                    
                    for result_ec in result_ec_number_filter:
                        local_result_ec = result_ec.split('\t')
                        local_result_ec[1] = local_result_ec[1].replace("ec:", "").strip()
                        file_ecs.write("{0}\n".format(local_result_ec[1]))
                        file_ecs_compound.write("{0};{1}\n".format(compound, local_result_ec[1]))
            
        file_ecs.close()
        file_ecs_compound.close()
        file_comp_not_found_from_product_sbml.close()
        file_reaction_not_found_from_product_sbml.close()

    # METHOD TO TEST EMAIL
    def sendMailTest(self):

        email = 'gianninimail@gmail.com'
        method = '0'
        nameFile = ""
        if method == "1":
            methodSelected = "FBA+FVA"
        else:
            methodSelected = "Only FBA"

        if nameFile == "":
            nameFile = "NAME NOT FOUND"

        mensagem = "FINAL REPORT - FIND TARGETS WEB \n\n\n" \
                   "" \
                   "FILE: " + str(nameFile) + " | METHOD: " + str(methodSelected) + "\n" \
                                                                                    "" \
                                                                                    "YOUR ANALYSIS HAS FINISHED SUCESSFULLY. THE FOLLOWING FILES HAVE BEEN GENERATED:\n" \
                                                                                    "" \
                                                                                    "08-filter_ECNumbers_drugbank.xls\n" \
                                                                                    "FIELDS: EC NUMBER, PRODUCT, ORGANISM NAME, UNIPROTID, DRUGBANKID\n\n" \
                                                                                    "" \
                                                                                    "11-hits_Uniprot.xls\n" \
                                                                                    "FIELDS: UNIPROTID, HIT_UNIPROTID, PERCENT_IDENT_BLAST, EVALUE, GENE_NAME, PATHWAY, FUNCTION, CATALYTIC ACTIVITY, LOCALIZATION, ID PDB\n\n" \
                                                                                    "" \
                                                                                    "13-list_inhibitors_per_target.xls AND 14-list_inhibitors_approved.xls\n" \
                                                                                    "FIELDS: UNIPROTID	DRUGBANKID	DRUGBANKDRUGID	DRUGNAME	DRUGGROUP	PHARMAACTION	ACTIONS\n\n" \
                                                                                    "" \
                                                                                    "model_data.xls\n" \
                                                                                    "ALL GENES, REACTIONS AND METABOLITES IN THE SBML FILE USED IN ANALYSIS\n\n" \
                                                                                    "" \
                                                                                    "summary_results.xls\n" \
                                                                                    "FILE WITH A SUMMARY OF DATA FOUND IN 08, 11 AND 13 XLS FILES.\n\n" \
 \
            # remetente    = 'findtargetsweb_fiocruz@hotmail.com'
        remetente = 'findtargetweb@gmail.com'
        # senha        = 'proccfiocruz1234'
        # senha = 'whzctyqjvqsxojcz'
        # senha = 'pjcgoufuudsazwgq'
        senha = 'ntucaeevfacubqyj'

        # Informacoes da mensagem
        destinatario = email
        assunto = 'REPORT THERAPEUTICS TARGETS FROM YOUR NETWORK MODEL'

        msg = MIMEMultipart()

        msg['From'] = remetente
        msg['To'] = destinatario
        msg['Subject'] = assunto

        # Preparando a mensagem
        '''
        mensagem = '\r\n'.join([
          'From: %s' % remetente,
          'To: %s' % destinatario,
          'Subject: %s' % assunto,
          '',
          '%s' % mensagem
          ])
        '''
        mensagem = '\r\n'.join([
            '%s' % mensagem
        ])

        mensagem = mensagem.encode("UTF-8")

        msg.attach(MIMEText(mensagem.decode("UTF-8"), 'plain'))

        # filename = zip_file_report.filename
        # attachment = open(filename, "rb")
        #
        # part = MIMEBase('application', 'zip')
        # part.set_payload((attachment).read())
        # encoders.encode_base64(part)
        # part.add_header('Content-Disposition', "attachment", filename=os.path.basename(filename))
        # msg.attach(part)

        # Enviando o email (USANDO O SMTP DO HOTMAIL PARA ENVIAR)
        # server = smtplib.SMTP("smtp.live.com: 587")
        server = smtplib.SMTP("smtp.gmail.com: 587")
        server.starttls()
        server.login(remetente, senha)
        text = msg.as_string()
        server.sendmail(remetente, destinatario, text)
        server.quit()