# from pipelineFindTargets import MyThread, FindTargets
# import front.pipelineFindTargets
# from pipelineFindTargets import FindTargets
# import pipelineFindTargets as find
import numpy as np
import cobra as cb
import pandas as pd
from pipelineFindTargets import MyThread, FindTargets


if __name__ == '__main__':
    #modelo = FindTargets()
    #model = modelo.readModel("iMO1056_MOPS_glycerol_integrated.xml")
    #roda = MyThread("Thiago", "Pseudomonas aeruginosa", "merigueti@gmail.com", model, "1")
    #roda.run()

    ## ERRO NO MOMENTO DE GERAR O XLS DO MODELO.
    #modelo = FindTargets()
    #model = modelo.readModel("iPAO1_MOPS_acetate_integrated.xml")
    #roda = MyThread("Thiago", "Pseudomonas aeruginosa", "merigueti@gmail.com", model, "1")
    #roda.run()

    # modelo = front.pipelineFindTargets.FindTargets()
    # model = modelo.readModel("iPAO1_MOPS_glycerol_integrated.xml")
    # roda = MyThread("Thiago", "Pseudomonas aeruginosa", "merigueti@gmail.com", model, "1")
    # roda.run()

    model = cb.io.load_model("textbook")
    print(len(model.reactions))
    print(len(model.metabolites))
    print(len(model.genes))
    cb.io.validate_sbml_model('/home/thiago/projetos/fiocruz/modelos/iMO1056_MOPS_glycerol_integrated.xml')
    model2 = cb.io.read_sbml_model('/home/thiago/projetos/fiocruz/modelos/iMO1056_MOPS_glycerol_integrated.xml')
    res = model.optimize()
    print(res.objective_value)
    print(res)

    res2 = model2.optimize()
    print(res2.objective_value)
    print(res2)
