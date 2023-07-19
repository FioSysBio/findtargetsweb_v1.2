# FindTargetsWEB v1.2

## Passos para deploy da aplicação pela Primeira vez:

tutorial de referência no no link: https://www.alura.com.br/artigos/fazendo-o-deploy-de-uma-aplicacao-django?gclid=CjwKCAjwp6CkBhB_EiwAlQVyxWHGgdsq70qLO77L_yPbeQ2_3U370eQS7FwOwvDT9ZcLtF8vt7laXBoCrBoQAvD_BwE

* instalar a versão de python 3.10
* instalar as libs do projeto (env.txt)
* instalar o gunicorn (pip install)
* fazer a clonagem do código no repositório no Github (ou git pull caso tenha já clonado o prjeto)
* fazer a inicialização do serviço com o comando na pasta do projeto /sistemas/findtargetsweb_v1.2:
  - source findtargetsweb/bin/activate
  - gunicorn -b 0.0.0.0:8184 findtargetsweb.wsgi


## Passos para deploy da aplicação em fase de atualização:

* fazer a parada do serviço buscando no gerenciador de tarefas "htop" a palavra "gunicorn" e fechar o serviço
* fazer a execução do comando git pull
* fazer a inicialização do serviço com o comando na pasta do projeto /sistemas/findtargetsweb_v1.2:
  - source findtargetsweb/bin/activate
  - gunicorn -b 0.0.0.0:8184 findtargetsweb.wsgi