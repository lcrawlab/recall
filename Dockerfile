FROM rocker/verse:4.0.5

RUN R -e "install.packages('.', type = 'source', repos = NULL)"