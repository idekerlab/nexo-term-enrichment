# Dockerfile for NeXO Term Enrichment service.
# 
FROM ipython/scipystack

MAINTAINER Keiichiro Ono <kono@ucsd.edu>

RUN pip install flask flask-restful

# Form a set of standard directories.
RUN mkdir -p /downloads
RUN mkdir -p /work

# Install python libraries

# qvalue library
RUN cd /downloads && git clone https://github.com/nfusi/qvalue.git && cd qvalue && python setup.py install

# Get actual application
RUN cd /downloads && git clone https://github.com/idekerlab/nexo-term-enrichment.git && cd nexo-term-enrichment

# EXPOSE PORT
EXPOSE 5000

WORKDIR /downloads/nexo-term-enrichment
CMD python nexo.py