FROM zakandrewking/me_dependencies:latest

RUN pip install git+git://github.com/zakandrewking/parallel_pandas.git

# install theseus and copy in models
RUN git clone https://github.com/zakandrewking/theseus.git
RUN cd theseus && pip install -e .
COPY data/prototype_*.pickle ./theseus/theseus/data/models/
COPY data/iML1515.xml ./theseus/theseus/data/models/

# install minime
COPY minime minime
RUN cd minime && pip install -e .

# install soplex
COPY lib lib
RUN pip install lib/soplex-0.0.5b2-cp27-none-linux_x86_64.whl

# install hindsight
COPY hindsight hindsight/hindsight
COPY setup.py hindsight/
RUN cd hindsight && pip install -e .

COPY bin bin
CMD python bin/3_run_secretion_tree.py
