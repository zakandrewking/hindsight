FROM zakandrewking/me_dependencies:latest

RUN pip install ipy_progressbar

COPY lib/soplex-0.0.5b2-cp27-none-linux_x86_64.whl .
RUN pip install soplex-0.0.5b2-cp27-none-linux_x86_64.whl

COPY minime minime
RUN cd minime && pip install -e .

COPY parallel_pandas parallel_pandas
RUN cd parallel_pandas && pip install -e .

COPY theseus theseus
RUN cd theseus && pip install -e .

COPY data/prototype_*.pickle ./theseus/theseus/data/models/
COPY data/iML1515.xml ./theseus/theseus/data/models/

COPY bin/3_run_secretion_tree.py ./
COPY data/sims_table.pickle ./

CMD python 3_run_secretion_tree.py
