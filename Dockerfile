FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

WORKDIR /kb/module/bin

RUN apt-get update && apt-get install -y wget zip gcc

RUN conda create --name py2 python=2.7

RUN /miniconda/envs/py2/bin/pip install numpy h5py && /miniconda/envs/py2/bin/pip install biom-format

RUN wget https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX/SECTION_Download/MODULE_Downloads/CLASS_Latest%20release/UNIT_FAPROTAX_1.2.1/FAPROTAX_1.2.1.zip && unzip FAPROTAX_1.2.1.zip && rm -rf FAPROTAX_1.2.1.zip

RUN cd FAPROTAX_1.2.1 && chmod +x collapse_table.py && sed -i 's/.usr.bin.python/\/miniconda\/envs\/py2\/bin\/python/' collapse_table.py

ENV PATH=/kb/module/bin/FAPROTAX_1.2.1:$PATH

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
