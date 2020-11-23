FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# TODO ......
WORKDIR /kb/module/bin

RUN apt-get update && \
apt-get install -y wget zip gcc vim tree

RUN cd /opt && \
wget https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX/SECTION_Download/MODULE_Downloads/CLASS_Old%20versions/UNIT_FAPROTAX_1.2.1/FAPROTAX_1.2.1.zip && \ 
unzip FAPROTAX_1.2.1.zip && \ 
rm FAPROTAX_1.2.1.zip && \
cd FAPROTAX_1.2.1 && \
sed -i '1s/\/usr\/bin\/python/\/miniconda\/envs\/py2\/bin\/python/' collapse_table.py && \
mv FAPROTAX.txt /kb/module/data

ENV PATH=/opt/FAPROTAX_1.2.1:$PATH

RUN conda create --yes --name py2 python=2.7 numpy

RUN pip install numpy==1.15.4 pandas
RUN pip install dotmap

RUN pip install numpy==1.19.2 # fewer warnings

ENV PYTHONUNBUFFERED=True

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
