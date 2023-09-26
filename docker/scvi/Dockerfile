FROM nciccbr/ccbr_ubuntu_base_20.04:v5

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}

# install dependencies
RUN python3 -m pip install --upgrade pip
RUN pip3 install click numpy scipy scikit-learn
RUN pip install scvi-tools
RUN pip install "scvi-tools[tutorials]"
RUN pip install .
RUN conda install -c conda-forge r-base r-essentials r-reticulate

# cleanup etc
# Save Dockerfile in the docker
COPY Dockerfile /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
RUN chmod a+r /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
#COPY argparse.bash /opt2
RUN chmod -R a+rx /opt2/argparse.bash
ENV PATH="/opt2/:$PATH"
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/
WORKDIR /data2