FROM rocker/r-ver:4.5.1

LABEL org.opencontainers.image.source="https://github.com/joookerz/sc-pcQTL"
LABEL org.opencontainers.image.description="Core R environment for sc-pcQTL"

RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential git libcurl4-openssl-dev libssl-dev libxml2-dev \
      libopenblas-dev plink1.9 ca-certificates \
    && ln -s /usr/bin/plink1.9 /usr/local/bin/plink \
    && rm -rf /var/lib/apt/lists/*

RUN R -q -e "install.packages(c('data.table','R.utils','Rcpp','RcppArmadillo','RcppEigen','BH','fastglm','roptim','testthat','remotes'), repos='https://cloud.r-project.org')" \
    && R -q -e "remotes::install_github('mkanai/fasthurdle@72c2cd7a91208d6285c52a80f198f045ba51b5d7', upgrade='never')" \
    && R -q -e "stopifnot(as.character(packageVersion('fasthurdle')) == '1.2.0')"
