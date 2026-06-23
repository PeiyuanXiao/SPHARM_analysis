# 1. Use the pre-built geospatial image

#    Pin the platform to linux/amd64. The conda env (analysis/scripts/
#    environment.yml) locks linux-64 package build hashes for numerical
#    reproducibility, so the image MUST be amd64. Without this flag, building
#    on an Apple Silicon Mac would target linux/arm64 and the conda solve
#    would fail (those build hashes do not exist for arm64). On amd64 hosts
#    (Linux, Windows/Intel, CI) this flag is the native default and is a no-op;
#    on Apple Silicon it forces emulation, giving bit-identical results.
FROM --platform=linux/amd64 rocker/geospatial:4.4.2

# 2. Install Python and Conda dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    libgl1 \
    libglu1-mesa \
    librsvg2-dev \
    libglpk-dev \
    libgit2-dev \
    curl \
    xvfb \
    && rm -rf /var/lib/apt/lists/*

# 3. Install Miniconda
#    conda 26.3.2 is required for `conda tos accept`.
#    conda goes at the end of the path, not the beginning
#.   coda for both amd64 and arm64 architectures
ARG TARGETARCH
RUN case "$TARGETARCH" in \
      "amd64")  url="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" ;; \
      "arm64")  url="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh" ;; \
    esac && \
    wget "$url" -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh
ENV PATH=$PATH:/opt/conda/bin

# --- RSTUDIO PROJECT AUTO-LOAD CONFIG ---
RUN mkdir -p /home/rstudio/.local/share/rstudio/projects_settings
RUN echo "/project/SPHARM_analysis.Rproj" > /home/rstudio/.local/share/rstudio/projects_settings/last-project-path
RUN mkdir -p /home/rstudio/.config/rstudio
RUN echo '{"initial_working_directory": "/project"}' > /home/rstudio/.config/rstudio/rstudio-prefs.json
RUN chown -R rstudio:rstudio /home/rstudio/.local /home/rstudio/.config

# --- TERMINAL CONFIG ---
RUN echo 'cd /project' >> /home/rstudio/.bashrc
RUN echo "source /opt/conda/etc/profile.d/conda.sh" >> /home/rstudio/.bashrc

# 4. Set up Project
WORKDIR /project
COPY . /project

# --- 5. RENV RESTORE ---
# A. Set RENV paths to location OUTSIDE /project
ENV RENV_PATHS_LIBRARY=/opt/renv/library
ENV RENV_PATHS_CACHE=/opt/renv/cache

# B. Create directories and give 'rstudio' user permission
RUN mkdir -p /opt/renv && chown -R rstudio:rstudio /opt/renv

# C. Remove ALL pre-installed system R packages to avoid renv conflicts
#    rocker/geospatial ships ~280 packages in the system library that
#    collide with renv's project library during restore.
RUN rm -rf /usr/local/lib/R/site-library/*

# D. Restore the renv project library.
#    Do NOT `install.packages('renv')` here: that pulls the LATEST renv from
#    CRAN (e.g. 1.2.3), newer than the version the lockfile pins (1.1.8). With
#    this lockfile a newer renv aborts the restore -- it installs 154 packages
#    then reports ~79 as "dependency failed" -- because Posit PPM (the lock's
#    moving "latest" repo) has rebuilt a pinned binary under a suffixed version
#    (lock pins Rcpp 1.1.1, PPM now serves 1.1.1-1; renv's "dependency tree was
#    repaired" message). The renv/activate.R sourced at R startup already
#    bootstraps the lockfile-pinned renv (1.1.8), which restores all 232
#    packages cleanly.
RUN R -e "options(renv.config.cache.symlinks = FALSE); renv::restore(prompt = FALSE)"

# --- 6. Build Python Env ---
#    environment.yml pins every package including BLAS/LAPACK build hashes.
#    A single conda solve is sufficient — do NOT add a second `conda install`
#    step after this, as it would trigger a re-solve and may alter numerical
#    library builds, breaking floating-point reproducibility.
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    conda env create -f analysis/scripts/environment.yml --solver=libmamba && \
    conda clean -afy

RUN git config --global --add safe.directory /project

# --- 7. Pre-create targets store directory with correct permissions ---
#    Prevents occasional tar_make() warning:
#    "cannot create file 'analysis/paper/_targets/meta/meta'"
RUN mkdir -p /project/analysis/paper/_targets/meta && \
    chmod -R 777 /project/analysis/paper/_targets