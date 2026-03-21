# 1. Use the pre-built geospatial image
FROM rocker/geospatial:4.4.2

# 2. Install Python and Conda dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    libgl1 \
    libglu1-mesa \
    xvfb \
    && rm -rf /var/lib/apt/lists/*

# 3. Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh \
    && bash miniconda.sh -b -p /opt/conda \
    && rm miniconda.sh
ENV PATH=/opt/conda/bin:$PATH

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

# C. Remove the pre-installed system ggplot2 (v3.5.1)
# This forces R to use the v4.0.1 version that renv installs below.
RUN rm -rf /usr/local/lib/R/site-library/ggplot2

# D. Install renv and restore
RUN R -e "install.packages('renv', repos='https://cloud.r-project.org')" && \
    R -e "options(renv.config.cache.symlinks = FALSE); renv::restore(prompt = FALSE)"

# --- 6. Build Python Env ---
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    conda env create -f analysis/scripts/environment.yml --solver=libmamba
  
RUN git config --global --add safe.directory /project
RUN /opt/conda/envs/spharm/bin/pip install pymeshfix networkx openpyxl
