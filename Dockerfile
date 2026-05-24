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

# C. Remove ALL pre-installed system R packages to avoid renv conflicts
#    rocker/geospatial ships ~280 packages in the system library that
#    collide with renv's project library during restore.
RUN rm -rf /usr/local/lib/R/site-library/*

# D. Install renv and restore
RUN R -e "install.packages('renv', repos='https://cloud.r-project.org')" && \
    R -e "options(renv.config.cache.symlinks = FALSE); renv::restore(prompt = FALSE)"

# --- 6. Build Python Env ---
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    conda env create -f analysis/scripts/environment.yml --solver=libmamba

RUN git config --global --add safe.directory /project

# --- 7. Pin binary Python geospatial/image stack inside the conda env ---
#    These packages must come from the same conda-forge solve to avoid
#    libtiff/libjpeg ABI mismatches when SPHARM imports PIL/VTK.
RUN conda install -n spharm -c conda-forge --override-channels \
    vtk libtiff libjpeg-turbo pillow pymeshfix networkx openpyxl \
    && conda clean -afy
    