FROM julia:1.11.5


# --- PART 1: StochasticBarrier.jl setup ---
COPY ./StochasticBarrierFunctions /StochasticBarrierFunctions
RUN chmod +x /StochasticBarrierFunctions/run_sos.bash
WORKDIR /StochasticBarrierFunctions

# Precompile Julia package
RUN julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); Pkg.precompile()'

# Alias to run experiment
RUN echo 'alias stochasticbarrier="/StochasticBarrierFunctions/run_sos.bash"' >> ~/.bashrc

# --- PART 2: PRoTECT setup ---
# Base Image
# -------------------------------
FROM ubuntu:22.04

# Set working directory
WORKDIR /app

# -------------------------------
# Install Python 3.10 and system dependencies
# -------------------------------
RUN apt-get update && apt-get install -y \
    software-properties-common \
    curl \
    git \
    unzip \
    nano \
    && add-apt-repository ppa:deadsnakes/ppa -y \
    && apt-get update \
    && apt-get install -y \
        python3.10 \
        python3.10-venv \
        python3.10-distutils \
        python3-pip \
    && rm -rf /var/lib/apt/lists/*

# -------------------------------
# Set up a Python virtual environment
# -------------------------------
RUN python3.10 -m venv /opt/protect-venv
ENV PATH="/opt/protect-venv/bin:$PATH"

# Upgrade pip inside venv
RUN pip install --upgrade pip

# -------------------------------
# Clone ProTECT repository
# -------------------------------
RUN git clone https://github.com/Kiguli/ProTECT.git /app/PRoTECT

# -------------------------------
# Install ProTECT Python dependencies
# -------------------------------
RUN pip install --no-cache-dir -r /app/PRoTECT/requirements.txt

# -------------------------------
# Set PYTHONPATH to include ProTECT
# -------------------------------
ENV PYTHONPATH="/app/PRoTECT:$PYTHONPATH"

# -------------------------------
# Entrypoint
# -------------------------------
# Allows interactive bash with environment ready
ENTRYPOINT ["/bin/bash"]