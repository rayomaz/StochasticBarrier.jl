FROM julia:1.11.5

# --- StochasticBarrier.jl setup ---
COPY ./StochasticBarrierFunctions /StochasticBarrierFunctions
RUN chmod +x /StochasticBarrierFunctions/run_sos.bash
WORKDIR /StochasticBarrierFunctions

# Precompile Julia package
RUN julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); Pkg.precompile()'

# Alias to run experiment
RUN echo 'alias stochasticbarrier="/StochasticBarrierFunctions/run_sos.bash"' >> ~/.bashrc

# -------------------------------
# Entrypoint
# -------------------------------
# Allows interactive bash with environment ready
ENTRYPOINT ["/bin/bash"]