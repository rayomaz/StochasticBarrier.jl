FROM julia:1.11.5

# --- StochasticBarrier.jl setup ---
COPY ./StochasticBarrierFunctions /StochasticBarrierFunctions
RUN chmod +x /StochasticBarrierFunctions/run_sos.bash
WORKDIR /StochasticBarrierFunctions

# Precompile Julia package
RUN julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); Pkg.precompile()'

# Aliases to run experiments
RUN echo 'alias stochasticbarrier="/StochasticBarrierFunctions/run_sos.bash"' >> ~/.bashrc \
    && echo 'alias stochasticpwc="/StochasticBarrierFunctions/run_pwc.bash"' >> ~/.bashrc

# -------------------------------
# Entrypoint
# -------------------------------
# Allows interactive bash with environment ready
ENTRYPOINT ["/bin/bash"]