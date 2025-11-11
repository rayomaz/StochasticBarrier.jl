FROM julia:1.11.5

# --- StochasticBarrier.jl setup ---
COPY ./StochasticBarrierFunctions /StochasticBarrierFunctions
RUN chmod +x /StochasticBarrierFunctions/run_sos.bash
WORKDIR /StochasticBarrierFunctions

# Precompile Julia package
RUN julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); Pkg.precompile()'

# Aliases to run experiments
RUN echo 'stochasticbarrier() {' >> ~/.bashrc \
    && echo '  case "$1" in' >> ~/.bashrc \
    && echo '    sos) shift; /StochasticBarrierFunctions/run_sos.bash "$@" ;;' >> ~/.bashrc \
    && echo '    pwc) shift; /StochasticBarrierFunctions/run_pwc.bash "$@" ;;' >> ~/.bashrc \
    && echo '    *) echo "Usage: stochasticbarrier {sos|pwc}" ;;' >> ~/.bashrc \
    && echo '  esac' >> ~/.bashrc \
    && echo '}' >> ~/.bashrc

# -------------------------------
# Entrypoint
# -------------------------------
# Allows interactive bash with environment ready
ENTRYPOINT ["/bin/bash"]