FROM julia:1.11.5

# --- StochasticBarrier.jl setup ---
COPY ./StochasticBarrierFunctions /StochasticBarrierFunctions
RUN chmod +x /StochasticBarrierFunctions/run_sos.bash \
    && chmod +x /StochasticBarrierFunctions/run_pwc.bash
WORKDIR /StochasticBarrierFunctions

# Precompile Julia package
RUN julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); Pkg.precompile()'

# Single alias that dispatches based on argument
RUN echo 'stochasticbarrier() {' >> ~/.bashrc \
    && echo '  case "$1" in' >> ~/.bashrc \
    && echo '    sos) /StochasticBarrierFunctions/run_sos.bash "$@" ;;' >> ~/.bashrc \
    && echo '    pwc) /StochasticBarrierFunctions/run_pwc.bash "$@" ;;' >> ~/.bashrc \
    && echo '    *) echo "Usage: stochasticbarrier {sos|pwc}" ;;' >> ~/.bashrc \
    && echo '  esac' >> ~/.bashrc \
    && echo '}' >> ~/.bashrc

# -------------------------------
# Entrypoint
# -------------------------------
ENTRYPOINT ["/bin/bash"]
