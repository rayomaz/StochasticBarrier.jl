FROM julia:1.11.5

# Install tools required by Mosek build
RUN apt-get update && \
    apt-get install -y tar bzip2 curl ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# --- StochasticBarrier.jl setup ---
COPY ./StochasticBarrierFunctions /StochasticBarrierFunctions
# Drop any host-side mosek.lic: the license must be mounted at runtime.
RUN rm -f /StochasticBarrierFunctions/benchmarks/mosek/mosek.lic
RUN chmod +x /StochasticBarrierFunctions/run_sos.bash \
    && chmod +x /StochasticBarrierFunctions/run_pwc.bash
WORKDIR /StochasticBarrierFunctions

# Mosek license is expected to be bind-mounted at /mosek/mosek.lic.
# Benchmark CSVs are written under $SB_RESULTS_DIR (bind-mount /results).
ENV MOSEKLM_LICENSE_FILE=/mosek/mosek.lic
ENV SB_RESULTS_DIR=/results

# Precompile Julia package
ENV JULIA_PROJECT='/StochasticBarrierFunctions/benchmarks'
RUN julia --project="$JULIA_PROJECT" -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

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
