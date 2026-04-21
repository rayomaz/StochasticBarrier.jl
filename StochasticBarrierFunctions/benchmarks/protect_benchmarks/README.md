Read the description below for repeatability of the PRoTECT experiments.

### Docker Image
A pre-built image tarball is distributed alongside this repository. Load it with
```sh
sudo docker load -i protect-bench.tar
```

If you prefer to build it yourself the Dockerfile at `StochasticBarrierFunctions/benchmarks/Dockerfile.protect` is still available:
```sh
sudo docker build -f Dockerfile.protect -t protect-bench .
```

### Run the container
The Mosek licence is bind-mounted into the Docker container at run time, along with the results directory the container will write its CSV to.

```sh
mkdir -p results
sudo docker run --rm -it \
    -v "$(pwd)/mosek/mosek.lic:/mosek/mosek.lic:ro" \
    -v "$(pwd)/results:/results" \
    protect-bench
```

Inside the container `MOSEKLM_LICENSE_FILE=/mosek/mosek.lic` and
`SB_RESULTS_DIR=/results` are pre-set.

## Run through bash

Use the provided alias to run every experiment:

```sh
protectbench
```

Each experiment appends a row to `/results/protect_results.csv`. See README.md in the repo root for how to combine with the results from StochasticBarrier.jl and StochasticBarrierFunc. for a combined table.