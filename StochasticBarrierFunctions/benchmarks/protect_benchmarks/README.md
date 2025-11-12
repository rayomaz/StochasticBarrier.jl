Read the description below for repeatability of the ProTECT experiments.

### Docker Image
The Dockerfile is provided in the **benchmarks** folder. Build this docker file to obtain all the required Python dependencies.
```sh
sudo docker build -f Dockerfile.protect -t protect-bench .
```

To start a container 

```sh
sudo docker run -it protect-bench
```

## Run through bash

Use the following alias command to run the benchmarks.

```sh
protectbench              
```
