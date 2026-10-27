# public_codes

## local_engine_codes

A local web GUI that runs MHC-I epitope predictors on your own machine:
**NetMHCpan 4.2** (binding affinity / eluted ligand) and **BigMHC** (presentation +
immunogenicity), with a Settings tab. It finds the models automatically.

The predictors themselves are not included: NetMHCpan is DTU academic-licence
software (download it from DTU), BigMHC is a separate open-source repository.

Quick start on Linux / WSL:

```bash
git clone https://github.com/ydv-amit-ydv/public_codes.git
mkdir -p ~/local_engines/runner_codes
cp -r public_codes/local_engine_codes/. ~/local_engines/runner_codes/
# put the NetMHCpan extraction in ~/local_engines/netmhc/ and the BigMHC clone in ~/local_engines/bigmhc/
cd ~/local_engines/runner_codes && bash start.sh      # http://localhost:5000
```

Full instructions: [`local_engine_codes/README.md`](local_engine_codes/README.md).
Update later with `git -C ~/public_codes pull` (or re-clone) and copy again.
