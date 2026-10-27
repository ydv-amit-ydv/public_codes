# MHC-I BA / EL: NetMHCpan 4.2 runner

A local Flask page around DTU's NetMHCpan 4.2 (`4.2d` static Linux build).
Paste a FASTA, protein or peptide list, choose alleles and lengths, and get a
sortable table of EL/BA scores, ranks, affinity (nM) and strong/weak binder calls.

## Install NetMHCpan (WSL / Ubuntu)

NetMHCpan is DTU academic-licence software: download it yourself and do not
commit it here.

```bash
sudo apt install -y tcsh gawk
mkdir -p ~/local_engines/netmhc && cd ~/local_engines/netmhc
tar -xzf netMHCpan-4.2*.Linux.tar.gz            # gives netMHCpan-4.2/ (with data/ inside)
# optional: the engine runs a patched copy, but for command-line use set NMHOME in netMHCpan-4.2/netMHCpan
echo -e "SIINFEKL\nGILGFVFTL\nNLVPMVATV" > ~/pep.txt
netMHCpan-4.2/netMHCpan -p ~/pep.txt -a HLA-A02:01 -BA
```

Expected: GILGFVFTL and NLVPMVATV strong binders (SB), SIINFEKL (a mouse Kb
epitope) a non-binder. `cannot execute: required file not found` means tcsh is
missing; `no binaries found ... /netMHCpan-4.2/...` means `NMHOME` was not edited.

## Run the GUI

Use the combined server (see the engine `README.md`): put this code in a folder beside
`netmhc/` and run `bash start.sh`. It finds NetMHCpan automatically, wherever the
tarball was extracted, and runs a patched copy of the script if `NMHOME` was never
edited. This tool alone: `python3 netmhcpan_app.py` (port 5000, `PORT` to change).

## Behaviour worth knowing

- **Two input modes.** *Scan* slides k-mers of the chosen lengths (8-14) over each
  FASTA record; *peptides as given* takes one peptide per line. Either way
  NetMHCpan runs in `-p` mode, one run per length, because `-p` takes a single
  length per file. Windows with non-standard residues are skipped and reported.
- **Allele presets** (filtered against `netMHCpan -listMHC`; hover a button to see
  any allele your install lacks): IEDB human reference 27, nine HLA supertype
  representatives, mouse H-2 (all, C57BL/6, BALB/c), plus a searchable list of
  every allele NetMHCpan supports.
- **Parameters exposed:** `-BA`, `-rankS`, `-rankW`, `-t`, and a raw-args box for
  anything else. Binders are SB at %Rank_EL < 0.5 and WB at < 2 by default.
- **Output columns** are prefixed `NMP_` (`NMP_Score_EL`, `NMP_%Rank_EL`,
  `NMP_Score_BA`, `NMP_%Rank_BA`, `NMP_Aff_nM`, `NMP_BindLevel`) alongside
  `Protein`, `Start`, `Peptide`, `Allele`, `Core`, so later models can be merged
  on peptide and allele without clashes.
- The server binds to 127.0.0.1 only.
- The parser reads each result table's own header line and has been checked
  against real NetMHCpan 4.2d output for the three-peptide test above.

## Not yet done

- Export to the `MHC_1_Binding_Elution` CSV schema the Workbench ingests.
- Context encoding (`-context`), XLS output and path overrides are not in the UI.
