# Installation

For convenience, we provide a statically linked `kd` Linux binary here: [https://github.com/dzif/kd-tree-overlapper/releases](https://github.com/dzif/kd-tree-overlapper/releases).

## Install from source

Building `kd` from source requires `make` and `cmake`. We have tested `kd` with FLANN v1.9.1 and SeqAn v2.3.2. 

#### Install FLANN:

```
git clone https://github.com/mariusmuja/flann.git
compile libs (If necessary install "cmake"):
cd flann;
cmake .;
make;
cd ..
```

#### Install SeqAn:

```
git clone https://github.com/seqan/seqan.git
```

#### Install kd-tree-overlapper:

```
git clone https://github.com/dzif/kd-tree-overlapper.git
cd kd-tree-overlapper;
make
```

This will place an executable named `kd` in the folder `kd-tree-overlapper`.

## KD-tree biobox

Following the assembler biobox specification, we provide an overlapper biobox to run `kd` in a Docker container.

#### Build the Biobox

```
docker build -t kdtree .
```

#### Run the biobox

```
docker run \
  --volume="$(pwd)/biobox.yaml:/bbx/input/biobox.yaml:ro"   \
  --volume="$(pwd)/input.fasta:/bbx/input/input.fasta:ro"   \
  --volume="$(pwd)/output:/bbx/output:rw" \
  --rm \
  kdtree \
  default
```

where

  * **input.fasta** is the fasta you would like to use

  * **output** is the directory where the overlap information will be placed

  * **biobox.yaml** is the yaml that describes your input, e.g.

```
---
version: "0.9.0"
arguments:
  - fasta:
    - id: "test_fasta"
      value: "/bbx/input/input.fasta"
```
