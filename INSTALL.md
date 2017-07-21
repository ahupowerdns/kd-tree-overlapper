# Installation

For convenience, we provide a statically linked `kd` binary here: [https://github.com/dzif/kd-tree-overlapper/releases](https://github.com/dzif/kd-tree-overlapper/releases).

Building `kd` from source requires `make` and `cmake` on your system. We have tested `kd` with FLANN v1.9.1 and SeqAn v2.3.2. 

### Install FLANN:

```
git clone https://github.com/mariusmuja/flann.git
libs (If necessary install "cmake"):
cd flann;
cmake .;
make;
cd ..
```

### Install SeqAn:

```
git clone https://github.com/seqan/seqan.git
```

### Install kd-tree-overlapper:

```
git clone https://github.com/dzif/kd-tree-overlapper.git
cd kd-tree-overlapper;
make
```

This will place an executable named `kd` in the folder `kd-tree-overlapper`.
