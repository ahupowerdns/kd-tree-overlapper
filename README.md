# kd-tree-overlapper

Fast and memory-efficient noisy read overlapping with KD-trees.

Software licensed under the [Apache License 2.0](LICENSE.txt), using the [FLANN](https://github.com/mariusmuja/flann) and [SeqAN](https://github.com/seqan/seqan) software libraries (and including [klib](https://github.com/attractivechaos/klib)'s [kseq.h](kseq.h) to parse FASTA/FASTQ input). Please see [INSTALL.md](INSTALL.md) for instructions how to install `kd`.

### Usage:

```
./kd [options] -i <input_file>

Options:
     -o outputfile with overlaps [overlaps.out]
     -k length of k-mers for geometric embedding [4]
     -r number of iterations for ANN search [600]
     -l length of tags (subreads placed on reads) [1200]
     -s average distance between tags [600]
     -n number of ANNs per tag, [40]
     -m minimum allowed distange between tags [200]
     -w window size for GC-peaks detection [100]
     -i input FASTA file
```
