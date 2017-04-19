# kd-tree-overlapper

Installation:

1. Get FLANN:
git clone https://github.com/mariusmuja/flann.git

2. Compile FLANN libs (If necessary install "cmake"):
cd flann;
cmake .;
make;
cd ..

3. Get SeqAn:
git clone https://github.com/seqan/seqan.git

4. Get kd-tree-overlapper:
git clone https://github.com/dzif/kd-tree-overlapper.git

5. Compile kd-tree-overlapper:
cd kd-tree-overlapper;
make

Usage: 
./kd [-o output=overlaps.out -k k-mer_len=4 -r kd-tree_iterations=600 
-l tags_len=1200 -s tags_spacing=600 -n num_NN=40 -m min_tags_space=200 
-w GC_window=100] [-i] input_file  

-o outputfile with overlaps  <br  />
-k length of k-mers for geometric embedding, default=4  
-r number of iterations for ANN search, default=600  
-l length of tags (subreads placed on reads), default=1200  
-s average distance between tags, default=600  
-n number of ANNs per tag, default=40  
-m minimum allowed distange between tags, default=200  
-w window size for GC-peaks detection, default=100  
-i input FASTA file  
