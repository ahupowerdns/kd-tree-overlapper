/**
 * Copyright 2016-2017 Dmitri Parkhomchuk, DZIF
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
**/

#include <flann/flann.hpp>
#include <stdio.h>
#include <map>
#include <tuple>
#include <string>

#include <seqan/modifier.h>

#include <sys/time.h>
#include <unistd.h>
#include "kseq.h"
#include <zlib.h>
#include <algorithm>

    #define REAL float
inline static REAL sqr(REAL x) {
        return x*x;
    }
int linreg(const std::vector<int> x, const std::vector<int> y, REAL* m, REAL* b, REAL* r)
    {
	int n = x.size();
        REAL   sumx = 0.0;                        /* sum of x                      */
        REAL   sumx2 = 0.0;                       /* sum of x**2                   */
        REAL   sumxy = 0.0;                       /* sum of x * y                  */
        REAL   sumy = 0.0;                        /* sum of y                      */
        REAL   sumy2 = 0.0;                       /* sum of y**2                   */

       for (int i=0;i<n;i++)   
          { 
          sumx  += x[i];       
          sumx2 += sqr(x[i]);  
          sumxy += x[i] * y[i];
          sumy  += y[i];      
          sumy2 += sqr(y[i]); 
          } 

       REAL denom = (n * sumx2 - sqr(sumx));
       if (denom == 0) {
           // singular matrix. can't solve the problem.
           *m = 0;
           *b = 0;
           *r = 0;
           return 1;
       }

       *m = (n * sumxy  -  sumx * sumy) / denom;
       *b = (sumy * sumx2  -  sumx * sumxy) / denom;
       if (r!=NULL) {
          *r = (sumxy - sumx * sumy / n) /          /* compute correlation coeff     */
                sqrt((sumx2 - sqr(sumx)/n) *
                (sumy2 - sqr(sumy)/n));
       }

       return 0; 
    }

timeval start_time_; //borrowed from FLANN examples
void start_timer(const std::string& message = "")
{
	if (!message.empty()) {
		printf("%s", message.c_str());
		fflush(stdout);
	}
    gettimeofday(&start_time_,NULL);
}

double stop_timer()
{
    timeval end_time;
    gettimeofday(&end_time,NULL);

	return double(end_time.tv_sec-start_time_.tv_sec)+ double(end_time.tv_usec-start_time_.tv_usec)/1000000;
}

using namespace seqan;

typedef unsigned short int flanns; // for tags with counts >255
//~ typedef unsigned char flanns; // for counts <256
typedef flann::Matrix<flanns> flannmat;

std::map<int, int> QIndex2AIndex; //mapping Q-gram index to frequency vector index

std::map<unsigned int, unsigned int> ReadId; //from flann index to fasta file index

int klen = 4; //q-gram length

int lf_len = 800; //length of tags

int tags_div = 400;//average space between tags (defines tags number)
int space = 200;//minimal space between tags

std::vector<flannmat> datasets;

std::string DNA = "ACGT";

std::map<std::string, int> Qmer2AIndex; //Qgram sequence to array index, revcoms mpped to the same
int Alen; //length of frequency vector

void populate_kmers_index()
{
	std::string nmer;
	int nt;
	int jj = 0;
	
	for (int i = 0; i < 1 << 2*klen; ++i)
	{
		nmer="";
		int ii=i;
		for (int j = klen - 1; j >= 0; --j)
		{
			nt=ii >> 2*j;
			nmer += DNA[nt];
			ii -= nt << 2*j;
		}
		if (Qmer2AIndex.count(nmer)) continue;
		else {
			Qmer2AIndex[nmer]=jj;
			seqan::reverseComplement(nmer);
			Qmer2AIndex[nmer]=jj;
			++jj;
			Alen = jj;
		}
		
	}
}

std::vector<unsigned int> len_loaded_seqs;

gzFile fp;
KSEQ_INIT(gzFile, gzread)
kseq_t *seq;

std::vector<unsigned int> istart;

std::vector<std::vector< unsigned int > * > coords;

using namespace std;
template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

  return idx;
}

int cg_win;//window for GC-content averaging
void tag(std::string &dna, std::vector<unsigned int> &rcoo) // tagging on GC-content peaks
{
	int wnd = cg_win;
	unsigned int ntags = dna.length()/tags_div;
	int at = 0, cg = 0;
	std::vector<int> prof; // (dna.length() - wnd, 2);
	for (unsigned int i = 0; i<dna.length() - wnd -1; ++i) {
		switch(dna.at(i)) {
			case 'A' : ++at;break;
			case 'T' : ++at;break;
			case 'C' : ++cg;break;
			case 'G' : ++cg;break;
		}
		switch(dna.at(i + wnd)) {
			case 'A' : --at;break;
			case 'T' : --at;break;
			case 'C' : --cg;break;
			case 'G' : --cg;break;
		}
		prof.push_back( abs(at-cg) );
		
					
	}
	
	int ofs = -wnd/2;// + lf_len/2;
	for (auto i: sort_indexes(prof)) {
          if (int(i) < ofs  or int(i) - ofs > int(dna.length()) - lf_len -1 ) continue;
		bool f = true;
		if (rcoo.size() ==0) rcoo.push_back(int(i) - ofs);
		else {
			if (rcoo.size() >= ntags) break;
			for (auto &ii : rcoo) if (abs(int(ii) + ofs - int(i)) < space) f = false;
				
			if (f) rcoo.push_back(int(i) - ofs);
		}
	}
}

void tag2(std::string &dna, std::vector<unsigned int> &rcoo) //tagging with regularly spaced tags
{
	int step = 20;// tagging with regular tags with step
	for (unsigned int i = 0; i<dna.length() - lf_len; i+=step)  {
		rcoo.push_back(int(i));
	}
}

int ns;
flannmat read_fastq(long int n)
{
	int l;
	int minl = lf_len + 100;
	std::vector<std::string> loaded_seqs;
	std::string nmer;
	std::deque<unsigned int> FlannId; //

	int fid = 1;

	unsigned int N_reads = 0;
	printf("Started tagging and counting\n");
	
	while ((l = kseq_read(seq)) >=0) {
		std::string ss(seq->seq.s);
		std::vector<unsigned int> * rcoo = new std::vector<unsigned int> ;
		loaded_seqs.push_back(ss);
		len_loaded_seqs.push_back(ss.length());
		if (l < minl ) {
			coords.push_back(NULL);
			continue;
		}
		tag(ss , *rcoo); //replace with tag2 for regular tags
		coords.push_back(rcoo);
		N_reads+=rcoo->size();
	}
	printf("N counted %i\n",N_reads);
	
	ns = N_reads;
	flannmat data(new flanns [N_reads*Alen],N_reads,Alen);
	fflush(stdout);

	int i = 0;

	for (auto ss : loaded_seqs) {
		l = ss.size();
		if (l < minl ) {
			++fid;
			istart.push_back(0);
			continue;
		}
		
		FlannId.clear();
		istart.push_back(i);
		
		for (auto &ii : *coords[fid-1]) {
			ReadId[i] = fid;
			for (int j = 0; j < lf_len; ++j) {
				nmer=ss.substr( int(ii) + j, klen);
				++data[i][Qmer2AIndex[nmer]];
			}
			++i;
		}
		++fid;
	}
	printf("num tags: %d\n", i);
	loaded_seqs.clear();
	return data;
}

template<typename A, typename B> // functions for map sorting from stackoverfow
std::pair<B,A> flip_pair(const std::pair<A,B> &p)
{
    return std::pair<B,A>(p.second, p.first);
}
template<typename A, typename B>
std::multimap<B,A> flip_map(const std::map<A,B> &src)
{
    std::multimap<B,A> dst;
    std::transform(src.begin(), src.end(), std::inserter(dst, dst.begin()), 
                   flip_pair<A,B>);
    return dst;
}

using namespace std;
int main(int argc, char** argv)
{  float version = 1.0;
	int nn,n_iter;
	string f_in, f_out;
	
  char * ivalue = NULL;
  char * ovalue = NULL;
  char * kvalue = NULL;
  char * rvalue = NULL;
  char * nvalue = NULL;
  char * lvalue = NULL;
  char * svalue = NULL;
  char * mvalue = NULL;
  char * wvalue = NULL;

  int c;

  for (int i = 0; i<argc; ++i)  if (std::string(argv[i]) == "-v") {return -1;}	

  opterr = 0;

  while ((c = getopt (argc, argv, "o:i:r:l:s:n:m:w:v:")) != -1)
    switch (c)
      {
      case 'o':
        ovalue = optarg;
        break;
      case 'i':
        ivalue = optarg;
        break;
      case 'n':
        nvalue = optarg;
        break;
      case 'l':
        lvalue = optarg;
        break;
      case 's':
        svalue = optarg;
        break;
      case 'r':
        rvalue = optarg;
        break;
      case 'm':
        mvalue = optarg;
        break;
      case 'w':
        wvalue = optarg;
        break;
      case '?':
	  if (optopt==118) {
		printf ("KD-tree overlapper, Version %g\n",version);
		return -1; 
	  }
          fprintf (stderr,
                   "Unknown option for -%c.\n",
                   optopt);
        return 1;
      default:
        abort ();
      }

      
     string fn;
      if (argc < 2) {
	      printf ("KD-tree overlapper, Version %g\n",version);
	      printf("usage: ./kd [-o output=overlaps.out -k k-mer_len=4 -r kd-tree_iterations=600 -l tags_len=1200 -s tags_spacing=600 -n num_NN=40 -m min_tags_space=200 -w GC_window=100 -v print version and quit] -i input_file  \n");
	      return -1;
	      }
     

	if (nvalue != NULL) nn = atoi(nvalue);
	else nn=40;
	if (kvalue != NULL) klen = atoi(kvalue);
	else klen=4;
	if (lvalue != NULL) lf_len = atoi(lvalue);
	else lf_len=1200;
	if (svalue != NULL) tags_div = atoi(svalue);
	else tags_div=600;
	if (rvalue != NULL) n_iter = atoi(rvalue);
	else n_iter=600;
	if (mvalue != NULL) space = atoi(mvalue);
	else space=200;
	if (wvalue != NULL) cg_win = atoi(wvalue);
	else cg_win=100;
	if (ovalue != NULL) f_out=ovalue;
	else f_out="overlaps.out";
	if (ivalue != NULL) {
		fn=ivalue;
		ifstream f(fn.c_str());
		if (!f.good()) {
			printf("Input file %s is not good\n",fn.c_str());
			return -1;
		}
	}
	else {
		printf("No input file provided.\n");
		return -1;
		
	}
   
	printf ("KD-tree overlapper, Version %g \ninput: %s, output: %s, k_len = %d , tags_len = %d, tags_dens = %d, n_iter = %d\n",
          version,fn.c_str(), f_out.c_str(), klen, lf_len, tags_div, n_iter);
    
	fp = gzopen(fn.c_str(), "r");
	seq = kseq_init(fp);

	populate_kmers_index();
	
	start_timer("Loading fastq...\n");
	datasets.push_back(read_fastq(ns));
	
	printf("Loaded %i tags (%g seconds)\n",ns, stop_timer());
	
	flann::SearchParams params(n_iter, 0.0f,true);
	params.cores = 16; // for multi-cores to work provide the option for GCC "-fopenmp" in the makefile
	
	start_timer("Creating index...\n");
	flann::Index<flann::L1<flanns> > kdindex(datasets.back(), flann::KDTreeIndexParams(4)); // number of trees
	kdindex.buildIndex();
	printf("Created index (%g seconds)\n", stop_timer());

	start_timer("Searching ANNs...\n");

	int ss = ns;
	flann::Matrix<int> indices(new int[ss*nn],ss,nn);
	flann::Matrix<float> dists(new float[ss*nn],ss,nn);
	//SEARCH ####################################
	kdindex.knnSearch(datasets.back(), indices, dists, nn, params);
	printf("Search done (%g seconds)\n", stop_timer());

	printf ("lf_len: %i tags_div %i \n", lf_len, tags_div);
	start_timer("Filtering overlapping pairs...\n");
	
	int hs = 0;
	
	unsigned int rf=1,rf1,rf2,r1,r2;
	
	std::map<int,  int> R2_count;
	std::map< int, double> R2_avg_weight;
	std::map<int, std::vector<std::pair<int,int> > > PairedFrags;

	FILE * mymhap;
	mymhap = fopen(f_out.c_str(),"w");
	
	std::map<int,int> rec;
	for ( int i = 0; i < ns; ++i) //for all tags
	{	
		for ( int j = 0; j < nn; ++j) //for all nearest neighbors of a tag
		{	
			if ( j == 0 and indices[i][j] == i and dists[i][j] == 0 ) continue; //skip hit on itself
			if (ReadId[i] == rf) {
				
				if (rec.count(indices[i][j]) > 0) continue;
				rec[indices[i][j]] = 1;
				R2_count[ReadId[indices[i][j]]]++;
				R2_avg_weight[ReadId[indices[i][j]]]+=dists[i][j];
				PairedFrags[ReadId[indices[i][j]]].push_back(make_pair(i,indices[i][j]));
			}
			else {
				rec.clear();
				std::multimap< int,  int> dst = flip_map(R2_count);
				
				int ii=0;
				r1 = rf;
				for (auto it = dst.rbegin(); it != dst.rend(); it++ ) { // going through paired reads
					
					r2 = it->second;
					int n_hits=PairedFrags[r2].size();
					if (PairedFrags[r2].size() <1) continue; //NUMBER OF HITS
					++ii;
					
					if (r1 == r2) continue;
					
					if (abs(r1 - r2)>100) hs++;
					
					std::vector<int> c1;
					std::vector<int> c2;
					int l1=len_loaded_seqs[r1-1];
					int l2=len_loaded_seqs[r2-1];
					
					for (auto it2 = PairedFrags[r2].begin(); it2 != PairedFrags[r2].end(); it2++) { //fragments of reads
						rf1 = it2->first;
						rf2 = it2->second;
						c1.push_back((*coords[r1-1])[rf1-istart[r1-1]]);
						c2.push_back((*coords[r2-1])[rf2-istart[r2-1]]);	
					}
					
					float m,b,r;
					
					if (linreg(c1,c2,&m,&b,&r) == 0) {
						if (abs(abs(m)-1.)*0<0.5) {
							float fit = m;
							int r1s,r2s,r1e,r2e;
							for (unsigned int i = 0; i < c1.size(); ++i) {
								if (m<0) {
									r1s=c1[i]+c2[i]+lf_len-l2;
									r1e=c1[i]+c2[i]+lf_len;
									r2s=c1[i]+c2[i]+lf_len-l1;
									r2e=c1[i]+c2[i]+lf_len;
								}
								else {
									r1s=c1[i]-c2[i];
									r1e=c1[i]-c2[i]+l2;
									r2s=c2[i]-c1[i];
									r2e=c2[i]-c1[i]+l1;
								}
								r1s = (r1s<0) ? 0 : r1s;
								r1e = (r1e>l1) ? l1 : r1e;
								r2s = (r2s<0) ? 0 : r2s;
								r2e = (r2e>l2) ? l2 : r2e;
								
								fprintf(mymhap,"%i %i %i %f 0 %i %i %i %i %i %i %i %i\n",
								r1,r2,n_hits,fit,r1s,r1e,l1,
								(m < 0) ? 1 : 0,r2s,r2e,l2,ii);
								break;
							}
						}
					}
					if (ii>200) break;
				}
				rf=ReadId[i];
				PairedFrags.clear();
				R2_count.clear();
				PairedFrags[ReadId[indices[i][j]]].push_back(make_pair(i,indices[i][j]));
				++R2_count[ReadId[indices[i][j]]];
			}	
		}
	}
	printf("done in (%g seconds)\n", stop_timer());
	fclose(mymhap);
	return 0;
}
