#ifndef GENE_BED_H
#define GENE_BED_H

#include "readdata.h"
#include <boost/format.hpp>

class tssdist {
 public:
  int all;
  int n1;
  int n5;
  int n10;
  int n100;
  int nover100;
 tssdist(): all(0), n1(0), n5(0), n10(0), n100(0), nover100(0){}
  void inc(int d, status st){
    all++;
    if(st != TSS && !d) nover100++;
    else if(abs(d) <= 1000) n1++;
    else if(abs(d) <= 5000) n5++;
    else if(abs(d) <= 10000) n10++;
    else if(abs(d) <= 100000) n100++;
    else nover100++;
  }
  void print(){
    BPRINT("# Input sites total: %1%, <1kbp: %2%, 1kbp~5kbp: %3%, 5kbp~10kbp: %4%, 10kbp~100kbp: %5%, 100kbp~: %6%\n") % all % n1 % n5 % n10 % n100 % nover100;
  }
};

class gdist {
 public:
  long genome;
  long up;
  long down;
  long genic;
  long exon;
  long intron;
  long inter;
  long conv;
  long div;
  long par;
 gdist(): genome(0), up(0), down(0), genic(0), exon(0), intron(0), inter(0), conv(0), div(0), par(0) {}
  void inc(const int val){
    genome++;
    switch(val) {
    case CONVERGENT: conv++;   break;
    case DIVERGENT:  div++;    break;
    case PARALLEL:   par++;    break;
    case UPSTREAM:   up++;     break;
    case DOWNSTREAM: down++;   break;
    case GENIC:      genic++;  break;
    case EXON:       exon++;   break;
    case INTRON:     intron++; break;
    case INTERGENIC: inter++;  break;
    default: break;
    }
  }
  double ratio(const int val){
    long n;
    switch(val) {
    case CONVERGENT: n=conv;   break;
    case DIVERGENT:  n=div;    break;
    case PARALLEL:   n=par;    break;
    case UPSTREAM:   n=up;     break;
    case DOWNSTREAM: n=down;   break;
    case GENIC:      n=genic;  break;
    case EXON:       n=exon;   break;
    case INTRON:     n=intron; break;
    case INTERGENIC: n=inter;  break;
    default: break;
    }
    return static_cast<double>(n*100/genome);
  }
};

template <class T>
void scanBedGene(T &x, const HashOfGeneDataMap &mp, int updist, int downdist)
{
  int summit = x.bed.summit;

  for(auto itr = mp.at(x.bed.chr).begin(); itr != mp.at(x.bed.chr).end(); ++itr) {
    std::string strand = itr->second.strand;
    int s = itr->second.txStart;
    int e = itr->second.txEnd;
    
    if((strand == "+" && s - updist <= summit && summit <= s) ||
       (strand == "-" && e          <= summit && summit <= e + updist)) {  // upstream
      x.update(UPSTREAM, itr->second);
      break;
    } else if((strand == "+" && e            <= summit && summit <= e + downdist) ||
	      (strand == "-" && s - downdist <= summit && summit <= s)) {  // downstream
      x.update(DOWNSTREAM, itr->second);
    } else if(s <= x.bed.summit && x.bed.summit <= e) {   // genic
      x.update(GENIC, itr->second);
    }
  }
  return;
}

template <class T>
int scanBedConv(T &x, const std::vector<convsite> &vconv)
{
  int on=0;
  int summit = x.bed.summit;
  for (auto conv: vconv) {
    if(conv.chr == x.bed.chr && conv.start <= summit && summit <= conv.end) {
      x.update(conv.st, conv.gene);
      on++;
    }
  }
  return on;
}

std::vector<convsite> gen_convergent(const int, const HashOfGeneDataMap &);

#endif  // GENE_BED_H
