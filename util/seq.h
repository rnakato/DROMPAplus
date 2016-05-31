#ifndef SEQ_H
#define SEQ_H

#include <iostream>
#include <vector>
#include <string>
#include <map>

using namespace std;

inline string addchr(string chr) {
  if(chr.find("chr") == string::npos) chr = "chr" + chr;
  return chr;
}

typedef enum{
  STRAND_PLUS,
  STRAND_MINUS,
  STRANDNUM
} Strand;

class range{
 public:
  int start;
  int end;
 range(int s, int e): start(s), end(e) {}
};

class bed{
 public:
  string chr;
  int start;
  int end;
  int summit;
 bed(): start(0), end(0), summit(0) {}
 bed(int s, int e, string c): start(s), end(e) {
   chr = addchr(c);
 }
 bed(vector<string> s): start(stoi(s[1])), end(stoi(s[2])), summit((start + end)/2) {
   chr = addchr(s[0]);
 }
 void print() const {
   cout << chr << "\t" << start << "\t" << end;
 }
 void printHead () const {
   cout << "chromosome\tstart\tend";
 }
 int length() {
   return abs(end - start);
 }
};

class bed12 : public bed {
 public:
  string name;
  int score;
  string strand;
  int thickStart;
  int thickEnd;
  string itemRgb;
  int blockCount;
  int blockSizes;
  int blockStarts;
 bed12(){}
 bed12(vector<string> s): bed(s) {
   int num = s.size();
   if(num > 3)  name        = s[3];
   if(num > 4)  score       = stoi(s[4]);
   if(num > 5)  strand      = s[5];
   if(num > 6)  thickStart  = stoi(s[6]);
   if(num > 7)  thickEnd    = stoi(s[7]);
   if(num > 8)  itemRgb     = s[8];
   if(num > 9)  blockCount  = stoi(s[9]);
   if(num > 10) blockSizes  = stoi(s[10]);
   if(num > 11) blockStarts = stoi(s[11]);
 }
 void print() const {
   cout << chr << "\t" << start << "\t" << end << "\t"
	     << name << "\t" << score << "\t" << strand << "\t"
	     << thickStart << "\t" << thickEnd << "\t" << itemRgb << "\t"
	     << blockCount << "\t" << blockSizes << "\t" << blockStarts;
 }
 void printHead () const {
   cout << "chromosome\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts";
 }
};

class macsxls : public bed {
 public:
  int len;
  int summit;
  double pileup;
  double p;
  double enrich;
  double q;
  string name;

 macsxls(){}
 macsxls(vector<string> s): bed(s) {
   len    = stoi(s[3]);
   summit = stoi(s[4]);
   pileup = stod(s[5]);
   p      = stod(s[6]);
   enrich = stod(s[7]);
   q      = stod(s[8]);
   name   = s[9];
 }
 void print() const {
   cout << chr << "\t" << start << "\t" << end << "\t"
	     << len << "\t" << summit << "\t" << pileup << "\t"
	     << p << "\t" << enrich << "\t" << q << "\t" << name;
 }
 void printHead () const {
   cout << "chromosome\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname";
 }
};

class genedata{
 public:
  string tname;
  string gname;
  string chr;
  int   txStart;   // "Transcription start position"
  int   txEnd;     // "Transcription end position"
  int   cdsStart;  // "Coding region start"
  int   cdsEnd;    // "Coding region end"
  int   exonCount; // "Number of exons"
  string strand;
  vector<range> exon;

  // for Ensembl
  string gsrc;  // gene source
  string tsrc;  // transcript source
  string gtype; // gene biotype
  string ttype; // transcript biotype

  genedata(): txStart(0), txEnd(0), cdsStart(0), cdsEnd(0), exonCount(0) {}

  genedata &operator=(const genedata &ob2){
    tname    = ob2.tname;
    gname    = ob2.gname;
    chr      = ob2.chr;
    txStart  = ob2.txStart;
    txEnd    = ob2.txEnd;
    cdsStart = ob2.cdsStart;
    cdsEnd   = ob2.cdsEnd;
    exonCount= ob2.exonCount;
    strand   = ob2.strand;
    exon     = ob2.exon;
    gsrc     = ob2.gsrc;
    tsrc     = ob2.tsrc;
    gtype    = ob2.gtype;
    ttype    = ob2.ttype;
    return *this;
  }

  int length() const {
    return (txEnd - txStart);
  }
  void printall() const {
    if(this){
      cout << tname << "\t" << gname << "\t" << chr << "\t" << strand << "\t" << txStart << "\t" << txEnd << "\t" << cdsStart << "\t" << cdsEnd << "\t" << exonCount << "\tgene source: " << gsrc << "\ttranscript source: "<< tsrc << "\tgene biotype: "<< gtype << "\ttranscript biotype: "<< ttype << "\t";
      for (auto x: exon) cout << x.start << "-" << x.end << ", ";
    }
  }
  void print() const {
    if(this){
      cout << tname << "\t" << gname << "\t" << strand << "\t" << txStart << "\t" << txEnd << "\t";
    }
  }
};

enum status {INTERGENIC, GENIC, INTRON, EXON, DOWNSTREAM, UPSTREAM, TSS, PARALLEL, DIVERGENT, CONVERGENT};

template <class T>
class bed_gene {
 public:
  status st;
  int d;
  T bed;
  const genedata *gene;
 bed_gene(): st(INTERGENIC), d(0), gene(nullptr) {}
 bed_gene(vector<string> s): st(INTERGENIC), d(0), bed(s), gene(nullptr) {}
  void print() { bed.print();}
  void printWithGene() {
    print();
    if(st == UPSTREAM)        cout << "\tupstream\t";
    else if(st == DOWNSTREAM) cout << "\tdownstream\t";
    else if(st == GENIC)      cout << "\tgenic\t";
    else if(st == INTERGENIC) cout << "\tintergenic\t";
    else if(st == CONVERGENT) cout << "\tconvergent\t";
    else if(st == DIVERGENT)  cout << "\tdivergent\t";
    else if(st == PARALLEL)   cout << "\tparallel\t";
    gene->print();
    cout << endl;
  }
  void printWithTss() {
    print();
    if(st == TSS) {
      cout << "\t" << d << "\t"; 	
      gene->print();
    }
    cout << endl;
  }
  
  void update(const status &pst, const genedata &pgene) {
    if(st < pst){
      st = pst;
      gene = &pgene;
    }
  }
  void update(const status &pst, const genedata *pgene) {
    if(st < pst){
      st = pst;
      gene = pgene;
    }
  }
  void printHead () const {
    bed.printHead();
  }
};

class strRange : public bed {
 public:
  string strand;
  const genedata *gene;
 strRange(): bed(), strand(0), gene(0) {}
 strRange(int s, int e, string c, string str, const genedata &pgene): bed(s,e,c), strand(str), gene(&pgene) {}
};

class convsite : public bed {
 public:
  status st;
  const genedata *gene;
 convsite(int s, int e, string c, const genedata *pgene): bed(s,e,c), gene(pgene) {}
};

class fasta {
 public:
  string name;
  long len, len_mpbl;
  int nbin;
  double p_mpbl;  /* mappability */
  double gcov;    /* genome coverage for bin */
  fasta (string str, int l=0): name(str), len(l), len_mpbl(0), nbin(0), p_mpbl(0), gcov(0) {}
  fasta (vector<string> &v): name(v[0]), len(stoi(v[1])), len_mpbl(0), p_mpbl(0), gcov(0) {}
  void print (){
    cout << name << "\t" << len << "\t" << nbin << "\t" << len_mpbl << "\t"<< p_mpbl << "\t" << gcov << endl;
  }
};

class RefGenome {
 public:
  fasta genome;
  vector<fasta> chr;
  int chrnum;
  int *GCdist;
  long sum_GCdist;
  RefGenome (): genome("Genome"), chrnum(0), sum_GCdist(0) {}
  RefGenome (map<string, int> gt): genome("Genome"), chrnum(0), sum_GCdist(0) {
    for(auto itr = gt.begin(); itr != gt.end(); ++itr){
      fasta v(itr->first, itr->second);
      chr.push_back(v);
      genome.len += itr->second;
    }
    chrnum = chr.size();
  }
  RefGenome (map<string, int> gt, int binsize): RefGenome(gt) {
    for(auto &x: chr) genome.nbin += x.nbin = x.len/binsize +1;
  }
  void print (){
    genome.print();
    for(auto x: chr) x.print();
    cout << "chrnum: " << chrnum << endl;
  }
};

#endif  // SEQ_H
