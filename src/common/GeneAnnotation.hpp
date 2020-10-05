/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _GENEANNOTATION_HPP_
#define _GENEANNOTATION_HPP_

#include "../../submodules/SSP/common/seq.hpp"
#include "../../submodules/SSP/common/inline.hpp"

enum status {INTERGENIC, GENIC, INTRON, EXON, DOWNSTREAM, UPSTREAM, TSS, PARALLEL, DIVERGENT, CONVERGENT};

class genedata {
 public:
  std::string tname;
  std::string gname;
  std::string tid;
  std::string gid;
  std::string chr;
  int32_t txStart;   // "Transcription start position"
  int32_t txEnd;     // "Transcription end position"
  int32_t cdsStart;  // "Coding region start"
  int32_t cdsEnd;    // "Coding region end"
  int32_t exonCount; // "Number of exons"
  std::string strand;
  std::vector<range> exon;

  // for Ensembl
  std::string gsrc;  // gene source
  std::string tsrc;  // transcript source
  std::string gtype; // gene biotype
  std::string ttype; // transcript biotype
  std::string ttag;  // Gencode tag

  genedata(): txStart(0), txEnd(0), cdsStart(0), cdsEnd(0), exonCount(0) {}

  int32_t length() const { return (txEnd - txStart); }
  void printall() const {
    // if(this){
      std::cout << tname << "\t" << gname << "\t" << tid << "\t" << gid << "\t" << chr << "\t" << strand << "\t" << txStart << "\t" << txEnd << "\t" << cdsStart << "\t" << cdsEnd << "\t" << exonCount << "\tgene source: " << gsrc << "\ttranscript source: "<< tsrc << "\tgene biotype: "<< gtype << "\ttranscript biotype: "<< ttype  << "\ttranscript tag: "<< ttag << "\t";
      for (auto &x: exon) std::cout << x.start << "-" << x.end << ", ";
      //    }
  }
  void print() const {
    std::cout << tname << "\t" << gname << "\t" << tid << "\t" << gid << "\t" << strand << "\t" << txStart << "\t" << txEnd << "\t";
  }
};

class hitGene {
 public:
  status st;
  int32_t d;
  const genedata *gene;
 hitGene(): st(INTERGENIC), d(0), gene(nullptr) {}
};

template <class T>
class bed_gene {
 public:
  T bed;
  hitGene gene;
  std::vector<hitGene> genelist;
 bed_gene(): gene() {}
 explicit bed_gene(std::vector<std::string> &s): bed(s), gene() {}
  void print() const { bed.print();}
  void printWithGene(bool redundant) const {
    if(redundant) {
      for(auto &x: genelist) {
	print();
	if(x.st == UPSTREAM)        std::cout << "\tupstream\t";
	else if(x.st == DOWNSTREAM) std::cout << "\tdownstream\t";
	else if(x.st == GENIC)      std::cout << "\tgenic\t";
	else if(x.st == INTERGENIC) std::cout << "\tintergenic\t";
	else if(x.st == CONVERGENT) std::cout << "\tconvergent\t";
	else if(x.st == DIVERGENT)  std::cout << "\tdivergent\t";
	else if(x.st == PARALLEL)   std::cout << "\tparallel\t";
	else PRINTERR_AND_EXIT("\nError: invalid status: " << gene.st);
	if(gene.st != INTERGENIC) x.gene->print();
	std::cout << std::endl;
      }
    } else {
      print();
      if(gene.st == UPSTREAM)        std::cout << "\tupstream\t";
      else if(gene.st == DOWNSTREAM) std::cout << "\tdownstream\t";
      else if(gene.st == GENIC)      std::cout << "\tgenic\t";
      else if(gene.st == INTERGENIC) std::cout << "\tintergenic\t";
      else if(gene.st == CONVERGENT) std::cout << "\tconvergent\t";
      else if(gene.st == DIVERGENT)  std::cout << "\tdivergent\t";
      else if(gene.st == PARALLEL)   std::cout << "\tparallel\t";
      else PRINTERR_AND_EXIT("\nError: invalid status: " << gene.st);
      if(gene.st != INTERGENIC) gene.gene->print();
      std::cout << std::endl;
    }
  }
  void printWithTss(bool redundant) const {
    if(redundant) {
      for(auto &x: genelist) {
	print();
	std::cout << "\t" << x.d << "\t";
	x.gene->print();
	std::cout << std::endl;
      }
    } else {
      print();
      if(gene.st == TSS) {
	std::cout << "\t" << gene.d << "\t";
	gene.gene->print();
      }
      std::cout << std::endl;
    }
  }

  void update(const status &pst, const genedata &pgene) {
    if(gene.st < pst){
      gene.st = pst;
      gene.gene = &pgene;
    }
    hitGene g;
    g.st = pst;
    g.gene = &pgene;
    genelist.push_back(g);
  }
  void update(const status &pst, const genedata *pgene) {
    if(gene.st < pst){
      gene.st = pst;
      gene.gene = pgene;
    }
    hitGene g;
    g.st = pst;
    g.gene = pgene;
    genelist.push_back(g);
  }
  void printHead () const { bed.printHead(); }
};

#endif  // _GENEANNOTATION_HPP_
