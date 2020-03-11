/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_PROFILE_H_
#define _DD_PROFILE_H_

#include "dd_gv.hpp"
#include "dd_readfile.hpp"

class ReadProfile {
  int32_t stype;
  std::string xlabel;
  std::string Rscriptname;
  std::string Rfigurename;

protected:
  int32_t binsize;
  int32_t nbin;
  int32_t width_from_center;
  int32_t binwidth_from_center;

  int32_t nsites;
  int32_t nsites_skipped;

  std::string RDataname;
  std::unordered_map<std::string, std::vector<double>> hprofile;

  std::vector<genedata> get_garray(const GeneDataMap &mp);

  int32_t isExceedRange(const int32_t posi, const int32_t chrlen) {
    return posi - width_from_center < 0 || posi + width_from_center >= chrlen;
  }

  void WriteValAroundPosi(std::ofstream &out,
			  const SamplePairOverlayed &pair, const vChrArray &vReadArray,
			  const int32_t posi, const std::string &strand);

  double getSumVal(const SamplePairOverlayed &pair,
		   const vChrArray &vReadArray,
		   const int32_t sbin, const int32_t ebin);

public:
  ReadProfile(const DROMPA::Global &p, const int32_t _nbin=0);
  void setOutputFilename(const DROMPA::Global &p);
  void MakeFigure(const DROMPA::Global &p);
  virtual void WriteTSV_EachChr(const DROMPA::Global &p, const chrsize &chr)=0;

  void printHead(const DROMPA::Global &p) {
    for (auto &x: p.samplepair) {
      std::string file(RDataname + "." + x.first.label + ".tsv");
      std::ofstream out(file);
      for (int32_t i=-binwidth_from_center; i<=binwidth_from_center; ++i) {
	out << "\t" << (i*binsize);
      }
      out << std::endl;
      out.close();
    }
  }

  void printNumOfSites() const {
    std::cout << "\n\nthe number of sites: " << nsites << std::endl;
    std::cout << "the number of skipped sites: " << nsites_skipped << std::endl;
  }
};


class ProfileTSS: public ReadProfile {

public:
  explicit ProfileTSS(const DROMPA::Global &p):
    ReadProfile(p) {}

  void WriteTSV_EachChr(const DROMPA::Global &p, const chrsize &chr);

};

class ProfileGene100: public ReadProfile {
  enum {GENEBLOCKNUM=100};

  void outputEachGene(std::ofstream &out, const SamplePairOverlayed &x,
		      const genedata &gene, const vChrArray &vReadArray, int32_t len);

public:
  explicit ProfileGene100(const DROMPA::Global &p):
    ReadProfile(p, GENEBLOCKNUM * 3)
  {}

  void WriteTSV_EachChr(const DROMPA::Global &p, const chrsize &chr);

  void printHead(const DROMPA::Global &p) {
    for (auto &x: p.samplepair) {
      std::string file(RDataname + "." + x.first.label + ".tsv");
      std::ofstream out(file);
      for (int32_t i=-GENEBLOCKNUM; i<GENEBLOCKNUM*2; ++i) {
	out << "\t" << i;
      }
      out << std::endl;
      out.close();
    }
  }
};

class ProfileBedSites: public ReadProfile {
public:
  explicit ProfileBedSites(const DROMPA::Global &p):
    ReadProfile(p) {}

  void WriteTSV_EachChr(const DROMPA::Global &p, const chrsize &chr);
};

#endif /* _DD_PROFILE_H_ */
