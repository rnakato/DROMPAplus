/* Copyright(c)  Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_PROFILE_H_
#define _DD_PROFILE_H_

#include "dd_gv.hpp"
#include "dd_readfile.hpp"

class ReadProfile {
  std::string xlabel;
  std::string Rscriptname;
  std::string Rfigurename;

protected:
  int32_t stype;
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
                          const SamplePairOverlayed &pair,
                          const vChrArray &vReadArray,
                          const int32_t posi,
                          const std::string &strand);

  double getAverageVal(const SamplePairOverlayed &pair,
                       const vChrArray &vReadArray,
                       const int32_t sbin, const int32_t ebin);
  double getMaxVal(const SamplePairOverlayed &pair,
                   const vChrArray &vReadArray,
                   const int32_t sbin, const int32_t ebin);

public:
  ReadProfile(const DROMPA::Global &p, const int32_t _nbin=0);
  void setOutputFilename(const DROMPA::Global &p, const std::string &);
  void MakeFigure(const DROMPA::Global &p);
  virtual void WriteTSV_EachChr(const DROMPA::Global &p, const chrsize &chr)=0;

  virtual void printHead(const DROMPA::Global &p) {
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

  void printNumOfSites(const int32_t nsample) const {
    std::cout << "\n\nthe number of sites: " << nsites/nsample << std::endl;
    std::cout << "the number of skipped sites: " << nsites_skipped/nsample << std::endl;
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
                      const genedata &gene, const vChrArray &vReadArray, const int32_t len);
  void outputEachGene_fixedlength(std::ofstream &out, const SamplePairOverlayed &x,
                                  const genedata &gene, const vChrArray &vReadArray, const int32_t len,
                                  const int32_t width_from_gene);


public:
  explicit ProfileGene100(const DROMPA::Global &p):
    ReadProfile(p, GENEBLOCKNUM * 3) {}

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

class ProfileMULTICI: public ReadProfile {
public:
  explicit ProfileMULTICI(const DROMPA::Global &p):
    ReadProfile(p) {}

  void WriteTSV_EachChr(const DROMPA::Global &p, const chrsize &chr);

  void printHead(const DROMPA::Global &p) {
    std::string file(RDataname + ".tsv");
    std::ofstream out(file);
    out << "chromosome\tstart\tend";
    if (p.isaddname())  out << "\tname";

    for (auto &x: p.samplepair) out << "\t" << x.first.label;
    out << std::endl;
    out.close();
  }
};


template <class T>
void exec_PROFILE_asType(DROMPA::Global &p)
{
  T profile(p);
  profile.setOutputFilename(p, "PROFILE");
  profile.printHead(p);

  for(auto &chr: p.gt) {
    if(!p.isincludeYM() && (chr.getname() == "Y" || chr.getname() == "M")) continue;
    std::cout << "\nchr" << chr.getname() << "..";

    profile.WriteTSV_EachChr(p, chr);
  }

  profile.printNumOfSites(p.samplepair.size());
  profile.MakeFigure(p);

  return;
}

#endif /* _DD_PROFILE_H_ */
