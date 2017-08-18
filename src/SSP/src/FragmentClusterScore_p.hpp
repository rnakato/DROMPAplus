/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _FRAGMENTCLUSTERSCORE_P_H_
#define _FRAGMENTCLUSTERSCORE_P_H_

class PropNeighborFrag {
  double sumOfvNeighborFrag;
  std::vector<int32_t> vNeighborFrag;
  std::vector<double> vcPNF;

 public:
 PropNeighborFrag(): sumOfvNeighborFrag(0),
		     vNeighborFrag(sizeOfvNeighborFrag, 0),
		     vcPNF(sizeOfvNeighborFrag, 0) {}

  void setNeighborFrag(const int32_t flen, const int32_t end, const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev);

  double getPNF(const int32_t i) const {
    if(i<0 || i>= sizeOfvNeighborFrag) {
      std::cerr << "error: invalid num " << i << "for getPNF max: " << sizeOfvNeighborFrag << std::endl;
      return -1;
    }
    else {
      return getratio(vNeighborFrag[i], sumOfvNeighborFrag);
    }
  }
  double getCumulativePNF(const int32_t i) const {
    if(i<0 || i>= sizeOfvNeighborFrag) {
      std::cerr << "error: invalid num " << i << "for getCumulativePNF max: " << sizeOfvNeighborFrag << std::endl;
      return -1;
    }
    else {
      return vcPNF[i];
    }
  }
};

class shiftFragVar {
  std::map<int32_t, double> mpFCS;
  //  std::map<int32_t, double> mpKLD;
  std::map<int32_t, PropNeighborFrag> pnf;
  std::map<int32_t, PropNeighborFrag> pnfbg;
  int32_t lenF3;
  int32_t flen;
  double r4cmp;
  int32_t numUsed4FCS;
  bool lackOfReads;
  int32_t ng_from_fcs;
  int32_t ng_to_fcs;
  int32_t ng_step_fcs;

  long nread;

 public:
  shiftFragVar(const FCSstats &fcsst, const SeqStatsGenome &genome);

  void execchr(const SeqStats &chr);
  uint32_t getnumUsed4FCS() const { return numUsed4FCS; }
  double getpnfbg(const int32_t i) const {
    double v(0);
    for(auto x: pnfbg) v += x.second.getCumulativePNF(i);
    return v/pnfbg.size();
  }
  /*  double getKLD(const int32_t len) const { // Kullback-Leibler divergence
    double e(0);
    for(size_t k=0; k<sizeOfvNeighborFrag; ++k) {
      double p(pnf.at(len).getCumulativePNF(k));
      double q(getpnfbg(k));
      if(p && q) e += p * log2(p/q);
    }
    return e;
    }*/
  double getFCS(const int32_t len) const {
    double diffMax(0);
    for(size_t k=0; k<sizeOfvNeighborFrag; ++k) {
      diffMax = std::max(diffMax, pnf.at(len).getCumulativePNF(k) - getpnfbg(k));
    }
    return diffMax;
  }
  void calcFCS() {
    std::cout << "Calculate FCS score..." << std::flush;
    for(auto x: pnf) {
      std::cout << x.first << "..." << std::flush;
      mpFCS[x.first] = getFCS(x.first);
      //      mpKLD[x.first] = getKLD(x.first);
    }
    std::cout << "done." << std::endl;
  }
  void outputPnf(const std::string &filename) const {
    std::ofstream out(filename);

    for(auto &x: v4pnf) out << "\tPNF len"  << x;
    for(auto &x: v4pnf) out << "\tCPNF len" << x;
    out << std::endl;

    for(size_t k=0; k<sizeOfvNeighborFrag-1; ++k) {
      out << k << "\t";
      for(auto &x: v4pnf) out << pnf.at(x).getPNF(k)           << "\t";
      for(auto &x: v4pnf) out << pnf.at(x).getCumulativePNF(k) << "\t";
      out << std::endl;
    }
  }

  void setFCSstats(FCSstats &p) {
    p.setfcsread(mpFCS.at(lenF3));
    p.setfcsflen(mpFCS.at(flen));
    p.setfcs1k(mpFCS.at(1000));
    p.setfcs10k(mpFCS.at(10000));
    p.setfcs100k(mpFCS.at(100000));
  }

  void outputFCS(const std::string &filename) const {
    if(!nread) std::cerr << filename << ": no read" << std::endl;

    std::ofstream out(filename);
    std::string str("");
    if(lackOfReads) str = " (read number is insufficient)";
    out << "Read length"     << str << "\t" << mpFCS.at(lenF3) << std::endl;
    out << "Fragment length" << str << "\t" << mpFCS.at(flen)  << std::endl;
    out << "Broad (1 kbp)"   << str << "\t" << mpFCS.at(1000)  << std::endl;
    out << "Broad (10 kbp)"  << str << "\t" << mpFCS.at(10000) << std::endl;
    out << "Strand shift\tFragment cluster score" << std::endl;
    for(auto x: mpFCS) {
      out << x.first << "\t" << mpFCS.at(x.first) << std::endl;
    }
  }
  /*  void outputKLD(const std::string &filename) const {
    if(!nread) std::cerr << filename << ": no read" << std::endl;

    std::ofstream out(filename);
    std::string str("");
    if(lackOfReads) str = " (read number is insufficient)";
    out << "Read length"     << str << "\t" << mpKLD.at(lenF3) << std::endl;
    out << "Fragment length" << str << "\t" << mpKLD.at(flen)  << std::endl;
    out << "Broad (1 kbp)"   << str << "\t" << mpKLD.at(1000)  << std::endl;
    out << "Broad (10 kbp)"  << str << "\t" << mpKLD.at(10000) << std::endl;
    out << "Strand shift\tFragment cluster score" << std::endl;
    for(auto x: mpKLD) {
      out << x.first << "\t" << mpKLD.at(x.first) << std::endl;
    }
    }*/
};

#endif /* _FRAGMENTCLUSTERSCORE_P_H_ */
