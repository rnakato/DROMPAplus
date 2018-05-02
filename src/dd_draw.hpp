/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_DRAW_H_
#define _DD_DRAW_H_

#include "dd_class.hpp"
#include "dd_readfile.hpp"
#include "WigStats.hpp"

class ChrArray {
public:
  int32_t binsize;
  int32_t nbin;
  WigArray array;
  WigStats stats;
  int32_t totalreadnum;
  std::unordered_map<std::string, int32_t> totalreadnum_chr;

  ChrArray(){}
  ChrArray(const DROMPA::Global &p,
	   const std::pair<const std::string, SampleInfo> &x,
	   const chrsize &chr):
    binsize(x.second.getbinsize()), nbin(chr.getlen()/binsize +1),
    array(loadWigData(x.first, x.second, chr)),
    stats(nbin, p.thre.pthre_inter),
    totalreadnum(x.second.gettotalreadnum()),
    totalreadnum_chr(x.second.gettotalreadnum_chr())
  {
    clock_t t1,t2;
    t1 = clock();
    if(p.getSmoothing()) array.Smoothing(p.getSmoothing());
    t2 = clock();
    PrintTime(t1, t2, "Smoothing");
    t1 = clock();
    stats.setWigStats(array);
    t2 = clock();
    PrintTime(t1, t2, "WigStats");
    t1 = clock();
    stats.peakcall(array, chr.getname());
    t2 = clock();
    PrintTime(t1, t2, "peakcall");
  }
};

class ChrArrayList {
  void loadSampleData(const DROMPA::Global &p, const chrsize &chr) {
    std::cout << "Loading data..";
    clock_t t1,t2;
    for (auto &x: p.vsinfo.getarray()) {
      t1 = clock();
      arrays[x.first] = ChrArray(p, x, chr);
      t2 = clock();
      PrintTime(t1, t2, "ChrArray new");
    }
  }

public:
  const chrsize &chr;
  std::unordered_map<std::string, ChrArray> arrays;

  ChrArrayList(DROMPA::Global &p, const chrsize &_chr): chr(_chr) {
    loadSampleData(p, chr);
    
#ifdef DEBUG    
    std::cout << "all WigArray:" << std::endl;
    for (auto x: arrays) {
      std::cout << x.first << ", binsize " << x.second.binsize << std::endl;
    }
    std::cout << "all SamplePair:" << std::endl;
    for (auto &x: p.samplepair) {
	std::cout << x.first.argvChIP << "," << x.first.argvInput
		  << ", binsize " << x.first.getbinsize() << std::endl;
    }
#endif
  }
};

class ReadProfile {
  int32_t stype;
  std::string Rscriptname;
  std::string RDataname;
  std::string Rfigurename;
  
protected:
  int32_t binsize;
  int32_t nbin;
  int32_t width_from_center;
  int32_t binwidth_from_center;

  std::unordered_map<std::string, std::vector<double>> hprofile;

  std::vector<genedata> get_garray(const GeneDataMap &mp)
  {
    std::vector<genedata> garray;
    for (auto &m: mp) {
      if (m.second.gtype == "nonsense_mediated_decay" ||
	  m.second.gtype == "processed_transcript" ||
	  m.second.gtype == "retained_intron") continue;
      garray.emplace_back(m.second);
    }
    return garray;
  }

  void addValAroundPosi(const SamplePairOverlayed &pair, const ChrArrayList &arrays,
			const int32_t posi, const std::string &strand,
			const std::string &name, const int32_t chrlen)
  {
    if(posi - width_from_center < 0 || posi + width_from_center >= chrlen) return;
    
    int32_t bincenter(posi/binsize);
    int32_t sbin(bincenter - binwidth_from_center);
    int32_t ebin(bincenter + binwidth_from_center);
    std::ofstream out(RDataname, std::ios::app);
    out << name;
    if (strand == "+") {
      //      int32_t n(0);
      for (int32_t i=sbin; i<=ebin; ++i) {
	//	printf("%d\t%f\n", i, getSumVal(pair, arrays, i, i));
	out << "\t" << getSumVal(pair, arrays, i, i);
      }
      //      for (int32_t i=sbin; i<=ebin; ++i) hprofile.at(pair.first.argvChIP)[n++] = getSumVal(pair, arrays, i, i);
    } else {
      // int32_t n(ebin);
      for (int32_t i=ebin; i>=sbin; --i) {
	//	printf("%d\t%f\n", i, getSumVal(pair, arrays, i, i));
	out << "\t" << getSumVal(pair, arrays, i, i);
      }
      //  for (int32_t i=sbin; i<=ebin; ++i) hprofile.at(pair.first.argvChIP)[n--] = getSumVal(pair, arrays, i, i);
    }
    out << std::endl;
  }

  double getSumVal(const SamplePairOverlayed &pair, const ChrArrayList &arrays, const int32_t sbin, const int32_t ebin) {
    double sum(0);
    for (int32_t i=sbin; i<=ebin; ++i) sum += arrays.arrays.at(pair.first.argvChIP).array[i];
    return sum;
  }
  
public:
  ReadProfile(const DROMPA::Global &p, const int32_t _nbin=0):
    stype(p.prof.stype),
    width_from_center(p.prof.width_from_center)
  {
    for (auto &x: p.samplepair) binsize = x.first.getbinsize();
    binwidth_from_center = width_from_center / binsize;
    
    if(_nbin) nbin = _nbin;
    else nbin = binwidth_from_center * 2 +1;

    std::vector<double> array(nbin, 0);
    for (auto &x: p.samplepair) hprofile[x.first.argvChIP] = array;
  }

  void setOutputFilename(const DROMPA::Global &p) {
    std::string prefix;
    if (stype==0)      prefix = p.getPrefixName() + ".ChIPread";
    else if (stype==1) prefix = p.getPrefixName() + ".Enrichment";
    Rscriptname = prefix + ".R";
    RDataname   = prefix + ".tsv";
    Rfigurename = prefix + ".pdf";
    for (auto &x: {Rscriptname, RDataname, Rfigurename}) unlink(x.c_str());
  }

  void output(const DROMPA::Global &p) {
    //    std::ofstream out();
    
    //      double color[]={CLR_RED, CLR_BLUE, CLR_GREEN, CLR_LIGHTCORAL, CLR_BLACK, CLR_PURPLE, CLR_GRAY3, CLR_OLIVE, CLR_YELLOW3, CLR_LIGHTCYAN, CLR_PINK, CLR_SALMON, CLR_GREEN2, CLR_BLUE3, CLR_PURPLE2, CLR_DARKORANGE};
    std::vector<std::string> posistr = {"",
					"Distance from TSS (bp)",
					"Distance from TTS (bp)",
					"% gene length from TSS",
					"Distance from the peak summit (bp)"};

    //    out << "# bsnum_allowed=%d, bsnum_skipped=%d\n", d->ntotal_profile, d->ntotal_skip) << std::endl;
    /*    for (int32_t i=0; i<p.samplepair.size(); ++i) {
      out << "p" << (i+1) << " <- c(" << std::endl;
      for (int32_t j=0; j<num; ++j) {
	out << sample[i].profile.IP[j];
	if(j != num-1) out << ","; else out << ")" << std::endl;
      }
      out << "p" << (i+1) << "_SE <- c(";
      for (int32_t j=0; j<num; ++j) out << (1.96 * sample[i].profile.SEarray[j]);  // 95% CI
      out << "p" << (i+1) << "_upper <- p" << (i+1) << " + p" << (i+1) << "_SE" << std::endl;
      out << "p" << (i+1) << "_lower <- p" << (i+1) << " + p" << (i+1) << "_SE" << std::endl;
      }*/
  }
  
  virtual void add(const DROMPA::Global &p, const ChrArrayList &arrays, const chrsize &chr)=0;
};

class ProfileBedSites: public ReadProfile {
public:
  ProfileBedSites(const DROMPA::Global &p):
    ReadProfile(p) {}

  void add(const DROMPA::Global &p, const ChrArrayList &arrays, const chrsize &chr) {
    for (auto &vbed: p.anno.vbedlist) {
      for (auto &bed: vbed.getvBed()) {
	if (bed.chr != chr.getname()) continue;
	for (auto &x: p.samplepair) addValAroundPosi(x, arrays, bed.summit, "+", bed.name, chr.getlen());
      }
    }
  }
};

class ProfileTSS: public ReadProfile {

public:
  ProfileTSS(const DROMPA::Global &p):
    ReadProfile(p) {}
  
  void add(const DROMPA::Global &p, const ChrArrayList &arrays, const chrsize &chr) {
    std::string chrname(rmchr(chr.getname()));
    if (p.anno.gmp.find(chrname) == p.anno.gmp.end()) return;

    auto gmp(get_garray(p.anno.gmp.at(chrname)));
    for (auto &gene: gmp) {
      std::cout << gene.gname << std::endl;
      for (auto &x: p.samplepair) addValAroundPosi(x, arrays, gene.txStart, gene.strand, gene.gname, chr.getlen());
    }
  }
};

class ProfileTES: public ReadProfile {
public:
  ProfileTES(const DROMPA::Global &p):
    ReadProfile(p) {}
  
  void add(const DROMPA::Global &p, const ChrArrayList &arrays, const chrsize &chr) {
    std::string chrname(rmchr(chr.getname()));
    if (p.anno.gmp.find(chrname) == p.anno.gmp.end()) return;
    
    auto gmp(get_garray(p.anno.gmp.at(chrname)));
    for (auto &gene: gmp) {
      for (auto &x: p.samplepair) addValAroundPosi(x, arrays, gene.txEnd, gene.strand, gene.gname, chr.getlen());
    }
  }
};

class ProfileGene100: public ReadProfile {
  enum {GENEBLOCKNUM=100};
public:
  ProfileGene100(const DROMPA::Global &p):
    ReadProfile(p, GENEBLOCKNUM * 3)
  {}

  void add(const DROMPA::Global &p, const ChrArrayList &arrays, const chrsize &chr) {
    std::string chrname(rmchr(chr.getname()));
    if (p.anno.gmp.find(chrname) == p.anno.gmp.end()) return;
    
    auto gmp(get_garray(p.anno.gmp.at(chrname)));
    for (auto &gene: gmp) {
      int32_t s,e;
      int32_t len(gene.length());
      double len100(len / (double)GENEBLOCKNUM);
      if(gene.txEnd + len >= chr.getlen() || gene.txStart - len < 0) continue;
      
      for(int32_t i=0; i<nbin; ++i) {
	if (gene.strand == "+") {
	  s = (gene.txStart - len + len100 *i)    / binsize;
	  e = (gene.txStart - len + len100 *(i+1))/ binsize;
	}else{
	  s = (gene.txEnd + len - len100 * (i+1))/ binsize;
	  e = (gene.txEnd + len - len100 * i)    / binsize;
	}
	for (auto &x: p.samplepair) hprofile.at(x.first.argvChIP)[i] = getSumVal(x, arrays, s, e);
      }
    }
  }
};


class Figure {
  const chrsize &chr;
  ChrArrayList arrays;
  std::vector<SamplePairOverlayed> &vsamplepairoverlayed;
  std::vector<bed> regionBed;
  
  void setSamplePairEachRatio(const DROMPA::Global &p, SamplePairEach &pair, const std::string &chrname)
  {
    DEBUGprint("setSamplePairEachRatio");
    if (pair.argvInput == "") return;
    pair.ratio = 1;
    
    switch (p.getNorm()) {
    case 0:
      pair.ratio = 1;
      break;
    case 1:
      pair.ratio = getratio(arrays.arrays.at(pair.argvChIP).totalreadnum, arrays.arrays.at(pair.argvInput).totalreadnum);
      break;
    case 2:
      pair.ratio = getratio(arrays.arrays.at(pair.argvChIP).totalreadnum_chr.at(chrname), arrays.arrays.at(pair.argvInput).totalreadnum_chr.at(chrname));
      break;
    case 3:
      pair.ratio = 1; // NCIS
      break;
    }
#ifdef DEBUG
    std::cout << "ChIP/Input Ratio for chr " << chrname << ": " << pair.ratio << std::endl;
#endif
  }
  
public:
  Figure(DROMPA::Global &p, const chrsize &_chr):
    chr(_chr),
    arrays(p, chr),
    vsamplepairoverlayed(p.samplepair),
    regionBed(p.drawregion.getRegionBedChr(chr.getname()))
  {

    for (auto &x: vsamplepairoverlayed) {
      setSamplePairEachRatio(p, x.first, chr.getname());
      if (x.OverlayExists()) setSamplePairEachRatio(p, x.second, chr.getname());
    }
  }

  int32_t Draw(DROMPA::Global &p) {
    if (p.drawregion.isRegionBed() && !regionBed.size()) return 0;
    std::cout << "Drawing.." << std::endl;
    DrawData(p);
    return 1;
  }

  void DrawData(DROMPA::Global &p);
};

#endif /* _DD_READFILE_H_ */
