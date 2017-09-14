/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_CLASS_H_
#define _DD_CLASS_H_

#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include "WigStats.hpp"
#include "SSP/common/BoostOptions.hpp"
#include "SSP/common/BedFormat.hpp"
#include "SSP/common/ReadAnnotation.hpp"
#include "SSP/common/util.hpp"

class chrsize;
class SamplePairChr;

enum class DrompaCommand {CHIP, NORM, THRE, ANNO_PC, ANNO_GV, DRAW, REGION, CG, PD, TR, PROF, OTHER};

class CommandParamSet {
public:
  int32_t sm;
  int32_t showctag;
  int32_t showratio;
  int32_t scaletag;
  int32_t scaleratio;
  int32_t scalepvalue;
  bool sigtest;

  CommandParamSet(const int32_t defsm,
		  const int32_t defshowctag,
		  const int32_t defshowratio,
		  const int32_t defscaletag,
		  const int32_t defscaleratio,
		  const int32_t defscalepvalue,
		  const bool defsigtest):
    sm(defsm),
    showctag(defshowctag),
    showratio(defshowratio),
    scaletag(defscaletag),
    scaleratio(defscaleratio),
    scalepvalue(defscalepvalue),
    sigtest(defsigtest)
  {}
};

class SampleFile {
  double lambda;
  double nb_p, nb_n, nb_p0;
  WigType iftype;
  int32_t binsize;
  int32_t totalreadnum;
  std::unordered_map<std::string, int32_t> totalreadnum_chr;
  std::string prefix;
public:
  std::vector<int> data;

  void setbinsize(std::string &v, const int32_t b) {
    if(b>0) binsize = b;
    else {
      try {
	binsize = stoi(v);
      }catch (...) {
	binsize = 0;
      }
    }
    if(binsize <= 0) PRINTERR("invalid binsize: " << v);
  }
  
 SampleFile() {}
  SampleFile(const std::string &filename, const std::vector<chrsize> gt,
	     const int32_t b, const WigType &type):
    binsize(0), totalreadnum(0), prefix("")
  {
   std::vector<std::string> v;
   boost::split(v, filename, boost::algorithm::is_any_of("."));
   int last(v.size()-1);

   if (type != WigType::NONE) iftype = type;
   else {
     if(v[last] == "wig") iftype = WigType::UNCOMPRESSWIG;    
     else if(v[last] == "gz" && v[last-1] == "wig") {
       iftype = WigType::COMPRESSWIG;
       --last;
     } else if(v[last] == "bedGraph") iftype = WigType::BEDGRAPH;
     else if(v[last] == "bw")         iftype = WigType::BIGWIG;
     else if(v[last] == "bin")        iftype = WigType::BINARY;
     else PRINTERR("invalid postfix: " << filename);
   }
   setbinsize(v[last-1], b);
   for (int32_t i=0; i<last; ++i) prefix += v[i] + ".";
   gettotalreadnum(filename, gt);
  }
  
  void scanStatsFile(const std::string &filename);
  void gettotalreadnum(const std::string &filename, const std::vector<chrsize> gt);
  int32_t getbinsize() const { return binsize; }
  WigType getiftype() const { return iftype; }

  int32_t gettotalreadnum() const { return totalreadnum; }
  const std::unordered_map<std::string, int32_t> & gettotalreadnum_chr() const { return totalreadnum_chr;}
  
};

class yScale {
  enum {TAG_DEFAULT=30, RATIO_DEFAULT=3, P_DEFAULT=5};
 public:
  double tag;
  double ratio;
  double pvalue;
 yScale(): tag(TAG_DEFAULT), ratio(RATIO_DEFAULT), pvalue(P_DEFAULT) {}
};

class SamplePairParam {
  int32_t binsize;
  yScale scale;
  std::unordered_map<std::string, std::vector<bed>> peaks;
  
public:
  std::string argvChIP, argvInput;
  std::string peak_argv;
  std::string label;
  double ratio;
  
  SamplePairParam():
    binsize(0), argvChIP(""), argvInput(""), peak_argv(""), label(""), ratio(1)
  {}
  SamplePairParam(const std::string &str, const int32_t b):
    binsize(0), argvChIP(""), argvInput(""), peak_argv(""), label(""), ratio(1)
  {
    std::vector<std::string> v;
    boost::split(v, str, boost::algorithm::is_any_of(","));

    if(v[0] != "") argvChIP  = v[0];
    if(v.size() >=2 && v[1] != "") argvInput = v[1];
    if(v.size() >=3 && v[2] != "") label     = v[2];
    if(v.size() >=4 && v[3] != "") peak_argv = v[3];
    if(peak_argv != "") peaks = parseBed_Hash<bed>(peak_argv);
    binsize = b;
    if(v.size() >=6 && v[5] != "") scale.tag = stod(v[5]);
    if(v.size() >=7 && v[6] != "") scale.ratio = stod(v[6]);
    if(v.size() >=8 && v[7] != "") scale.pvalue = stod(v[7]);

    //    printBed_Hash(peaks);
  }
  std::vector<bed> getpeaksChr(const std::string &chrname) const {
    if (peak_argv != "") return peaks.at(rmchr(chrname));
    else return std::vector<bed>();
  }
  void print() const {
    std::cout << boost::format("ChIP: %1% label: %2% peaklist: %3%\n") % argvChIP % label % peak_argv;
    std::cout << boost::format("   Input: %1%\n") % argvInput;
    std::cout << boost::format("   binsize: %1%\n") % binsize;
  }
  int32_t getbinsize() const { return binsize; }
};
  
class SamplePair {
 public:
  bool overlay;
  SamplePairParam first;
  SamplePairParam second;

  SamplePair(const std::string &str, const int32_t b):
    overlay(false), first(str, b)
  {}
  void print() const {
    first.print();
    if (overlay) {
      printf("Overlay ");
      second.print();
    }
  }
};

class pdSample {
 public:
  std::string argv;
  std::string name;
  pdSample(){}
};

class GraphFile {
  std::string filename;
  int32_t binsize;

public:
  GraphFile(): filename(""), binsize(0) {}

  void setValue(const std::string &f, const int32_t s) {
    filename = f;
    binsize = s;
  }
  const std::string & getfilename() const { return filename; }
  int32_t getbinsize() const { return binsize; }
  bool isOn() const {
    if(filename != "") return true;
    else return false;
  }
};

namespace DROMPA {
  class Global;
  
  class Annotation {

  public:
    std::string genefile;
    int32_t gftype;
    HashOfGeneDataMap gmp;
    bool showtranscriptname;
    std::string arsfile;
    bool showars;
    std::string terfile;
    std::vector<std::vector<bed>> vbedlist;
    std::vector<InteractionSet> vinterlist;
    std::string repeatfile;
    std::string mpfile;
    double mpthre;
    std::string gapfile;
    GraphFile GC;
    GraphFile GD;

    Annotation():
      genefile(""), arsfile(""), terfile(""),
      repeatfile(""), mpfile(""), gapfile("")
    {}

    void setOptsPC(MyOpt::Opts &allopts);
    void setOptsGV(MyOpt::Opts &allopts);

    HashOfGeneDataMap getGMP() {
      HashOfGeneDataMap tmp;
      if(!gftype)        tmp = parseRefFlat(genefile);
      else if(gftype==1) tmp = parseGtf(genefile);
      else if(gftype==2) tmp = parseSGD(genefile);
      else PRINTERR("invalid --gftype: " << gftype);
      
      //      printMap(tmp);
      
	return tmp; // hash for transcripts
    }
    
    void setValuesPC(const MyOpt::Variables &values);
    void setValuesGV(const MyOpt::Variables &values);
    void InitDumpPC(const MyOpt::Variables &values) const;
    void InitDumpGV(const MyOpt::Variables &values) const;
    
    int32_t getgftype() const { return gftype; }
  };
  
  class Threshold {
  public:
    double pthre_inter;
    double pthre_enrich;
    double qthre;
    double ethre;
    double ipm;
    bool sigtest;
    //    int32_t width4lmd;
    
    Threshold(): sigtest(false) {}

    void setOpts(MyOpt::Opts &allopts, const bool sig);
    void setValues(const MyOpt::Variables &values);
    void InitDump() const;
  };
  
  class DrawRegion {
    bool isRegion;
    std::vector<bed> regionBed;

    std::string chr;
    std::string genelocifile;
    int32_t len_geneloci;

  public:
    DrawRegion():
      isRegion(false), chr(""), genelocifile(""), len_geneloci(0)
    {}

    void setOpts(MyOpt::Opts &allopts);
    void setValues(const MyOpt::Variables &values);
    void InitDump(const MyOpt::Variables &values) const;
    
    std::vector<bed> getRegionBedChr(const std::string &chrname) {
      std::vector<bed> vbed;
      for(auto &x: regionBed) {
	if (x.chr == chrname || x.chr == "chr" + chrname) vbed.emplace_back(x);
      }
      return vbed;
    }
    const std::string & getchr() const { return chr; }
    bool isRegionBed() const { return isRegion; }
    void isRegionOff() { isRegion=false; }
  };
  
  class DrawParam {
    int32_t linenum_per_page;
    int32_t barnum;
    double ystep;
    bool showymem;
    bool showylab;

    int32_t samplenum;

  public:
    int32_t width_per_line;
    int32_t showctag;
    int32_t showitag;
    int32_t showratio;
    int32_t showpinter;
    int32_t showpenrich;
    int32_t viz;
    double scale_tag;
    double scale_ratio;
    double scale_pvalue;

    DrawParam(): showymem(true), showylab(true) {}

    void setOpts(MyOpt::Opts &allopts, const CommandParamSet &cps);
    void setValues(const MyOpt::Variables &values, const int32_t n);
    void InitDump() const;
    bool isshowymem() const { return showymem; };
    bool isshowylab() const { return showylab; };

    double getlineheight() const { return ystep * barnum; }

    int32_t getNumLine(const int32_t s, const int32_t e) const{
      return (e-s)/width_per_line +1;
    }
    int32_t getNumPage(const int32_t s, const int32_t e) const {
      return (getNumLine(s,e) -1) / linenum_per_page +1;
    }

    int32_t getlpp() const { return linenum_per_page; }
    int32_t getHeightEachSample(const SamplePairParam &pair) const;
    int32_t getHeightAllSample(const Global &p, const std::vector<SamplePairChr> &pairs) const;
    int32_t getPageHeight(const Global &p, const std::vector<SamplePairChr> &pairs) const;
  };
  
  class Global {
    bool ispng;
    bool showchr;
    WigType iftype;
    std::string oprefix;
    bool includeYM;
    int32_t norm;
    int32_t smoothing;

  public:
    MyOpt::Opts opts;
    DrawParam drawparam;
    DrawRegion drawregion;
    Threshold thre;
    Annotation anno;

    std::vector<chrsize> gt;
    std::unordered_map<std::string, SampleFile> sample;
    std::vector<SamplePair> samplepair;
    std::vector<pdSample> pd;

    bool isGV;
    
    Global():
      ispng(false), showchr(false), iftype(WigType::NONE),
      oprefix(""), includeYM(false),
      opts("Options"), isGV(false)
    {}

    void setOpts(const std::vector<DrompaCommand> &st, const CommandParamSet &cps);
    void setValues(const std::vector<DrompaCommand> &vopts, const MyOpt::Variables &values);
    void setOptsNorm(MyOpt::Opts &allopts, const int32_t defsm);
    void setValuesNorm(const MyOpt::Variables &values);
    void setOptsOther(MyOpt::Opts &allopts);
    void setValuesOther(const MyOpt::Variables &values);
    void InitDumpChIP() const;
    void InitDumpNorm() const;
    void InitDumpOther() const;

    WigType getIfType() const { return iftype; }

    int32_t getSmoothing() const { return smoothing; }
    int32_t getNorm() const { return norm; }
    const std::string getFigFileName() const
    {
      return oprefix + ".pdf";
    }
    const std::string getFigFileNameChr(const std::string &chr) const
    {
      return oprefix + "_" + chr + ".pdf";
    }
    bool isincludeYM() const { return includeYM; }
    bool isshowchr() const { return showchr; }
  };

}

#endif /* _DD_CLASS_H_ */
