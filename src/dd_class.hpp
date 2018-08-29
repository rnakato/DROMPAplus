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

enum class DrompaCommand {CHIP, NORM, THRE, ANNO_PC, ANNO_GV, DRAW, REGION, CG, PD, TR, PROF, OTHER};

class CommandParamSet {
public:
  int32_t sm;
  int32_t showctag;
  int32_t showratio;
  double scaletag;
  double scaleratio;
  double scalepvalue;
  bool sigtest;

  CommandParamSet(const int32_t defsm,
		  const int32_t defshowctag,
		  const int32_t defshowratio,
		  const double defscaletag,
		  const double defscaleratio,
		  const double defscalepvalue,
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

class SampleInfo {
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
  
  SampleInfo() {}
  SampleInfo(const std::string &filename,
	     const std::vector<chrsize> gt,
	     const int32_t b,
	     const WigType &type):
    binsize(0), totalreadnum(0), prefix("")
  {
   std::vector<std::string> v;
   ParseLine(v, filename, '.');
   int last(v.size()-1);

   if (type != WigType::NONE) iftype = type;
   else {
     if(v[last] == "wig") iftype = WigType::UNCOMPRESSWIG;    
     else if(v[last] == "gz" && v[last-1] == "wig") {
       iftype = WigType::COMPRESSWIG;
       --last;
     } else if(v[last] == "bedGraph") iftype = WigType::BEDGRAPH;
     else if(v[last] == "bw")         iftype = WigType::BIGWIG;
     //     else if(v[last] == "bin")        iftype = WigType::BINARY;
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
  const std::unordered_map<std::string, int32_t>& gettotalreadnum_chr() const& { return totalreadnum_chr; }
};


class SampleInfoList {
  std::unordered_map<std::string, SampleInfo> vsinfo;

public:
  SampleInfoList(){}
  
  void addSampleInfo(const std::string &str, const std::vector<chrsize> &gt, const WigType iftype) {
    int32_t binsize(0);  
    std::vector<std::string> v;
    ParseLine(v, str, ',');
      
    if(v.size() >8) {
      std::cerr << "error: sample std::string has ',' more than 8: " << str << std::endl;
      exit(1);
    }
    if(v[0] == "") {
      std::cerr << "please specify ChIP sample: " << str << std::endl;
      exit(1);
    }
    isFile(v[0]);
      
    if(v.size() >4 && v[4] != "") {
      try { binsize = stoi(v[4]); }
      catch (...) { std::cerr << "Warning: invalid binsize " << v[4] << "." << std::endl; }
    }
      
    // ChIP sample
    if(!Exists(v[0])) vsinfo[v[0]] = SampleInfo(v[0], gt, binsize, iftype);
    if(vsinfo[v[0]].getbinsize() <= 0) PRINTERR("please specify binsize.\n");
      
    // Input sample
    if(v.size() >=2 && v[1] != "") {
      if(!Exists(v[1])) vsinfo[v[1]] = SampleInfo(v[1], gt, binsize, iftype);
      if(vsinfo[v[0]].getbinsize() != vsinfo[v[1]].getbinsize()) PRINTERR("binsize of ChIP and Input should be same. " << str);
    }
  }

  bool Exists(const std::string &str) const { return vsinfo.find(str) != vsinfo.end(); }
  int32_t getbinsize(const std::string &str) const { return vsinfo.at(str).getbinsize(); }
  
  const SampleInfo& operator[](const std::string &str) const& {
    return vsinfo.at(str);
  }
  const std::unordered_map<std::string, SampleInfo> &getarray() const & { return vsinfo; }
};
  

class yScale {
 public:
  double tag;
  double ratio;
  double pvalue;
   yScale(): tag(0), ratio(0), pvalue(0) {}
};

class SamplePairEach {
  int32_t binsize;
  std::unordered_map<std::string, std::vector<bed>> peaks;
  
public:
  std::string argvChIP, argvInput;
  std::string peak_argv;
  std::string label;
  double ratio;
  yScale scale;
  
  SamplePairEach():
    binsize(0), argvChIP(""), argvInput(""), peak_argv(""), label(""), ratio(1)
  {}
  SamplePairEach(const std::string &str,  const SampleInfoList &vsinfo):
    binsize(0), argvChIP(""), argvInput(""), peak_argv(""), label(""), ratio(1)
  {
    std::vector<std::string> v;
    ParseLine(v, str, ',');

    /* 1:ChIP   2:Input   3:name   4:peaklist   5:binsize
       6:scale_tag   7:scale_ratio   8:scale_pvalue */
    if(v[0] != "") argvChIP = v[0];
    if(v.size() >=2 && v[1] != "") argvInput = v[1];
    if(v.size() >=3 && v[2] != "") label     = v[2];
    if(v.size() >=4 && v[3] != "") peak_argv = v[3];
    if(peak_argv != "") peaks = parseBed_Hash<bed>(peak_argv);
    binsize = vsinfo.getbinsize(argvChIP);
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
  bool InputExists() const { return argvInput != "";}
};
  
class SamplePairOverlayed {
  bool overlay;
  
 public:
  SamplePairEach first;
  SamplePairEach second;

  SamplePairOverlayed(const std::string &str, const SampleInfoList &vsinfo):
    overlay(false), first(str, vsinfo)
  {}

  void setSecondSample(const std::string &str, const SampleInfoList &vsinfo) {
    second = SamplePairEach(str, vsinfo);
    overlay = true;
  }

  void print() const {
    first.print();
    if (overlay) {
      printf("Overlay ");
      second.print();
    }
  }
  bool OverlayExists() const { return overlay; }
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

    template <class T>
    void readBedFile(const std::vector<std::string> &v) {
      auto vbed = parseBed<T>(v[0]);
      //    printBed(vbed);
      if(v.size()>1) vbedlist.emplace_back(vbed, v[1]);
      else vbedlist.emplace_back(vbed, "Bed");
    }
    
  public:
    std::string genefile;
    int32_t gftype;
    HashOfGeneDataMap gmp;
    bool showtranscriptname;
    std::string arsfile;
    bool showars;
    std::string terfile;
    std::vector<vbed<bed12>> vbedlist;
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
  
  class Profile {
    enum {TSS, TTS, GENE100, BEDSITES, PTYPENUM};

  public:
    int32_t ptype;
    int32_t stype;
    int32_t ntype;
    int32_t width_from_center;
    double maxval;
    int32_t hmsort;
    
    Profile(){}
    
    void setOpts(MyOpt::Opts &allopts);
    void setValues(const MyOpt::Variables &values);
    void InitDump() const;
    bool isPtypeTSS() const { return ptype == TSS; }
    bool isPtypeTTS() const { return ptype == TTS; }
    bool isPtypeGene100() const { return ptype == GENE100; }
    bool isPtypeBed() const { return ptype == BEDSITES; }
  };
  
  class Threshold {
  public:
    double pthre_inter;
    double pthre_enrich;
    double qthre;
    double ethre;
    double ipm;
    bool sigtest;
    
    Threshold(): sigtest(false) {}

    void setOpts(MyOpt::Opts &allopts);
    void setValues(const MyOpt::Variables &values);
    void InitDump() const;
  };
  
  class DrawRegion {
    bool isRegion;
    std::vector<bed> regionBed;

    std::string chr;
    std::unordered_map<std::string, int32_t> geneloci;
    int32_t len_geneloci;

    void getGeneLoci(const std::string &genelocifile) {
        std::ifstream in(genelocifile);
	if (!in) PRINTERR("cannot open " << genelocifile);
	
	std::string lineStr;
	while (!in.eof()) {
	  getline(in, lineStr);
	  if (lineStr.empty()) continue;
	  std::vector<std::string> v;
	  ParseLine(v, lineStr, '\t');
	  geneloci[v[0]] = 1;
	}
    }

  public:
    DrawRegion():
      isRegion(false), chr(""), len_geneloci(0)
    {}

    void setOpts(MyOpt::Opts &allopts);
    void setValues(const MyOpt::Variables &values);
    void InitDump(const MyOpt::Variables &values) const;
    
    std::vector<bed> getRegionBedChr(const std::string &chrname) const {
      std::vector<bed> vbed;
      for(auto &x: regionBed) {
	if (x.chr == chrname || x.chr == "chr" + chrname) vbed.emplace_back(x);
      }
      return vbed;
    }
    const std::string & getchr() const { return chr; }
    bool isRegionBed() const { return isRegion; }
    bool isRegionLociFile() const { return geneloci.size() != 0; }
    int32_t getLenGeneLoci() const { return len_geneloci; }
    void isRegionOff() { isRegion=false; }
    bool ExistGeneLociFile(const std::string &genename) const {
      return geneloci.find(genename) != geneloci.end();
    }
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
    double scale_tag;
    double scale_ratio;
    double scale_pvalue;
    double alpha;

    DrawParam(): showymem(true), showylab(true) {}

    void setOpts(MyOpt::Opts &allopts, const CommandParamSet &cps);
    void setValues(const MyOpt::Variables &values, const int32_t n);
    void InitDump() const;
    bool isshowymem() const { return showymem; };
    bool isshowylab() const { return showylab; };

    double getlineheight() const { return ystep * barnum; }

    int32_t getNumLine(const int32_t s, const int32_t e) const{
      int32_t nline = (e-s -1) / width_per_line +1;
      return nline;
    }
    int32_t getNumPage(const int32_t s, const int32_t e) const {
      int32_t npage = (getNumLine(s,e) -1) / linenum_per_page +1;
      return npage;
    }

    int32_t getlpp() const { return linenum_per_page; }
    int32_t getHeightEachSample(const SamplePairEach &pair) const;
    int32_t getHeightAllSample(const Global &p, const std::vector<SamplePairOverlayed> &pairs) const;
    int32_t getPageHeight(const Global &p, const std::vector<SamplePairOverlayed> &pairs) const;
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
    Profile prof;

    std::vector<chrsize> gt;
    SampleInfoList vsinfo;
    std::vector<SamplePairOverlayed> samplepair;
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
    const std::string getPrefixName() const
    {
      return oprefix;
    }
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
