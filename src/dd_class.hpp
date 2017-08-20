/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DD_CLASS_H_
#define _DD_CLASS_H_

#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include "WigStats.hpp"
#include "SSP/common/BoostOptions.hpp"

class chrsize;

enum class DrompaCommand {CHIP, NORM, THRE, ANNO_PC, ANNO_GV, DRAW, REGION, SCALE, CG, PD, TR, PROF, OVERLAY, OTHER};

class SampleFile {
  double lambda;
  double nb_p, nb_n, nb_p0;
  WigType iftype;
  int32_t binsize;
 public:
  std::vector<int> data;

  void setbinsize(std::string &v, const int32_t b) {
    if(b>0) binsize = b;
    else {
      try {
	binsize = std::stoi(v);
      }catch (...) {
	binsize = 0;
      }
    }
    if(binsize <= 0) PRINTERR("invalid binsize: " << v);
  }
  
 SampleFile() {}
  SampleFile(const std::string &str, const int32_t b, const WigType &type): binsize(0) {
   std::vector<std::string> v;
   boost::split(v, str, boost::algorithm::is_any_of("."));
   int last(v.size()-1);

   if(type != WigType::NONE) {
     iftype = type;
     binsize = b;
   } else {
     if(v[last] == "wig") iftype = WigType::UNCOMPRESSWIG;    
     else if(v[last] == "gz" && v[last-1] == "wig") {
       iftype = WigType::COMPRESSWIG;
       --last;
     } else if(v[last] == "bedGraph") iftype = WigType::BEDGRAPH;
     else if(v[last] == "bw")         iftype = WigType::BIGWIG;
     else if(v[last] == "bin")        iftype = WigType::BINARY;
     else PRINTERR("invalid postfix: " << str);
   }
   setbinsize(v[last-1], b);
  }
 int getbinsize()    const { return binsize; }
 WigType getiftype() const { return iftype; }
};

class yScale {
  enum {TAG_DEFAULT=30, RATIO_DEFAULT=3, P_DEFAULT=5};
 public:
  double tag;
  double ratio;
  double pvalue;
 yScale(): tag(TAG_DEFAULT), ratio(RATIO_DEFAULT), pvalue(P_DEFAULT) {}
};

class SamplePair {
  int overlay;
  //  double fc; //comp

  std::string peak_argv;
  /*  Peak *peak;
  char *peak_argv;
  char *peakarray;
  char *linename;
  int *binnum;*/
  int binsize;
  
  yScale scale;
  
 public:
  std::string argvChIP, argvInput;
  std::string label;

  SamplePair(const std::vector<std::string> &v, const int32_t b): peak_argv(""), argvInput(""), label("") {
    if(v[0] != "") argvChIP  = v[0];
    if(v.size() >=2 && v[1] != "") argvInput = v[1];
    if(v.size() >=3 && v[2] != "") label     = v[2];
    if(v.size() >=4 && v[3] != "") peak_argv = v[3];
    binsize = b;
    if(v.size() >=6 && v[5] != "") scale.tag = stof(v[5]);
    if(v.size() >=7 && v[6] != "") scale.ratio  = stof(v[6]);
    if(v.size() >=8 && v[7] != "") scale.pvalue = stof(v[7]);
  }
  void print() {
    std::cout << boost::format("ChIP: %1% label: %2% peaklist: %3%\n") % argvChIP % label % peak_argv;
    std::cout << boost::format("   Input: %1%\n") % argvInput;
    std::cout << boost::format("   binsize: %1%\n") % binsize;
  }
  void printall() {
    std::cout << boost::format("ChIP:%1% Input:%2% name:%3% peak:%4% scale_tag %5% scale_ratio %6% scale_pvalue %7%\n")
      % argvChIP % argvInput % label % peak_argv % scale.tag % scale.ratio % scale.pvalue;
  }
};

class pdSample {
 public:
  std::string argv;
  std::string name;
  pdSample(){}
};


namespace DROMPA {

  class Scale {
    MyOpt::Opts opt;
    int32_t barnum;
    double ystep;
    
  public:
    double scale_tag;
    double scale_ratio;
    double scale_pvalue;
    Scale(): opt("Scale for Y axis",100)
    {
      opt.add_options()
	("scale_tag",    boost::program_options::value<double>()->default_value(30), "Scale for read line")
	("scale_ratio",  boost::program_options::value<double>()->default_value(5),  "Scale for fold enrichment")
	("scale_pvalue", boost::program_options::value<double>()->default_value(5),  "Scale for -log10(p)")
	("bn",           boost::program_options::value<int32_t>()->default_value(2), "Number of memories of y-axis")
	("ystep",        boost::program_options::value<double>()->default_value(15), "Height of read line")
	;
    }

    void setOpts(MyOpt::Opts &allopts) {
      allopts.add(opt);
    }
  
    void setValues(const MyOpt::Variables &values) {
      DEBUGprint("Scale setValues...");

      scale_tag    = MyOpt::getVal<double>(values, "scale_tag");
      scale_ratio  = MyOpt::getVal<double>(values, "scale_ratio");
      scale_pvalue = MyOpt::getVal<double>(values, "scale_pvalue");
      barnum = MyOpt::getVal<int32_t>(values, "bn");
      ystep  = MyOpt::getVal<double>(values, "ystep");

      DEBUGprint("Scale setValues done.");
    }
    double getlineheight() const { return ystep * barnum; }
  };
  
  class DrawRegion {
    MyOpt::Opts opt;

    std::string chr;
    std::string region;
    std::string generegion;
    int32_t len_genefile;

  public:
    DrawRegion():
      opt("Region to draw",100),
      chr(""), region(""),
      generegion(""), len_genefile(0)
    {
      opt.add_options()
	("chr",         boost::program_options::value<std::string>(), "Output the specified chromosome only")
	("region,r",    boost::program_options::value<std::string>(), "Specify genomic regions for drawing")
	("genefile",    boost::program_options::value<std::string>(), "Specify gene loci to visualize")  
	("len_genefile",boost::program_options::value<int32_t>()->default_value(50000), "extended length for each gene locus")
	;
    }

    void setOpts(MyOpt::Opts &allopts) {
      allopts.add(opt);
    }
  
    void setValues(const MyOpt::Variables &values) {
      DEBUGprint("DrawRegion setValues...");

      if(values.count("chr"))          chr          = MyOpt::getVal<std::string>(values, "chr");
      if(values.count("region"))       region       = MyOpt::getVal<std::string>(values, "region");
      if(values.count("genefile"))     generegion   = MyOpt::getVal<std::string>(values, "genefile");
      if(values.count("len_genefile")) len_genefile = MyOpt::getVal<int32_t>(values, "len_genefile");
      
      DEBUGprint("DrawRegion setValues done.");
    }

    std::string getchr() const { return chr; } 
  };
  
  class DrawParam {
    
    MyOpt::Opts opt;
    bool showars;
    int32_t linenum_per_page;
    bool showymem;
    bool showylab;
    int32_t viz;

    int32_t lineheight;
    int32_t samplenum;

  public:
    int32_t width_per_line;
    int32_t showctag;
    int32_t showitag;
    int32_t showratio;
    int32_t showpinter;
    int32_t showpenrich;
    
    DrawParam():
      opt("Drawing",100),
      showars(false), showymem(true), showylab(true),
      lineheight(0)
    {
      opt.add_options()
	("showctag",
	 boost::program_options::value<int32_t>()->default_value(1)->notifier(boost::bind(&MyOpt::range<int32_t>, _1, 0, 1, "--showctag")),
	 "Display ChIP read lines")
	("showitag",
	 boost::program_options::value<int32_t>()->default_value(0)->notifier(boost::bind(&MyOpt::range<int32_t>, _1, 0, 2, "--showitag")),
	 "Display Input read lines (0:off 1:all 2:first one)")
	("showratio",
	 boost::program_options::value<int32_t>()->default_value(0)->notifier(boost::bind(&MyOpt::range<int32_t>, _1, 0, 2, "--showratio")),
	 "Display ChIP/Input ratio (0:off 1:liner scale 2:logscale)")
	("showpinter",
	 boost::program_options::value<int32_t>()->default_value(0)->notifier(boost::bind(&MyOpt::range<int32_t>, _1, 0, 1, "--showpinter")),
	 "Display -log10(p) lines for ChIP internal")
	("showpenrich",
	 boost::program_options::value<int32_t>()->default_value(0)->notifier(boost::bind(&MyOpt::range<int32_t>, _1, 0, 1, "--showpenrich")),
	 "Display -log10(p) lines for ChIP/Input enrichment")
	("showars", "(For S.servisiae) Display ARS and do not display genes)")
	("ls",
	 boost::program_options::value<int32_t>()->default_value(1000)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--ls")),
	 "Width for each line (kp)")
	("lpp",
	 boost::program_options::value<int32_t>()->default_value(1)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--lpp")),
	 "Line number per page")
	("offymem", "Omit Y memory")
	("offylabel", "Omit Y label")
	("viz",
	 boost::program_options::value<int32_t>()->default_value(0)->notifier(boost::bind(&MyOpt::range<int32_t>, _1, 0, 2, "--viz")),
	 "Color of read profile\n     0: normal color\n     1: semitransparent color\n")
	;
    }

    void setOpts(MyOpt::Opts &allopts) {
      allopts.add(opt);
    }
  
    void setValues(const MyOpt::Variables &values, const int32_t l, const int32_t n) {
      DEBUGprint("DrawParam setValues...");

      showctag = MyOpt::getVal<int32_t>(values, "showctag");
      showitag = MyOpt::getVal<int32_t>(values, "showctag");
      showratio = MyOpt::getVal<int32_t>(values, "showctag");
      showpinter = MyOpt::getVal<int32_t>(values, "showctag");
      showpenrich = MyOpt::getVal<int32_t>(values, "showctag");
      showars  = values.count("showars");
      width_per_line = 1000 * MyOpt::getVal<int32_t>(values, "ls");
      linenum_per_page = MyOpt::getVal<int32_t>(values, "lpp");
      showymem = !values.count("offymem");
      showylab = !values.count("offylabel");
      viz      = MyOpt::getVal<int32_t>(values, "viz");

      lineheight = l;
      samplenum = n;
      
      DEBUGprint("DrawParam setValues done.");
    }

    bool isshowymem() const { return showymem; };
    bool isshowylab() const { return showylab; };
    int32_t getNumLine(const int32_t s, const int32_t e) const{
      return (e-s)/width_per_line +1;
    }
    int32_t getNumPage(const int32_t s, const int32_t e) const {
      return (getNumLine(s,e) -1) / linenum_per_page +1;
    }

    int32_t getlpp() const { return linenum_per_page; }
    int32_t getHeightEachSample() const;
    int32_t getHeightAllSample() const;
    int32_t getPageHeight() const;
  };
  
  class Global {
    bool ispng;
    WigType iftype;
    std::string oprefix;
    bool includeYM;

    // smoothing
  public:
    MyOpt::Opts opts;
    DrawParam drawparam;
    DrawRegion drawregion;
    Scale scale;

    std::vector<chrsize> gt;
    std::unordered_map<std::string, SampleFile> sample;
    std::vector<SamplePair> samplepair;
    std::vector<pdSample> pd;
    
    Global(): ispng(false), iftype(WigType::NONE), oprefix(""), includeYM(false), opts("Options") {}
    
    void setOpts(std::vector<DrompaCommand> &st);
    void setValues(const std::vector<DrompaCommand> &vopts, const MyOpt::Variables &values);

    WigType getIfType() const {return iftype;}

    const std::string getFigFileNameChr(const std::string &chr) const
    {
      return oprefix + "_chr" + chr + ".pdf";
    }
    bool isincludeYM() const { return includeYM; }
  };

}

#endif /* _DD_CLASS_H_ */
