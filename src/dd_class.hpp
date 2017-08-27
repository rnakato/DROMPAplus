/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
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

enum class DrompaCommand {CHIP, NORM, THRE, ANNO_PC, ANNO_GV, DRAW, REGION, CG, PD, TR, PROF, OVERLAY, OTHER};

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
	binsize = std::stoi(v);
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

class SamplePair {
  int32_t overlay;

  std::string peak_argv;
  std::unordered_map<std::string, std::vector<bed>> peaks;
  int32_t binsize;
  
  yScale scale;
  
 public:
  std::string argvChIP, argvInput;
  std::string label;

  SamplePair(const std::vector<std::string> &v, const int32_t b):
    peak_argv(""), argvInput(""), label("")
  {
    if(v[0] != "") argvChIP  = v[0];
    if(v.size() >=2 && v[1] != "") argvInput = v[1];
    if(v.size() >=3 && v[2] != "") label     = v[2];
    if(v.size() >=4 && v[3] != "") peak_argv = v[3];
    peaks = parseBed_Hash<bed>(peak_argv);
    binsize = b;
    if(v.size() >=6 && v[5] != "") scale.tag = stof(v[5]);
    if(v.size() >=7 && v[6] != "") scale.ratio  = stof(v[6]);
    if(v.size() >=8 && v[7] != "") scale.pvalue = stof(v[7]);

    printBed_Hash(peaks);
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
  const std::vector<bed> & getpeaksChr(const std::string &chrname) const { return peaks.at(rmchr(chrname)); }
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

    MyOpt::Opts optPC;
    MyOpt::Opts optGV;

  public:
    std::string genefile;
    int32_t gftype;
    //    bool showgene;
    HashOfGeneDataMap gmp;
    bool showtranscriptname;
    std::string arsfile;
    bool showars;
    std::string terfile;
    std::vector<std::vector <bed>> vbedlist;
    std::string repeatfile;
    std::string interfile;
    std::string mpfile;
    double mpthre;
    std::string gapfile;
    GraphFile GC;
    GraphFile GD;
    
    Annotation(): optPC("Annotation",100), optGV("Optional data",100),
		  genefile(""), arsfile(""), terfile(""), repeatfile(""),
		  interfile(""), mpfile(""), gapfile("")
    {
      optPC.add_options()
	("genefile,g", boost::program_options::value<std::string>(), "Gene annotation file")
	("gftype",
	 boost::program_options::value<int32_t>()->default_value(0)->notifier(boost::bind(&MyOpt::range<int32_t>, _1, 0, 2, "--gftype")),
	 "Format of gene annotation\n     0: RefFlat\n     1: GTF\n     2: SGD (for S. cerevisiae)\n")
	//	("showasgene", "Show one representative for each gene (default: all isoforms)")
	("showtranscriptname", "Show transcript name (default: gene name)")
	("ars",    boost::program_options::value<std::string>(), "ARS list (for yeast)")
	("ter",    boost::program_options::value<std::string>(), "TER list (for S.cerevisiae)")
	("showars", "Display ARS and TER and do not display genes")
	("bed",    boost::program_options::value<std::vector<std::string>>(), "<bedfile>,<label>: Specify bed file and name (<label> can be omited)")
	("repeat", boost::program_options::value<std::string>(), "Display repeat annotation (RepeatMasker format)")
	;
      optGV.add_options()
	("inter",  boost::program_options::value<std::vector<std::string>>(), "<interaction file>,<label>: Specify interaction file and name (<label> can be omited)")  // FDR de iro kaery
	("mp",     boost::program_options::value<std::string>(), "Mappability file")
	("mpthre", boost::program_options::value<double>()->default_value(0.3), "Low mappability threshold")
	("gap",    boost::program_options::value<std::string>(), "Specify gapped regions to be shaded")
	("GC",     boost::program_options::value<std::string>(), "Visualize GC contents graph")
	("gcsize",
	 boost::program_options::value<int32_t>()->default_value(100000)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--gcsize")),
	 "Window size for GC contents")
	("GD",     boost::program_options::value<std::string>(), "Visualize gene density (number of genes for each window)")
	("gdsize",
	 boost::program_options::value<int32_t>()->default_value(100000)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--gdsize")),
	 "Window size for gene density")
	;
    }

    void setOptsPC(MyOpt::Opts &allopts) { allopts.add(optPC); }
    void setOptsGV(MyOpt::Opts &allopts) { allopts.add(optGV); }

    HashOfGeneDataMap getGMP() {
      HashOfGeneDataMap tmp;
      if(!gftype)        tmp = parseRefFlat(genefile);
      else if(gftype==1) tmp = parseGtf(genefile);
      else if(gftype==2) tmp = parseSGD(genefile);
      else PRINTERR("invalid --gftype: " << gftype);
      
      //      printMap(tmp);
      
      //      if(showgene && gftype!=2) return construct_gmp(tmp); // hash for genes
	return tmp; // hash for transcripts
    }
    
    void setValuesPC(const MyOpt::Variables &values) {
      DEBUGprint("AnnoPC setValues...");

      try {
	for (auto x: {"genefile", "ars", "ter", "repeat"}) if (values.count(x)) isFile(MyOpt::getVal<std::string>(values, x));
	if (values.count("genefile")) {
	  genefile = MyOpt::getVal<std::string>(values, "genefile");
	  gftype   = MyOpt::getVal<int32_t>(values, "gftype");
	  //	  showgene = values.count("showasgene");
	  gmp = getGMP();
	}
	if (values.count("ars")) {
	  arsfile = MyOpt::getVal<std::string>(values, "ars");
	  parseARSOriDB(arsfile, gmp);
	}
	showars = values.count("showars");
	if (values.count("ter")) {
	  terfile = MyOpt::getVal<std::string>(values, "ter");
	  parseTER(terfile, gmp);
	}
	showtranscriptname = values.count("showtranscriptname");
	//	printMap(gmp);
	if (values.count("repeat")) repeatfile = MyOpt::getVal<std::string>(values, "repeat");

	if (values.count("bed")) { 
	  for(auto x: MyOpt::getVal<std::vector<std::string>>(values, "bed")) {
	    //    std::cout << x << std::endl;
	    auto vbed = parseBed<bed>(x);
	    //    printBed(vbed);
	    vbedlist.emplace_back(vbed);
	  }
	}
      } catch(const boost::bad_any_cast& e) {
	std::cout << e.what() << std::endl;
	exit(0);
      }

      DEBUGprint("AnnoPC setValues done.");
    }
    void setValuesGV(const MyOpt::Variables &values) {
      DEBUGprint("AnnoGV setValues...");
      
      try {
	if (values.count("inter")) {
	  interfile = MyOpt::getVal<std::string>(values, "inter");
	  isFile(interfile);
	}
	if (values.count("mp")) mpfile = MyOpt::getVal<std::string>(values, "mp");
	mpthre = MyOpt::getVal<double>(values, "mpthre");
	if (values.count("gap")) gapfile = MyOpt::getVal<std::string>(values, "gap");
	if (values.count("GC")) GC.setValue(MyOpt::getVal<std::string>(values, "GC"), MyOpt::getVal<int32_t>(values, "gcsize"));
	if (values.count("GD")) GD.setValue(MyOpt::getVal<std::string>(values, "GD"), MyOpt::getVal<int32_t>(values, "gdsize"));

      } catch(const boost::bad_any_cast& e) {
	std::cout << e.what() << std::endl;
	exit(0);
      }
      DEBUGprint("AnnoGV setValues done.");
    }

    void InitDumpPC(const MyOpt::Variables &values) const {
      std::vector<std::string> str_gftype = {"refFlat", "GTF", "SGD"};
      std::vector<std::string> str_gtype = {"Transcript", "Gene"};
      DEBUGprint("INITDUMP:DrompaCommand::ANNO_PC");
      std::cout << boost::format("\nAnnotations:\n");
      if(genefile != "") std::cout << boost::format("   Gene file: %1%, Format: %2%\n") % genefile % str_gftype[gftype]; //% str_gtype[showgene];
      //      if(arsfile != "")  std::cout << boost::format("   ARS file: %1%\n") % arsfile;
      //if(terfile != "")  std::cout << boost::format("   TER file: %1%\n") % terfile;
      //if(repeatfile != "") std::cout << boost::format("   Repeat file: %1%\n") % repeatfile;
      MyOpt::printOpt<std::string>(values, "ars",    "   ARS file");
      MyOpt::printOpt<std::string>(values, "ter",    "   TER file");
      if(showars) std::cout << "Display ARS and TER only." << std::endl;
      MyOpt::printOpt<std::string>(values, "repeat", "   Repeat file");
      MyOpt::printOpt<std::string>(values, "region", "   Region file");
      MyOpt::printVOpt<std::string>(values, "bed", "   Bed file");
    }
    void InitDumpGV(const MyOpt::Variables &values) const {
      DEBUGprint("INITDUMP:DrompaCommand::ANNO_GV");
      std::cout << boost::format("\nAnnotations:\n");
      MyOpt::printVOpt<std::string>(values, "inter", "   Interaction file");
      if (mpfile != "") {
	std::cout << boost::format("Mappability file directory: %1%\n") % mpfile;
	std::cout << boost::format("\tLow mappablitiy threshold: %1$2f\n") % mpthre;
      }
      MyOpt::printOpt<std::string>(values, "gc", "   GCcontents file");
      MyOpt::printOpt<std::string>(values, "gd", "   Gene density file");
      /*	if(d->GC.argv)     std::cout << boost::format("   GCcontents file: %1%\n")   % values["gc"].as<std::string>();
		if(d->GD.argv)     std::cout << boost::format("   Gene density file: %1%\n") % values["gd"].as<std::string>();*/
    }
    
    int32_t getgftype() const { return gftype; }
  };
  
  class Threshold {
    MyOpt::Opts opt;
    
  public:
    double pthre_inter;
    double pthre_enrich;
    double qthre;
    double ethre;
    double ipm;
    bool sigtest;
    int32_t width4lmd;
    
    Threshold(): opt("Threshold",100)
    {
      opt.add_options()
	("pthre_internal", boost::program_options::value<double>()->default_value(1e-4), "p-value for ChIP internal")
	("pthre_enrich",   boost::program_options::value<double>()->default_value(1e-4), "p-value for ChIP/Input enrichment")
	("qthre",          boost::program_options::value<double>()->default_value(1),    "FDR")
	("ethre",          boost::program_options::value<double>()->default_value(2),    "IP/Input fold enrichment")
	("ipm",            boost::program_options::value<double>()->default_value(0),    "Read intensity of peak summit")
	("nosig", "Omit highlighting peak regions")
	("width4lmd",
	 boost::program_options::value<int32_t>()->default_value(100000)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 0, "--width4lmd")),
	 "Width for calculating local lambda")
	;
    }

    void setOpts(MyOpt::Opts &allopts) {
      allopts.add(opt);
    }
    void setValues(const MyOpt::Variables &values) {
      DEBUGprint("Threshold setValues...");
      try {
	pthre_inter  = MyOpt::getVal<double>(values, "pthre_internal");
	pthre_enrich = MyOpt::getVal<double>(values, "pthre_enrich");
	qthre        = MyOpt::getVal<double>(values, "qthre");
	ethre        = MyOpt::getVal<double>(values, "ethre");
	ipm          = MyOpt::getVal<double>(values, "ipm");
	sigtest      = !values.count("nosig");
	width4lmd    = MyOpt::getVal<int32_t>(values, "width4lmd");
      } catch(const boost::bad_any_cast& e) {
	std::cout << e.what() << std::endl;
	exit(0);
      }
      DEBUGprint("Threshold setValues done.");
    }
    void InitDump() const {
      std::vector<std::string> str_bool = {"OFF", "ON"};
      DEBUGprint("INITDUMP:DrompaCommand::THRE");
      
      std::cout << boost::format("   significance test: %1%\n") % str_bool[sigtest];
      if (sigtest) {
	std::cout << boost::format("   p-value threshold (internal, -log10): %1$.2e\n")   % pthre_inter;
	std::cout << boost::format("   p-value threshold (internal, -log10): %1$.2e\n")   % pthre_inter;
	std::cout << boost::format("   p-value threshold (enrichment, -log10): %1$.2e\n") % pthre_enrich;
	std::cout << boost::format("   FDR threshold: %1$.2e\n")                          % qthre;
	std::cout << boost::format("   Peak intensity threshold: %1$.2f\n")               % ipm;
	std::cout << boost::format("   Enrichment threshold: %1$.2f\n")                   % ethre;
      }
    }
  };
  
  class DrawRegion {
    MyOpt::Opts opt;
    bool isRegion;
    std::vector<bed> regionBed;

    std::string chr;
    std::string genelocifile;
    int32_t len_geneloci;

  public:
    DrawRegion():
      opt("Region to draw",100),
      isRegion(false), chr(""), genelocifile(""), len_geneloci(0)
    {
      opt.add_options()
	("chr",
         boost::program_options::value<std::string>(),
	 "Output the specified chromosome only")
	("region,r",
	 boost::program_options::value<std::string>(),
	 "Specify genomic regions for drawing")
	("genelocifile",
	 boost::program_options::value<std::string>(),
	 "Specify of a file gene namaes to visualize")  
	("len_geneloci",
	 boost::program_options::value<int32_t>()->default_value(50000)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--len_geneloci")),
	 "extended length for each gene locus")
	;
    }

    void setOpts(MyOpt::Opts &allopts) {
      allopts.add(opt);
    }
    void setValues(const MyOpt::Variables &values) {
      DEBUGprint("DrawRegion setValues...");
      try {
	for (auto x: {"region", "genelocifile"}) if (values.count(x)) isFile(MyOpt::getVal<std::string>(values, x));
	if (values.count("chr")) chr = MyOpt::getVal<std::string>(values, "chr");
	if (values.count("region")) {
	  isRegion = true;
	  regionBed = parseBed<bed>(MyOpt::getVal<std::string>(values, "region"));
	  if(!regionBed.size()) PRINTERR("Error no bed regions in " << MyOpt::getVal<std::string>(values, "region"));
	  printBed(regionBed);
	}
	if (values.count("genelocifile")) genelocifile = MyOpt::getVal<std::string>(values, "genelocifile");
	len_geneloci = MyOpt::getVal<int32_t>(values, "len_geneloci");
      
      } catch(const boost::bad_any_cast& e) {
	std::cout << e.what() << std::endl;
	exit(0);
      }
      DEBUGprint("DrawRegion setValues done.");
    }
    void InitDump(const MyOpt::Variables &values) const {
      std::vector<std::string> str_bool = {"OFF", "ON"};
      
      DEBUGprint("INITDUMP:DrompaCommand::DRAWREGION");
      MyOpt::printOpt<std::string>(values, "region",    "   Region file");
      if (chr != "") std::cout << boost::format("   output chr%1% only.\n") % chr;
      if (genelocifile != "") {
	std::cout << boost::format("    Geneloci file: %1%, around %2% bp\n") % genelocifile % len_geneloci;
      }
    }
    
    std::vector<bed> getRegionBedChr(const std::string &chrname) {
      std::vector<bed> vbed;
      for(auto &x: regionBed) {
	if (x.chr == chrname || x.chr == "chr" + chrname) vbed.emplace_back(x);
      }
      return vbed;
    }
    const std::string & getchr() const { return chr; }
    bool isRegionBed() const { return isRegion; }
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

    void setOpts(MyOpt::Opts &allopts) {
      MyOpt::Opts opt("Drawing",100);
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
	("scale_tag",
	 boost::program_options::value<double>()->default_value(30)->notifier(boost::bind(&MyOpt::over<double>, _1, 0, "--scale_tag")),
	 "Scale for read line")
	("scale_ratio",
	 boost::program_options::value<double>()->default_value(5)->notifier(boost::bind(&MyOpt::over<double>, _1, 0, "--scale_ratio")),
	 "Scale for fold enrichment")
	("scale_pvalue",
	 boost::program_options::value<double>()->default_value(5)->notifier(boost::bind(&MyOpt::over<double>, _1, 0, "--scale_pvalue")),
	 "Scale for -log10(p)")
	("ls",
	 boost::program_options::value<int32_t>()->default_value(1000)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--ls")),
	 "Width for each line (kp)")
	("lpp",
	 boost::program_options::value<int32_t>()->default_value(1)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--lpp")),
	 "Line number per page")
	("bn",
	 boost::program_options::value<int32_t>()->default_value(2)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--bn")),
	 "Number of memories of y-axis")
	("ystep",
	 boost::program_options::value<double>()->default_value(15)->notifier(boost::bind(&MyOpt::over<double>, _1, 1, "--ystep")),
	 "Height of read line")
	("offymem", "Omit Y memory")
	("offylabel", "Omit Y label")
	("viz",
	 boost::program_options::value<int32_t>()->default_value(0)->notifier(boost::bind(&MyOpt::range<int32_t>, _1, 0, 2, "--viz")),
	 "Color of read profile\n     0: normal color\n     1: semitransparent color\n")
	;
      allopts.add(opt);
    }
  
    void setValues(const MyOpt::Variables &values, const int32_t n) {
      DEBUGprint("DrawParam setValues...");
      try {
	showctag = MyOpt::getVal<int32_t>(values, "showctag");
	showitag = MyOpt::getVal<int32_t>(values, "showitag");
	showratio = MyOpt::getVal<int32_t>(values, "showratio");
	showpinter = MyOpt::getVal<int32_t>(values, "showpinter");
	showpenrich = MyOpt::getVal<int32_t>(values, "showpenrich");
	width_per_line = 1000 * MyOpt::getVal<int32_t>(values, "ls");
	linenum_per_page = MyOpt::getVal<int32_t>(values, "lpp");
	barnum = MyOpt::getVal<int32_t>(values, "bn");
	ystep  = MyOpt::getVal<double>(values, "ystep");
	showymem = !values.count("offymem");
	showylab = !values.count("offylabel");
	viz      = MyOpt::getVal<int32_t>(values, "viz");

	scale_tag    = MyOpt::getVal<double>(values, "scale_tag");
	scale_ratio  = MyOpt::getVal<double>(values, "scale_ratio");
	scale_pvalue = MyOpt::getVal<double>(values, "scale_pvalue");
      
	samplenum = n;
      } catch(const boost::bad_any_cast& e) {
	std::cout << e.what() << std::endl;
	exit(0);
      }
      DEBUGprint("DrawParam setValues done.");
    }
    void InitDump() const {
      std::vector<std::string> str_bool = {"OFF", "ON"};
      std::vector<std::string> str_input = {"OFF", "ALL", "FIRST"};
      std::vector<std::string> str_ratio = {"OFF", "Linear", "Logratio"};
      
      DEBUGprint("INITDUMP:DrompaCommand::DRAW");
      std::cout << boost::format("\nFigure parameter:\n");
      std::cout << boost::format("   Display read: ChIP %1%, Input %2%, y-axis scale: %3%\n") % str_bool[showctag] % str_input[showitag] % scale_tag; 
      std::cout << boost::format("   Display enrichment: %1%, y-axis scale: %2%\n")           % str_ratio[showratio] % scale_ratio;
      std::cout << boost::format("   Display pvalue (internal): %1%, y-axis scale: %2%\n")    % str_bool[showpinter] % scale_pvalue;
      std::cout << boost::format("   Display pvalue (ChIP/Input): %1%, y-axis scale: %2%\n")  % str_bool[showpenrich] % scale_pvalue;
      std::cout << boost::format("   Width per line: %1% kbp\n")           % (width_per_line/1000);
      std::cout << boost::format("   Y-axis label: %1%\n")                 % str_bool[showylab];
      std::cout << boost::format("   Y-axis memory: %1%\n")                % str_bool[showymem];
    }
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
    int32_t getHeightEachSample(const SamplePairChr &pair) const;
    int32_t getHeightAllSample(const Global &p, const std::vector<SamplePairChr> &pairs) const;
    int32_t getPageHeight(const Global &p, const std::vector<SamplePairChr> &pairs) const;
  };
  
  class Global {
    bool ispng;
    bool rmchr;
    WigType iftype;
    std::string oprefix;
    bool includeYM;
    int32_t norm;
    int32_t sm;

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
    
    Global():
      ispng(false), rmchr(false), iftype(WigType::NONE),
      oprefix(""), includeYM(false),
      opts("Options")
    {}

    void setOpts(std::vector<DrompaCommand> &st);
    void setValues(const std::vector<DrompaCommand> &vopts, const MyOpt::Variables &values);
    
    void setOptsNorm(MyOpt::Opts &allopts) {
      MyOpt::Opts o("Normalization",100);
      o.add_options()
	("norm",
	 boost::program_options::value<int32_t>()->default_value(1)->notifier(boost::bind(&MyOpt::range<int32_t>, _1, 0, 2, "--norm")),
	 "Normalization between ChIP and Input\n      0: not normalize\n      1: with total read number (genome)\n      2: with total read number (each chr)\n      3: with NCIS method\n")
	("sm",
	 boost::program_options::value<int32_t>()->default_value(0)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 0, "--sm")),
	 "Smoothing width") // gausian ??
	;
      allopts.add(o);
    }
    void setValuesNorm(const MyOpt::Variables &values) {
      DEBUGprint("Norm setValues...");
      try {
	norm = MyOpt::getVal<int32_t>(values, "norm");
	sm   = MyOpt::getVal<int32_t>(values, "sm");
      } catch(const boost::bad_any_cast& e) {
	std::cout << e.what() << std::endl;
	exit(0);
      }
      DEBUGprint("Norm setValues done.");
    }
    void setOptsOther(MyOpt::Opts &allopts) {
      MyOpt::Opts o("Others",100);
      o.add_options()
	("includeYM", "output peaks of chromosome Y and M")
	("rmchr",   "Remove chromosome-separated pdf files")
	("png",     "Output with png format (Note: output each page separately)")
	("threads,p",
	 boost::program_options::value<int32_t>()->default_value(1)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--thread")),
	 "number of threads to launch")
	("help,h", "show help message")
	;
      allopts.add(o);
    }  
    void setValuesOther(const MyOpt::Variables &values) {
      DEBUGprint("Other setValues...");
      try {
	includeYM = values.count("includeYM");
	ispng = values.count("png");
	rmchr = values.count("rmchr");
      } catch(const boost::bad_any_cast& e) {
	std::cout << e.what() << std::endl;
	exit(0);
      }
      DEBUGprint("Other setValues done.");
    }
    void InitDumpNorm() const {
      std::vector<std::string> str_norm = { "OFF", "TOTALREAD GENOME", "TOTALREAD CHR", "NCIS" };
      
      DEBUGprint("INITDUMP:DrompaCommand::NORM");
      std::cout << boost::format("   ChIP/Input normalization: %1%\n") % str_norm[norm];
      if (sm) std::cout << boost::format("   smoothing width: %1% bp\n") % sm;
    }
    void InitDumpOther() const {
      std::vector<std::string> str_bool = {"OFF", "ON"};
      std::vector<std::string> str_format = {"PDF", "PNG"};
      
      DEBUGprint("INITDUMP:DrompaCommand::OTHER");
      std::cout << boost::format("   Output format: %1%\n") % str_format[ispng];
      if (includeYM) std::cout << boost::format("   include chromosome Y and M\n");
      if(rmchr) std::cout << boost::format("   remove chr pdfs\n");
    }

    WigType getIfType() const { return iftype; }

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
    bool isrmchr() const { return rmchr; }
  };

}

#endif /* _DD_CLASS_H_ */
