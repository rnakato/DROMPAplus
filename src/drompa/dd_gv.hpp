/* Copyright(c) Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_GV_H_
#define _DD_GV_H_

#include <boost/algorithm/string.hpp>
#include "dd_sample_definition.hpp"
#include "ReadAnnotation.hpp"
#include "../submodules/SSP/common/util.hpp"
#include "../submodules/SSP/common/BoostOptions.hpp"

class chrsize;

enum class DrompaCommand {
  CHIP, NORM, THRE, ANNO_PC, ANNO_GV,
  DRAW, REGION, GENWIG, PROF, MULTICI,
  CG, TR, OTHER
};

class CommandParamSet {
public:
  int32_t sm;
  int32_t showctag;
  int32_t showratio;
  double scaletag;
  double scaleratio;
  double scalepvalue;
  double thre_pinter;
  double thre_penrich;
  double thre_ethre;
  double thre_ipm;

  CommandParamSet(const int32_t defsm,
                  const int32_t defshowctag,
                  const int32_t defshowratio,
                  const double defscaletag,
                  const double defscaleratio,
                  const double defscalepvalue,
                  const double defthre_pinter,
                  const double defthre_penrich,
                  const double defthre_ethre,
                  const double defthre_ipm):
    sm(defsm),
    showctag(defshowctag),
    showratio(defshowratio),
    scaletag(defscaletag),
    scaleratio(defscaleratio),
    scalepvalue(defscalepvalue),
    thre_pinter(defthre_pinter),
    thre_penrich(defthre_penrich),
    thre_ethre(defthre_ethre),
    thre_ipm(defthre_ipm)
  {}
};

namespace DROMPA {
  class Global;

  /// GraphDataFileName (used in Annotation)
  class GraphDataFileName {
    std::string filename;
    int32_t binsize;

  public:
    GraphDataFileName(): filename(""), binsize(0) {}

    void setValue(const std::string &f, const int32_t s) {
      filename = f;
      binsize = s;
    }
    const std::string & getfilename() const { return filename; }
    int32_t getbinsize() const { return binsize; }
    bool isOn() const {
      if (filename != "") return true;
      else return false;
    }
  };


  //// Annotation
  class Annotation {
    bool isUCSC;
    bool isIdeogram;
    bool isChIADrop;

    template <class T>
    std::vector<vbed<T>> readBedFile(const std::vector<std::string> &vfilename) {
      std::vector<vbed<T>> list;

      for (auto &x: vfilename) {
        std::vector<std::string> v;
        ParseLine(v, x, ',');
        auto vbed = parseBed<T>(v[0]);
        //    printBed(vbed);
        if (v.size()>1) list.emplace_back(vbed, v[1]);
        else            list.emplace_back(vbed, "Bed");
      }
      return list;
    }

    void parse_ChIADropData(const std::string &fileName);

  public:
    std::string genefile;
    int32_t gftype;
    HashOfGeneDataMap gmp;
    HashOfGeneDataMap arsgmp;
    bool showtranscriptname;
    std::string arsfile;
    std::string terfile;
    std::vector<vbed<bed>> vbedlist;
    std::vector<vbed<bed12>> vbed12list;
    std::vector<InteractionSet> vinterlist;
    std::vector<cytoband> vcytoband;
    std::unordered_map<std::string,
                       std::unordered_map<std::string, std::vector<int32_t>>> mp_ChIADrop;
    int32_t chia_distance_thre;
    std::string repeatfile;
    std::string mpfile;
    double mpthre;
    std::string gapfile;
    GraphDataFileName GC;
    GraphDataFileName GD;

    Annotation():
      isUCSC(false), //showars(false),
      isIdeogram(false),
      isChIADrop(false),
      genefile(""), gftype(0),
      showtranscriptname(false),
      arsfile(""),
      terfile(""),
      chia_distance_thre(100000),
      repeatfile(""),
      mpfile(""), mpthre(0),
      gapfile("")
    {}

    void setOptsPC(MyOpt::Opts &allopts);
    void setOptsGV(MyOpt::Opts &allopts);

    HashOfGeneDataMap getGMP() {
      HashOfGeneDataMap tmp;
      if (!gftype)        tmp = parseRefFlat(genefile);
      else if (gftype==1) tmp = parseGtf(genefile);
      else if (gftype==2) tmp = parseSGD(genefile);
      else PRINTERR_AND_EXIT("invalid --gftype: " << gftype);

      isUCSC = isGeneUCSC(tmp);
      //      printMap(tmp);

      return tmp; // hash for transcripts
    }

    void setValuesPC(const MyOpt::Variables &values);
    void setValuesGV(const MyOpt::Variables &values);
    void InitDumpPC(const MyOpt::Variables &values) const;
    void InitDumpGV(const MyOpt::Variables &values) const;

    int32_t getgftype() const { return gftype; }
    bool is_Anno_UCSC() const {return isUCSC; }
    //    bool isshowars() const { return showars; }
    bool existChIADrop() const { return isChIADrop; }
    bool showIdeogram() const { return isIdeogram; }
  };


  ////// Profile
  class Profile {
    enum {TSS, TTS, GENE100, BEDSITES, PTYPENUM};
    int32_t width_from_center;
    int32_t ptype;
    int32_t stype;
    bool usefixedlength;
    //   int32_t ntype;

  public:

    Profile(): width_from_center(0),
               ptype(0), stype(0),
               usefixedlength(false)
    {} //, ntype(0)

    void setOpts(MyOpt::Opts &allopts);
    void setValues(const MyOpt::Variables &values);
    void InitDump() const;
    bool isPtypeTSS()     const { return ptype == TSS; }
    bool isPtypeTTS()     const { return ptype == TTS; }
    bool isPtypeGene100() const { return ptype == GENE100; }
    bool isPtypeBed()     const { return ptype == BEDSITES; }
    bool isusefixedlength() const { return usefixedlength; }
    int32_t get_width_from_center() const { return width_from_center; }
    int32_t is_distribution_type() const { return stype; }
//    int32_t is_normalization_type() const { return ntype; }
  };


  //// Threshold
  class Threshold {
  public:
    double pthre_inter;
    double pthre_enrich;
    double ethre;
    double ipm;
    bool sigtest;

    Threshold(): pthre_inter(0), pthre_enrich(0),
                 ethre(0), ipm(0), sigtest(false) {}

    void setOpts(MyOpt::Opts &allopts, const CommandParamSet &cps);
    void setValues(const MyOpt::Variables &values);
    void InitDump() const;
  };


  //// DrawRegion
  class DrawRegion {
    bool isRegion;
    std::vector<bed> regionBed;

    std::string chr;
    std::unordered_map<std::string, int32_t> geneloci;
    int32_t len_geneloci;

    void getGeneLoci(const std::string &genelocifile) {
      std::ifstream in(genelocifile);
      if (!in) PRINTERR_AND_EXIT("cannot open " << genelocifile);

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
      for (auto &x: regionBed) {
        if (x.chr == chrname || x.chr == "chr" + chrname) vbed.emplace_back(x);
      }
      return vbed;
    }
    const std::string & getchr() const { return chr; }
    bool isRegionBed() const { return isRegion; }
    bool isGeneLociFile() const { return geneloci.size() != 0; }
    int32_t getLenGeneLoci() const { return len_geneloci; }
    void isRegionOff() { isRegion=false; }
    bool ExistInGeneLociFile(const std::string &genename) const {
      return geneloci.find(genename) != geneloci.end();
    }
  };


  /// Drawparam
  class DrawParam {
    int32_t linenum_per_page;
    int32_t barnum;
    double ystep;
    bool shownegative;
    bool showymem;
    bool showylab;
    bool showbedname;
    bool showpdf;

    int32_t samplenum;

  public:
    int32_t width_page_pixel;
    int32_t width_draw_pixel;
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

    DrawParam(): linenum_per_page(0), barnum(0), ystep(0),
                 showymem(true), showylab(true), showbedname(true), showpdf(true),
                 samplenum(0),
                 width_page_pixel(0), width_draw_pixel(0), width_per_line(0),
                 showctag(0), showitag(0), showratio(0), showpinter(0), showpenrich(0),
                 scale_tag(0), scale_ratio(0), scale_pvalue(0),
                 alpha(0) {}

    void setOpts(MyOpt::Opts &allopts, const CommandParamSet &cps);
    void setValues(const MyOpt::Variables &values, const int32_t n);
    void InitDump() const;
    bool isshownegative() const { return shownegative; };
    bool isshowymem() const { return showymem; };
    bool isshowylab() const { return showylab; };
    bool isshowbedname() const { return showbedname; };
    bool isshowpdf() const { return showpdf; };

    int32_t getbarnum() const { return barnum; }
    double getystep() const { return ystep; }
    double getHeightDf() const { return ystep * barnum; }

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


  /// Global
  class Global {
    bool ispng;
    bool showchr;
    WigType iftype;
    std::string oprefix;
    bool includeYM;
    int32_t norm;
    int32_t smoothing;

    WigType genwig_oftype;
    int32_t genwig_ofvalue;
    std::string genometablefilename;

    // MULTICI
    bool getmaxval;
    bool addname;

  public:
    MyOpt::Opts opts;
    DrawParam drawparam;
    DrawRegion drawregion;
    Threshold thre;
    Annotation anno;
    Profile prof;

    std::vector<chrsize> gt;
    vSampleInfo vsinfo;
    std::vector<SamplePairOverlayed> samplepair;

    bool isGV;

    Global():
      ispng(false), showchr(false), iftype(WigType::NONE),
      oprefix(""), includeYM(false), norm(0), smoothing(0),
      genwig_ofvalue(0), getmaxval(false), addname(false),
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

    bool ischr_in_gt(const std::string &query) const {
      for(auto &chr: gt) {
        if (chr.getname() == query) return true;
      }
      return false;
    }

    int32_t getSmoothing() const { return smoothing; }
    int32_t getChIPInputNormType() const { return norm; }
    const std::string getPrefixName() const { return oprefix; }
    const std::string getFigFileName() const { return oprefix + ".pdf"; }
    const std::string getGenomeTableFileName() const { return genometablefilename; }
    const std::string getFigFileNameChr(const std::string &chr) const
    {
      return oprefix + "_" + chr + ".pdf";
    }
    const std::string genwig_getOutputFileTypeStr() const {
      std::vector<std::string> strType = {"COMPRESSED WIG", "WIG", "BEDGRAPH", "BIGWIG"};
      return strType[static_cast<int32_t>(genwig_oftype)];
    }
    bool isincludeYM() const { return includeYM; }
    bool isshowchr() const { return showchr; }
    bool isgetmaxval() const { return getmaxval; }
    bool isaddname() const { return addname; }

    void genwig_openfilestream() {
      for (auto &x: samplepair) x.first.genwig_openfilestream(getPrefixName(), genwig_oftype, genwig_ofvalue);
    }
    void genwig_closefilestream() {
      for (auto &x: samplepair) x.first.genwig_closefilestream(getGenomeTableFileName());
    }
  };
}

#endif /* _DD_GV_H_ */
