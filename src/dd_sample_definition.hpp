/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_SAMPLE_DEFINITION_H_
#define _DD_SAMPLE_DEFINITION_H_

#include <unordered_map>
#include "WigStats.hpp"
#include "extendBedFormat.hpp"

class vChrArray;

class SampleInfo {
  WigType iftype;
  int32_t binsize;
  int32_t totalreadnum;
  std::unordered_map<std::string, int32_t> totalreadnum_chr;
  std::string prefix;

  void setbinsize(std::string &v, const int32_t b);

public:
  std::vector<int> data;

  SampleInfo() {}
  SampleInfo(const std::string &filename,
	     const std::vector<chrsize> &gt,
	     const int32_t b,
	     const WigType &type);

  void scanStatsFile(const std::string &filename);
  void gettotalreadnum(const std::string &filename, const std::vector<chrsize> &gt);
  int32_t getbinsize() const { return binsize; }
  WigType getiftype() const { return iftype; }

  int32_t gettotalreadnum() const { return totalreadnum; }
  const std::unordered_map<std::string, int32_t>& gettotalreadnum_chr() const & {
    return totalreadnum_chr;
  }
};

class vSampleInfo {
  std::unordered_map<std::string, SampleInfo> vsinfo;

public:
  vSampleInfo(){}

  void addSampleInfo(const std::string &str, const std::vector<chrsize> &gt, const WigType iftype);
  bool Exists(const std::string &str) const {
    return vsinfo.find(str) != vsinfo.end();
  }
  int32_t getbinsize(const std::string &str) const {
    return vsinfo.at(str).getbinsize();
  }
  const std::unordered_map<std::string, SampleInfo> &getarray() const & {
    return vsinfo;
  }
  const SampleInfo& operator[](const std::string &str) const& {
    return vsinfo.at(str);
  }
};

class SamplePairEach {
  FILE* File;
  std::string genwig_filename;
  std::string genwig_tmpfile;

  int32_t binsize;

  std::unordered_map<std::string, std::vector<bed>> vbedregions;
  std::unordered_map<std::string, std::vector<Peak>> vPeak;

  class yScale {
  public:
    double tag;
    double ratio;
    double pvalue;
    yScale(): tag(0), ratio(0), pvalue(0) {}
  };

public:
  std::string argvChIP, argvInput;
  std::string peak_argv;
  std::string label;
  double ratio;
  yScale scale;

  SamplePairEach():
    genwig_filename(""), binsize(0), argvChIP(""), argvInput(""), peak_argv(""), label(""), ratio(1)
  {}
  SamplePairEach(const std::string &str, const vSampleInfo &vsinfo);

  void setScalingFactor(const int32_t normtype, const vChrArray &vReadArray, const std::string &chrname);

  void genEnrichWig(const vChrArray &vReadArray, const std::string &chrname);

  void peakcall_withInput(const vChrArray &vReadArray, const std::string &chrname,
			  const double pthre_inter, const double pthre_enrich,
			  const double ethre, const double ipm);
  void peakcall_onlyChIP(const vChrArray &vReadArray, const std::string &chrname, const double pthre_inter);

  void printPeak(const std::string &prefix) const {
    std::string filename = prefix + "." + label + ".peak.tsv";

    DEBUGprint("Peakfile name: " << filename);
    std::ofstream out(filename);

    Peak v;
    v.printHead(out);
    int32_t num(0);
    for (auto &x: vPeak) {
      for (auto &peak: x.second) peak.print(out, num++);
    }
  }

  std::vector<bed> getBedChr(const std::string &chrname) const {
    if (vbedregions.find(rmchr(chrname)) != vbedregions.end()) return vbedregions.at(rmchr(chrname));
    else return std::vector<bed>();
  }
  std::vector<Peak> getPeakChr(const std::string &chrname) const {
    if (vPeak.find(chrname) != vPeak.end()) return vPeak.at(chrname);
    else return std::vector<Peak>();
  }
  void print() const;
  int32_t getbinsize() const { return binsize; }
  bool BedExists() const { return peak_argv != ""; }
  bool InputExists() const { return argvInput != ""; }

  void genwig_openfilestream(const std::string &prefix, WigType oftype) {
    genwig_filename = prefix + "." + label + ".enrich";
    if (oftype==WigType::COMPRESSWIG || oftype==WigType::UNCOMPRESSWIG) {
      genwig_filename += ".wig";
      File = fopen(genwig_filename.c_str(), "w");
      fprintf(File, "track type=wiggle_0\tname=\"%s\"\tdescription=\"Merged tag counts for every %d bp\"\n",
	      genwig_filename.c_str(), binsize);

    } else if (oftype==WigType::BEDGRAPH) {
      genwig_filename += ".bedGraph";
      File = fopen(genwig_filename.c_str(), "w");
    } else if (oftype==WigType::BIGWIG) {
      genwig_filename += ".bw";
      int32_t fd(0);
      char tmpfile[] = "/tmp/drompa+_bedGraph_XXXXXX";
      if ((fd = mkstemp(tmpfile)) < 0){
	perror("mkstemp");
      }
      File = fopen(tmpfile, "w");
      genwig_tmpfile = std::string(tmpfile);
    }
  }

  void genwig_closefilestream(WigType oftype) {
    if (oftype==WigType::COMPRESSWIG || oftype==WigType::UNCOMPRESSWIG) {
      fclose(File);
      if (oftype==WigType::COMPRESSWIG) {
	std::string command = "gzip -f " + genwig_filename;
	if (system(command.c_str())) PRINTERR_AND_EXIT("gzip .wig failed.");
      }
    } else if (oftype==WigType::BEDGRAPH) {
      fclose(File);
    } else if (oftype==WigType::BIGWIG) {
      std::string command = "bedGraphToBigWig " + genwig_tmpfile + " " + genwig_filename;
      if (system(command.c_str())) {
	unlink(genwig_tmpfile.c_str());
	std::cerr << "Error: command " << command << "return nonzero status. "
		  << "Add the PATH to 'DROMPAplus/otherbins'." << std::endl;
      }
      unlink(genwig_tmpfile.c_str());
    }
  }

};

class SamplePairOverlayed {
  bool overlay;

 public:
  SamplePairEach first;
  SamplePairEach second;

  SamplePairOverlayed(const std::string &str, const vSampleInfo &vsinfo):
    overlay(false), first(str, vsinfo)
  {}

  void setSecondSample(const std::string &str, const vSampleInfo &vsinfo) {
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

#endif // _DD_SAMPLE_DEFINITION_H_
