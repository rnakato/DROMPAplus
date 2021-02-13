/* Copyright(c)  Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "../submodules/SSP/common/gzstream.h"
#include "dd_readfile.hpp"

namespace {
  void SplitBedGraphLine(std::vector<std::string> &v, const std::string &str)
  {
    size_t current(0), found;
    while((found = str.find_first_of(" \t", current)) != std::string::npos) {
      v.emplace_back(std::string(str, current, found - current));
      current = found + 1;
    }
    v.emplace_back(std::string(str, current, str.size() - current));
    return;
  }

  template <class T>
  void readWig(T &in, WigArray &array, const std::string &chrname, const int binsize)
  {
    std::string head("chrom="+ chrname);
    int32_t on(0);

    std::string lineStr;
    while (!in.eof()) {
      getline(in, lineStr);
      if(lineStr.empty() || !lineStr.find("track")) continue;
      if(on && isStr(lineStr, "chrom=")) break;
      if(isStr(lineStr, head)) {
        on=1;
        continue;
      }
      if(!on) continue;
      std::vector<std::string> v;
      SplitBedGraphLine(v, lineStr);
      array.setval((stoi(v[0])-1)/binsize, stod(v[1]));
    }
    return;
  }

  void readBedGraph(WigArray &array, const std::string &filename,
                    const std::string &chrname, const int32_t binsize)
  {
    std::ifstream in(filename);
    if (!in) PRINTERR_AND_EXIT("cannot open " << filename);

    DEBUGprint_FUNCStart();

    int32_t on(0);
    std::string lineStr;
    while (!in.eof()) {
      getline(in, lineStr);
      if (lineStr.empty()) continue;
      std::vector<std::string> v;
      SplitBedGraphLine(v, lineStr);
      if (v[0] != chrname) {
        if (!on) continue;
        else break;
      }
      if (!on) on=1;

      if(v.size() < 4) {
        std::cerr << "\nError: invalid delimitar in BedGraph file?: " << lineStr << std::endl;
        exit(1);
      }

      double val(0);
      if(v[3] == "") val = 0; else val = stod(v[3]);
      //std::cout << chrname << "\t" << v[0] << "\t" << binsize << "\t" << val << "\t" << v[2] << "\t" << std::endl;

      try {
        int32_t start(stoi(v[1]));
        int32_t end(stoi(v[2])-1);
        if (start%binsize) PRINTERR_AND_EXIT("ERROR: invalid start position: " << start << " for binsize " << binsize);
        int32_t s(start/binsize);
        int32_t e(end/binsize);
        //    std::cout << s << "\t " << e << "\t " << array.size() << "\t" << stod(v[3]) << std::endl;
        for(int32_t i=s; i<=e; ++i) array.setval(i, val);
      } catch (const boost::bad_any_cast& e) {
        PRINTERR_AND_EXIT("Error: invalid value in BedGraph. " + lineStr + ": :" + std::string(e.what()));
      }

    }

    in.close();

    DEBUGprint_FUNCend();
    return;
  }

  void funcWig(WigArray &array, const std::string &filename,
               const int32_t binsize, const std::string &chrname)
  {
    DEBUGprint_FUNCStart();

    std::ifstream in(filename);
    if (!in) PRINTERR_AND_EXIT("cannot open " << filename);
    readWig(in, array, chrname, binsize);
    in.close();

    DEBUGprint_FUNCend();
  }

  void funcCompressWig(WigArray &array, const std::string &filename,
                       const int32_t binsize, const std::string &chrname)
  {
    DEBUGprint_FUNCStart();

    igzstream in(filename.c_str());
    readWig(in, array, chrname, binsize);

    DEBUGprint_FUNCend();
  }

  void funcBigWig(WigArray &array, const std::string &filename,
                  const int32_t binsize, const std::string &chrname)
  {
    DEBUGprint_FUNCStart();

    int32_t fd(0);
    char tmpfile[] = "/tmp/DROMPAplus_bigWigToBedGraph_XXXXXXXX";
    if ((fd = mkstemp(tmpfile)) < 0) {
      perror("mkstemp");
      //      fd = EXIT_FAILURE;
    }
    std::string command = "bigWigToBedGraph -chrom=chr" + rmchr(chrname) + " " + filename + " " + std::string(tmpfile);
    int32_t return_code = system(command.c_str());
    if (WEXITSTATUS(return_code)) {
      PRINTERR_AND_EXIT("Error: command " << command
                        << "return nonzero status. "
                        << "Add the PATH to 'DROMPAplus/otherbins'.");
    }
    readBedGraph(array, std::string(tmpfile), chrname, binsize);
    unlink(tmpfile);
    close(fd);

    DEBUGprint_FUNCend();
  }

  void funcBedGraph(WigArray &array, const std::string &filename,
                    const int32_t binsize, const std::string &chrname)
  {
    DEBUGprint_FUNCStart();

    readBedGraph(array, filename, chrname, binsize);

    DEBUGprint_FUNCend();
  }
}

WigArray loadWigData(const std::string &filename, const SampleInfo &x, const chrsize &chr)
{
  int32_t binsize(x.getbinsize());
  int32_t nbin(chr.getlen()/binsize +1);

  WigArray array(nbin, 0);
  std::string chrname(chr.getrefname());
  WigType iftype(x.getiftype());

  if (iftype == WigType::NONE) PRINTERR_AND_EXIT("Suffix error of "<< filename <<". please specify --iftype option.");
  else if (iftype == WigType::UNCOMPRESSWIG) funcWig(array, filename, binsize, chrname);
  else if (iftype == WigType::COMPRESSWIG)   funcCompressWig(array, filename, binsize, chrname);
  else if (iftype == WigType::BIGWIG)        funcBigWig(array, filename, binsize, chrname);
  else if (iftype == WigType::BEDGRAPH)      funcBedGraph(array, filename, binsize, chrname);

  //array.dump();

  return array;
}


vChrArray::vChrArray(const DROMPA::Global &p, const chrsize &_chr):
  chr(_chr)
{
  std::cout << "Load sample data..";
  for (auto &x: p.vsinfo.getarray()) {
    clock_t t1,t2;
    t1 = clock();
    arrays[x.first] = ChrArray(p, x, chr);
    t2 = clock();
    PrintTime(t1, t2, "ChrArray new");
  }

#ifdef DEBUG
  std::cout << "all WigArray:" << std::endl;
  for (auto &x: arrays) {
    std::cout << x.first << ", binsize " << x.second.binsize << std::endl;
  }
  std::cout << "all SamplePair:" << std::endl;
  for (auto &x: p.samplepair) {
    std::cout << x.first.argvChIP << "," << x.first.argvInput
              << ", binsize " << x.first.getbinsize() << std::endl;
  }
#endif
}
