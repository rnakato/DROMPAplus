/* Copyright(c)  Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <boost/filesystem.hpp>
#include "ReadMpbldata.hpp"
#include "../submodules/SSP/common/seq.hpp"
#include "../submodules/SSP/src/SeqStats.hpp"
#include "../submodules/SSP/common/gzstream.h"

std::vector<int32_t> readMpblWigArray(const std::string &mpfile,
				      const std::string &chrname,
				      const int32_t binsize,
				      const int32_t nbin)
{
  DEBUGprint_FUNCStart();
  std::string filename = mpfile + "/map_" + chrname + "." + std::to_string(binsize) + ".wig.gz";
  std::vector<int32_t> mparray(nbin, 0);

  DEBUGprint("mpfile: " << filename);

  igzstream in(filename.c_str());

  std::string lineStr;
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty()) continue;
    std::vector<std::string> v;
    ParseLine(v, lineStr, '\t');
    mparray[(stoi(v[0]))/binsize] = stod(v[1]);
  }

  DEBUGprint_FUNCend();
  return mparray;
}

namespace {
  void generateMpblWigData(const std::string &filename, std::vector<BpStatus> &mparray, const int32_t binsize)
  {
    std::cout << filename << ".gz not found. Generating.." << std::endl;

    FILE* File = fopen(filename.c_str(), "w");
    int32_t nbin(mparray.size()/binsize +1);
    std::vector<int32_t> wigarray(nbin, 0);

    for (size_t i=0; i<mparray.size(); ++i) {
      if (mparray[i] == BpStatus::MAPPABLE) ++wigarray[i/binsize];
    }
    for (int32_t i=0; i<nbin; ++i) {
      fprintf(File, "%d\t%.4f\n", i*binsize, wigarray[i]/(double)binsize);
    }
    fclose(File);

    // compression
    std::string command = "gzip -f " + filename;
    if (system(command.c_str())) PRINTERR_AND_EXIT("gzip .wig failed.");
  }
}


std::vector<BpStatus> readMpblBpArray(const std::string &mpfile,
				      const std::string &chrname,
				      const int32_t chrlen,
				      const int32_t binsize)
{
  static int32_t on(0);

  if(mpfile == "") {
    if(!on) {
      std::cout << "Mappability file is not specified. All genomeic regions are considered as mappable." << std::endl;
      on=1;
    }
    return std::vector<BpStatus>(chrlen, BpStatus::MAPPABLE);
  }

  if(!on) {
    std::cout << "Reading binary mappability file.." << std::flush;
    on=1;
  }
  std::vector<BpStatus> mparray(chrlen, BpStatus::UNMAPPABLE);

  std::string filename = mpfile + "/map_" + chrname + "_binary.txt.gz";
  isFile(filename);

  igzstream in(filename.c_str());

  int32_t n(0);
  int8_t c;
  while (!in.eof()) {
    c = in.get();
    if(c==' ') continue;
    if(c=='1') mparray[n] = BpStatus::MAPPABLE;
    ++n;
    if(n >= chrlen-1) break;
  }

  std::string mpblwigfile = mpfile + "/map_" + chrname + "." + std::to_string(binsize) + ".wig";
  boost::filesystem::path const file(mpblwigfile + ".gz");
  if(!boost::filesystem::exists(file)) {
    generateMpblWigData(mpblwigfile, mparray, binsize);
  }

  return mparray;
}

void setPeak_to_MpblBpArray(std::vector<BpStatus> &array,
			    const std::string &chrname,
			    const std::vector<bed> &vbed)
{
  int32_t chrlen(array.size());

  for(auto &bed: vbed) {
//    std::cout << bed.chr << "?" << chrname << std::endl;
    if(bed.chr == chrname) {
      size_t s(std::max(0, bed.start));
      size_t e(std::min(bed.end, chrlen-1));
//      printf("s=%d, e=%d\n", s,e);
      for(size_t i=s; i<=e; ++i) array[i] = BpStatus::INBED;
    }
  }

/*  int32_t n(0);
  for(size_t i=0; i<chrlen; ++i) {
    if(array[i]==BpStatus::INBED) ++n;
  }
  printf("test m %d\n", n);*/

 return;
}
