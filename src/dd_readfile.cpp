/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <ext/stdio_filebuf.h>
#include "dd_readfile.hpp"
#include "version.hpp"

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
      //    std::cout << chrname << "\t" << v[0] << "\t" << binsize << "\t" << stod(v[3]) << "\t" << v[2] << "\t" << std::endl;

      int32_t start(stoi(v[1]));
      int32_t end(stoi(v[2])-1);
      if (start%binsize) PRINTERR_AND_EXIT("ERROR: invalid start position: " << start << " for binsize " << binsize);
      int32_t s(start/binsize);
      int32_t e(end/binsize);
      //    std::cout << s << "\t " << e << "\t " << array.size() << "\t" << stod(v[3]) << std::endl;
      for(int32_t i=s; i<=e; ++i) array.setval(i, stod(v[3]));
    }

    in.close();

    DEBUGprint_FUNCend();
    return;
  }

  /*void readBinary(WigArray &array, const std::string &filename, const int32_t nbin)
    {
    static int nbinsum(0);
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    if (!in) PRINTERR_AND_EXIT("cannot open " << filename);

    in.seekg(nbinsum * sizeof(int32_t));

    array.readBinary(in, nbin);
    // array.printArray();

    nbinsum += nbin;
    return;
    }*/

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

    std::string command = "zcat " + filename;
    FILE *fp = popen(command.c_str(), "r");
    __gnu_cxx::stdio_filebuf<char> *p_fb = new __gnu_cxx::stdio_filebuf<char>(fp, std::ios_base::in);
    std::istream in(static_cast<std::streambuf *>(p_fb));
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
      std::cerr << "Error: command " << command << "return nonzero status. "
		<< "Add the PATH to 'DROMPAplus/otherbins'." << std::endl;
      exit(0);
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

  /*void funcBinary(WigArray &array, const std::string &filename, const int32_t nbin)
    {
    DEBUGprint("WigType::BINARY");
    readBinary(array, filename, nbin);
    }*/
  int32_t getNcolReadNum(std::string &lineStr)
  {
    int32_t ncol_readnum(0);
    std::vector<std::string> v;
    ParseLine(v, lineStr, '\t');
    for (size_t i=0; i<v.size(); ++i) {
      if(isStr(v[i], "normalized read number")) {
	ncol_readnum = i;
	//      std::cout << ncol_readnum << std::endl;
	break;
      }
    }
    return ncol_readnum;
  }

  void OutputStatsfileForOtherData(const std::string &filename, const std::string &statsfile,
				   const std::vector<chrsize> &gt,
				   const std::unordered_map<std::string, int32_t> &totalreadnum_chr,
				   int32_t totalreadnum)
  {
    std::cout << statsfile << " not found. Generating...";
    std::ofstream out(statsfile);

    out << "Generated by drompa+ version " << VERSION << std::endl;
    out << "Input file: \"" << filename << "\"" << std::endl;

    // Global stats
    out << "\ttotal reads\tread depth" << std::endl;

    int32_t flen(150); // Predefined fragment length
    uint64_t lengenome(0);
    for(auto &x: gt) lengenome += x.getlen();

    out << "Genome\t"
	<< totalreadnum << "\t"
	<< getratio(totalreadnum * flen, lengenome)
	<< std::endl;
    for(auto &chr: gt) {
      out << chr.getrefname() << "\t"
	  << totalreadnum_chr.at(chr.getname()) << "\t"
	  << getratio(totalreadnum_chr.at(chr.getname()) * flen, chr.getlen())
	  << std::endl;
    }

    std::cout << "done." << std::endl;
  }
}

WigArray loadWigData(const std::string &filename, const SampleInfo &x, const chrsize &chr)
{
  //  std::cout << chr.getname() << std::endl;

  int32_t binsize(x.getbinsize());
  int32_t nbin(chr.getlen()/binsize +1);

  WigArray array(nbin, 0);
  //  std::string filename(x.first);
  std::string chrname(chr.getrefname());
  WigType iftype(x.getiftype());

  if (iftype == WigType::NONE) PRINTERR_AND_EXIT("Suffix error of "<< filename <<". please specify --iftype option.");
  else if (iftype == WigType::UNCOMPRESSWIG) funcWig(array, filename, binsize, chrname);
  else if (iftype == WigType::COMPRESSWIG)   funcCompressWig(array, filename, binsize, chrname);
  else if (iftype == WigType::BIGWIG)        funcBigWig(array, filename, binsize, chrname);
  else if (iftype == WigType::BEDGRAPH)      funcBedGraph(array, filename, binsize, chrname);
  //  else if (iftype == WigType::BINARY)        funcBinary(array, filename, nbin);

  //array.dump();

  return array;
}

void SampleInfo::scanStatsFile(const std::string &filename)
{
  DEBUGprint_FUNCStart();

  enum {NONE, PARSE2WIG, OTHER, OTHER_CHR};
  std::string lineStr;
  int32_t status(NONE);
  int32_t ncol_readnum(0);

  std::ifstream in(filename);
  if (!in) PRINTERR_AND_EXIT("cannot open " << filename);
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || isStr(lineStr, "% genome")) continue;

    std::vector<std::string> v;
    switch(status) {
    case NONE:
      if(isStr(lineStr, "normalized read number")) ncol_readnum = getNcolReadNum(lineStr);
      else if (isStr(lineStr, "Generated by drompa+")) {
	status = OTHER;
      } else if (isStr(lineStr, "Genome")) {
	std::vector<std::string> v;
	ParseLine(v, lineStr, '\t');
	totalreadnum = stoi(v[ncol_readnum]);
	status = PARSE2WIG;
      }
      break;
    case PARSE2WIG:
      ParseLine(v, lineStr, '\t');
      totalreadnum_chr[v[0]] = stoi(v[ncol_readnum]);
      break;
    case OTHER:
      if (isStr(lineStr, "Genome")) {
	ParseLine(v, lineStr, '\t');
	totalreadnum = stoi(v[1]);
	status = OTHER_CHR;
      }
      break;
    case OTHER_CHR:
      ParseLine(v, lineStr, '\t');
      totalreadnum_chr[rmchr(v[0])] = stoi(v[1]);
      break;
    }
  }

  DEBUGprint_FUNCend();
}

void SampleInfo::gettotalreadnum(const std::string &filename, const std::vector<chrsize> &gt)
{
  std::string statsfile(prefix + "tsv");
  if (checkFile(statsfile)) scanStatsFile(statsfile);
  else {
    DEBUGprint("loadWigData: noStatsFile...");
    for(auto &chr: gt) {
      WigArray array(loadWigData(filename, *this, chr));
      totalreadnum_chr[chr.getname()] = array.getArraySum();
      totalreadnum += totalreadnum_chr[chr.getname()];
    }
    OutputStatsfileForOtherData(filename, statsfile, gt, totalreadnum_chr, totalreadnum);
  }

#ifdef DEBUG
    std::cout << "Total read number:" << std::endl;
    std::cout << "Whole genome: " << totalreadnum << std::endl;
    for(auto &chr: gt) {
      std::cout << chr.getname() << ": " << totalreadnum_chr[chr.getname()] << std::endl;
    }
#endif
}

void SamplePairEach::setScalingFactor(const int32_t normtype, const vChrArray &vReadArray, const std::string &chrname)
{
  DEBUGprint_FUNCStart();
  if (argvInput == "") return; // ratio = 1;

  switch (normtype) {
  case 0:  // not normalize
    ratio = 1;
    break;
  case 1:  // total read for genome
    ratio = getratio(vReadArray.getArray(argvChIP).totalreadnum,
		     vReadArray.getArray(argvInput).totalreadnum);
    break;
  case 2:  // total read for each chromosome
    ratio = getratio(vReadArray.getArray(argvChIP).totalreadnum_chr.at(chrname),
		     vReadArray.getArray(argvInput).totalreadnum_chr.at(chrname));
    break;
  case 3:  // NCIS
    ratio = 1;
    break;
  }
#ifdef DEBUG
  std::cout << "ChIP/Input Ratio for chr " << chrname << ": " << ratio << std::endl;
#endif

  DEBUGprint_FUNCend();
}
