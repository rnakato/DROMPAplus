/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <ext/stdio_filebuf.h>
#include "dd_readfile.hpp"

pdSample scan_pdstr(const std::string &str)
{
  std::vector<std::string> v;
  ParseLine(v, str, ',');
  //  boost::split(v, str, boost::algorithm::is_any_of(","));

  if(v.size() > 2) {
    std::cerr << "error: sample std::string has ',' more than 2: " << str << std::endl;
    exit(1);
  }

  pdSample pd;
  if(v[0] == "") {
      std::cerr << "please specify file: " << str << std::endl;
      exit(1);
  } else {
    pd.argv = v[0];
  }
  if(v[1] != "") pd.name = v[1];
  else pd.name = v[0];

  return pd;
}

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

void readBedGraph(WigArray &array, const std::string &filename, const std::string &chrname, const int32_t binsize)
{
  std::ifstream in(filename);
  if (!in) PRINTERR("cannot open " << filename);
    
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
    
    int32_t start(stoi(v[1]));
    int32_t end(stoi(v[2])-1);
    if (start%binsize) PRINTERR("ERROR: invalid start position: " << start << " for binsize " << binsize);
    for(int32_t i=start; i<=end; ++i) array.setval(i/binsize, stol(v[3]));
  }
}

/*void readBinary(WigArray &array, const std::string &filename, const int32_t nbin)
{
  static int nbinsum(0);
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (!in) PRINTERR("cannot open " << filename);

  in.seekg(nbinsum * sizeof(int32_t));

  array.readBinary(in, nbin);
  // array.printArray();

  nbinsum += nbin;
  return;
  }*/

void funcWig(WigArray &array, const std::string &filename, const int32_t binsize, const std::string &chrname)
{
  DEBUGprint("WigType::UNCOMPRESSWIG");
  std::ifstream in(filename);
  if (!in) PRINTERR("cannot open " << filename);
  readWig(in, array, chrname, binsize);
}

void funcCompressWig(WigArray &array, const std::string &filename, const int32_t binsize, const std::string &chrname)
{
  DEBUGprint("WigType::COMPRESSWIG");
  std::string command = "zcat " + filename;
  FILE *fp = popen(command.c_str(), "r");
  __gnu_cxx::stdio_filebuf<char> *p_fb = new __gnu_cxx::stdio_filebuf<char>(fp, std::ios_base::in);
  std::istream in(static_cast<std::streambuf *>(p_fb));
  readWig(in, array, chrname, binsize);
}

void funcBigWig(WigArray &array, const std::string &filename, const int32_t binsize, const std::string &chrname)
{
  DEBUGprint("WigType::BIGWIG");
  int32_t fd(0);
  char tmpfile[] = "/tmp/DROMPAplus_bigWigToBedGraph_XXXXXXXX";
  if ((fd = mkstemp(tmpfile)) < 0){
    perror("mkstemp");
    //      fd = EXIT_FAILURE;
  }
  std::string command = "bigWigToBedGraph -chrom=chr" + rmchr(chrname) + " " + filename + " " + std::string(tmpfile);
  int32_t return_code = system(command.c_str());
  if(WEXITSTATUS(return_code)) {
    std::cerr << "Error: command " << command << "return nonzero status." << std::endl;
    exit(0);
  }
  readBedGraph(array, std::string(tmpfile), chrname, binsize);
  unlink(tmpfile);
}
void funcBedGraph(WigArray &array, const std::string &filename, const int32_t binsize, const std::string &chrname)
{
  DEBUGprint("WigType::BEDGRAPH");
  readBedGraph(array, filename, chrname, binsize);
}

/*void funcBinary(WigArray &array, const std::string &filename, const int32_t nbin)
{
  DEBUGprint("WigType::BINARY");
  readBinary(array, filename, nbin);
  }*/


WigArray loadWigData(const std::string &filename, const SampleInfo &x, const chrsize &chr)
{
  //  std::cout << chr.getname() << std::endl;

  int32_t binsize(x.getbinsize());
  int32_t nbin(chr.getlen()/binsize +1);

  WigArray array(nbin, 0);
  //  std::string filename(x.first);
  std::string chrname(chr.getrefname());
  WigType iftype(x.getiftype());

  if (iftype == WigType::NONE) PRINTERR("Suffix error of "<< filename <<". please specify --iftype option.");
  else if (iftype == WigType::UNCOMPRESSWIG) funcWig(array, filename, binsize, chrname);
  else if (iftype == WigType::COMPRESSWIG)   funcCompressWig(array, filename, binsize, chrname);
  else if (iftype == WigType::BIGWIG)        funcBigWig(array, filename, binsize, chrname);
  else if (iftype == WigType::BEDGRAPH)      funcBedGraph(array, filename, binsize, chrname);
  //  else if (iftype == WigType::BINARY)        funcBinary(array, filename, nbin);
  
  //array.dump();

  return array;
}

int32_t getNcolReadNum(std::string &lineStr)
{
  int32_t ncol_readnum(0);
  std::vector<std::string> v;
  ParseLine(v, lineStr, '\t');
  //  boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
  for (size_t i=0; i<v.size(); ++i) {
    if(isStr(v[i], "normalized read number")) {
      ncol_readnum = i;
      //      std::cout << ncol_readnum << std::endl;
      break;
    }
  }
  return ncol_readnum;
}

void SampleInfo::scanStatsFile(const std::string &filename)
{
  DEBUGprint("scanStatsFile...");
   
  std::ifstream in(filename);
  if (!in) PRINTERR("cannot open " << filename);
    
  std::string lineStr;
  int32_t on(0);
  int32_t ncol_readnum(0);
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || isStr(lineStr, "% genome")) continue;

    if(!on){
      if(isStr(lineStr, "normalized read number")) ncol_readnum = getNcolReadNum(lineStr);
      else if(isStr(lineStr, "Genome")) {
	std::vector<std::string> v;
	ParseLine(v, lineStr, '\t');
	//	boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
	totalreadnum = stoi(v[ncol_readnum]);
	on=1;
      }
    } else {
      std::vector<std::string> v;
      ParseLine(v, lineStr, '\t');
      //      boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
      totalreadnum_chr[v[0]] = stoi(v[ncol_readnum]);
    }
  }
}

void SampleInfo::gettotalreadnum(const std::string &filename, const std::vector<chrsize> gt)
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
  }
  
#ifdef DEBUG
    std::cout << "Total read number:" << std::endl;
    std::cout << "Whole genome: " << totalreadnum << std::endl;
    for(auto &chr: gt) {
      std::cout << chr.getname() << ": " << totalreadnum_chr[chr.getname()] << std::endl;
    }
#endif
}
