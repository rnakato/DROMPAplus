/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <unordered_map>
#include <ext/stdio_filebuf.h>
#include "dd_readfile.hpp"

/* 1:ChIP   2:Input   3:name   4:peaklist   5:binsize
   6:scale_tag   7:scale_ratio   8:scale_pvalue */
//void scan_samplestr(const std::string &str, std::unordered_map<std::string, SampleFile> &sample, std::vector<SamplePair> &samplepair)
void scan_samplestr(const std::string &str, const std::vector<chrsize> gt,
		    std::unordered_map<std::string, SampleFile> &sample,
		    std::vector<SamplePair> &samplepair,
		    WigType iftype)
{
  std::vector<std::string> v;
  boost::split(v, str, boost::algorithm::is_any_of(","));
  int32_t binsize(0);

  if(v.size() >8) {
    std::cerr << "error: sample std::string has ',' more than 8: " << str << std::endl;
    exit(1);
  }
  if(v[0] == "") {
      std::cerr << "please specify ChIP sample: " << str << std::endl;
      exit(1);
  }
  
  if(v.size() >4) binsize = stoi(v[4]);
  
  if(sample.find(v[0]) == sample.end()) sample[v[0]] = SampleFile(v[0], gt, binsize, iftype);
  if(sample[v[0]].getbinsize() <= 0) PRINTERR("please specify binsize.\n");
  
  if(v.size() >=2 && v[1] != "") {
    if(sample.find(v[1]) == sample.end()) sample[v[1]] = SampleFile(v[1], gt, binsize, iftype);
    if(sample[v[0]].getbinsize() != sample[v[1]].getbinsize()) PRINTERR("binsize of ChIP and Input should be same. " << str);
  }

  samplepair.emplace_back(v, sample[v[0]].getbinsize());
  
  return;
}

pdSample scan_pdstr(const std::string &str)
{
  std::vector<std::string> v;
  boost::split(v, str, boost::algorithm::is_any_of(","));

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
    boost::split(v, lineStr, boost::algorithm::is_any_of(" \t"), boost::algorithm::token_compress_on);
    array.setval((stoi(v[0])-1)/binsize, stol(v[1]));
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
    boost::split(v, lineStr, boost::algorithm::is_any_of(" \t"), boost::algorithm::token_compress_on);
    if (v[0] != chrname) {
      if (!on) continue;
      else break;
    }
    if (!on) on=1;
    
    int32_t start(stoi(v[1]));
    if (start%binsize) PRINTERR("ERROR: invalid start position: " << start << " for binsize " << binsize);
    array.setval(start/binsize, stol(v[3]));
  }
}
    
void readBinary(WigArray &array, const std::string &filename, const int32_t nbin)
{
  static int nbinsum(0);
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (!in) PRINTERR("cannot open " << filename);

  in.seekg(nbinsum * sizeof(int32_t));

  array.readBinary(in, nbin);
  // array.printArray();

  nbinsum += nbin;
  return;
}

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
  char tmpfile[] = "/tmp/DROMPAplus_bigWigToBedGraph_XXXXXX";
  if ((fd = mkstemp(tmpfile)) < 0){
    perror("mkstemp");
    //      fd = EXIT_FAILURE;
  }
  std::string command = "bigWigToBedGraph -chrom=chr" + rmchr(chrname) + " " + filename + " " + std::string(tmpfile);
  int32_t return_code = system(command.c_str());
  if(WEXITSTATUS(return_code)) {
    std::cerr << "Warning: command " << command << "return nonzero status." << std::endl;
  }
  readBedGraph(array, std::string(tmpfile), chrname, binsize);
  unlink(tmpfile);
}
void funcBedGraph(WigArray &array, const std::string &filename, const int32_t binsize, const std::string &chrname)
{
  DEBUGprint("WigType::BEDGRAPH");
  readBedGraph(array, filename, chrname, binsize);
}

void funcBinary(WigArray &array, const std::string &filename, const int32_t nbin)
{
  DEBUGprint("WigType::BINARY");
  readBinary(array, filename, nbin);
}


WigArray loadWigData(const std::string &filename, const SampleFile &x, const chrsize &chr)
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
  else if (iftype == WigType::BINARY)        funcBinary(array, filename, nbin);
  
  /*  if (iftype == WigType::NONE) {
    if (isStr(filename, ".bin"))           funcBinary(array, filename, nbin);
    else if (isStr(filename, ".bedGraph")) funcBedGraph(array, filename, binsize, chrname);
    else if (isStr(filename, ".bw"))       funcBigWig(array, filename, binsize, chrname);
    else if (isStr(filename, ".wig.gz"))   funcCompressWig(array, filename, binsize, chrname);
    else if (isStr(filename, ".wig"))      funcWig(array, filename, binsize, chrname);
    else PRINTERR("Suffix error of "<< filename <<". please specify --iftype option.");
  } else if (iftype == WigType::UNCOMPRESSWIG) funcWig(array, filename, binsize, chrname);
  else if (iftype == WigType::COMPRESSWIG)     funcCompressWig(array, filename, binsize, chrname);
  else if (iftype == WigType::BIGWIG)          funcBigWig(array, filename, binsize, chrname);
  else if (iftype == WigType::BEDGRAPH)        funcBedGraph(array, filename, binsize, chrname);
  else if (iftype == WigType::BINARY)          funcBinary(array, filename, nbin);
  */
  //array.dump();

  // NCIS hakokode

  
  
  return array;
}


int32_t getNcolReadNum(std::string &lineStr)
{
  int32_t ncol_readnum(0);
  std::vector<std::string> v;
  boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
  for (size_t i=0; i<v.size(); ++i) {
    if(isStr(v[i], "normalized read number")) {
      ncol_readnum = i;
      //      std::cout << ncol_readnum << std::endl;
      break;
    }
  }
  return ncol_readnum;
}

void SampleFile::scanStatsFile(const std::string &filename)
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
	boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
	totalreadnum = stoi(v[ncol_readnum]);
	on=1;
      }
    } else {
      std::vector<std::string> v;
      boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
      totalreadnum_chr[v[0]] = stoi(v[ncol_readnum]);
    }
  }
}

void SampleFile::gettotalreadnum(const std::string &filename, const std::vector<chrsize> gt)
{
  std::string statsfile(prefix + "csv");
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
