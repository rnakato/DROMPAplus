/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "Mapfile.hpp"

void FragmentLengthDist::outputDistFile(const std::string &prefix, const uint64_t nread) {
  std::string outputfile = prefix + ".ReadLengthDist.csv";
  std::ofstream out(outputfile);
  printVector(out, vlenF3, "F3", nread);
  if(pairedend) {
    out << std::endl;
    printVector(out, vlenF5, "F5", nread);
  }
  out.close();
  
  if(pairedend) {
    outputfile = prefix + ".FragmentLengthDist.csv";
    std::ofstream out(outputfile);
    printVector(out, vflen, "Fragmemt", nread);
  }
}

void SeqStatsGenome::setValues(const MyOpt::Variables &values) {
  DEBUGprint("SeqStatsGenome setValues...");
  
  inputfilename = MyOpt::getVal<std::string>(values, "input");
  pairedend     = values.count("pair");
  maxins        = MyOpt::getVal<int32_t>(values, "maxins");
  specifyFtype  = values.count("ftype");
  if(onFtype()) {
    ftype = MyOpt::getVal<std::string>(values, "ftype");
    if(ftype != "SAM" && ftype != "BAM" && ftype != "BOWTIE" && ftype != "TAGALIGN") PRINTERR("invalid --ftype.\n");
  }
  
  dflen.setValues(values);
  genometable = MyOpt::getVal<std::string>(values, "gt");
  readGenomeTable(genometable);
  
  if(values.count("mptable")) {
    for(auto &x: chr) x.getMptable(MyOpt::getVal<std::string>(values, "mptable"));
  }
  
  // Greekchr
  for(auto &x: chr) {
    if(x.getname() == "I") {
      for(auto &x:chr) x.Greekchron();
      break;
    }
  }
  
  // sepchr
  vsepchr = MyMthread::getVsepchr(getlen(), chr, MyOpt::getVal<int32_t>(values, "threads"));
  
  DEBUGprint("SeqStatsGenome setValues done.");
#ifdef DEBUG
  std::cout << "chr\tautosome" << std::endl;
  for(auto &x: chr) printList(x.getname(), x.isautosome());
  for(uint32_t i=0; i<vsepchr.size(); i++)
    std::cout << "thread " << (i+1) << ": " << vsepchr[i].s << "-" << vsepchr[i].e << std::endl;
  printReadstats();
#endif
}


void SeqStatsGenome::readGenomeTable(const std::string &gt)
{
  std::vector<std::string> v;
  std::string lineStr;
  std::ifstream in(gt);
  if(!in) PRINTERR("Could nome open " << gt << ".");
  
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || lineStr[0] == '#') continue;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    chr.emplace_back(v[0], stoi(v[1]));
  }
  return;
}

int32_t setIdLongestChr(SeqStatsGenome &genome)
{
  int32_t id(0);
  uint64_t lenmax(0);
  for(size_t i=0; i<genome.chr.size(); ++i) {
    if(lenmax < genome.chr[i].getlenmpbl()) {
      lenmax = genome.chr[i].getlenmpbl();
      id = i;
    }
  }
  return id;
}
