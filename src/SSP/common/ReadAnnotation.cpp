/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <boost/algorithm/string.hpp>
#include "ReadAnnotation.hpp"
#include "inline.hpp"
#include "util.hpp"

int32_t countmp(HashOfGeneDataMap &mp)
{
  int32_t n(0);
  for(auto pair: mp) {
    //  for(auto itr = mp.begin(); itr != mp.end(); ++itr) {
    //    for(auto itr2 = mp.at(itr->first).begin(); itr2 != mp.at(itr->first).end(); ++itr2) n++;
    for(auto x: mp.at(pair.first)) n++;
  }
  return n;
}

std::vector<std::string> scanGeneName(const HashOfGeneDataMap &mp)
{
  std::vector<std::string> vgname;
  for(auto pair: mp) {
    for(auto x: mp.at(pair.first)) vgname.push_back(x.first);
    //  for(auto itr = mp.begin(); itr != mp.end(); ++itr) {
    //    for(auto itr2 = mp.at(itr->first).begin(); itr2 != mp.at(itr->first).end(); ++itr2) {
      //    }
  }
  return vgname;
}

HashOfGeneDataMap extract_mp(const HashOfGeneDataMap &tmp, const std::vector<std::string> glist)
{
  HashOfGeneDataMap mp;

  for(auto &x: glist) {
    for(auto pair: tmp) {
    //    for(auto itr = tmp.begin(); itr != tmp.end(); ++itr) {
      std::string chr(pair.first);
      if (tmp.at(chr).find(x) != tmp.at(chr).end()) {
	mp[chr][x] = tmp.at(chr).at(x);
	break;
      }
    }
  }
  return mp;
}

std::vector<std::string> readGeneList(const std::string& fileName)
{
  std::ifstream in(fileName);
  if(!in) PRINTERR("genelist file does not exist.");

  std::vector<std::string> glist;
  std::string lineStr;
  
  while (!in.eof()) {
    getline(in, lineStr);
    if(!lineStr.empty()) glist.push_back(lineStr);
  }

  return glist;
}

HashOfGeneDataMap parseRefFlat(const std::string& fileName)
{
  if(isStr(fileName, ".gtf")) {
    std::cerr << "Warning: gene file seems to be gtf format but is parsed as refFlat." << std::endl;
  }

  std::ifstream in(fileName);
  if(!in) PRINTERR("refFlat file does not exist.");

  HashOfGeneDataMap tmp;
  std::string lineStr;
  
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty()) continue;

    std::vector<std::string> v, exonStarts, exonEnds;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));

    std::string tname(v[1]);
    std::string chr = rmchr(v[2]);
    
    tmp[chr][tname].tname   = tname;
    tmp[chr][tname].gname   = v[0];
    tmp[chr][tname].chr     = chr;
    tmp[chr][tname].strand  = v[3];
    tmp[chr][tname].txStart = stoi(v[4]);
    tmp[chr][tname].txEnd   = stoi(v[5]);
    tmp[chr][tname].cdsStart = stoi(v[6]);
    tmp[chr][tname].cdsEnd   = stoi(v[7]);
    tmp[chr][tname].exonCount = stoi(v[8]);

    boost::split(exonStarts, v[9], boost::algorithm::is_any_of(","));
    boost::split(exonEnds,  v[10], boost::algorithm::is_any_of(","));

    for(int32_t i=0; i<tmp[chr][tname].exonCount; i++) {
      tmp[chr][tname].exon.emplace_back(stoi(exonStarts[i]), stoi(exonEnds[i]));
    }
  }
  return tmp;
}

HashOfGeneDataMap parseGtf(const std::string& fileName)
{
  if(!isStr(fileName, ".gtf")) {
    std::cerr << "Warning: gene file may not be gtf format but is parsed as gtf." << std::endl;
  }

  std::ifstream in(fileName);
  if(!in) PRINTERR("gtf file does not exist.");

  HashOfGeneDataMap tmp;
  std::string lineStr;
  
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || lineStr[0] == '#') continue;
    std::vector<std::string> v;
    
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    std::string feat = v[2];
    if(feat == "gene" || feat == "transcript" || feat == "three_prime_utr" || feat == "five_prime_utr") continue;
    
    std::string chr = rmchr(v[0]);
    int32_t start = stoi(v[3]);
    int32_t end   = stoi(v[4]);
    std::string strand = v[6];
    std::string annotation = v[8];

    std::string gname, tname, gid, tid, gsrc, gtype, tsrc, ttype, ttag="";
    std::vector<std::string> idtab, vc;
    boost::split(idtab, annotation, boost::algorithm::is_any_of(";"));
    for (auto term: idtab) {
      boost::split(vc, term, boost::algorithm::is_any_of("\""));
      if(isStr(term, "gene_source"))             gsrc  = vc[1];
      else if(isStr(term, "gene_biotype"))       gtype = vc[1];
      else if(isStr(term, "transcript_source"))  tsrc  = vc[1];
      else if(isStr(term, "transcript_biotype")) ttype = vc[1];
      else if(isStr(term, "tag")) {
	if(ttag=="" || vc[1]=="CCDS")           ttag = vc[1];
	else if(vc[1]=="basic" && ttag!="CCDS") ttag = vc[1];
	else if( ttag!="basic" && ttag!="CCDS") ttag = vc[1];
      }
      else if(isStr(term, "transcript_name")) tname = vc[1];
      else if(isStr(term, "gene_name"))       gname = vc[1];
      else if(isStr(term, "transcript_id"))   tid = vc[1];
      else if(isStr(term, "gene_id"))         gid = vc[1];
    }

    tmp[chr][tname].tname  = tname;
    tmp[chr][tname].gname  = gname;
    tmp[chr][tname].tid  = tid;
    tmp[chr][tname].gid  = gid;
    tmp[chr][tname].chr  = chr;
    tmp[chr][tname].strand = strand;
    if(feat == "start_codon") {
      if(strand == "+") tmp[chr][tname].cdsStart = start;
      else tmp[chr][tname].cdsEnd = end;
    } else if(feat == "stop_codon") {
      if(strand == "+") tmp[chr][tname].cdsEnd = end;
      else tmp[chr][tname].cdsStart = start;
    } else if(feat == "exon") {
      tmp[chr][tname].exonCount++;
      if(!tmp[chr][tname].txStart || start < tmp[chr][tname].txStart) tmp[chr][tname].txStart = start;
      if(end > tmp[chr][tname].txEnd) tmp[chr][tname].txEnd = end;
      tmp[chr][tname].exon.emplace_back(start, end);
    }

    tmp[chr][tname].gsrc = gsrc;
    tmp[chr][tname].tsrc = tsrc;
    tmp[chr][tname].gtype = gtype;
    tmp[chr][tname].ttype = ttype;
    tmp[chr][tname].ttag = ttag;
  }

  return tmp;
}

HashOfGeneDataMap construct_gmp(const HashOfGeneDataMap &tmp)
{
  HashOfGeneDataMap gmp;

  for(auto pair: tmp) {
    std::string chr(pair.first);
    for(auto x: pair.second) {
      std::string gname(x.second.gname);
      if(gmp.find(chr) == gmp.end() || gmp[chr].find(gname) == gmp[chr].end()) gmp[chr][gname] = x.second;
      else if(x.second.ttag == "CCDS") {
	if(gmp[chr][gname].ttag != "CCDS" || gmp[chr][gname].length() < x.second.length()) gmp[chr][gname] = x.second;
      }else if(x.second.ttag == "basic") {
	if(gmp[chr][gname].ttag != "CCDS" &&
	   ((gmp[chr][gname].ttag == "basic" && gmp[chr][gname].length() < x.second.length()) ||
	    gmp[chr][gname].ttag != "basic")) gmp[chr][gname] = x.second;
      }
      else if(gmp[chr][gname].ttag != "basic" && (x.second.ttag == "basic" || gmp[chr][gname].length() < x.second.length())) gmp[chr][gname] = x.second;
      else if(x.second.ttag == "basic" || gmp[chr][gname].length() < x.second.length()) gmp[chr][gname] = x.second;
    }
  }
  return gmp;
}

void printMap(const HashOfGeneDataMap &mp)
{
  for(auto pair: mp) {
    for(auto x: pair.second) {
      x.second.printall();
      std::cout << std::endl;
    }
  }
  return;
}

void printRefFlat(const HashOfGeneDataMap &mp, const int32_t nameflag)
{
  for(auto pair: mp) {
    for(auto x: pair.second) {
      if(nameflag) std::cout << x.second.gname << "\t" << x.second.tname << "\t";
      else std::cout << x.second.gid << "\t" << x.second.tid << "\t";
      std::cout << pair.first << "\t"
	   << x.second.strand << "\t"
	   << x.second.txStart << "\t"
	   << x.second.txEnd << "\t";
      if(x.second.cdsStart) {
	std::cout << x.second.cdsStart << "\t"
	     << x.second.cdsEnd   << "\t";
      } else {
	std::cout << x.second.txEnd << "\t"
	     << x.second.txEnd << "\t";
      }
      std::cout << x.second.exonCount << "\t";
      for (auto &ex: x.second.exon) std::cout << ex.start << ",";
      std::cout << "\t";
      for (auto &ex: x.second.exon) std::cout << ex.end   << ",";
      std::cout << std::endl;
    }
  }
  return;
}

std::vector<chrsize> read_genometable(const std::string& fileName)
{
  std::ifstream in(fileName);
  if(!in) PRINTERR("genometable file does not exist.");

  std::vector<chrsize> gt;
  std::string lineStr;
  
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || lineStr[0] == '#') continue;
    std::vector<std::string> v;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    gt.emplace_back(rmchr(v[0]), stoi(v[1]));
  }
  return gt;
}

