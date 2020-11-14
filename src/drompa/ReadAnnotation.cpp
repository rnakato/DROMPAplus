/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <boost/algorithm/string.hpp>
#include "ReadAnnotation.hpp"
#include "../submodules/SSP/common/inline.hpp"
#include "../submodules/SSP/common/util.hpp"

namespace {
  const std::string changeIntToGreek(const std::string& name)
  {
    if (name=="1")       return "I";
    else if (name=="2")  return "II";
    else if (name=="3")  return "III";
    else if (name=="4")  return "IV";
    else if (name=="5")  return "V";
    else if (name=="6")  return "VI";
    else if (name=="7")  return "VII";
    else if (name=="8")  return "VIII";
    else if (name=="9")  return "IX";
    else if (name=="10") return "X";
    else if (name=="11") return "XI";
    else if (name=="12") return "XII";
    else if (name=="13") return "XIII";
    else if (name=="14") return "XIV";
    else if (name=="15") return "XV";
    else if (name=="16") return "XVI";
    else return name;
  }
}

int32_t countmp(HashOfGeneDataMap &mp)
{
  int32_t n(0);
  for (auto pair: mp) {
    for (auto &x: mp.at(pair.first)) ++n;
  }
  return n;
}

std::vector<std::string> scanGeneName(const HashOfGeneDataMap &mp)
{
  std::vector<std::string> vgname;
  for (auto &pair: mp) {
    for (auto &x: mp.at(pair.first)) vgname.push_back(x.first);
  }
  return vgname;
}

HashOfGeneDataMap extract_mp(const HashOfGeneDataMap &tmp, const std::vector<std::string> glist)
{
  HashOfGeneDataMap mp;

  for (auto &x: glist) {
    for (auto &pair: tmp) {
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
  if (!in) PRINTERR_AND_EXIT("genelist file does not exist.");

  std::vector<std::string> glist;
  std::string lineStr;

  while (!in.eof()) {
    getline(in, lineStr);
    if (!lineStr.empty()) glist.push_back(lineStr);
  }

  return glist;
}

void parseARSOriDB(const std::string& fileName, HashOfGeneDataMap &mp)
{
  std::ifstream in(fileName);
  if (!in) PRINTERR_AND_EXIT("ARS file does not exist.");

  std::string lineStr;

  while (!in.eof()) {
    getline(in, lineStr);
    if (lineStr.empty() || isStr(lineStr, "Status")) continue;

    std::vector<std::string> v;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));

    if (v.size() < 7) {
      std::cout << lineStr << std::endl;
      PRINTERR_AND_EXIT("Error: the ARS file does not have 7 columns.");
    }

    std::string tname;
    if (isStr(v[2], "ARS")) tname = v[2];
    else tname = "ARS_" + v[2];
    std::string chr(changeIntToGreek(rmchr(v[4])));

    mp[chr][tname].chr     = chr;
    mp[chr][tname].gname   = tname;
    mp[chr][tname].strand  = "";
    mp[chr][tname].txStart = stoi(v[5]);
    mp[chr][tname].txEnd   = stoi(v[6]);
    mp[chr][tname].gtype   = "ARS";
  }
  return;
}

void parseTER(const std::string& fileName, HashOfGeneDataMap &mp)
{
  std::ifstream in(fileName);
  if (!in) PRINTERR_AND_EXIT("ARS file does not exist.");

  std::string lineStr;

  while (!in.eof()) {
    getline(in, lineStr);
    if (lineStr.empty() || isStr(lineStr, "#") || isStr(lineStr, "Name")) continue;

    std::vector<std::string> v;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));

    std::string tname(v[0]);
    std::string chr(changeIntToGreek(rmchr(v[1])));

    mp[chr][tname].chr     = chr;
    mp[chr][tname].gname   = tname;
    mp[chr][tname].strand  = "";
    mp[chr][tname].txStart = stoi(v[2]);
    mp[chr][tname].txEnd   = stoi(v[3]);
    mp[chr][tname].gtype   = "TER";
  }
  return;
}

HashOfGeneDataMap parseSGD(const std::string& fileName)
{
  if (isStr(fileName, ".gtf")) {
    std::cerr << "Warning: gene file seems to be gtf format but is parsed as refFlat." << std::endl;
  }

  std::ifstream in(fileName);
  if (!in) PRINTERR_AND_EXIT(fileName << " does not exist.");

  HashOfGeneDataMap tmp;
  std::string lineStr;

  while (!in.eof()) {
    getline(in, lineStr);
    if (lineStr.empty()) continue;

    std::vector<std::string> v;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));

    std::string chr(changeIntToGreek(rmchr(v[8])));
    std::string type(v[1]);
    std::string tname(v[0]);

    if (type == "ARS") tmp[chr][tname].gname = v[3];
    else if (type == "centromere") tmp[chr][tname].gname = "CEN_chr" + chr;
    else if (type == "teromere")   tmp[chr][tname].gname = v[3];
    else if (type == "ORF") {
      if (v[4] != "") tmp[chr][tname].gname = v[4];
      else tmp[chr][tname].gname = v[3];
    }
    else if (type == "tRNA") {
      if (v[4] != "") tmp[chr][tname].gname = v[4];
      else tmp[chr][tname].gname = v[3];
    }
    else if (type == "rRNA")   tmp[chr][tname].gname = "rRNA";
    else if (type == "snoRNA") tmp[chr][tname].gname = "snoRNA";
    else if (type == "long_terminal_repeat") tmp[chr][tname].gname = "LTR";
    else if (type == "repeat_region") tmp[chr][tname].gname = v[6];
    else if (type == "retrotransposon") tmp[chr][tname].gname = v[3];
    else continue;

    tmp[chr][tname].chr = chr;
    tmp[chr][tname].gtype = type;
    if (v[11] == "C") {
      tmp[chr][tname].txStart = stoi(v[10]);
      tmp[chr][tname].txEnd   = stoi(v[9]);
      tmp[chr][tname].strand  = "-";
    } else {
      tmp[chr][tname].txStart = stoi(v[9]);
      tmp[chr][tname].txEnd   = stoi(v[10]);
      tmp[chr][tname].strand  = "+";
    }
    if (type == "ARS" || type == "centromere"|| type == "teromere") tmp[chr][tname].strand = "";
  }
  return tmp;
}

HashOfGeneDataMap parseRefFlat(const std::string& fileName)
{
  if (isStr(fileName, ".gtf")) {
    std::cerr << "Warning: gene file seems to be gtf format but is parsed as refFlat." << std::endl;
  }

  bool UCSC(false);
  std::ifstream in(fileName);
  if (!in) PRINTERR_AND_EXIT("refFlat file does not exist.");

  HashOfGeneDataMap tmp;
  std::string lineStr;

  while (!in.eof()) {
    getline(in, lineStr);
    if (lineStr.empty() || lineStr[0] == '#') continue;

    std::vector<std::string> v, exonStarts, exonEnds;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));

    std::string tname(v[1]);
    std::string chr = rmchr(v[2]);
    if (isStr(tname, "NM") || isStr(tname, "NR")) UCSC = true;

    try {
      tmp[chr][tname].gname   = v[0];
      tmp[chr][tname].tname   = tname;
      tmp[chr][tname].chr     = chr;
      tmp[chr][tname].strand  = v[3];
      tmp[chr][tname].txStart = stoi(v[4]);
      tmp[chr][tname].txEnd   = stoi(v[5]);
      tmp[chr][tname].cdsStart = stoi(v[6]);
      tmp[chr][tname].cdsEnd   = stoi(v[7]);
      tmp[chr][tname].exonCount = stoi(v[8]);
      if (v.size() >= 13) tmp[chr][tname].gtype = v[12];
      else if (v.size() >= 12) tmp[chr][tname].gtype = v[11];
      else if (UCSC) {
        if (isStr(tname, "NM")) tmp[chr][tname].gtype = "protein_coding";
        else tmp[chr][tname].gtype = "noncoding RNA";
      }

      boost::split(exonStarts, v[9], boost::algorithm::is_any_of(","));
      boost::split(exonEnds,  v[10], boost::algorithm::is_any_of(","));

      for (int32_t i=0; i<tmp[chr][tname].exonCount; i++) {
        tmp[chr][tname].exon.emplace_back(stoi(exonStarts[i]), stoi(exonEnds[i]));
      }
    } catch (std::exception &e) {
      PRINTERR_AND_EXIT("invalid columns in refFlat format. " + std::string(e.what()));
    }
  }
  return tmp;
}

HashOfGeneDataMap parseGtf(const std::string& fileName)
{
  if (!isStr(fileName, ".gtf")) {
    std::cerr << "Warning: gene file may not be gtf format but is parsed as gtf." << std::endl;
  }

  std::ifstream in(fileName);
  if (!in) PRINTERR_AND_EXIT("gtf file does not exist.");

  HashOfGeneDataMap tmp;
  std::string lineStr;

  while (!in.eof()) {
    getline(in, lineStr);
    if (lineStr.empty() || lineStr[0] == '#') continue;
    std::vector<std::string> v;

    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    std::string feat(v[2]);
    if (feat == "gene" || feat == "transcript" || feat == "three_prime_utr" || feat == "five_prime_utr") continue;

    std::string chr(rmchr(v[0]));
    int32_t start(stoi(v[3]));
    int32_t end(stoi(v[4]));
    std::string strand(v[6]);
    std::string annotation(v[8]);

    std::string gname, tname, gid, tid, gsrc, gtype, tsrc, ttype, ttag="";
    std::vector<std::string> idtab, vc;
    boost::split(idtab, annotation, boost::algorithm::is_any_of(";"));
    for (auto &term: idtab) {
      boost::split(vc, term, boost::algorithm::is_any_of("\""));
      if (isStr(vc[0], "gene_source"))             gsrc  = vc[1];
      else if (isStr(vc[0], "gene_biotype"))       gtype = vc[1];
      else if (isStr(vc[0], "transcript_source"))  tsrc  = vc[1];
      else if (isStr(vc[0], "transcript_biotype")) ttype = vc[1];
      else if (isStr(vc[0], "transcript_name")) tname = vc[1];
      else if (isStr(vc[0], "gene_name"))       gname = vc[1];
      else if (isStr(vc[0], "transcript_id"))   tid = vc[1];
      else if (isStr(vc[0], "gene_id"))         gid = vc[1];
      else if (isStr(vc[0], "tag")) {
        if (ttag=="" || vc[1]=="CCDS")           ttag = vc[1];
        else if (vc[1] == "basic" && ttag != "CCDS") ttag = vc[1];
        else if (ttag  != "basic" && ttag != "CCDS") ttag = vc[1];
      }
    }

    if (tname =="") continue;

    tmp[chr][tid].tname  = tname;
    tmp[chr][tid].gname  = gname;
    tmp[chr][tid].tid    = tid;
    tmp[chr][tid].gid    = gid;
    tmp[chr][tid].chr    = chr;
    tmp[chr][tid].strand = strand;
    if (feat == "start_codon") {
      if (strand == "+") tmp[chr][tid].cdsStart = start;
      else tmp[chr][tid].cdsEnd = end;
    } else if (feat == "stop_codon") {
      if (strand == "+") tmp[chr][tid].cdsEnd = end;
      else tmp[chr][tid].cdsStart = start;
    } else if (feat == "exon") {
      ++tmp[chr][tid].exonCount;
      if (!tmp[chr][tid].txStart || start < tmp[chr][tid].txStart) tmp[chr][tid].txStart = start;
      if (end > tmp[chr][tid].txEnd) tmp[chr][tid].txEnd = end;
      tmp[chr][tid].exon.emplace_back(start, end);
    }

    tmp[chr][tid].gsrc = gsrc;
    tmp[chr][tid].tsrc = tsrc;
    tmp[chr][tid].gtype = gtype;
    tmp[chr][tid].ttype = ttype;
    tmp[chr][tid].ttag = ttag;
  }

  // "start_codon", "stop_codon"がないとcdsStart, cdsEndが0になる
  for (auto &pair: tmp) {
    for (auto &x: pair.second) {
      if (!x.second.cdsStart) x.second.cdsStart = x.second.txStart;
      if (!x.second.cdsEnd)   x.second.cdsEnd   = x.second.txEnd;
    }
  }

  return tmp;
}

HashOfGeneDataMap construct_gmp(const HashOfGeneDataMap &tmp)
{
  HashOfGeneDataMap gmp;

  for (auto pair: tmp) {
    std::string chr(pair.first);
    for (auto x: pair.second) {
      std::string gname(x.second.gid);
      if (gmp.find(chr) == gmp.end() || gmp[chr].find(gname) == gmp[chr].end()) {
        gmp[chr][gname] = x.second;
      }
      else if (x.second.ttag == "CCDS") {
        if (gmp[chr][gname].ttag != "CCDS" || gmp[chr][gname].length() < x.second.length()) gmp[chr][gname] = x.second;
      }
      else if (x.second.ttag == "basic") {
        if (gmp[chr][gname].ttag != "CCDS" &&
            ((gmp[chr][gname].ttag == "basic" && gmp[chr][gname].length() < x.second.length())
             || gmp[chr][gname].ttag != "basic")) gmp[chr][gname] = x.second;
      }
      else if (gmp[chr][gname].ttag != "basic"
               && (x.second.ttag == "basic" || gmp[chr][gname].length() < x.second.length())) {
        gmp[chr][gname] = x.second;
      }
      else if (x.second.ttag == "basic" || gmp[chr][gname].length() < x.second.length()) gmp[chr][gname] = x.second;
    }
  }
  return gmp;
}

void printMap(const HashOfGeneDataMap &mp)
{
  for (auto &pair: mp) {
    for (auto &x: pair.second) {
      x.second.printall();
      std::cout << std::endl;
    }
  }
  return;
}

bool isGeneUCSC(const HashOfGeneDataMap &mp)
{
  for (auto &pair: mp) {
    for (auto &x: pair.second) {
      if (isStr(x.second.tname, "ENS")) return false;
    }
  }
  return true;
}

std::vector<chrsize> readGenomeTable(const std::string& fileName)
{
  std::ifstream in(fileName);
  if (!in) PRINTERR_AND_EXIT("genometable file does not exist.");

  std::vector<chrsize> gt;
  std::string lineStr;

  while (!in.eof()) {
    getline(in, lineStr);
    if (lineStr.empty() || lineStr[0] == '#') continue;
    std::vector<std::string> v;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    gt.emplace_back(v[0], stoi(v[1]));
  }
  // Greekchr
  /*  for (auto &x: gt) {
    if (x.getname() == "I") {
      for (auto &x:gt) x.Greekchron();
      break;
    }
    }*/
  return gt;
}
