#include "readdata.h"
#include "warn.h"
#include <boost/algorithm/string.hpp>

using namespace std;

int countmp(unordered_map<string, unordered_map<string, genedata>> &mp)
{
  int n(0);
  for(auto itr = mp.begin(); itr != mp.end(); ++itr) {
    for(auto itr2 = mp.at(itr->first).begin(); itr2 != mp.at(itr->first).end(); ++itr2) n++;
  }
  return n;
}

vector<string> scanGeneName(const unordered_map<string, unordered_map<string, genedata>> &mp)
{
  vector<string> vgname;
  for(auto itr = mp.begin(); itr != mp.end(); ++itr) {
    for(auto itr2 = mp.at(itr->first).begin(); itr2 != mp.at(itr->first).end(); ++itr2) {
      vgname.push_back(itr2->first);
    }
  }
  return vgname;
}

unordered_map<string, unordered_map<string, genedata>>
  extract_mp(const unordered_map<string, unordered_map<string, genedata>> &tmp,
	     const vector<string> glist)
{
  unordered_map<string, unordered_map<string, genedata>> mp;

  for(auto x: glist) {
    for(auto itr = tmp.begin(); itr != tmp.end(); ++itr) {
      string chr = itr->first;
      if(tmp.at(chr).find(x) != tmp.at(chr).end()){
	mp[chr][x] = tmp.at(chr).at(x);
	break;
      }
    }
  }

  return mp;
}

vector<string> readGeneList(const string& fileName)
{
  ifstream in(fileName);
  if(!in) printerr("genelist file does not exist.");

  vector<string> glist;
  string lineStr;
  
  while (!in.eof()) {
    getline(in, lineStr);
    if(!lineStr.empty()) glist.push_back(lineStr);
  }

  return glist;
}

unordered_map<string, unordered_map<string, genedata>> parseRefFlat(const string& fileName)
{
  ifstream in(fileName);
  if(!in) printerr("refFlat file does not exist.");

  unordered_map<string, unordered_map<string, genedata>> tmp;
  string lineStr;
  
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty()) continue;

    vector<string> v, exonStarts, exonEnds;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));

    string tname(v[1]);
    string chr = addchr(v[2]);
    
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

    for(int i=0; i<tmp[chr][tname].exonCount; i++){
      range exon(stoi(exonStarts[i]), stoi(exonEnds[i]));
      tmp[chr][tname].exon.push_back(exon);
    }
  }

  return tmp;
}

unordered_map<string, unordered_map<string, genedata>> parseGtf(const string& fileName, const int nameflag)
{
  ifstream in(fileName);
  if(!in) printerr("gtf file does not exist.");

  unordered_map<string, unordered_map<string, genedata>> tmp;
  string lineStr;
  
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || lineStr[0] == '#') continue;
    string gname, tname;
    vector<string> v;
    
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    string feat = v[2];
    if(feat == "gene" || feat == "transcript") continue;
    if(feat == "three_prime_utr" || feat == "five_prime_utr") continue;
    
    string chr  = addchr(v[0]);
    int start   = stoi(v[3]);
    int end     = stoi(v[4]);
    string strand = v[6];
    string id   = v[8];
    string gsrc, gtype, tsrc, ttype;

    vector<string> idtab, vc;
    boost::split(idtab, id, boost::algorithm::is_any_of(";"));
    for (auto term: idtab) {
      boost::split(vc, term, boost::algorithm::is_any_of("\""));
      if(term.find("gene_source") != string::npos)             gsrc  = vc[1];
      else if(term.find("gene_biotype") != string::npos)       gtype = vc[1];
      else if(term.find("transcript_source") != string::npos)  tsrc  = vc[1];
      else if(term.find("transcript_biotype") != string::npos) ttype = vc[1];
      else{
	if(nameflag) {
	  if(term.find("transcript_name") != string::npos) tname = vc[1];
	  else if(term.find("gene_name")  != string::npos) gname = vc[1];
	} else {
	  if(term.find("transcript_id") != string::npos) tname = vc[1];
	  else if(term.find("gene_id")  != string::npos) gname = vc[1];
	}
      }
    }

    tmp[chr][tname].tname  = tname;
    tmp[chr][tname].gname  = gname;
    tmp[chr][tname].chr    = chr;
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
      if(end > tmp[chr][tname].txEnd) tmp[chr][tname].txEnd   = end;
      range exon(start, end);
      tmp[chr][tname].exon.push_back(exon);
    }

    tmp[chr][tname].gsrc = gsrc;
    tmp[chr][tname].tsrc = tsrc;
    tmp[chr][tname].gtype = gtype;
    tmp[chr][tname].ttype = ttype;
  }

  return tmp;
}

unordered_map<string, unordered_map<string, genedata>> construct_gmp(const unordered_map<string, unordered_map<string, genedata>> &tmp)
{
  unordered_map<string, unordered_map<string, genedata>> gmp;

  for(auto itr = tmp.begin(); itr != tmp.end(); ++itr) {
    string chr = itr->first;
    for(auto itr2 = itr->second.begin(); itr2 != itr->second.end(); ++itr2) {
      string gname = itr2->second.gname;
      if(gmp.find(chr) == gmp.end() || gmp[chr].find(gname) == gmp[chr].end()) gmp[chr][gname] = itr2->second;
      else if((itr2->second.tsrc.empty() || itr2->second.tsrc == "ensembl_havana") && (gmp[chr][gname].length() < itr2->second.length())) {
	 gmp[chr][gname] = itr2->second;
      }
    }
  }
  return gmp;
}

void printMap(const unordered_map<string, unordered_map<string, genedata>> &mp)
{
  for(auto itr = mp.begin(); itr != mp.end(); ++itr) {
    for(auto itr2 = itr->second.begin(); itr2 != itr->second.end(); ++itr2) {
      itr2->second.printall();
      cout << endl;
    }
  }
}

void printRefFlat(const unordered_map<string, unordered_map<string, genedata>> &mp)
{
  for(auto itr = mp.begin(); itr != mp.end(); ++itr) {
    for(auto itr2 = itr->second.begin(); itr2 != itr->second.end(); ++itr2) {
      cout << itr2->second.gname << "\t"
	   << itr2->first << "\t"
	   << itr->first << "\t"
	   << itr2->second.strand << "\t"
	   << itr2->second.txStart << "\t"
	   << itr2->second.txEnd << "\t";
      if(itr2->second.cdsStart) {
	cout << itr2->second.cdsStart << "\t"
	     << itr2->second.cdsEnd   << "\t";
      } else {
	cout << itr2->second.txEnd << "\t"
	     << itr2->second.txEnd << "\t";
      }
      cout << itr2->second.exonCount << "\t";
      for (auto x: itr2->second.exon) cout << x.start << ",";
      cout << "\t";
      for (auto x: itr2->second.exon) cout << x.end   << ",";
      cout << endl;
    }
  }
  return;
}

map<string, int> read_genometable(const string& fileName)
{
  ifstream in(fileName);
  if(!in) printerr("genometable file does not exist.");

  map<string, int> gt;
  string lineStr;
  
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || lineStr[0] == '#') continue;
    vector<string> v;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    string chr = addchr(v[0]);
    gt[chr] = stoi(v[1]);
  }
  return gt;
}