/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include "readdata.h"
#include "macro.h"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <sstream>

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
  if(!in) PRINTERR("genelist file does not exist.");

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
  if(fileName.find(".gtf") != string::npos) {
    cerr << "Warning: gene file seems to be gtf format but is parsed as refFlat." << endl;
  }

  ifstream in(fileName);
  if(!in) PRINTERR("refFlat file does not exist.");

  unordered_map<string, unordered_map<string, genedata>> tmp;
  string lineStr;
  
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty()) continue;

    vector<string> v, exonStarts, exonEnds;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));

    string tname(v[1]);
    string chr = rmchr(v[2]);
    
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
  if(fileName.find(".gtf") == string::npos) {
    cerr << "Warning: gene file may not be gtf format but is parsed as gtf." << endl;
  }

  ifstream in(fileName);
  if(!in) PRINTERR("gtf file does not exist.");

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
    
    string chr  = rmchr(v[0]);
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
  return;
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

vector<chrsize> read_genometable(const string& fileName)
{
  ifstream in(fileName);
  if(!in) PRINTERR("genometable file does not exist.");

  vector<chrsize> gt;
  string lineStr;
  
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || lineStr[0] == '#') continue;
    vector<string> v;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    chrsize temp;
    temp.name = rmchr(v[0]);
    temp.len = stoi(v[1]);
    gt.push_back(temp);
  }
  return gt;
}

vector<int> readMpbl(string mpfile, string chrname, int binsize, int nbin)
{
  string filename = mpfile + "/map_fragL150_" + chrname + "_bin" + IntToString(binsize) +".txt";
  vector<int> mparray(nbin, 0);

  isFile(filename);
  ifstream in(filename);

  string lineStr;
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty()) continue;
    vector<string> v;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));

    int n(stoi(v[0])/binsize);
    double val(stof(v[1])*binsize);
    mparray[n] = val;
  }

  return mparray;
}

vector<char> readMpbl_binary(string mpfile, string chrname, int chrlen)
{
  string filename = mpfile + "/map_" + chrname + "_binary.txt";
  vector<char> mparray(chrlen, UNMAPPABLE);

  isFile(filename);
  int n(0);
  char c;
  ifstream in(filename);
  while (!in.eof()) {
    c = in.get();
    if(c==' ') continue;
    if(c=='1') mparray[n]=MAPPABLE;
    ++n;
    if(n >= chrlen-1) break;
  }

  return mparray;
}

vector<char> readMpbl_binary(int chrlen)
{
  vector<char> mparray(chrlen, MAPPABLE);
  return mparray;
}

vector<char> arraySetBed(vector<char> &array, string chrname, vector<bed> vbed)
{
  for(auto bed: vbed) {
    if(bed.chr == chrname) {
      int s(bed.start);
      int e(bed.end);
      if(e>=(int)array.size()) {
	cerr << "Warning: bedfile" << bed.start <<"-"<<bed.end << " > array size " << array.size()<< endl;
	e = array.size()-1;
      }
      for(int i=s; i<=e; ++i) array[i] = INBED;
    } 
  }

  return array;
}

void isFile(string str)
{
  boost::filesystem::path const file(str);
  if(!boost::filesystem::exists(file)) PRINTERR(str << " does not exist.");
}

string IntToString(int n)
{
  ostringstream stream;
  stream << n;
  return stream.str();
}

std::string rmchr(const std::string &chr)
{
  std::string s;
  if(!chr.find("chr")) s = chr.substr(3);
  else s = chr;
  return s;
}
