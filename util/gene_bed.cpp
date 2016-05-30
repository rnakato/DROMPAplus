#include "gene_bed.h"
#include "macro.h"

void print_gdist(const cmdline::parser p, gdist n, string str)
{
  cout << str << "\t";
  if(p.exist("intron")) cout << boost::format("%1%\t%2%\t%3%\t%4%\t%5%\t%6%") % n.genome % n.up % n.down % n.exon % n.intron % n.inter;
  else cout << boost::format("%1%\t%2%\t%3%\t%4%\t%5%") % n.genome % n.up % n.genic % n.intron % n.inter;
  if(p.exist("conv")) cout << boost::format("\t%1%\t%2%\t%3%\n") % n.conv % n.div % n.par;
  else cout << endl;

  cout << str << " (%)\t";
  if(p.exist("intron")) printf("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f", 100.0, n.ratio(UPSTREAM), n.ratio(DOWNSTREAM), n.ratio(EXON), n.ratio(INTRON), n.ratio(INTERGENIC));
  else printf("%.2f\t%.2f\t%.2f\t%.2f\t%.2f", 100.0, n.ratio(UPSTREAM), n.ratio(DOWNSTREAM), n.ratio(GENIC), n.ratio(INTERGENIC));
  if(p.exist("conv")) printf("\t%.2f\t%.2f\t%.2f\n", n.ratio(CONVERGENT), n.ratio(DIVERGENT), n.ratio(PARALLEL));
  else cout << endl;
  return;
}

void mergeArray(vector<strRange> &array)
{
  auto itr = array.begin();
  while (itr != array.end()) {
    auto next = itr+1;
    if(next != array.end() && next->strand == itr->strand &&
       overlap(itr->start, itr->end, next->start, next->end)) {
      next->start = min(itr->start, next->start);
      next->end   = max(itr->end, next->end);
      itr = array.erase(itr);
    } else {
      itr++;
    }
  }

  return;
}

void scanConvergent(vector<convsite> &vconv, const vector<strRange> array, int limit)
{
  for(auto itr = array.begin(); itr != array.end(); ++itr) {
    for(auto itr2= itr+1; itr2 != array.end(); ++itr2) {
      int dist = itr2->start - itr->end;
      if(0 < dist && dist <= limit){
	convsite conv(itr->end, itr2->start, itr->chr, itr->gene);
	if(itr->strand == itr2->strand) conv.st = PARALLEL;
	else if(itr->strand == "+" && itr2->strand == "-") conv.st = CONVERGENT;
	else if(itr->strand == "-" && itr2->strand == "+") conv.st = DIVERGENT;
	vconv.push_back(conv);
      }else if(dist > limit) break;
    }
  }
  return;
}

vector<convsite> gen_convergent(const cmdline::parser p, const unordered_map<string, unordered_map<string, genedata>> &mp)
{
  vector<convsite> vconv;
  for(auto itr = mp.begin(); itr != mp.end(); ++itr){
    string chr(itr->first);
    if(mp.find(chr) == mp.end()) continue;

    vector<strRange> array;
    for(auto pg = mp.at(chr).begin(); pg != mp.at(chr).end(); ++pg) {
      //      if(pg->second.txStart==59206450 || pg->second.txStart==59210641 || pg->second.txStart==59206450 || pg->second.txStart==59213948) cout << pg->second.txStart << "." <<  pg->second.txEnd << ","<< chr << "." << pg->second.strand << endl;
      strRange r(pg->second.txStart, pg->second.txEnd, chr, pg->second.strand, pg->second);
      array.push_back(r);
    }
    sort(array.begin(), array.end(),
	 [](const strRange& x, const strRange& y) { return x.start < y.start; });
    mergeArray(array);
    scanConvergent(vconv, array, p.get<int>("limconv"));
  }
  return vconv;
}
