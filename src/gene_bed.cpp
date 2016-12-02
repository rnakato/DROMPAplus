#include "gene_bed.h"
#include "macro.h"

void mergeArray(std::vector<strRange> &array)
{
  auto itr = array.begin();
  while (itr != array.end()) {
    auto next = itr+1;
    if(next != array.end() && next->strand == itr->strand &&
       overlap(itr->start, itr->end, next->start, next->end)) {
      next->start = std::min(itr->start, next->start);
      next->end   = std::max(itr->end, next->end);
      itr = array.erase(itr);
    } else {
      itr++;
    }
  }

  return;
}

void scanConvergent(std::vector<convsite> &vconv, const std::vector<strRange> array, int limit)
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

std::vector<convsite> gen_convergent(const int limconv, const HashOfGeneDataMap &mp)
{
  std::vector<convsite> vconv;
  for(auto itr = mp.begin(); itr != mp.end(); ++itr){
    std::string chr(itr->first);
    if(mp.find(chr) == mp.end()) continue;

    std::vector<strRange> array;
    for(auto pg = mp.at(chr).begin(); pg != mp.at(chr).end(); ++pg) {
      strRange r(pg->second.txStart, pg->second.txEnd, chr, pg->second.strand, pg->second);
      array.push_back(r);
    }
    sort(array.begin(), array.end(),
	 [](const strRange& x, const strRange& y) { return x.start < y.start; });
    mergeArray(array);
    scanConvergent(vconv, array, limconv);
  }
  return vconv;
}
