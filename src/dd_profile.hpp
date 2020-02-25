/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_PROFILE_H_
#define _DD_PROFILE_H_

#include "dd_gv.hpp"

class ReadProfile {
  int32_t stype;
  std::string xlabel;
  std::string Rscriptname;
  std::string Rfigurename;

protected:
  int32_t binsize;
  int32_t nbin;
  int32_t width_from_center;
  int32_t binwidth_from_center;

  int32_t nsites;
  int32_t nsites_skipped;

  std::string RDataname;
  std::unordered_map<std::string, std::vector<double>> hprofile;

  std::vector<genedata> get_garray(const GeneDataMap &mp)
  {
    std::vector<genedata> garray;
    for (auto &m: mp) {
      if (m.second.gtype == "nonsense_mediated_decay" ||
	  m.second.gtype == "processed_transcript"    ||
	  m.second.gtype == "retained_intron") continue;
      garray.emplace_back(m.second);
    }
    return garray;
  }

  int32_t isExceedRange(const int32_t posi, const int32_t chrlen) {
    return posi - width_from_center < 0 || posi + width_from_center >= chrlen;
  }

  void WriteValAroundPosi(std::ofstream &out, const SamplePairOverlayed &pair, const vChrArray &vReadArray,
			  const int32_t posi, const std::string &strand)
  {
    int32_t bincenter(posi/binsize);
    int32_t sbin(bincenter - binwidth_from_center);
    int32_t ebin(bincenter + binwidth_from_center);

    if (strand == "+") {
      for (int32_t i=sbin; i<=ebin; ++i) out << "\t" << getSumVal(pair, vReadArray, i, i);
    } else {
      for (int32_t i=ebin; i>=sbin; --i) out << "\t" << getSumVal(pair, vReadArray, i, i);
    }
  }

  double getSumVal(const SamplePairOverlayed &pair, const vChrArray &vReadArray,
		   const int32_t sbin, const int32_t ebin) {
    double value(-1);
    if (!stype) { // ChIP read
      double sumIP(0);
      for (int32_t i=sbin; i<=ebin; ++i) sumIP += vReadArray.getArray(pair.first.argvChIP).array[i];
      value = getratio(sumIP, (ebin - sbin + 1));
    } else if (stype == 1) { // ChIP/Input enrichment
      double sumRatio(0);
      for (int32_t i=sbin; i<=ebin; ++i) {
	double r = getratio(vReadArray.getArray(pair.first.argvChIP).array[i],
			     vReadArray.getArray(pair.first.argvInput).array[i]);
//	printf("%f, %f, r=%f\n",vReadArray.getArray(pair.first.argvChIP).array[i],
//	       vReadArray.getArray(pair.first.argvInput).array[i],
//	       r);
	sumRatio += r;
      }
      value = getratio(sumRatio, (ebin - sbin + 1));
    }
    return value;
  }

public:
  ReadProfile(const DROMPA::Global &p, const int32_t _nbin=0):
    stype(p.prof.stype),
    width_from_center(p.prof.width_from_center),
    nsites(0), nsites_skipped(0)
  {

    if(p.prof.isPtypeTSS())          xlabel = "Distance from TSS (bp)";
    else if(p.prof.isPtypeTTS())     xlabel = "Distance from TES (bp)";
    else if(p.prof.isPtypeGene100()) xlabel = "% gene length from TSS";
    else if(p.prof.isPtypeBed())     xlabel = "Distance from the peak summit (bp)";

    for (auto &x: p.samplepair) binsize = x.first.getbinsize();
    binwidth_from_center = width_from_center / binsize;
    if (binwidth_from_center <= 1) {
      PRINTERR_AND_EXIT("please specify larger size for --widthfromcenter:" << width_from_center
	       << " than binsize:" << binsize << ".");
    }

    if (_nbin) nbin = _nbin;
    else nbin = binwidth_from_center * 2 +1;

    std::vector<double> array(nbin, 0);
    for (auto &x: p.samplepair) hprofile[x.first.argvChIP] = array;
  }

  void setOutputFilename(const DROMPA::Global &p) {
    std::string prefix;
    if (stype==0)      prefix = p.getPrefixName() + ".ChIPread";
    else if (stype==1) prefix = p.getPrefixName() + ".Enrichment";
    Rscriptname = prefix + ".R";
    RDataname   = prefix;    //  prefix + ".tsv";
    Rfigurename = prefix + ".pdf";
    for (auto &x: {Rscriptname, RDataname, Rfigurename}) unlink(x.c_str());
    for (auto &x: p.samplepair) {
      std::string file(RDataname + "." + x.first.label + ".tsv");
//      std::cout << file << std::endl;
      unlink(file.c_str());
    }
    for (auto &x: {Rscriptname, Rfigurename}) unlink(x.c_str());
  }

  void MakeFigure(const DROMPA::Global &p) {
    std::cout << "\nMake figure.." << std::endl;
    std::ofstream out(Rscriptname);

    std::vector<double> vcol({CLR_RED, CLR_BLUE, CLR_GREEN, CLR_LIGHTCORAL, CLR_BLACK, CLR_PURPLE, CLR_GRAY3, CLR_OLIVE, CLR_YELLOW3, CLR_SLATEGRAY, CLR_PINK, CLR_SALMON, CLR_GREEN2, CLR_BLUE3, CLR_PURPLE2, CLR_DARKORANGE});

    for (size_t i=1; i<=p.samplepair.size(); ++i) {
      out << "t <- read.table('"
	  << RDataname << "." << p.samplepair[i-1].first.label << ".tsv"
	  << "', header=F, sep='\\t', quote='')" << std::endl;
      out << "t <- t[,-1]" << std::endl;
      out << "n <- nrow(t)" << std::endl;


      out << boost::format("t%1% <- t\n") % i;
      out << boost::format("colnames(t%1%) <- t%1%[1,]\n") % i;
      out << boost::format("t%1% <- t%1%[-1,]\n") % i;
      out << boost::format("p%1% <- apply(t%1%,2,mean)\n") % i;
      out << boost::format("sd%1% <- apply(t%1%,2,sd)\n") % i;
      out << boost::format("SE%1% <- 1.96 * sd%1%/sqrt(n)\n") % i;
      out << boost::format("p%1%_upper <- p%1% + SE%1%\n") % i;
      out << boost::format("p%1%_lower <- p%1% - SE%1%\n") % i;
      out << boost::format("x <- as.integer(colnames(t%1%))\n") % i;
    }
    out << "pdf('" << Rfigurename << "',6,6)" << std::endl;
    out << "plot(x,p1,type='l',col=rgb(" << vcol[0] << "," << vcol[1] << "," << vcol[2] << ")";
    if (p.prof.isPtypeGene100()) out << ",log='y'";
    if (!stype) {
      out << ",xlab='" << xlabel << "',ylab='Read density')" << std::endl;
    } else if (stype==1) {
      out << ",xlab='" << xlabel << "',ylab='Read enrichment')" << std::endl;
    }
    out << "polygon(c(x, rev(x)), c(p1_lower, rev(p1_upper)), col=rgb("
	<< vcol[0] << "," << vcol[1] << "," << vcol[2] << ",0.3), border=NA)" << std::endl;

    for (size_t i=1; i<p.samplepair.size(); ++i) {
      std::string colstr =
	std::to_string(vcol[i*3]) + ","
	+ std::to_string(vcol[i*3 +1]) + ","
	+ std::to_string(vcol[i*3 +2]) + "";
      out << boost::format("lines(x,p%1%,col=rgb(%2%))\n") % (i+1) % colstr;
      out << boost::format("polygon(c(x, rev(x)), c(p%1%_lower, rev(p%1%_upper)), col=rgb(%2%,0.3), border=NA)\n") % (i+1) % colstr;
    }

    out << "legend('bottomleft',c(";
    for (size_t i=0; i<p.samplepair.size(); ++i) {
      out << "'" << p.samplepair[i].first.label << "'";
      if (i < p.samplepair.size()-1) out << ",";
    }
    out << "),col=c(rgb(" << vcol[0] << "," << vcol[1] << "," << vcol[2] <<")";
    for (size_t i=1; i<p.samplepair.size(); ++i) {
      out << ",rgb(" << vcol[i*3] << "," << vcol[i*3 +1] << "," << vcol[i*3 +2] << ")";
    }
    out << "), lwd=1.5)" << std::endl;

    std::string command("R --vanilla < " + Rscriptname);
    if(system(command.c_str())) PRINTERR_AND_EXIT("Rscript execution failed.");
  }

  virtual void WriteTSV_EachChr(const DROMPA::Global &p, const chrsize &chr)=0;

  void printHead(const DROMPA::Global &p) {
    for (auto &x: p.samplepair) {
      std::string file(RDataname + "." + x.first.label + ".tsv");
      std::ofstream out(file);
      for (int32_t i=-binwidth_from_center; i<=binwidth_from_center; ++i) {
	out << "\t" << (i*binsize);
      }
      out << std::endl;
      out.close();
    }
  }

  void printNumOfSites() const {
    std::cout << "\n\nthe number of sites: " << nsites << std::endl;
    std::cout << "the number of skipped sites: " << nsites_skipped << std::endl;
  }
};


class ProfileTSS: public ReadProfile {

public:
  explicit ProfileTSS(const DROMPA::Global &p):
    ReadProfile(p) {}

  void WriteTSV_EachChr(const DROMPA::Global &p, const chrsize &chr) {
    DEBUGprint("ProfileTSS:WriteTSV_EachChr");
    if(p.anno.genefile == "") PRINTERR_AND_EXIT("Please specify --gene.");

    std::string chrname(rmchr(chr.getname()));
    if (p.anno.gmp.find(chrname) == p.anno.gmp.end()) return;

    vChrArray vReadArray(p, chr);

    for (auto &x: p.samplepair) {
      std::string file(RDataname + "." + x.first.label + ".tsv");
      std::ofstream out(file, std::ios::app);

      auto gmp(get_garray(p.anno.gmp.at(chrname)));

      for (auto &gene: gmp) {
	++nsites;

	int32_t position(0);
	if (p.prof.isPtypeTSS()) {
	  if (gene.strand == "+") position = gene.txStart;
	  else                    position = gene.txEnd;
	} else if (p.prof.isPtypeTTS()) {
	  if (gene.strand == "+") position = gene.txEnd;
	  else                    position = gene.txStart;
	}
	if (isExceedRange(position, chr.getlen())) {
	  ++nsites_skipped;
	  continue;
	}
	out << gene.tname;
	WriteValAroundPosi(out, x, vReadArray, position, gene.strand);
	out << std::endl;
      }

      out.close();
    }
  }
};

class ProfileGene100: public ReadProfile {
  enum {GENEBLOCKNUM=100};

  void outputEachGene(std::ofstream &out, const SamplePairOverlayed &x,
		      const genedata &gene, const vChrArray &vReadArray, int32_t len) {
    int32_t s,e;
    double len100(len / (double)GENEBLOCKNUM);

    for (int32_t i=0; i<nbin; ++i) {
      if (gene.strand == "+") {
	s = (gene.txStart - len + len100 *i)       / binsize;
	e = (gene.txStart - len + len100 *(i+1) -1)/ binsize;
      }else{
	s = (gene.txEnd + len - len100 * (i+1))/ binsize;
	e = (gene.txEnd + len - len100 * i -1) / binsize;
      }
      out << "\t" << getSumVal(x, vReadArray, s, e);
    }
  }

public:
  explicit ProfileGene100(const DROMPA::Global &p):
    ReadProfile(p, GENEBLOCKNUM * 3)
  {}

  void WriteTSV_EachChr(const DROMPA::Global &p, const chrsize &chr) {
    DEBUGprint("ProfileGene100:WriteTSV_EachChr");
    if(p.anno.genefile == "") PRINTERR_AND_EXIT("Please specify --gene.");

    std::string chrname(rmchr(chr.getname()));
    if (p.anno.gmp.find(chrname) == p.anno.gmp.end()) return;

    vChrArray vReadArray(p, chr);
//    std::ofstream out(RDataname, std::ios::app);


    for (auto &x: p.samplepair) {
      std::string file(RDataname + "." + x.first.label + ".tsv");
      std::ofstream out(file, std::ios::app);

      auto gmp(get_garray(p.anno.gmp.at(chrname)));
      for (auto &gene: gmp) {
	++nsites;

	int32_t len(gene.length());
	if (len < 1000 || gene.txEnd + len >= chr.getlen() || gene.txStart - len < 0) {
	  ++nsites_skipped;
	  continue;
	}

	out << gene.gname;
	outputEachGene(out, x, gene, vReadArray, len);
	out << std::endl;
      }
    }
  }

  void printHead(const DROMPA::Global &p) {
    for (auto &x: p.samplepair) {
      std::string file(RDataname + "." + x.first.label + ".tsv");
      std::ofstream out(file);
      for (int32_t i=-GENEBLOCKNUM; i<GENEBLOCKNUM*2; ++i) {
	out << "\t" << i;
      }
      out << std::endl;
      out.close();
    }
  }
};


class ProfileBedSites: public ReadProfile {
public:
  explicit ProfileBedSites(const DROMPA::Global &p):
    ReadProfile(p) {}

  void WriteTSV_EachChr(const DROMPA::Global &p, const chrsize &chr) {
    DEBUGprint("ProfileBedSites:WriteTSV_EachChr");

    if(!p.anno.vbedlist.size()) PRINTERR_AND_EXIT("Please specify --bed.");

    vChrArray vReadArray(p, chr);

    for (auto &x: p.samplepair) {
      std::string file(RDataname + "." + x.first.label + ".tsv");
      std::ofstream out(file);

      for (auto &vbed: p.anno.vbedlist) {
	for (auto &bed: vbed.getvBed()) {
	  if (bed.chr != chr.getname()) continue;
	  ++nsites;
	  if (isExceedRange(bed.summit, chr.getlen())) {
	    ++nsites_skipped;
	    continue;
	  }
	  out << bed.getSiteStr();
	  WriteValAroundPosi(out, x, vReadArray, bed.summit, "+");
	  out << std::endl;
	}
      }

      out.close();
    }

  }
};

#endif /* _DD_PROFILE_H_ */
