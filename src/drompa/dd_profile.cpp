/* Copyright(c)  Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "dd_profile.hpp"
#include "color.hpp"

namespace {
  double getReadVal(const vChrArray &vReadArray, const SamplePairOverlayed &pair,
                    const int32_t i, const int32_t stype)
  {
    double val(0);
    if (!stype) {   // ChIP read
      val = vReadArray.getArray(pair.first.argvChIP).array[i];
    } else if (stype == 1) {   // ChIP/Input enrichment
      val = getratio(vReadArray.getArray(pair.first.argvChIP).array[i],
                     vReadArray.getArray(pair.first.argvInput).array[i]);
    }
    return val;
  }
}

std::vector<genedata> ReadProfile::get_garray(const GeneDataMap &mp)
{
  std::vector<genedata> garray;
  for (auto &m: mp) {
    if (m.second.gtype == "nonsense_mediated_decay"
        || m.second.gtype == "processed_transcript"
        || m.second.gtype == "retained_intron") continue;
    garray.emplace_back(m.second);
  }
  return garray;
}

void ReadProfile::WriteValAroundPosi(std::ofstream &out,
                                     const SamplePairOverlayed &pair,
                                     const vChrArray &vReadArray,
                                     const int32_t posi, const std::string &strand)
{
  int32_t bincenter(posi/binsize);
  int32_t sbin(bincenter - binwidth_from_center);
  int32_t ebin(bincenter + binwidth_from_center);

  if (strand == "+") {
    for (int32_t i=sbin; i<=ebin; ++i) out << "\t" << getReadVal(vReadArray, pair, i, stype);
  } else {
    for (int32_t i=ebin; i>=sbin; --i) out << "\t" << getReadVal(vReadArray, pair, i, stype);
  }
}

double ReadProfile::getAverageVal(const SamplePairOverlayed &pair,
                                  const vChrArray &vReadArray,
                                  const int32_t sbin,
                                  const int32_t ebin)
{
  double sumIP(0);
  for (int32_t i=sbin; i<=ebin; ++i) sumIP += getReadVal(vReadArray, pair, i, stype);
  return getratio(sumIP, (ebin - sbin + 1));
}

double ReadProfile::getMaxVal(const SamplePairOverlayed &pair,
                              const vChrArray &vReadArray,
                              const int32_t sbin,
                              const int32_t ebin)
{
  double maxIP(0);
  for (int32_t i=sbin; i<=ebin; ++i) {
    maxIP = std::max(maxIP, getReadVal(vReadArray, pair, i, stype));
  }
  return maxIP;
}

ReadProfile::ReadProfile(const DROMPA::Global &p, const int32_t _nbin):
  stype(p.prof.is_distribution_type()),
  width_from_center(p.prof.get_width_from_center()),
  nsites(0), nsites_skipped(0)
{

  if(p.prof.isPtypeTSS())          xlabel = "Distance from TSS (bp)";
  else if(p.prof.isPtypeTTS())     xlabel = "Distance from TES (bp)";
  else if(p.prof.isPtypeGene100()) xlabel = "% gene length from TSS";
  else if(p.prof.isPtypeBed())     xlabel = "Distance from the peak summit (bp)";

  for (auto &x: p.samplepair) binsize = x.first.getbinsize();
  binwidth_from_center = width_from_center / binsize;
  if (binwidth_from_center <= 1) {
    PRINTERR_AND_EXIT("please specify larger size for --widthfromcenter:"
		      << width_from_center
		      << " than binsize:"
		      << binsize << ".");
  }

  if (_nbin) nbin = _nbin;
  else nbin = binwidth_from_center * 2 +1;

  std::vector<double> array(nbin, 0);
  for (auto &x: p.samplepair) hprofile[x.first.argvChIP] = array;
}

void ReadProfile::setOutputFilename(const DROMPA::Global &p, const std::string &commandname)
{
  std::string prefix(p.getPrefixName() + "." + commandname);

  if (p.isgetmaxval()) prefix += ".maxvalue";
  else                 prefix += ".averaged";
  if      (stype==0)   prefix += ".ChIPread";
  else if (stype==1)   prefix += ".Enrichment";

  Rscriptname = prefix + ".R";
  RDataname   = prefix;
  Rfigurename = prefix + ".pdf";
  for (auto &x: {Rscriptname, RDataname, Rfigurename}) unlink(x.c_str());
  for (auto &x: p.samplepair) {
    std::string file(RDataname + "." + x.first.label + ".tsv");
    //      std::cout << file << std::endl;
    unlink(file.c_str());
  }
  for (auto &x: {Rscriptname, Rfigurename}) unlink(x.c_str());
}

void ReadProfile::MakeFigure(const DROMPA::Global &p)
{
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

  out << "ymax <- ceiling(max(c(";
  for (size_t i=1; i<=p.samplepair.size(); ++i) {
    out << boost::format("max(p%1%)") % i;
    if(i<p.samplepair.size()) out << ",";
    else out << ")))\n";
  }
  out << "ymin <- ceiling(min(c(";
  for (size_t i=1; i<=p.samplepair.size(); ++i) {
    out << boost::format("min(p%1%)") % i;
    if(i<p.samplepair.size()) out << ",";
    else out << ")))\n";
  }

  out << "pdf('" << Rfigurename << "',6,6)" << std::endl;
  out << "plot(x,p1,type='l',col=rgb(" << vcol[0] << "," << vcol[1] << "," << vcol[2] << ")";

/*  double min(1);
  for (size_t i=1; i<=p.samplepair.size(); ++i) {
    for(j=0; j<num; j++){
      if(min > sample[i].profile.IP[j]) min = sample[i].profile.IP[j];
    }
  }
  min = (int)(min*10)/10.0;
  if(!min) min = 0.1;*/

  if (p.prof.isPtypeGene100()) out << ", log='y', ylim=c(ymin, ymax)";
  else out << ", ylim=c(0, ymax)";

  if (!stype) {
    out << ",xlab='" << xlabel << "',ylab='Read density')" << std::endl;
  } else if (stype==1) {
    out << ",xlab='" << xlabel << "',ylab='ChIP/Input enrichment')" << std::endl;
  }
  out << "polygon(c(x, rev(x)), c(p1_lower, rev(p1_upper)), col=rgb("
      << vcol[0] << "," << vcol[1] << "," << vcol[2] << ",0.3), border=NA)" << std::endl;

  for (size_t i=1; i<p.samplepair.size(); ++i) {
    std::string colstr = std::to_string(vcol[i*3]) + ","
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


void ProfileTSS::WriteTSV_EachChr(const DROMPA::Global &p, const chrsize &chr)
{
  DEBUGprint_FUNCStart();

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

  DEBUGprint_FUNCend();
}

void ProfileGene100::outputEachGene_fixedlength(std::ofstream &out, const SamplePairOverlayed &x,
                                                const genedata &gene, const vChrArray &vReadArray, int32_t len,
                                                const int32_t width_from_gene)
{
  double len100(len / (double)GENEBLOCKNUM);
  double div100(width_from_gene/ (double)GENEBLOCKNUM);

//  printf("genestart %d, geneend %d, strand %s, GENEBLOCKNUM %d, width_from_gene %d, len100 %1f, div100 %1f\n",
//         gene.txStart, gene.txEnd, gene.strand.c_str(), GENEBLOCKNUM, width_from_gene, len100, div100);

  for (int32_t i=0; i<nbin; ++i) {
    int32_t s(0), e(0);
    if (i < GENEBLOCKNUM) {
      if (gene.strand == "+") {
        s = (gene.txStart - width_from_gene + div100 *i)         / binsize;
        e = (gene.txStart - width_from_gene + div100 *(i+1) -1)  / binsize;
      } else {
        s = (gene.txEnd   + width_from_gene - div100 * (i+1)) / binsize;
        e = (gene.txEnd   + width_from_gene - div100 * i -1)  / binsize;
      }
    } else if (i < GENEBLOCKNUM*2) {
      if (gene.strand == "+") {
        s = (gene.txStart + len100 * (i - GENEBLOCKNUM))      / binsize;
        e = (gene.txStart + len100 * (i+1 - GENEBLOCKNUM) -1) / binsize;
      } else {
        s = (gene.txEnd   - len100 * (i+1 - GENEBLOCKNUM))  / binsize;
        e = (gene.txEnd   - len100 * (i - GENEBLOCKNUM) -1) / binsize;
      }
    } else {
      if (gene.strand == "+") {
        s = (gene.txEnd   + div100 * (i - 2*GENEBLOCKNUM))      / binsize;
        e = (gene.txEnd   + div100 * (i+1 - 2*GENEBLOCKNUM) -1) / binsize;
      } else {
        s = (gene.txStart - div100 * (i+1 - 2*GENEBLOCKNUM))  / binsize;
        e = (gene.txStart - div100 * (i - 2*GENEBLOCKNUM) -1) / binsize;
      }
    }
  //  printf("i %d, s %d  e %d\n", i,s,e);

    out << "\t" << getAverageVal(x, vReadArray, s, e);
  }
//  if (gene.strand == "-") exit(0);

}


void ProfileGene100::outputEachGene(std::ofstream &out, const SamplePairOverlayed &x,
                                    const genedata &gene, const vChrArray &vReadArray, int32_t len)
{
  double len100(len / (double)GENEBLOCKNUM);

  for (int32_t i=0; i<nbin; ++i) {
    int32_t s(0), e(0);
    if (gene.strand == "+") {
      s = (gene.txStart - len + len100 *i)       / binsize;
      e = (gene.txStart - len + len100 *(i+1) -1)/ binsize;
    }else{
      s = (gene.txEnd + len - len100 * (i+1))/ binsize;
      e = (gene.txEnd + len - len100 * i -1) / binsize;
    }
    out << "\t" << getAverageVal(x, vReadArray, s, e);
  }
}

void ProfileGene100::WriteTSV_EachChr(const DROMPA::Global &p, const chrsize &chr)
{
  DEBUGprint_FUNCStart();

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

      if (p.prof.isusefixedlength()) outputEachGene_fixedlength(out, x, gene, vReadArray, len, p.prof.get_width_from_center());
      else outputEachGene(out, x, gene, vReadArray, len);
      out << std::endl;
    }
  }

  DEBUGprint_FUNCend();
}

void ProfileBedSites::WriteTSV_EachChr(const DROMPA::Global &p, const chrsize &chr)
{
  DEBUGprint_FUNCStart();

  if(!p.anno.vbedlist.size()) PRINTERR_AND_EXIT("Please specify --bed.");

  vChrArray vReadArray(p, chr);

  for (auto &x: p.samplepair) {
    std::string file(RDataname + "." + x.first.label + ".tsv");
    std::ofstream out(file, std::ios::app);

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

  DEBUGprint_FUNCend();
}

void ProfileMULTICI::WriteTSV_EachChr(const DROMPA::Global &p, const chrsize &chr)
{
  DEBUGprint_FUNCStart();

  if(!p.anno.vbedlist.size()) PRINTERR_AND_EXIT("Please specify --bed.");

  vChrArray vReadArray(p, chr);

  std::string file(RDataname + ".tsv");
  std::ofstream out(file, std::ios::app);

  for (auto &vbed: p.anno.vbedlist) {
    for (auto &bed: vbed.getvBed()) {
      if (bed.chr != rmchr(chr.getname())) continue;
      ++nsites;

      if (bed.start < 0 || bed.end >= chr.getlen()) {
        ++nsites_skipped;
        continue;
      }

      int32_t sbin(bed.start/binsize);
      int32_t ebin((bed.end-1)/binsize);

      if (p.isaddname()) out << bed.getSiteStrTABwithNAME();
      else out << bed.getSiteStrTAB();

      for (auto &x: p.samplepair) {
        if (p.isgetmaxval()) out << "\t" << getMaxVal(x, vReadArray, sbin, ebin);
        else                 out << "\t" << getAverageVal(x, vReadArray, sbin, ebin);
      }
      out << std::endl;
    }
  }
  out.close();

  DEBUGprint_FUNCend();
}
