/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "FragmentClusterScore.hpp"
#include "Mapfile.hpp"
#include "FragmentClusterScore_p.hpp"

namespace {
  std::vector<int8_t> genVector4FixedReadsNum(const SeqStats &chr, const double r4cmp, int32_t &numUsed4FCS, Strand::Strand strand)
  {
    int32_t chrlen(chr.getlen());
    std::vector<int8_t> array(chrlen, 0);
    for (auto &x: chr.getvReadref(strand)) {
      if(!x.duplicate && my_range(x.F3, 0, chrlen-1)){
	if(rand() >= r4cmp) continue;
	++array[x.F3];
	++numUsed4FCS;
      }
    }
    return array;
  }

  void makeRscript(const std::string prefix)
  {
    std::string Rscript(prefix + ".FCS.R");
    std::ofstream out(Rscript);
    out << "data <- read.csv('" << prefix << ".pnf.csv', header=TRUE, row.names=1, sep='\\t', quote='')" << std::endl;
    out << "colnames(data) <- colnames(data)[-1]" << std::endl;
    out << "data <- data[,-ncol(data)]" << std::endl;
    out << "nrow <- nrow(data)" << std::endl;
    out << "array <- c('50','150', '500', '1000', '5000', '10000', '100000', '1000000')" << std::endl;

    out << "ncol <- length(array)" << std::endl;
    out << "cols <- rainbow(ncol)" << std::endl;
    out << "cpnf <- data[,paste('CPNF.len', array, sep='')]" << std::endl;
    out << "colnames(cpnf) <- paste('len', array, sep='')" << std::endl;
    out << "x <- seq(1, nrow, 10)" << std::endl;
    out << "cpnf10 <- cpnf[x,]" << std::endl;
    out << "pnf <- rbind(cpnf10[-1,],cpnf[nrow,]) - cpnf10" << std::endl;
    out << std::endl;
    out << "pdf('" << prefix << ".FCS.pdf', height=5, width=15)" << std::endl;
    out << "par(mfrow=c(1,3))" << std::endl;
    // Proportion of NN fragments
    out << "plot(0, 0, type = 'n', xlim = range(1:nrow), ylim = c(0,max(c(max(pnf),0.2))), xlab = 'Neighboring distance (bp)', ylab = 'Proportion of nearest neibor fragments')" << std::endl;
    out << "for (i in 1:ncol) { lines(x, pnf[,i], col=cols[i])}" << std::endl;
    out << "legend('topright', legend = colnames(pnf), lty = 1, col = cols)" << std::endl;
    // Cumurative proportion
    out << "plot(0, 0, type = 'n', xlim = range(1:nrow), ylim = c(0,max(c(max(cpnf),0.4))), xlab = 'Neighboring distance (bp)', ylab = 'Cumulative proportion')" << std::endl;
    out << "for (i in 1:ncol) { lines(1:nrow, cpnf[,i], col=cols[i])}" << std::endl;
    out << "legend('topright', legend = colnames(cpnf), lty = 1, col = cols)" << std::endl;
    // FCS
    out << "data <- read.csv('" << prefix << ".fcs.csv', header=TRUE, skip=4, sep='\t', quote='')" << std::endl;
    out << "plot(data[,1],data[,2], log='x', type='l', ylim = c(0,max(c(data[,2],0.2))), xlab = 'Read-pair distance (bp)', ylab = 'Fragment cluster score')" << std::endl;
    ///    out << "data <- read.csv('" << prefix << ".KLD.csv', header=TRUE, skip=4, sep='\t', quote='')" << std::endl;
    //    out << "plot(data[,1],data[,2], log='x', type='l', ylim = c(0,max(c(data[,2],0.2))), xlab = 'Read-pair distance (bp)', ylab = 'Kullback-Leibler divergence')" << std::endl;
    out << "dev.off()" << std::endl;

    std::string command = "R --vanilla < " + Rscript + " > " + Rscript + ".log 2>&1";
    std::cout << command << std::endl;
    int32_t return_code = system(command.c_str());
    if(WEXITSTATUS(return_code)) {
      std::cerr << "Warning: command " << command << "return nonzero status." << std::endl;
    }
    return;
  }
}

void makeFCSProfile(FCSstats &fcsst, const SeqStatsGenome &genome, const std::string &head, const std::string &typestr)
{
  shiftFragVar dist(fcsst, genome);
  std::cout << "Making FCS profile..." << std::flush;

  for(size_t i=0; i<genome.chr.size(); ++i) {
    if(!genome.chr[i].isautosome()) continue;
    std::cout << genome.chr[i].getname() << ".." << std::flush;
    dist.execchr(genome.chr[i]);
  }
  std::cout << "\nread number for calculating FCS: " << dist.getnumUsed4FCS() << std::endl;

  dist.calcFCS();
  
  std::string filename1 = head + ".pnf.csv";
  dist.outputPnf(filename1);
  std::string filename2 = head + "." + typestr + ".csv";
  dist.outputFCS(filename2);
  //  filename2 = head + ".KLD.csv";
  // dist.outputKLD(filename2);

  makeRscript(head);
  
  dist.setFCSstats(fcsst);

  return;
}

shiftFragVar::shiftFragVar(const FCSstats &fcsst, const SeqStatsGenome &genome):
   lenF3(genome.dflen.getlenF3()), flen(genome.dflen.getflen()),
   r4cmp(0), numUsed4FCS(0), lackOfReads(false),
   ng_from_fcs(fcsst.getNgFromFCS()),
   ng_to_fcs(fcsst.getNgToFCS()),
   ng_step_fcs(fcsst.getNgStepFCS()),
   nread(0)
  {
    std::vector<int32_t> v{lenF3, flen};
    std::copy(v4pnf.begin(), v4pnf.end(), std::back_inserter(v));
    for(auto len: v) pnf[len] = PropNeighborFrag();
    
    for(auto &x: genome.chr) {
      if(x.isautosome()) nread += x.getnread_nonred(Strand::BOTH);
    }
	
    //double r = (getratio(sspst.getnum4ssp(), nread)) / (NUM_100M/static_cast<double>(dist.getlen()));
    double r = getratio(fcsst.getnum4fcs(), nread);
#ifdef DEBUG
    std::cout << "\nr for FCS\t" << r << "\t reads: " << nread << std::endl;
#endif
    if(r>1){
      std::cerr << "\nWarning: number of reads (" << nread << ") is less than num4ssp ("<< fcsst.getnum4fcs() <<").\n";
      lackOfReads=true;
    }
    r4cmp = r*RAND_MAX;
  }

void PropNeighborFrag::setNeighborFrag(const int32_t flen, const int32_t end,
		      const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev)
{
  if(flen < 0) {
    std::cerr << "error: invalid flen " << flen << "for setNeighborFrag." << std::endl;
    return;
  }
  
  int32_t distance(0);
  for(int32_t i=0; i<end-flen; ++i) {
    if(fwd[i] && rev[i+flen]) {
      if(distance < sizeOfvNeighborFrag-1) ++vNeighborFrag[distance];
      else ++vNeighborFrag[sizeOfvNeighborFrag-1];
      distance = 0;
    }
    ++distance;
  }

  sumOfvNeighborFrag = 0;
  for(auto x: vNeighborFrag) sumOfvNeighborFrag += x;

  
  for(size_t i=0; i<sizeOfvNeighborFrag; ++i) {
    if(!i) vcPNF[i] = getPNF(0);
    else   vcPNF[i] = vcPNF[i-1] + getPNF(i);
  }
}

void shiftFragVar::execchr(const SeqStats &chr)
{
  auto fwd = genVector4FixedReadsNum(chr, r4cmp, numUsed4FCS, Strand::FWD);
  auto rev = genVector4FixedReadsNum(chr, r4cmp, numUsed4FCS, Strand::REV);
  
  for(int32_t flen=ng_from_fcs; flen<ng_to_fcs; flen += ng_step_fcs) {
    pnfbg[flen].setNeighborFrag(flen, chr.getlen(), fwd, rev);
  }
  for(auto x: pnf) {
    pnf[x.first].setNeighborFrag(x.first, chr.getlen(), fwd, rev);
  }
}

