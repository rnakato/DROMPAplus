/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_DRAW_DATAFRAME_H_
#define _DD_DRAW_DATAFRAME_H_

#include "dd_draw_pdfpage.hpp"
#include "color.hpp"

class DataFrame {
  virtual void StrokeYmem(const int32_t nlayer);

  void StrokeFrame() {
    cr->set_line_width(0.4);
    cr->set_source_rgba(CLR_BLACK, 1);
    if(bothdirection || shownegative) rel_xline(cr, OFFSET_X, par.yaxis_now - len_minus, par.getXaxisLen());
    else rel_xline(cr, OFFSET_X, par.yaxis_now - height_df, par.getXaxisLen());
    rel_yline(cr, OFFSET_X, par.yaxis_now - height_df, height_df);
    cr->stroke();
  }

  void StrokeBins(const SamplePairEach &pair, const vChrArray &vReadArray, const int32_t nlayer);
  virtual void StrokeEachBin(const SamplePairEach &pair, const vChrArray &vReadArray,
			     const int32_t i, const double xcen, const int32_t yaxis,
			     const int32_t nlayer);

protected:
  enum { POSI_ASSAYLABEL=6, POSI_XLABEL=66, LEN_EDGE=2};
  const Cairo::RefPtr<Cairo::Context> cr;
  const DParam &par;
  double scale, scale2nd;
  std::string label, label2nd;
  double width_df;
  double height_df;

  int32_t barnum_minus;
  int32_t barnum_plus;
  double len_minus, len_plus;

  bool sigtest;
  double threshold;
  std::string chrname;
  bool shownegative;

  bool bothdirection;
  uint32_t ndigit;

  double getEthre(const DROMPA::Global &p) const {
    if (p.isGV) return p.drawparam.scale_ratio;
    else return p.thre.ethre;
  }

 public:
  DataFrame(const Cairo::RefPtr<Cairo::Context> cr_,
	    const std::string &l, const std::string &l2,
	    const double sglobal, const double slocal1, const double slocal2,
	    const DParam &refparam, const bool sig, const double thre,
	    const std::string &_chrname, const int32_t width_draw):
    cr(cr_), par(refparam), label(l), label2nd(l2),
    width_df(width_draw), height_df(refparam.get_height_df()),
    barnum_minus(refparam.barnum/2),
    barnum_plus(refparam.barnum - barnum_minus),
    len_minus(barnum_minus*par.ystep),
    len_plus(barnum_plus*par.ystep),
    sigtest(sig), threshold(thre), chrname(_chrname),
    shownegative(false), bothdirection(false), ndigit(1)
  {
    if (slocal1) scale    = slocal1; else scale    = sglobal;
    if (slocal2) scale2nd = slocal2; else scale2nd = sglobal;
  }

  void Draw(const DROMPA::Global &p, const SamplePairOverlayed &pair, const vChrArray &vReadArray);

  virtual void HighlightPeaks(const SamplePairOverlayed &pair){ (void)(pair); return; }

  virtual void getColor1st(const double alpha)=0;
  virtual void getColor2nd(const double alpha)=0;
  virtual void setColor(const double value, const int32_t nlayer, const double alpha);
  virtual double getVal(const SamplePairEach &pair, const vChrArray &vReadArray, const int32_t i)=0;
  virtual const std::string getAssayName() const =0;
  virtual double get_yscale_num(int32_t i, double scale) const {
    if (shownegative) return (i - barnum_minus) * scale;
    else return i*scale;
  }

  void StrokeSampleLabel(const SamplePairOverlayed &pair);
};


class ChIPDataFrame : public DataFrame {
  double getVal(const SamplePairEach &pair, const vChrArray &vReadArray, const int32_t i) {
    return vReadArray.getArray(pair.argvChIP).array[i];
  }

  void getColor1st(const double alpha) { cr->set_source_rgba(CLR_GREEN3, alpha); }
  void getColor2nd(const double alpha) { cr->set_source_rgba(CLR_ORANGE, alpha); }
  const std::string getAssayName() const { return "Reads"; }

 public:
  ChIPDataFrame(const Cairo::RefPtr<Cairo::Context> cr_,
		const DROMPA::Global &p,
		const SamplePairOverlayed &pair,
		const DParam &refparam,
		const std::string &chrname):
    DataFrame(cr_, pair.first.label, pair.second.label,
	      p.drawparam.scale_tag, pair.first.scale.tag, pair.second.scale.tag,
	      refparam, p.thre.sigtest, p.thre.ipm, chrname,
	      p.drawparam.width_draw_pixel)
  {
    shownegative = p.drawparam.isshownegative();
  }

  template <class T>
  void strokePeaks(T &bed) {
    if (!my_overlap(bed.start, bed.end, par.xstart, par.xend)) return;

    int32_t s(std::max(bed.start, par.xstart) - par.xstart);
    int32_t e(std::min(bed.end, par.xend)     - par.xstart);
    double x(BP2PIXEL(s));
    double len((e-s)* par.dot_per_bp);
/*    DEBUGprint("bed.start " << bed.start << " bed.end "  << bed.end
	       << " bed.len "  << bed.length()
	       << " par.xstart " << par.xstart << " par.xend "  << par.xend
	       << " s " << s << " e " << e << " x " <<x << " len " << len);*/
    rel_xline(cr, x, par.yaxis_now - height_df/2, len);
    cr->stroke();
  }

  void HighlightPeaks(const SamplePairOverlayed &pair) {
    cr->set_source_rgba(CLR_ORANGE, 0.4);
    cr->set_line_width(height_df + 6);

    if (pair.first.BedExists()) { // specified BED
      for (auto &peak: pair.first.getBedChr(chrname)) {
	strokePeaks<bed>(peak);
      }
    } else {     // peak calling by DROMPA+
      for (auto &peak: pair.first.getPeakChr(chrname)) {
//	peak.print();
	strokePeaks<Peak>(peak);
      }
    }

    cr->stroke();
  }
};

class InputDataFrame : public DataFrame {
  double getVal(const SamplePairEach &pair, const vChrArray &vReadArray, const int32_t i) {
    return vReadArray.getArray(pair.argvInput).array[i];
  }

  void getColor1st(const double alpha) { cr->set_source_rgba(CLR_BLUE, alpha); }
  void getColor2nd(const double alpha) { cr->set_source_rgba(CLR_OLIVE, alpha); }
  const std::string getAssayName() const { return "Input"; }

 public:
  InputDataFrame(const Cairo::RefPtr<Cairo::Context> cr_, const DROMPA::Global &p,
		 const SamplePairOverlayed &pair, const DParam &refparam,
		 const std::string &chrname):
    DataFrame(cr_, pair.first.label, pair.second.label,
	      p.drawparam.scale_tag, pair.first.scale.tag, pair.second.scale.tag,
	      refparam, false, 0, chrname,
	      p.drawparam.width_draw_pixel)
  {
    (void)(pair);
    shownegative = p.drawparam.isshownegative();
  }
};


class RatioDataFrame : public DataFrame {
  bool isGV;

  void setColor(const double value, const int32_t nlayer, const double alpha);

  double getVal(const SamplePairEach &pair, const vChrArray &vReadArray, const int32_t i)
  {
    return CalcRatio(vReadArray.getArray(pair.argvChIP).array[i],
		     vReadArray.getArray(pair.argvInput).array[i],
		     pair.ratio);
  }

  void getColor1st(const double alpha) { cr->set_source_rgba(CLR_ORANGE, alpha); }
  void getColor2nd(const double alpha) { cr->set_source_rgba(CLR_DEEPSKYBLUE, alpha); }
  const std::string getAssayName() const { return "ChIP/Input"; }

public:
  RatioDataFrame(const Cairo::RefPtr<Cairo::Context> cr_, const DROMPA::Global &p,
		 const SamplePairOverlayed &pair, const DParam &refparam,
		 const std::string &chrname):
    DataFrame(cr_, pair.first.label, pair.second.label,
	      p.drawparam.scale_ratio, pair.first.scale.ratio, pair.second.scale.ratio,
	      refparam, p.thre.sigtest, getEthre(p), chrname,
	      p.drawparam.width_draw_pixel),
    isGV(p.isGV)
  {
    shownegative = p.drawparam.isshownegative();
  }
};

class LogRatioDataFrame : public DataFrame {
  bool isGV;

  void setColor(const double value, const int32_t nlayer, const double alpha);

  void getColor1st(const double alpha) { cr->set_source_rgba(CLR_ORANGE, alpha); }
  void getColor2nd(const double alpha) { cr->set_source_rgba(CLR_DEEPSKYBLUE, alpha); }
  const std::string getAssayName() const { return "Enrichment"; }

  double getVal(const SamplePairEach &pair, const vChrArray &vReadArray, const int32_t i)
  {
    return CalcRatio(vReadArray.getArray(pair.argvChIP).array[i],
		     vReadArray.getArray(pair.argvInput).array[i],
		     pair.ratio);
  }

  double get_yscale_num(int32_t i, double scale) const {
    return std::pow(scale, i - barnum_minus);
  }

  void StrokeEachBin(const SamplePairEach &pair, const vChrArray &vReadArray,
		     const int32_t i, const double xcen, const int32_t yaxis,
		     const int32_t nlayer);

 public:
  LogRatioDataFrame(const Cairo::RefPtr<Cairo::Context> cr_, const DROMPA::Global &p,
		    const SamplePairOverlayed &pair, const DParam &refparam,
		    const std::string &chrname):
    DataFrame(cr_, pair.first.label, pair.second.label,
	      p.drawparam.scale_ratio, pair.first.scale.ratio, pair.second.scale.ratio,
	      refparam, p.thre.sigtest, getEthre(p), chrname,
	      p.drawparam.width_draw_pixel),
    isGV(p.isGV)
  {
    bothdirection = true;
    ndigit = 2;
  }
};

class PvalueDataFrame : public DataFrame {
  void setColor(const double value, const int32_t nlayer, const double alpha);
  virtual double getVal(const SamplePairEach &pair, const vChrArray &vReadArray, const int32_t i)=0;
  void getColor1st(const double alpha) { cr->set_source_rgba(CLR_RED, alpha); }
  void getColor2nd(const double alpha) { cr->set_source_rgba(CLR_BLUE, alpha); }
  virtual const std::string getAssayName() const =0;

public:
  PvalueDataFrame(const Cairo::RefPtr<Cairo::Context> cr_, const DROMPA::Global &p,
		  const SamplePairOverlayed &pair, const DParam &refparam,
		  const std::string &chrname,
		  double pvalue):
    DataFrame(cr_,
	      pair.first.label, pair.second.label,
	      p.drawparam.scale_pvalue,
	      pair.first.scale.pvalue,
	      pair.second.scale.pvalue,
	      refparam, p.thre.sigtest,
	      pvalue,
	      chrname,
	      p.drawparam.width_draw_pixel)
  {}
};

class PinterDataFrame : public PvalueDataFrame {
  double getVal(const SamplePairEach &pair, const vChrArray &vReadArray, const int32_t i) {
    const ChrArray &a = vReadArray.getArray(pair.argvChIP);
    double myu(a.array.getLocalAverage(i, a.binsize));

    return getlogp_Poisson(a.array[i], myu);
  }
  const std::string getAssayName() const { return "logp(ChIP)"; }

public:
  PinterDataFrame(const Cairo::RefPtr<Cairo::Context> cr_, const DROMPA::Global &p,
		  const SamplePairOverlayed &pair, const DParam &refparam,
		  const std::string &chrname):
    PvalueDataFrame(cr_, p, pair, refparam, chrname, p.thre.pthre_inter)
  {}
};

class PenrichDataFrame : public PvalueDataFrame {
  double getVal(const SamplePairEach &pair, const vChrArray &vReadArray, const int32_t i) {
    return getlogp_BinomialTest(vReadArray.getArray(pair.argvChIP).array[i],
				vReadArray.getArray(pair.argvInput).array[i],
				pair.ratio);
  }
  const std::string getAssayName() const { return "logp(Enrich)"; }

 public:
  PenrichDataFrame(const Cairo::RefPtr<Cairo::Context> cr_,
		   const DROMPA::Global &p, const SamplePairOverlayed &pair,
		   const DParam &refparam, const std::string &chrname):
    PvalueDataFrame(cr_, p, pair, refparam, chrname, p.thre.pthre_enrich)
  {}
};

#endif /* _DD_DRAW_DATAFRAME_H_ */
