/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "dd_draw_dataframe.hpp"

void DataFrame::StrokeYmem(const int32_t nlayer) {
  cr->set_source_rgba(CLR_BLACK, 0.5);
  for (int32_t i=0; i<par.barnum; ++i) rel_xline(cr, OFFSET_X, par.yaxis_now - i*par.ystep, par.getXaxisLen());
  cr->stroke();

  cr->set_source_rgba(CLR_BLACK, 1);
  double x(0);
  if (!nlayer) x = OFFSET_X + par.getXaxisLen() + 7;
  else x = OFFSET_X - 20;

  for(int32_t i=1; i<=par.barnum; ++i) {
    std::string y_value;
    if (!nlayer) y_value = float2string(i*scale, 1);
    else         y_value = float2string(i*scale2nd, 1);
    showtext_cr(cr, x, par.yaxis_now - i*(par.ystep - 1.5), y_value, 9);
  }
}

void DataFrame::StrokeBins(const SamplePairEach &pair, const vChrArray &vReadArray, const int32_t nlayer) {
  int32_t binsize(pair.getbinsize());
  int32_t sbin(par.xstart/binsize);
  int32_t ebin(par.xend/binsize);
  double dot_per_bin(binsize * par.dot_per_bp);
  int32_t yaxis(par.yaxis_now);  // convert to int

  int32_t thin(std::min(par.width_per_line/(1000*binsize), 20));

  double xcen(BP2PIXEL(binsize/2)); // initial position
  if (thin > 1) cr->set_line_width(dot_per_bin*thin);
  else cr->set_line_width(dot_per_bin);

  for (int32_t i=sbin; i<ebin; ++i, xcen += dot_per_bin) {
    if (thin > 1 && i%thin) continue;
    StrokeEachBin(pair, vReadArray, i, xcen, yaxis, nlayer);
  }
  cr->stroke();
}

void DataFrame::StrokeEachBin(const SamplePairEach &pair,
			      const vChrArray &vReadArray,
			      const int32_t i, const double xcen,
			      const int32_t yaxis, const int32_t nlayer) {
  double value(getVal(pair, vReadArray, i));
  if (!value) return;

  int32_t len;
  if (!nlayer) len = getbinlen(value / scale);
  else len = getbinlen(value / scale2nd);

  // waku
  setColor(value, nlayer, 1);
  if (len <0) rel_yline(cr, xcen, yaxis + len, len_binedge);
  cr->stroke();

  // nakami
  if (par.alpha) {
    setColor(value, nlayer, par.alpha);
    rel_yline(cr, xcen, yaxis, len);
    cr->stroke();
  }
}

void DataFrame::Draw(const DROMPA::Global &p, const SamplePairOverlayed &pair, const vChrArray &vReadArray) {
  int32_t nlayer(0);
  StrokeBins(pair.first, vReadArray, nlayer);
  if (p.drawparam.isshowymem()) StrokeYmem(nlayer);

  StrokeFrame();
  HighlightPeaks(pair);

  // Overlayed
  if (pair.OverlayExists()) {
    nlayer = 1;
    StrokeBins(pair.second, vReadArray, nlayer);
    if (p.drawparam.isshowymem()) StrokeYmem(nlayer);
  }

  // Sample Label
  if (p.drawparam.isshowylab()) StrokeYlab(pair);
}

void RatioDataFrame::setColor(const double value, const int32_t nlayer, const double alpha)
{
  if (!nlayer) { // first layer
    if(isGV || sigtest) {
      if (value > threshold) cr->set_source_rgba(CLR_RED, alpha);
      else cr->set_source_rgba(CLR_GRAY, alpha);
    } else {
      cr->set_source_rgba(CLR_ORANGE, alpha);
    }
  } else {    // second layer
    if(isGV || sigtest) {
      if (value > threshold) cr->set_source_rgba(CLR_DARKORANGE, alpha);
      else cr->set_source_rgba(CLR_GRAY2, alpha);
    } else {
      cr->set_source_rgba(CLR_DEEPSKYBLUE, alpha);
    }
  }
}

void LogRatioDataFrame::setColor(const double value, const int32_t nlayer, const double alpha){
  if (!nlayer) {
    if (isGV || sigtest) {
      if (value > 0) cr->set_source_rgba(CLR_SALMON, 1);
      else cr->set_source_rgba(CLR_BLUEPURPLE, 1);
    } else cr->set_source_rgba(CLR_ORANGE, 1);
  } else {
    if (isGV || sigtest) {
      if (value > 0) cr->set_source_rgba(CLR_DEEPSKYBLUE, alpha);
      else cr->set_source_rgba(CLR_GRAY2, alpha);
    } else cr->set_source_rgba(CLR_PURPLE, alpha);
  }
}

void LogRatioDataFrame::StrokeEachBin(const SamplePairEach &pair,
				      const vChrArray &vReadArray,
				      const int32_t i,
				      const double xcen,
				      const int32_t yaxis,
				      const int32_t nlayer) {
  double value(CalcRatio(vReadArray.getArray(pair.argvChIP).array[i],
			 vReadArray.getArray(pair.argvInput).array[i],
			 pair.ratio));
  if (!value) return;

  if (!nlayer) value = value ? log(value) / log(scale): 0;
  else         value = value ? log(value) / log(scale2nd): 0;

  int32_t len(0);
  if (value > 0)      len = -std::min(par.ystep*value, len_plus);
  else if (value < 0) len = -std::max(par.ystep*value, -len_minus);

  // waku
  setColor(value, nlayer, par.alpha);
  rel_yline(cr, xcen, yaxis - len_minus + len, len_binedge);
  cr->stroke();
  // nakami
  if (par.alpha) {
    setColor(value, nlayer, par.alpha);
    rel_yline(cr, xcen, yaxis - len_minus, len);
    cr->stroke();
  }
}

// StrokeYlab

void ChIPDataFrame::StrokeYlab(const SamplePairOverlayed &pair) {
  if (pair.OverlayExists()) {
    cr->set_source_rgba(CLR_GREEN3, 1);
    showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2 - 7, label, 12);
    cr->set_source_rgba(CLR_ORANGE, 1);
    showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2 + 7, label2nd, 12);
  } else {
    cr->set_source_rgba(CLR_BLACK, 1);
    showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2, label, 12);
  }
}

void InputDataFrame::StrokeYlab(const SamplePairOverlayed &pair)
{
  if (pair.OverlayExists()) {
    cr->set_source_rgba(CLR_BLACK, 1);
    showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2 - 7, label, 12);
    showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2 + 7, "1: blue, 2: olive", 12);
  } else {
    cr->set_source_rgba(CLR_BLACK, 1);
    showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2, label, 12);
  }
}

void RatioDataFrame::StrokeYlab(const SamplePairOverlayed &pair)
{
  if (pair.OverlayExists()) {
    if (label2nd != "") {
      cr->set_source_rgba(CLR_ORANGE, 1);
      showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2 - 7, label, 12);
      cr->set_source_rgba(CLR_DEEPSKYBLUE, 1);
      showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2 + 7, label2nd, 12);
    } else {
      cr->set_source_rgba(CLR_BLACK, 1);
      showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2 - 7, label, 12);
      showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2 + 7, "1: orange, 2: purple", 12);
    }
  } else {
    cr->set_source_rgba(CLR_BLACK, 1);
    showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2, label, 12);
  }
}

void LogRatioDataFrame::StrokeYlab(const SamplePairOverlayed &pair) {
  if (pair.OverlayExists()) {
    if (label2nd != "") {
      cr->set_source_rgba(CLR_ORANGE, 1);
      showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2 - 7, label, 12);
      cr->set_source_rgba(CLR_PURPLE, 1);
      showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2 + 7, label2nd, 12);
    } else {
      cr->set_source_rgba(CLR_BLACK, 1);
      showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2 - 7, label, 12);
      showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2 + 7, "1: orange, 2: purple", 12);
    }
  } else {
    cr->set_source_rgba(CLR_BLACK, 1);
    showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2, label, 12);
  }
}

void PinterDataFrame::StrokeYlab(const SamplePairOverlayed &pair)
{
  if (pair.OverlayExists()) {
    if (label2nd != "") {
      cr->set_source_rgba(CLR_RED, 1);
      showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2 - 7, label, 12);
      cr->set_source_rgba(CLR_YELLOW2, 1);
      showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2 + 7, label2nd, 12);
    } else {
      cr->set_source_rgba(CLR_BLACK, 1);
      showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2 - 7, label, 12);
      showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2 + 7, "1: red, 2: yellow", 12);
    }
  } else {
    cr->set_source_rgba(CLR_BLACK, 1);
    showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2, label, 12);
  }
}

void PenrichDataFrame::StrokeYlab(const SamplePairOverlayed &pair)
{
  if (pair.OverlayExists()) {
    if (label2nd != "") {
      cr->set_source_rgba(CLR_RED, 1);
      showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2 - 7, label, 12);
      cr->set_source_rgba(CLR_YELLOW2, 1);
      showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2 + 7, label2nd, 12);
    } else {
      cr->set_source_rgba(CLR_BLACK, 1);
      showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2 - 7, label, 12);
      showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2 + 7, "1: red, 2: yellow", 12);
    }
  } else {
    cr->set_source_rgba(CLR_BLACK, 1);
    showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2, label, 12);
  }
}
