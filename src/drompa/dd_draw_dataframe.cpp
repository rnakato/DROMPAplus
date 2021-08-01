/* Copyright(c)  Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "dd_draw_dataframe.hpp"

void DataFrame::Draw(const DROMPA::Global &p,
                     const SamplePairOverlayed &pair,
                     const vChrArray &vReadArray)
{
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

  if (p.drawparam.isshowylab()) StrokeSampleLabel(pair);
}

void DataFrame::StrokeYmem(const int32_t nlayer)
{
  cr->set_line_width(0.4);
  cr->set_source_rgba(CLR_BLACK, 0.5);
  for (int32_t i=0; i<par.barnum; ++i) rel_xline(cr, OFFSET_X, par.yaxis_now - i*par.ystep, par.getXaxisLen());
  cr->stroke();

  cr->set_source_rgba(CLR_BLACK, 1);
  double x(0);
  if (!nlayer) x = OFFSET_X + par.getXaxisLen() + 7;
  else x = OFFSET_X - 20;

  for(int32_t i=0; i<=par.barnum; ++i) {
    if(!i && !bothdirection && !shownegative) continue;
    double y_value;
    if (!nlayer) y_value = get_yscale_num(i, scale);
    else         y_value = get_yscale_num(i, scale2nd);
    showtext_cr(cr, x, par.yaxis_now - i*(par.ystep - 1.5), float2string(y_value, ndigit), 9);
  }
}

void DataFrame::StrokeSampleLabel(const SamplePairOverlayed &pair) {
  cr->set_source_rgba(CLR_BLACK, 1);
  showtext_cr(cr, POSI_ASSAYLABEL, par.yaxis_now - height_df/2, getAssayName(), 10);

  if (pair.OverlayExists()) {
    getColor1st(1);
    showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2 - 7, label, 12);
    getColor2nd(1);
    showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2 + 7, label2nd, 12);
  } else {
    cr->set_source_rgba(CLR_BLACK, 1);
    showtext_cr(cr, POSI_XLABEL, par.yaxis_now - height_df/2, label, 12);
  }
}

double DataFrame::getmax(const SamplePairEach &pair,
              const vChrArray &vReadArray,
              const int32_t i, const int32_t thin)
{
  double value(0);
  double v;
  for (int32_t j=0; j<thin; ++j) {
    v = getVal(pair, vReadArray, i+j);
    if (abs(v) > abs(value)) value = v;
  }
  return value;
}

void DataFrame::StrokeBins(const SamplePairEach &pair,
                           const vChrArray &vReadArray,
                           const int32_t nlayer)
{
  int32_t binsize(pair.getbinsize());
  int32_t sbin(par.xstart/binsize);
  int32_t ebin(par.xend/binsize);
  double dot_per_bin(binsize * par.dot_per_bp);
  int32_t yaxis(par.yaxis_now);  // convert to int

  int32_t thin(std::min(par.width_per_line/(1000*binsize), 20));
  if(!thin) thin = 1;

/*  printf("\n\nthin %d\n", thin);
  printf("\n\ndot_per_bin %f\n", dot_per_bin);
  printf("\n\nthin %d\n", par.width_per_line/(1000*binsize));*/

  double xcen(BP2PIXEL(binsize/2)); // initial position
  if (thin > 1) cr->set_line_width(dot_per_bin*thin);
  else cr->set_line_width(dot_per_bin);

  for (int32_t i=sbin; i<ebin; i += thin, xcen += dot_per_bin*thin) {
    double value(0);
    if (thin > 1) value = getmax(pair, vReadArray, i, thin);
    else value = getVal(pair, vReadArray, i);

    StrokeEachBin(value, xcen, yaxis, nlayer);
 }

  /*for (int32_t i=sbin; i<ebin; ++i, xcen += dot_per_bin) {
    if (thin > 1 && i%thin) continue;
    double value(getVal(pair, vReadArray, i));
    StrokeEachBin(value, xcen, yaxis, nlayer);
  }*/
  cr->stroke();
}

void DataFrame::StrokeEachBin(const double value, const double xcen,
                              const int32_t yaxis, const int32_t nlayer)
{
  if (!value) return;

  double len(0);
  if (!nlayer) len = par.ystep * value / scale;
  else         len = par.ystep * value / scale2nd;
  if (shownegative) {
    if (value > 0)      len = -std::min(len,  len_plus);
    else if (value < 0) len = -std::max(len, -len_minus);
  } else                len = -std::min(len,  height_df);
  if (!len) return;

  /*  DEBUGprint("value " << value
      << " ystep:" << par.ystep
      << " scale:" << scale
      << " len:" << len
      << " height_df " << height_df
      << " yaxis" << yaxis
      << " len_binedge" << len_binedge
      << " len_minus " << len_minus);*/

  // waku
  setColor(value, nlayer, 1);
  rel_yline(cr, xcen, yaxis -2 + len - len_minus*shownegative, LEN_EDGE);
  cr->stroke();

  // nakami
  if (par.alpha) {
    setColor(value, nlayer, par.alpha);
    rel_yline(cr, xcen, yaxis - len_minus*shownegative, len);
    cr->stroke();
  }
}

void LogRatioDataFrame::StrokeEachBin(const double _value, const double xcen,
                                      const int32_t yaxis, const int32_t nlayer)
{
  double value(0);
  if (!nlayer) value = _value ? log(_value) / log(scale): 0;
  else         value = _value ? log(_value) / log(scale2nd): 0;

  int32_t len(0);
  if (value > 0)      len = -std::min(par.ystep*value,  len_plus);
  else if (value < 0) len = -std::max(par.ystep*value, -len_minus);
  if (!len) return;

  // waku
  setColor(value, nlayer, 1);
  rel_yline(cr, xcen, yaxis + len - len_minus, LEN_EDGE);
  cr->stroke();

  // nakami
  if (par.alpha) {
    setColor(value, nlayer, par.alpha);
    rel_yline(cr, xcen, yaxis - len_minus, len);
    cr->stroke();
  }
}

// SetColor
void DataFrame::setColor(const double value, const int32_t nlayer, const double alpha)
{
  (void)(value);
  if (!nlayer) getColor1st(alpha);
  else getColor2nd(alpha);
}

void RatioDataFrame::setColor(const double value, const int32_t nlayer, const double alpha)
{
  if (!nlayer) { // first layer
    if(isGV || sigtest) {
      if (value >= threshold) cr->set_source_rgba(CLR_RED, alpha);
      else cr->set_source_rgba(CLR_GRAY, alpha);
    } else getColor1st(alpha);
  } else {    // second layer
    if(isGV || sigtest) {
      if (value >= threshold) cr->set_source_rgba(CLR_DARKORANGE, alpha);
      else cr->set_source_rgba(CLR_GRAY2, alpha);
    } else getColor2nd(alpha);
  }
}

void LogRatioDataFrame::setColor(const double value, const int32_t nlayer, const double alpha)
{
  if (!nlayer) {
    if (isGV || sigtest) {
      if (value >= 0) cr->set_source_rgba(CLR_SALMON, 1);
      else cr->set_source_rgba(CLR_BLUEPURPLE, 1);
    } else getColor1st(alpha);
  } else {
    if (isGV || sigtest) {
      if (value >= 0) cr->set_source_rgba(CLR_DEEPSKYBLUE, alpha);
      else cr->set_source_rgba(CLR_GRAY2, alpha);
    } else getColor2nd(alpha);
  }
}

void PvalueDataFrame::setColor(const double value, const int32_t nlayer, const double alpha)
{
  if (!nlayer) { // first layer
    if (value >= threshold) getColor1st(alpha);
    else cr->set_source_rgba(CLR_GRAY, alpha);
  } else {    // second layer
    if (value >= threshold) getColor2nd(alpha);
    else cr->set_source_rgba(CLR_GRAY, alpha);
  }
}
