/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_DRAW_P_H_
#define _DD_DRAW_P_H_

#include <cairommconfig.h>
#include <cairomm/context.h>
#include <cairomm/surface.h>
#include "color.hpp"


#define rel_xline(cr, x1, y1, xlen) do{		\
    cr->move_to(x1,   (int32_t)y1);		\
    cr->line_to(x1+xlen, (int32_t)y1); }while(0)
#define rel_yline(cr, x1, y1, ylen) do{		\
    cr->move_to(x1, (int32_t)y1);		\
    cr->line_to(x1, (int32_t)(y1+ylen)); }while(0)

#define CALCRATIO(c,i,r) ((i) ? ((c)/static_cast<double>((i)*(r))): 0)


namespace {
  enum class LineType { CHIP, INPUT, RATIO, RATIO_GV, PVALUE_INTER, PVALUE_ENRICH};

  enum {OFFSET_X=190, OFFSET_Y=50, MERGIN_BETWEEN_DATA=6, MERGIN_BETWEEN_LINE=30};
  int32_t pagewidth(1088);
  int32_t width_draw(820);
  int32_t mergin_between_graph_data(15);
  int32_t memnum_GC(10);
  int32_t boxheight_graph(80);
}

inline void showtext_cr(Cairo::RefPtr<Cairo::Context> cr, const double x, const double y, const std::string &str, const int32_t fontsize)
{
  cr->move_to(x, y);
  cr->set_font_size(fontsize);
  cr->show_text(str);
  cr->stroke();
  return;
}

inline std::string float2string(const double f, const int32_t digits)
{
  std::ostringstream oss;
  oss << std::setprecision(digits) << std::setiosflags(std::ios::fixed) << f;
  return oss.str();
}


class DParam {
public:
  int32_t pstart;
  int32_t pend;

  int32_t start;
  int32_t end;
  int32_t num_line;
  int32_t num_page;
  int32_t width_per_line;

  double dot_per_bp;
  double yaxis_now;
  int32_t xstart;
  int32_t xend;

  double ystep;
  int32_t barnum;

  DParam(const int32_t s, const int32_t e, const DROMPA::Global &p):
    start(s), end(e),
    num_line(p.drawparam.getNumLine(start, end)),
    num_page(p.drawparam.getNumPage(start, end)),
    width_per_line(p.drawparam.width_per_line),
    dot_per_bp(getratio(width_draw, width_per_line)),
    yaxis_now(0), xstart(0), xend(0), ystep(12), barnum(2)
  {}

  void set_xstart_xend(const int32_t i) {
    xstart = start + i * width_per_line;
    if(i==num_line-1) xend = end;
    else xend = start + (i+1) * width_per_line -1;
  }
};

class Page {
  const std::unordered_map<std::string, ChrArray> &arrays;
  const std::vector<SamplePairChr> &pairs;

  DParam par;
  Cairo::RefPtr<Cairo::Context> cr;
  
  public:
  
  Page(const DROMPA::Global &p,
       const std::unordered_map<std::string, ChrArray> &refarrays,
       const std::vector<SamplePairChr> &refpairs,
       Cairo::RefPtr<Cairo::PdfSurface> surface, const int32_t s, const int32_t e):
    arrays(refarrays), pairs(refpairs),
    par(s, e, p),
    cr(Cairo::Context::create(surface))
  {
    cr->select_font_face( "Arial", Cairo::FONT_SLANT_NORMAL, Cairo::FONT_WEIGHT_NORMAL );
  }
  
  //  void stroke_keys_dataframe(const DROMPA::Global &p, const SamplePairChr &pair, const LineType type);
  // void stroke_frame_dataframe();
  //void stroke_ymem_dataframe(const DROMPA::Global &p, const LineType type, const int32_t nlayer);
  void stroke_each_layer(const DROMPA::Global &p, const SamplePairChr &pair, const int32_t nlayer);

  template <class T>
  void stroke_readdist(const DROMPA::Global &p, const SamplePairChr &pair, const int32_t nlayer)
  {
    par.yaxis_now += getHeightDf() + MERGIN_BETWEEN_DATA;

    T df(cr, p, pair, par, getWidthDf(), getHeightDf());
    df.stroke_bindata(p, pair, arrays, nlayer);
    df.stroke_dataframe(p, nlayer);
    // if(!nlayer) stroke_xaxis(d, cr, xstart, xend);*/
    return;
  }

  
  //  void stroke_dataframe(const DROMPA::Global &p, const SamplePairChr &pair, const LineType type, const int32_t nlayer);
  //void stroke_bindata(const DROMPA::Global &p, const SamplePairChr &pair, const LineType type, const int32_t nlayer);
  void draw(const DROMPA::Global &p, const int32_t page_curr);

  void set_xstart_xend(const int32_t i) {
    par.set_xstart_xend(i);
  }

  std::tuple<int32_t, int32_t> get_start_end_linenum(const int32_t page, const int32_t linenum_per_page) const {
    int32_t start(0), end(0);
    start = page * linenum_per_page;
    if(page == par.num_page-1) end = par.num_line;
    else end = (page+1) * linenum_per_page;
    return std::forward_as_tuple(start, end);
  }

  double getWidthDf() const { return (par.xend - par.xstart +1) * par.dot_per_bp; }
  double getHeightDf() const { return par.ystep * par.barnum; }
};


class DataFrame {

protected:
  const DParam &par;
  const Cairo::RefPtr<Cairo::Context> cr;
  double scale;
  std::string label;
  double width_df;
  double height_df;
  
 public:
  DataFrame(const Cairo::RefPtr<Cairo::Context> cr_, const std::string &l, const double s,
	    const DParam &refparam, const double wdf, const double hdf):
    cr(cr_), par(refparam), scale(s), label(l), width_df(wdf), height_df(hdf)
  {}

  void stroke_frame()
  {
    cr->set_line_width(0.4);
    cr->set_source_rgba(CLR_BLACK, 1);
    rel_xline(cr, OFFSET_X, par.yaxis_now, width_df);
    rel_yline(cr, OFFSET_X, par.yaxis_now - height_df, height_df);
    cr->stroke();
  }
  void stroke_ymem(const int32_t nlayer)
  {
    cr->set_source_rgba(CLR_BLACK, 1);
  
    double x(0);
    if (!nlayer) x = OFFSET_X + width_df + 7; else x = OFFSET_X - 20;
    for(int32_t i=1; i<=par.barnum; ++i) {
      std::string str(float2string(i*scale, 1));
      showtext_cr(cr, x, par.yaxis_now - i*(par.ystep - 1.5), str, 9);
    }
    return;
  }
  void stroke_dataframe(const DROMPA::Global &p, const int32_t nlayer) {
    stroke_frame();
    /* y memory */
    cr->set_line_width(0.4);
    cr->set_source_rgba(CLR_BLACK, 0.5);
    for(int32_t i=0; i<par.barnum; ++i) rel_xline(cr, OFFSET_X, par.yaxis_now - i*par.ystep, width_df);
    cr->stroke();

    if (p.drawparam.isshowymem()) stroke_ymem(nlayer);
    if (!nlayer && p.drawparam.isshowylab()) {
      cr->set_source_rgba(CLR_BLACK, 1);
      showtext_cr(cr, 50, par.yaxis_now - height_df/2, label, 12);
    }
    return;
  }
  
  void stroke_bindata(const DROMPA::Global &p, const SamplePairChr &pair,
		      const std::unordered_map<std::string, ChrArray> &arrays, const int32_t nlayer);

  double bp2axis(const double v) const { return v * par.dot_per_bp + OFFSET_X; }
  int32_t getbinlen(const double value) const { return -std::min(par.ystep*value, height_df); }
  
  virtual void stroke_bin(const SamplePairChr &pair,
			  const std::unordered_map<std::string, ChrArray> &arrays,
			  const int32_t i, const double xcen, const int32_t yaxis, const int32_t viz)=0;
};

#endif /* _DD_READFILE_P_H_ */
