/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_DRAW_PDFPAGE_H_
#define _DD_DRAW_PDFPAGE_H_

#include <cairommconfig.h>
#include <cairomm/context.h>
#include <cairomm/surface.h>
#include "dd_gv.hpp"
#include "dd_readfile.hpp"
#include "significancetest.hpp"
#include "dd_draw_environment_variable.hpp"
#include "dd_draw_myfunc.hpp"

class RGB;

class DParam {
public:
  int32_t pstart;
  int32_t pend;

  int32_t start;
  int32_t end;
  int32_t num_line;
  int32_t num_page;
  int32_t width_per_line;

  double yaxis_now;
  int32_t xstart;
  int32_t xend;

  double ystep;
  int32_t barnum;
  double height_df;

  int32_t width_draw;
  double dot_per_bp;

  double alpha;

  DParam(const int32_t s, const int32_t e, const DROMPA::Global &p):
    pstart(0), pend(0), start(s), end(e),
    num_line(p.drawparam.getNumLine(start, end)),
    num_page(p.drawparam.getNumPage(start, end)),
    width_per_line(p.drawparam.width_per_line),
    yaxis_now(0), xstart(0), xend(0), ystep(p.drawparam.getystep()), barnum(p.drawparam.getbarnum()),
    height_df(p.drawparam.getHeightDf()),
    width_draw(p.drawparam.width_draw_pixel),
    dot_per_bp(getratio(width_draw, width_per_line)),
    alpha(p.drawparam.alpha)
  {}

  void set_xstart_xend(const int32_t i) {
    xstart = start + i * width_per_line;
    if(i==num_line-1) xend = end;
    else xend = start + (i+1) * width_per_line -1;
  }

  double getXaxisLen() const { return (xend - xstart) * dot_per_bp; }
  double get_height_df() const { return height_df; }
};

class GraphData {
  enum {MEM_MAX_DEFAULT=20, MEMNUM_GC=10};
public:
  int32_t binsize;
  std::vector<double> array;
  std::string label;
  int32_t memnum;
  int32_t boxheight;

  double mmin;
  double mmax;
  double mwid;
  GraphData(): binsize(0), memnum(0), boxheight(0), mmin(0), mmax(0), mwid(0){}

  void setValue(const DROMPA::GraphDataFileName &g,
		const std::string &chr, const int32_t chrlen,
		const std::string &l,	const double ymin, const double ymax);

  double getylen(const int32_t i) const {
    return boxheight * (array[i] - mmin)/mwid;
  }
  double getBoxHeight4mem() const { return boxheight/memnum; }

  const std::string getmemory(const int32_t i) const {
    std::string str;
    double mem(mwid/memnum);
    if (mem <1) str = float2string(mmin + i*mem, 2);
    else        str = float2string(mmin + i*mem, 1);
    return str;
  }
};

class PDFPage {
  enum {GFTYPE_REFFLAT=0, GFTYPE_GTF=1, GFTYPE_SGD=2};

  const vChrArray &vReadArray;
  std::string chrname;
  const std::vector<SamplePairOverlayed> &vsamplepairoverlayed;
  GraphData GC, GD;

  Cairo::RefPtr<Cairo::Context> cr;

  int32_t setInterval() const {
    int32_t interval;
    if (par.width_per_line > 100*NUM_1M) interval = 10*NUM_1M;       // 10Mbp
    else if (par.width_per_line > 10*NUM_1M) interval = 100*NUM_1K;  // 100kbp
    else interval = par.width_per_line/10;
    return interval;
  }

  void StrokeEachLayer(const DROMPA::Global &p);
  void StrokeReadLines(const DROMPA::Global &p);
  void StrokeGraph(const GraphData &graph);
  void DrawIdeogram(const DROMPA::Global &p);
  void DrawGeneAnnotation(const DROMPA::Global &p);
  void strokeARS(const HashOfGeneDataMap &mp, const double ycenter);
  void strokeGeneSGD(const DROMPA::Global &p, const double ycenter);
  void strokeGene(const DROMPA::Global &p, const double ycenter);

  template <class T> void Draw_vbedlist(const std::vector<vbed<T>> &vlist);
  void drawBedAnnotation(const vbed<auto> &vbed);
  void drawInteraction(const InteractionSet &vinter);
  void StrokeChIADrop(const DROMPA::Global &p);
  void strokeChIADropBarcode(const std::vector<int32_t> &v, const std::string &nbarcode, const double _ywidth, const double yaxis, const RGB &color);

  public:
  DParam par;

  PDFPage(const DROMPA::Global &p,
	  const vChrArray &_vReadArray,
	  const std::vector<SamplePairOverlayed> &pair,
	  const Cairo::RefPtr<Cairo::PdfSurface> surface,
	  const int32_t s, const int32_t e):
    vReadArray(_vReadArray),
    chrname(vReadArray.getchr().getrefname()),
    vsamplepairoverlayed(pair),
    cr(Cairo::Context::create(surface)),
    par(s, e, p)
  {
    cr->select_font_face( "Arial", Cairo::FONT_SLANT_NORMAL, Cairo::FONT_WEIGHT_NORMAL);
    if(p.anno.GC.isOn()) GC.setValue(p.anno.GC, chrname, vReadArray.getchrlen(), "GC%",         20, 70);
    if(p.anno.GD.isOn()) GD.setValue(p.anno.GD, chrname, vReadArray.getchrlen(), "Num of genes", 0, 40);
  }

  template <class T>
  void StrokeDataFrame(const DROMPA::Global &p, const SamplePairOverlayed &pair)
  {
    par.yaxis_now += par.get_height_df() + MERGIN_BETWEEN_DATA;
    T df(cr, p, pair, par, chrname);
    df.Draw(p, pair, vReadArray);
    stroke_xaxis(par.yaxis_now);

    return;
  }

  void MakePage(const DROMPA::Global &p, const int32_t page_no, const std::string &pagelabel);

  void set_xstart_xend(const int32_t i) { par.set_xstart_xend(i); }
  void stroke_xaxis(const double y);
  void stroke_xaxis_num(const double y, const int32_t fontsize);
  void StrokeWidthOfInteractionSite(const bed &site, const double y);
  void drawArc_from_to(const Interaction &inter, const int32_t start, const int32_t end, const int32_t ref_height, const double ref_ytop);
  void drawArc_from_none(const Interaction &inter, const int32_t start, const int32_t end, const int32_t ref_height, const double ref_ytop);
  void drawArc_none_to(const Interaction &inter, const int32_t start, const int32_t end, const int32_t ref_height, const double ref_ytop);

  std::tuple<int32_t, int32_t> get_start_end_linenum(const int32_t page, const int32_t linenum_per_page) const {
    int32_t start(0), end(0);
    start = page * linenum_per_page;
    if(page == par.num_page-1) end = par.num_line;
    else end = (page+1) * linenum_per_page;
    return std::forward_as_tuple(start, end);
  }
};

#endif /* _DD_DRAW_PDFPAGE_H_ */
