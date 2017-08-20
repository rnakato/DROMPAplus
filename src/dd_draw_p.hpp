/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_DRAW_P_H_
#define _DD_DRAW_P_H_

namespace {
  enum class LineType { CHIP, INPUT, RATIO, RATIO_GV, PVALUE_INTER, PVALUE_ENRICH};

  enum {OFFSET_X=190, OFFSET_Y=50, MERGIN_BETWEEN_DATA=6, MERGIN_BETWEEN_LINE=30};
  int32_t pagewidth(1088);
  int32_t width_draw(820);
  int32_t mergin_between_graph_data(15);
  int32_t memnum_GC(10);
  int32_t boxheight_graph(80);
}

class Page {
  const std::unordered_map<std::string, ChrArray> &arrays;
  
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
  Cairo::RefPtr<Cairo::Context> cr;
  
  public:
  
  Page(const DROMPA::Global &p, 
       const std::unordered_map<std::string, ChrArray> &refarrays,
       Cairo::RefPtr<Cairo::PdfSurface> surface, const int32_t s, const int32_t e):
    arrays(refarrays),
    start(s), end(e),
    num_line(p.drawparam.getNumLine(start, end)),
    num_page(p.drawparam.getNumPage(start, end)),
    width_per_line(p.drawparam.width_per_line),
    dot_per_bp(getratio(width_draw, width_per_line)),
    yaxis_now(0), xstart(0), xend(0),
    ystep(12), barnum(2),
    cr(Cairo::Context::create(surface))
  {
    cr->select_font_face( "Serif", Cairo::FONT_SLANT_NORMAL, Cairo::FONT_WEIGHT_NORMAL );
  }

  void stroke_each_layer(const DROMPA::Global &p, const SamplePairChr &pair, const int32_t nlayer);
  void stroke_readdist(const DROMPA::Global &p, const SamplePairChr &pair, const LineType type, const int32_t nlayer);
  void stroke_dataframe(const DROMPA::Global &p, const SamplePairChr &pair, const LineType type, const int32_t nlayer);
  void stroke_bindata(const DROMPA::Global &p, const SamplePairChr &pair, const LineType type, const int32_t nlayer);
  void draw(const DROMPA::Global &p, const std::vector<SamplePairChr> &pairs, const int32_t page_curr);

  void set_xstart_xend(const int32_t i) {
    xstart = start + i * width_per_line;
    if(i==num_line-1) xend = end;
    else xend = start + (i+1) * width_per_line -1;
  }

  std::tuple<int32_t, int32_t> get_start_end_linenum(const int32_t page, const int32_t linenum_per_page) const {
    int32_t start(0), end(0);
    start = page * linenum_per_page;
    if(page == num_page-1) end = num_line;
    else end = (page+1) * linenum_per_page;
    return std::forward_as_tuple(start, end);
  }
  
  double bp2axis(const double v) const {
    return v * dot_per_bp + OFFSET_X;
  }

  double getWidthDf() const { return (xend - xstart +1) * dot_per_bp; }
  double getHeightDf() const { return ystep * barnum; }
};



#endif /* _DD_READFILE_P_H_ */
