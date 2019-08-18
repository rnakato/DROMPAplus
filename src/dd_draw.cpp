/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */

#include <sstream>
#include <iomanip>
#include "dd_draw.hpp"
#include "dd_draw_pdfpage.hpp"
#include "dd_draw_dataframe.hpp"
#include "../submodules/SSP/common/inline.hpp"
#include "../submodules/SSP/common/util.hpp"

/*    if (p->gapfile && d->gaparray[i] >= GAP_THRE){
      cr->set_source_rgba(CLR_BLACK, 0.3);
      rel_yline(cr, xcen, yaxis - height_df, height_df);
      cr->stroke();
      } else if(p->mpfile && d->mparray[i] < p->mpthre){
      cr->set_source_rgba(CLR_BLACK, 0.3);
      rel_yline(cr, xcen, yaxis - height_df, height_df);
      cr->stroke();
      }*/

namespace {
  class posivector {
    int32_t start;
  public:
    std::vector<int32_t> v;
    posivector(std::vector<int32_t> _v, int32_t xstart): start(xstart), v(_v) {}
    virtual ~posivector(){}

    bool operator < (const posivector &another) const {
      for (size_t i=0; i<std::min(v.size(), another.v.size()); ++i) {
	if (v[i] < start && another.v[i] < start) continue;
	if (v[i] < another.v[i]) return 1;
	if (v[i] == another.v[i] && i==v.size()-1) return 1;
	if (v[i] == another.v[i] && i==another.v.size()-1) return 0;
	if (v[i] > another.v[i]) return 0;
      }
      return 1;
    };
    bool operator == (const posivector &another) const {
      if(v.size() != another.v.size()) return 0;

      for (size_t i=0; i<v.size(); ++i) {
	if (v[i] != another.v[i]) return 0;
      }
      return 1;
    };
  };

  void strokeGraph4EachWindow(const Cairo::RefPtr<Cairo::Context> cr,
			      double x_pre, double y_pre,
			      double x_cen, double y_cen,
			      int32_t bottom)
  {
    if(y_pre > bottom && y_cen > bottom) return;

    double x1(x_pre);
    double x2(x_cen);
    double y1(y_pre);
    double y2(y_cen);
    if (y_pre > bottom || y_cen > bottom) {
      double xlen = abs(x2 - x1);
      double ylen = abs(y2 - y1);
      if (y_pre > bottom) {
	double ydiff(abs(y_cen - bottom));
	x1 = x2 - xlen*(ydiff/ylen);
	y1 = bottom;
      } else {
	double ydiff(abs(y_pre - bottom));
	x2 = x1 - xlen*(ydiff/ylen);
	y2 = bottom;
      }
    }

    cr->move_to(x1, y1);
    cr->line_to(x2, y2);
    cr->stroke();
    return;
  }

  RGB getInterRGB(double val)
  {
    val = val*0.4 + 0.6;
    HSV color(val, 1.0, 1.0);

    return HSVtoRGB(color);
  }

  void showColorBar(const Cairo::RefPtr<Cairo::Context> cr, int32_t x, const int32_t y, const double maxval)
  {
    int32_t barwidth(1);
    cr->set_line_width(barwidth+0.1);

    cr->set_source_rgba(CLR_BLACK, 1);
    showtext_cr(cr, x-36, y+20, "-log(p)", 10);

    for (int32_t i=0; i<=50; ++i) {
      RGB color(getInterRGB(i*0.02));
      cr->set_source_rgba(color.r, color.g, color.b, 0.8);
      rel_yline(cr, x, y, 12);
      cr->stroke();
      cr->set_source_rgba(CLR_BLACK, 1);
      if (!i)    showtext_cr(cr, x-3, y+20, "0", 10);
      if (i==40) showtext_cr(cr, x-3, y+20, float2string(maxval/3, 1), 10);
      x += barwidth;
      cr->stroke();
    }
    return;
  }
  void showColorBar_ChIADrop(const Cairo::RefPtr<Cairo::Context> cr, int32_t x, const int32_t y, const int32_t maxval)
  {
    int32_t barwidth(1);
    cr->set_line_width(barwidth+0.1);

    cr->set_source_rgba(CLR_BLACK, 1);
    showtext_cr(cr, x-15, y+30, "num of barcodes", 10);

    for (int32_t i=0; i<=50; ++i) {
      RGB color(getInterRGB(i*0.02));
      cr->set_source_rgba(color.r, color.g, color.b, 0.8);
      rel_yline(cr, x, y, 12);
      cr->stroke();
      cr->set_source_rgba(CLR_BLACK, 1);
      if (!i)    showtext_cr(cr, x-3, y+20, "1", 10);
      if (i==48) showtext_cr(cr, x-3, y+20, std::to_string(maxval), 10);
      x += barwidth;
      cr->stroke();
    }
    return;
  }
}

void PDFPage::StrokeWidthOfInteractionSite(const bed site, const double y)
{
  cr->set_line_width(2);
  cr->set_source_rgba(CLR_DARKORANGE, 0.8);
  double s = par.bp2xaxis(site.start - par.xstart);
  double e = par.bp2xaxis(site.end - par.xstart);
  rel_xline(cr, s, y, e-s);
  cr->stroke();
}

// cr->arc(中心x, 中心y, 半径, start角度, end角度) 角度はラジアン
void PDFPage::drawArc_from_to(const Interaction &inter, const int32_t start, const int32_t end, const int32_t ref_height, const double ref_ytop)
{
  double ytop = ref_ytop + 10;
  int32_t height = ref_height - 20;
  double radius((end - start)/2.0 * par.dot_per_bp); // 半径
  double r = std::min(0.4, height/radius);
  //    printf("r %f %f %d %d %d\n", r, radius, height, start, end);

  cr->set_line_width(3);
  cr->scale(1, r);
  cr->arc(par.bp2xaxis((start + end) /2), ytop/r, radius, 0, M_PI);
  cr->stroke();
  cr->scale(1, 1/r);

  // Highlight each site of interaction
  StrokeWidthOfInteractionSite(inter.first, ytop);
  StrokeWidthOfInteractionSite(inter.second, ytop);
}

void PDFPage::drawArc_from_none(const Interaction &inter, const int32_t start, const int32_t end, const int32_t ref_height, const double ref_ytop)
{
  double ytop = ref_ytop + 10;
  int32_t height = ref_height;
  double radius(height*3);
  double r(1/3.0);

  double bp_s(par.bp2xaxis(start));
  double bp_e(par.bp2xaxis(end));
  double bp_x(bp_s + radius);
  double bp_y(ytop/r);

  cr->set_line_width(4);
  cr->scale(1, r);
  cr->arc(bp_x, bp_y, radius, 0.5*M_PI, M_PI);
  if (bp_e - bp_x > 0) rel_xline(cr, bp_x, bp_y + radius, bp_e - bp_x);
  cr->stroke();
  cr->scale(1, 1/r);

  // bin of interaction
  StrokeWidthOfInteractionSite(inter.first, ytop);
}

void PDFPage::drawArc_none_to(const Interaction &inter, const int32_t start, const int32_t end, const int32_t ref_height, const double ref_ytop)
{
  double ytop = ref_ytop + 10;
  int32_t height = ref_height;
  double radius(height*3);
  double r(1/3.0);

  double bp_s(par.bp2xaxis(start));
  double bp_e(par.bp2xaxis(end));
  double bp_x(bp_e - radius);
  double bp_y(ytop/r);

  cr->set_line_width(4);
  cr->scale(1, r);
  cr->arc(bp_x, bp_y, radius, 0, 0.5*M_PI);
  if (bp_x - bp_s > 0) rel_xline(cr, bp_x, bp_y + radius, -(bp_x - bp_s));
  cr->stroke();
  cr->scale(1, 1/r);

  // bin of interaction
  StrokeWidthOfInteractionSite(inter.second, ytop);
}

void PDFPage::drawInteraction(const InteractionSet &vinter)
{
  DEBUGprint("drawInteraction");
  int32_t boxheight(BOXHEIGHT_INTERACTION);
  std::string chr(rmchr(chrname));
  double ytop(par.yaxis_now);
  double ycenter(par.yaxis_now + boxheight/2);

  // colorbar
  showColorBar(cr, 70, ycenter-2, vinter.getmaxval());

  // label
  cr->set_source_rgba(CLR_BLACK, 1);
  showtext_cr(cr, 70, ycenter-6, vinter.getlabel(), 12);

  for (auto &x: vinter.getvinter()) {
    //    printList(x.first.chr, x.second.chr, chr);
    if (x.first.chr != chr && x.second.chr != chr) continue;

    // interchromosomalは描画しない
    if (x.first.chr != chr || x.second.chr != chr) continue;

    RGB color(getInterRGB(x.getval()/vinter.getmaxval() *3)); // maxval の 1/3 を色のmax値に設定
    cr->set_source_rgba(color.r, color.g, color.b, 0.8);
    /*    else {   // inter-chromosomal
	  RGB color(getInterRGB(x.getval()/vinter.getmaxval() *3)); // maxval の 1/3 を色のmax値に設定
	  cr->set_source_rgba(color.r, color.g, color.b, 0.8);
	  }*/

    int32_t xcen_head(-1);
    int32_t xcen_tail(-1);
    if (par.xstart <= x.first.summit  && x.first.summit  <= par.xend) xcen_head = x.first.summit  - par.xstart;
    if (par.xstart <= x.second.summit && x.second.summit <= par.xend) xcen_tail = x.second.summit - par.xstart;

    if (xcen_head < 0 && xcen_tail < 0) continue;

    //    printf("%d, %d, %d, %d, %d, %d\n", x.first.start, x.first.summit, x.first.end, x.second.start, x.second.summit, x.second.end);
    if (xcen_head >= 0 && xcen_tail >= 0) drawArc_from_to(x, xcen_head, xcen_tail, boxheight, ytop);
    if (xcen_head > 0 && xcen_tail < 0)   drawArc_from_none(x, xcen_head, par.xend - par.xstart, boxheight, ytop);
    if (xcen_head < 0 && xcen_tail > 0)   drawArc_none_to(x, xcen_head, xcen_tail, boxheight, ytop);
  }
  cr->stroke();
  par.yaxis_now += boxheight;

  return;
}

void PDFPage::stroke_xaxis(const double y)
{
  double x;
  int32_t interval_large(setInterval());
  int32_t interval(interval_large/10);

  cr->set_source_rgba(CLR_BLACK, 1);
  for(int32_t i=setline(par.xstart, interval); i<=par.xend; i+=interval) {
    x = par.bp2xaxis(i - par.xstart);
    if (!(i%interval_large)) {
      cr->set_line_width(1);
      rel_yline(cr, x, y-4, 8);
    } else {
      cr->set_line_width(0.5);
      rel_yline(cr, x, y-1.5, 3);
    }
    cr->stroke();
  }
  return;
}

void PDFPage::StrokeGraph(const GraphData &graph)
{
  DEBUGprint("PDFPage::DrawGraph");
  int32_t s(par.xstart/graph.binsize);
  int32_t e(par.xend/graph.binsize +1);
  double diff = graph.binsize * par.dot_per_bp;

  double ytop(par.yaxis_now);
  double ycenter(par.yaxis_now + graph.boxheight/2);
  double ybottom(par.yaxis_now + graph.boxheight);

  // graph line
  cr->set_line_width(0.6);
  cr->set_source_rgba(CLR_GREEN, 1);
  double xpre(OFFSET_X);
  double xcen(par.bp2xaxis(0.5*graph.binsize));
  double ypre(ybottom - graph.getylen(s));
  for (int32_t i=s; i<e; ++i, xcen += diff) {
    double ycen(ybottom - graph.getylen(i));
    strokeGraph4EachWindow(cr, xpre, ypre, xcen, ycen, ybottom + 10);
    xpre = xcen;
    ypre = ycen;
  }

  cr->set_source_rgba(CLR_BLACK, 1);

  // label
  showtext_cr(cr, OFFSET_X - 5*graph.label.length() -55, ycenter, graph.label, 13);

  // x-axis
  stroke_xaxis(ybottom);

  // y-axis
  cr->set_line_width(0.4);
  rel_yline(cr, OFFSET_X, ytop, graph.boxheight);
  cr->stroke();
  cr->set_line_width(1.5);
  rel_xline(cr, OFFSET_X, ybottom, (par.xend - par.xstart+1) * par.dot_per_bp);
  cr->stroke();

  // y memory
  cr->set_line_width(0.5);
  for (int32_t i=0; i<=graph.memnum; ++i) {
    std::string str(graph.getmemory(i));
    double x(OFFSET_X - 5*str.length() - 7);
    double y(ybottom - i*graph.getBoxHeight4mem());
    showtext_cr(cr, x, y+2, str, 9);
    rel_xline(cr, OFFSET_X-2, y, 2);
    cr->stroke();
  }
  par.yaxis_now += graph.boxheight + MERGIN_BETWEEN_GRAPH_DATA;

  return;
}

void PDFPage::drawBedAnnotation(const vbed<bed12> &vbed)
{
  DEBUGprint("drawBedAnnotation");
  int32_t boxheight(BOXHEIGHT_BEDANNOTATION);
  std::string chr(rmchr(chrname));
  double ycenter(par.yaxis_now + boxheight/2);

  // label
  cr->set_source_rgba(CLR_BLACK, 1);
  showtext_cr(cr, 90, ycenter+4, vbed.getlabel(), 12);

  cr->set_line_width(0.5);
  rel_xline(cr, OFFSET_X, ycenter, par.getXaxisLen());
  cr->stroke();

  // bed
  cr->set_line_width(boxheight/2);
  int32_t on(0);
  for (auto &x: vbed.getvBed()) {
    if (x.chr != chr) continue;

    if (x.rgb_r != -1) {
      cr->set_source_rgba(x.rgb_r/(double)255, x.rgb_g/(double)255, x.rgb_b/(double)255, 0.6);
    } else {
      if(!on) {
	cr->set_source_rgba(CLR_GRAY4, 1);
	on=1;
      } else {
	cr->set_source_rgba(CLR_GREEN, 1);
	on=0;
      }
    }
    if (par.xstart <= x.end && x.start <= par.xend) {
      double x1 = par.bp2xaxis(x.start - par.xstart);
      double len = (x.end - x.start) * par.dot_per_bp;
      rel_xline(cr, x1, ycenter, len);
      cr->stroke();
    }
  }
  cr->stroke();
  par.yaxis_now += boxheight;

  return;
}

void PDFPage::stroke_xaxis_num(const double y, const int32_t fontsize)
{
  int32_t mega, kilo;
  double x;
  int32_t interval(setInterval());

  cr->set_source_rgba(CLR_BLACK, 1);
  for(int32_t i=setline(par.xstart, interval); i<=par.xend; i+=interval) {
    std::string str;
    x = par.bp2xaxis(i - par.xstart);
    if (par.width_per_line > 100*NUM_1M)     str = float2string(i/static_cast<double>(NUM_1M), 1) + "M";
    else if (par.width_per_line > 10*NUM_1M) str = float2string(i/static_cast<double>(NUM_1K), 1) + "k";
    else {
      mega = i/NUM_1M;
      kilo = (i%NUM_1M)/NUM_1K;
      if (par.width_per_line > 10*NUM_1K) str = float2string(i/static_cast<double>(NUM_1M), 3) + "M";
      else if (par.width_per_line > 10) {
	if (mega) str = std::to_string(mega) + "," + float2string((i%NUM_1M)/static_cast<double>(NUM_1K), 3) + "K";
	else str = float2string((i%NUM_1M)/static_cast<double>(NUM_1K), 3) + "K";
      } else {
	if (mega) str = std::to_string(mega) + "," + std::to_string(kilo) + "," + std::to_string(i%NUM_1K);
	else if (kilo) str = std::to_string(kilo) + "," + std::to_string(i%NUM_1K);
	else str = std::to_string(i%NUM_1K);
      }
    }
    showtext_cr(cr, x - 3*str.length(), y+10, str, fontsize);
  }
  return;
}

void PDFPage::StrokeReadLines(const DROMPA::Global &p)
{
  auto &d = p.drawparam;

  for (auto &pair: vsamplepairoverlayed) {
    if (d.showpinter) StrokeDataFrame<PinterDataFrame>(p, pair);
    if (d.showpenrich && pair.first.InputExists()) StrokeDataFrame<PenrichDataFrame>(p, pair);
    if (d.showratio && pair.first.InputExists()) {
      if (d.showratio == 1)      StrokeDataFrame<RatioDataFrame>(p, pair);
      else if (d.showratio == 2) StrokeDataFrame<LogRatioDataFrame>(p, pair);
    }
    if (d.showctag) StrokeDataFrame<ChIPDataFrame>(p, pair);
    if (d.showitag==1 && pair.first.InputExists()) StrokeDataFrame<InputDataFrame>(p, pair);
  }
  if (d.showitag==2 && vsamplepairoverlayed[0].first.InputExists()) StrokeDataFrame<InputDataFrame>(p, vsamplepairoverlayed[0]);

  stroke_xaxis_num(par.yaxis_now, 9);
  return;
}

void PDFPage::strokeChIADropBarcode(const std::vector<int32_t> &v, const std::string &nbarcode, const double _ywidth, const RGB &color)
{
  double ywidth = std::min(_ywidth, 2.0);
  double ycenter(par.yaxis_now + ywidth/2);

  int32_t s = std::max(v[0], par.xstart);
  int32_t e = std::min(v[v.size()-1], par.xend);
//  cr->set_source_rgba(CLR_GRAY2, 1);
  cr->set_source_rgba(color.r, color.g, color.b, 0.8);
  cr->set_line_width(ywidth*0.1);
  rel_xline_double(cr, par.bp2xaxis(s - par.xstart), ycenter, (e-s) * par.dot_per_bp);
  cr->stroke();

  // barcode number
  cr->set_source_rgba(CLR_BLACK, 1);
  showtext_cr(cr, par.bp2xaxis(e - par.xstart) + 0.5, par.yaxis_now + ywidth, nbarcode, 1.0);
  cr->stroke();

  cr->set_line_width(ywidth*0.8);
  cr->set_source_rgba(color.r, color.g, color.b, 0.8);
  for (auto &posi: v) {
    if(posi >= par.xstart && posi <= par.xend) {
      double x1 = par.bp2xaxis(posi - par.xstart);
      double len = std::max(1000 * par.dot_per_bp, 0.02);
      rel_xline_double(cr, x1, ycenter, len);
      cr->stroke();
    }
  }

  par.yaxis_now += ywidth;
}

void PDFPage::StrokeChIADrop(const DROMPA::Global &p)
{
  DEBUGprint("StrokeChIADrop");
  int32_t boxheight(BOXHEIGHT_ChIADROP);
  std::string chr(rmchr(chrname));

  /* frame */
  cr->set_source_rgba(CLR_GRAY4, 1);
  cr->rectangle(OFFSET_X, par.yaxis_now, (par.xend - par.xstart+1) * par.dot_per_bp, boxheight);
  cr->stroke();

  std::vector<posivector> vv;

  for (auto &x: p.anno.mp_ChIADrop.at(chr)) {
    const std::vector<int32_t> &v = x.second;
    if (v.size() ==1) continue;
    if (v[v.size()-1] < par.xstart || v[0] > par.xend) continue;
    if (v[0] < par.xstart && v[v.size()-1] > par.xend) continue;
    vv.emplace_back(v, par.xstart);
  }

  std::sort(vv.begin(), vv.end());

  int32_t max(0);
  int32_t num_line(0);
  for (size_t i=0; i<vv.size(); ++i) {
    int32_t n(1);
    while(i < vv.size()-1 && vv[i].v == vv[i+1].v) { ++i; ++n; }
    max = std::max(max, n);
    ++num_line;
  }

  // colorbar
  showColorBar_ChIADrop(cr, 80, par.yaxis_now + 10, max);

  for (size_t i=0; i<vv.size(); ++i) {
    int32_t n(1);
    while(i < vv.size()-1 && vv[i].v == vv[i+1].v) { ++i; ++n; }
    RGB color(getInterRGB((n-1)/(double)max));
    strokeChIADropBarcode(vv[i-n+1].v, std::to_string(n), boxheight/(double)num_line, color);
  }

  return;
}

void PDFPage::DrawIdeogram(const DROMPA::Global &p)
{
  DEBUGprint("PDFPage::DrawIdeogram");
  int32_t boxheight(BOXHEIGHT_IDEOGRAM);
  int32_t on(0);
  int32_t acen_once(0);

  cr->set_line_width(1);
  // frame
  cr->set_source_rgba(CLR_BLACK, 1);
  cr->rectangle(OFFSET_X, par.yaxis_now, (par.xend - par.xstart+1) * par.dot_per_bp, boxheight);
  cr->stroke();

  for (auto &x: p.anno.vcytoband) {
    if (rmchr(chrname) != x.chr) continue;
    //    x.print();
    if (x.stain == "acen") cr->set_source_rgba(CLR_RED, 1);
    else if (x.stain == "gneg") cr->set_source_rgba(CLR_GRAY0, 1);
    else if (x.stain == "gpos25" || x.stain == "stalk") cr->set_source_rgba(CLR_GRAY, 1);
    else if (x.stain == "gpos50") cr->set_source_rgba(CLR_GRAY2, 1);
    else if (x.stain == "gpos75") cr->set_source_rgba(CLR_GRAY3, 1);
    else if (x.stain == "gpos100" || x.stain == "gvar") cr->set_source_rgba(CLR_GRAY4, 1);
    else { std::cout << "Warning: stain " << x.stain << " is not annotated." << std::endl; }

    double s(OFFSET_X + x.start * par.dot_per_bp);
    double len((x.end - x.start +1) * par.dot_per_bp);

    if (x.stain == "acen") {
      if(!acen_once) {
	mytriangle(s,     par.yaxis_now,
		   s,     par.yaxis_now + boxheight,
		   s+len, par.yaxis_now + boxheight/2);
	++acen_once;
      } else {
	mytriangle(s,     par.yaxis_now + boxheight/2,
		   s+len, par.yaxis_now,
		   s+len, par.yaxis_now + boxheight);
      }
    } else {
      cr->rectangle(s, par.yaxis_now, len, boxheight);
    }
    cr->fill();
    cr->stroke();

    cr->set_source_rgba(CLR_BLACK, 1);
    double y;
    if(on) y = par.yaxis_now + boxheight + 5;
    else y = par.yaxis_now - 3;
    showtext_cr(cr, s+1, y, x.name, 5);
    if(on) on=0; else { ++on; }
  }

  par.yaxis_now += boxheight + MERGIN_BETWEEN_GRAPH_DATA;

  return;
}


void PDFPage::StrokeEachLayer(const DROMPA::Global &p)
{
  if (p.anno.showIdeogram()) DrawIdeogram(p);

  if (p.anno.GC.isOn()) StrokeGraph(GC);
  if (p.anno.GD.isOn()) StrokeGraph(GD);

  // Gene
  if (p.anno.genefile != "" || p.anno.arsfile != "" || p.anno.terfile != "")DrawGeneAnnotation(p);

  // Read
  StrokeReadLines(p);

  // Bed file
  if(p.anno.vbedlist.size()) {
    par.yaxis_now += MERGIN_BETWEEN_READ_BED;
    for (auto &x: p.anno.vbedlist) {
      drawBedAnnotation(x);
      par.yaxis_now += +2;
    }
  }

  // ChIA-Drop
  if (p.anno.existChIADrop()) {
    par.yaxis_now += MERGIN_BETWEEN_READ_BED;
    StrokeChIADrop(p);
  }

  // Interaction
  if (p.anno.vinterlist.size()) {
    par.yaxis_now += MERGIN_BETWEEN_READ_BED;
    for (auto &x: p.anno.vinterlist) {
      drawInteraction(x);
      par.yaxis_now += +5;
    }
  }

  //if(d->repeat.argv) draw_repeat(d, cr, xstart, xend);

  return;
}

void PDFPage::MakePage(const DROMPA::Global &p,
		       const int32_t page_no,
		       const std::string &pagelabel)
{
  DEBUGprint("PDFPage::MakePage");
  int32_t line_start, line_end;
  std::tie(line_start, line_end) = get_start_end_linenum(page_no, p.drawparam.getlpp());

  par.yaxis_now = OFFSET_Y;
  cr->set_source_rgba(CLR_WHITE, 1);
  cr->paint();

  // Stroke each layer
  for (int32_t i=line_start; i<line_end; ++i) {
    set_xstart_xend(i);
    if (par.xstart >= par.xend) continue;
    StrokeEachLayer(p);
    if (i != line_end-1) par.yaxis_now += MERGIN_BETWEEN_EACHLAYER;
  }

  // Page title
  std::string title;
  if(par.num_page>1) title = chrname + "_" + pagelabel + "_" + std::to_string(page_no+1);
  else               title = chrname + "_" + pagelabel;
  cr->set_source_rgba(CLR_BLACK, 1);
  showtext_cr(cr, 50, 30, title, 16);
  cr->stroke();

  cr->show_page();
  return;
}

void Figure::Draw_Region(DROMPA::Global &p,
			 std::string &pdffilename,
			 int32_t width,
			 int32_t height)
{
  const auto surface = Cairo::PdfSurface::create(pdffilename, width, height);
  int32_t region_no(1);
  for (auto &x: regionBed) {
    int32_t num_page = p.drawparam.getNumPage(x.start, x.end);
    for(int32_t i=0; i<num_page; ++i) {
      std::cout << boost::format("   page %5d/%5d/%5d\r")
	% (i+1) % num_page % region_no << std::flush;
      PDFPage page(p, vReadArray, vsamplepairoverlayed, surface, x.start, x.end);
      page.MakePage(p, i, std::to_string(region_no));
    }
    ++region_no;
    printf("\n");
  }
}

void Figure::Draw_GeneLoci(DROMPA::Global &p,
			   std::string &pdffilename,
			   int32_t width,
			   int32_t height)
{
  const auto surface = Cairo::PdfSurface::create(pdffilename, width, height);
  int32_t len(p.drawregion.getLenGeneLoci());
  auto &gmp_chr = p.anno.gmp.at(vReadArray.getchr().getname());
  for (auto &m: gmp_chr) {
    if(!p.drawregion.ExistInGeneLociFile(m.second.gname)) continue;

    int32_t start = std::max(0, m.second.txStart - len);
    int32_t end   = std::min(m.second.txEnd + len, vReadArray.getchrlen() -1);
    int32_t num_page(p.drawparam.getNumPage(start, end));
    for(int32_t i=0; i<num_page; ++i) {
      std::cout << boost::format("   page %5d/%5d/%s\r") % (i+1) % num_page % m.second.gname << std::flush;
      PDFPage page(p, vReadArray, vsamplepairoverlayed, surface, start, end);
      page.MakePage(p, i, m.second.gname);
    }
    printf("\n");
  }
}

void Figure::Draw_Whole(DROMPA::Global &p,
			std::string &pdffilename,
			int32_t width,
			int32_t height)
{
#ifdef CAIRO_HAS_PDF_SURFACE
  const auto surface = Cairo::PdfSurface::create(pdffilename, width, height);
  int32_t num_page = p.drawparam.getNumPage(0, vReadArray.getchrlen());
  for (int32_t i=0; i<num_page; ++i) {
    std::cout << boost::format("   page %5d/%5d\r") % (i+1) % num_page << std::flush;
    PDFPage page(p, vReadArray, vsamplepairoverlayed, surface, 0, vReadArray.getchrlen());
    page.MakePage(p, i, "1");
  }
  printf("\n");
#else
  std::cout << "You must compile cairo with PDF support for DROMPA+." << std::endl;
  return;
#endif
}
