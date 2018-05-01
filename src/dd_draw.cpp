/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */

#include <sstream>
#include <iomanip>
#include "dd_draw.hpp"
#include "dd_draw_p.hpp"
#include "SSP/common/inline.hpp"
#include "SSP/common/util.hpp"

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
  void strokeGraph4EachWindow(const Cairo::RefPtr<Cairo::Context> cr,
			      double x_pre, double y_pre, double x_cen, double y_cen, int32_t bottom)
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

}

void DataFrame::StrokeEachBin(const SamplePairParam &pair,
			       const ChrArrayMap &arrays,
			       const int32_t i, const double xcen, const int32_t yaxis,
			       const int32_t nlayer)
{
  double value(getVal(pair, arrays, i));
  if (!value) return;

  int32_t len(getbinlen(value / scale));

  setColor(value, nlayer, 1);
  if(len <0) rel_yline(cr, xcen, yaxis + len, len_binedge);
  cr->stroke();

  if(alpha) {
    setColor(value, nlayer, alpha);
    rel_yline(cr, xcen, yaxis, len);
    cr->stroke();
  }

  return;
}

void LogRatioDataFrame::StrokeEachBin(const SamplePairParam &pair, const ChrArrayMap &arrays,
			       const int32_t i, const double xcen, const int32_t yaxis,
			       const int32_t nlayer)
{
  double value(getVal(pair, arrays, i) / scale);
  if (!value) return;

  int32_t len(0);
  if(value>0) len = -std::min(par.ystep*value, height_df/2);
  else        len = -std::max(par.ystep*value, -height_df/2);

  setColor(value, nlayer, 1);
  rel_yline(cr, xcen, yaxis + len, len_binedge);
  cr->stroke();

  if(alpha) {
    setColor(value, nlayer, alpha);
    rel_yline(cr, xcen, yaxis-height_df/2, len);
    cr->stroke();
  }

  return;
}

void ChIPDataFrame::stroke_peakregion(const SamplePairChr &pair)
{
  //  DEBUGprint("stroke_peakregion");
  cr->set_source_rgba(CLR_RED, 0.4);

  for (auto &bed: pair.peaks1st) {
    if (!my_overlap(bed.start, bed.end, par.xstart, par.xend)) continue;
    cr->set_line_width(bed.length() * par.dot_per_bp);
    double x(par.bp2xaxis(bed.summit - par.xstart));
    rel_yline(cr, x, par.yaxis_now - height_df -5, height_df +10);
  }
  cr->stroke();
  return;
}

void DataFrame::StrokeBinsInLine(const SamplePairParam &pair, const ChrArrayMap &arrays,
			       const int32_t nlayer)
{
  int32_t binsize(pair.getbinsize());
  int32_t sbin(par.xstart/binsize);
  int32_t ebin(par.xend/binsize);
  double dot_per_bin(binsize * par.dot_per_bp);
  int32_t yaxis(par.yaxis_now);  // intに変換

  int32_t thin(std::min(par.width_per_line/(1000*binsize), 20));

  double xcen(par.bp2xaxis(binsize/2)); // initial position
  if (thin > 1) cr->set_line_width(dot_per_bin*thin);
  else cr->set_line_width(dot_per_bin);

  for (int32_t i=sbin; i<ebin; ++i, xcen += dot_per_bin) {
    if (thin > 1 && i%thin) continue;
    StrokeEachBin(pair, arrays, i, xcen, yaxis, nlayer);
  }
  cr->stroke();
  return;
}

void Page::drawInteraction(const InteractionSet &vinter)
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
  
  par.yaxis_now += boxheight +5;

  return;
}

void Page::StrokeReadLines(const DROMPA::Global &p, const SamplePairChr &pair)
{
  if (p.drawparam.showpinter) stroke_readdist<PinterDataFrame>(p, pair);
  if (p.drawparam.showpenrich && pair.pair.first.argvInput!="") stroke_readdist<PenrichDataFrame>(p, pair);
  if (p.drawparam.showratio && pair.pair.first.argvInput!="") {
    if (p.drawparam.showratio == 1) stroke_readdist<RatioDataFrame>(p, pair);
    else if (p.drawparam.showratio == 2) stroke_readdist<LogRatioDataFrame>(p, pair);
  }
  if (p.drawparam.showctag) stroke_readdist<ChIPDataFrame>(p, pair);
  if (p.drawparam.showitag==1 && pair.pair.first.argvInput!="") stroke_readdist<InputDataFrame>(p, pair);
  return;
}

void Page::StrokeGraph(const GraphData &graph)
{
  DEBUGprint("Page::DrawGraph");
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

void Page::StrokeEachLayer(const DROMPA::Global &p)
{  
  if (p.anno.GC.isOn()) StrokeGraph(GC);
  if (p.anno.GD.isOn()) StrokeGraph(GD);
  
  if (p.anno.genefile != "" || p.anno.arsfile != "" || p.anno.terfile != "") DrawGeneAnnotation(p);

  for (size_t j=0; j<pairs.size(); ++j) StrokeReadLines(p, pairs[j]);
  stroke_xaxis_num(par.yaxis_now, 9);
  
  if (p.anno.vinterlist.size()) {
    par.yaxis_now += 5;
    for (auto &x: p.anno.vinterlist) drawInteraction(x);
  }

  //    if(d->bednum) draw_bedfile(d, cr, xstart, xend, chr);
  //if(d->repeat.argv) draw_repeat(d, cr, xstart, xend);
  
  //    if(d->visualize_itag==2) stroke_readdist(p, d, cr, &(sample[0]), xstart, xend, LTYPE_INPUT, 0);
  
  par.yaxis_now += MERGIN_BETWEEN_LINE;
  return;
}

void Page::MakePage(const DROMPA::Global &p, const int32_t page_no, const int32_t region_no)
{
  DEBUGprint("Page::MakePage");
  int32_t line_start, line_end;
  std::tie(line_start, line_end) = get_start_end_linenum(page_no, p.drawparam.getlpp());
 
  par.yaxis_now = OFFSET_Y;
  cr->set_source_rgba(CLR_WHITE, 1);
  cr->paint();

  // Stroke each layer
  for (int i=line_start; i<line_end; ++i) {
    set_xstart_xend(i);
    if (par.xstart >= par.xend) continue;
    StrokeEachLayer(p);
  }

  // Page title
  std::string title;
  if(par.num_page>1) title = chrname + "_" + std::to_string(region_no) + "_" + std::to_string(page_no+1);
  else               title = chrname + "_" + std::to_string(region_no);
  cr->set_source_rgba(CLR_BLACK, 1);
  showtext_cr(cr, 50, 30, title, 16);
  
  cr->show_page();
  return;
}

void Figure::DrawData(DROMPA::Global &p, const chrsize &chr)
{
  DEBUGprint("Figure::DrawData");
  int32_t width(pagewidth);
  int32_t height(p.drawparam.getPageHeight(p, pairs));
  std::string pdffilename(p.getFigFileNameChr(chr.getrefname()));
  //  std::cout << chr.getrefname() << std::endl;
    
#ifdef CAIRO_HAS_PDF_SURFACE
  const auto surface = Cairo::PdfSurface::create(pdffilename, width, height);

  if (!p.drawregion.isRegionBed()){  // whole chromosome
    int32_t num_page(p.drawparam.getNumPage(0, chr.getlen()));
    for (int32_t i=0; i<num_page; ++i) {
      std::cout << boost::format("   page %5d/%5d\r") % (i+1) % num_page << std::flush;
      Page page(p, arrays, pairs, surface, chr, 0, chr.getlen());
      page.MakePage(p, i, 1);
    }
  } else {
    int32_t region_no(1);
    for (auto &x: regionBed) {
      int32_t num_page(p.drawparam.getNumPage(x.start, x.end));
      for(int32_t i=0; i<num_page; ++i) {
	std::cout << boost::format("   page %5d/%5d/%5d\r") % (i+1) % num_page % region_no << std::flush;
	Page page(p, arrays, pairs, surface, chr, x.start, x.end);
	page.MakePage(p, i, region_no);
      }
      ++region_no;
    }
  } 
  std::cout << "Wrote PDF file \"" << pdffilename << "\"" << std::endl;

#else
  std::cout << "You must compile cairo with PDF support for DROMPA." << std::endl;
  return;
#endif
}

int32_t DROMPA::DrawParam::getHeightEachSample(const SamplePairParam &pair) const {
  int32_t height(0);
  int32_t n(0);
  if (showctag)                            { height += getlineheight(); ++n; }
  if (showitag==1 && pair.argvInput != "") { height += getlineheight(); ++n; }
  if (showratio   && pair.argvInput != "") { height += getlineheight(); ++n; }
  if (showpinter)                          { height += getlineheight(); ++n; }
  if (showpenrich && pair.argvInput != "") { height += getlineheight(); ++n; }
  height += MERGIN_BETWEEN_DATA * (n-1);
  
#ifdef DEBUG
  std::cout << "LineHeight: " << getlineheight() << ",n: " << n << std::endl;
  std::cout << "HeightEachSample: " << height << std::endl;
#endif
  return height;
}
    
int32_t DROMPA::DrawParam::getHeightAllSample(const DROMPA::Global &p, const std::vector<SamplePairChr> &pairs) const {
  int32_t height(0);
  for (auto x: pairs) height += getHeightEachSample(x.pair.first);
  height += MERGIN_BETWEEN_DATA * (samplenum-1);
  if (showitag==2) height += getlineheight() + MERGIN_BETWEEN_DATA;

  if (p.anno.GC.isOn()) height += BOXHEIGHT_GRAPH + MERGIN_BETWEEN_GRAPH_DATA;
  if (p.anno.GD.isOn()) height += BOXHEIGHT_GRAPH + MERGIN_BETWEEN_GRAPH_DATA;
  
  if (p.anno.genefile != "" || p.anno.arsfile != "" || p.anno.terfile != "") {
    height += BOXHEIGHT_GENEBOX_EXON + MERGIN_BETWEEN_DATA;
  }
  height += (BOXHEIGHT_INTERACTION +5) * p.anno.vinterlist.size();

#ifdef DEBUG
  std::cout << "HeightAllSample; " << height << std::endl;
#endif
  return height;
}
    

int32_t DROMPA::DrawParam::getPageHeight(const DROMPA::Global &p, const std::vector<SamplePairChr> &pairs) const {
  int32_t height(OFFSET_Y*2);
  height += getHeightAllSample(p, pairs) * linenum_per_page;
  height += MERGIN_BETWEEN_LINE * (linenum_per_page-1);
  height += 50; // 保険

#ifdef DEBUG
  std::cout << "PageHeight: " << height << std::endl;
#endif
  return height;
}

