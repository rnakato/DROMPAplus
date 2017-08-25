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

enum {GFTYPE_REFFLAT=0, GFTYPE_GTF=1, GFTYPE_SGD=2};

void ChIPDataFrame::stroke_bin(const SamplePairChr &pair,
			       const std::unordered_map<std::string, ChrArray> &arrays,
			       const int32_t i, const double xcen, const int32_t yaxis, const int32_t viz)
{
  double value((arrays.at(pair.argvChIP)).array[i] / scale);
  if (!value) return;

  int32_t len(getbinlen(value));
  cr->set_source_rgba(CLR_GREEN3, 1);
  if (!viz) {
    rel_yline(cr, xcen, yaxis, len);
  } else {
    if(len <=-1 && value <= par.barnum) rel_yline(cr, xcen, yaxis + len, 1);
    cr->stroke();
    cr->set_source_rgba(CLR_GREEN3, 0.5);
    rel_yline(cr, xcen, yaxis, len);     
    cr->stroke();
  }
  return;
}

void InputDataFrame::stroke_bin(const SamplePairChr &pair,
				const std::unordered_map<std::string, ChrArray> &arrays,
				const int32_t i, const double xcen, const int32_t yaxis, const int32_t viz)
{
  double value((arrays.at(pair.argvInput)).array[i] / scale);
  if (!value) return;

  int32_t len(getbinlen(value));

  cr->set_source_rgba(CLR_BLUE, 1);
  rel_yline(cr, xcen, yaxis, len);
  return;
}

void RatioDataFrame::stroke_bin(const SamplePairChr &pair,
				const std::unordered_map<std::string, ChrArray> &arrays,
				const int32_t i, const double xcen, const int32_t yaxis, const int32_t viz)
{
  double value(CALCRATIO(arrays.at(pair.argvChIP).array[i], arrays.at(pair.argvInput).array[i], 1) / scale);  // TOTAL READ 正規化入れる
  if (!value) return;
 
  int32_t len(getbinlen(value));

  cr->set_source_rgba(CLR_ORANGE, 1);
  rel_yline(cr, xcen, yaxis, len);
  return;
}

void LogRatioDataFrame::stroke_bin(const SamplePairChr &pair,
				   const std::unordered_map<std::string, ChrArray> &arrays,
				   const int32_t i, const double xcen, const int32_t yaxis, const int32_t viz)
{
  double data(CALCRATIO(arrays.at(pair.argvChIP).array[i], arrays.at(pair.argvInput).array[i], 1));  // TOTAL READ 正規化入れる
  if (!data) return;

  double value(log10(data)/scale);
  int32_t len(0);
    
  if(value>0) len = -std::min(par.ystep*value, height_df/2);
  else        len = -std::max(par.ystep*value, -height_df/2);
    
  cr->set_source_rgba(CLR_ORANGE, 1);
  rel_yline(cr, xcen, yaxis-height_df/2, len);
  return;
}

void PinterDataFrame::stroke_bin(const SamplePairChr &pair,
				 const std::unordered_map<std::string, ChrArray> &arrays,
				 const int32_t i, const double xcen, const int32_t yaxis, const int32_t viz)
{
  // double data = zero_inflated_binomial_test(WIGARRAY2VALUE(sample->ChIP->data[i]), sample->ChIP->nb_p, sample->ChIP->nb_n);
  double data = CALCRATIO(arrays.at(pair.argvChIP).array[i], arrays.at(pair.argvInput).array[i], 1);
  double value(data/scale);
  if (!value) return;
 
  int32_t len(getbinlen(value));

  //    if(data > p->pthre_internal) cr->set_source_rgba(CLR_PINK, 1);
  //else cr->set_source_rgba(CLR_CLR_GRAY, 1);
  cr->set_source_rgba(CLR_PINK, 1);
  rel_yline(cr, xcen, yaxis, len);
  return;
}

void PenrichDataFrame::stroke_bin(const SamplePairChr &pair,
				  const std::unordered_map<std::string, ChrArray> &arrays,
				  const int32_t i, const double xcen, const int32_t yaxis, const int32_t viz)
{
  //double data = binomial_test(WIGARRAY2VALUE(sample->ChIP->data[i]), WIGARRAY2VALUE(sample->Input->data[i]), sample->comp->genome->ratio);
  double data = CALCRATIO(arrays.at(pair.argvChIP).array[i], arrays.at(pair.argvInput).array[i], 1);
  double value(data/scale);
  if (!value) return;
 
  int32_t len(getbinlen(value));

  //    if(data > p->pthre_enrich) cr->set_source_rgba(CLR_PINK, 1);
  //else cr->set_source_rgba(CLR_CLR_GRAY, 1);
  cr->set_source_rgba(CLR_PINK, 1);
  rel_yline(cr, xcen, yaxis, len);
  return;
}


void DataFrame::stroke_dataframe(const DROMPA::Global &p, const int32_t nlayer) {
  DEBUGprint("stroke_dataframe");
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

void DataFrame::stroke_bindata(const DROMPA::Global &p, const SamplePairChr &pair,
			       const std::unordered_map<std::string, ChrArray> &arrays, const int32_t nlayer)
{
  DEBUGprint("stroke_bindata");
  int32_t binsize(pair.getbinsize());
  int32_t sbin(par.xstart/binsize);
  int32_t ebin(par.xend/binsize);
  double dot_per_bin(binsize * dot_per_bp);
  int32_t yaxis(par.yaxis_now);  // intに変換

  int32_t thin(std::min(par.width_per_line/(1000*binsize), 20));

  // peak region 色を変えるのではなく、その領域のバックを光らせる
  double xcen(bp2xaxis(binsize/2)); // initial position
  if(thin > 1) cr->set_line_width(dot_per_bin*thin);
  else cr->set_line_width(dot_per_bin);

  for (int32_t i=sbin; i<ebin; ++i, xcen += dot_per_bin) {
    if (thin > 1 && i%thin) continue;
    stroke_bin(pair, arrays, i, xcen, yaxis, p.drawparam.viz);
  }
  cr->stroke();
  return;
}

void Page::stroke_each_layer(const DROMPA::Global &p, const SamplePairChr &pair, const int32_t nlayer)
{
  DEBUGprint("stroke_each_layer");
  if (p.drawparam.showpinter) stroke_readdist<PinterDataFrame>(p, pair, nlayer);
  if (p.drawparam.showpenrich && pair.argvInput!="") stroke_readdist<PenrichDataFrame>(p, pair, nlayer);
  if (p.drawparam.showratio && pair.argvInput!="") {
    //    if(p->ftype==FTYPE_GV)  stroke_readdist(p, d, cr, sample, xstart, xend, LTYPE_RATIO_GV, nlayer);
    //else
    if (p.drawparam.showratio == 1) stroke_readdist<RatioDataFrame>(p, pair, nlayer);
    else if (p.drawparam.showratio == 2) stroke_readdist<LogRatioDataFrame>(p, pair, nlayer);
  }
  if (p.drawparam.showctag) stroke_readdist<ChIPDataFrame>(p, pair, nlayer);
  if (p.drawparam.showitag==1 && pair.argvInput!="") stroke_readdist<InputDataFrame>(p, pair, nlayer);
  return;
}


/*void Page::draw_graph(DDParam *d, cairo_t *cr, Graph *graph, gint memnum, gint boxheight, bool color, bool xaxis)
{
  gint i;
  gchar str[16];
  gint s = xstart/graph->wsize, e = xend/graph->wsize +1;
  gdouble xcen, xpre, ycen, ypre;
  gdouble width = (xend-xstart+1) * dot_per_bp;
  gdouble diff = graph->wsize * dot_per_bp;

  yaxis_now += boxheight/2;

  // graph line
  if(color==false) cairo_set_source_rgba(cr, CLR_GREEN, 1); else cairo_set_source_rgba(cr, CLR_BLUE, 1);
  cairo_set_line_width(cr, 0.6);
  xpre = BP2XAXIS(0);
  xcen = BP2XAXIS(0.5*graph->wsize);
  ypre = defy_graph(graph->array[s], boxheight, graph->mmin, graph->mmax);
  for(i=s; i<e; i++, xcen += diff){
    ycen = defy_graph(graph->array[i], boxheight, graph->mmin, graph->mmax);
    check_bottom(cr, xpre, ypre, xcen, ycen, yaxis_now + boxheight/2 + 10);
    xpre = xcen;
    ypre = ycen;
  }
  // keys
  cairo_set_source_rgba(cr, CLR_BLACK, 1);
  showtext_cr(cr, OFFSET_X - 5*strlen(graph->name)-55, yaxis_now, graph->name, 13);

  yaxis_now += boxheight/2;

  // axis
  cairo_set_source_rgba(cr, CLR_BLACK, 1);
  cairo_set_line_width(cr, 0.4);
  rel_yline(cr, OFFSET_X, yaxis_now - boxheight, boxheight);
  cairo_stroke(cr);
  cairo_set_line_width(cr, 1.5);
  rel_xline(cr, OFFSET_X, yaxis_now, width);
  cairo_stroke(cr);
  stroke_xaxis(d, cr, xstart, xend);
  
  // memory 
  gdouble x, y;
  gdouble mem = (graph->mmax - graph->mmin)/memnum;
  //  gdouble memnum = (graph->mmax - graph->mmin)/graph->mem;
  cairo_set_source_rgba(cr, CLR_BLACK, 1);
  cairo_set_line_width(cr, 0.5);
  for(i=0; i<=memnum; i++){
    if(mem <1) sprintf(str, "%.2f", graph->mmin + i*mem);
    else       sprintf(str, "%d", (int)(graph->mmin + i*mem));
    x = OFFSET_X - 5*strlen(str) - 7;
    y = yaxis_now - i*boxheight/memnum;
    showtext_cr(cr, x, y+2, str, 9);
    rel_xline(cr, OFFSET_X-2, y, 2);
    cairo_stroke(cr);
  }
  if(xaxis==true) stroke_xaxis_num(d, cr, xstart, xend, yaxis_now, 9);

  return;
  }*/

void show_colors(const Cairo::RefPtr<Cairo::Context> cr, const int32_t x, int32_t &ycen,
		 const std::string &label, const double r, const double g, const double b)
{
  cr->set_source_rgba(r,g,b, 1);
  rel_xline(cr, x, ycen, 20);
  cr->stroke();
  cr->set_source_rgba(CLR_BLACK, 1);
  showtext_cr(cr, x+24, ycen+4, label, 10);
  ycen += 15;
  return;
}

void ShowGeneAnnotation(const Cairo::RefPtr<Cairo::Context> cr, const double yaxis, const int32_t gftype)
{
  int32_t x(50);
  int32_t y(yaxis-20);

  cr->set_line_width(2.5);
  show_colors(cr, x, y, "Coding", CLR_BLUE);
  show_colors(cr, x, y, "Noncoding", CLR_GREEN);
  if (gftype == GFTYPE_SGD) {
    show_colors(cr, x, y, "rRNA", CLR_BLACK);
    show_colors(cr, x, y, "LTR", CLR_PURPLE);
  } else {
    show_colors(cr, x, y, "Processed transcript", CLR_ORANGE);
    show_colors(cr, x, y, "MicroRNA", CLR_PINK);
    show_colors(cr, x, y, "Pseudo", CLR_GRAY2);
    show_colors(cr, x, y, "Others", CLR_BLACK);
  }
  return;
}

void Page::strokeGeneSGD(const DROMPA::Global &p, const double ycenter)
{
  DEBUGprint("strokeGeneSGD");
  
  ShowGeneAnnotation(cr, ycenter, p.anno.getgftype());

  try {
    int32_t ars_on(0);
    int32_t on_plus(0);
    int32_t on_minus(0);
    for(auto &m: p.anno.gmp.at(rmchr(chrname))) {
      if(!my_overlap(m.second.txStart, m.second.txEnd, par.xstart, par.xend)) continue;

      cr->set_line_width(0.3);
      GeneElement g(m, par.xstart, ycenter, 0, on_plus, on_minus);

      if (m.second.gtype=="centromere" || m.second.gtype=="teromere") {
	cr->set_source_rgba(CLR_GREEN, 1);
	rel_yline(cr, g.xcen, ycenter -2, g.ylen);
	showtext_cr(cr, g.x_name, g.y_name-6, m.second.gname, 8);
      }
      else if (m.second.gtype=="ARS") {
	cr->set_source_rgba(CLR_RED, 1);
	rel_yline(cr, g.xcen, ycenter -2, g.ylen +10 - ars_on*8);
	showtext_cr(cr, g.x_name, g.y_name +4 - ars_on*8, m.second.gname, 7);
	if(ars_on==2) ars_on=0; else ++ars_on;
      }
      else if (m.second.gtype=="TER") {
	cr->set_source_rgba(CLR_OLIVE, 1);
	rel_yline(cr, g.xcen, ycenter -2, g.ylen -5);
	showtext_cr(cr, g.x_name, g.y_name-11, m.second.gname, 7);
      }
      else if (m.second.gtype=="rRNA" || m.second.gtype=="snoRNA") {
	cr->set_source_rgba(CLR_BLACK, 1);
	showtext_cr(cr, g.x_name, g.y_name, m.second.gname, 6);
	cr->set_line_width(1.5);
	rel_xline(cr, g.x1, g.ybar, g.xwid);
	cr->stroke();
      }
      else if (m.second.gtype=="LTR" || m.second.gtype=="retrotransposon" || isStr(m.second.gtype, "repeat")) {
	cr->set_source_rgba(CLR_PURPLE, 1);
	showtext_cr(cr, g.x_name, g.y_name, m.second.gname, 6);
	cr->set_line_width(1.5);
	rel_xline(cr, g.x1, g.ybar, g.xwid);
	cr->stroke();
      }
      else if (m.second.gtype=="tRNA") {
	cr->set_source_rgba(CLR_GREEN, 1);
	showtext_cr(cr, g.x_name, g.y_name, m.second.gname, 6);
	cr->set_line_width(1.5);
	rel_xline(cr, g.x1, g.ybar, g.xwid);
	cr->stroke();
      }
      else {
	cr->set_source_rgba(CLR_BLUE, 1);
	showtext_cr(cr, g.x_name, g.y_name, m.second.gname, 6);
	cr->set_line_width(1.5);
	rel_xline(cr, g.x1, g.ybar, g.xwid);
	cr->stroke();
      }
    }
  } catch (...) {
    std::cerr << "Warning: " << chrname  << " has no gene." << std::endl;
  }
  
  return;
}

void Page::strokeARS(const HashOfGeneDataMap &mp, const double ycenter)
{
  int32_t ars_on(0);
  cr->set_line_width(0.3);
  int32_t on_plus(0);
  int32_t on_minus(0);

  try {
    for(auto &m: mp.at(rmchr(chrname))) {
      if(!my_overlap(m.second.txStart, m.second.txEnd, par.xstart, par.xend)) continue;

      GeneElement g(m, par.xstart, ycenter, 0, on_plus, on_minus);

      if (m.second.gtype=="ARS") {
	cr->set_source_rgba(CLR_RED, 1);
	rel_yline(cr, g.xcen, ycenter -2, g.ylen +14 - 8 * ars_on);
	showtext_cr(cr, g.x_name, g.y_name +8 - ars_on*8, m.second.gname, 8);
	if(ars_on==2) ars_on=0; else ++ars_on;
      }
      else if (m.second.gtype=="centromere") {
	cr->set_source_rgba(CLR_GREEN, 1);
	rel_yline(cr, g.xcen, ycenter -2, g.ylen -2);
	showtext_cr(cr, g.x_name, g.y_name -8, m.second.gname, 8);
      }
      else if (m.second.gtype=="teromere") {
	cr->set_source_rgba(CLR_OLIVE, 1);
	rel_yline(cr, g.xcen, ycenter -2, g.ylen -5);
	showtext_cr(cr, g.x_name, g.y_name -11, m.second.gname, 7);
      }else continue;
    }
  } catch (...) {
    std::cerr << "Warning: " << chrname  << " has no gene." << std::endl;
  }
  
  return;
}

void Page::strokeGene(const DROMPA::Global &p, const double ycenter)
{
  DEBUGprint("strokeGene");
  
  ShowGeneAnnotation(cr, ycenter, p.anno.getgftype());

  try {
    int32_t on_plus(0);
    int32_t on_minus(0);
    
    for(auto &m: p.anno.gmp.at(rmchr(chrname))) {
      if(!my_overlap(m.second.txStart, m.second.txEnd, par.xstart, par.xend)) continue;

      GeneElement g(m, par.xstart, ycenter, 1, on_plus, on_minus);

      if(isStr(m.second.gtype, "coding"))         cr->set_source_rgba(CLR_BLUE, 1);  
      else if(isStr(m.second.gtype, "noncoding")) cr->set_source_rgba(CLR_GREEN, 1);
      else if(isStr(m.second.gtype, "miRNA"))     cr->set_source_rgba(CLR_PINK, 1);
      else if(isStr(m.second.gtype, "process"))   cr->set_source_rgba(CLR_ORANGE, 1);
      else if(isStr(m.second.gtype, "pseudo"))    cr->set_source_rgba(CLR_GRAY2, 1);
      else                                        cr->set_source_rgba(CLR_BLACK, 1);

      // gene body 
      cr->set_line_width(1.5);
      rel_yline(cr, g.x1, g.ybar-4, 8);
      rel_yline(cr, g.x2, g.ybar-4, 8);
      cr->stroke();
      cr->set_line_width(3);
      rel_xline(cr, g.x1, g.ybar, g.xwid);
      cr->stroke();
      // exon
      cr->set_line_width(6);
      for (int32_t i=0; i<m.second.exonCount; ++i) {
	double x(bp2xaxis(m.second.exon[i].start - par.xstart +1));
	double xlen(m.second.exon[i].getlen() * dot_per_bp);
	rel_xline(cr, x, g.ybar, xlen);
	cr->stroke();
      }
      // name
      cr->set_source_rgba(CLR_BLACK, 1);
      showtext_cr(cr, g.x_name, g.y_name, m.second.gname, 10);

    }
  } catch (...) {
    std::cerr << "Warning: " << chrname  << " has no gene." << std::endl;
  }
  return;
}

void Page::DrawAnnotation(const DROMPA::Global &p)
{
  DEBUGprint("DrawAnnotation");
  int32_t boxheight;
  if(p.anno.getgftype() == GFTYPE_SGD) boxheight = BOXHEIGHT_GENEBOX_NOEXON;
  else boxheight = BOXHEIGHT_GENEBOX_EXON;

  double ytop(par.yaxis_now);
  double ycenter(ytop + boxheight/2);

  if (p.anno.showars) {
    strokeARS(p.anno.gmp, ycenter);
    showtext_cr(cr, 70, ycenter, "ARS", 12);
  } else {
    if(p.anno.getgftype() == GFTYPE_SGD) strokeGeneSGD(p, ycenter);
    else strokeGene(p, ycenter);
  }
  /* frame */
  cr->rectangle(par.xstart, par.xend, ytop, boxheight);
  //  if(d->backcolors) cr->fill(); // fill_rectangle(cr, xstart, xend, ytop, boxheight, CLR_YELLOW2, 0.1);

  /* genome line */
  cr->set_source_rgba(CLR_BLACK, 1);
  cr->set_line_width(1.5);
  rel_xline(cr, OFFSET_X, ycenter, par.getXaxisLen());
  cr->stroke();

  /* memory */
  stroke_xaxis(ycenter);
  stroke_xaxis_num(ycenter, 9);

  par.yaxis_now += boxheight + MERGIN_BETWEEN_DATA;
  return;
}

void Page::Draw(const DROMPA::Global &p, const int32_t page_curr, const int32_t region_no)
{
  DEBUGprint("Page::Draw");
  int32_t line_start, line_end;
  std::tie(line_start, line_end) = get_start_end_linenum(page_curr, p.drawparam.getlpp());
 
  par.yaxis_now = OFFSET_Y;
  cr->set_source_rgba(CLR_WHITE, 1);
  cr->paint();

  for(int i=line_start; i<line_end; ++i) {
    set_xstart_xend(i);
    if(par.xstart >= par.xend) continue;
      /*   if(d->GC.argv){
	    draw_graph(d, cr, &(d->GC), xstart, xend, memnum_GC, boxheight_graph, false, true);
	    yaxis_now += mergin_between_graph_data;
	    }
	    if(d->GD.argv){
	    draw_graph(d, cr, &(d->GD), xstart, xend, memnum_GC, boxheight_graph, false, true);
	    yaxis_now += mergin_between_graph_data;
	    }*/

    if(p.anno.genefile != "" || p.anno.arsfile != "" || p.anno.terfile != "") DrawAnnotation(p);
     /*   if(d->internum){
      for(j=0; j<d->internum; ++j) draw_interaction(cr, &(d->inter[j]), xstart, xend, chr);
      }*/
    //    double ytemp = par.yaxis_now;
    int32_t nlayer = 0;
    for(size_t j=0; j<pairs.size(); ++j) stroke_each_layer(p, pairs[j], nlayer);
    stroke_xaxis_num(par.yaxis_now, 9);

    //    if(d->bednum) draw_bedfile(d, cr, xstart, xend, chr);
    //if(d->repeat.argv) draw_repeat(d, cr, xstart, xend);
      
    /*    if(p->samplenum_overlay){
      nlayer = 1;
      yaxis_now = ytemp;
      for(; j<p->samplenum; ++j) stroke_each_layer(p, d, &(sample[j]), cr, xstart, xend, nlayer);
      }*/
    //    if(d->visualize_itag==2) stroke_readdist(p, d, cr, &(sample[0]), xstart, xend, LTYPE_INPUT, 0);
      
    par.yaxis_now += MERGIN_BETWEEN_LINE;
  }

  // page titile
  std::string title;
  if(par.num_page>1) title = chrname + "_" + std::to_string(region_no) + "_" + std::to_string(page_curr+1);
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
  std::string pdffile(p.getFigFileNameChr(chr.getrefname()));
  std::cout << chr.getrefname() << std::endl;
    
#ifdef CAIRO_HAS_PDF_SURFACE
  const auto surface = Cairo::PdfSurface::create(pdffile, width, height);

  if(!p.drawregion.isRegionBed()){  // whole chromosome
    int32_t num_page(p.drawparam.getNumPage(0, chr.getlen()));
    for(int32_t i=0; i<num_page; ++i) {
      std::cout << boost::format("   page %5d/%5d\r") % (i+1) % num_page << std::flush;
      Page page(p, arrays, pairs, surface, chr.getrefname(), 0, chr.getlen());
      page.Draw(p, i, 1);
    }
  }else{
    int32_t region_no(1);
    for (auto &x: regionBed) {
      int32_t num_page(p.drawparam.getNumPage(x.start, x.end));
      for(int32_t i=0; i<num_page; ++i) {
	std::cout << boost::format("   page %5d/%5d/%5d\r") % (i+1) % num_page % region_no << std::flush;
	Page page(p, arrays, pairs, surface, chr.getrefname(), x.start, x.end);
	page.Draw(p, i, region_no);
      }
      ++region_no;
    }
  } 
  std::cout << "Wrote PDF file \"" << pdffile << "\"" << std::endl;

#else
  std::cout << "You must compile cairo with PDF support for DROMPA." << std::endl;
  return;
#endif
}

int32_t DROMPA::DrawParam::getHeightEachSample(const SamplePairChr &pair) const {
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
  for(auto x: pairs) height += getHeightEachSample(x);
  height += MERGIN_BETWEEN_DATA * (samplenum-1);
  if(showitag==2) height += getlineheight() + MERGIN_BETWEEN_DATA;

  if(p.anno.genefile != "" || p.anno.arsfile != "" || p.anno.terfile != "") {
    if(p.anno.getgftype() == GFTYPE_SGD) height += BOXHEIGHT_GENEBOX_NOEXON;
    else height += BOXHEIGHT_GENEBOX_EXON;
    height += MERGIN_BETWEEN_DATA;
  }
  //  height_lpp += BOXHEIGHT_INTERACTION * d->internum;

#ifdef DEBUG
  std::cout << "HeightAllSample; " << height << std::endl;
#endif
  return height;
}
    

int32_t DROMPA::DrawParam::getPageHeight(const DROMPA::Global &p, const std::vector<SamplePairChr> &pairs) const {
  int32_t height(OFFSET_Y*2);
  height += getHeightAllSample(p, pairs) * linenum_per_page;
  height += MERGIN_BETWEEN_LINE * (linenum_per_page-1);

#ifdef DEBUG
  std::cout << "PageHeight: " << height << std::endl;
#endif
  return height;
}

