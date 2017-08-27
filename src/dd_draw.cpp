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
  double value(CALCRATIO(arrays.at(pair.argvChIP).array[i], arrays.at(pair.argvInput).array[i], pair.ratio) / scale);  // TOTAL READ 正規化入れる
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
  double data(CALCRATIO(arrays.at(pair.argvChIP).array[i], arrays.at(pair.argvInput).array[i], pair.ratio));  // TOTAL READ 正規化入れる
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
  //  DEBUGprint("stroke_dataframe");
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
  //  DEBUGprint("stroke_bindata");
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
  //  DEBUGprint("stroke_each_layer");
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

std::vector<genedata> get_garray(const GeneDataMap &mp, const int32_t xstart, const int32_t xend)
{
  std::vector<genedata> garray;
  for (auto &m: mp) {
    if (!my_overlap(m.second.txStart, m.second.txEnd, xstart, xend)) continue;
    if (m.second.gtype == "nonsense_mediated_decay" ||
	m.second.gtype == "processed_transcript" ||
	m.second.gtype == "retained_intron") continue;
    garray.emplace_back(m.second);
  }
  sort(garray.begin(), garray.end(), 
       [](const genedata &x, const genedata &y) { return x.txStart < y.txStart;});
  return garray;
}

void Page::strokeGeneSGD(const DROMPA::Global &p, const double ycenter)
{
  DEBUGprint("strokeGeneSGD");

  cr->set_line_width(2.5);
  int32_t ycen(ycenter-30);
  show_colors(cr, 50, ycen, "Coding",    CLR_BLUE);
  show_colors(cr, 50, ycen, "Noncoding", CLR_GREEN);
  show_colors(cr, 50, ycen, "rRNA",      CLR_BLACK);
  show_colors(cr, 50, ycen, "LTR",       CLR_PURPLE);

  try {
    std::vector<genedata> garray(get_garray(p.anno.gmp.at(rmchr(chrname)), par.xstart, par.xend));

    int32_t ars_on(0);
    int32_t on_plus(0);
    int32_t on_minus(0);
    for (auto &m: garray) {
      GeneElement g(m, par.xstart, ycenter, 0, on_plus, on_minus);
 
      cr->set_line_width(0.3);
      if (m.gtype=="centromere" || m.gtype=="teromere") {
	cr->set_source_rgba(CLR_GREEN, 1);
	rel_yline(cr, g.xcen, ycenter -2, g.ylen);
	showtext_cr(cr, g.x_name, g.y_name-6, m.gname, 8);
      }
      else if (m.gtype=="ARS") {
	cr->set_source_rgba(CLR_RED, 1);
	rel_yline(cr, g.xcen, ycenter -2, g.ylen +10 - ars_on*8);
	showtext_cr(cr, g.x_name, g.y_name +4 - ars_on*8, m.gname, 7);
	if(ars_on==2) ars_on=0; else ++ars_on;
      }
      else if (m.gtype=="TER") {
	cr->set_source_rgba(CLR_OLIVE, 1);
	rel_yline(cr, g.xcen, ycenter -2, g.ylen -5);
	showtext_cr(cr, g.x_name, g.y_name-11, m.gname, 7);
      }
      else if (m.gtype=="rRNA" || m.gtype=="snoRNA") {
	cr->set_source_rgba(CLR_BLACK, 1);
	showtext_cr(cr, g.x_name, g.y_name, m.gname, 6);
	cr->set_line_width(1.5);
	rel_xline(cr, g.x1, g.ybar, g.xwid);
	cr->stroke();
      }
      else if (m.gtype=="LTR" || m.gtype=="retrotransposon" || isStr(m.gtype, "repeat")) {
	cr->set_source_rgba(CLR_PURPLE, 1);
	showtext_cr(cr, g.x_name, g.y_name, m.gname, 6);
	cr->set_line_width(1.5);
	rel_xline(cr, g.x1, g.ybar, g.xwid);
	cr->stroke();
      }
      else if (m.gtype=="tRNA") {
	cr->set_source_rgba(CLR_GREEN, 1);
	showtext_cr(cr, g.x_name, g.y_name, m.gname, 6);
	cr->set_line_width(1.5);
	rel_xline(cr, g.x1, g.ybar, g.xwid);
	cr->stroke();
      }
      else {
	cr->set_source_rgba(CLR_BLUE, 1);
	showtext_cr(cr, g.x_name, g.y_name, m.gname, 6);
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
  cr->set_line_width(0.3);
  try {
    std::vector<genedata> garray(get_garray(mp.at(rmchr(chrname)), par.xstart, par.xend));

    int32_t ars_on(0);
    int32_t on_plus(0);
    int32_t on_minus(0);
    for (auto &m: garray) {
      GeneElement g(m, par.xstart, ycenter, 0, on_plus, on_minus);

      if (m.gtype=="ARS") {
	cr->set_source_rgba(CLR_RED, 1);
	rel_yline(cr, g.xcen, ycenter -2, g.ylen +14 - 8 * ars_on);
	showtext_cr(cr, g.x_name, g.y_name +8 - ars_on*8, m.gname, 8);
	if(ars_on==2) ars_on=0; else ++ars_on;
      }
      else if (m.gtype=="centromere") {
	cr->set_source_rgba(CLR_GREEN, 1);
	rel_yline(cr, g.xcen, ycenter -2, g.ylen -2);
	showtext_cr(cr, g.x_name, g.y_name -8, m.gname, 8);
      }
      else if (m.gtype=="teromere") {
	cr->set_source_rgba(CLR_OLIVE, 1);
	rel_yline(cr, g.xcen, ycenter -2, g.ylen -5);
	showtext_cr(cr, g.x_name, g.y_name -11, m.gname, 7);
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
  
  cr->set_line_width(2.5);
  int32_t ycen(ycenter-20);
  show_colors(cr, 50, ycen, "CodingRNA", CLR_BLUE);
  show_colors(cr, 50, ycen, "lincRNA", CLR_PINK);
  show_colors(cr, 50, ycen, "Antisence", CLR_GREEN);
  show_colors(cr, 50, ycen, "other RNA", CLR_ORANGE);
  show_colors(cr, 50, ycen, "Pseudo", CLR_GRAY2);
  show_colors(cr, 50, ycen, "Others", CLR_BLACK);

  try {
    std::vector<genedata> garray(get_garray(p.anno.gmp.at(rmchr(chrname)), par.xstart, par.xend));

    double llimit(150);
    double rlimit(OFFSET_X + width_draw + 60);
    int32_t on_plus(0);
    int32_t on_minus(0);
    for (auto &m: garray) {
      GeneElement g(m, par.xstart, ycenter, 1, on_plus, on_minus);

      if (isStr(m.gtype, "protein_coding")) cr->set_source_rgba(CLR_BLUE, 1);
      else if (isStr(m.gtype, "lincRNA"))   cr->set_source_rgba(CLR_PINK, 1);
      else if (isStr(m.gtype, "antisense")) cr->set_source_rgba(CLR_GREEN, 1);
      else if (isStr(m.gtype, "RNA"))       cr->set_source_rgba(CLR_ORANGE, 1);
      else if (isStr(m.gtype, "pseudo"))    cr->set_source_rgba(CLR_GRAY2, 1);
      else                                  cr->set_source_rgba(CLR_BLACK, 1);

      // gene body
      cr->set_line_width(1.5);
      if(g.x1 >= llimit) rel_yline(cr, g.x1, g.ybar-4, 8);
      if(g.x2 <= rlimit) rel_yline(cr, g.x2, g.ybar-4, 8);
      cr->stroke();
      cr->set_line_width(3);
      cr->move_to(std::max(g.x1,llimit), g.ybar);
      cr->line_to(std::min(g.x2,rlimit), g.ybar);
      cr->stroke();
      
      // exon
      cr->set_line_width(6);
      for (int32_t i=0; i<m.exonCount; ++i) {
	double x(bp2xaxis(m.exon[i].start - par.xstart +1));
	double xlen(std::max(1.0,m.exon[i].getlen() * dot_per_bp));
	if (x >= llimit && x <= rlimit) {
	  rel_xline(cr, x, g.ybar, xlen);
	  cr->stroke();
	}
      }
      
      // name
      cr->set_source_rgba(CLR_BLACK, 1);
      if(p.anno.showtranscriptname) showtext_cr(cr, g.x_name, g.y_name, m.tname, 7);
      else showtext_cr(cr, g.x_name, g.y_name, m.gname, 7);
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
  //  stroke_xaxis_num(ycenter, 9);

  par.yaxis_now += boxheight + MERGIN_BETWEEN_DATA;
  return;
}

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

void Page::drawGraph(const GraphData &graph)
{
  DEBUGprint("Page::DrawGraph");
  int32_t s(par.xstart/graph.binsize);
  int32_t e(par.xend/graph.binsize +1);
  double diff = graph.binsize * dot_per_bp;

  double ytop(par.yaxis_now);
  double ycenter(par.yaxis_now + graph.boxheight/2);
  double ybottom(par.yaxis_now + graph.boxheight);

  // graph line
  cr->set_line_width(0.6);
  cr->set_source_rgba(CLR_GREEN, 1);
  double xpre(OFFSET_X);
  double xcen(bp2xaxis(0.5*graph.binsize));
  double ypre(ybottom - graph.getylen(s));
  for (int32_t i=s; i<e; ++i, xcen += diff) {
    double ycen(ybottom - graph.getylen(i));
    strokeGraph4EachWindow(cr, xpre, ypre, xcen, ycen, ybottom + 10);
    xpre = xcen;
    ypre = ycen;
  }
  // label
  cr->set_source_rgba(CLR_BLACK, 1);
  showtext_cr(cr, OFFSET_X - 5*graph.label.length() -55, ycenter, graph.label, 13);

  // y-axis
  cr->set_source_rgba(CLR_BLACK, 1);
  cr->set_line_width(0.4);
  rel_yline(cr, OFFSET_X, ytop, graph.boxheight);
  cr->stroke();
  cr->set_line_width(1.5);
  rel_xline(cr, OFFSET_X, ybottom, (par.xend - par.xstart+1) * dot_per_bp);
  cr->stroke();
  // x-axis
  stroke_xaxis(ybottom);
  
  // y memory 
  cr->set_source_rgba(CLR_BLACK, 1);
  cr->set_line_width(0.5);
  for (int32_t i=0; i<=graph.memnum; ++i) {
    std::string str(graph.getmemory(i));
    double x(OFFSET_X - 5*str.length() - 7);
    double y(ybottom - i*graph.getBoxHeight4mem());
    showtext_cr(cr, x, y+2, str, 9);
    rel_xline(cr, OFFSET_X-2, y, 2);
    cr->stroke();
  }
  //  if(xaxis==true) stroke_xaxis_num(d, cr, xstart, xend, yaxis_now, 9);

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

    if(p.anno.GC.isOn()){
      drawGraph(GC);
      par.yaxis_now += GC.boxheight + mergin_between_graph_data;
    }
    if(p.anno.GD.isOn()){
      drawGraph(GD);
      par.yaxis_now += GD.boxheight + mergin_between_graph_data;
    }

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
      Page page(p, arrays, pairs, surface, chr, 0, chr.getlen());
      page.Draw(p, i, 1);
    }
  }else{
    int32_t region_no(1);
    for (auto &x: regionBed) {
      int32_t num_page(p.drawparam.getNumPage(x.start, x.end));
      for(int32_t i=0; i<num_page; ++i) {
	std::cout << boost::format("   page %5d/%5d/%5d\r") % (i+1) % num_page % region_no << std::flush;
	Page page(p, arrays, pairs, surface, chr, x.start, x.end);
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

  if(p.anno.GC.isOn()) height += BOXHEIGHT_GRAPH + mergin_between_graph_data;
  if(p.anno.GD.isOn()) height += BOXHEIGHT_GRAPH + mergin_between_graph_data;
  
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

