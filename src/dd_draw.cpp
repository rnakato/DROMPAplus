/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */

#include <sstream>
#include <iomanip>
#include "dd_draw.hpp"
#include "dd_draw_p.hpp"
#include "SSP/common/inline.hpp"

      /*    if (p->gapfile && d->gaparray[i] >= GAP_THRE){
	    cr->set_source_rgba(CLR_BLACK, 0.3);
	    rel_yline(cr, xcen, yaxis - height_df, height_df);
	    cr->stroke();
	    } else if(p->mpfile && d->mparray[i] < p->mpthre){
	    cr->set_source_rgba(CLR_BLACK, 0.3);
	    rel_yline(cr, xcen, yaxis - height_df, height_df);
	    cr->stroke();
	    }*/


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
  int32_t binsize(pair.getbinsize());
  int32_t sbin(par.xstart/binsize);
  int32_t ebin(par.xend/binsize);
  double dot_per_bin(binsize * par.dot_per_bp);
  int32_t yaxis(par.yaxis_now);  // intに変換

  int32_t thin(std::min(par.width_per_line/(1000*binsize), 20));

  // peak region 色を変えるのではなく、その領域のバックを光らせる
  double xcen(par.bp2xaxis(binsize/2)); // initial position
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


void Page::draw(const DROMPA::Global &p, const int32_t page_curr, const std::string &chrname, const int32_t region_no)
{
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

    /*    if(d->gene.argv || d->arsfile || d->terfile) draw_annotation(d, cr, xstart, xend);
    if(d->internum){
      for(j=0; j<d->internum; ++j) draw_interaction(cr, &(d->inter[j]), xstart, xend, chr);
      }*/
    //    double ytemp = par.yaxis_now;
    int32_t nlayer = 0;
    for(size_t j=0; j<pairs.size(); ++j) stroke_each_layer(p, pairs[j], nlayer);
    stroke_xaxis_num(9);

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
  if(par.num_page>1) title = chrname + "_" + std::to_string(region_no+1) + "_" + std::to_string(page_curr+1);
  else               title = chrname + "_" + std::to_string(region_no+1);
  cr->set_source_rgba(CLR_BLACK, 1);
  showtext_cr(cr, 50, 30, title, 16);
  
  cr->show_page();
  return;
}

void Figure::DrawData(DROMPA::Global &p)
{
  int32_t width(pagewidth);
  int32_t height(p.drawparam.getPageHeight(pairs));

#ifdef CAIRO_HAS_PDF_SURFACE
  const auto surface = Cairo::PdfSurface::create(p.getFigFileNameChr(chrname), width, height);

  std::cout << "Wrote PDF file \"" << p.getFigFileNameChr(chrname) << "\"" << std::endl;

  int32_t start(0);
  int32_t end(chrlen);
  
  int32_t num_page(p.drawparam.getNumPage(start, end));
  int32_t region_no(0);
  for(int32_t i=0; i<num_page; ++i) {
    std::cout << boost::format("   page %5d/%5d/%5d\r") % (i+1) % num_page % (region_no+1) << std::flush;
    Page page(p, arrays, pairs, surface, start, end);
    page.draw(p, i, "chr" + chrname, region_no);
  }
#else
  std::cout << "You must compile cairo with PDF support for DROMPA." << std::endl;
  return;
#endif
}


int32_t DROMPA::DrawParam::getHeightEachSample(const SamplePairChr &pair) const {
  int32_t height(0);
  int32_t n(0);
  if (showctag)                            { height += lineheight; ++n; }
  if (showitag==1 && pair.argvInput != "") { height += lineheight; ++n; }
  if (showratio   && pair.argvInput != "") { height += lineheight; ++n; }
  if (showpinter)                          { height += lineheight; ++n; }
  if (showpenrich && pair.argvInput != "") { height += lineheight; ++n; }
  height += MERGIN_BETWEEN_DATA * (n-1);
  
#ifdef DEBUG
  std::cout << "LineHeight: " << lineheight << ",n: " << n << std::endl;
  std::cout << "HeightEachSample: " << height << std::endl;
#endif
  return height;
}
    
int32_t DROMPA::DrawParam::getHeightAllSample(const std::vector<SamplePairChr> &pairs) const {
  int32_t height(0);
  for(auto x: pairs) height += getHeightEachSample(x);
  height += MERGIN_BETWEEN_DATA * (samplenum-1);
  if(showitag==2) height += lineheight + MERGIN_BETWEEN_DATA;

#ifdef DEBUG
  std::cout << "HeightAllSample; " << height << std::endl;
#endif
  return height;
}
    

int32_t DROMPA::DrawParam::getPageHeight(const std::vector<SamplePairChr> &pairs) const {
  int32_t height(OFFSET_Y*2);
  height += getHeightAllSample(pairs) * linenum_per_page;
  height += MERGIN_BETWEEN_LINE * (linenum_per_page-1);

#ifdef DEBUG
  std::cout << "PageHeight: " << height << std::endl;
#endif
  return height;
}

