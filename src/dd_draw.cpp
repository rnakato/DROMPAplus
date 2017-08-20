/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */

#include <cairommconfig.h>
#include <cairomm/context.h>
#include <cairomm/surface.h>
#include <sstream>
#include <iomanip>
#include "dd_draw.hpp"
#include "dd_draw_p.hpp"
#include "color.hpp"
#include "SSP/common/inline.hpp"

#define rel_xline(cr, x1, y1, xlen) do{		\
    cr->move_to(x1,   (int32_t)y1);		\
    cr->line_to(x1+xlen, (int32_t)y1); }while(0)
#define rel_yline(cr, x1, y1, ylen) do{		\
    cr->move_to(x1, (int32_t)y1);		\
    cr->line_to(x1, (int32_t)(y1+ylen)); }while(0)

void showtext_cr(Cairo::RefPtr<Cairo::Context> cr, const double x, const double y, const std::string &str, const int32_t fontsize)
{
  cr->move_to(x, y);
  cr->select_font_face("Serif", Cairo::FONT_SLANT_NORMAL, Cairo::FONT_WEIGHT_NORMAL);
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

/*
double define_value_and_color(const DROMPA::Global &p, const SamplePairChr &pair, const int32_t i, const LineType type, const int32_t nlayer)
{
  gdouble value(0), ratio, data;
  switch(type) {
  case LineType::CHIP:
    value = arrays[pair.argvChIP][i] / p.scale.scale_tag;
    if(!nlayer) cairo_set_source_rgba(cr, CLR_GREEN3, r_trans);
    else cairo_set_source_rgba(cr, CLR_ORANGE, r_trans);
    if(sample->peakarray && sample->peakarray[i] && d->do_peakcall) cairo_set_source_rgba(cr, CLR_PINK, r_trans);
    if(!sample->peakarray && d->do_peakcall) peakcall(p, cr, sample, i);
    break;
  case LineType::INPUT:
    value = WIGARRAY2VALUE(sample->Input->data[i]) /sample->scale_tag;
    cairo_set_source_rgba(cr, CLR_BLUE, r_trans); 
    break;
  case LineType::RATIO_GV: // same as LTYPE_RATIO 
  case LineType::RATIO:
    ratio = CALCRATIO(sample->ChIP->data[i], sample->Input->data[i], sample->comp->genome->ratio);
    if(!ratio) value = 0;
    else{
      if(d->visualize_ratio==2) value = log10(ratio) / sample->scale_ratio;
      else value = ratio / sample->scale_ratio; 
    }
    if(type==LineType::RATIO_GV){
      if(d->visualize_ratio==1 && value > 1) cairo_set_source_rgba(cr, CLR_ORANGE, r_trans);
      else if(d->visualize_ratio==2 && value > 0) cairo_set_source_rgba(cr, CLR_ORANGE, r_trans);
      else cairo_set_source_rgba(cr, CLR_SLATEGRAY, r_trans);
      break;
    }
    if(d->do_peakcall){
      if(ratio > p->enrichthre) cairo_set_source_rgba(cr, CLR_PINK, r_trans);
      else cairo_set_source_rgba(cr, CLR_GRAY, r_trans);
    }else cairo_set_source_rgba(cr, CLR_ORANGE, r_trans);
    if(nlayer) cairo_set_source_rgba(cr, CLR_BLUE, r_trans);
    if(d->visualize_ratio==2 && value < 0) cairo_set_source_rgba(cr, CLR_GRAY3, r_trans);
    break;
  case LineType::PVALUE_INTER:
    data = zero_inflated_binomial_test(WIGARRAY2VALUE(sample->ChIP->data[i]), sample->ChIP->nb_p, sample->ChIP->nb_n);
    value = data /sample->scale_pvalue;
    if(data > p->pthre_internal) cairo_set_source_rgba(cr, CLR_PINK, r_trans);
    else cairo_set_source_rgba(cr, CLR_GRAY, r_trans);
    break;
  case LineType::PVALUE_ENRICH:
    data = binomial_test(WIGARRAY2VALUE(sample->ChIP->data[i]), WIGARRAY2VALUE(sample->Input->data[i]), sample->comp->genome->ratio);
    value = data /sample->scale_pvalue;
    if(data > p->pthre_enrich) cairo_set_source_rgba(cr, CLR_PINK, r_trans);
    else cairo_set_source_rgba(cr, CLR_GRAY, r_trans);
    break;
  default:
    fprintf(stderr, "error: invalid stroke value.\n");
    exit(1);
  }
  return value;
  }*/

  
void Page::stroke_bindata(const DROMPA::Global &p, const SamplePairChr &pair, const LineType type, const int32_t nlayer)
{
  int32_t binsize(pair.getbinsize());
  int32_t sbin(xstart/binsize);
  int32_t ebin(xend/binsize);
  double dot_per_bin(binsize * dot_per_bp);
  int32_t yaxis(yaxis_now);
  int32_t dfHeight(getHeightDf());

  double value(0);
  int32_t thin(std::min(width_per_line/(1000*binsize), 20));
  int32_t barnum_minus(barnum/2);
  int32_t barnum_plus(barnum - barnum_minus);

  // peak region 色を変えるのではなく、その領域のバックを光らせる
  double xcen(bp2axis(binsize/2));  // initial position
  if(thin > 1) cr->set_line_width(dot_per_bin*thin);
  else cr->set_line_width(dot_per_bin);

  for (int32_t i=sbin; i<ebin; ++i, xcen += dot_per_bin) {
    if (thin > 1 && i%thin) continue;

    /*    if (p->gapfile && d->gaparray[i] >= GAP_THRE){
      cr->set_source_rgba(CLR_BLACK, 0.3);
      rel_yline(cr, xcen, yaxis - dfHeight, dfHeight);
      cr->stroke();
    } else if(p->mpfile && d->mparray[i] < p->mpthre){
      cr->set_source_rgba(CLR_BLACK, 0.3);
      rel_yline(cr, xcen, yaxis - dfHeight, dfHeight);
      cr->stroke();
      }*/

    value = (arrays.at(pair.argvChIP)).array[i] / p.scale.scale_tag;
    cr->set_source_rgba(CLR_ORANGE, 1);
    
    //    value = define_value_and_color(p, pair, i, type, nlayer);
    if (!value) continue;

    /*if (d->visualize_ratio==2) {
      if(value>0) len = -d->ystep * (min(value, barnum_plus));
      else        len = -d->ystep * (max(value, -barnum_minus));
      rel_yline(cr, xcen, (gint)yaxis_now - d->ystep * barnum_minus, (gint)len);
      } else {*/
    double len(std::max(-ystep*value, -getHeightDf()));

    //      if (!d->viz) {
    //	rel_yline(cr, xcen, (gint)yaxis_now, (gint)len);
    //  } else { // d->viz ==1
    if(len <-1 && value <= barnum) rel_yline(cr, xcen, (int)(yaxis_now + len +1), 1);
    cr->stroke();
    //    define_value_and_color(p, d, cr, sample, i, type, 0.5, nlayer);
    rel_yline(cr, xcen, (int)yaxis_now, (int)len);
    // }
      //    }
    cr->stroke();
  }

  return;
}


void Page::stroke_dataframe(const DROMPA::Global &p, const SamplePairChr &pair, const LineType type, const int32_t nlayer)
{
  double width_df(getWidthDf());
  double height_df(getHeightDf()); 
  double scale=0;

  /* y memory */
  cr->set_line_width(0.4);
  cr->set_source_rgba(CLR_BLACK, 0.5);
  for(int32_t i=0; i<barnum; i++) rel_xline(cr, OFFSET_X, yaxis_now - i*ystep, width_df);
  cr->stroke();

  /* y key */
  if(p.drawparam.isshowymem()) {
    int32_t barnum_minus = barnum/2;
    if(type == LineType::CHIP               || type == LineType::INPUT)        scale = p.scale.scale_tag;
    else if(type == LineType::RATIO         || type == LineType::RATIO_GV)     scale = p.scale.scale_ratio;
    else if(type == LineType::PVALUE_ENRICH || type == LineType::PVALUE_INTER) scale = p.scale.scale_pvalue;
    else{ PRINTERR("error: invalid LineType."); }

    cr->set_source_rgba(CLR_BLACK, 1);
    double x(0);
    if(!nlayer) x = OFFSET_X + width_df + 7; else x = OFFSET_X - 20;
    for(int32_t i=1; i<=barnum; i++){
      double y = yaxis_now - i*(ystep - 1.5);
      std::string str;
      if(p.drawparam.showratio==2){
	if(i < barnum_minus) str = "1/" + std::to_string(static_cast<int>(pow(2, (barnum_minus-i) * scale)));
	else str = std::to_string(static_cast<int>(pow(2, (i-barnum_minus) * scale)));
      }else str = float2string(i*scale, 1);
      showtext_cr(cr, x, y, str, 9);
    }
  }
  if(nlayer) return;

  /* frame */
  cr->set_line_width(0.4);
  cr->set_source_rgba(CLR_BLACK, 1);
  rel_xline(cr, OFFSET_X, yaxis_now, width_df);
  rel_yline(cr, OFFSET_X, yaxis_now - height_df, height_df);
  cr->stroke();
  
  /* keys */
  if(p.drawparam.isshowylab()){
  cr->set_source_rgba(CLR_BLACK, 1);
  std::string str;
    switch(type){
    case LineType::CHIP: str = pair.label; break;
    case LineType::INPUT: str = "Input"; break;
    case LineType::RATIO_GV: // same as LineType::RATIO
    case LineType::RATIO:
      if(p.drawparam.showctag) str = "IP/Input";
      else str = pair.label;
      break;
    case LineType::PVALUE_INTER:
      if(p.drawparam.showctag || p.drawparam.showratio) str = "pval (ChIP internal)";
      else str = pair.label;
      break;
    case LineType::PVALUE_ENRICH: 
      if(p.drawparam.showctag || p.drawparam.showratio) str = "pval (IP/Input)";
      else str = pair.label;
      break;
    }
    showtext_cr(cr, 50, yaxis_now - ystep * barnum/2, str, 12);
  }

  return;
}

void Page::stroke_readdist(const DROMPA::Global &p, const SamplePairChr &pair, const LineType type, const int32_t nlayer)
{
  yaxis_now += getHeightDf() + MERGIN_BETWEEN_DATA;
  stroke_bindata(p, pair, type, nlayer);
  stroke_dataframe(p, pair, type, nlayer);
  // if(!nlayer) stroke_xaxis(d, cr, xstart, xend);*/
  return;
}


void Page::stroke_each_layer(const DROMPA::Global &p, const SamplePairChr &pair, const int32_t nlayer)
{
  //  if(d->visualize_p_inter)  stroke_readdist(p, d, cr, sample, xstart, xend, LTYPE_PVALUE_INTER,  nlayer);
  //if(d->visualize_p_enrich) stroke_readdist(p, d, cr, sample, xstart, xend, LTYPE_PVALUE_ENRICH, nlayer);
  /* if(d->visualize_ratio){
    if(p->ftype==FTYPE_GV)  stroke_readdist(p, d, cr, sample, xstart, xend, LTYPE_RATIO_GV, nlayer);
    else                    stroke_readdist(p, d, cr, sample, xstart, xend, LTYPE_RATIO,    nlayer);
    }*/
  if(p.drawparam.showctag) stroke_readdist(p, pair, LineType::CHIP, nlayer);
  //if(d->visualize_itag==1)  stroke_readdist(p, d, cr, sample, xstart, xend, LTYPE_INPUT,    nlayer);
  return;
}

void Page::draw(const DROMPA::Global &p, const std::vector<SamplePairChr> &pairs, const int32_t page_curr)
{
  int32_t line_start, line_end;
  std::tie(line_start, line_end) = get_start_end_linenum(page_curr, p.drawparam.getlpp());
 
  yaxis_now = OFFSET_Y;
  cr->set_source_rgba(CLR_WHITE, 1);
  cr->paint();

  for(int i=line_start; i<line_end; i++){
    set_xstart_xend(i);
    if(xstart >= xend) continue;
    
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
      for(j=0; j<d->internum; j++) draw_interaction(cr, &(d->inter[j]), xstart, xend, chr);
      }*/
    double ytemp = yaxis_now;
    int32_t nlayer = 0;
    for(size_t j=0; j<pairs.size(); j++) stroke_each_layer(p, pairs[j], nlayer);
    //stroke_xaxis_num(d, cr, xstart, xend, yaxis_now, 9);

    //    if(d->bednum) draw_bedfile(d, cr, xstart, xend, chr);
    //if(d->repeat.argv) draw_repeat(d, cr, xstart, xend);
      
    /*    if(p->samplenum_overlay){
      nlayer = 1;
      yaxis_now = ytemp;
      for(; j<p->samplenum; j++) stroke_each_layer(p, d, &(sample[j]), cr, xstart, xend, nlayer);
      }*/
    //    if(d->visualize_itag==2) stroke_readdist(p, d, cr, &(sample[0]), xstart, xend, LTYPE_INPUT, 0);
      
    yaxis_now += MERGIN_BETWEEN_LINE;
  }
  //  print_pagetitle(cr, g->chr[chr].name, page_curr, num_page, region_no);
  cr->show_page();
  return;
}

void Figure::DrawData(DROMPA::Global &p)
{
  int32_t width(pagewidth);
  int32_t height(p.drawparam.getPageHeight());

#ifdef CAIRO_HAS_PDF_SURFACE
  const auto surface = Cairo::PdfSurface::create(p.getFigFileNameChr(chrname), width, height);

  std::cout << "Wrote PDF file \"" << p.getFigFileNameChr(chrname) << "\"" << std::endl;

  int32_t start(0);
  int32_t end(chrlen);
  
  int32_t num_page(p.drawparam.getNumPage(start, end));
  int32_t region_no(0);
  for(int32_t i=0; i<num_page; ++i) {
    std::cout << boost::format("   page %5d/%5d/%5d\r") % (i+1) % num_page % (region_no+1) << std::flush;
    Page page(p, arrays, surface, start, end);
    page.draw(p, pairs, i);
  }
#else
  std::cout << "You must compile cairo with PDF support for DROMPA." << std::endl;
  return;
#endif
}


int32_t DROMPA::DrawParam::getHeightEachSample() const {
  int32_t height(0);
  int32_t n(0);
  if (showpinter)  { height += lineheight; ++n; }
  if (showpenrich) { height += lineheight; ++n; }
  if (showratio)   { height += lineheight; ++n; }
  if (showctag)    { height += lineheight; ++n; }
  if (showitag==1) { height += lineheight; ++n; }
  height += MERGIN_BETWEEN_DATA * (n-1);
  
#ifdef DEBUG
  std::cout << "HeightEachSample; " << height << std::endl;
#endif
  return height;
}
    
int32_t DROMPA::DrawParam::getHeightAllSample() const {
  int32_t height(0);
  height += getHeightEachSample() * samplenum;
  height += MERGIN_BETWEEN_DATA * (samplenum-1);
  if(showitag==2) height += lineheight + MERGIN_BETWEEN_DATA;

#ifdef DEBUG
  std::cout << "HeightAllSample; " << height << std::endl;
#endif
  return height;
}
    

int32_t DROMPA::DrawParam::getPageHeight() const {
  int32_t height(OFFSET_Y*2);
  height += getHeightAllSample() * linenum_per_page;
  height += MERGIN_BETWEEN_LINE * (linenum_per_page-1);

#ifdef DEBUG
  std::cout << "PageHeight: " << height << std::endl;
#endif
  return height;
}


//int32_t Figure::calc_pageheight(const DROMPA::Global &p){
//  int32_t height(OFFSET_Y*2), height_lpp(0);

  /* GC and GD */
  //  if(d->GC.argv) height += boxheight_graph + mergin_between_graph_data;
  //if(d->GD.argv) height += boxheight_graph + mergin_between_graph_data;
  //if(p->ftype==FTYPE_PD){  // PD
  //  height += BOXHEIGHT_PD *d->pdnum + MERGIN_BETWEEN_GRAPH_DATA *(d->pdnum-1);
  //}else{
    /* annotation */
  //    if(d->gftype == GFTYPE_REFFLAT || d->gftype == GFTYPE_ENSEMBL) height_lpp += BOXHEIGHT_GENEBOX_EXON;
  //  else height_lpp += BOXHEIGHT_GENEBOX_NOEXON;
  //    height_lpp += MERGIN_BETWEEN_DATA;
  //height_lpp += BOXHEIGHT_INTERACTION * d->internum;
    /* dataline */

  ///  height_lpp = p.drawparam.getHeightAllSample();
    /* bedline */
    /*if(d->bednum){
      height_lpp += MERGIN_FOR_BED;
      height_lpp += LINEHEIGHT_BED * d->bednum;
      }*/
    //    if(d->repeat.argv) height_lpp += MERGIN_FOR_REPEAT + HEIGHT_REPEAT * (REPEATTYPENUM-3);
  // height += height_lpp * p.drawparam.linenum_per_page + mergin_between_line * (p.drawparam.linenum_per_page-1);
    //}

  //  printf("%d %d %d %d %d\n", height, height_lpp, height_dataline, height_sample, lineheight);

// return height;}

    /*   
    cr->save();
    // draw a border around the image
    cr->set_line_width(20.0); // make the line wider
    cr->rectangle(0.0, 0.0, cairo_image_surface_get_width(surface->cobj()), height);
    cr->stroke();

    cr->set_source_rgba(0.0, 0.0, 0.0, 0.7);
    // draw a circle in the center of the image
    cr->arc(width / 2.0, height / 2.0,             height / 4.0, 0.0, 2.0 * M_PI);
    cr->stroke();

    // draw a diagonal line
    cr->move_to(width / 4.0, height / 4.0);
    cr->line_to(width * 3.0 / 4.0, height * 3.0 / 4.0);
    cr->stroke();
    cr->restore();

    cr->show_page();

    */
