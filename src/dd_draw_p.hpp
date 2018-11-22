/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_DRAW_P_H_
#define _DD_DRAW_P_H_

#include <cairommconfig.h>
#include <cairomm/context.h>
#include <cairomm/surface.h>
#include "dd_peakcall.hpp"

#define rel_xline(cr, x1, y1, xlen) do{		\
    cr->move_to(x1,   (int32_t)y1);		\
    cr->line_to(x1+xlen, (int32_t)y1); }while(0)
#define rel_yline(cr, x1, y1, ylen) do{		\
    cr->move_to(x1, (int32_t)y1);		\
    cr->line_to(x1, (int32_t)(y1+ylen)); }while(0)

using ChrArrayMap = std::unordered_map<std::string, ChrArray>;

inline double CalcRatio(const double c, const double i, const double r)
{
  return i ? c/i*r: 0;
}

namespace {
  enum class LineType{CHIP, INPUT, RATIO, RATIO_GV, PVALUE_INTER, PVALUE_ENRICH};

  enum {OFFSET_X=190, OFFSET_Y=50, MERGIN_BETWEEN_DATA=10, MERGIN_BETWEEN_LINE=30};
  enum {BOXHEIGHT_GENEBOX_EXON=140, BOXHEIGHT_GENEBOX_NOEXON=60};
  enum {BOXHEIGHT_GRAPH=80, MEMNUM_GC=10, MERGIN_BETWEEN_GRAPH_DATA=15, BOXHEIGHT_INTERACTION=45, BOXHEIGHT_BEDANNOTATION=12};
  int32_t pagewidth(1088);
  int32_t width_draw(750);
}

inline int32_t setline(const int32_t start, const int32_t interval)
{
  int32_t posi(start-1);
  if(!posi%interval) return posi;
  else return (posi/interval +1) * interval;
}

inline void showtext_cr(const Cairo::RefPtr<Cairo::Context> cr, const double x, const double y, const std::string &str, const int32_t fontsize)
{
  cr->move_to(x, y);
  cr->set_font_size(fontsize);
  cr->show_text(str);
  cr->stroke();
  return;
}

inline const std::string float2string(const double f, const int32_t digits)
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

  double yaxis_now;
  int32_t xstart;
  int32_t xend;

  double ystep;
  int32_t barnum;
  
  double dot_per_bp;

  double alpha;
  
  DParam(const int32_t s, const int32_t e, const DROMPA::Global &p):
    start(s), end(e),
    num_line(p.drawparam.getNumLine(start, end)),
    num_page(p.drawparam.getNumPage(start, end)),
    width_per_line(p.drawparam.width_per_line),
    yaxis_now(0), xstart(0), xend(0), ystep(12), barnum(2),
    dot_per_bp(getratio(width_draw, width_per_line)),
    alpha(p.drawparam.alpha)
  {}

  void set_xstart_xend(const int32_t i) {
    xstart = start + i * width_per_line;
    if(i==num_line-1) xend = end;
    else xend = start + (i+1) * width_per_line -1;
  }
  
  double getXaxisLen() const { return (xend - xstart) * dot_per_bp; }
  
  double bp2xaxis(const int32_t bp) const { return bp * dot_per_bp + OFFSET_X; }
  
  double getHeightDf() const { return ystep * barnum; }
};


class GraphData {
  enum {MEM_MAX_DEFAULT=20};
public:
  int32_t binsize;
  std::vector<double> array;
  std::string label;
  int32_t memnum;
  int32_t boxheight;
  
  double mmin;
  double mmax;
  double mwid;
  GraphData(){}

  void setValue(const GraphFile &g, const std::string &chr, const int32_t chrlen, const std::string &l,
		const double ymin, const double ymax,
		const int32_t m, const int32_t bh)
  {
    binsize = g.getbinsize();
    label = l;
    memnum = m;
    boxheight = bh;
    std::string filename(g.getfilename() + "/" + chr + "-bs" + std::to_string(binsize));
    setArray(filename, chrlen, ymin, ymax);
  }
  void setArray(const std::string &filename, const int32_t chrlen, const double ymin, const double ymax) {
    std::ifstream in(filename);
    if (!in) PRINTERR("cannot open " << filename);

    array = std::vector<double>(chrlen/binsize +1);
    
    double maxtemp(0);
    std::string lineStr;
    while (!in.eof()) {
      getline(in, lineStr);
      if (lineStr.empty()) continue;
      
      std::vector<std::string> v;
      boost::split(v, lineStr, boost::algorithm::is_any_of(" \t"), boost::algorithm::token_compress_on);

      int32_t start(stoi(v[0]));
      if(start % binsize){
	printf("%d %d\n", start, binsize);
	PRINTERR("[E]graph: invalid start position or binsize:  " << filename);
      }
      double val(stod(v[1]));
      array[start/binsize] = val;
      if(maxtemp < val) maxtemp = val;
    }
    if (ymax) {
      mmin = ymin;
      mmax = ymax;
    } else {
      mmin = 0;
      if(MEM_MAX_DEFAULT > maxtemp) mmax = MEM_MAX_DEFAULT;
      else mmax = maxtemp;
    }
    mwid = mmax - mmin;
  }
  double getylen(const int32_t i) const {
    return boxheight * (array[i] - mmin)/mwid;
  }
  double getBoxHeight4mem() const {
    return boxheight/memnum;
  }
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

  const ChrArrayMap &arrays;
  const std::vector<SamplePairOverlayed> &vsamplepairoverlayed;
  GraphData GC;
  GraphData GD;

  Cairo::RefPtr<Cairo::Context> cr;

  std::string chrname;
  
  int32_t setInterval() const {
    int32_t interval;
    if (par.width_per_line > 100*NUM_1M) interval = 10*NUM_1M;       // 10Mbp
    else if (par.width_per_line > 10*NUM_1M) interval = 100*NUM_1K;  // 100kbp
    else interval = par.width_per_line/10;
    return interval;
  }
  
  void StrokeEachLayer(const DROMPA::Global &p);
  void StrokeReadLines(const DROMPA::Global &p, const SamplePairOverlayed &pair);
  void StrokeGraph(const GraphData &graph);
  void DrawGeneAnnotation(const DROMPA::Global &p);
  void strokeARS(const HashOfGeneDataMap &mp, const double ycenter);
  void strokeGeneSGD(const DROMPA::Global &p, const double ycenter);
  void strokeGene(const DROMPA::Global &p, const double ycenter);

  void drawBedAnnotation(const vbed<bed12> &vbed);
  void drawInteraction(const InteractionSet &vinter);

  public:
  DParam par;
  
  PDFPage(const DROMPA::Global &p,
       const ChrArrayMap &refarrays,
       const std::vector<SamplePairOverlayed> &pair,
       const Cairo::RefPtr<Cairo::PdfSurface> surface,
       const chrsize &chr, const int32_t s, const int32_t e):
    arrays(refarrays),
    vsamplepairoverlayed(pair),
    cr(Cairo::Context::create(surface)),
    chrname(chr.getrefname()),
    par(s, e, p)
  {
    cr->select_font_face( "Arial", Cairo::FONT_SLANT_NORMAL, Cairo::FONT_WEIGHT_NORMAL);
    if(p.anno.GC.isOn()) GC.setValue(p.anno.GC, chrname, chr.getlen(), "GC%",         20, 70, MEMNUM_GC, BOXHEIGHT_GRAPH);
    if(p.anno.GD.isOn()) GD.setValue(p.anno.GD, chrname, chr.getlen(), "Num of genes", 0, 40, MEMNUM_GC, BOXHEIGHT_GRAPH);
  }

  template <class T>
  void StrokeDataFrame(const DROMPA::Global &p, const SamplePairOverlayed &pair)
  {
    par.yaxis_now += par.getHeightDf() + MERGIN_BETWEEN_DATA;
    T df(cr, p, pair, par, chrname);

    df.Draw(p, pair, arrays);
    stroke_xaxis(par.yaxis_now);

    return;
  }

  void MakePage(const DROMPA::Global &p, const int32_t page_no, const std::string &pagelabel);

  void set_xstart_xend(const int32_t i) { par.set_xstart_xend(i); }

  void stroke_xaxis(const double y) {
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

  void stroke_xaxis_num(const double y, const int32_t fontsize) {
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

  void StrokeWidthOfInteractionSite(const bed site, const double y) {
    cr->set_line_width(2);
    cr->set_source_rgba(CLR_DARKORANGE, 0.8);
    double s = par.bp2xaxis(site.start - par.xstart);
    double e = par.bp2xaxis(site.end - par.xstart);
    rel_xline(cr, s, y, e-s);
    cr->stroke();
  }
  
  // cr->arc(中心x, 中心y, 半径, start角度, end角度) 角度はラジアン
  void drawArc_from_to(const Interaction &inter, const int32_t start, const int32_t end, const int32_t ref_height, const double ref_ytop) {
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
  
  void drawArc_from_none(const Interaction &inter, const int32_t start, const int32_t end, const int32_t ref_height, const double ref_ytop) {
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
  void drawArc_none_to(const Interaction &inter, const int32_t start, const int32_t end, const int32_t ref_height, const double ref_ytop) {
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

  std::tuple<int32_t, int32_t> get_start_end_linenum(const int32_t page, const int32_t linenum_per_page) const {
    int32_t start(0), end(0);
    start = page * linenum_per_page;
    if(page == par.num_page-1) end = par.num_line;
    else end = (page+1) * linenum_per_page;
    return std::forward_as_tuple(start, end);
  }
};


#endif /* _DD_READFILE_P_H_ */
