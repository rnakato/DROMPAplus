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
  enum {BOXHEIGHT_GENEBOX_EXON=130, BOXHEIGHT_GENEBOX_NOEXON=80};
  enum {BOXHEIGHT_GRAPH=80, MEMNUM_GC=10};
  int32_t pagewidth(1088);
  int32_t width_draw(820);
  double dot_per_bp(0);
  int32_t mergin_between_graph_data(15);
}

inline double bp2xaxis(const int32_t bp) { return bp * dot_per_bp + OFFSET_X; }

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

  DParam(const int32_t s, const int32_t e, const DROMPA::Global &p):
    start(s), end(e),
    num_line(p.drawparam.getNumLine(start, end)),
    num_page(p.drawparam.getNumPage(start, end)),
    width_per_line(p.drawparam.width_per_line),
    yaxis_now(0), xstart(0), xend(0), ystep(12), barnum(2)
  {
    dot_per_bp = getratio(width_draw, width_per_line);
  }

  void set_xstart_xend(const int32_t i) {
    xstart = start + i * width_per_line;
    if(i==num_line-1) xend = end;
    else xend = start + (i+1) * width_per_line -1;
  }
  
  double getXaxisLen() const { return (xend - xstart) * dot_per_bp; }
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
      double val(stof(v[1]));
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

class Page {
  const std::unordered_map<std::string, ChrArray> &arrays;
  const std::vector<SamplePairChr> &pairs;
  GraphData GC;
  GraphData GD;

  DParam par;
  Cairo::RefPtr<Cairo::Context> cr;

  std::string chrname;
  
  public:
  
  Page(const DROMPA::Global &p,
       const std::unordered_map<std::string, ChrArray> &refarrays,
       const std::vector<SamplePairChr> &refpairs,
       const Cairo::RefPtr<Cairo::PdfSurface> surface,
       const chrsize &chr, const int32_t s, const int32_t e):
    arrays(refarrays), pairs(refpairs),
    par(s, e, p), cr(Cairo::Context::create(surface)),
    chrname(chr.getrefname())
  {
    cr->select_font_face( "Arial", Cairo::FONT_SLANT_NORMAL, Cairo::FONT_WEIGHT_NORMAL);
    if(p.anno.GC.isOn()) GC.setValue(p.anno.GC, chrname, chr.getlen(), "GC%",         20, 70, MEMNUM_GC, BOXHEIGHT_GRAPH);
    if(p.anno.GD.isOn()) GD.setValue(p.anno.GD, chrname, chr.getlen(), "Num of genes", 0, 40, MEMNUM_GC, BOXHEIGHT_GRAPH);
  }

  void stroke_each_layer(const DROMPA::Global &p, const SamplePairChr &pair, const int32_t nlayer);

  template <class T>
  void stroke_readdist(const DROMPA::Global &p, const SamplePairChr &pair, const int32_t nlayer)
  {
    par.yaxis_now += getHeightDf() + MERGIN_BETWEEN_DATA;
    T df(cr, p, pair, par, getWidthDf(), getHeightDf());
    df.stroke_bindata(p, pair, arrays, nlayer);
    df.stroke_dataframe(p, nlayer);
    if(!nlayer) stroke_xaxis();
    return;
  }
  void DrawAnnotation(const DROMPA::Global &p);
  void strokeARS(const HashOfGeneDataMap &mp, const double ycenter);
  void strokeGeneSGD(const DROMPA::Global &p, const double ycenter);
  void strokeGene(const DROMPA::Global &p, const double ycenter);
  void drawGraph(const GraphData &graph);

  void Draw(const DROMPA::Global &p, const int32_t page_curr, const int32_t region_no);

  void set_xstart_xend(const int32_t i) {
    par.set_xstart_xend(i);
  }

  /*  
  int32_t set_interval_large() {
    int32_t interval;
    if(d->gw==true){
      if(d->large_genome==false) interval = 100*NUM_1K;  // 100kbp
      else interval = 10*NUM_1M;  // 10Mbp
    }else{
      interval = width_per_line/10;
    }
    return interval;
    }*/
  
  void stroke_xaxis(const double y) {
    double x;
    int32_t interval_large(par.width_per_line/10); //set_interval_large(d);
    int32_t interval = interval_large/10;
    
    cr->set_source_rgba(CLR_BLACK, 1);
    for(int32_t i=setline(par.xstart, interval); i<=par.xend; i+=interval) {
      x = bp2xaxis(i - par.xstart);
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
    int32_t interval(par.width_per_line/10);
    
    cr->set_source_rgba(CLR_BLACK, 1);
    for(int32_t i=setline(par.xstart, interval); i<=par.xend; i+=interval) {
      std::string str;
      x = bp2xaxis(i - par.xstart);
      /*  if(d->gw==true){
	if(d->large_genome==false) sprintf(str, "%dk", i/NUM_1K);
	else sprintf(str, "%dM", i/NUM_1M);
	showtext_cr(cr, x - 3*strlen(str), y+10, str, fontsize);
	}else{*/
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
	showtext_cr(cr, x - 3*str.length(), y+10, str, fontsize);
	//      }
    }
    return;
  }

  void stroke_xaxis(){
    double x;
    //  int32_t interval_large = set_interval_large(d);
    int32_t interval_large(par.width_per_line/10);
    int32_t interval(interval_large/10);

    cr->set_source_rgba(CLR_BLACK, 1);
    for(int32_t i=setline(par.xstart, interval); i<=par.xend; i+=interval) {
      x = bp2xaxis(i-par.xstart);
    
      if (!(i%interval_large)) {
	cr->set_line_width(1);
	rel_yline(cr, x, par.yaxis_now-4, 8);
      } else {
	cr->set_line_width(0.5);
	rel_yline(cr, x, par.yaxis_now-1.5, 3);
      }
      cr->stroke();
    }
    return;
  }

  std::tuple<int32_t, int32_t> get_start_end_linenum(const int32_t page, const int32_t linenum_per_page) const {
    int32_t start(0), end(0);
    start = page * linenum_per_page;
    if(page == par.num_page-1) end = par.num_line;
    else end = (page+1) * linenum_per_page;
    return std::forward_as_tuple(start, end);
  }

  double getWidthDf() const { return (par.xend - par.xstart +1) * dot_per_bp; }
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
  void stroke_dataframe(const DROMPA::Global &p, const int32_t nlayer);
  void stroke_bindata(const DROMPA::Global &p, const SamplePairChr &pair,
		      const std::unordered_map<std::string, ChrArray> &arrays, const int32_t nlayer);

  int32_t getbinlen(const double value) const { return -std::min(par.ystep*value, height_df); }
  
  virtual void stroke_bin(const SamplePairChr &pair,
			  const std::unordered_map<std::string, ChrArray> &arrays,
			  const int32_t i, const double xcen, const int32_t yaxis, const int32_t viz)=0;
};


class ChIPDataFrame : public DataFrame {
 public:
  ChIPDataFrame(const Cairo::RefPtr<Cairo::Context> cr_, const DROMPA::Global &p, const SamplePairChr &pair,
		const DParam &refparam, const double wdf, const double hdf):
    DataFrame(cr_, pair.label, p.drawparam.scale_tag, refparam, wdf, hdf)
  {}

  void stroke_bin(const SamplePairChr &pair,
		  const std::unordered_map<std::string, ChrArray> &arrays,
		  const int32_t i, const double xcen, const int32_t yaxis, const int32_t viz);
};


class InputDataFrame : public DataFrame {
 public:
  InputDataFrame(const Cairo::RefPtr<Cairo::Context> cr_, const DROMPA::Global &p, const SamplePairChr &pair,
		const DParam &refparam, const double wdf, const double hdf):
    DataFrame(cr_, "Input", p.drawparam.scale_tag, refparam, wdf, hdf)
  {}

  void stroke_bin(const SamplePairChr &pair,
		  const std::unordered_map<std::string, ChrArray> &arrays,
		  const int32_t i, const double xcen, const int32_t yaxis, const int32_t viz);
};

class RatioDataFrame : public DataFrame {
 public:
  RatioDataFrame(const Cairo::RefPtr<Cairo::Context> cr_, const DROMPA::Global &p, const SamplePairChr &pair,
		const DParam &refparam, const double wdf, const double hdf):
    DataFrame(cr_, getlabel(p, pair), p.drawparam.scale_ratio, refparam, wdf, hdf)
  {}

  const std::string getlabel(const DROMPA::Global &p, const SamplePairChr &pair) const {
    if(p.drawparam.showctag) return "IP/Input";
    else return pair.label;
  }

  void stroke_bin(const SamplePairChr &pair,
		  const std::unordered_map<std::string, ChrArray> &arrays,
		  const int32_t i, const double xcen, const int32_t yaxis, const int32_t viz);
};

class LogRatioDataFrame : public DataFrame { // log10(ratio)
  int32_t barnum_minus;
  int32_t barnum_plus;
 public:
  LogRatioDataFrame(const Cairo::RefPtr<Cairo::Context> cr_, const DROMPA::Global &p, const SamplePairChr &pair,
		const DParam &refparam, const double wdf, const double hdf):
    DataFrame(cr_, getlabel(p, pair), p.drawparam.scale_ratio, refparam, wdf, hdf),
    barnum_minus(refparam.barnum/2), barnum_plus(refparam.barnum - barnum_minus)
  {}

  const std::string getlabel(const DROMPA::Global &p, const SamplePairChr &pair) const {
    if(p.drawparam.showctag) return "IP/Input";
    else return pair.label;
  }

  void stroke_ymem(const int32_t nlayer)
  {
    cr->set_source_rgba(CLR_BLACK, 1);

    int32_t barnum_minus = par.barnum/2;
    double x(0);
    if (!nlayer) x = OFFSET_X + width_df + 7; else x = OFFSET_X - 20;
    for(int32_t i=1; i<=par.barnum; ++i) {
      std::string str;
      if (i < barnum_minus) str = "1/" + std::to_string(static_cast<int>(pow(2, (barnum_minus-i) * scale)));
      else str = std::to_string(static_cast<int>(pow(2, (i-barnum_minus) * scale)));
      
      showtext_cr(cr, x, par.yaxis_now - i*(par.ystep - 1.5), str, 9);
    }
    return;
  }
  void stroke_bin(const SamplePairChr &pair,
		  const std::unordered_map<std::string, ChrArray> &arrays,
		  const int32_t i, const double xcen, const int32_t yaxis, const int32_t viz);
};

class PinterDataFrame : public DataFrame {
 public:
  PinterDataFrame(const Cairo::RefPtr<Cairo::Context> cr_, const DROMPA::Global &p, const SamplePairChr &pair,
		const DParam &refparam, const double wdf, const double hdf):
    DataFrame(cr_, getlabel(p, pair), p.drawparam.scale_ratio, refparam, wdf, hdf)
  {}

  const std::string getlabel(const DROMPA::Global &p, const SamplePairChr &pair) const {
    if(p.drawparam.showctag || p.drawparam.showratio) return "pval (ChIP internal)";
    else return pair.label;
  }

  void stroke_bin(const SamplePairChr &pair,
		  const std::unordered_map<std::string, ChrArray> &arrays,
		  const int32_t i, const double xcen, const int32_t yaxis, const int32_t viz);
};

class PenrichDataFrame : public DataFrame {
 public:
  PenrichDataFrame(const Cairo::RefPtr<Cairo::Context> cr_, const DROMPA::Global &p, const SamplePairChr &pair,
		   const DParam &refparam, const double wdf, const double hdf):
    DataFrame(cr_, getlabel(p, pair), p.drawparam.scale_ratio, refparam, wdf, hdf)
  {}

  const std::string getlabel(const DROMPA::Global &p, const SamplePairChr &pair) const {
    if(p.drawparam.showctag || p.drawparam.showratio) return "pval (IP/Input)";
    else return pair.label;
  }

  void stroke_bin(const SamplePairChr &pair,
		  const std::unordered_map<std::string, ChrArray> &arrays,
		  const int32_t i, const double xcen, const int32_t yaxis, const int32_t viz);
};


class GeneElement{
  int32_t dif;
  int32_t cnt;
  
 public:
  double x1, x2, xcen, xwid;
  double x_name;
  int32_t y_name, ylen;
  int32_t ybar;

  GeneElement(const std::pair<const std::string, genedata> &m,
	      const int32_t xstart, const double ycenter,
	      const int32_t ty,
	      int32_t &on_plus, int32_t &on_minus):
    dif(6), cnt(1),
    x1(bp2xaxis(m.second.txStart - xstart +1)),
    x2(bp2xaxis(m.second.txEnd   - xstart +1)),
    xcen((x1+x2)/2), xwid(x2 - x1),
    ybar(0)
  {
    x_name = xcen - 3.25 * m.second.gname.length() + 6;
    if (x_name < 0) x_name = 10;
    else if (x_name > pagewidth) x_name = pagewidth - 50;

    if(!ty) { // SGD
      if (m.second.strand == "+") {
	ybar   = ycenter - dif;
	y_name = ybar -5 - on_minus*6;
	if(on_minus == cnt) on_minus=0; else ++on_minus;
      }
      else if (m.second.strand == "-") {
	ybar   = ycenter + dif;
	y_name = ybar +9 + on_plus*6;
	if (on_plus == cnt) on_plus=0; else ++on_plus;
      }
      else y_name = ycenter -22;
      ylen = y_name - ycenter;
    } else {  // Others
      if (m.second.strand == "+") {
	ybar = ycenter -6 - on_minus * 13;
	y_name = ybar -6;
	if(on_minus==3) on_minus=0; else ++on_minus;
      } else{
	ybar = ycenter +12 + on_plus * 13;
	y_name = ybar +11;
	if(on_plus==3) on_plus=0; else ++on_plus;
      }
    }
    
  }
};

#endif /* _DD_READFILE_P_H_ */
