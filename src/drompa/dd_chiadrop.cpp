/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */

#include "dd_draw.hpp"
#include "dd_draw_pdfpage.hpp"
#include "color.hpp"
#include "../submodules/SSP/common/inline.hpp"
#include "../submodules/SSP/common/util.hpp"

namespace {
  class posivector {
    int32_t start;
  public:
    std::vector<int32_t> v;
    posivector(std::vector<int32_t> _v, int32_t start): start(start), v(_v) {}
    virtual ~posivector(){}

    bool operator < (const posivector &another) const {
      if (v[0] < start && another.v[0] < start) {
	if (v[v.size()-1] < another.v[another.v.size()-1]) return 1;
	else return 0;
      }
      if (v[0] < start ) return 1;
      if (another.v[0] < start ) return 0;
      if (v[0] < another.v[0]) return 1;
      if (v[0] == another.v[0]) {
	if (v[v.size()-1] < another.v[another.v.size()-1]) return 1;
	else return 0;
      }
/*      printf("v: ");
      for (auto &x: v) printf("%d ", x);
      printf("\n");
      for (auto &x: another.v) printf("%d ", x);
      printf("\n");*/
      return 0;
    };
/*    bool operator < (const posivector &another) const {
      for (size_t i=0; i<std::min(v.size(), another.v.size()); ++i) {
	if ((v[i] < start && another.v[i] < start)
	    || v[i] == another.v[i]) {
	  if (v[v.size()-1] < another.v[another.v.size()-1]) return 1;
	  else return 0;
	}
	if (v[i] < start ) return 1;
	if (another.v[i] < start ) return 0;
	if (v[i] < another.v[i]) return 1;
      }
      std::cout << v.size() << "ttt" << another.v.size() << std::endl;
      return 0;
    };*/
    bool operator == (const posivector &another) {
      if (v[0] == another.v[0] && v[v.size()-1] == another.v[another.v.size()-1]) return 1;
      else return 0;
    }
/*    bool operator == (const posivector &another) {
      if(v.size() != another.v.size()) return 0;

      for (size_t i=0; i<v.size(); ++i) {
	if (v[i] != another.v[i]) return 0;
      }
      return 1;
    }*/
    void printvnum() const {
      printf("vnum: %lu\n", v.size());
    }
  };

  RGB getInterRGB(double val)
  {
    val = val*0.4 + 0.6;
    HSV color(val, 1.0, 1.0);

    return HSVtoRGB(color);
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

void PDFPage::strokeChIADropBarcode(const std::vector<int32_t> &v, const std::string &nbarcode, const double _ywidth, const double yaxis, const RGB &color)
{
  double ywidth = std::min(_ywidth, 0.4);
  double ycenter(yaxis + ywidth/2);

  int32_t s = std::max(v[0], par.xstart);
  int32_t e = std::min(v[v.size()-1], par.xend);

  //  cr->set_source_rgba(CLR_GRAY2, 1);
  cr->set_source_rgba(color.r, color.g, color.b, 0.4);
  cr->set_line_width(ywidth*0.1);
  rel_xline_double(cr, BP2PIXEL(s - par.xstart), ycenter, (e-s) * par.dot_per_bp);
  cr->stroke();

  // barcode number
  cr->set_source_rgba(CLR_BLACK, 1);
  showtext_cr(cr, BP2PIXEL(s - par.xstart) - 3.5, yaxis + ywidth, nbarcode, 1.0);
  cr->stroke();

  // barcode
  cr->set_line_width(ywidth * 2);
  cr->set_source_rgba(color.r, color.g, color.b, 1);
  for (auto &posi: v) {
    if(posi >= par.xstart && posi <= par.xend) {
      double x1 = BP2PIXEL(posi - par.xstart);
//      double len = std::max(1000 * par.dot_per_bp, 0.05);
      double len = std::max(1000 * par.dot_per_bp, 4.0);
      rel_xline_double(cr, x1 - len/2, ycenter, len);
      cr->stroke();
    }
  }
}

void splitbarcode(const std::vector<int32_t> &v,
		  std::vector<std::vector<int32_t>> &vv,
		  int32_t distance_thre)
{
  for (size_t i=1; i<v.size(); ++i) {
    if (v[i] - v[i-1] > distance_thre) {
      std::vector<int> one(v.begin(), v.begin() + i);
      std::vector<int> rest(v.begin() + i, v.end());
      vv.emplace_back(one);
      splitbarcode(rest, vv, distance_thre);
      return;
    }
  }
  vv.emplace_back(v);
  return;
}

void add_barcode(std::vector<posivector> &vChIA, const std::vector<int32_t> &v, int32_t start, int32_t end)
{
  if (v.size() ==1) return;
  if (v[v.size()-1] < start || v[0] > end) return;
  if (v[0] < start && v[v.size()-1] > end) return;

  vChIA.emplace_back(v, start);
}

void PDFPage::StrokeChIADrop(const DROMPA::Global &p)
{
  DEBUGprint_FUNCStart();

  int32_t boxheight(BOXHEIGHT_ChIADROP);
  std::string chr(rmchr(chrname));

  /* frame */
  cr->set_source_rgba(CLR_GRAY4, 1);
  cr->rectangle(OFFSET_X, par.yaxis_now, par.getXaxisLen(), boxheight);
  cr->stroke();

  std::vector<posivector> vChIA;

  // chia dataが0でない場合
  if (p.anno.mp_ChIADrop.find(chr) != p.anno.mp_ChIADrop.end()) {
    for (auto &x: p.anno.mp_ChIADrop.at(chr)) {
      const std::vector<int32_t> &v = x.second;
//      std::vector<int32_t> v{100, 200, 300, 100000, 300000, 400000, 450000, 700000, 10000000};
      if (v.size() ==1) continue;
      if (v[v.size()-1] < par.xstart || v[0] > par.xend) continue;
      if (v[0] < par.xstart && v[v.size()-1] > par.xend) continue;

      std::vector<std::vector<int32_t>> vv;
      splitbarcode(v, vv, p.anno.chia_distance_thre);
#ifdef DEBUG
      for (auto &x: vv) {
	printf("xxx:  ");
	for (auto &y:x) {
	  std::cout << y << " ";
	}
	printf("\n");
      }
#endif
      for (auto &x: vv) {
	add_barcode(vChIA, x, par.xstart, par.xend);
      }
    }

#ifdef DEBUG
    std::cout << "ChIA size: " << vChIA.size() << std::endl;
    for (auto &x: vChIA) x.printvnum();
#endif

    std::sort(vChIA.begin(), vChIA.end());

#ifdef DEBUG
    printf("sortdone\n");
#endif

    int32_t max(0);
    int32_t num_line(0);
    for (size_t i=0; i<vChIA.size(); ++i) {
      int32_t n(1);
      while(i < vChIA.size()-1 && vChIA[i].v == vChIA[i+1].v) { ++i; ++n; }
      max = std::max(max, n);
      ++num_line;
    }

#ifdef DEBUG
    printf("num line %d\n", num_line);
#endif

    showColorBar_ChIADrop(cr, 80, par.yaxis_now + 10, max);

    double ywidth = std::min(boxheight/(double)num_line, 2.0);

    int32_t nbarcode(1);
    for (size_t i=0; i<vChIA.size(); ++i) {
      int32_t n(1);
      while(i < vChIA.size()-1 && vChIA[i].v == vChIA[i+1].v) { ++i; ++n; }
      RGB color(getInterRGB((n-1)/(double)max));
      strokeChIADropBarcode(vChIA[i-n+1].v, std::to_string(n), ywidth, par.yaxis_now + (nbarcode++)*ywidth, color);
    }
  }
  par.yaxis_now += boxheight + MERGIN_BETWEEN_READ_BED;

  DEBUGprint_FUNCStart();
  return;
}
