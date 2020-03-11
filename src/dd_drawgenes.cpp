/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */

#include <sstream>
#include "dd_draw.hpp"
#include "dd_draw_pdfpage.hpp"
#include "color.hpp"
#include "../submodules/SSP/common/inline.hpp"
#include "../submodules/SSP/common/util.hpp"

namespace {
  class GeneElement{
    int32_t dif;
    int32_t cnt;

  public:
    double x1, x2, xcen, xwid;
    double x_name;
    int32_t y_name, ylen;
    int32_t ybar;

    GeneElement(const genedata &m, const DParam &par,
		const double ycenter, const int32_t ty,
		int32_t &on_plus, int32_t &on_minus):
      dif (6), cnt(1),
      x1(BP2PIXEL(m.txStart - par.xstart +1)),
      x2(BP2PIXEL(m.txEnd   - par.xstart +1)),
      xcen((x1+x2)/2), xwid(x2 - x1),
      ylen(0), ybar(0)
    {
      //      std::cout << m.txStart << " "  << m.txEnd << " "  << x1 << " "<< xstart << " " << x2 << " "  << xwid << std::endl;

      if (!ty) { // SGD
	x_name = xcen - 3.25 * m.gname.length() + 6;
	if (m.strand == "+") {
	  ybar   = ycenter - dif;
	  y_name = ybar -5 - on_minus*6;
	  if (on_minus == cnt) on_minus=0; else ++on_minus;
	}
	else if (m.strand == "-") {
	  ybar   = ycenter + dif;
	  y_name = ybar +9 + on_plus*6;
	  if (on_plus == cnt) on_plus=0; else ++on_plus;
	}
	else y_name = ycenter -22;
	ylen = y_name - ycenter;
      } else {  // Others
	if (m.strand == "+") {
	  ybar = ycenter - 8 - on_minus * 8;
	  if (on_minus==7) on_minus=0; else ++on_minus;
	} else{
	  ybar = ycenter + 8 + on_plus * 8;
	  if (on_plus==7) on_plus=0; else ++on_plus;
	}
	x_name = x2 + 2;
	if (x_name < 150) x_name = 150;
	else if (x_name > OFFSET_X + par.width_draw + 60) x_name = OFFSET_X + par.width_draw + 60;
	y_name = ybar +3;
      }

    }
  };

  void ShowColorAnnotation(const Cairo::RefPtr<Cairo::Context> cr, const int32_t x, int32_t &ycen,
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
}

void PDFPage::strokeARS(const HashOfGeneDataMap &mp, const double ycenter)
{
  cr->set_line_width(0.3);
  try {
    std::vector<genedata> garray(get_garray(mp.at(rmchr(chrname)), par.xstart, par.xend));

    int32_t ars_on(0);
    int32_t on_plus(0);
    int32_t on_minus(0);
    for (auto &m: garray) {
      GeneElement g(m, par, ycenter, 0, on_plus, on_minus);

      if (m.gtype=="ARS") {
	cr->set_source_rgba(CLR_RED, 1);
	rel_yline(cr, g.xcen, ycenter -2, g.ylen +14 - 8 * ars_on);
	showtext_cr(cr, g.x_name, g.y_name +8 - ars_on*8, m.gname, 8);
	if (ars_on==2) ars_on=0; else ++ars_on;
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
    std::cerr << "Warning: " << chrname << " has no gene." << std::endl;
  }

  return;
}

void PDFPage::strokeGeneSGD(const DROMPA::Global &p, const double ycenter)
{
  DEBUGprint_FUNCStart();

  cr->set_line_width(2.5);
  int32_t ycen(ycenter-30);
  ShowColorAnnotation(cr, 50, ycen, "Coding",    CLR_BLUE);
  ShowColorAnnotation(cr, 50, ycen, "Noncoding", CLR_GREEN);
  ShowColorAnnotation(cr, 50, ycen, "rRNA",      CLR_BLACK);
  ShowColorAnnotation(cr, 50, ycen, "LTR",       CLR_PURPLE);

  try {
    std::vector<genedata> garray(get_garray(p.anno.gmp.at(rmchr(chrname)), par.xstart, par.xend));

    int32_t ars_on(0);
    int32_t on_plus(0);
    int32_t on_minus(0);
    for (auto &m: garray) {
      GeneElement g(m, par, ycenter, 0, on_plus, on_minus);

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
	if (ars_on==2) ars_on=0; else ++ars_on;
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

  DEBUGprint_FUNCend();
  return;
}

void PDFPage::strokeGene(const DROMPA::Global &p, const double ycenter)
{
  DEBUGprint_FUNCStart();

  cr->set_line_width(2.5);
  int32_t ycen(ycenter-20);

  if (p.anno.is_Anno_UCSC()){
    ShowColorAnnotation(cr, 50, ycen, "CodingRNA", CLR_BLUE);
    ShowColorAnnotation(cr, 50, ycen, "ncRNA", CLR_GREEN);
  } else {
    ShowColorAnnotation(cr, 50, ycen, "CodingRNA", CLR_BLUE);
    ShowColorAnnotation(cr, 50, ycen, "lincRNA", CLR_PINK);
    ShowColorAnnotation(cr, 50, ycen, "Antisence", CLR_GREEN);
    ShowColorAnnotation(cr, 50, ycen, "other RNA", CLR_ORANGE);
    ShowColorAnnotation(cr, 50, ycen, "Pseudo", CLR_GRAY2);
    ShowColorAnnotation(cr, 50, ycen, "Others", CLR_BLACK);
  }

  try {
    std::vector<genedata> garray(get_garray(p.anno.gmp.at(rmchr(chrname)), par.xstart, par.xend));

    double llimit(150);
    double rlimit(OFFSET_X + p.drawparam.width_draw_pixel + 60);
    int32_t on_plus(0);
    int32_t on_minus(0);
    for (auto &m: garray) {
      GeneElement g(m, par, ycenter, 1, on_plus, on_minus);

      if (isStr(m.gtype, "protein_coding")) cr->set_source_rgba(CLR_BLUE, 1);
      else if (isStr(m.gtype, "noncoding RNA")) cr->set_source_rgba(CLR_GREEN, 1); // UCSC
      else if (isStr(m.gtype, "lincRNA"))   cr->set_source_rgba(CLR_PINK, 1);
      else if (isStr(m.gtype, "antisense")) cr->set_source_rgba(CLR_GREEN, 1);
      else if (isStr(m.gtype, "RNA"))       cr->set_source_rgba(CLR_ORANGE, 1);
      else if (isStr(m.gtype, "pseudo"))    cr->set_source_rgba(CLR_GRAY2, 1);
      else                                  cr->set_source_rgba(CLR_BLACK, 1);

      // Gene body
      cr->set_line_width(1.5);
      if (g.x1 >= llimit) rel_yline(cr, g.x1, g.ybar-4, 8);
      if (g.x2 <= rlimit) rel_yline(cr, g.x2, g.ybar-4, 8);
      cr->stroke();
      cr->set_line_width(3);
      cr->move_to(std::max(g.x1, llimit), g.ybar);
      cr->line_to(std::min(g.x2, rlimit), g.ybar);
      cr->stroke();

      // Exon
      cr->set_line_width(6);
      for (int32_t i=0; i<m.exonCount; ++i) {
	double x(BP2PIXEL(m.exon[i].start - par.xstart +1));
	double xlen(std::max(1.0,m.exon[i].getlen() * par.dot_per_bp));
	if (x >= llimit && x <= rlimit) {
	  rel_xline(cr, x, g.ybar, xlen);
	  cr->stroke();
	}
      }

      // Name
      cr->set_source_rgba(CLR_BLACK, 1);
      if (p.anno.showtranscriptname) showtext_cr(cr, g.x_name, g.y_name, m.tname, 8);
      else showtext_cr(cr, g.x_name, g.y_name, m.gname, 8);
    }
  } catch (...) {
    std::cerr << "Warning: " << chrname  << " has no gene." << std::endl;
  }

  DEBUGprint_FUNCend();
  return;
}

void PDFPage::DrawGeneAnnotation(const DROMPA::Global &p)
{
  DEBUGprint_FUNCStart();

  int32_t boxheight;
  if (p.anno.getgftype() == GFTYPE_SGD) boxheight = BOXHEIGHT_GENEBOX_NOEXON;
  else boxheight = BOXHEIGHT_GENEBOX_EXON;

  double ytop(par.yaxis_now);
  double ycenter(ytop + boxheight/2);

  if (p.anno.arsfile != "") {
    DEBUGprint("DrawARS");
    strokeARS(p.anno.arsgmp, ycenter);
  }
  if (p.anno.genefile != "") {
    DEBUGprint("DrawGene");
    if (p.anno.getgftype() == GFTYPE_SGD) strokeGeneSGD(p, ycenter);
    else strokeGene(p, ycenter);
  }
  /* frame */
/*  cr->set_line_width(0.4);
  cr->rectangle(OFFSET_X, ytop, par.getXaxisLen(), boxheight);
  cr->stroke();*/

  /* genome line */
  cr->set_source_rgba(CLR_BLACK, 1);
  cr->set_line_width(1.5);
  rel_xline(cr, OFFSET_X, ycenter, par.getXaxisLen());
  cr->stroke();

  /* memory */
  stroke_xaxis(ycenter);
  par.yaxis_now += boxheight + MERGIN_BETWEEN_DATA;

  DEBUGprint_FUNCend();
  return;
}
