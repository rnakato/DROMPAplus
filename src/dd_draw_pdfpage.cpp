/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */

#include "dd_draw_pdfpage.hpp"

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
