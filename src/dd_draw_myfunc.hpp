/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_DRAW_MYFUNC_H_
#define _DD_DRAW_MYFUNC_H_

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cairommconfig.h>
#include <cairomm/context.h>
#include <cairomm/surface.h>

#define rel_xline_double(cr, x1, y1, xlen) do{ \
    cr->move_to(x1,   y1); \
    cr->line_to(x1+xlen, y1); } while(0)
#define rel_xline(cr, x1, y1, xlen) do{ \
    cr->move_to(x1,   (int32_t)y1); \
    cr->line_to(x1+xlen, (int32_t)y1); } while(0)
#define rel_yline(cr, x1, y1, ylen) do{ \
    cr->move_to(x1, (int32_t)y1); \
    cr->line_to(x1, (int32_t)(y1+ylen)); } while(0)
#define mytriangle(x1, y1, x2, y2, x3, y3) do{ \
    cr->move_to(x1, y1);		       \
    cr->rel_line_to((x2) - (x1), (y2) - (y1));	       \
    cr->rel_line_to((x3) - (x2), (y3) - (y2));	       \
    cr->rel_line_to((x1) - (x3), (y1) - (y3));	       \
    cr->close_path(); } while(0)

inline double CalcRatio(const double c, const double i, const double r)
{
  return i ? c/i*r: 0;
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

#endif /* _DD_DRAW_MYFUNC_H_ */
