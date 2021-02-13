/* Copyright(c)  Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <cmath>
#include <algorithm>
#include "color.hpp"

double Fraction(const double v) {
  return v - floor(v);
}

RGB HSVtoRGB(const HSV &hsv)
{
  const double h = Fraction(hsv.h);
  const double s = hsv.s;
  const double v = hsv.v;
  const double hueF = h * 6.0;
  const int hueI = static_cast<int>(hueF);
  const double fr = hueF - hueI;
  const double m = v * (1.0-s);
  const double n = v * (1.0-s*fr);
  const double p = v * (1.0-s*(1.0-fr));

  RGB rgb;

  switch (hueI) {
  case 0:  rgb.r = v; rgb.g = p; rgb.b = m; break;
  case 1:  rgb.r = n; rgb.g = v; rgb.b = m; break;
  case 2:  rgb.r = m; rgb.g = v; rgb.b = p; break;
  case 3:  rgb.r = m; rgb.g = n; rgb.b = v; break;
  case 4:  rgb.r = p; rgb.g = m; rgb.b = v; break;
  default: rgb.r = v; rgb.g = m; rgb.b = n; break;
  }

  return rgb;
}

HSV RGBtoHSV(const RGB &rgb)
{
  const double min(std::min(std::min(rgb.r, rgb.g), rgb.b));
  const double max(std::max(std::max(rgb.r, rgb.g), rgb.b));

  HSV hsv(0.0, 0.0, max);	

  const double delta(max - min);

  if (delta != 0.0) {
    hsv.s = delta / max;
    if (rgb.r == max) hsv.h = (rgb.g-rgb.b) / delta;
    else if (rgb.g == max) hsv.h = 2.0 + (rgb.b-rgb.r) / delta;
    else hsv.h = 4.0 + (rgb.r-rgb.g) / delta;
    
    hsv.h /= 6.0;

    if (hsv.h < 0.0) hsv.h += 1.0;
  }

  return hsv;
}
