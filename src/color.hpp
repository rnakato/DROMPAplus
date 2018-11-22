/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _COLOR_H_
#define _COLOR_H_

#define CLR_BLACK 0, 0, 0
#define CLR_WHITE 1, 1, 1
#define CLR_GRAY 200.0/255.0, 200.0/255.0, 200.0/255.0
#define CLR_GRAY2 150.0/255.0, 150.0/255.0, 150.0/255.0
#define CLR_GRAY3 100.0/255.0, 100.0/255.0, 100.0/255.0
#define CLR_GRAY4 50.0/255.0, 50.0/255.0, 50.0/255.0
#define CLR_GREEN 50.0/255.0, 150.0/255.0, 50.0/255.0
#define CLR_GREEN2 0, 1, 0
#define CLR_GREEN3 60.0/255.0, 179.0/255.0, 113.0/255.0
#define CLR_GREEN4 140.0/255.0, 199.0/255.0, 0/255.0
#define CLR_GREEN5 0, 1, 64/255.0
#define CLR_OLIVE 128.0/255.0, 128.0/255.0, 0
#define CLR_BLUE 50.0/255.0, 50.0/255.0, 150.0/255.0
#define CLR_BLUE2 150/255, 1, 1
#define CLR_BLUE3 0, 0, 1
#define CLR_SLATEGRAY 112/255.0, 128.0/255.0, 144.0/255.0
#define CLR_RED 1, 0, 0
#define CLR_PINK 1, 50.0/255.0, 50.0/255.0
#define CLR_PINK2 1, 0/255.0, 64.0/255.0
#define CLR_BROWN 165/255.0, 42/255.0, 42/255.0
#define CLR_LIGHTCORAL 240.0/255.0, 128.0/255.0, 128.0/255.0
#define CLR_SALMON 250.0/255.0, 128.0/255.0, 114.0/255.0
#define CLR_ORANGE 1, 165.0/255.0, 0
#define CLR_DARKORANGE 1, 140.0/255.0, 0
#define CLR_YELLOW 240.0/255.0, 240.0/255.0, 168.0/255.0
#define CLR_YELLOW2 0.8, 1, 0
#define CLR_YELLOW3 1, 1, 150/255
#define CLR_PURPLE 150.0/255.0, 50.0/255.0, 150.0/255.0
#define CLR_PURPLE2 100.0/255.0, 10.0/255.0, 150.0/255.0
#define CLR_LAKEBLUE 79/255.0, 255/255.0, 246/255.0
#define CLR_DEEPSKYBLUE 0/255.0, 191/255.0, 246/255.0
#define CLR_BLUEPURPLE 110/255.0, 69/255.0, 247/255.0
#define CLR_BLUEGRAY 43/255.0, 199/255.0, 242/255.0

class RGB {
public:
  double r,g,b;
  RGB(): r(0), g(0), b(0) {}
  RGB(double r_, double g_, double b_ ): r(r_), g(g_), b(b_) {}
};

class HSV {
public:
  double h;
  double s,v;
  HSV(): h(0), s(0), v(0) {}
  HSV(double h_, double s_, double v_ ): h(h_), s(s_), v(v_) {}
};

RGB HSVtoRGB(const HSV &hsv);
HSV RGBtoHSV(const RGB &rgb);

#endif
