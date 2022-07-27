/* Copyright(c)  Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "extendBedFormat.hpp"

bed12::bed12(std::vector<std::string> &s):
  bed(s), rgb_r(-1), rgb_g(-1), rgb_b(-1)
{
  try {
    int32_t num = s.size();
    if(num > 3)  name        = s[3];
    if(num > 4)  score       = stoi(s[4]);
    if(num > 5)  strand      = s[5];
    if(num > 6)  thickStart  = stoi(s[6]);
    if(num > 7)  thickEnd    = stoi(s[7]);
    if(num > 8) {
      std::vector<std::string> v;
      ParseLine_NoDelimCheck(v, s[8], ',');
      if(v.size() >= 3) {
//        std::cout << v[0] << "\t" << v[1] << "\t" << v[2] << "\n"
        rgb_r = stoi(v[0]);
        rgb_g = stoi(v[1]);
        rgb_b = stoi(v[2]);
 //       printf("%d, %d, %d\n", rgb_r, rgb_g, rgb_b);
      }
    }
    if(num > 9)  blockCount  = stoi(s[9]);
    if(num > 10) blockSizes  = stoi(s[10]);
    if(num > 11) blockStarts = stoi(s[11]);
  } catch (std::exception &e) {
    PRINTERR_AND_EXIT("invalid columns in bed12 format. " + std::string(e.what()));
  }
}

void InteractionSet::setAsMango(const std::string &lineStr)
{
  if (isStr(lineStr, "color")) {
    std::cerr << "Warning: Interaction does not look mango format." << std::endl;
    return;
  }
  if (isStr(lineStr, "chrom1")) return;

  std::vector<std::string> v;
  ParseLine_NoDelimCheck(v, lineStr, '\t');
  if(v.size() < 8) {
    std::cerr << "Warning: " << lineStr << " does not contain 8 columns." << std::endl;
    return;
  }
  try {
    double val(0), val_tmp(0);
    if(v.size() > 8) val_tmp = stod(v[15]); else val_tmp = stod(v[7]); // P
    if(val_tmp) val = -log10(val_tmp); else val = -log10(1e-12);

    vinter.emplace_back(bed({v[0], v[1], v[2]}),
                        bed({v[3], v[4], v[5]}),
                        val);
    maxval = std::max(val, maxval);
  } catch (std::exception &e) {
    PRINTERR_AND_EXIT(e.what());
  }
}

void InteractionSet::setAsHICCUPS(const std::string &lineStr)
{
  if (isStr(lineStr, "color")) return;
  std::vector<std::string> v;
  ParseLine_NoDelimCheck(v, lineStr, '\t');
   if(v.size() < 19) {
    std::cerr << "Warning: " << lineStr << " does not contain 19 columns." << std::endl;
    return;
  }

  try {
    double val(-log10(stod(v[17]))); // fdrDonut
    vinter.emplace_back(bed({v[0], v[1], v[2]}),
                        bed({v[3], v[4], v[5]}),
                        val);
    if(std::isfinite(val)) maxval = std::max(val, maxval);
  } catch (std::exception &e) {
    PRINTERR_AND_EXIT(e.what());
  }
}


void InteractionSet::setAsBEDPE(const std::string &lineStr)
{
  std::vector<std::string> v;
  ParseLine_NoDelimCheck(v, lineStr, '\t');
  if(v.size() < 6) {
    std::cerr << "Warning: " << lineStr << " does not contain 6 columns." << std::endl;
    return;
  }

  try {
    vinter.emplace_back(bed({v[0], v[1], v[2]}),
                        bed({v[3], v[4], v[5]}),
                        1);
    maxval = 1;
  } catch (std::exception &e) {
    PRINTERR_AND_EXIT(e.what());
  }
}

void InteractionSet::compare_bed_loop(const std::vector<bed> &bed1,
                                      const std::vector<bed> &bed2,
                                      const bool nobs)
{
  int32_t aa(0), bb(0), ab(0), an(0), bn(0), nn(0), hit_bed1(0), hit_bed2(0);

  for (auto &x: bed1) {
    if (isoverlap_asBed(x, vinter)) ++hit_bed1;
  }
  for (auto &x: bed2) {
    if (isoverlap_asBed(x, vinter)) ++hit_bed2;
  }
  for (auto &x: vinter) {
    int32_t on(0);
    x.ofirst.peakovrlpd1  = isoverlap_asloop(x.first, bed1);
    x.ofirst.peakovrlpd2  = isoverlap_asloop(x.first, bed2);
    x.osecond.peakovrlpd1 = isoverlap_asloop(x.second, bed1);
    x.osecond.peakovrlpd2 = isoverlap_asloop(x.second, bed2);

    if ((x.ofirst.peakovrlpd1 && x.osecond.peakovrlpd2) || (x.ofirst.peakovrlpd2 && x.osecond.peakovrlpd1)) {
      ++ab;
      ++on;
    }
    if (x.ofirst.peakovrlpd1 && x.osecond.peakovrlpd1) {
      ++aa;
      ++on;
    }
    if (x.ofirst.peakovrlpd2 && x.osecond.peakovrlpd2) {
      ++bb;
      ++on;
    }
    if (on) continue;
    if (x.ofirst.peakovrlpd1 || x.osecond.peakovrlpd1) ++an;
    else if (x.ofirst.peakovrlpd2 || x.osecond.peakovrlpd2) ++bn;
    else ++nn;

  }

  std::cout << "# Number: " << vinter.size() << std::endl;
  printf("# bed1-bed2: %d (%.1f%%)\n", ab, getpercent(ab, vinter.size()));
  printf("# bed1-bed1: %d (%.1f%%)\n", aa, getpercent(aa, vinter.size()));
  printf("# bed2-bed2: %d (%.1f%%)\n", bb, getpercent(bb, vinter.size()));
  printf("# bed1-none: %d (%.1f%%)\n", an, getpercent(an, vinter.size()));
  printf("# bed2-none: %d (%.1f%%)\n", bn, getpercent(bn, vinter.size()));
  printf("# none: %d (%.1f%%)\n",      nn, getpercent(nn, vinter.size()));
  printf("# bed1: %d (%.1f%%)\n", hit_bed1, getpercent(hit_bed1, bed1.size()));
  printf("# bed2: %d (%.1f%%)\n", hit_bed2, getpercent(hit_bed2, bed2.size()));

  if(!nobs) {
    for (auto &x: vinter) {
      if((x.ofirst.peakovrlpd1 && x.osecond.peakovrlpd2) || (x.ofirst.peakovrlpd2 && x.osecond.peakovrlpd1)) x.print();
    }
  }
}
