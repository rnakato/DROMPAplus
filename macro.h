/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _MACRO_H_
#define _MACRO_H_

#define VERSION "3.3.0+"

enum {NUM_1K=1000,
      NUM_1M=1000000,
      NUM_10M=10000000,
      NUM_100M=100000000};

#define VALUE2WIGARRAY(v) ((v) * 1000.0)
#define WIGARRAY2VALUE(v) ((v) * (1.0/1000.0))

enum PWfile_Type {
  TYPE_BINARY,
  TYPE_COMPRESSWIG,
  TYPE_UNCOMPRESSWIG,
  TYPE_BEDGRAPH,
  TYPE_BIGWIG,
  PWFILETYPENUM
};

#define BPRINT std::cout << boost::format
#define RANGE(i, min, max) (((i) >=(min)) && ((i) <=(max)) ? 1: 0)
#define overlap(s1,e1,s2,e2) ((e1 >= s2) && (e2 >= s1))
#define PRINTERR(...) do{ std::cerr << "Error: " << __VA_ARGS__ << std::endl; std::exit(1); }while(0)

#endif /* _MACRO_H_ */
