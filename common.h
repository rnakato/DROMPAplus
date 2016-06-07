/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _COMMON_H_
#define _COMMON_H_

#include <iostream>
#include <string>

#define VERSION "3.3.0+"
#define NUM_1K 1000
#define NUM_1M 1000000
#define NUM_10M 10000000
#define NUM_100M 100000000

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

template <class T>
class var {
  std::string name;
 public:
  T val;
 var(std::string str, T v=0): name(str), val(v) {}
  friend std::ostream &operator<<(std::ostream &stream, var v) {
    return stream << v.name << "\t" << v.val;
  }
  var operator++() {
    val++;
    return *this;
  }
};

#endif /* _COMMON_H_ */
