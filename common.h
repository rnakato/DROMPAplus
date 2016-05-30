/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _COMMON_H_
#define _COMMON_H_

#define VERSION "3.3.0+"

#define NUM_1K 1000
#define NUM_1M 1000000
#define NUM_10M 10000000
#define NUM_100M 100000000

typedef enum{
  TYPE_BINARY,
  TYPE_COMPRESSWIG,
  TYPE_UNCOMPRESSWIG,
  TYPE_BEDGRAPH,
  TYPE_BIGWIG,
  PWFILETYPENUM
} PWfile_Type;

#endif /* _COMMON_H_ */
