/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _MACRO_H_
#define _MACRO_H_

#include <boost/format.hpp>

#define BPRINT std::cout << boost::format
#define RANGE(i, min, max) (((i) >=(min)) && ((i) <=(max)) ? 1: 0)
#define overlap(s1,e1,s2,e2) ((e1 >= s2) && (e2 >= s1))
#define JOIN(a,b) (a ## b)

//#define CALCRATIO(c,i,r) ((i) ? ((c)/(double)((i)*(r))): 0)
//#define MKSTRING(str) (# str)
//#define FLUSH(args...) do{printf(args); fflush(stdout); }while(0)

#endif /* _MACRO_H_ */
