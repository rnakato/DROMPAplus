/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DD_GV_H_
#define _DD_GV_H_

class ddparam {
  string cmd;
 public:
  int nosig;
  int chr;
  int showars;
  int png;
  int rmchr;
  int vctag;
  int vitag;
  int vratio;
  int vpi;
  int vpe;
  int barnum;
  
  void setcmd(const string &str){ cmd = str;}
  string getcmd(){ return cmd; }
  ddparam () {}
};

#endif /* _DD_GV_H_ */
