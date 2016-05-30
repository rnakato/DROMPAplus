#include <cstdlib>
#include "warn.h"

void printerr(string str){
  cerr << "Error: " << str << endl;
  exit(1);
}
void printwrn(string str){
  cerr << "Warning: " << str << endl;
}
