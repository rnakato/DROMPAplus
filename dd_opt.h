/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DD_OPT_H_
#define _DD_OPT_H_

using namespace boost::program_options;
enum optstatus {OPTREQ, OPTIO, OPTTHRE, OPTANNO, OPTDRAW, OPTSCALE, OPTCG, OPTPD, OPTTR, OPTPROF, OPTOTHER};

class opt {
public:
  options_description opts;
  opt(const string str): opts(str){}
  void add(vector<optstatus> st);
};

void opt::add(vector<optstatus> st)
{
  for(auto x: st) {
    switch(x) {
    case OPTREQ:
      {
	options_description o("Required");
	o.add_options()
	  ("input,i",   value<vector<string> >(), "specify ChIP data, Input data and name of ChIP sample (separated by ',', <Input> and <name> can be omitted)")
	  ("output,p",  value<string>(),	  "output prefix")
	  ("gt",        value<string>(),	  "genome table")
	  ;
	opts.add(o);
      break;
      }
    case OPTIO:
      {
	options_description o("Input/Output");
	o.add_options()
	  ("ioverlay",  value<vector<string> >(),	      "input file to overlay")
	  ("binsize,b", value<int>()->default_value(100), "bin size")
	  ("binsize2",  value<int>()->default_value(100), "bin size for overlay")
	  ("if",        value<int>(),		      "itype")
	  ("sm",        value<int>(),		      "smoothing width") // gausian ??
	  ("width4lmd", value<int>(),		      "width4lambda")
	  ("norm",      value<int>(),	              "ntype")
	  ;
	opts.add(o);
	break;
      }
    case OPTTHRE: 
      {
	options_description o("Threshold");
	o.add_options()
	  ("pthre_internal", value<double>(), "pthre_internal")
	  ("pthre_enrich",   value<double>(), "pthre_enrich")
	  ("qthre",     value<double>(),	  "FDR")
	  ("ipm",       value<double>(),	  "IPmaxthre")
	  ("ethre,e",   value<double>(),	  "enrichthre")
	  ;
	opts.add(o);
	break;
      }
    case OPTANNO:
      {
	options_description o("Annotation");
	o.add_options()
	  ("gene,g", value<string>(),	  "gene")
	  ("gftype", value<int>(),		  "gftype")
	  ("ars",    value<string>(),	  "ARS")  
	  ("ter",    value<string>(),	  "TER")
	  ("mp",     value<string>(),  	  "mappability")
	  ("mpthre", value<double>()->default_value(0.3), "low mappability threshold")
	  ("gap",    value<string>(),	  "gap")
	  ("bed,b",  value<vector<string> >(),"bed")
	  ("inter",  value<vector<string> >(),"interaction")  // FDRde iro kaeru
	  ("gc",     value<string>(), 	  "GC")
	  ("gcsize", value<int>(),	 	  "gcsize")
	  ("gd",     value<string>(), 	  "GD")
	  ("gdsize", value<int>(),		  "gdsize")
	  ("repeat", value<string>(),	  "repeat") 
	  ;
	opts.add(o);
	break;
      }
    case OPTDRAW:
      {
	options_description o("Drawing");
	o.add_options()
	  ("nosig",       value<int>(),     "nosig")
	  ("chr",         value<int>(),     "chronly")
	  ("region,r",    value<string>(),	"regions for drawing")
	  ("genefile",    value<string>(),	"genefile")  
	  ("pdsize",      value<int>(),	"pdsize")
	  ("len_genefile",value<int>(),	"length")
	  ("showars",     value<int>(),     "nosig")
	  ("png",         value<int>(),     "png format")
	  ("rmchr",       value<int>(),     "rmchr")
	  ("ls",          value<int>()->default_value(1000),"line size")
	  ("lpp",         value<int>()->default_value(1),	"line num per page")
	  ("offbg",       value<int>(),     "png format")
	  ("offymem",     value<int>(),     "png format")
	  ("offylab",     value<int>(),     "png format")
	  ("viz",         value<int>(),     "png format")
	  ;
	opts.add(o);
	break;
      }
    case OPTSCALE:
      {
	options_description o("Scaling");
	o.add_options()
	  ("scale_tag",    value<double>(), "scale_tag")
	  ("scale_ratio",  value<double>(), "scale_ratio")
	  ("scale_pvalue", value<double>(),	"scale_pvalue")
	  ("scale_tag2",   value<double>(),	"scale_tag2")
	  ("scale_ratio2", value<double>(),	"scale_ratio2")
	  ("scale_pvalue2",value<double>(),	"scale_pvalue2")
	  ("showctag",     value<int>(),    "ctag")
	  ("showitag",     value<int>(),    "itag")
	  ("showratio",    value<int>(),    "ratop")
	  ("showpinter",   value<int>(),    "itag")
	  ("showpenrich",  value<int>(),    "itag")
	  ("bn",           value<int>(),    "itag")
	  ("ystep",        value<double>(),	"ystep")
	  ;
	opts.add(o);
	break;
      }
    case OPTCG: 
      {
	options_description o("CG");
	o.add_options()
	  ("cgthre",    value<double>(), "(for CG) ")
	  ;
	opts.add(o);
	break;
      }
    case OPTTR: 
      {
	options_description o("TR");
	o.add_options()
	  ("tssthre",    value<double>(), "(for CG) ")
	  ;
	opts.add(o);
	break;
      }
    case OPTPD:
      {
	options_description o("PD");
	o.add_options()
	  ("prop",   value<double>(),  "scale_tag")
	  ;
	opts.add(o);
	break;
      }
    case OPTPROF:
      {
	options_description o("PROFILE AND HEATMAP");
	o.add_options()
	  ("ptype",   value<double>(),  "")
	  ("stype",   value<double>(),  "")
	  ("ntype",   value<double>(),  "")
	  ("cw",      value<double>(),  "")
	  ("offse",   value<double>(),  "")
	  ("hmsort",   value<double>(),  "")
	  ("sortgbody",   value<double>(),  "")
	  ("pdetail",   value<double>(),  "")
	  ;
	opts.add(o);
	break;
      }
      
    case OPTOTHER:
      {
	options_description o("Others");
	o.add_options()
	  ("help,h", "show help message")
	  ;
	opts.add(o);
	break;
      }
    }
  }
}

//             PC_BROAD    peak-calling (for broad mode)
//             FRIP        accumulate read counts in bed regions specified
//             3DMAP       accumulate read counts in bed regions specified

int func_PCSHARP(int argc, char* argv[]);
int func_ENRICH(int argc, char* argv[]);
int func_GV(int argc, char* argv[]);
int func_PD(int argc, char* argv[]);
int func_CI(int argc, char* argv[]);
int func_CG(int argc, char* argv[]);
int func_GOVERLOOK(int argc, char* argv[]);
int func_PROFILE(int argc, char* argv[]);
int func_HEATMAP(int argc, char* argv[]);
int func_TR(int argc, char* argv[]);


#endif /* _DD_OPT_H_ */
