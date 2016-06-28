/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DD_OPT_H_
#define _DD_OPT_H_

#include "util.h"
#include "common.h"
#include "macro.h"

using namespace boost::program_options;
enum optstatus {OPTREQ, OPTIO, OPTTHRE, OPTANNO, OPTDRAW, OPTSCALE, OPTCG, OPTPD, OPTTR, OPTPROF, OPTOTHER};

class opt {
public:
  options_description opts;
  opt(const string str): opts(str){}
  void add(vector<optstatus> st);
};

class Command {
  opt opts;
  public:
  string name;
  string desc;
  vector<optstatus> vopts;
  variables_map values;
  
 Command(string n, string d, vector<optstatus> v): opts("Options"), name(n), desc(d), vopts(v) {
    opts.add(v);
  };
  void print() {
    cout << setw(8) << " " << left << setw(12) << name
	 << left << setw(40) << desc << endl;
  }
  void printhelp() {
    BPRINT("%1%:  %2%\n") % name % desc;
    BPRINT("Usage: drompa %1% [options] -o <output> -gt <genometable> -i <ChIP>,<Input>,<name> [-i <ChIP>,<Input>,<name> ...]\n\n") % name;
    cout << opts.opts << endl;
  }
  void checkParam();
  void InitDump();
  void getOpts(int argc, char* argv[]) {
    if (argc ==1) {
      printhelp();
      exit(0);
    }
    try {
      store(parse_command_line(argc, argv, opts.opts), values);
      
      if (values.count("help")) {
	printhelp();
	exit(0);
      }
      
      notify(values);
      checkParam();
      InitDump();
    } catch (exception &e) {
      cout << e.what() << endl;
    }
  }
};

void opt::add(vector<optstatus> st)
{
  for(auto x: st) {
    switch(x) {
    case OPTREQ:
      {
	options_description o("Required");
	o.add_options()
	  ("input,i",   value<vector<string>>(), "specify ChIP data, Input data and name of ChIP sample (separated by ',', <Input> and <name> can be omitted)")
	  ("output,o",  value<string>(),	  "output prefix")
	  ("gt",        value<string>(),	  "genome table")
	  ;
	opts.add(o);
      break;
      }
    case OPTIO:
      {
	options_description o("Input/Output");
	o.add_options()
	  ("ioverlay",  value<vector<string>>(),	      "input file to overlay")
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
	  ("qthre",     value<double>(),      "FDR")
	  ("ipm",       value<double>(),      "IPmaxthre")
	  ("ethre,e",   value<double>(),      "enrichthre")
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
	  ("sortgbody",value<double>(),  "")
	  ("pdetail",   value<double>(),  "")
	  ;
	opts.add(o);
	break;
      }
      
    case OPTOTHER:
      {
	options_description o("Others");
	o.add_options()
	  ("threads,p", value<int>()->default_value(1), "number of threads to launch")
	  ("help,h", "show help message")
	  ;
	opts.add(o);
	break;
      }
    }
  }
}

void Command::checkParam() {
  for(auto x: vopts) {
    switch(x) {
    case OPTREQ:
      {
	for (auto x: {"input", "output", "gt"}) if (!values.count(x)) PRINTERR("specify --" << x << " option.");
	break;
      }
    case OPTIO:
      {
	for (auto x: {"binsize", "binsize2", "sm", "width4lmd"}) chkminus<int>(values, x, 0);
	//	if(!RANGE(values["if"].as<int>(), 0, )) printerr("invalid wigfile type.\n");
	//if(!RANGE(values["norm"].as<int>(), 0, )) printerr("invalid wigfile type.\n");
	break;
      }
    case OPTTHRE: 
      {
	for (auto x: {"pthre_internal", "pthre_enrich", "qthre", "ipm", "ethre"}) chkminus<int>(values, x, -1);
	break;
      }
    case OPTANNO:
      {
	break;
      }
    case OPTDRAW:
      {
	for (auto x: {"pdsize", "len_genefile", "ls", "lpp"}) chkminus<int>(values, x, 0);
	break;
      }
    case OPTSCALE:
      {
	for (auto x: {"scale_tag","scale_ratio","scale_pvalue","scale_tag2","scale_ratio2","scale_pvalue2","bn","ystep"}) chkminus<int>(values, x, 0);
	break;
      }
    case OPTCG: 
      {
	for (auto x: {"cgthre"}) chkminus<int>(values, x, -1);
	break;
      }
    case OPTTR: 
      {
	for (auto x: {"tssthre"}) chkminus<int>(values, x, -1);
	break;
      }
    case OPTPD:
      {
	break;
      }
    case OPTPROF:
      {
	break;
      }
      
    case OPTOTHER:
      {
	for (auto x: {"threads"}) chkminus<int>(values, x, 0);
	break;
      }
    }
    return;
  }
}

void Command::InitDump()
{
  uint i;
  vector<string> str_wigfiletype = {"BINARY", "COMPRESSED WIG", "WIG", "BEDGRAPH", "BIGWIG"};
  vector<string> str_norm  = { "OFF", "TOTALREAD", "NCIS" };
  vector<string> str_stype = { "ChIP read", "Enrichment ratio", "Enrichment P-value" };
  vector<string> str_ptype = { "NONE", "TSS", "TTS", "GENE100", "SPECIFIEDSITES" };
  vector<string> str_ntype = { "WHOLE GENOME", "TARGET REGIONS ONLY" };

  BPRINT("\n======================================\n");
  BPRINT("drompa version %1%: %2%\n\n") % VERSION % name;
  BPRINT("output prefix: %1%\n")     % values["output"].as<string>();
  BPRINT("Genome-table file: %1%\n") % values["gt"].as<string>();
  BPRINT("\tInput format: %1%\n")    % str_wigfiletype[values["if"].as<int>()];
  if(values["sm"].as<int>()) BPRINT("   smoothing width: %1% bp\n") % values["sm"].as<int>();
  if (values.count("mp")) {
    BPRINT("Mappability file directory: %1%\n") % values["mp"].as<string>();
    BPRINT("\tLow mappablitiy threshold: %1%\n") % values["mpthre"].as<double>();
  }
  BPRINT("   ChIP/Input normalization: %s\n") %str_norm[values["norm"].as<int>()];

  
  for (auto x:) {
    BPRINT("Input file %1%\n") % x;
  }

  if(p->ftype != FTYPE_GOVERLOOK){
    BPRINT("\nSamples\n");
    vector<string> vsample(values["input"].as<vector<string>>());
    for(i=0; i<vsample.size(); ++i) {
      BPRINT("   ChIP%1%: %2%\tname: %3%\tbinsize:%4%\n") % i+1 % sample[i].ChIP->argv % sample[i].linename % sample[i].binsize;
      if(sample[i].Input->argv) BPRINT("   Input%1%: %2%\n") % i+1 % sample[i].Input->argv;
      if(sample[i].peak_argv) BPRINT("      peak list: %1%\n") % sample[i].peak_argv;
    }
    vector<string> vsample(values["ioverlay"].as<vector<string>>());
    for(; i<p->samplenum; i++){
      BPRINT("   Overlayed ChIP%1%: %2%\tname: %3%\tbinsize:%4%\n") % i+1 % sample[i].ChIP->argv % sample[i].linename % sample[i].binsize;
      if(sample[i].Input->argv) BPRINT("   Input%1%: %2%\n") % i+1 % sample[i].Input->argv;
      if(sample[i].peak_argv) BPRINT("      peak list: %1%\n") % sample[i].peak_argv;
    }
  }
  for(auto x: vopts) {
    switch(x) {
    case OPTREQ:
      {
	break;
      }
    case OPTIO:
      {
	break;
      }
    case OPTTHRE: 
      {
	BPRINT("   Peak intensity threshold: %1$.2f\n") % values["IPmaxthre"].as<double>();
	BPRINT("   Enrichment threshold: %1$.2f\n")     % values["enrichthre"].as<double>();
	BPRINT("   p-value threshold (internal, -log10): %1$.2e\n")   % values["pthre_internal"].as<double>();
	BPRINT("   p-value threshold (enrichment, -log10): %1$.2e\n") % values["pthre_enrich"].as<double>();
	BPRINT("   FDR threshold: %1$.2e\n") % values["FDR"].as<double>();
  	break;
      }
    case OPTANNO:
      {
	BPRINT("\nAnnotations:\n");
	if(d->gene.argv)   BPRINT("   gene file: %s, format: %1%\n") %d->gene.argv, str_gftype[d->gftype];
	if(d->arsfile)     BPRINT("   ARS file: %1%\n") % d->arsfile;
	if(d->terfile)     BPRINT("   TER file: %1%\n") % d->terfile;
	if(d->repeat.argv) BPRINT("   repeat file: %1%\n") %d->repeat.argv;
	if(d->GC.argv)     BPRINT("   GCcontents file: %1%\n") %d->GC.argv;
	if(d->GD.argv)     BPRINT("   gene density file: %1%\n") %d->GD.argv;
	if(d->internum){
	  for(i=0; i<d->internum; i++) BPRINT("   interaction file %d: %1%\n",i+1, d->inter[i].argv;
					      }
	  if(d->drawregion_argv) BPRINT("   region file: %1%\n") %d->drawregion_argv;
	  for(i=0; i<d->bednum; i++){
	    if(p->ftype==FTYPE_PROFILE || p->ftype==FTYPE_HEATMAP) BPRINT("   bedfile%d: %1%\n") %i+1, d->peak[i]->argv;
	    else BPRINT("   bedfile%d: %s, name: %1%\n") %i+1, d->bed[i]->argv, d->bed[i]->name;
	  }
	  break;
      }
    case OPTDRAW:
      {
	BPRINT("\nFigure parameter:\n");
	BPRINT("   display read: ChIP %s, Input %1%\n") %str_bool[d->visualize_ctag], str_Inputread[d->visualize_itag];
	BPRINT("   display enrichment: %1%\n") %str_ratio[d->visualize_ratio];
	BPRINT("   display pvalue (internal): %1%\n") %str_bool[d->visualize_p_inter];
	BPRINT("   display pvalue (ChIP/Input): %1%\n") %str_bool[d->visualize_p_enrich];
	BPRINT("   background color: %1%\n") %str_bool[d->backcolors];
	BPRINT("   Y label: %1%\n") %str_bool[d->stroke_ylab];
	BPRINT("   Y memory: %1%\n") %str_bool[d->stroke_ymem];
	break;
      }
    case OPTSCALE:
      {
	for (auto x: {"scale_tag","scale_ratio","scale_pvalue","scale_tag2","scale_ratio2","scale_pvalue2","bn","ystep"}) chkminus<int>(values, x, 0);
	break;
      }
    case OPTCG: 
      {
	for (auto x: {"cgthre"}) chkminus<int>(values, x, -1);
	break;
      }
    case OPTTR: 
      {
	for (auto x: {"tssthre"}) chkminus<int>(values, x, -1);
	break;
      }
    case OPTPD:
      {
	//	for(i=0; i<d->pdnum; i++) BPRINT("   IP%d: %s\tname: %s\n") %i+1, d->PD[i].argv, d->PD[i].name;
	//	BPRINT("   pdsize: %d\n") %d->pdsize;
	break;
      }
    case OPTPROF:
      {
	BPRINT("   show type: %1$\n")             %str_stype[values["stype"].as<int>()];
	BPRINT("   profile type: %1$\n")          %str_ptype[values["ptype"].as<int>()];
	BPRINT("   profile normalization: %1$\n") %str_ntype[values["ntype"].as<int>()];
	break;
      }
      
    case OPTOTHER:
      {
	for (auto x: {"threads"}) chkminus<int>(values, x, 0);
	break;
      }
    }
    return;
  }

  }
  
  if(d->chronly) BPRINT("output %1% only.\n") % values["chr"].as<int>();

  printf("======================================\n");
  return;
}

#endif /* _DD_OPT_H_ */
