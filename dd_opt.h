/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DD_OPT_H_
#define _DD_OPT_H_

#include "util.h"
#include "dd_gv.h"
#include "dd_readfile.h"

using namespace boost::program_options;
enum optstatus {OPTCHIP, OPTNORM, OPTTHRE, OPTANNO_PC, OPTANNO_GV, OPTDRAW, OPTREGION, OPTSCALE, OPTCG, OPTPD, OPTTR, OPTPROF, OPTOVERLAY, OPTOTHER};

class opt {
public:
  options_description opts;
  opt(const string str): opts(str) {}
  void add(vector<optstatus> st, int binsize);
};

class Command {
  opt opts;
  string desc;
  string requiredstr;
  vector<optstatus> vopts;
  variables_map values;
  Param p;
  function<void(variables_map &, Param &)> func;
  
  public:
  string name;

 Command(string n, string d, string r, function<void(variables_map &, Param &)> _func, vector<optstatus> v, int binsize=50): opts("Options"), desc(d), requiredstr(r), vopts(v), func(_func), name(n) {
    opts.add(v, binsize);
  };
  void print() const {
    cout << setw(8) << " " << left << setw(12) << name
	 << left << setw(40) << desc << endl;
  }
  void printhelp() const {
    BPRINT("%1%:  %2%\n") % name % desc;
    BPRINT("Usage: drompa %1% [options] -o <output> -gt <genometable> %2%\n\n") % name % requiredstr;
    cout << opts.opts << endl;
  }
  void checkParam();
  void InitDump();
  void execute(int argc, char* argv[]) {
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
      p.gt = read_genometable(values["gt"].as<string>(), values["binsize"].as<int>());

      func(values, p);

    } catch (exception &e) {
      cout << e.what() << endl;
    }
  }
};

void opt::add(vector<optstatus> st, int binsize)
{
  options_description o("Required",100);
  o.add_options()
    ("output,o",  value<string>(),	 "Output prefix")
    ("gt",        value<string>(),	 "Genome table")
    ;
  opts.add(o);

  for(auto x: st) {
    switch(x) {
    case OPTCHIP:
      {
	options_description o("Input",100);
	o.add_options()
	  ("input,i",   value<vector<string>>(), "Specify ChIP data, Input data and name of ChIP sample\n     (separated by ',', values except for 1 can be omitted)\n     1:ChIP   2:Input   3:name   4:peaklist   5:binsize\n     6:scale_tag   7:scale_ratio   8:scale_pvalue\n")
	  ("if",        value<int>()->default_value(0),   "Input file format\n     0: Binary (.bin)\n     1: Compressed wig (.wig.gz)\n     2: Uncompressed wig (.wig)\n     3: bedGraph (.bedGraph)\n")
	  ("binsize,b", value<int>()->default_value(binsize), "Bin size")
	  ;
	opts.add(o);
	break;
      }
    case OPTNORM:
      {
	options_description o("",100);
	o.add_options()
	  ("norm",      value<int>()->default_value(1),	     "Normalization between ChIP and Input\n      0: not normalize\n      1: with total read number\n      2: with NCIS method\n")
	  ("sm",        value<int>()->default_value(0),      "Smoothing width") // gausian ??
	  ;
	opts.add(o);
	break;
      }
    case OPTTHRE: 
      {
	options_description o("Threshold",100);
	o.add_options()
	  ("pthre_internal", value<double>()->default_value(1e-4), "p-value for ChIP internal")
	  ("pthre_enrich",   value<double>()->default_value(1e-4), "p-value for ChIP/Input enrichment")
	  ("qthre",          value<double>()->default_value(1),    "FDR")
	  ("ethre,e",        value<double>()->default_value(2),    "IP/Input fold enrichment")
	  ("ipm",            value<double>()->default_value(0),    "Read intensity of peak summit")
	  ("nosig", "Omit highlighting peak regions")
	  ("width4lmd", value<int>()->default_value(100000), "Width for calculating local lambda")
	  ;
	opts.add(o);
	break;
      }
    case OPTANNO_PC:
      {
	options_description o("Annotation",100);
	o.add_options()
	  ("gene,g", value<string>(),	  "Gene annotation file")
	  ("gftype", value<int>()->default_value(1), "Format of gene annotation\n     0: RefFlat (default)\n     1: Ensembl\n     2: GTF (for S. pombe)\n     3: SGD (for S. cerevisiae)\n")
	  ("ars",    value<string>(),	  "ARS list (for yeast)")
	  ("ter",    value<string>(),	  "TER list (for S.cerevisiae)")  
	  ("bed",    value<vector<string>>(), "<bedfile>,<label>: Specify bed file and name (<label> can be omited)")
	  ("repeat", value<string>(),	  "Display repeat annotation (RepeatMasker format)") 
	  ;
	opts.add(o);
	break;
      }
    case OPTANNO_GV:
      {
	options_description o("Optional data",100);
	o.add_options()
	  ("mp",     value<string>(),  	  "Mappability file")
	  ("mpthre", value<double>()->default_value(0.3), "Low mappability threshold")
	  ("gap",    value<string>(),	  "Specify gapped regions to be shaded")
	  ("inter",  value<vector<string>>(), "<interaction file>,<label>: Specify interaction file and name (<label> can be omited)")  // FDRde iro kaeru
	  ("gc",     value<string>(), 	  "Visualize GC contents graph")
	  ("gcsize", value<int>()->default_value(100000), "Window size for GC contents")
	  ("gd",     value<string>(), 	  "Visualize gene density (number of genes for each window)")
	  ("gdsize", value<int>()->default_value(100000), "Window size for gene density")
	  ;
	opts.add(o);
	break;
      }
    case OPTDRAW:
      {
	options_description o("Drawing",100);
	o.add_options()
	  ("showctag",     value<int>(),    "Display ChIP read lines")
	  ("showitag",     value<int>(),    "Display Input read lines (0:off 1:all 2:first one)")
	  ("showratio",    value<int>(),    "Display ChIP/Input ratio (0:off 1:liner scale 2:logscale)")
	  ("showpinter",   value<int>(),    "Display -log10(p) lines for ChIP internal")
	  ("showpenrich",  value<int>(),    "Display -log10(p) lines for ChIP/Input enrichment")
	  ("showars",     value<int>(),     "Display ARS only (do not display genes)")
	  ("ls",          value<int>()->default_value(1000), "Width for each line (kb)")
	  ("lpp",         value<int>()->default_value(1),    "Line number per page")
	  ("offbg",       value<int>(),     "Omit background color of read lines")
	  ("offymem",     value<int>(),     "Omit Y memory")
	  ("offylab",     value<int>(),     "Omit Y label")
	  ("viz",         value<int>()->default_value(0), "Color of read profile\n     0: normal color\n     1: semitransparent color\n")
	  ;
	opts.add(o);
	break;
      }
    case OPTREGION:
      {
	options_description o("Region to draw",100);
	o.add_options()
	  ("chr",         value<int>(),     "Output the specified chromosome only")
	  ("region,r",    value<string>(),  "Specify genomic regions for drawing")
	  ("genefile",    value<string>(),  "Specify gene loci to visualize")  
	  ("len_genefile",value<int>()->default_value(50000), "extended length for each gene locus")
	  ;
	opts.add(o);
	break;
      }
    
    case OPTSCALE:
      {
	options_description o("Scale for Y axis",100);
	o.add_options()
	  ("scale_tag",    value<double>(), "Scale for read line")
	  ("scale_ratio",  value<double>(), "Scale for fold enrichment")
	  ("scale_pvalue", value<double>(), "Scale for -log10(p)")
	  ("bn",           value<int>()->default_value(2),     "Number of memories of y-axis")
	  ("ystep",        value<double>()->default_value(20), "Height of read line")
	  ;
	opts.add(o);
	break;
      }
    case OPTOVERLAY:
      {
	options_description o("For overlay",100);
	o.add_options()
	  ("ioverlay",  value<vector<string>>(),	  "Input file")
	  ("binsize2",  value<int>()->default_value(binsize), "Bin size")
	  ("scale_tag2",   value<double>(), "Scale for read line")
	  ("scale_ratio2", value<double>(), "Scale for fold enrichment")
	  ("scale_pvalue2",value<double>(), "Scale for -log10(p)")
	  ;
	opts.add(o);
	break;
      }
    case OPTCG: 
      {
	options_description o("CG",100);
	o.add_options()
	  ("cgthre",    value<double>(), "Minimum threshold per kbp")
	  ;
	opts.add(o);
	break;
      }
    case OPTTR: 
      {
	options_description o("TR",100);
	o.add_options()
	  ("tssthre",    value<double>(), "")
	  ;
	opts.add(o);
	break;
      }
    case OPTPD:
      {
	options_description o("PD",100);
	o.add_options()
	  ("pd",   value<vector<string>>(), "Peak density file and name\n(separated by ',' <name> can be omited)")
	  ("prop",   value<double>(),  "scale_tag")
	  ("pdsize", value<int>()->default_value(100000), "windowsize for peak density")
	  ;
	opts.add(o);
	break;
      }
    case OPTPROF:
      {
	options_description o("PROFILE AND HEATMAP",100);
	o.add_options()
	  ("ptype",   value<int>(),  "Region type: 1; around TSS, 2; around TES, 3; divide gene into 100 subregions 4; around peak sites")
	  ("stype",   value<int>(),  "Show type: 0; ChIP read (default) 1; ChIP/Input enrichment")
	  ("ntype",   value<int>(),  "Normalization type: 0; total read 1; target regions only")
	  ("cw",      value<double>()->default_value(2500), "width from the center")
	  ("maxval",   value<double>(),  "Upper limit for heatmap")
	  ("offse",  "Omit the standard error in profile")
	  ("hmsort",   value<int>()->default_value(1),  "Column number for sorting sites")
	  ("sortgbody",  "Sort sites by read number of gene body (default: TSS)")
	  ("pdetail",  "")
	  ;
	opts.add(o);
	break;
      }
      
    case OPTOTHER:
      {
	options_description o("Others",100);
	o.add_options()
	  ("rmchr",   "Remove chromosome-separated pdf files")
	  ("png",     "Output with png format (Note: output each page separately)")
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
  for (auto x: {"output", "gt"}) if (!values.count(x)) PRINTERR("specify --" << x << " option.");

  for(auto x: vopts) {
    switch(x) {
    case OPTCHIP:
      {
	for (auto x: {"input"}) if (!values.count(x)) PRINTERR("specify --" << x << " option.");
	
	chkrange<int>(values, "if", 0, 3);
	chkminus<int>(values, "binsize", 0);

	vector<string> v(values["input"].as<vector<string>>());
	for(auto x:v) scan_samplestr(x, p.sample, p.samplepair);
	
	break;
      }
    case OPTNORM:
      {
	chkminus<int>(values, "sm", 0);
	chkrange<int>(values, "norm", 0, 1);
	break;
      }
    case OPTTHRE: 
      {
	for (auto x: {"pthre_internal", "pthre_enrich", "qthre", "ipm", "ethre"}) chkminus<int>(values, x, -1);
	chkminus<int>(values, "width4lmd", 0);
	break;
      }
    case OPTANNO_PC:
      {
	chkrange<int>(values, "gftype", 0, 3);
	for (auto x: {"gene", "ars", "ter"}) if (values.count(x)) isFile(values[x].as<string>());
	break;
      }
    case OPTANNO_GV:
      {
	chkminus<int>(values, "mpthre", -1);
	for (auto x: {"gcsize", "gdsize"}) chkminus<int>(values, x, 0);
	break;
      }
    case OPTDRAW:
      {
	for (auto x: {"ls", "lpp"}) chkminus<int>(values, x, 0);
	break;
      }
    case OPTREGION:
      {
	for (auto x: {"region", "genefile"}) if (values.count(x)) isFile(values[x].as<string>());
	chkminus<int>(values, "len_genefile", -1);
	break;
      }
    case OPTSCALE:
      {
	for (auto x: {"scale_tag","scale_ratio","scale_pvalue","bn","ystep"}) chkminus<int>(values, x, 0);
	break;
      }
    case OPTOVERLAY:
      {
	for (auto x: {"scale_tag2","scale_ratio2","scale_pvalue2","binsize2"}) chkminus<int>(values, x, 0);
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
	for (auto x: {"pd"}) if (!values.count(x)) PRINTERR("specify --" << x << " option.");
	for (auto x: {"pdsize"}) chkminus<int>(values, x, 0);
	
	vector<string> v(values["pd"].as<vector<string>>());
	for(auto x:v) p.pd.push_back(scan_pdstr(x));
	break;
      }
    case OPTPROF:
      {
	chkrange<int>(values, "ptype", 0, 4);
	chkrange<int>(values, "stype", 0, 1);
	chkrange<int>(values, "ntype", 0, 1);
	for (auto x: {"cw", "maxval", "hmsort"}) chkminus<int>(values, x, 0);
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
  vector<string> str_bool = {"ON", "OFF"};
  vector<string> str_gftype = {"refFlat", "Ensembl", "gtf", "SGD"};
  vector<string> str_wigfiletype = {"BINARY", "COMPRESSED WIG", "WIG", "BEDGRAPH", "BIGWIG"};
  vector<string> str_norm  = { "OFF", "TOTALREAD", "NCIS" };
  vector<string> str_stype = { "ChIP read", "Enrichment ratio", "Enrichment P-value" };
  vector<string> str_ptype = { "NONE", "TSS", "TTS", "GENE100", "SPECIFIEDSITES" };
  vector<string> str_ntype = { "WHOLE GENOME", "TARGET REGIONS ONLY" };

  BPRINT("\n======================================\n");
  BPRINT("drompa version %1%: %2%\n\n") % VERSION % name;
  BPRINT("output prefix: %1%\n")     % values["output"].as<string>();
  BPRINT("Genome-table file: %1%\n") % values["gt"].as<string>();

  for(auto x: vopts) {
    switch(x) {
    case OPTCHIP:
      {
	BPRINT("\nSamples\n");
	for(uint i=0; i<p.samplepair.size(); ++i) {
	  cout << (i+1) << ": ";
	  p.samplepair[i].print();
	}
	BPRINT("   Input format: %1%\n")    % str_wigfiletype[values["if"].as<int>()];
	break;
      }
    case OPTNORM:
      {
	BPRINT("   ChIP/Input normalization: %s\n") % str_norm[values["norm"].as<int>()];
	if(values["sm"].as<int>()) BPRINT("   smoothing width: %1% bp\n") % values["sm"].as<int>();
	break;
      }
    case OPTTHRE: 
      {
	BPRINT("   Peak intensity threshold: %1$.2f\n")               % values["IPmaxthre"].as<double>();
	BPRINT("   Enrichment threshold: %1$.2f\n")                   % values["enrichthre"].as<double>();
	BPRINT("   p-value threshold (internal, -log10): %1$.2e\n")   % values["pthre_internal"].as<double>();
	BPRINT("   p-value threshold (enrichment, -log10): %1$.2e\n") % values["pthre_enrich"].as<double>();
	BPRINT("   FDR threshold: %1$.2e\n")                          % values["FDR"].as<double>();
  	break;
      }
    case OPTANNO_PC:
      {
	BPRINT("\nAnnotations:\n");
	if(values.count("gene")) BPRINT("   Gene file: %1%, Format: %2%\n")
			     % values["gene"].as<string>() % str_gftype[values["gftype"].as<int>()];
	printOpt<string>(values, "ars",    "   ARS file");
	printOpt<string>(values, "ter",    "   TER file");
	printOpt<string>(values, "repeat", "   Repeat file");
	printOpt<string>(values, "gc", "   GCcontents file");
	printOpt<string>(values, "gd", "   Gene density file");
	/*	if(d->arsfile)     BPRINT("   ARS file: %1%\n")          % values["ars"].as<string>();
	if(d->terfile)     BPRINT("   TER file: %1%\n")          % values["ter"].as<string>();
	if(d->repeat.argv) BPRINT("   Repeat file: %1%\n")       % values["repeat"].as<string>();*/
	printOpt<string>(values, "region", "   Region file");
	printVOpt<string>(values, "bed", "   Bed file");
	//	if(name != "PROFILE" || name != "HEATMAP") BPRINT("   name: %1%\n") % d->bed[i]->name;
	break;
      }
    case OPTANNO_GV:
      {
	BPRINT("\nAnnotations:\n");
	if(values.count("gene")) BPRINT("   Gene file: %1%, Format: %2%\n")
			     % values["gene"].as<string>() % str_gftype[values["gftype"].as<int>()];
	printOpt<string>(values, "gc", "   GCcontents file");
	printOpt<string>(values, "gd", "   Gene density file");
	/*	
	if(d->GC.argv)     BPRINT("   GCcontents file: %1%\n")   % values["gc"].as<string>();
	if(d->GD.argv)     BPRINT("   Gene density file: %1%\n") % values["gd"].as<string>();*/
	printVOpt<string>(values, "inter", "   Interaction file");
	if (values.count("mp")) {
	  BPRINT("Mappability file directory: %1%\n") % values["mp"].as<string>();
	  BPRINT("\tLow mappablitiy threshold: %1%\n") % values["mpthre"].as<double>();
	}
	break;
      }
    case OPTDRAW:
      {
	BPRINT("\nFigure parameter:\n");
	BPRINT("   Display read: ChIP %1%, Input %2%\n") % str_bool[values["showctag"].as<int>()] % str_bool[values["showitag"].as<int>()];
	BPRINT("   Display enrichment: %1%\n")           % str_bool[values["showratio"].as<int>()];
	BPRINT("   Display pvalue (internal): %1%\n")    % str_bool[values["showpinter"].as<int>()];
	BPRINT("   Display pvalue (ChIP/Input): %1%\n")  % str_bool[values["showpenrich"].as<int>()];
	BPRINT("   Background color: %1%\n")             % str_bool[!values["offbg"].as<int>()];
	BPRINT("   Y label: %1%\n")                      % str_bool[!values["offylab"].as<int>()];
	BPRINT("   Y memory: %1%\n")                     % str_bool[!values["offymem"].as<int>()];
	break;
      }
    
    case OPTREGION:
      {
	break;
      }
    case OPTSCALE:
      {
	for (auto x: {"scale_tag","scale_ratio","scale_pvalue","scale_tag2","scale_ratio2","scale_pvalue2","bn","ystep"}) chkminus<int>(values, x, 0);
	break;
      }
    case OPTOVERLAY:
      {
	if(values.count("ioverlay")) { 
	  BPRINT("\nOverlayed samples\n");
	  for(uint i=0; i<p.samplepair.size(); ++i) {
	    cout << (i+1) << ": ";
	    p.samplepair[i].print();
	  }
	}
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
	BPRINT("\nSamples\n");
	for(uint i=0; i<p.pd.size(); ++i) {
	  BPRINT("   IP%1%: %2%\tname: %3%\n") % (i+1) % p.pd[i].argv % p.pd[i].name;
	}
	break;
      }
    case OPTPROF:
      {
	BPRINT("   show type: %1$\n")             % str_stype[values["stype"].as<int>()];
	BPRINT("   profile type: %1$\n")          % str_ptype[values["ptype"].as<int>()];
	BPRINT("   profile normalization: %1$\n") % str_ntype[values["ntype"].as<int>()];
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
  
  if(values.count("chr")) BPRINT("output %1% only.\n") % values["chr"].as<int>();

  printf("======================================\n");
  return;
}

#endif /* _DD_OPT_H_ */
