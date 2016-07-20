/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DD_OPT_H_
#define _DD_OPT_H_

#include "util.h"
#include "dd_gv.h"
#include "dd_readfile.h"

enum optstatus {OPTCHIP, OPTNORM, OPTTHRE, OPTANNO_PC, OPTANNO_GV, OPTDRAW, OPTREGION, OPTSCALE, OPTCG, OPTPD, OPTTR, OPTPROF, OPTOVERLAY, OPTOTHER};

class opt {
public:
  boost::program_options::options_description opts;
  opt(const string str): opts(str) {}
  void add(vector<optstatus> st);
};

class Command {
  opt opts;
  string desc;
  string requiredstr;
  vector<optstatus> vopts;
  boost::program_options::variables_map values;
  Param p;
  function<void(boost::program_options::variables_map &, Param &)> func;
  
  public:
  string name;

 Command(string n, string d, string r, function<void(boost::program_options::variables_map &, Param &)> _func, vector<optstatus> v): opts("Options"), desc(d), requiredstr(r), vopts(v), func(_func), name(n) {
    opts.add(v);
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
      p.gt = read_genometable(values["gt"].as<string>());

      func(values, p);

    } catch (exception &e) {
      cout << e.what() << endl;
    }
  }
};

void opt::add(vector<optstatus> st)
{
  boost::program_options::options_description o("Required",100);
  o.add_options()
    ("output,o",  boost::program_options::value<string>(),	 "Output prefix")
    ("gt",        boost::program_options::value<string>(),	 "Genome table")
    ;
  opts.add(o);

  for(auto x: st) {
    switch(x) {
    case OPTCHIP:
      {
	boost::program_options::options_description o("Input",100);
	o.add_options()
	  ("input,i",   boost::program_options::value<vector<string>>(), "Specify ChIP data, Input data and name of ChIP sample\n     (separated by ',', values except for 1 can be omitted)\n     1:ChIP   2:Input   3:name   4:peaklist   5:binsize\n     6:scale_tag   7:scale_ratio   8:scale_pvalue\n")
	  //	  ("binsize,b", boost::program_options::value<int>()->default_value(binsize), "Bin size")
	  ;
	opts.add(o);
	break;
      }
    case OPTNORM:
      {
	boost::program_options::options_description o("",100);
	o.add_options()
	  ("norm",      boost::program_options::value<int>()->default_value(1),	     "Normalization between ChIP and Input\n      0: not normalize\n      1: with total read number\n      2: with NCIS method\n")
	  ("sm",        boost::program_options::value<int>()->default_value(0),      "Smoothing width") // gausian ??
	  ;
	opts.add(o);
	break;
      }
    case OPTTHRE: 
      {
	boost::program_options::options_description o("Threshold",100);
	o.add_options()
	  ("pthre_internal", boost::program_options::value<double>()->default_value(1e-4), "p-value for ChIP internal")
	  ("pthre_enrich",   boost::program_options::value<double>()->default_value(1e-4), "p-value for ChIP/Input enrichment")
	  ("qthre",          boost::program_options::value<double>()->default_value(1),    "FDR")
	  ("ethre,e",        boost::program_options::value<double>()->default_value(2),    "IP/Input fold enrichment")
	  ("ipm",            boost::program_options::value<double>()->default_value(0),    "Read intensity of peak summit")
	  ("nosig", "Omit highlighting peak regions")
	  ("width4lmd", boost::program_options::value<int>()->default_value(100000), "Width for calculating local lambda")
	  ;
	opts.add(o);
	break;
      }
    case OPTANNO_PC:
      {
	boost::program_options::options_description o("Annotation",100);
	o.add_options()
	  ("gene,g", boost::program_options::value<string>(),	  "Gene annotation file")
	  ("gftype", boost::program_options::value<int>()->default_value(1), "Format of gene annotation\n     0: RefFlat (default)\n     1: Ensembl\n     2: GTF (for S. pombe)\n     3: SGD (for S. cerevisiae)\n")
	  ("ars",    boost::program_options::value<string>(),	  "ARS list (for yeast)")
	  ("ter",    boost::program_options::value<string>(),	  "TER list (for S.cerevisiae)")  
	  ("bed",    boost::program_options::value<vector<string>>(), "<bedfile>,<label>: Specify bed file and name (<label> can be omited)")
	  ("repeat", boost::program_options::value<string>(),	  "Display repeat annotation (RepeatMasker format)") 
	  ;
	opts.add(o);
	break;
      }
    case OPTANNO_GV:
      {
	boost::program_options::options_description o("Optional data",100);
	o.add_options()
	  ("mp",     boost::program_options::value<string>(),  	  "Mappability file")
	  ("mpthre", boost::program_options::value<double>()->default_value(0.3), "Low mappability threshold")
	  ("gap",    boost::program_options::value<string>(),	  "Specify gapped regions to be shaded")
	  ("inter",  boost::program_options::value<vector<string>>(), "<interaction file>,<label>: Specify interaction file and name (<label> can be omited)")  // FDRde iro kaeru
	  ("gc",     boost::program_options::value<string>(), 	  "Visualize GC contents graph")
	  ("gcsize", boost::program_options::value<int>()->default_value(100000), "Window size for GC contents")
	  ("gd",     boost::program_options::value<string>(), 	  "Visualize gene density (number of genes for each window)")
	  ("gdsize", boost::program_options::value<int>()->default_value(100000), "Window size for gene density")
	  ;
	opts.add(o);
	break;
      }
    case OPTDRAW:
      {
	boost::program_options::options_description o("Drawing",100);
	o.add_options()
	  ("showctag",     boost::program_options::value<int>(),    "Display ChIP read lines")
	  ("showitag",     boost::program_options::value<int>(),    "Display Input read lines (0:off 1:all 2:first one)")
	  ("showratio",    boost::program_options::value<int>(),    "Display ChIP/Input ratio (0:off 1:liner scale 2:logscale)")
	  ("showpinter",   boost::program_options::value<int>(),    "Display -log10(p) lines for ChIP internal")
	  ("showpenrich",  boost::program_options::value<int>(),    "Display -log10(p) lines for ChIP/Input enrichment")
	  ("showars",     boost::program_options::value<int>(),     "Display ARS only (do not display genes)")
	  ("ls",          boost::program_options::value<int>()->default_value(1000), "Width for each line (kb)")
	  ("lpp",         boost::program_options::value<int>()->default_value(1),    "Line number per page")
	  ("offbg",       boost::program_options::value<int>(),     "Omit background color of read lines")
	  ("offymem",     boost::program_options::value<int>(),     "Omit Y memory")
	  ("offylab",     boost::program_options::value<int>(),     "Omit Y label")
	  ("viz",         boost::program_options::value<int>()->default_value(0), "Color of read profile\n     0: normal color\n     1: semitransparent color\n")
	  ;
	opts.add(o);
	break;
      }
    case OPTREGION:
      {
	boost::program_options::options_description o("Region to draw",100);
	o.add_options()
	  ("chr",         boost::program_options::value<int>(),     "Output the specified chromosome only")
	  ("region,r",    boost::program_options::value<string>(),  "Specify genomic regions for drawing")
	  ("genefile",    boost::program_options::value<string>(),  "Specify gene loci to visualize")  
	  ("len_genefile",boost::program_options::value<int>()->default_value(50000), "extended length for each gene locus")
	  ;
	opts.add(o);
	break;
      }
    
    case OPTSCALE:
      {
	boost::program_options::options_description o("Scale for Y axis",100);
	o.add_options()
	  ("scale_tag",    boost::program_options::value<double>(), "Scale for read line")
	  ("scale_ratio",  boost::program_options::value<double>(), "Scale for fold enrichment")
	  ("scale_pvalue", boost::program_options::value<double>(), "Scale for -log10(p)")
	  ("bn",           boost::program_options::value<int>()->default_value(2),     "Number of memories of y-axis")
	  ("ystep",        boost::program_options::value<double>()->default_value(20), "Height of read line")
	  ;
	opts.add(o);
	break;
      }
    case OPTOVERLAY:
      {
	boost::program_options::options_description o("For overlay",100);
	o.add_options()
	  ("ioverlay",  boost::program_options::value<vector<string>>(),	  "Input file")
	  //	  ("binsize2",  boost::program_options::value<int>()->default_value(binsize), "Bin size")
	  ("scale_tag2",   boost::program_options::value<double>(), "Scale for read line")
	  ("scale_ratio2", boost::program_options::value<double>(), "Scale for fold enrichment")
	  ("scale_pvalue2",boost::program_options::value<double>(), "Scale for -log10(p)")
	  ;
	opts.add(o);
	break;
      }
    case OPTCG: 
      {
	boost::program_options::options_description o("CG",100);
	o.add_options()
	  ("cgthre",    boost::program_options::value<double>(), "Minimum threshold per kbp")
	  ;
	opts.add(o);
	break;
      }
    case OPTTR: 
      {
	boost::program_options::options_description o("TR",100);
	o.add_options()
	  ("tssthre",    boost::program_options::value<double>(), "")
	  ;
	opts.add(o);
	break;
      }
    case OPTPD:
      {
	boost::program_options::options_description o("PD",100);
	o.add_options()
	  ("pd",   boost::program_options::value<vector<string>>(), "Peak density file and name\n(separated by ',' <name> can be omited)")
	  ("prop",   boost::program_options::value<double>(),  "scale_tag")
	  ("pdsize", boost::program_options::value<int>()->default_value(100000), "windowsize for peak density")
	  ;
	opts.add(o);
	break;
      }
    case OPTPROF:
      {
	boost::program_options::options_description o("PROFILE AND HEATMAP",100);
	o.add_options()
	  ("ptype",   boost::program_options::value<int>(),  "Region type: 1; around TSS, 2; around TES, 3; divide gene into 100 subregions 4; around peak sites")
	  ("stype",   boost::program_options::value<int>(),  "Show type: 0; ChIP read (default) 1; ChIP/Input enrichment")
	  ("ntype",   boost::program_options::value<int>(),  "Normalization type: 0; total read 1; target regions only")
	  ("cw",      boost::program_options::value<double>()->default_value(2500), "width from the center")
	  ("maxval",   boost::program_options::value<double>(),  "Upper limit for heatmap")
	  ("offse",  "Omit the standard error in profile")
	  ("hmsort",   boost::program_options::value<int>()->default_value(1),  "Column number for sorting sites")
	  ("sortgbody",  "Sort sites by read number of gene body (default: TSS)")
	  ("pdetail",  "")
	  ;
	opts.add(o);
	break;
      }
      
    case OPTOTHER:
      {
	boost::program_options::options_description o("Others",100);
	o.add_options()
	  ("rmchr",   "Remove chromosome-separated pdf files")
	  ("png",     "Output with png format (Note: output each page separately)")
	  ("threads,p", boost::program_options::value<int>()->default_value(1), "number of threads to launch")
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

	//	chkminus<int>(values, "binsize", 0);

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
	for (auto x: {"scale_tag2","scale_ratio2","scale_pvalue2"}) chkminus<int>(values, x, 0);
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
  //  vector<string> str_wigfiletype = {"BINARY", "COMPRESSED WIG", "WIG", "BEDGRAPH", "BIGWIG"};
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
	//	BPRINT("   Input format: %1%\n")    % str_wigfiletype[values["if"].as<int>()];
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
