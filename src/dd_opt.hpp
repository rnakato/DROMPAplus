/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _DD_OPT_HPP_
#define _DD_OPT_HPP_

#include "SSP/common/util.hpp"
#include "dd_gv.hpp"
#include "dd_readfile.hpp"
#include "version.hpp"

std::vector<chrsize> read_genometable(const std::string&);

enum class DrompaOpt {CHIP, NORM, THRE, ANNO_PC, ANNO_GV, DRAW, REGION, SCALE, CG, PD, TR, PROF, OVERLAY, OTHER};

class opt {
public:
  boost::program_options::options_description opts;
  opt(const std::string str): opts(str) {}
  void add(std::vector<DrompaOpt> st);
};

class Command {
  opt opts;
  std::string desc;
  std::string requiredstr;
  std::vector<DrompaOpt> vopts;
  boost::program_options::variables_map values;
  Param p;
  std::function<void(boost::program_options::variables_map &, Param &)> func;
  
  public:
  std::string name;

 Command(std::string n, std::string d, std::string r, std::function<void(boost::program_options::variables_map &, Param &)> _func, std::vector<DrompaOpt> v): opts("Options"), desc(d), requiredstr(r), vopts(v), func(_func), name(n) {
    opts.add(v);
  };
  void print() const {
    std::cout << std::setw(8) << " " << std::left << std::setw(12) << name
	      << std::left << std::setw(40) << desc << std::endl;
  }
  void printhelp() const {
    std::cout << boost::format("%1%:  %2%\n") % name % desc;
    std::cout << boost::format("Usage: drompa %1% [options] -o <output> --gt <genometable> %2%\n\n") % name % requiredstr;
    std::cout << opts.opts << std::endl;
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
      p.gt = read_genometable(values["gt"].as<std::string>());

      func(values, p);

    } catch (std::exception &e) {
      std::cout << e.what() << std::endl;
    }
  }
};

void opt::add(std::vector<DrompaOpt> st)
{
  boost::program_options::options_description o("Required",100);
  o.add_options()
    ("output,o",  boost::program_options::value<std::string>(),	 "Output prefix")
    ("gt",        boost::program_options::value<std::string>(),	 "Genome table")
    ;
  opts.add(o);

  for(auto x: st) {
    switch(x) {
    case DrompaOpt::CHIP:
      {
	boost::program_options::options_description o("Input",100);
	o.add_options()
	  ("input,i",   boost::program_options::value<std::vector<std::string>>(), "Specify ChIP data, Input data and name of ChIP sample\n     (separated by ',', values except for 1 can be omitted)\n     1:ChIP   2:Input   3:name   4:peaklist   5:binsize\n     6:scale_tag   7:scale_ratio   8:scale_pvalue\n")
	  //	  ("binsize,b", boost::program_options::value<int>()->default_value(binsize), "Bin size")
	  ;
	opts.add(o);
	break;
      }
    case DrompaOpt::NORM:
      {
	boost::program_options::options_description o("",100);
	o.add_options()
	  ("norm",      boost::program_options::value<int>()->default_value(1),	     "Normalization between ChIP and Input\n      0: not normalize\n      1: with total read number\n      2: with NCIS method\n")
	  ("sm",        boost::program_options::value<int>()->default_value(0),      "Smoothing width") // gausian ??
	  ;
	opts.add(o);
	break;
      }
    case DrompaOpt::THRE: 
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
    case DrompaOpt::ANNO_PC:
      {
	boost::program_options::options_description o("Annotation",100);
	o.add_options()
	  ("gene,g", boost::program_options::value<std::string>(),	  "Gene annotation file")
	  ("gftype", boost::program_options::value<int>()->default_value(1), "Format of gene annotation\n     0: RefFlat (default)\n     1: Ensembl\n     2: GTF (for S. pombe)\n     3: SGD (for S. cerevisiae)\n")
	  ("ars",    boost::program_options::value<std::string>(),	  "ARS list (for yeast)")
	  ("ter",    boost::program_options::value<std::string>(),	  "TER list (for S.cerevisiae)")  
	  ("bed",    boost::program_options::value<std::vector<std::string>>(), "<bedfile>,<label>: Specify bed file and name (<label> can be omited)")
	  ("repeat", boost::program_options::value<std::string>(),	  "Display repeat annotation (RepeatMasker format)") 
	  ;
	opts.add(o);
	break;
      }
    case DrompaOpt::ANNO_GV:
      {
	boost::program_options::options_description o("Optional data",100);
	o.add_options()
	  ("mp",     boost::program_options::value<std::string>(),  	  "Mappability file")
	  ("mpthre", boost::program_options::value<double>()->default_value(0.3), "Low mappability threshold")
	  ("gap",    boost::program_options::value<std::string>(),	  "Specify gapped regions to be shaded")
	  ("inter",  boost::program_options::value<std::vector<std::string>>(), "<interaction file>,<label>: Specify interaction file and name (<label> can be omited)")  // FDRde iro kaeru
	  ("gc",     boost::program_options::value<std::string>(), 	  "Visualize GC contents graph")
	  ("gcsize", boost::program_options::value<int>()->default_value(100000), "Window size for GC contents")
	  ("gd",     boost::program_options::value<std::string>(), 	  "Visualize gene density (number of genes for each window)")
	  ("gdsize", boost::program_options::value<int>()->default_value(100000), "Window size for gene density")
	  ;
	opts.add(o);
	break;
      }
    case DrompaOpt::DRAW:
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
    case DrompaOpt::REGION:
      {
	boost::program_options::options_description o("Region to draw",100);
	o.add_options()
	  ("chr",         boost::program_options::value<int>(),     "Output the specified chromosome only")
	  ("region,r",    boost::program_options::value<std::string>(),  "Specify genomic regions for drawing")
	  ("genefile",    boost::program_options::value<std::string>(),  "Specify gene loci to visualize")  
	  ("len_genefile",boost::program_options::value<int>()->default_value(50000), "extended length for each gene locus")
	  ;
	opts.add(o);
	break;
      }
    
    case DrompaOpt::SCALE:
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
    case DrompaOpt::OVERLAY:
      {
	boost::program_options::options_description o("For overlay",100);
	o.add_options()
	  ("ioverlay",  boost::program_options::value<std::vector<std::string>>(),	  "Input file")
	  //	  ("binsize2",  boost::program_options::value<int>()->default_value(binsize), "Bin size")
	  ("scale_tag2",   boost::program_options::value<double>(), "Scale for read line")
	  ("scale_ratio2", boost::program_options::value<double>(), "Scale for fold enrichment")
	  ("scale_pvalue2",boost::program_options::value<double>(), "Scale for -log10(p)")
	  ;
	opts.add(o);
	break;
      }
    case DrompaOpt::CG: 
      {
	boost::program_options::options_description o("CG",100);
	o.add_options()
	  ("cgthre",    boost::program_options::value<double>(), "Minimum threshold per kbp")
	  ;
	opts.add(o);
	break;
      }
    case DrompaOpt::TR: 
      {
	boost::program_options::options_description o("TR",100);
	o.add_options()
	  ("tssthre",    boost::program_options::value<double>(), "")
	  ;
	opts.add(o);
	break;
      }
    case DrompaOpt::PD:
      {
	boost::program_options::options_description o("PD",100);
	o.add_options()
	  ("pd",   boost::program_options::value<std::vector<std::string>>(), "Peak density file and name\n(separated by ',' <name> can be omited)")
	  ("prop",   boost::program_options::value<double>(),  "scale_tag")
	  ("pdsize", boost::program_options::value<int>()->default_value(100000), "windowsize for peak density")
	  ;
	opts.add(o);
	break;
      }
    case DrompaOpt::PROF:
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
      
    case DrompaOpt::OTHER:
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
    case DrompaOpt::CHIP:
      {
	for (auto x: {"input"}) if (!values.count(x)) PRINTERR("specify --" << x << " option.");

	//	chkminus<int>(values, "binsize", 0);

	std::vector<std::string> v(values["input"].as<std::vector<std::string>>());
	for(auto x:v) scan_samplestr(x, p.sample, p.samplepair);
	
	break;
      }
    case DrompaOpt::NORM:
      {
	chkminus<int>(values, "sm", 0);
	chkrange<int>(values, "norm", 0, 1);
	break;
      }
    case DrompaOpt::THRE: 
      {
	for (auto x: {"pthre_internal", "pthre_enrich", "qthre", "ipm", "ethre"}) chkminus<int>(values, x, -1);
	chkminus<int>(values, "width4lmd", 0);
	break;
      }
    case DrompaOpt::ANNO_PC:
      {
	chkrange<int>(values, "gftype", 0, 3);
	for (auto x: {"gene", "ars", "ter"}) if (values.count(x)) isFile(values[x].as<std::string>());
	break;
      }
    case DrompaOpt::ANNO_GV:
      {
	chkminus<int>(values, "mpthre", -1);
	for (auto x: {"gcsize", "gdsize"}) chkminus<int>(values, x, 0);
	break;
      }
    case DrompaOpt::DRAW:
      {
	for (auto x: {"ls", "lpp"}) chkminus<int>(values, x, 0);
	break;
      }
    case DrompaOpt::REGION:
      {
	for (auto x: {"region", "genefile"}) if (values.count(x)) isFile(values[x].as<std::string>());
	chkminus<int>(values, "len_genefile", -1);
	break;
      }
    case DrompaOpt::SCALE:
      {
	for (auto x: {"scale_tag","scale_ratio","scale_pvalue","bn","ystep"}) chkminus<int>(values, x, 0);
	break;
      }
    case DrompaOpt::OVERLAY:
      {
	for (auto x: {"scale_tag2","scale_ratio2","scale_pvalue2"}) chkminus<int>(values, x, 0);
	break;
      }
    case DrompaOpt::CG: 
      {
	for (auto x: {"cgthre"}) chkminus<int>(values, x, -1);
	break;
      }
    case DrompaOpt::TR: 
      {
	for (auto x: {"tssthre"}) chkminus<int>(values, x, -1);
	break;
      }
    case DrompaOpt::PD:
      {
	for (auto x: {"pd"}) if (!values.count(x)) PRINTERR("specify --" << x << " option.");
	for (auto x: {"pdsize"}) chkminus<int>(values, x, 0);
	
	std::vector<std::string> v(values["pd"].as<std::vector<std::string>>());
	for(auto &x: v) p.pd.push_back(scan_pdstr(x));
	break;
      }
    case DrompaOpt::PROF:
      {
	chkrange<int>(values, "ptype", 0, 4);
	chkrange<int>(values, "stype", 0, 1);
	chkrange<int>(values, "ntype", 0, 1);
	for (auto x: {"cw", "maxval", "hmsort"}) chkminus<int>(values, x, 0);
	break;
      }
      
    case DrompaOpt::OTHER:
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
  std::vector<std::string> str_bool = {"ON", "OFF"};
  std::vector<std::string> str_gftype = {"refFlat", "Ensembl", "gtf", "SGD"};
  //  std::vector<std::string> str_wigfiletype = {"BINARY", "COMPRESSED WIG", "WIG", "BEDGRAPH", "BIGWIG"};
  std::vector<std::string> str_norm  = { "OFF", "TOTALREAD", "NCIS" };
  std::vector<std::string> str_stype = { "ChIP read", "Enrichment ratio", "Enrichment P-value" };
  std::vector<std::string> str_ptype = { "NONE", "TSS", "TTS", "GENE100", "SPECIFIEDSITES" };
  std::vector<std::string> str_ntype = { "WHOLE GENOME", "TARGET REGIONS ONLY" };

  std::cout << boost::format("\n======================================\n");
  std::cout << boost::format("drompa version %1%: %2%\n\n") % VERSION % name;
  std::cout << boost::format("output prefix: %1%\n")     % values["output"].as<std::string>();
  std::cout << boost::format("Genome-table file: %1%\n") % values["gt"].as<std::string>();

  for(auto x: vopts) {
    switch(x) {
    case DrompaOpt::CHIP:
      {
	std::cout << boost::format("\nSamples\n");
	for(uint i=0; i<p.samplepair.size(); ++i) {
	  std::cout << (i+1) << ": ";
	  p.samplepair[i].print();
	}
	//	std::cout << boost::format("   Input format: %1%\n")    % str_wigfiletype[values["if"].as<int>()];
	break;
      }
    case DrompaOpt::NORM:
      {
	std::cout << boost::format("   ChIP/Input normalization: %s\n") % str_norm[values["norm"].as<int>()];
	if(values["sm"].as<int>()) std::cout << boost::format("   smoothing width: %1% bp\n") % values["sm"].as<int>();
	break;
      }
    case DrompaOpt::THRE: 
      {
	std::cout << boost::format("   Peak intensity threshold: %1$.2f\n")               % values["IPmaxthre"].as<double>();
	std::cout << boost::format("   Enrichment threshold: %1$.2f\n")                   % values["enrichthre"].as<double>();
	std::cout << boost::format("   p-value threshold (internal, -log10): %1$.2e\n")   % values["pthre_internal"].as<double>();
	std::cout << boost::format("   p-value threshold (enrichment, -log10): %1$.2e\n") % values["pthre_enrich"].as<double>();
	std::cout << boost::format("   FDR threshold: %1$.2e\n")                          % values["FDR"].as<double>();
  	break;
      }
    case DrompaOpt::ANNO_PC:
      {
	std::cout << boost::format("\nAnnotations:\n");
	if(values.count("gene")) std::cout << boost::format("   Gene file: %1%, Format: %2%\n")
			     % values["gene"].as<std::string>() % str_gftype[values["gftype"].as<int>()];
	MyOpt::printOpt<std::string>(values, "ars",    "   ARS file");
	MyOpt::printOpt<std::string>(values, "ter",    "   TER file");
	MyOpt::printOpt<std::string>(values, "repeat", "   Repeat file");
	MyOpt::printOpt<std::string>(values, "gc", "   GCcontents file");
	MyOpt::printOpt<std::string>(values, "gd", "   Gene density file");
	/*	if(d->arsfile)     std::cout << boost::format("   ARS file: %1%\n")          % values["ars"].as<std::string>();
	if(d->terfile)     std::cout << boost::format("   TER file: %1%\n")          % values["ter"].as<std::string>();
	if(d->repeat.argv) std::cout << boost::format("   Repeat file: %1%\n")       % values["repeat"].as<std::string>();*/
	MyOpt::printOpt<std::string>(values, "region", "   Region file");
	MyOpt::printVOpt<std::string>(values, "bed", "   Bed file");
	//	if(name != "PROFILE" || name != "HEATMAP") std::cout << boost::format("   name: %1%\n") % d->bed[i]->name;
	break;
      }
    case DrompaOpt::ANNO_GV:
      {
	std::cout << boost::format("\nAnnotations:\n");
	if(values.count("gene")) std::cout << boost::format("   Gene file: %1%, Format: %2%\n")
			     % values["gene"].as<std::string>() % str_gftype[values["gftype"].as<int>()];
	MyOpt::printOpt<std::string>(values, "gc", "   GCcontents file");
	MyOpt::printOpt<std::string>(values, "gd", "   Gene density file");
	/*	
	if(d->GC.argv)     std::cout << boost::format("   GCcontents file: %1%\n")   % values["gc"].as<std::string>();
	if(d->GD.argv)     std::cout << boost::format("   Gene density file: %1%\n") % values["gd"].as<std::string>();*/
	MyOpt::printVOpt<std::string>(values, "inter", "   Interaction file");
	if (values.count("mp")) {
	  std::cout << boost::format("Mappability file directory: %1%\n") % values["mp"].as<std::string>();
	  std::cout << boost::format("\tLow mappablitiy threshold: %1%\n") % values["mpthre"].as<double>();
	}
	break;
      }
    case DrompaOpt::DRAW:
      {
	std::cout << boost::format("\nFigure parameter:\n");
	std::cout << boost::format("   Display read: ChIP %1%, Input %2%\n") % str_bool[values["showctag"].as<int>()] % str_bool[values["showitag"].as<int>()];
	std::cout << boost::format("   Display enrichment: %1%\n")           % str_bool[values["showratio"].as<int>()];
	std::cout << boost::format("   Display pvalue (internal): %1%\n")    % str_bool[values["showpinter"].as<int>()];
	std::cout << boost::format("   Display pvalue (ChIP/Input): %1%\n")  % str_bool[values["showpenrich"].as<int>()];
	std::cout << boost::format("   Background color: %1%\n")             % str_bool[!values["offbg"].as<int>()];
	std::cout << boost::format("   Y label: %1%\n")                      % str_bool[!values["offylab"].as<int>()];
	std::cout << boost::format("   Y memory: %1%\n")                     % str_bool[!values["offymem"].as<int>()];
	break;
      }
    
    case DrompaOpt::REGION:
      {
	break;
      }
    case DrompaOpt::SCALE:
      {
	for (auto x: {"scale_tag","scale_ratio","scale_pvalue","scale_tag2","scale_ratio2","scale_pvalue2","bn","ystep"}) chkminus<int>(values, x, 0);
	break;
      }
    case DrompaOpt::OVERLAY:
      {
	if(values.count("ioverlay")) { 
	  std::cout << boost::format("\nOverlayed samples\n");
	  for(uint i=0; i<p.samplepair.size(); ++i) {
	    std::cout << (i+1) << ": ";
	    p.samplepair[i].print();
	  }
	}
	break;
      }
    case DrompaOpt::CG: 
      {
	for (auto x: {"cgthre"}) chkminus<int>(values, x, -1);
	break;
      }
    case DrompaOpt::TR: 
      {
	for (auto x: {"tssthre"}) chkminus<int>(values, x, -1);
	break;
      }
    case DrompaOpt::PD:
      {
	std::cout << boost::format("\nSamples\n");
	for(uint i=0; i<p.pd.size(); ++i) {
	  std::cout << boost::format("   IP%1%: %2%\tname: %3%\n") % (i+1) % p.pd[i].argv % p.pd[i].name;
	}
	break;
      }
    case DrompaOpt::PROF:
      {
	std::cout << boost::format("   show type: %1$\n")             % str_stype[values["stype"].as<int>()];
	std::cout << boost::format("   profile type: %1$\n")          % str_ptype[values["ptype"].as<int>()];
	std::cout << boost::format("   profile normalization: %1$\n") % str_ntype[values["ntype"].as<int>()];
	break;
      }
      
    case DrompaOpt::OTHER:
      {
	for (auto x: {"threads"}) chkminus<int>(values, x, 0);
	break;
      }
    }
    return;
  }
  
  if(values.count("chr")) std::cout << boost::format("output %1% only.\n") % values["chr"].as<int>();

  printf("======================================\n");
  return;
}

#endif /* _DD_OPT_HPP_ */
