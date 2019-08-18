/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "dd_gv.hpp"
#include "dd_readfile.hpp"
#include "../submodules/SSP/common/BedFormat.hpp"

using namespace boost::program_options;
using namespace MyOpt;
using namespace DROMPA;

#define SETOPT_OVER(name,type,defval,val) name, value<type>()->default_value(defval)->NOTIFY_OVER(type, val, name)
#define SETOPT_RANGE(name,type,defval,min,max) name, value<type>()->default_value(defval)->NOTIFY_RANGE(type, min, max, name)
#define NOTIFY_OVER(type,val,name)      notifier(boost::bind(&MyOpt::over<type>, _1, val, name))
#define NOTIFY_RANGE(type,min,max,name) notifier(boost::bind(&MyOpt::range<type>, _1, min, max, name))

inline std::string addDoubleQuotes (const std::string &str) { return "\"" + str + "\"";}

void vSampleInfo::addSampleInfo(const std::string &str, const std::vector<chrsize> &gt, const WigType iftype) {
  int32_t binsize(0);
  std::vector<std::string> v;
  ParseLine(v, str, ',');

  if(v.size() >8) {
    std::cerr << "error: sample std::string has ',' more than 8: " << str << std::endl;
    exit(1);
  }
  if(v[0] == "") {
    std::cerr << "please specify ChIP sample: " << str << std::endl;
    exit(1);
  }
  isFile(v[0]);

  if(v.size() >4 && v[4] != "") {
    try { binsize = stoi(v[4]); }
    catch (...) { std::cerr << "Warning: invalid binsize " << v[4] << "." << std::endl; }
  }

  // ChIP sample
  if(!Exists(v[0])) vsinfo[v[0]] = SampleInfo(v[0], gt, binsize, iftype);
  if(vsinfo[v[0]].getbinsize() <= 0) PRINTERR("please specify binsize.\n");

  // Input sample
  if(v.size() >=2 && v[1] != "") {
    if(!Exists(v[1])) vsinfo[v[1]] = SampleInfo(v[1], gt, binsize, iftype);
    if(vsinfo[v[0]].getbinsize() != vsinfo[v[1]].getbinsize()) PRINTERR("binsize of ChIP and Input should be same. " << str);
  }
}

void Annotation::setOptsPC(MyOpt::Opts &allopts)
{
  MyOpt::Opts opt("Annotation",100);
  opt.add_options()
    ("gene,g", value<std::string>(), "Gene annotation file")
    (SETOPT_RANGE("gftype", int32_t, 0, 0, 2),
     "Format of gene annotation\n     0: RefFlat\n     1: GTF\n     2: SGD (for S. cerevisiae)\n")
    ("showtranscriptname", "Show transcript name (default: gene name)")
    ("ars",    value<std::string>(), "ARS list (for yeast)")
    ("ter",    value<std::string>(), "TER list (for S.cerevisiae)")
    ("showars", "Display ARS and TER and do not display genes")
    ("bed",    value<std::vector<std::string>>(), "<bedfile>,<label>: BED file (<label> can be omited)")
    ("repeat", value<std::string>(), "Display repeat annotation (RepeatMasker format)")
    ;
  allopts.add(opt);
}

void Annotation::setOptsGV(MyOpt::Opts &allopts) {
  MyOpt::Opts opt("Genome view",100);
  opt.add_options()
    ("inter",  value<std::vector<std::string>>(), "<interaction file>,<label>: Interaction file and label (<label> can be omited)")
    ("chiadrop", value<std::string>(), "ChIADrop file (single-cell ChIA-PET)")
    ("mp",     value<std::string>(), "Mappability file")
    ("mpthre", value<double>()->default_value(0.3), "Low mappability threshold")
    ("gap",    value<std::string>(), "Specify gapped regions to be shaded")
    ("GC",     value<std::string>(), "Visualize GC contents graph")
    (SETOPT_OVER("gcsize", int32_t, 100000, 1), "Window size for GC contents")
    ("GD",     value<std::string>(), "Visualize gene density (number of genes for each window)")
    (SETOPT_OVER("gdsize", int32_t, 500000, 1), "Window size for gene density")
    ;
  allopts.add(opt);
}

void Annotation::setValuesPC(const Variables &values) {
  DEBUGprint("AnnoPC setValues...");

  try {
    for (auto x: {"gene", "ars", "ter", "repeat"}) if (values.count(x)) isFile(getVal<std::string>(values, x));
    if (values.count("gene")) {
      genefile = getVal<std::string>(values, "gene");
      gftype   = getVal<int32_t>(values, "gftype");
      gmp = getGMP();
    }
    if (values.count("ars")) {
      arsfile = getVal<std::string>(values, "ars");
      parseARSOriDB(arsfile, gmp);
    }
    showars = values.count("showars");
    if (values.count("ter")) {
      terfile = getVal<std::string>(values, "ter");
      parseTER(terfile, gmp);
    }
    showtranscriptname = values.count("showtranscriptname");
    //	printMap(gmp);
    if (values.count("repeat")) repeatfile = getVal<std::string>(values, "repeat");

    if (values.count("bed")) {
      for(auto &x: getVal<std::vector<std::string>>(values, "bed")) {
	std::vector<std::string> v;
	ParseLine(v, x, ',');
	readBedFile<bed12>(v);
      }
    }
  } catch (const boost::bad_any_cast& e) {
    std::cout << e.what() << std::endl;
    exit(0);
  }

  DEBUGprint("AnnoPC setValues done.");
}

void Annotation::setValuesGV(const Variables &values) {
  DEBUGprint("AnnoGV setValues...");

  try {
    if (values.count("inter")) {
      for(auto &x: getVal<std::vector<std::string>>(values, "inter")) {
	//    std::cout << x << std::endl;
	std::vector<std::string> v;
	boost::split(v, x, boost::algorithm::is_any_of(","));
	std::string label(v[0]);
	std::string tool("mango");
	if(v.size() >= 2) label = v[1];
	if(v.size() >= 3) tool  = v[2];
	vinterlist.emplace_back(InteractionSet(v[0], label, tool));
      }
    }
    if (values.count("chiadrop")) parse_ChIADropData(getVal<std::string>(values, "chiadrop"));
    if (values.count("mp")) mpfile = getVal<std::string>(values, "mp");
    mpthre = getVal<double>(values, "mpthre");
    if (values.count("gap")) gapfile = getVal<std::string>(values, "gap");
    if (values.count("GC")) GC.setValue(getVal<std::string>(values, "GC"), getVal<int32_t>(values, "gcsize"));
    if (values.count("GD")) GD.setValue(getVal<std::string>(values, "GD"), getVal<int32_t>(values, "gdsize"));

  } catch (const boost::bad_any_cast& e) {
    std::cout << e.what() << std::endl;
    exit(0);
  }
  DEBUGprint("AnnoGV setValues done.");
}

void Annotation::InitDumpPC(const Variables &values) const {
  std::vector<std::string> str_gftype = {"refFlat", "GTF", "SGD"};
  std::vector<std::string> str_gtype = {"Transcript", "Gene"};
  DEBUGprint("INITDUMP:DrompaCommand::ANNO_PC");
  std::cout << boost::format("\nAnnotations:\n");
  if(genefile != "") std::cout << boost::format("   Gene file: %1%, Format: %2%\n") % genefile % str_gftype[gftype];
  //      if(arsfile != "")  std::cout << boost::format("   ARS file: %1%\n") % arsfile;
  //if(terfile != "")  std::cout << boost::format("   TER file: %1%\n") % terfile;
  //if(repeatfile != "") std::cout << boost::format("   Repeat file: %1%\n") % repeatfile;
  printOpt<std::string>(values, "ars",    "   ARS file");
  printOpt<std::string>(values, "ter",    "   TER file");
  if(showars) std::cout << "Display ARS and TER only." << std::endl;
  printOpt<std::string>(values, "repeat", "   Repeat file");
  printOpt<std::string>(values, "region", "   Region file");
  printVOpt<std::string>(values, "bed", "   Bed file");
}

void Annotation::InitDumpGV(const Variables &values) const {
  DEBUGprint("INITDUMP:DrompaCommand::ANNO_GV");
  std::cout << boost::format("\nAnnotations:\n");
  printVOpt<std::string>(values, "inter", "   Interaction file");
  if (mpfile != "") {
    std::cout << boost::format("Mappability file directory: %1%\n") % mpfile;
    std::cout << boost::format("\tLow mappablitiy threshold: %1$2f\n") % mpthre;
  }
  printOpt<std::string>(values, "gc", "   GCcontents file");
  printOpt<std::string>(values, "gd", "   Gene density file");
  /*	if(d->GC.argv)     std::cout << boost::format("   GCcontents file: %1%\n")   % values["gc"].as<std::string>();
	if(d->GD.argv)     std::cout << boost::format("   Gene density file: %1%\n") % values["gd"].as<std::string>();*/
}


void Profile::setOpts(MyOpt::Opts &allopts) {
  MyOpt::Opts opt("PROFILE AND HEATMAP",100);
  opt.add_options()
    (SETOPT_RANGE("ptype", int32_t, TSS, TSS, PTYPENUM-1),
     "Region: 0; around TSS, 1; around TES, 2; divide gene into 100 subregions 3; around peak sites")
    (SETOPT_RANGE("stype", int32_t, 0, 0, 1),
     "Distribution: 0; ChIP read, 1; ChIP/Input enrichment")
    (SETOPT_RANGE("ntype", int32_t, 0, 0, 1),
     "Normalization: 0; total read 1; target regions only")
    (SETOPT_OVER("widthfromcenter", int32_t, 2500, 1), "width from the center")
    (SETOPT_OVER("maxval", double, 0, 0), "Upper limit for heatmap")
    (SETOPT_OVER("hmsort", int32_t, 1, 1),  "Column number for sorting sites")
    ;
  allopts.add(opt);
}

void Profile::setValues(const Variables &values) {
  DEBUGprint("PROF setValues...");

  try {
    ptype = getVal<int32_t>(values, "ptype");
    stype = getVal<int32_t>(values, "stype");
    ntype = getVal<int32_t>(values, "ntype");
    width_from_center = getVal<int32_t>(values, "widthfromcenter");
    maxval = getVal<double>(values, "maxval");
    hmsort = getVal<int32_t>(values, "hmsort");
  } catch (const boost::bad_any_cast& e) {
    std::cout << e.what() << std::endl;
    exit(0);
  }
  DEBUGprint("PROF setValues done.");
}

void Profile::InitDump() const {
  DEBUGprint("INITDUMP:DrompaCommand::PROF");

  std::vector<std::string> str_ptype = { "TSS", "TTS", "GENE100", "BEDSITES" };
  std::vector<std::string> str_stype = { "ChIP read", "Enrichment ratio", "Enrichment P-value" };
  std::vector<std::string> str_ntype = { "WHOLE GENOME", "TARGET REGIONS ONLY" };

  std::cout << boost::format("   Profile type: %1%\n")  % str_ptype[ptype];
  std::cout << boost::format("   Show type: %1%\n")     % str_stype[stype];
  std::cout << boost::format("   Normalization: %1%\n") % str_ntype[ntype];
  std::cout << boost::format("   Width from center: %1%\n") % width_from_center;
}


void Threshold::setOpts(MyOpt::Opts &allopts) {
  sigtest = false;

  MyOpt::Opts opt("Threshold",100);
  opt.add_options()
    ("pthre_internal", value<double>()->default_value(1e-3), "p-value for ChIP internal")
    ("pthre_enrich",   value<double>()->default_value(1e-3), "p-value for ChIP/Input enrichment")
    ("qthre",          value<double>()->default_value(1),    "FDR")
    ("ethre",          value<double>()->default_value(2),    "IP/Input fold enrichment")
    ("ipm",            value<double>()->default_value(0),    "Read intensity of peak summit")
    ("callpeak",       "Call peaks by DROMPA")
    //    (SETOPT_OVER("width4lmd", int32_t, 100000, 0), "Width for calculating local lambda")
    ;
  allopts.add(opt);
}

void Threshold::setValues(const Variables &values) {
  DEBUGprint("Threshold setValues...");
  try {
    pthre_inter  = getVal<double>(values, "pthre_internal");
    pthre_enrich = getVal<double>(values, "pthre_enrich");
    qthre        = getVal<double>(values, "qthre");
    ethre        = getVal<double>(values, "ethre");
    ipm          = getVal<double>(values, "ipm");
    if(values.count("callpeak")) sigtest = true;
    //    width4lmd    = getVal<int32_t>(values, "width4lmd");
  } catch (const boost::bad_any_cast& e) {
    std::cout << e.what() << std::endl;
    exit(0);
  }
  DEBUGprint("Threshold setValues done.");
}

void Threshold::InitDump() const {
  std::vector<std::string> str_bool = {"OFF", "ON"};
  DEBUGprint("INITDUMP:DrompaCommand::THRE");

  //  std::cout << boost::format("\nPeak calling by DROMPA+: %1%\n") % str_bool[sigtest];
  if (sigtest) {
    std::cout << boost::format("   p-value threshold (internal, -log10): %1$.2e\n")   % pthre_inter;
    std::cout << boost::format("   p-value threshold (internal, -log10): %1$.2e\n")   % pthre_inter;
    std::cout << boost::format("   p-value threshold (enrichment, -log10): %1$.2e\n") % pthre_enrich;
    std::cout << boost::format("   FDR threshold: %1$.2e\n")                          % qthre;
    std::cout << boost::format("   Peak intensity threshold: %1$.2f\n")               % ipm;
    std::cout << boost::format("   Enrichment threshold: %1$.2f\n")                   % ethre;
  }
}

void DrawRegion::setOpts(MyOpt::Opts &allopts) {
  MyOpt::Opts opt("Region to draw",100);
  opt.add_options()
    ("chr", value<std::string>(),
     "Output the specified chromosome only")
    ("region,r", value<std::string>(),
     "Specify genomic regions for drawing")
    ("genelocifile", value<std::string>(),
     "Specify of a file gene namaes to visualize")
    (SETOPT_OVER("len_geneloci", int32_t, 50000, 1), "extended length for each gene locus")
    ;
  allopts.add(opt);
}

void DrawRegion::setValues(const Variables &values) {
  DEBUGprint("DrawRegion setValues...");
  try {
    for (auto x: {"region", "genelocifile"}) if (values.count(x)) isFile(getVal<std::string>(values, x));
    if (values.count("chr")) chr = getVal<std::string>(values, "chr");
    if (values.count("region")) {
      isRegion = true;
      regionBed = parseBed<bed>(getVal<std::string>(values, "region"));
      if(!regionBed.size()) PRINTERR("Error no bed regions in " << getVal<std::string>(values, "region"));
      //      printBed(regionBed);
    }
    if (values.count("genelocifile")) {
      getGeneLoci(getVal<std::string>(values, "genelocifile"));
      if (!values.count("gene")) PRINTERR("Please specify --gene option when supplying --genelocifile.");
    }
    len_geneloci = getVal<int32_t>(values, "len_geneloci");

  } catch (const boost::bad_any_cast& e) {
    std::cout << e.what() << std::endl;
    exit(0);
  }
  DEBUGprint("DrawRegion setValues done.");
}

void DrawRegion::InitDump(const Variables &values) const {
  std::vector<std::string> str_bool = {"OFF", "ON"};

  DEBUGprint("INITDUMP:DrompaCommand::DRAWREGION");
  printOpt<std::string>(values, "region",    "   Region file");
  if (chr != "") std::cout << boost::format("   output chr%1% only.\n") % chr;
  if (values.count("genelocifile")) {
    std::cout << boost::format("   Geneloci file: %1%, around %2% bp\n") % getVal<std::string>(values, "genelocifile") % len_geneloci;
  }
}

void Global::setOpts(const std::vector<DrompaCommand> &st, const CommandParamSet &cps) {
  MyOpt::Opts o("Required",100);
  o.add_options()
    ("output,o",  value<std::string>(), "Output prefix")
    ("gt",        value<std::string>(), "Genome table")
    ;
  opts.add(o);
  for(auto &x: st) {
    switch(x) {
    case DrompaCommand::CHIP:
      {
	options_description o("Input",100);
	o.add_options()
	  ("input,i", value<std::vector<std::string>>(),
	   "Specify ChIP-Input pair and label\n     (separated by ',', values except for 1 can be omitted)\n     1:ChIP   2:Input   3:label   4:peaklist   5:binsize\n     6:scale_tag   7:scale_ratio   8:scale_pvalue")
	  ("ioverlay", value<std::vector<std::string>>(), "Specify sample pairs to overlay (same manner as -i)")
	  (SETOPT_RANGE("if", int32_t, static_cast<int32_t>(WigType::NONE), 0, static_cast<int32_t>(WigType::WIGTYPENUM) -1),
	   "Input file format\n   0: compressed wig (.wig.gz)\n   1: uncompressed wig (.wig)\n   2: bedGraph (.bedGraph)\n   3: bigWig (.bw)")
	  ;
	opts.add(o);
	break;
      }
    case DrompaCommand::NORM:    setOptsNorm(opts, cps.sm); break;
    case DrompaCommand::THRE:    thre.setOpts(opts); break;
    case DrompaCommand::ANNO_PC: anno.setOptsPC(opts); break;
    case DrompaCommand::ANNO_GV: anno.setOptsGV(opts); break;
    case DrompaCommand::DRAW:    drawparam.setOpts(opts, cps); break;
    case DrompaCommand::REGION:  drawregion.setOpts(opts); break;
    case DrompaCommand::PROF:    prof.setOpts(opts); break;
    case DrompaCommand::CG:
      {
	options_description o("CG",100);
	o.add_options()
	  ("cgthre",    value<double>(), "Minimum threshold per kbp")
	  ;
	opts.add(o);
	break;
      }
    case DrompaCommand::TR:
      {
	options_description o("TR",100);
	o.add_options()
	  ("tssthre",    value<double>(), "")
	  ;
	opts.add(o);
	break;
      }
    case DrompaCommand::PD:
      {
	options_description o("PD",100);
	o.add_options()
	  ("pd",   value<std::vector<std::string>>(), "Peak density file and name\n(separated by ',' <name> can be omited)")
	  ("prop",   value<double>(),  "scale_tag")
	  (SETOPT_OVER("pdsize", int32_t, 100000, 1), "windowsize for peak density")
	  ;
	opts.add(o);
	break;
      }
    case DrompaCommand::OTHER: setOptsOther(opts); break;
    }
  }
}

void Global::setValues(const std::vector<DrompaCommand> &vopts, const Variables &values) {
  for (auto x: {"output", "gt"}) if (!values.count(x)) PRINTERR("specify --" << x << " option.");

  oprefix = getVal<std::string>(values, "output");
  gt = read_genometable(getVal<std::string>(values, "gt"));

  for (auto op: vopts) {
    switch(op) {
    case DrompaCommand::CHIP:
      {
	DEBUGprint("ChIP setValues...");
	try {
//	  if (!values.count("input")) PRINTERR("specify --input option.");

	  // SamplePairOverlayed first
	  if (values.count("input")) {
	    std::vector<std::string> v(getVal<std::vector<std::string>>(values, "input"));
	    for (auto &x:v) {
	      vsinfo.addSampleInfo(x, gt, iftype);
	      samplepair.emplace_back(x, vsinfo);
	    }
	  }

	  // SamplePairOverlayed second
	  if (values.count("ioverlay")) {
	    std::vector<std::string> v(getVal<std::vector<std::string>>(values, "ioverlay"));
	    int32_t num(0);
	    for (auto &x:v) {
	      vsinfo.addSampleInfo(x, gt, iftype);
	      samplepair[num].setSecondSample(x, vsinfo);
	      ++num;
	    }
	  }

	  if (values.count("if")) iftype = static_cast<WigType>(getVal<int32_t>(values, "if"));

	} catch(const boost::bad_any_cast& e) {
	  std::cout << e.what() << std::endl;
	  exit(0);
	}
	DEBUGprint("ChIP setValues done.");

	break;
      }
    case DrompaCommand::NORM: setValuesNorm(values); break;
    case DrompaCommand::THRE: thre.setValues(values); break;
    case DrompaCommand::ANNO_PC: anno.setValuesPC(values); break;
    case DrompaCommand::ANNO_GV: anno.setValuesGV(values); break;
    case DrompaCommand::PROF: prof.setValues(values); break;
    case DrompaCommand::DRAW: drawparam.setValues(values, samplepair.size()); break;
    case DrompaCommand::REGION: drawregion.setValues(values); break;
    case DrompaCommand::CG:
      {
	DEBUGprint("Global::setValues::CG");
	//	for (auto x: {"cgthre"}) chkminus<int>(values, x, -1);
	break;
      }
    case DrompaCommand::TR:
      {
	DEBUGprint("Global::setValues::TR");
	//	for (auto x: {"tssthre"}) chkminus<int>(values, x, -1);
	break;
      }
    case DrompaCommand::PD:
      {
	DEBUGprint("Global::setValues::PD");
	for (auto x: {"pd"}) if (!values.count(x)) PRINTERR("specify --" << x << " option.");
	//	for (auto x: {"pdsize"}) chkminus<int>(values, x, 0);

	std::vector<std::string> v(getVal<std::vector<std::string>>(values, "pd"));
	for(auto &x: v) pd.emplace_back(scan_pdstr(x));
	break;
      }
    case DrompaCommand::OTHER: setValuesOther(values); break;
    }
  }

  return;
}

void Global::setOptsNorm(MyOpt::Opts &allopts, const int32_t defsm) {
  MyOpt::Opts o("Normalization",100);
  o.add_options()
    (SETOPT_RANGE("norm", int32_t, 1, 0, 2),
     "Normalization between ChIP and Input\n      0: not normalize\n      1: with total read number (genome)\n      2: with total read number (each chr)\n      3: with NCIS method\n")
    (SETOPT_OVER("sm", int32_t, defsm, 0), "# of bins for Gausian smoothing")
    ;
  allopts.add(o);
}

void Global::setValuesNorm(const Variables &values) {
  DEBUGprint("Norm setValues...");
  try {
    norm = getVal<int32_t>(values, "norm");
    smoothing = getVal<int32_t>(values, "sm");
  } catch (const boost::bad_any_cast& e) {
    std::cout << e.what() << std::endl;
    exit(0);
  }
  DEBUGprint("Norm setValues done.");
}

void Global::setOptsOther(MyOpt::Opts &allopts) {
  MyOpt::Opts o("Others",100);
  o.add_options()
    ("includeYM", "output peaks of chromosome Y and M")
    ("showchr",   "Output chromosome-separated pdf files")
    ("png",     "Output with png format (Note: output each page separately)")
    (SETOPT_OVER("threads,p", int32_t, 1, 1), "number of threads to launch")
    ("help,h", "show help message")
    ;
  allopts.add(o);
}

void Global::setValuesOther(const Variables &values) {
  DEBUGprint("Other setValues...");
  try {
    includeYM = values.count("includeYM");
    ispng = values.count("png");
    showchr = values.count("showchr");
  } catch(const boost::bad_any_cast& e) {
    std::cout << e.what() << std::endl;
    exit(0);
  }
  DEBUGprint("Other setValues done.");
}

void Global::InitDumpChIP() const {
  std::vector<std::string> str_wigfiletype = {"BINARY", "COMPRESSED WIG", "WIG", "BEDGRAPH", "BIGWIG"};
  DEBUGprint("INITDUMP:DrompaCommand::CHIP");
  printf("\nSamples\n");
  for (size_t i=0; i<samplepair.size(); ++i) {
    std::cout << (i+1) << ": ";
    samplepair[i].print();
  }
  if (iftype < WigType::NONE) std::cout << boost::format("Input format: %1%\n") % str_wigfiletype[static_cast<int32_t>(iftype)];
}

void Global::InitDumpNorm() const {
  std::vector<std::string> str_norm = { "OFF", "TOTALREAD GENOME", "TOTALREAD CHR", "NCIS" };

  DEBUGprint("INITDUMP:DrompaCommand::NORM");
  std::cout << boost::format("   ChIP/Input normalization: %1%\n") % str_norm[norm];
  if (smoothing) std::cout << boost::format("   smoothing width: %1% bins\n") % smoothing;
}

void Global::InitDumpOther() const {
  std::vector<std::string> str_bool = {"OFF", "ON"};
  std::vector<std::string> str_format = {"PDF", "PNG"};

  DEBUGprint("INITDUMP:DrompaCommand::OTHER");
  std::cout << boost::format("   Output format: %1%\n") % str_format[ispng];
  if (includeYM) std::cout << boost::format("   include chromosome Y and M\n");
  //  if (!showchr) std::cout << boost::format("   remove chr pdfs\n");
}


void DrawParam::setOpts(MyOpt::Opts &allopts, const CommandParamSet &cps) {
  MyOpt::Opts opt("Drawing",100);
  opt.add_options()
    (SETOPT_RANGE("showctag",   int32_t, cps.showctag,  0, 1), "Display ChIP read lines")
    (SETOPT_RANGE("showitag",   int32_t, 0,             0, 2), "Display Input read lines (0:off 1:all 2:first one)")
    (SETOPT_RANGE("showratio",  int32_t, cps.showratio, 0, 2), "Display ChIP/Input ratio (0:off 1:liner scale 2:logscale)")
    (SETOPT_RANGE("showpinter", int32_t, 0,             0, 1), "Display -log10(p) lines for ChIP internal")
    (SETOPT_RANGE("showpenrich",int32_t, 0,             0, 1), "Display -log10(p) lines for ChIP/Input enrichment")
    (SETOPT_OVER("scale_tag",    double, cps.scaletag,     0), "Scale for read line")
    (SETOPT_OVER("scale_ratio",  double, cps.scaleratio,   0), "Scale for fold enrichment")
    (SETOPT_OVER("scale_pvalue", double, cps.scalepvalue,  0), "Scale for -log10(p)")
    (SETOPT_OVER("ls",  int32_t, 500, 1), "Width for each line (kp)")
    (SETOPT_OVER("lpp", int32_t, 1, 1), "Line number per page")
    (SETOPT_OVER("bn",  int32_t, 2, 1), "Number of memories of y-axis")
    (SETOPT_OVER("ystep", double, 15, 1), "Height of read line")
    (SETOPT_OVER("width_page", int32_t, 1088, 1), "Width(pixel) of pdf page")
    (SETOPT_OVER("width_draw", int32_t, 750, 1), "Width(pixel) of read line")
    (SETOPT_RANGE("alpha", double, 1,  0, 1), "Transparency of read distribution")
    ("offymem", "Omit Y memory")
    ("offylabel", "Omit Y label")
    ;
  allopts.add(opt);
}

void DrawParam::setValues(const Variables &values, const int32_t n) {
  DEBUGprint("DrawParam setValues...");
  try {
    showctag = getVal<int32_t>(values, "showctag");
    showitag = getVal<int32_t>(values, "showitag");
    showratio = getVal<int32_t>(values, "showratio");
    showpinter = getVal<int32_t>(values, "showpinter");
    showpenrich = getVal<int32_t>(values, "showpenrich");
    width_page_pixel = getVal<int32_t>(values, "width_page");
    width_draw_pixel = getVal<int32_t>(values, "width_draw");
    width_per_line = 1000 * getVal<int32_t>(values, "ls");
    linenum_per_page = getVal<int32_t>(values, "lpp");
    barnum = getVal<int32_t>(values, "bn");
    ystep  = getVal<double>(values, "ystep");
    showymem = !values.count("offymem");
    showylab = !values.count("offylabel");
    alpha = getVal<double>(values, "alpha");

    scale_tag    = getVal<double>(values, "scale_tag");
    scale_ratio  = getVal<double>(values, "scale_ratio");
    scale_pvalue = getVal<double>(values, "scale_pvalue");

    samplenum = n;
  } catch (const boost::bad_any_cast& e) {
    std::cout << e.what() << std::endl;
    exit(0);
  }
  DEBUGprint("DrawParam setValues done.");
}

void DrawParam::InitDump() const {
  std::vector<std::string> str_bool = {"OFF", "ON"};
  std::vector<std::string> str_input = {"OFF", "ALL", "FIRST"};
  std::vector<std::string> str_ratio = {"OFF", "Linear", "Logratio"};

  DEBUGprint("INITDUMP:DrompaCommand::DRAW");
  std::cout << boost::format("\nFigure parameter:\n");
  std::cout << boost::format("   Display read: ChIP %1%, Input %2%, y-axis scale: %3%\n") % str_bool[showctag] % str_input[showitag] % scale_tag;
  std::cout << boost::format("   Display enrichment: %1%, y-axis scale: %2%\n")           % str_ratio[showratio] % scale_ratio;
  std::cout << boost::format("   Display pvalue (internal): %1%, y-axis scale: %2%\n")    % str_bool[showpinter] % scale_pvalue;
  std::cout << boost::format("   Display pvalue (ChIP/Input): %1%, y-axis scale: %2%\n")  % str_bool[showpenrich] % scale_pvalue;
  std::cout << boost::format("   Width per line: %1% kbp\n")           % (width_per_line/1000);
  std::cout << boost::format("   Y-axis label: %1%\n")                 % str_bool[showylab];
  std::cout << boost::format("   Y-axis memory: %1%\n")                % str_bool[showymem];
}

Global::pdSample Global::scan_pdstr(const std::string &str){
  std::vector<std::string> v;
  ParseLine(v, str, ',');

  if(v.size() > 2) {
    std::cerr << "error: sample std::string has ',' more than 2: " << str << std::endl;
    exit(1);
  }

  pdSample pd;
  if(v[0] == "") {
    std::cerr << "please specify file: " << str << std::endl;
    exit(1);
  } else {
    pd.argv = v[0];
  }
  if(v[1] != "") pd.name = v[1];
  else pd.name = v[0];

  return pd;
}
