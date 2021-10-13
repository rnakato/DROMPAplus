/* Copyright(c) Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "dd_gv.hpp"
#include "dd_readfile.hpp"
#include "extendBedFormat.hpp"

using namespace boost::program_options;
using namespace MyOpt;
using namespace DROMPA;

#define SETOPT_OVER(name,type,defval,val) name, value<type>()->default_value(defval)->NOTIFY_OVER(type, val, name)
#define SETOPT_RANGE(name,type,defval,min,max) name, value<type>()->default_value(defval)->NOTIFY_RANGE(type, min, max, name)
#define NOTIFY_OVER(type,val,name)      notifier(std::bind(&MyOpt::over<type>, std::placeholders::_1, val, name))
#define NOTIFY_RANGE(type,min,max,name) notifier(std::bind(&MyOpt::range<type>, std::placeholders::_1, min, max, name))

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
    ("bed",    value<std::vector<std::string>>(), "<bedfile>,<label>: BED file (<label> can be omited)")
    ("bed12",  value<std::vector<std::string>>(), "<bed12file>,<label>: BED12 file (<label> can be omited)")
    ("repeat", value<std::string>(), "Display repeat annotation (RepeatMasker format)")
    ;
  allopts.add(opt);
}

void Annotation::setOptsGV(MyOpt::Opts &allopts)
{
  MyOpt::Opts opt("Genome view",100);
  opt.add_options()
    ("ideogram", value<std::string>(), "Cytoband file for drawing ideogram")
    ("inter",  value<std::vector<std::string>>(),
     "<interaction file>,<label>,<tool>: Interaction file\n     (label and tool can be omitted)\n     <tool>:mango (Mango format, default)\n            hiccups (HICCUPS format)\n            bedpe (BEDPE format)\n")
    ("chiadrop", value<std::string>(), "ChIADrop file (single-cell ChIA-PET)")
    ("chia_distance_thre", value<int32_t>()->default_value(100000), "(with --chiadrip) max distance of neighboring reads in each GEM")
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

void Annotation::setValuesPC(const Variables &values)
{
  DEBUGprint_FUNCStart();

  try {
    for (auto x: {"gene", "ars", "ter", "repeat"}) {
      if (values.count(x)) isFile(getVal<std::string>(values, x));
    }
    if (values.count("gene")) {
      genefile = getVal<std::string>(values, "gene");
      gftype   = getVal<int32_t>(values, "gftype");
      gmp = getGMP();
    }
    if (values.count("ars")) {
      arsfile = getVal<std::string>(values, "ars");
      parseARSOriDB(arsfile, arsgmp);
    }
    if (values.count("ter")) {
      terfile = getVal<std::string>(values, "ter");
      parseTER(terfile, arsgmp);
    }
    showtranscriptname = values.count("showtranscriptname");
    //	printMap(gmp);
    if (values.count("repeat")) repeatfile = getVal<std::string>(values, "repeat");

    if (values.count("bed")) {
      vbedlist = readBedFile<bed>(getVal<std::vector<std::string>>(values, "bed"));
    }
    if (values.count("bed12")) {
      vbed12list = readBedFile<bed12>(getVal<std::vector<std::string>>(values, "bed12"));
    }
  } catch (const boost::bad_any_cast& e) {
    PRINTERR_AND_EXIT(e.what());
  }
  DEBUGprint_FUNCend();
}

void Annotation::setValuesGV(const Variables &values)
{
  DEBUGprint_FUNCStart();

  try {
    if (values.count("inter")) {
      for (auto &x: getVal<std::vector<std::string>>(values, "inter")) {
        std::vector<std::string> v;
        boost::split(v, x, boost::algorithm::is_any_of(","));
        std::string label(v[0]);
        std::string tool("mango");
        if (v.size() >= 2) label = v[1];
        if (v.size() >= 3) tool  = v[2];
        vinterlist.emplace_back(InteractionSet(v[0], label, tool));
      }
    }
    if (values.count("ideogram")) {
      std::ifstream in(getVal<std::string>(values, "ideogram"));
      if (!in) {
        PRINTERR_AND_EXIT("Error: ideogram file " << getVal<std::string>(values, "ideogram") << " does not exist.");
      }
      while (!in.eof()) {
        std::string lineStr;
        getline(in, lineStr);
        if (lineStr.empty() || lineStr[0] == '#') continue;
        std::vector<std::string> v;
        ParseLine(v, lineStr, '\t');
        vcytoband.emplace_back(v);
      }
      isIdeogram = true;
    }
    if (values.count("chiadrop")) parse_ChIADropData(getVal<std::string>(values, "chiadrop"));
    chia_distance_thre = getVal<int32_t>(values, "chia_distance_thre");

    if (values.count("mp")) mpfile = getVal<std::string>(values, "mp");
    mpthre = getVal<double>(values, "mpthre");
    if (values.count("gap")) gapfile = getVal<std::string>(values, "gap");
    if (values.count("GC")) GC.setValue(getVal<std::string>(values, "GC"), getVal<int32_t>(values, "gcsize"));
    if (values.count("GD")) GD.setValue(getVal<std::string>(values, "GD"), getVal<int32_t>(values, "gdsize"));

  } catch (const boost::bad_any_cast& e) {
    PRINTERR_AND_EXIT(e.what());
  }
  DEBUGprint_FUNCend();
}

void Annotation::InitDumpPC(const Variables &values) const {
  DEBUGprint_FUNCStart();

  std::vector<std::string> str_gftype = {"refFlat", "GTF", "SGD"};
  std::cout << boost::format("\nAnnotations:\n");
  if (genefile != "") std::cout << boost::format("   Gene file: %1%, Format: %2%\n") % genefile % str_gftype[gftype];
  printOpt<std::string>(values, "ars",    "   ARS file");
  printOpt<std::string>(values, "ter",    "   TER file");
  printOpt<std::string>(values, "repeat", "   Repeat file");
  printOpt<std::string>(values, "region", "   Region file");
  printVOpt<std::string>(values, "bed",   "   Bed file");
  DEBUGprint_FUNCend();
}

void Annotation::InitDumpGV(const Variables &values) const {
  DEBUGprint_FUNCStart();

  printVOpt<std::string>(values, "inter", "   Interaction file");
  if (mpfile != "") {
    std::cout << boost::format("Mappability file directory: %1%\n") % mpfile;
    std::cout << boost::format("\tLow mappablitiy threshold: %1$2f\n") % mpthre;
  }
  if (isChIADrop) {
    std::cout << boost::format("ChIA-Drop file: %1%\n") % getVal<std::string>(values, "chiadrop");
    std::cout << boost::format("\tMax distance threshold: %1%\n") % chia_distance_thre;
  }
  printOpt<std::string>(values, "gc", "   GCcontents file");
  printOpt<std::string>(values, "gd", "   Gene density file");
  DEBUGprint_FUNCend();
}

void Profile::setOpts(MyOpt::Opts &allopts)
{
  MyOpt::Opts opt("PROFILE AND HEATMAP",100);
  opt.add_options()
    (SETOPT_RANGE("ptype", int32_t, TSS, TSS, PTYPENUM-1),
     "Region: 0; around TSS, 1; around TES, 2; divide gene into 100 subregions 3; around peak sites")
    (SETOPT_RANGE("stype", int32_t, 0, 0, 1),
     "Distribution: 0; ChIP read, 1; ChIP/Input enrichment")
//    (SETOPT_RANGE("ntype", int32_t, 0, 0, 1),
//     "Normalization: 0; total read 1; target regions only")
    (SETOPT_OVER("widthfromcenter", int32_t, 2500, 1), "width from the center")
    ("fixedlengthfromgene", "(for --ptype 2) use fixed length (specied by --widthfromcenter) from a gene")
    ;
  allopts.add(opt);
}

void Profile::setValues(const Variables &values)
{
  DEBUGprint_FUNCStart();

  try {
    ptype = getVal<int32_t>(values, "ptype");
    stype = getVal<int32_t>(values, "stype");
 //   ntype = getVal<int32_t>(values, "ntype");
    width_from_center = getVal<int32_t>(values, "widthfromcenter");
    usefixedlength = values.count("fixedlengthfromgene");
  } catch (const boost::bad_any_cast& e) {
    PRINTERR_AND_EXIT(e.what());
  }

  DEBUGprint_FUNCend();
}

void Profile::InitDump() const {
  DEBUGprint_FUNCStart();

  std::vector<std::string> str_ptype = { "TSS", "TTS", "GENE100", "BEDSITES" };
  std::vector<std::string> str_stype = { "ChIP read", "Enrichment ratio", "Enrichment P-value" };
//  std::vector<std::string> str_ntype = { "WHOLE GENOME", "TARGET REGIONS ONLY" };

  std::cout << boost::format("   Profile type: %1%\n")  % str_ptype[ptype];
  std::cout << boost::format("   Show type: %1%\n")     % str_stype[stype];
//  std::cout << boost::format("   Normalization: %1%\n") % str_ntype[ntype];
//  std::cout << boost::format("   Width from center: %1%\n") % width_from_center;

  DEBUGprint_FUNCend();
}


void Threshold::setOpts(MyOpt::Opts &allopts, const CommandParamSet &cps)
{
  sigtest = false;

  MyOpt::Opts opt("Threshold",100);
  opt.add_options()
    ("callpeak",       "Implement peak-calling (default: skip)")
    ("pthre_internal", value<double>()->default_value(cps.thre_pinter), "p-value threshold (-log10) for ChIP internal ")
    ("pthre_enrich",   value<double>()->default_value(cps.thre_penrich), "p-value threshold (-log10) for ChIP/Input enrichment")
    ("ethre",          value<double>()->default_value(cps.thre_ethre), "IP/Input fold enrichment")
    ("ipm",            value<double>()->default_value(cps.thre_ipm), "Read intensity of peak summit")
    ;
  allopts.add(opt);
}

void Threshold::setValues(const Variables &values)
{
  DEBUGprint_FUNCStart();
  try {
    pthre_inter  = getVal<double>(values, "pthre_internal");
    pthre_enrich = getVal<double>(values, "pthre_enrich");
    ethre        = getVal<double>(values, "ethre");
    ipm          = getVal<double>(values, "ipm");
    sigtest      = values.count("callpeak");
  } catch (const boost::bad_any_cast& e) {
    PRINTERR_AND_EXIT(e.what());
  }
  DEBUGprint_FUNCend();
}

void Threshold::InitDump() const {
  DEBUGprint_FUNCStart();

  if (sigtest) {
    std::cout << boost::format("\nThreshold:\n");
    std::cout << boost::format("   p-value threshold (internal, -log10): %1$.1f\n")   % pthre_inter;
    std::cout << boost::format("   p-value threshold (enrichment, -log10): %1$.1f\n") % pthre_enrich;
    std::cout << boost::format("   Peak intensity threshold: %1$.2f\n")               % ipm;
    std::cout << boost::format("   Enrichment threshold: %1$.2f\n")                   % ethre;
  }
  DEBUGprint_FUNCend();
}

void DrawRegion::setOpts(MyOpt::Opts &allopts)
{
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

void DrawRegion::setValues(const Variables &values)
{
  DEBUGprint_FUNCStart();

  try {
    for (auto x: {"region", "genelocifile"}) if (values.count(x)) isFile(getVal<std::string>(values, x));
    if (values.count("chr")) chr = getVal<std::string>(values, "chr");
    if (values.count("region")) {
      isRegion = true;
      regionBed = parseBed<bed>(getVal<std::string>(values, "region"));
      if (!regionBed.size()) PRINTERR_AND_EXIT("Error no bed regions in " << getVal<std::string>(values, "region"));
      //      printBed(regionBed);
    }
    if (values.count("genelocifile")) {
      getGeneLoci(getVal<std::string>(values, "genelocifile"));
      if (!values.count("gene")) PRINTERR_AND_EXIT("Please specify --gene option when supplying --genelocifile.");
    }
    len_geneloci = getVal<int32_t>(values, "len_geneloci");

  } catch (const boost::bad_any_cast& e) {
    PRINTERR_AND_EXIT(e.what());
  }
  DEBUGprint_FUNCend();
}

void DrawRegion::InitDump(const Variables &values) const {
  DEBUGprint_FUNCStart();

  printOpt<std::string>(values, "region",    "   Region file");
  //  if (chr != "") std::cout << boost::format("   output chr%1% only.\n") % chr;
  if (values.count("genelocifile")) {
    std::cout << boost::format("   Geneloci file: %1%, around %2% bp\n") % getVal<std::string>(values, "genelocifile") % len_geneloci;
  }
  DEBUGprint_FUNCend();
}

void Global::setOpts(const std::vector<DrompaCommand> &st, const CommandParamSet &cps)
{
  MyOpt::Opts o("Required",100);
  o.add_options()
    ("output,o",  value<std::string>(), "Output prefix")
    ("gt",        value<std::string>(), "Genome table")
    ;
  opts.add(o);
  for (auto &x: st) {
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
    case DrompaCommand::NORM:    setOptsNorm(opts, cps.sm);    break;
    case DrompaCommand::THRE:    thre.setOpts(opts, cps);      break;
    case DrompaCommand::ANNO_PC: anno.setOptsPC(opts);         break;
    case DrompaCommand::ANNO_GV: anno.setOptsGV(opts);         break;
    case DrompaCommand::DRAW:    drawparam.setOpts(opts, cps); break;
    case DrompaCommand::REGION:  drawregion.setOpts(opts);     break;
    case DrompaCommand::PROF:    prof.setOpts(opts);           break;
    case DrompaCommand::MULTICI:
      {
        options_description o("MULTICI",100);
        o.add_options()
          ("maxvalue",  "output the max bin value (default: averaged bin value) ")
          ("addname",  "output peakname in addition to genomic position (the 4th column required)")
          ;
        opts.add(o);
        break;
      }
    case DrompaCommand::GENWIG:
      {
        options_description o("GENWIG",100);
        o.add_options()
          ("outputformat",
           boost::program_options::value<int32_t>()->default_value(3)->notifier(std::bind(&MyOpt::range<int32_t>, std::placeholders::_1, 0, static_cast<int>(WigType::WIGTYPENUM) -2, "--outputformat")),
           "Output format\n   0: compressed wig (.wig.gz)\n   1: uncompressed wig (.wig)\n   2: bedGraph (.bedGraph)\n   3: bigWig (.bw)")
          ("outputvalue",
           boost::program_options::value<int32_t>()->default_value(0)->notifier(std::bind(&MyOpt::range<int32_t>, std::placeholders::_1, 0, 2, "--outputvalue")),
           "Output value\n   0: ChIP/Input enrichment\n   1: P-value (ChIP internal)\n   2: P-value (ChIP/Input enrichment)")
          ;
        opts.add(o);
        break;
      }
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
    case DrompaCommand::OTHER: setOptsOther(opts); break;
    }
  }
}

void Global::setValues(const std::vector<DrompaCommand> &vopts, const Variables &values)
{
  for (auto x: {"output", "gt"}) if (!values.count(x)) PRINTERR_AND_EXIT("specify --" << x << " option.");

  oprefix = getVal<std::string>(values, "output");
  genometablefilename = getVal<std::string>(values, "gt");
  gt = readGenomeTable(genometablefilename);

  for (auto op: vopts) {
    switch(op) {
    case DrompaCommand::CHIP:
      {
        DEBUGprint("ChIP setValues...");
        try {
          // SamplePairOverlayed first
          if (values.count("input")) {
            auto v(getVal<std::vector<std::string>>(values, "input"));
            for (auto &x: v) {
              vsinfo.addSampleInfo(x, gt, iftype);
              samplepair.emplace_back(x, vsinfo);
            }
          }

          // SamplePairOverlayed second
          if (values.count("ioverlay")) {
            auto v(getVal<std::vector<std::string>>(values, "ioverlay"));
            int32_t num(0);
            for (auto &x: v) {
              vsinfo.addSampleInfo(x, gt, iftype);
              samplepair[num].setSecondSample(x, vsinfo);
              ++num;
            }
          }

          iftype = static_cast<WigType>(getVal<int32_t>(values, "if"));

        } catch(const boost::bad_any_cast& e) {
          PRINTERR_AND_EXIT(e.what());
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
    case DrompaCommand::MULTICI:
      {
        DEBUGprint("Global::setValues::MULTICI");
        getmaxval = values.count("maxvalue");
        addname = values.count("addname");
        break;
      }
    case DrompaCommand::GENWIG:
      {
        DEBUGprint("Global::setValues::GENWIG");
        genwig_oftype = static_cast<WigType>(MyOpt::getVal<int32_t>(values, "outputformat"));
        genwig_ofvalue = MyOpt::getVal<int32_t>(values, "outputvalue");
        break;
      }
    case DrompaCommand::CG:
      {
        DEBUGprint("Global::setValues::CG");
        break;
      }
    case DrompaCommand::TR:
      {
        DEBUGprint("Global::setValues::TR");
        break;
      }
    case DrompaCommand::OTHER: setValuesOther(values); break;
    }
  }

  return;
}

void Global::setOptsNorm(MyOpt::Opts &allopts, const int32_t defsm)
{
  MyOpt::Opts o("Normalization",100);
  o.add_options()
    (SETOPT_RANGE("norm", int32_t, 0, 0, 2),
     "Normalization between ChIP and Input\n      0: not normalize\n      1: with total read number (genome)\n      2: with total read number (each chr)\n      3: with NCIS method\n")
    (SETOPT_OVER("sm", int32_t, defsm, 0), "# of bins for Gausian smoothing")
    ;
  allopts.add(o);
}

void Global::setValuesNorm(const Variables &values) {
  DEBUGprint_FUNCStart();
  try {
    norm = getVal<int32_t>(values, "norm");
    smoothing = getVal<int32_t>(values, "sm");
  } catch (const boost::bad_any_cast& e) {
    PRINTERR_AND_EXIT(e.what());
  }
  DEBUGprint_FUNCend();
}

void Global::setOptsOther(MyOpt::Opts &allopts)
{
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

void Global::setValuesOther(const Variables &values)
{
  DEBUGprint_FUNCStart();
  try {
    includeYM = values.count("includeYM");
    ispng = values.count("png");
    showchr = values.count("showchr");
  } catch(const boost::bad_any_cast& e) {
    PRINTERR_AND_EXIT(e.what());
  }
  DEBUGprint_FUNCend();
}

void Global::InitDumpChIP() const {
  DEBUGprint_FUNCStart();

  std::vector<std::string> str_wigfiletype = {"BINARY", "COMPRESSED WIG", "WIG", "BEDGRAPH", "BIGWIG"};
  printf("\nSamples\n");
  for (size_t i=0; i<samplepair.size(); ++i) {
    std::cout << (i+1) << ": ";
    samplepair[i].print();
  }
  if (iftype < WigType::NONE) {
    std::cout << boost::format("Input format: %1%\n") % str_wigfiletype[static_cast<int32_t>(iftype)];
  }
  DEBUGprint_FUNCend();
}

void Global::InitDumpNorm() const {
  DEBUGprint_FUNCStart();

  std::vector<std::string> str_norm = { "OFF", "TOTALREAD GENOME", "TOTALREAD CHR", "NCIS" };
  std::cout << boost::format("   ChIP/Input normalization: %1%\n") % str_norm[norm];
  if (smoothing) std::cout << boost::format("   smoothing width: %1% bins\n") % smoothing;
  DEBUGprint_FUNCend();
}

void Global::InitDumpOther() const {
  DEBUGprint_FUNCStart();

  std::vector<std::string> str_format = {"PDF", "PNG"};
  if (drawparam.isshowpdf()) std::cout << boost::format("   Output format: %1%\n") % str_format[ispng];
  else std::cout << boost::format("   Output format: do not depict figure files.\n");
  if (includeYM) std::cout << boost::format("   include chromosome Y and M\n");
  DEBUGprint_FUNCend();
}

void DrawParam::setOpts(MyOpt::Opts &allopts, const CommandParamSet &cps)
{
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
    ("shownegative", "allow negative values in historgram")
    ("offymem", "Omit Y memory")
    ("offylabel", "Omit Y label")
    ("offbedname", "Omit name in BED12")
    ("offpdf", "Omit pdf generation (for peak calling)")
    ;
  allopts.add(opt);
}

void DrawParam::setValues(const Variables &values, const int32_t n)
{
  DEBUGprint_FUNCStart();

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
    shownegative = values.count("shownegative");
    showymem = !values.count("offymem");
    showylab = !values.count("offylabel");
    showbedname = !values.count("offbedname");
    showpdf  = !values.count("offpdf");
    alpha = getVal<double>(values, "alpha");

    scale_tag    = getVal<double>(values, "scale_tag");
    scale_ratio  = getVal<double>(values, "scale_ratio");
    scale_pvalue = getVal<double>(values, "scale_pvalue");

    samplenum = n;
  } catch (const boost::bad_any_cast& e) {
    PRINTERR_AND_EXIT(e.what());
  }
  DEBUGprint_FUNCend();
}

void DrawParam::InitDump() const {
  DEBUGprint_FUNCStart();

  std::vector<std::string> str_bool = {"OFF", "ON"};
  std::vector<std::string> str_input = {"OFF", "ALL", "FIRST"};
  std::vector<std::string> str_ratio = {"OFF", "Linear", "Logratio"};

  std::cout << boost::format("\nFigure parameter:\n");
  std::cout << boost::format("   Display read: ChIP %1%, Input %2%, y-axis scale: %3%\n") % str_bool[showctag] % str_input[showitag] % scale_tag;
  std::cout << boost::format("   Display enrichment: %1%, y-axis scale: %2%\n")           % str_ratio[showratio] % scale_ratio;
  std::cout << boost::format("   Display pvalue (internal): %1%, y-axis scale: %2%\n")    % str_bool[showpinter] % scale_pvalue;
  std::cout << boost::format("   Display pvalue (ChIP/Input): %1%, y-axis scale: %2%\n")  % str_bool[showpenrich] % scale_pvalue;
  std::cout << boost::format("   Width per line: %1% kbp\n")           % (width_per_line/1000);
  std::cout << boost::format("   Y-axis label: %1%\n")                 % str_bool[showylab];
  std::cout << boost::format("   Y-axis memory: %1%\n")                % str_bool[showymem];

  DEBUGprint("barnum " << barnum);
  DEBUGprint("ystep " << ystep);
  DEBUGprint_FUNCend();
}
