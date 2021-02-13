/* Copyright(c)  Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _READANNOTATION_HPP_
#define _READANNOTATION_HPP_

#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include "../common/GeneAnnotation.hpp"

using GeneDataMap = std::unordered_map<std::string, genedata>;
using HashOfGeneDataMap = std::unordered_map<std::string, GeneDataMap>;

int32_t countmp(HashOfGeneDataMap &);
std::vector<std::string> scanGeneName(const HashOfGeneDataMap &);
HashOfGeneDataMap extract_mp(const HashOfGeneDataMap &, const std::vector<std::string>);
std::vector<std::string> readGeneList(const std::string&);
void parseARSOriDB(const std::string& fileName, HashOfGeneDataMap &mp);
void parseTER(const std::string& fileName, HashOfGeneDataMap &mp);
HashOfGeneDataMap parseSGD(const std::string&);
HashOfGeneDataMap parseRefFlat(const std::string&);
HashOfGeneDataMap parseGtf(const std::string&);
HashOfGeneDataMap construct_gmp(const HashOfGeneDataMap &);
void printMap(const HashOfGeneDataMap &);
bool isGeneUCSC(const HashOfGeneDataMap &);
std::vector<chrsize> readGenomeTable(const std::string&);

#endif  // _READANNOTATION_HPP_
