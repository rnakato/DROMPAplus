/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _READANNOTATION_HPP_
#define _READANNOTATION_HPP_

#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include "SSP/common/seq.hpp"

using GeneDataMap = std::unordered_map<std::string, genedata>;
using HashOfGeneDataMap = std::unordered_map<std::string, GeneDataMap>;

int32_t countmp(HashOfGeneDataMap &);
std::vector<std::string> scanGeneName(const HashOfGeneDataMap &);
HashOfGeneDataMap extract_mp(const HashOfGeneDataMap &, const std::vector<std::string>);
std::vector<std::string> readGeneList(const std::string&);
HashOfGeneDataMap parseRefFlat(const std::string&);
HashOfGeneDataMap parseGtf(const std::string&);
HashOfGeneDataMap construct_gmp(const HashOfGeneDataMap &);
void printMap(const HashOfGeneDataMap &);
void printRefFlat(const HashOfGeneDataMap &, const int32_t nameflag);
std::vector<chrsize> read_genometable(const std::string&);

#endif  // _READANNOTATION_HPP_
