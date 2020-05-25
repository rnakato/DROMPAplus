/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "util.hpp"
std::string basename(const std::string &path)
{
    return path.substr(path.find_last_of('/') + 1);
}
