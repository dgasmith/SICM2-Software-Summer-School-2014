//
// author: Daniel Smith (dsmith@auburn.edu)
// date  : July 14, 2014
// the use of this software is permitted under the conditions GNU General Public License (GPL) version 2
//

#ifndef __s2i2_core_blockkernel_h_DEFINED
#define __s2i2_core_blockkernel_h_DEFINED

// standard C++ headers
#include <cstddef>

void dgemm_blocked(double* a, double* b, double* c, size_t n, size_t bsize);

#endif // __s2i2_core_blockkernel_h_DEFINED
