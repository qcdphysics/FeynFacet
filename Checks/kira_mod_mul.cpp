// Release/CI diagnostic: link against the actual Kira build's libpyred.a.
// The reference uses repeated addition, independently of the patched __int128 path.
#include <cfloat>
#include <cstdint>
#include <iostream>
#include <random>
#include "pyred/coeff_helper.h"

static uint64_t reference(uint64_t a, uint64_t b, uint64_t m) {
  uint64_t result = 0;
  while (b) {
    if (b & 1) result = (result + a) % m;
    a = (a + a) % m; // m < 2^63, so additions cannot overflow uint64_t.
    b >>= 1;
  }
  return result;
}

int main() {
  std::mt19937_64 random(6172026);
  const uint64_t primes[] = {9223372036854775783ULL, 9223372036854775643ULL};
  unsigned errors = 0, count = 0;
  for (uint64_t m : primes) {
    const uint64_t edges[] = {0, 1, 2, m / 2, m - 2, m - 1};
    for (uint64_t a : edges) for (uint64_t b : edges) {
      ++count;
      errors += pyred::CoeffHelper::mod_mul(a, b, m) != reference(a, b, m);
    }
    for (unsigned i = 0; i < 50000; ++i) {
      const uint64_t a = random() % m, b = random() % m;
      ++count;
      errors += pyred::CoeffHelper::mod_mul(a, b, m) != reference(a, b, m);
    }
  }
  std::cout << "LDBL_MANT_DIG=" << LDBL_MANT_DIG << "; mismatches="
            << errors << "/" << count << "\n";
  return errors ? 1 : 0;
}
