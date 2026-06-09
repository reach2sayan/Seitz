#include <spglib/spglib.hpp>

#include <cstdio>

int main() {
  std::printf("CppCrystal %s — modern C++ port of spglib "
              "(targeting reference v%d.%d.%d)\n",
              spglib::version_string(), spglib::kReferenceSpglibVersion.major,
              spglib::kReferenceSpglibVersion.minor,
              spglib::kReferenceSpglibVersion.patch);
  return 0;
}
