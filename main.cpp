#include <cppcrystal/cppcrystal.hpp>

#include <cstdio>

int main() {
  std::printf("CppCrystal %s — modern C++23 crystallography library "
              "(validated against reference spglib v%d.%d.%d)\n",
              cppcrystal::version_string(), cppcrystal::kReferenceSpglibVersion.major,
              cppcrystal::kReferenceSpglibVersion.minor,
              cppcrystal::kReferenceSpglibVersion.patch);
  return 0;
}
