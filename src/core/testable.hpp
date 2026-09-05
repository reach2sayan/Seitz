#pragma once

// The library is built with hidden visibility: only the declarations inside a
// public header's `#pragma GCC visibility push(default)` window reach the
// shared object's symbol table. The unit tests link that shared object and
// deliberately exercise the pipeline stages behind the public API, so the
// internal entry points they call are exported explicitly with this.
//
// It is not an API promise -- nothing under src/ is installed, so no consumer
// can even declare these. It says only "the test binary links this".
//
// Empty on MSVC, where the library is a static archive (see CMakeLists.txt):
#if defined(_MSC_VER)
#define CPPCRYSTAL_TESTABLE
#else
#define CPPCRYSTAL_TESTABLE __attribute__((visibility("default")))
#endif
