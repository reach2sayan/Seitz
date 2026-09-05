#pragma once

// The library is built with hidden visibility: only declarations inside a
// public header's `#pragma GCC visibility push(default)` window reach the
// symbol table. The unit tests link that shared object and exercise pipeline
// stages behind the public API, so the internal entry points they call are
// exported with this.
//
// Not an API promise -- nothing under src/ is installed, so no consumer can
// declare these. It says only "the test binary links this".
//
// Empty on MSVC, where the library is a static archive (see CMakeLists.txt):
#if defined(_MSC_VER)
#define SEITZ_TESTABLE
#else
#define SEITZ_TESTABLE __attribute__((visibility("default")))
#endif
