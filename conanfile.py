import os
import re

from conan import ConanFile
from conan.errors import ConanInvalidConfiguration
from conan.tools.build import check_min_cppstd
from conan.tools.cmake import CMake, CMakeDeps, CMakeToolchain, cmake_layout
from conan.tools.files import copy, load, rmdir

# Builds this checkout: `conan create .` from the repository root. Eigen and
# Boost come from conan (SEITZ_EXTERNAL_EIGEN_AND_BOOST, the same hook the vcpkg
# port uses). The symmetry tables are transcribed at configure time from
# spglib's v2.7.0 sources, which FetchContent clones -- an installed spglib
# carries no sources -- by numpy/pandas programs, run with a venv built in the
# build folder. Those two steps need the network at build time.


class SeitzConan(ConanFile):
    name = "seitz"
    description = "Crystal symmetry: a modern C++23 port of spglib"
    license = "BSD-3-Clause"
    url = "https://github.com/reach2sayan/Seitz"
    homepage = url
    topics = ("crystallography", "symmetry", "spglib", "space-group")
    package_type = "library"
    settings = "os", "arch", "compiler", "build_type"
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {
        "shared": False,
        "fPIC": True,
        # container and graph are the only compiled Boost libraries seitz links
        # (graph pulls in math, random, regex, serialization); leaf, parser,
        # flyweight and algorithm are header-only. Build nothing else.
        **{f"boost/*:without_{lib}": True for lib in (
            "atomic", "charconv", "chrono", "cobalt", "context", "contract",
            "coroutine", "date_time", "exception", "fiber", "filesystem",
            "iostreams", "json", "locale", "log", "nowide", "process",
            "program_options", "stacktrace", "test", "thread", "timer",
            "type_erasure", "url", "wave")},
    }
    exports_sources = ("CMakeLists.txt", "LICENSE", "cmake/*", "include/*",
                       "src/*", "tools/*", "pyproject.toml")

    def set_version(self):
        cmake = load(self, os.path.join(self.recipe_folder, "CMakeLists.txt"))
        self.version = re.search(r"project\(Seitz VERSION (\S+)", cmake)[1]

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def configure(self):
        if self.options.shared:
            self.options.rm_safe("fPIC")

    def layout(self):
        cmake_layout(self)

    def requirements(self):
        # Public headers include both, so both are transitive.
        self.requires("eigen/5.0.0", transitive_headers=True)
        self.requires("boost/1.88.0", transitive_headers=True)

    def validate(self):
        check_min_cppstd(self, 23)
        # Boost.LEAF keeps error payloads in thread_local statics with no usable
        # visibility attribute off GCC, so across a DLL boundary every typed
        # error would arrive as the generic base. CMakeLists.txt forces STATIC
        # on Windows for that reason.
        if self.settings.os == "Windows" and self.options.shared:
            raise ConanInvalidConfiguration("seitz is static-only on Windows")

    def _venv_python(self):
        venv = os.path.join(self.build_folder, "transcribe-venv")
        python = os.path.join(venv, "Scripts" if self.settings.os == "Windows" else "bin",
                              "python")
        if not os.path.exists(venv):
            self.run(f"python3 -m venv {venv}")
            self.run(f"{python} -m pip install --quiet --disable-pip-version-check numpy pandas")
        return python

    def generate(self):
        tc = CMakeToolchain(self)
        tc.cache_variables.update({
            "SEITZ_EXTERNAL_EIGEN_AND_BOOST": True,
            "SEITZ_STATIC_LIBRARY": not self.options.shared,
            "SEITZ_BUILD_TESTS": False,
            "SEITZ_BUILD_DEMO": False,
            "SEITZ_BUILD_ORACLE_TESTS": False,
            "Python_EXECUTABLE": self._venv_python().replace("\\", "/"),
        })
        tc.generate()
        CMakeDeps(self).generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        copy(self, "LICENSE", self.source_folder, os.path.join(self.package_folder, "licenses"))
        CMake(self).install()
        rmdir(self, os.path.join(self.package_folder, "lib", "cmake"))

    def package_info(self):
        self.cpp_info.set_property("cmake_file_name", "Seitz")
        self.cpp_info.set_property("cmake_target_name", "seitz::seitz")
        self.cpp_info.libs = ["seitz"]
        self.cpp_info.requires = ["eigen::eigen", "boost::headers", "boost::container",
                                  "boost::graph"]
        if self.settings.os in ("Linux", "FreeBSD"):
            self.cpp_info.system_libs = ["pthread", "m"]
