# MIT License
#
# Copyright (c) 2022 Ignacio Vizzo, Tiziano Guadagnino, Benedikt Mersch, Cyrill
# Stachniss.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
include(ExternalProject)

ExternalProject_Add(
  external_pmc
  PREFIX pmc
  URL https://github.com/LimHyungTae/pmc/archive/refs/tags/libpmc.tar.gz
  URL_HASH SHA256=64ea6e628fe0920df771d21a243df6a771ebde9221042fe203390ff1520e969d
  UPDATE_COMMAND ""
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
)

ExternalProject_Get_Property(external_pmc INSTALL_DIR)
add_library(PMCHelper INTERFACE)
add_dependencies(PMCHelper external_pmc)
target_include_directories(PMCHelper INTERFACE ${INSTALL_DIR}/include)
target_link_directories(PMCHelper INTERFACE ${INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR})
target_link_libraries(PMCHelper INTERFACE pmc)
set_property(TARGET PMCHelper PROPERTY EXPORT_NAME pmc::pmc)
add_library(pmc::pmc ALIAS PMCHelper)
