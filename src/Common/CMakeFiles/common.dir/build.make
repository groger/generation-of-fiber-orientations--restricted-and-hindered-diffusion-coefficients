# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /biomed-resimg/crews_rodent/devel/linux/cmake-2.8.2/cmake-2.8.2/bin/cmake

# The command to remove a file.
RM = /biomed-resimg/crews_rodent/devel/linux/cmake-2.8.2/cmake-2.8.2/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /biomed-resimg/crews_rodent/devel/linux/cmake-2.8.2/cmake-2.8.2/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig

# Include any dependencies generated for this target.
include Common/CMakeFiles/common.dir/depend.make

# Include the progress variables for this target.
include Common/CMakeFiles/common.dir/progress.make

# Include the compile flags for this target's objects.
include Common/CMakeFiles/common.dir/flags.make

Common/CMakeFiles/common.dir/Exception.cc.o: Common/CMakeFiles/common.dir/flags.make
Common/CMakeFiles/common.dir/Exception.cc.o: Common/Exception.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object Common/CMakeFiles/common.dir/Exception.cc.o"
	cd /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig/Common && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/common.dir/Exception.cc.o -c /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig/Common/Exception.cc

Common/CMakeFiles/common.dir/Exception.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/common.dir/Exception.cc.i"
	cd /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig/Common && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig/Common/Exception.cc > CMakeFiles/common.dir/Exception.cc.i

Common/CMakeFiles/common.dir/Exception.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/common.dir/Exception.cc.s"
	cd /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig/Common && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig/Common/Exception.cc -o CMakeFiles/common.dir/Exception.cc.s

Common/CMakeFiles/common.dir/Exception.cc.o.requires:
.PHONY : Common/CMakeFiles/common.dir/Exception.cc.o.requires

Common/CMakeFiles/common.dir/Exception.cc.o.provides: Common/CMakeFiles/common.dir/Exception.cc.o.requires
	$(MAKE) -f Common/CMakeFiles/common.dir/build.make Common/CMakeFiles/common.dir/Exception.cc.o.provides.build
.PHONY : Common/CMakeFiles/common.dir/Exception.cc.o.provides

Common/CMakeFiles/common.dir/Exception.cc.o.provides.build: Common/CMakeFiles/common.dir/Exception.cc.o
.PHONY : Common/CMakeFiles/common.dir/Exception.cc.o.provides.build

Common/CMakeFiles/common.dir/IndexValuePair.cc.o: Common/CMakeFiles/common.dir/flags.make
Common/CMakeFiles/common.dir/IndexValuePair.cc.o: Common/IndexValuePair.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object Common/CMakeFiles/common.dir/IndexValuePair.cc.o"
	cd /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig/Common && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/common.dir/IndexValuePair.cc.o -c /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig/Common/IndexValuePair.cc

Common/CMakeFiles/common.dir/IndexValuePair.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/common.dir/IndexValuePair.cc.i"
	cd /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig/Common && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig/Common/IndexValuePair.cc > CMakeFiles/common.dir/IndexValuePair.cc.i

Common/CMakeFiles/common.dir/IndexValuePair.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/common.dir/IndexValuePair.cc.s"
	cd /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig/Common && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig/Common/IndexValuePair.cc -o CMakeFiles/common.dir/IndexValuePair.cc.s

Common/CMakeFiles/common.dir/IndexValuePair.cc.o.requires:
.PHONY : Common/CMakeFiles/common.dir/IndexValuePair.cc.o.requires

Common/CMakeFiles/common.dir/IndexValuePair.cc.o.provides: Common/CMakeFiles/common.dir/IndexValuePair.cc.o.requires
	$(MAKE) -f Common/CMakeFiles/common.dir/build.make Common/CMakeFiles/common.dir/IndexValuePair.cc.o.provides.build
.PHONY : Common/CMakeFiles/common.dir/IndexValuePair.cc.o.provides

Common/CMakeFiles/common.dir/IndexValuePair.cc.o.provides.build: Common/CMakeFiles/common.dir/IndexValuePair.cc.o
.PHONY : Common/CMakeFiles/common.dir/IndexValuePair.cc.o.provides.build

# Object files for target common
common_OBJECTS = \
"CMakeFiles/common.dir/Exception.cc.o" \
"CMakeFiles/common.dir/IndexValuePair.cc.o"

# External object files for target common
common_EXTERNAL_OBJECTS =

Common/libcommon.a: Common/CMakeFiles/common.dir/Exception.cc.o
Common/libcommon.a: Common/CMakeFiles/common.dir/IndexValuePair.cc.o
Common/libcommon.a: Common/CMakeFiles/common.dir/build.make
Common/libcommon.a: Common/CMakeFiles/common.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libcommon.a"
	cd /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig/Common && $(CMAKE_COMMAND) -P CMakeFiles/common.dir/cmake_clean_target.cmake
	cd /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig/Common && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/common.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Common/CMakeFiles/common.dir/build: Common/libcommon.a
.PHONY : Common/CMakeFiles/common.dir/build

Common/CMakeFiles/common.dir/requires: Common/CMakeFiles/common.dir/Exception.cc.o.requires
Common/CMakeFiles/common.dir/requires: Common/CMakeFiles/common.dir/IndexValuePair.cc.o.requires
.PHONY : Common/CMakeFiles/common.dir/requires

Common/CMakeFiles/common.dir/clean:
	cd /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig/Common && $(CMAKE_COMMAND) -P CMakeFiles/common.dir/cmake_clean.cmake
.PHONY : Common/CMakeFiles/common.dir/clean

Common/CMakeFiles/common.dir/depend:
	cd /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig/Common /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig/Common /biomed-resimg/NAMIC/DTITractographyPhantom/fiber-sig/Common/CMakeFiles/common.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Common/CMakeFiles/common.dir/depend

