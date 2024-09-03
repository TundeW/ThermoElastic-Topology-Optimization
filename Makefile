# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/root/DEV/FEA Based Solver/FEA Based Solver"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/root/DEV/FEA Based Solver/FEA Based Solver"

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/usr/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start "/root/DEV/FEA Based Solver/FEA Based Solver/CMakeFiles" "/root/DEV/FEA Based Solver/FEA Based Solver//CMakeFiles/progress.marks"
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start "/root/DEV/FEA Based Solver/FEA Based Solver/CMakeFiles" 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named HeatExchanger

# Build rule for target.
HeatExchanger: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 HeatExchanger
.PHONY : HeatExchanger

# fast build rule for target.
HeatExchanger/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/build
.PHONY : HeatExchanger/fast

src/ElasticFeaFunction.o: src/ElasticFeaFunction.cpp.o
.PHONY : src/ElasticFeaFunction.o

# target to build an object file
src/ElasticFeaFunction.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/ElasticFeaFunction.cpp.o
.PHONY : src/ElasticFeaFunction.cpp.o

src/ElasticFeaFunction.i: src/ElasticFeaFunction.cpp.i
.PHONY : src/ElasticFeaFunction.i

# target to preprocess a source file
src/ElasticFeaFunction.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/ElasticFeaFunction.cpp.i
.PHONY : src/ElasticFeaFunction.cpp.i

src/ElasticFeaFunction.s: src/ElasticFeaFunction.cpp.s
.PHONY : src/ElasticFeaFunction.s

# target to generate assembly for a file
src/ElasticFeaFunction.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/ElasticFeaFunction.cpp.s
.PHONY : src/ElasticFeaFunction.cpp.s

src/FeaHelperFunctions.o: src/FeaHelperFunctions.cpp.o
.PHONY : src/FeaHelperFunctions.o

# target to build an object file
src/FeaHelperFunctions.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/FeaHelperFunctions.cpp.o
.PHONY : src/FeaHelperFunctions.cpp.o

src/FeaHelperFunctions.i: src/FeaHelperFunctions.cpp.i
.PHONY : src/FeaHelperFunctions.i

# target to preprocess a source file
src/FeaHelperFunctions.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/FeaHelperFunctions.cpp.i
.PHONY : src/FeaHelperFunctions.cpp.i

src/FeaHelperFunctions.s: src/FeaHelperFunctions.cpp.s
.PHONY : src/FeaHelperFunctions.s

# target to generate assembly for a file
src/FeaHelperFunctions.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/FeaHelperFunctions.cpp.s
.PHONY : src/FeaHelperFunctions.cpp.s

src/FowardAnalysisFunctions.o: src/FowardAnalysisFunctions.cpp.o
.PHONY : src/FowardAnalysisFunctions.o

# target to build an object file
src/FowardAnalysisFunctions.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/FowardAnalysisFunctions.cpp.o
.PHONY : src/FowardAnalysisFunctions.cpp.o

src/FowardAnalysisFunctions.i: src/FowardAnalysisFunctions.cpp.i
.PHONY : src/FowardAnalysisFunctions.i

# target to preprocess a source file
src/FowardAnalysisFunctions.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/FowardAnalysisFunctions.cpp.i
.PHONY : src/FowardAnalysisFunctions.cpp.i

src/FowardAnalysisFunctions.s: src/FowardAnalysisFunctions.cpp.s
.PHONY : src/FowardAnalysisFunctions.s

# target to generate assembly for a file
src/FowardAnalysisFunctions.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/FowardAnalysisFunctions.cpp.s
.PHONY : src/FowardAnalysisFunctions.cpp.s

src/MeshHelperFunctions.o: src/MeshHelperFunctions.cpp.o
.PHONY : src/MeshHelperFunctions.o

# target to build an object file
src/MeshHelperFunctions.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/MeshHelperFunctions.cpp.o
.PHONY : src/MeshHelperFunctions.cpp.o

src/MeshHelperFunctions.i: src/MeshHelperFunctions.cpp.i
.PHONY : src/MeshHelperFunctions.i

# target to preprocess a source file
src/MeshHelperFunctions.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/MeshHelperFunctions.cpp.i
.PHONY : src/MeshHelperFunctions.cpp.i

src/MeshHelperFunctions.s: src/MeshHelperFunctions.cpp.s
.PHONY : src/MeshHelperFunctions.s

# target to generate assembly for a file
src/MeshHelperFunctions.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/MeshHelperFunctions.cpp.s
.PHONY : src/MeshHelperFunctions.cpp.s

src/ProjectionHelperFunctions.o: src/ProjectionHelperFunctions.cpp.o
.PHONY : src/ProjectionHelperFunctions.o

# target to build an object file
src/ProjectionHelperFunctions.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/ProjectionHelperFunctions.cpp.o
.PHONY : src/ProjectionHelperFunctions.cpp.o

src/ProjectionHelperFunctions.i: src/ProjectionHelperFunctions.cpp.i
.PHONY : src/ProjectionHelperFunctions.i

# target to preprocess a source file
src/ProjectionHelperFunctions.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/ProjectionHelperFunctions.cpp.i
.PHONY : src/ProjectionHelperFunctions.cpp.i

src/ProjectionHelperFunctions.s: src/ProjectionHelperFunctions.cpp.s
.PHONY : src/ProjectionHelperFunctions.s

# target to generate assembly for a file
src/ProjectionHelperFunctions.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/ProjectionHelperFunctions.cpp.s
.PHONY : src/ProjectionHelperFunctions.cpp.s

src/ThermalFeaFunction.o: src/ThermalFeaFunction.cpp.o
.PHONY : src/ThermalFeaFunction.o

# target to build an object file
src/ThermalFeaFunction.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/ThermalFeaFunction.cpp.o
.PHONY : src/ThermalFeaFunction.cpp.o

src/ThermalFeaFunction.i: src/ThermalFeaFunction.cpp.i
.PHONY : src/ThermalFeaFunction.i

# target to preprocess a source file
src/ThermalFeaFunction.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/ThermalFeaFunction.cpp.i
.PHONY : src/ThermalFeaFunction.cpp.i

src/ThermalFeaFunction.s: src/ThermalFeaFunction.cpp.s
.PHONY : src/ThermalFeaFunction.s

# target to generate assembly for a file
src/ThermalFeaFunction.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/ThermalFeaFunction.cpp.s
.PHONY : src/ThermalFeaFunction.cpp.s

src/WriteOutputFunctions.o: src/WriteOutputFunctions.cpp.o
.PHONY : src/WriteOutputFunctions.o

# target to build an object file
src/WriteOutputFunctions.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/WriteOutputFunctions.cpp.o
.PHONY : src/WriteOutputFunctions.cpp.o

src/WriteOutputFunctions.i: src/WriteOutputFunctions.cpp.i
.PHONY : src/WriteOutputFunctions.i

# target to preprocess a source file
src/WriteOutputFunctions.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/WriteOutputFunctions.cpp.i
.PHONY : src/WriteOutputFunctions.cpp.i

src/WriteOutputFunctions.s: src/WriteOutputFunctions.cpp.s
.PHONY : src/WriteOutputFunctions.s

# target to generate assembly for a file
src/WriteOutputFunctions.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/WriteOutputFunctions.cpp.s
.PHONY : src/WriteOutputFunctions.cpp.s

src/WriteVtuFunction.o: src/WriteVtuFunction.cpp.o
.PHONY : src/WriteVtuFunction.o

# target to build an object file
src/WriteVtuFunction.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/WriteVtuFunction.cpp.o
.PHONY : src/WriteVtuFunction.cpp.o

src/WriteVtuFunction.i: src/WriteVtuFunction.cpp.i
.PHONY : src/WriteVtuFunction.i

# target to preprocess a source file
src/WriteVtuFunction.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/WriteVtuFunction.cpp.i
.PHONY : src/WriteVtuFunction.cpp.i

src/WriteVtuFunction.s: src/WriteVtuFunction.cpp.s
.PHONY : src/WriteVtuFunction.s

# target to generate assembly for a file
src/WriteVtuFunction.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/WriteVtuFunction.cpp.s
.PHONY : src/WriteVtuFunction.cpp.s

src/main.o: src/main.cpp.o
.PHONY : src/main.o

# target to build an object file
src/main.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/main.cpp.o
.PHONY : src/main.cpp.o

src/main.i: src/main.cpp.i
.PHONY : src/main.i

# target to preprocess a source file
src/main.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/main.cpp.i
.PHONY : src/main.cpp.i

src/main.s: src/main.cpp.s
.PHONY : src/main.s

# target to generate assembly for a file
src/main.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/HeatExchanger.dir/build.make CMakeFiles/HeatExchanger.dir/src/main.cpp.s
.PHONY : src/main.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... HeatExchanger"
	@echo "... src/ElasticFeaFunction.o"
	@echo "... src/ElasticFeaFunction.i"
	@echo "... src/ElasticFeaFunction.s"
	@echo "... src/FeaHelperFunctions.o"
	@echo "... src/FeaHelperFunctions.i"
	@echo "... src/FeaHelperFunctions.s"
	@echo "... src/FowardAnalysisFunctions.o"
	@echo "... src/FowardAnalysisFunctions.i"
	@echo "... src/FowardAnalysisFunctions.s"
	@echo "... src/MeshHelperFunctions.o"
	@echo "... src/MeshHelperFunctions.i"
	@echo "... src/MeshHelperFunctions.s"
	@echo "... src/ProjectionHelperFunctions.o"
	@echo "... src/ProjectionHelperFunctions.i"
	@echo "... src/ProjectionHelperFunctions.s"
	@echo "... src/ThermalFeaFunction.o"
	@echo "... src/ThermalFeaFunction.i"
	@echo "... src/ThermalFeaFunction.s"
	@echo "... src/WriteOutputFunctions.o"
	@echo "... src/WriteOutputFunctions.i"
	@echo "... src/WriteOutputFunctions.s"
	@echo "... src/WriteVtuFunction.o"
	@echo "... src/WriteVtuFunction.i"
	@echo "... src/WriteVtuFunction.s"
	@echo "... src/main.o"
	@echo "... src/main.i"
	@echo "... src/main.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

