# -*- Makefile -*-

#arch = ../setup/Make.Linux_MPI
setup_file = ../setup/Make.Linux_MPI

#include setup/Make.$(arch)
include $(setup_file)


HPCG_DEPS = src/CG.o \
	    src/CG_ref.o \
	    src/TestCG.o \
	    src/ComputeResidual.o \
	    src/ExchangeHalo.o \
	    src/GenerateGeometry.o \
	    src/GenerateProblem.o \
	    src/OptimizeProblem.o \
	    src/ReadHpcgDat.o \
	    src/ReportResults.o \
	    src/SetupHalo.o \
	    src/TestSymmetry.o \
	    src/TestNorms.o \
	    src/WriteProblem.o \
	    src/YAML_Doc.o \
	    src/YAML_Element.o \
	    src/ComputeDotProduct.o \
	    src/ComputeDotProduct_ref.o \
	    src/mytimer.o \
	    src/ComputeSPMV.o \
	    src/ComputeSPMV_ref.o \
	    src/ComputeSYMGS.o \
	    src/ComputeSYMGS_ref.o \
	    src/ComputeWAXPBY.o \
	    src/ComputeWAXPBY_ref.o \
	    src/ComputeMG_ref.o \
	    src/ComputeMG.o \
	    src/ComputeProlongation_ref.o \
	    src/ComputeRestriction_ref.o \
	    src/GenerateCoarseProblem.o \
	    src/init.o \
	    src/finalize.o

all: bin/xhpcg

bin/xhpcg: testing/main.o $(HPCG_DEPS)
	$(LINKER) $(LINKFLAGS) testing/main.o $(HPCG_DEPS) $(HPCG_LIBS) -o bin/xhpcg

clean:
	rm -f src/*.o testing/*.o bin/xhpcg

.PHONY: all clean

testing/main.o: ../testing/main.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/CG.o: ../src/CG.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/CG_ref.o: ../src/CG_ref.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/TestCG.o: ../src/TestCG.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/ComputeResidual.o: ../src/ComputeResidual.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/ExchangeHalo.o: ../src/ExchangeHalo.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/GenerateGeometry.o: ../src/GenerateGeometry.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/GenerateProblem.o: ../src/GenerateProblem.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/OptimizeProblem.o: ../src/OptimizeProblem.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/ReadHpcgDat.o: ../src/ReadHpcgDat.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/ReportResults.o: ../src/ReportResults.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/SetupHalo.o: ../src/SetupHalo.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/TestSymmetry.o: ../src/TestSymmetry.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/TestNorms.o: ../src/TestNorms.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/WriteProblem.o: ../src/WriteProblem.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/YAML_Doc.o: ../src/YAML_Doc.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/YAML_Element.o: ../src/YAML_Element.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/ComputeDotProduct.o: ../src/ComputeDotProduct.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/ComputeDotProduct_ref.o: ../src/ComputeDotProduct_ref.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/finalize.o: ../src/finalize.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/init.o: ../src/init.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/mytimer.o: ../src/mytimer.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/ComputeSPMV.o: ../src/ComputeSPMV.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/ComputeSPMV_ref.o: ../src/ComputeSPMV_ref.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/ComputeSYMGS.o: ../src/ComputeSYMGS.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/ComputeSYMGS_ref.o: ../src/ComputeSYMGS_ref.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/ComputeWAXPBY.o: ../src/ComputeWAXPBY.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/ComputeWAXPBY_ref.o: ../src/ComputeWAXPBY_ref.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/ComputeMG_ref.o: ../src/ComputeMG_ref.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/ComputeMG.o: ../src/ComputeMG.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/ComputeProlongation_ref.o: ../src/ComputeProlongation_ref.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/ComputeRestriction_ref.o: ../src/ComputeRestriction_ref.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

src/GenerateCoarseProblem.o: ../src/GenerateCoarseProblem.cpp
	$(CXX) -c $(CXXFLAGS) -I../src $< -o $@

