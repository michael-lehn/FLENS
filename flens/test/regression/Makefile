#
#  Regression tests.  Each test is a self-checking program that returns a
#  non-zero exit status on failure.
#
#    make check                        # build and run all tests
#    make check CXX=clang++ CXXSTD=c++20     (see toolchain.mk)
#    make check REFLAPACK="-L/opt/homebrew/opt/openblas/lib -lopenblas"
#
#  With REFLAPACK set, results are additionally compared with an external
#  (reference) LAPACK, and the tests of the Fortran interface are built.
#

FLENS_ROOT ?= ../../..
include $(FLENS_ROOT)/toolchain.mk

BUILD     ?= build
REFLAPACK ?=

TESTS         = sterf tiny-copy tiny-gemv lantr matrixmarket lsame-interface
TESTS_LAPACK  = zheev-interface lqf-ilaenv gesvj-interface rot-interface

ifneq ($(strip $(REFLAPACK)),)
    TESTS    += $(TESTS_LAPACK)
    REFFLAGS  = -DWITH_REFERENCE_LAPACK
endif

ALL_CXXFLAGS = -std=$(CXXSTD) $(CXXFLAGS) -I$(FLENS_ROOT) $(REFFLAGS) -MMD -MP
PROGRAMS     = $(addprefix $(BUILD)/,$(TESTS))

all : $(PROGRAMS)

check : $(PROGRAMS)
	@failed=0; \
	for t in $(TESTS); do \
	    if (cd $(BUILD) && ./$$t > $$t.log 2>&1); then echo "PASS: $$t"; \
	    else echo "FAIL: $$t (see $(BUILD)/$$t.log)"; \
	         failed=$$((failed+1)); fi; \
	done; \
	if [ -z "$(strip $(REFLAPACK))" ]; then \
	    echo "note: REFLAPACK not set, reference comparisons skipped"; \
	fi; \
	test $$failed -eq 0

$(BUILD)/% : %.cc | $(BUILD)
	$(CXX) $(ALL_CXXFLAGS) -o $@ $< $(REFLAPACK)

$(BUILD) :
	mkdir -p $@

clean :
	$(RM) -r $(BUILD)

.PHONY : all check clean

-include $(PROGRAMS:%=%.d)
