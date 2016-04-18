Place the public header files in this directory. They will be
available to your code (and other modules) with

     #include <IMP/fuzzyrestraint/myheader.h>

All headers should include `IMP/fuzzyrestraint/fuzzyrestraint_config.h` as their
first include and surround all code with `IMPFUZZYRESTRAINT_BEGIN_NAMESPACE`
and `IMPFUZZYRESTRAINT_END_NAMESPACE` to put it in the IMP::fuzzyrestraint namespace
and manage compiler warnings.

Headers should also be exposed to SWIG in the `pyext/swig.i-in` file.
