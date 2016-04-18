Place the private header files in this directory. They will be
available to your code with

     #include <IMP/fuzzyrestraint/internal/myheader.h>

All headers should include `IMP/fuzzyrestraint/fuzzyrestraint_config.h` as their
first include and surround all code with `IMPFUZZYRESTRAINT_BEGIN_INTERNAL_NAMESPACE`
and `IMPFUZZYRESTRAINT_END_INTERNAL_NAMESPACE` to put it in the
IMP::fuzzyrestraint::internal namespace and manage compiler warnings.
