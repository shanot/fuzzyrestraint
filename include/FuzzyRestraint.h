/**
 *  \file IMP/fuzzyrestraint/FuzzyRestraint.h
 *  \brief A restraint for ambiguous contact restraints
 *
 *  Copyright 2007-2016 IMP Inventors. All rights reserved.
 */

#ifndef IMPFUZZYRESTRAINT_FUZZY_RESTRAINT_H
#define IMPFUZZYRESTRAINT_FUZZY_RESTRAINT_H
#include "fuzzyrestraint_config.h"
#include <IMP/Restraint.h>
#include <IMP/Particle.h>
#include <Python.h>

IMPFUZZYRESTRAINT_BEGIN_NAMESPACE
//! Fuzzy Restraint
/** Fast and general implementation of the restraint provided by pmi/restraints/proteomics.
 ** Allows to combine multiple restraint into one using fuzzy logic operations.
 */

class IMPFUZZYRESTRAINTEXPORT FuzzyRestraint;
class IMPFUZZYRESTRAINTEXPORT FuzzyOr;
class IMPFUZZYRESTRAINTEXPORT FuzzyAnd;
class IMPFUZZYRESTRAINTEXPORT FuzzyRestraint : public Restraint {
public:
  FuzzyRestraint(Model *m, std::string name) : Restraint(m, name) {}
  virtual double unprotected_evaluate(DerivativeAccumulator *) const {
    return -log(get_probability());
  };
  virtual double get_probability() const = 0;
  virtual ModelObjectsTemp do_get_inputs() const = 0;
  FuzzyOr *operator|(FuzzyRestraint *r);
  FuzzyAnd *operator&(FuzzyRestraint *r);
	PyObject* get_output(){
		PyObject* out(PyDict_New());
		PyDict_SetItem(out, PyString_FromString(get_name().c_str()), PyFloat_FromDouble(get_probability()));
		return out;
	}
  IMP_OBJECT_METHODS(FuzzyRestraint);
};

class IMPFUZZYRESTRAINTEXPORT FuzzyOr : public FuzzyRestraint {
  IMP::Pointer<FuzzyRestraint> r1_;
  IMP::Pointer<FuzzyRestraint> r2_;

public:
  FuzzyOr(Model *m, FuzzyRestraint *r1, FuzzyRestraint *r2)
      : FuzzyRestraint(m, "FuzzyOr %1%"), r1_(r1), r2_(r2) {}
  virtual double get_probability() const IMP_OVERRIDE {
    Model *m(get_model());
    double const a(r1_->get_probability());
    double const b(r2_->get_probability());
    return b+a-a*b;
  }
  virtual ModelObjectsTemp do_get_inputs() const IMP_OVERRIDE {
    ModelObjectsTemp ret(r1_->get_inputs());
    ModelObjectsTemp const tmp(r2_->get_inputs());
    ret.insert(ret.end(), tmp.begin(), tmp.end());
    return ret;
  }
  IMP_OBJECT_METHODS(FuzzyOr);
};

class IMPFUZZYRESTRAINTEXPORT FuzzyAnd : public FuzzyRestraint {
  IMP::Pointer<FuzzyRestraint> r1_;
  IMP::Pointer<FuzzyRestraint> r2_;

public:
  FuzzyAnd(Model *m, FuzzyRestraint *r1, FuzzyRestraint *r2)
      : FuzzyRestraint(m, "FuzzyAnd %1%"), r1_(r1), r2_(r2) {}
  virtual double get_probability() const IMP_OVERRIDE {
    Model *m(get_model());
    double const a(r1_->get_probability());
    double const b(r2_->get_probability());
    return a * b;
  }
  virtual ModelObjectsTemp do_get_inputs() const IMP_OVERRIDE {
    ModelObjectsTemp ret(r1_->get_inputs());
    ModelObjectsTemp const tmp(r2_->get_inputs());
    ret.insert(ret.end(), tmp.begin(), tmp.end());
    return ret;
  }
  IMP_OBJECT_METHODS(FuzzyAnd);
};

class IMPFUZZYRESTRAINTEXPORT FuzzyLinearRestraint : public FuzzyRestraint {
  IMP::Pointer<Particle> p1_;
  IMP::Pointer<Particle> p2_;

  double start_;
  double end_;
  double innerslope_;

private:
  double slope_;

public:
  FuzzyLinearRestraint(Model *m, Particle *p1, Particle *p2, double start = 1,
                       double end = 10, double innerslope = 0.01)
      : FuzzyRestraint(m, "FuzzyLinearRestraint %1%"), p1_(p1), p2_(p2),
        start_(start), end_(end), innerslope_(innerslope) {
    slope_ = 1.0 / (end - start);
  }
  virtual double get_probability() const
      IMP_OVERRIDE {
    Model *m(get_model());
    double const d(core::get_distance(core::XYZ(m, p1_->get_index()),
                                      core::XYZ(m, p2_->get_index())));
    double ret(1e-14);
    if (d <= start_) {
      ret = 1.0;
    } else if (d < end_) {
      ret = 1.0 - (d - start_) * slope_;
    }
    return ret*exp(-d * innerslope_);
  }
  virtual ModelObjectsTemp do_get_inputs() const IMP_OVERRIDE {
    ModelObjectsTemp ret;
    Model *m(get_model());
    ret.push_back(m->get_particle(p1_->get_index()));
    ret.push_back(m->get_particle(p2_->get_index()));
    return ret;
  }
  IMP_OBJECT_METHODS(FuzzyLinearRestraint);
};

class IMPFUZZYRESTRAINTEXPORT FuzzySigmoidRestraint : public FuzzyRestraint {
  IMP::Pointer<Particle> p1_;
  IMP::Pointer<Particle> p2_;

  double theta_;
  double slope_;
  double plateau_;
  double innerslope_;

public:
  FuzzySigmoidRestraint(Model *m, Particle *p1, Particle *p2,
                        double theta = 5.0, double slope = 2.0,
                        double plateau = 1e-14, double innerslope = 0.01)
      : FuzzyRestraint(m, "FuzzySigmoidRestraint %1%"), p1_(p1), p2_(p2),
        theta_(theta), slope_(slope), plateau_(plateau),
        innerslope_(innerslope) {}
  virtual double get_probability() const IMP_OVERRIDE {
    Model *m(get_model());
    double const d(core::get_distance(core::XYZ(m, p1_->get_index()),
                                      core::XYZ(m, p2_->get_index())));
    double const argvalue((d - theta_) / slope_);
    return ((1.0 - (1.0 - plateau_) / (1.0 + exp(-argvalue)))*exp(-innerslope_ * d));
  };
  virtual ModelObjectsTemp do_get_inputs() const IMP_OVERRIDE {
    ModelObjectsTemp ret;
    Model *m(get_model());
    ret.push_back(m->get_particle(p1_->get_index()));
    ret.push_back(m->get_particle(p2_->get_index()));
    return ret;
  }
  IMP_OBJECT_METHODS(FuzzySigmoidRestraint);
};

FuzzyOr *FuzzyRestraint::operator|(FuzzyRestraint *r) {
  IMP_NEW(FuzzyOr, ret, (get_model(), this, r));
  return ret.release();
}
FuzzyAnd *FuzzyRestraint::operator&(FuzzyRestraint *r) {
  IMP_NEW(FuzzyAnd, ret, (get_model(), this, r));
  return ret.release();
}
IMPFUZZYRESTRAINT_END_NAMESPACE
#endif /* IMPFUZZYRESTRAINT_FUZZY_RESTRAINT_H */
