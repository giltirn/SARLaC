#ifndef _GLOBALFIT_DEF_BINDING_H_
#define _GLOBALFIT_DEF_BINDING_H_

#define _GLOBALFIT_BIND_SUB_TO_SUPER(r, dummy, elem) out.BOOST_PP_TUPLE_ELEM(0,elem) = in.BOOST_PP_TUPLE_ELEM(1,elem);
#define _GLOBALFIT_BIND_SUPER_TO_SUB(r, dummy, elem) out.BOOST_PP_TUPLE_ELEM(1,elem) = in.BOOST_PP_TUPLE_ELEM(0,elem);

//Automatically generate the binding structure for a sub-fitfunc
#define DEF_GLOBALFIT_BINDING(BINDING_NAME, SUBSET_TYPE, SUPERSET_TYPE, SUPERSET_SUBSET_BINDINGS)	\
  struct BINDING_NAME{							\
  typedef SUBSET_TYPE subsetType;					\
  typedef SUPERSET_TYPE supersetType;					\
									\
  static void supersetToSubset(subsetType &out, const supersetType &in){ \
    TUPLE_SEQUENCE_FOR_EACH(_GLOBALFIT_BIND_SUPER_TO_SUB, , SUPERSET_SUBSET_BINDINGS); \
  }									\
  static void subsetToSuperset(supersetType &out, const subsetType &in){ \
    TUPLE_SEQUENCE_FOR_EACH(_GLOBALFIT_BIND_SUB_TO_SUPER, , SUPERSET_SUBSET_BINDINGS); \
  }									\
  };

#endif
