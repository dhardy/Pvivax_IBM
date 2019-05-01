/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///                                                                       ///
///  Individual-based Plasmodium vivax transmission model.                ///
///                                                                       ///
///  Dr Michael White                                                     ///
///  Institut Pasteur                                                     ///
///  michael.white@pasteur.fr                                             ///
///                                                                       ///
///                                                                       ///
///  Random number generator support code                                 ///
///                                                                       ///
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


#include "randlib.h"

// Support function to sample a number between 0 and 1
inline double gen_u01() {
    return genunf(0, 1);
}

// Support function to sample a boolean with given probability
inline bool gen_bool(double p) {
    return gen_u01() < p;
}
