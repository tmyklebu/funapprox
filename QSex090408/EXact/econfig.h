#ifndef __econfig__
#define __econfig__
#include "../EG/eg_config.h"
#include "../EG/eg_macros.h"
#include "../EG/eg_nummacros.h"
#include "../EG/eg_mem.h"
#include "../EG/eg_lpnum.h"
#include "../EG/eg_lpnum.dbl.h"
#include "../EG/eg_lpnum.int.h"
#include "../EG/dbl_eg_numutil.h"
#include "../EG/mpf_eg_numutil.h"
#include "../EG/mpq_eg_numutil.h"
#include "../EG/eg_lpnum.mpf.h"
#include "../EG/eg_lpnum.mpq.h"
#include "../EG/eg_lpnum.mpz.h"
#ifdef USE_SUNOS
extern int strcasecmp();
extern int strncasecmp();
extern int snprintf(char *, size_t, const char *, ...);
#endif
#endif
