/*
* defs-env.h
* Portability definitions
*
* Copyright Dec. 1994, Jeff Fessler, The University of Michigan
*/

#ifndef DefsEnv
#define DefsEnv

/*
* set the following flag for compiling on PC/NT/x86
*/
#if defined(Is_pc)
#define	Userand
#define	Need_uint
#define Need_proto_lgamma
#define Need_proto_erf
#define Need_time_h
#endif

#ifdef Use_thread
#include "def,thread.h"
#endif

#include <stdio.h>
#include <string.h>
#if defined(Need_assert)
#include <assert.h>
#endif
#include <sys/types.h>
#include <ctype.h>
#if defined(__osf__) | defined(Need_time_h)
#include <time.h>
#endif
#if defined(Need_ieeefp)
#include <ieeefp.h>	/* needed for solaris for finite() */
#endif
#if defined(Need_dynlib)
#include <dlfcn.h>	/* needed for dynamic libraries	*/
#endif

#ifdef Use_mpi
#include <mpi.h>
#endif

#if defined(hphp)
#define _INCLUDE_XOPEN_SOURCE
#define Userand48
	extern double drand48(void);
#if !defined(Mmex)
#	define Need_uint
#endif
#endif

#if defined(Need_uint)
	typedef unsigned int	uint;
	typedef unsigned short	ushort;
#endif

#if defined(titan)
#define Isfinite isfinite	/* no finite() on titan */
#define Provide_strstr		/* no strstr() on titan */
#define Need_proto_strstr
#include <vmath.h>
#else
#include <math.h>
#endif

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2	0.70710678118654752440
#endif
#ifndef M_SQRT2
#define M_SQRT2		1.4142135623730950488
#endif

#ifndef CLOCKS_PER_SEC		/* stupid SunOS doesn't have it */
#	define CLOCKS_PER_SEC 1000000
#endif

#ifndef Clock
	extern clock_t GlobalClock0;
#define	Clock0	GlobalClock0 = clock();
#define	Clock	( (clock() - GlobalClock0) / (double) CLOCKS_PER_SEC )
#endif

/*
* No stdlib on titan!
*/
#if defined(titan)
	extern void	exit();
	extern char	*calloc(), *malloc();
	extern char	*memcpy(), *memset();
	extern int	free();
#else
#	include <stdlib.h>
#endif

/*
* Other externs
*/
#ifdef Need_protos_sol2		/* needed for sol2 */
#	define Need_proto_erf
#	define Need_proto_lgamma
#	define Need_proto_random
#endif

#if defined(titan) | defined(sun)
	extern int	bzero(), bcopy();
	/*#	define Free(p)	free((char *) (p)) */
#endif


/*
* And of course this varies too
*/
#if defined(_SIZE_T_)
#	define Size_t	size_t
#else
#	define Size_t	unsigned
#endif

#if !defined(noConst)
#	define Const	const
#else
#	define Const
#endif

#if !defined(MAXint)
#	define MAXchar		127
#	define MAXbyte		255
#	define MAXshort		32767
#	define MAXint		2147483647
#	define MAXfloat		1e30
#	define MAXdouble	1e300
#endif

/*
* others
*/
#include "def,macro.h"

#ifndef No_defs_type
#include "def,type.h"
#endif

#include "def,inline.h"
#include "def,alloc.h"
#include "def,proto.h"

#ifndef SHRT_MAX
#define SHRT_MAX	32767
#endif
#ifndef SHRT_MIN
#define SHRT_MIN	(-32768)
#endif

/*
* define all the globals that are needed by the macros
*/
#ifndef GlobalsDefinedHere
#	define GlobalsDefinedHere	\
		FILE	*Global_fp_exist;	\
		clock_t GlobalClock0;	\
		int Mpi_rank=0, Mpi_numprocs=1;
#endif

#endif /* DefsEnv */
