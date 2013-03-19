/*
* def,macro.h
* favorite macro definitions
*
* Copyright 1994-12, Jeff Fessler, The University of Michigan
*/

#ifndef Defs_macro
#define Defs_macro

#define	Chat		((chat > 0) ? (chat-1) : 0)

#define	Calloc(n,s)	calloc((unsigned) (n), (unsigned) (s))

#define	Swab(s,d,n)	{ \
	if ((s) == (d)) Exit("swab problem") \
	swab((char *)(s), (char *)(d), n); }

#define	Nonull(str)	((str) ? (str) : "null")

#define Argflt(n,def)	(argc > n ? atof(argv[n]) : def)
#define Argint(n,def)	(argc > n ? atoi(argv[n]) : def)
#define Argstr(n,def)	(argc > n ? argv[n] : def)
#define Argstn(n)	((argc > n && strcmp(argv[n], "-")) ? argv[n] : NULL)

#define Isdigit(c)	(('0' <= (c)) && ((c) <= '9'))

/* use -DIsfinite=isfinite on lnx86 */
#if !defined(Isfinite)
#define Isfinite(x)	finite(x)
#endif


#ifndef Notes
#	define Notes(msg)	{ if (chat) { \
		printf("Note %s %d: ", __FILE__, __LINE__); \
		printf msg; (void)fflush(stdout); } }
#endif

#ifndef Note
#	define Note(msg)	\
	{ (void)fprintf(stdout, "Note %s %d: %s\n", __FILE__, __LINE__, msg); \
	(void)fflush(stdout); }
#endif

#ifndef Note1
#	define Note1(msg, arg)	\
	{ char Zstr[1024]; (void)sprintf(Zstr, msg, arg); Note(Zstr) }
#endif

#ifndef Note2
#	define Note2(msg, a1, a2)	\
	{ char Zstr[1024]; (void)sprintf(Zstr, msg, a1, a2); Note(Zstr) }
#endif

#ifndef Note3
#	define Note3(msg, a1, a2, a3)	\
	{ char Zstr[1024]; (void)sprintf(Zstr, msg, a1, a2, a3); Note(Zstr) }
#endif

#ifndef Note4
#	define Note4(msg, a1, a2, a3, a4)	\
	{ char Zstr[1024]; (void)sprintf(Zstr, msg, a1, a2, a3, a4); Note(Zstr) }
#endif

#ifndef Note5
#	define Note5(msg, a1, a2, a3, a4, a5)	\
	{ char Zstr[1024]; (void)sprintf(Zstr, msg, a1, a2, a3, a4, a5); Note(Zstr) }
#endif

#ifndef Note6
#	define Note6(msg, a1, a2, a3, a4, a5, a6)	\
	{ char Zstr[1024]; (void)sprintf(Zstr, msg, a1, a2, a3, a4, a5, a6); Note(Zstr) }
#endif

#ifndef Warn
#	define Warn(msg)	{ \
	(void)fprintf(stderr, "WARN %s %d: %s\n", __FILE__, __LINE__, msg); \
	(void)fflush(stderr); }
#endif

#ifndef Warn1
#	define Warn1(msg, arg)	\
	{ char Zstr[1024]; (void)sprintf(Zstr, msg, arg); Warn(Zstr) }
#endif

#ifndef Warn2
#	define Warn2(msg, a1, a2)	\
	{ char Zstr[1024]; (void)sprintf(Zstr, msg, a1, a2); Warn(Zstr) }
#endif

#ifndef Warn3
#	define Warn3(msg, a1, a2, a3)	\
	{ char Zstr[1024]; (void)sprintf(Zstr, msg, a1, a2, a3); Warn(Zstr) }
#endif

#ifndef Warn4
#	define Warn4(msg, a1, a2, a3, a4)	\
	{ char Zstr[1024]; (void)sprintf(Zstr, msg, a1, a2, a3, a4); Warn(Zstr) }
#endif

#ifndef Fail
#	define Fail(msg)	{ \
	(void)fprintf(stderr, "FAIL %s %d: %s\n", __FILE__, __LINE__, msg); \
	(void)fflush(stderr); \
	return Failure; }
#endif

#ifndef Fail1
#	define Fail1(msg, arg)	\
	{ char Zstr[1024]; (void)sprintf(Zstr, msg, arg); Fail(Zstr) }
#endif

#ifndef Fail2
#	define Fail2(msg, a1, a2)	\
	{ char Zstr[1024]; (void)sprintf(Zstr, msg, a1, a2); Fail(Zstr) }
#endif

#ifndef Fail3
#	define Fail3(msg, a1, a2, a3)	\
	{ char Zstr[1024]; (void)sprintf(Zstr, msg, a1, a2, a3); Fail(Zstr) }
#endif

#ifndef Fail4
#	define Fail4(msg, a1, a2, a3, a4)	\
	{ char Zstr[1024]; (void)sprintf(Zstr, msg, a1, a2, a3, a4); Fail(Zstr) }
#endif

#ifndef Fail5
#	define Fail5(msg, a1, a2, a3, a4, a5)	\
	{ char Zstr[1024]; (void)sprintf(Zstr, msg, a1, a2, a3, a4, a5); Fail(Zstr) }
#endif

#ifndef Fail6
#	define Fail6(msg, a1, a2, a3, a4, a5, a6)	\
	{ char Zstr[1024]; (void)sprintf(Zstr, msg, a1, a2, a3, a4, a5, a6); Fail(Zstr) }
#endif

#ifndef Exit
#	define Exit(msg)	{ \
	(void)fprintf(stderr, "EXIT %s %d: %s\n", __FILE__, __LINE__, msg); \
	exit(-1); }
#endif

#ifndef Call
#define Call(fun, arg)		{ if (!(fun arg)) Fail(#fun) }
#define Call1(fun, arg, str)	{ if (!(fun arg)) Fail1(#fun" %s", str) }
#endif

/*
* And those pesky IO functions....
*/
#if defined(titan)
#	define Fread(p,s,n,f)	fread ((char *)(p), (Size_t)(s), (int)(n), f)
#	define Fwrite(p,s,n,f)	fwrite((char *)(p), (Size_t)(s), (int)(n), f)
#elif defined(sun)
#	define Fread(p,s,n,f)	fread ((char *)(p), (int)(s), (int)(n), f)
#	define Fwrite(p,s,n,f)	fwrite((char *)(p), (int)(s), (int)(n), f)
#else
#	define Fread(p,s,n,f)	fread ((void *)(p), (Size_t)(s), (Size_t)(n), f)
#	define Fwrite(p,s,n,f)	fwrite((cvoid *)(p), (Size_t)(s), (Size_t)(n), f)
#	define Strncat(s1,s2,n)	strncat(s1, s2, (int) (n))
#endif

#if defined(titan) | defined(sun)
#define Memset(p,c,n)		memset((char *)(p), (int)(c), (int)(n));
#define Memcpy(t,f,n)		memcpy((char *)(t), (char *)(f), (int)(n));
#define Strncat(s1,s2,n)	strncat(s1, s2, (Size_t) (n))
#define Bzero(p,n)	bzero((char *) (p), (int) ((n) * sizeof(*(p))));
#define Bcopy(s,d,n)	bcopy((char *) (s), (char *) (d), (int) ((n) * sizeof(*(s))));
#else
#define Memset(p,v,n)	(void)memset((void *) (p), v, (Size_t) (n));
#define Memcpy(d,s,n)	{ if ((d) == (s)) Exit("Identical memcpy pointers") \
		else (void)memcpy((void *) (d), (cvoid *) (s), (Size_t) (n)); }
#define Bzero(p,n)	Memset(p, 0, (n) * sizeof(*(p)))
#define Bcopy(s,d,n)	Memcpy(d, s, (n) * sizeof(*(s)))
#endif

#if 0 /* superceded by versions in io library */
#ifndef FileExist
	/* (!fclose(fopen(name, "rb"))) did not work on linux */
	extern FILE *Global_fp_exist;
#define FileExist(name)		((Global_fp_exist = fopen(name, "rb"), \
				Global_fp_exist && !fclose(Global_fp_exist)))
#endif

/*
* macro to check to see if we can write a file to this name.
* use this with caution since if the file exists, it will be truncated.
*/
#ifndef FileWriteable
#define FileWriteable(s)	((Global_fp_exist = fopen(s, "wb"), \
				Global_fp_exist && fclose(Global_fp_exist))) 
#endif
#endif

/*
* The suffix "0" means return 0 on error.
*/
#define	Fopen0(fp,n,t)		{if ( !(fp = fopen(n, t)) )	\
					Fail1("fopen '%s'", n)}
#define	Fopen0r(fp,n)		Fopen0(fp,n,"rb")
#define	Fopen0w(fp,n)		Fopen0(fp,n,"wb")
#define	Fopen0e(fp,n)		Fopen0(fp,n,"rb+")	/* for editing */
#define	Fflush0(fp)		{if (fflush(fp))	Fail("fflush")}
/* seek from start of file */
#define	Fseek0(fp,n)		{if (fseek(fp, (long) (n), 0))	\
					Fail1("fseek %d", n)}
#define	Fskip0(fp,n)		{if (fseek(fp, (long) (n), 1))	\
					Fail1("fseek1 %d", n)}
#define	Fclose0(fp)		{if (fclose(fp))	Fail("fclose")}
#define	Fread0(p,s,n,fp)	{uint ZF = Fread(p, s, n, fp); \
				if (ZF != (unsigned) (n)) \
				Fail3("fread found %d bytes, not %d at %d", ZF, (n), (int) ftell(fp))}
#define	Fwrite0(p,s,n,fp)	{if (((unsigned) (n)) != (unsigned) Fwrite(p, s, n, fp)) Fail("fwrite")}


#if defined(Userand48)
#	define Rand01		drand48()
#	define Seeder(s)	srand48(s)
#elif defined(Userand)				/* for PC */
#	define Rand01		(rand() / RAND_MAX)
#	define Seeder(s)	srand(s)
#else
#	define Rand01		(random() / 2147483647.0)
#	define Seeder(s)	srandom(s)
#endif

/*
* The remaining macro definitions should be environment independent
*/

#define Sqr(x)		((x) * (x))
#define Odd(n)		((n)/2 != ((n)+1)/2)
#define Even(n)		((n)/2 == ((n)+1)/2)
#ifndef Abs
#	define Abs(x)	(((x) >= 0) ? (x) : (-x))
#endif
#ifndef Max
#	define Max(x,y)	(((x) > (y)) ? (x) : (y))
#endif
#ifndef Min
#	define Min(x,y)	(((x) < (y)) ? (x) : (y))
#endif
#ifndef Rint	/* round to nearest integer, ala rint() */
#	define Rint(x)	( (int) ((x) >= 0 ? ((x) + 0.5) : ((x) - 0.5)) )
	/* formerly ((int) ((x) + 0.5)), which is wrong for negatives! */
#endif

#define Floor(x)	((int) floor(x))
#define Ceil(x)		((int) ceil(x))

#ifndef True
#	define	True	((bool) 1)
#endif
#ifndef False
#	define	False	((bool) 0)
#endif

#define Streq(a,b)	(!strcmp(a,b))
#define Streqn(a,b,n)	(!strncmp(a,b,n))

#define	Success	((bool) 1)
#define	Failure	0
#define	Ok	return Success;

#ifdef CountAlloc

#define	Alloc0(ptr,type,n,s)	{ \
	if ( !((ptr) = (type *) io_mem_alloc((uint) (n), (uint) (s), \
		__FILE__, __LINE__, #ptr)) ) \
		Fail2("io_mem_alloc(%d,%d)", (int) (n), (int) (s)) }
#define Free(p)		{ \
		if (!io_mem_free((void *) (p), __FILE__, __LINE__, #p)) { \
			Warn2("bug free %s %d", __FILE__, __LINE__) \
			Exit("bug") } }
#define Free0(p)	{ \
		Call(io_mem_free, ((void *) (p), __FILE__, __LINE__, #p)) \
		p = NULL; }

#define PrintAlloc	io_mem_print(__FILE__, __LINE__, 1);
#define MemUsed		io_mem_usage(__FILE__, __LINE__, Chat);
#define MemInfo(p)	io_mem_info((cvoid *) p, __FILE__, __LINE__, #p);

#else

#define	Alloc0(ptr,type,n,s)	{ \
	if ( !((ptr) = (type *) calloc((unsigned) (n), (unsigned) (s))) ) \
		Fail2("alloc %d %d", (int) (n), (int) (s)) }
#define Free(p)		free((void *) (p))
#define Free0(p)	Free(p);

#define PrintAlloc	{}
#define MemUsed		{}
#define MemInfo(p)	{}

#endif

#define Mem0(p,n)	Alloc0(p, void, n, sizeof(*(p)))

#endif /* Defs_macro */
