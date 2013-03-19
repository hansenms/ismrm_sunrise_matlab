/*
* def,inline.h
*
* Simple vector operations done in-line.
*
* Copyright 1994-12, Jeff Fessler, The University of Michigan
*/

#define	whileZn	for (; Zn; --Zn)

/*
* v op c(onstant)
*/
#define Inline1vector(type1, v1, op, type2, c, n)	\
{				\
register type1	*Z1 = v1;	\
register type2	Zc = (type2) c;	\
register int	Zn = n;		\
	whileZn			\
		*Z1++ op Zc;	\
}

/*
* v1 op v2
*/
#define Inline2vector(type1, v1, op, type2, v2, n)	\
{				\
register type1	*Z1 = v1;	\
register type2	*Z2 = v2;	\
register int	Zn = n;		\
	whileZn {		\
		*Z1++ op (*Z2);	\
	++Z2; }			\
}




/*
* Vector with Scalar
*/

/*
* v *= c(onstant)
*/
#define VectScale(type, v, c, n)	\
	Inline1vector(type, v, *=, cdouble, c, n)

/*
* v += c(onstant)
*/
#define VectInc(type, v, c, n)	\
	Inline1vector(type, v, +=, cdouble, c, n)

/*
* v = c(onstant)
*/
#define VectSet(type, v, c, n)	\
	Inline1vector(type, v, =, Const type, c, n)

/*
* v = Min(v,c)
*/
#define VectMin(type, v, c, n)	\
{				\
register type	*Zv = v;	\
register int	Zn = n;		\
	whileZn {		\
		if (*Zv > c) *Zv = c;	\
		++Zv; }		\
}


/*
* v = Max(v,c)
*/
#define VectMax(type, v, c, n)	\
{				\
register type	*Zv = v;	\
register int	Zn = n;		\
	whileZn {		\
		if (*Zv < c) *Zv = c;	\
		++Zv; }		\
}

/*
* v = Max(v,0)
*/
#define VectNonneg(type, v, n)	VectMax(type, v, 0, n)


/*
* s = sum(v)
*/
#define VectAccum(type, s, v, n)	\
{					\
register double	Zsum = 0.;		\
register Const type *Zv = v;		\
register int	Zn = n;			\
	whileZn				\
		Zsum += *(Zv++);	\
	s = Zsum;			\
}


/*
* s = sum(|v|)
*/
#define VectNorm1(type, s, v, n)	\
{					\
register Const type *Zv = v;		\
register int	Zn = n;			\
register double	Zsum = 0;		\
	whileZn {			\
		Zsum += Abs(*Zv);	\
		++Zv;			\
	} s = Zsum;			\
}


/*
* s = <v,v>
*/
#define VectNorm2(type, s, v, n)	\
{					\
register Const type *Zv = v;	\
register int	Zn = n;			\
register double	Zsum = 0;		\
	whileZn {			\
		Zsum += *Zv * *Zv;	\
		++Zv;			\
	} s = Zsum;			\
}


/*
* s = <v,v>_w
*/
#define VectNorm2wtd(type, s, v, w, n)	\
{					\
register Const type *Zv = v;	\
register Const type *Zw = w;	\
register int	Zn = n;				\
register double	Zsum = 0;			\
	whileZn {				\
		Zsum += *Zv * *Zv * *Zw;	\
		++Zv; ++Zw;			\
	} s = Zsum;				\
}


/*
* [v; v; ...; v] repeated nrep times
*/
#define VectRepeat(type, v, npoint, nrep)	\
{					\
register type *Zv = v;			\
Const int	Zn = npoint;		\
Const int	Zr = nrep;		\
int	Zi;				\
	for (Zi=1; Zi < Zr; Zi++)	\
		Bcopy(Zv, Zv+Zi*Zn, Zn)	\
}


/*
* Vector with Vector
*/

/*
* v1 = v2
*/
#define VectCopy(type1, v1, type2, v2, n)	\
	Inline2vector(type1, v1, =, Const type2, v2, n)

/*
* r *= v
*/
#define VectMul(type, r, v, n)	\
	Inline2vector(type, r, *=, Const type, v, n)

/*
* r /= v
*/
#define VectDiv(type, r, v, n)	\
	Inline2vector(type, r, /=, Const type, v, n)

/*
* r += v
*/
#define VectAdd(type, r, v, n)	\
	Inline2vector(type, r, +=, Const type, v, n)

/*
* r += const * v
*/
#define VectAddScale(type, r, c, v, n)	\
	Inline2vector(type, r, += c *, Const type, v, n)

/*
* r -= v
*/
#define VectSub(type, r, v, n)	\
	Inline2vector(type, r, -=, Const type, v, n)

/*
* r += v1 + v2
*/
#define VectInc2(type, r, v1, v2, n)	\
{					\
register type	*Zr = r;		\
register Const type *Z1 = v1;		\
register Const type *Z2 = v2;		\
register int	Zn = n;			\
	whileZn				\
		*Zr++ += *(Z1++) + *(Z2++);	\
}

/*
* r = r / v or 0 if v=0
*/
#define VectDiv0(type, r, v, n)	\
{				\
register type *Zr = r;		\
register Const type *Zv = v;	\
register int	Zn = n;		\
	whileZn {		\
		if (*Zv)	\
			*Zr++ /= *Zv;	\
		else {			\
			*Zr++ = 0;	\
		}			\
		++Zv; }		\
}


/*
* s = <v1, v2>
*/
#define VectInprod3d(s, type1, v1, type2, v2, n1, n2, n3)	\
{					\
register Const type1 *Z1 = v1;		\
register Const type2 *Z2 = v2;		\
int	Zn = n3;		\
cint	Zn2 = n2;		\
register cint	Zn1 = n1;	\
double	Zsum = 0.;		\
	whileZn	{		\
		int Zi2;		\
		double	Zsum2 = 0.;	\
		for (Zi2=0; Zi2 < Zn2; ++Zi2) {			\
			register int Zi1;			\
			register double	Zsum1 = 0.;		\
			for (Zi1=0; Zi1 < Zn1; ++Zi1)		\
				Zsum1 += *(Z1++) * *(Z2++);	\
			Zsum2 += Zsum1;				\
		}		\
		Zsum += Zsum2;	\
	}			\
	s = Zsum;		\
}

#define VectInprod2d(s, type1, v1, type2, v2, n1, n2)	\
	VectInprod3d(s, type1, v1, type2, v2, n1, n2, 1)

/*
* s = <v1, v2>
*/
#define VectInprod(type, s, v1, v2, n)	\
{					\
register Const type *Z1 = v1;		\
register Const type *Z2 = v2;		\
register int	Zn = n;			\
register double	Zsum = 0.;		\
	whileZn				\
		Zsum += *(Z1++) * *(Z2++);	\
	s = Zsum;				\
}


/*
* s = <v1, v2.^2>
*/
#define VectInprod2(type, s, v1, v2, n)	\
{					\
register Const type *Z1 = v1;		\
register Const type *Z2 = v2;		\
register int	Zn = n;			\
register double	Zsum = 0.;		\
	whileZn {			\
		Zsum += *(Z1++) * *Z2 * *Z2;	\
		++Z2;			\
	} s = Zsum;			\
}


/*
* s = <v1, v2.^2>
* Supposes v2 is mostly zeros
*/
#define VectInprod2_0(type, s, v1, v2, n)	\
{					\
register Const type *Z1 = v1;		\
register Const type *Z2 = v2;		\
register int	Zn = n;			\
register double	Zsum = 0.;		\
	whileZn {			\
		if (*Z2)		\
			Zsum += *Z1 * *Z2 * *Z2;	\
		++Z1; ++Z2;		\
	} s = Zsum;			\
}


/*
* min and max of vector
*/
#define VectMinMax(type1, vmin, vmax, type2, v, n)	\
{						\
register Const type2	*Zv = v;		\
register type1	Zmin = (type1) *Zv;		\
register type1	Zmax = (type1) *Zv++;		\
register int	Zn = n-1;			\
	whileZn {				\
		if (*Zv > Zmax)	Zmax = *Zv;	\
		else				\
		if (*Zv < Zmin)	Zmin = *Zv;	\
		++Zv;				\
	}					\
	vmin = Zmin;				\
	vmax = Zmax;				\
}


/*
* min and max of vector, and min > 0
*/
#define VectMinMax0(type1, this_type_max, vmin, vmin0, vmax, type2, v, n)	\
{						\
register Const type2	*Zv = v;		\
register type1	Zmin = (type1) *Zv;		\
register type1	Zmin0 = (type1) (Zmin > 0 ? Zmin : this_type_max);	\
register type1	Zmax = (type1) *Zv++;			\
register int	Zn = n-1;			\
	whileZn {				\
		if (*Zv > Zmax)	Zmax = *Zv;	\
		else				\
		if (*Zv < Zmin)	Zmin = *Zv;	\
		if (*Zv > 0 && *Zv < Zmin0)	Zmin0 = *Zv;	\
		++Zv;				\
	}					\
	vmin = Zmin;				\
	vmin0 = Zmin0;				\
	vmax = Zmax;				\
}


/*
* debugging
*/

#define VectInfo(type, v, n, msg) { \
	double Ztmin, Ztmax, Ztmin0, ZZsum; \
	VectMinMax0(type, MAX##type, Ztmin, Ztmin0, Ztmax, type, v, n) \
	VectAccum(type, ZZsum, v, n) \
	Note6("%s [%d] min=%g min0=%g max=%g sum=%g", \
		msg, n, Ztmin, Ztmin0, Ztmax, ZZsum) }
