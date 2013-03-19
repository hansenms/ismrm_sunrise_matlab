// cbct,def.h
// cone-beam CT definitions
// Copyright 2008-10-09, Jeff Fessler, University of Michigan

#ifndef jf_cbct_def_h
#define jf_cbct_def_h

#include "defs-env.h"

// initial version used "double" for scalars
#ifndef rscalar
#define rscalar float
#endif
#define crscalar Const rscalar

// image geometry
typedef struct {
	int nx; // image dimensions
	int ny;
	int nz;
	rscalar dx; // pixel size (can be negative)
	rscalar dy; // can be negative to cause flip
	rscalar dz; // cannot be negative
	rscalar offset_x; // center offset in pixels (usually 0)
	rscalar offset_y;
	rscalar offset_z;
	byte *mask2; // [nx ny] 2D support mask: 0 or 1 ... nthread
} cbct_ig;

// cone-beam geometry
typedef struct {
	rscalar dso; // distance from source to isocenter (infinity for ||)
	rscalar dsd; // distance from source to detector (infinity for ||)
	rscalar dfs; // distance from detector focal point to source
			// 0 for 3rd-gen CT, infinity for flat detector
	int ns; // # detector channels per row
	int nt; // # detector rows (along axial direction)
	rscalar ds; // horizontal (transaxial) ray spacing
	rscalar dt; // vertical (axial) ray spacing
	rscalar offset_s; // detector channel offset [pixels]
	rscalar offset_t; // detector vertical offset [pixels]
} cbct_cg;

#define cbct_view_pad 3 // extra space to support 4-wide simd

// work space (one for each thread)
typedef struct {
	float *view; // ([ns nt] or [nt ns]) + cbct_view_pad

	// for sf1:
	float *gamma; // [ns] azimuthal angles between each ray and center one
	float *vec; // [Max(ns,nt)]
	float *weight_s; // [ns]
#if 0
	float *uppers; // [nx+1]
#endif
	// for sf4 and sf5 only
	float dl; // [1], [0 fovz/2), dividing line for rectangle and trapezoid approximation in t
	// for sf2 only
	float *weight_inner; // [nt] backprojector: sum_k(weight_s(k) * projection_t(k)) for each it
				// forward projector: sum_k(weight_t(k) * image(k)) for each it
} cbct_work;


// type of CT geometry
typedef enum {
	cbct_geom_unknown,
	cbct_geom_par, // parallel beam (dso = infinity)
	cbct_geom_flat, // flat detector (dfs = infinity)
	cbct_geom_arc3 // dfs = 0
} cbct_geom;


// cbct,mask2.c

extern sof cbct_mask_init(
byte *mask_int, // [nx ny] out: integers 0 or 1 ... nthread
cbyte *mask_bin, // [nx ny] in: binary, possibly null
cint nx, cint ny,
cint nthread,
cint chat);

typedef enum {
	cbct_back_error,
	cbct_back_zero, // zero before back-projection
	cbct_back_inc, // incrementing back-projection (no zeroing first)
} cbct_back_init;


// cbct,*,back.c

typedef sof cbct_any_back1_type(
float *image, // [nz nx ny] <- trick!
cint nx,
cint ny,
cint nz,
crscalar dx,
crscalar dy, // can be negative to cause flip
crscalar dz,
crscalar offset_x, // center offset in pixels (usually 0)
crscalar offset_y,
crscalar offset_z_shift, // offset_z + source_z / dz
cbyte *mask2, // [nx ny] 2D support mask: 0, 1, ..., nthread
cbyte mask_id, // 1 ... nthread
crscalar dso, // distance from source to isocenter
crscalar dsd, // distance from source to detector
crscalar dfs, // distance from focal point to source (0 or inf)
cint ns,
cint nt,
crscalar ds, // horizontal ray spacing
crscalar dt, // vertical ray spacing
crscalar beta, // source angle [radians]
crscalar offset_s, // channel offset [pixels]
crscalar offset_t, // vertical offset on detector [pixels]
cfloat *proj_in, // [ns*nt] <- trick! projection view at angle beta
ctruf i_is_ns_nt, // 1 if input is [ns nt] or 0 if [nt ns]
cbct_work *cw, // work space
cfloat scale, // scale input projection view by this factor before backprojecting
cint iz_start, // usually 0
cint iz_end); // do [iz_start, iz_end)

extern cbct_any_back1_type cbct_nn1_back1; // cbct,nn1,back.c
extern cbct_any_back1_type cbct_pd1_back1; // cbct,pd1,back.c
extern cbct_any_back1_type cbct_sf1_back1; // cbct,sf1,back.c
extern cbct_any_back1_type cbct_sf2_back1; // cbct,sf2,back.c
extern cbct_any_back1_type cbct_sf3_back1; // cbct,sf3,back.c
extern cbct_any_back1_type cbct_sf4_back1; // cbct,sf4,back.c
extern cbct_any_back1_type cbct_sf5_back1; // cbct,sf5,back.c

// cbct,any,back,t.c

extern sof cbct_any_back_t(
float *image, // [nz nx ny] <- trick!
const cbct_ig *ig,
const cbct_cg *cg,
cbct_work *cw,
cint na, // # of views
cfloat *proj, // [ns*nt na] <- trick! projection views
ctruf i_is_ns_nt, // 1 if [ns nt] or 0 if [nt ns]
cchar *systype, // nn1 pd1 sf1 ...
cfloat *beta, // [na] source angles [radians]
cfloat *source_z, // [na] source positions, possibly null
cfloat *offset_s, // [na] possibly null
cfloat *offset_t, // [na] possibly null
cint nthread, // # of threads
cint iblock, // for OS, which block
cint nblock, // for OS, # blocks
ctruf compact, // 0 (for now)
cbct_back_init,
cfloat scale,
cint iz_start, // usually 0
cint iz_end, // do [iz_start, iz_end)
cint chat);


// cbct,*,proj.c

typedef sof cbct_any_proj1_type(
cfloat *image, // [nz nx ny] <- trick!
cint nx,
cint ny,
cint nz,
crscalar dx,
crscalar dy, // can be negative to cause flip
crscalar dz,
crscalar offset_x, // center offset in pixels (usually 0)
crscalar offset_y,
crscalar offset_z_shift, // offset_z + source_z / dz
cbyte *mask2, // [nx ny] 2D support mask: 0, 1, ..., nthread
crscalar dso, // distance from source to isocenter
crscalar dsd, // distance from source to detector
crscalar dfs, // distance from focal point to source (0 or inf)
cint ns,
cint nt,
crscalar ds, // horizontal ray spacing
crscalar dt, // vertical ray spacing
crscalar beta, // source angle [radians]
crscalar offset_s, // channel offset [pixels]
crscalar offset_t, // vertical offset on detector [pixels]
float *proj_out, // [ns*nt] <- trick! projection view at angle beta
ctruf o_is_ns_nt, // 1 if output should be [ns nt], 0 if [nt ns]
cbct_work *cw, // work space
cfloat scale); // scale output projections by this factor

extern cbct_any_proj1_type cbct_nn1_proj1; // cbct,nn1,proj.c
extern cbct_any_proj1_type cbct_pd1_proj1; // cbct,pd1,proj.c
extern cbct_any_proj1_type cbct_sf1_proj1; // cbct,sf1,proj.c
extern cbct_any_proj1_type cbct_sf2_proj1; // cbct,sf2,proj.c
extern cbct_any_proj1_type cbct_sf3_proj1; // cbct,sf3,proj.c
extern cbct_any_proj1_type cbct_sf4_proj1; // cbct,sf4,proj.c
extern cbct_any_proj1_type cbct_sf5_proj1; // cbct,sf5,proj.c
// cbct,any,proj,t.c

extern sof cbct_any_proj_t(
cfloat *image, // [nz nx ny] <- trick!
const cbct_ig *ig,
const cbct_cg *cg,
cbct_work *cw,
cint na, // # of views
float *proj, // [ns*nt na] <- trick! projection views
ctruf o_is_ns_nt, // 1 for [ns nt] or 0 for [nt ns]
cchar *systype, // nn1 pd1 sf1 ...
cfloat *beta, // [na] source angles [radians]
cfloat *source_z, // [na] source positions, possibly null
cfloat *offset_s, // [na] possibly null
cfloat *offset_t, // [na] possibly null
cint nthread, // # of threads
cint iblock, // for OS, which block
cint nblock, // for OS, # blocks
ctruf compact, // 0 (for now)
cfloat scale,
cint chat);


// cbct,work.c

extern truf isnan_vect_f(cfloat *v, cint nn);

extern float cbct_inprodf(
register cfloat *v1, // [nn]
register cfloat *v2, // [nn]
cint nn);

extern cbct_geom cbct_geom_parse(cfloat dso, cfloat dsd, cfloat dfs);

extern cbct_work *cbct_work_alloc(
const cbct_ig *ig,
const cbct_cg *cg,
cint nthread,
cfloat dl);

extern sof cbct_work_free(cbct_work *cw, cint nthread);

extern void cbct_view_transpose(
float *po, // [n2 n1]
cfloat *pi, // [n1 n2]
cint n1, cint n2);

extern cfloat *cbct_back_view_prep(
cfloat *proj, // [ns nt] or [nt ns]
ctruf i_ns_nt, // [ns nt] input?
ctruf o_ns_nt, // [ns nt] output?
float *work, // [ns*nt] work space
cint ns,
cint nt,
cfloat scale);

extern float cbct_pd1_scale(
const cbct_cg *cg,
const cbct_ig *ig,
cint chat);

extern sof cbct_footprint_size(
float *p_foot_s,
const cbct_cg *cg,
const cbct_ig *ig,
cint chat);


// cbct,incs.c

extern sof cbct_sf_incs(
cfloat *inprod_t, // [nt]
register cfloat *weight_s, // [nbins]
cint nbins,
cint ismin,
cint it_start, // typically 0
cint it_stop, // typically nt
float *proj, // [ns*nt] projection data to be incremented
ctruf is_ns_nt, // 1 for [ns nt] or 0 for [nt ns]
cint ns,
cint nt);


// cbct,sf1,misc.c

extern void cbct_sf1_proj_val(
register cfloat inprod_t,
register cfloat *pws, // [nbins] weight_s
register float *pps, // [>= nbins] proj
cint nbins);

extern void cbct_sf_inner_prods(
register float *inprod, // [nt] resulting inner products
register cfloat *weight_s, // [nbins]
cint nbins,
cint it_start, // typically 0
cint it_stop, // typically nt
register cfloat *pp, // [ns nt] input data
cint ns);

extern void cbct_sf1_sort_taus(float *taus);
extern float cbct_sf1_weight_value_s(cfloat *taus, cfloat bc);

#endif // jf_cbct_def_h
