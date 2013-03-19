// cbct,mask2.c
// mask operations

#include "cbct,def.h"

//
// cbct_mask_init()
// initialize mask with values for each thread in it
// input and output masks can be the same pointer!
//
sof cbct_mask_init(
byte *mask_int,	// [nx ny] out: integers 0 or 1 ... nthread
cbyte *mask_bin, // [nx ny] in: binary, possibly null
cint nx, cint ny,
cint nthread,
ctruf chat)
{
	if (nthread > 250) Fail("too many threads")

	long sum = 0;
	if (!mask_bin)
		sum = nx * ny;
	else
		for (int jj=0; jj < nx * ny; ++jj)
			sum += mask_bin[jj] != 0;

	const int nj_per_thread = (int) Ceilf(sum / (float) nthread);
	if (chat > 99)
		Note3("%ld voxels / %d threads = %d each",
			sum, nthread, nj_per_thread)

	int jj=0;
	for (int it=1; it <= nthread; ++it) {
		int nj_here = Min(sum, it * nj_per_thread) - (it-1) * nj_per_thread;
		if (chat > 99) Note2("thread %d nj %d", it-1, nj_here)

		while (nj_here) {
			if (!mask_bin || mask_bin[jj]) {
				mask_int[jj] = it;
				--nj_here;
			}
			else
				mask_int[jj] = 0;
			++jj;
		}
	}
	// if (jj != nx * ny) Fail("bug")
	// Iwrite2byte("mask-int.fld", mask_int, nx, ny)
	Ok
}
