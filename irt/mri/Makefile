# mri/Makefile

# compile MRI-related .mex* files

include $(HOME)/l/src/make/Makedef.$(Arch)

#suf_i686 = mexglx
#suf_ppc = mexmac
#suf_ = suf_i686; # default is linux
#suf = $(suf_$(Arch))
suf = mex$(MatlabSuf)
# from Makedef.*

dir = ../mex/v7
mex = $(MatlabCom)
#mex = mex

goal1 = $(dir)/mri_exp_mult_mex.$(suf)
goal2 = $(dir)/exp_xform_mex.$(suf)

all:
	$(MAKE) $(goal1) # date=-DDate='\"date='`date +%Y-%m-%d,%T`'\"'
	$(MAKE) $(goal2) # date=-DDate='\"date='`date +%Y-%m-%d,%T`'\"'

aflag_x86_64 = -fPIC
aflag_mi64 = -arch x86_64
aflag = $(aflag_$(Arch))

CFLAGS = -D_GNU_SOURCE $(date) \
	-std=c99 -Wall -Wshadow -Wpointer-arith -W -Werror $(aflag)

com = $(mex) CFLAGS='$(CFLAGS) -fexceptions'

$(goal1):	mri_exp_mult_mex.c mri_exp_mult_template.c
	echo $(date)
	$(com) mri_exp_mult_mex.c
	mv mri_exp_mult_mex.$(suf) $@
	chmod 755 $@

$(goal2):	exp_xform_mex.c exp_xform_template.c
	$(com) exp_xform_mex.c
	mv exp_xform_mex.$(suf) $@
	chmod 755 $@

reset:
	@echo suf = $(suf)
	rm -i $(goal1) $(goal2)
