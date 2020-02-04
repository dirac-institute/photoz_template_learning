

DESCRIBING COSMOS DATA:

zCOSMOS 10k-Bright reference:
https://iopscience.iop.org/article/10.1088/0067-0049/184/2/218/meta
https://www.eso.org/sci/observing/phase3/data_releases/zcosmos_dr3_b2.pdf

photometry source:
https://iopscience.iop.org/article/10.1088/0004-637X/690/2/1236/pdf
https://arxiv.org/abs/1604.02350

COSMOS filters:
http://cosmos.astro.caltech.edu/page/filterset

obtaining data from SQL module at http://cesam.lam.fr/hstcosmos/search/sql-module


Only keeping galaxies with zflag between 3 and 4.5. This is to keep all the best redshifts, per the Lilly paper, and includes all of the various decimal flags described there as well.

Only keeping galaxies with flag_maskB, flag_maskV, flag_maski, flag_maskz, flag_maskD all equal 0 (per the photometry README).

Type = 0 (galaxy)

NOTE I NEED TO INCLUDE LINKS TO THE READMEs.

Saving the following bands: