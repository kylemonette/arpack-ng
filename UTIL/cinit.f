c-----------------------------------------------------------------------
c\Name: cinit
c
c\Description:
c  Initializes the ARPACK debug.h / stat.h COMMON block variables that
c  mydsapps reads (logfil, ndigit, msapps, tsapps). Normally an ARPACK
c  driver program (e.g. dsaupd's caller) sets these once at start-up;
c  since mydsapps is being called directly from a MEX gateway with no
c  such driver, this routine stands in for that step. Call it once
c  before the first call to mydsapps.
c
c\Usage:
c  call cinit()
c
c-----------------------------------------------------------------------
c
	subroutine cinit()
c
	include   'debug.h'
	include   'stat.h'
c
c     %--------------------------------------------------------%
c     | logfil: Fortran unit number for ARPACK debug printing. |
c     | ndigit: negative value = number of digits per line for |
c     |         the ivout/dvout debug printers. Non-positive   |
c     |         msapps below means those printers are never    |
c     |         actually invoked, so this is mostly a formality|
c     %--------------------------------------------------------%
c
	logfil = 6
	ndigit = -3
c
c     %----------------------------------------------------%
c     | msapps: message/debug level for the *sapps routines |
c     | 0 = silent. Raise this only if you want mydsapps to |
c     | print diagnostic output via ivout/dvout.            |
c     %----------------------------------------------------%
c
	msapps = 0
c
c     %------------------------------------------%
c     | tsapps: cumulative timing accumulator for |
c     | mydsapps, updated in-place by the routine |
c     | itself via arscnd. Must start at zero.    |
c     %------------------------------------------%
c
	tsapps = 0.0D+0
c
	return
c
c     %---------------%
c     | End of cinit  |
c     %---------------%
c
	end
