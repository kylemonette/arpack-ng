c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: mydnaupd
c
c\Description:
c
c  Reverse communication interface for the Implicitly Restarted Arnoldi
c  Iteration, for real nonsymmetric problems. This is a local fork of
c  dnaupd: the public interface (arguments) is UNCHANGED, but internally
c  every restart is performed via mydnaup2 -- a fork of dnaup2 that
c  applies the NONSYMMETRIC ARROWHEAD restart (mydnapps, via
c  narrowgivens) instead of classical implicit-shift bulge chasing. See
c
c    "A Unified View of Arrowhead Matrix Transformations and Lanczos
c     Restarts" (and nonsymmetric companion), James Baglama,
c     Kyle Monette, Vasilije Perovic (2026)
c
c\Usage:
c  call mydnaupd
c     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
c       IPNTR, WORKD, WORKL, LWORKL, INFO )
c
c\Arguments
c  Identical to dnaupd -- see dnaupd.f for the full description of
c  every argument (IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV,
c  IPARAM, IPNTR, WORKD, WORKL, LWORKL, INFO all mean exactly what they
c  mean there). WHICH must be one of 'LM','SM','LR','SR','LI','SI'.
c  LWORKL must be at least 3*NCV**2 + 6*NCV.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. A Unified View of Arrowhead Matrix Transformations and Lanczos
c     Restarts (and nonsymmetric companion), James Baglama,
c     Kyle Monette, Vasilije Perovic, (2026).
c
c\Routines called:
c     mydnaup2 Local fork of dnaup2 that additionally computes the real
c              Schur decomposition of the current KEV+NP Hessenberg
c              matrix, reorders it, and passes it into mydnapps.
c     dstats   ARPACK routine that initializes timing and other
c              statistics variables.
c     ivout    ARPACK utility routine that prints integers.
c     arscnd   ARPACK utility routine for timing.
c     dvout    ARPACK utility routine that prints vectors.
c     dlamch   LAPACK routine that determines machine constants.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\SCCS Information: @(#)
c FILE: mynaupd.F   SID: 2.8   DATE OF SID: 04/10/01   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine mydnaupd
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam,
     &     ipntr, workd, workl, lworkl, info )
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Double precision
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(11), ipntr(14)
      Double precision
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0 , zero = 0.0D+0 )
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    bounds, ierr, ih, iq, ishift, iupd, iw,
     &           ldh, ldq, msglvl, mxiter, mode, nb,
     &           nev0, next, np, ritzr, ritzi, j
      save       bounds, ierr, ih, iq, ishift, iupd, iw,
     &           ldh, ldq, msglvl, mxiter, mode, nb,
     &           nev0, next, np, ritzr, ritzi
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   mydnaup2 ,  dvout , ivout, arscnd, dstats
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dlamch
      external   dlamch
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (ido .eq. 0) then
c
c        %-------------------------------%
c        | Initialize timing statistics  |
c        | & message level for debugging |
c        %-------------------------------%
c
         call dstats
         call arscnd (t0)
         msglvl = mnaupd
c
         ierr   = 0
         ishift = iparam(1)
         mxiter = iparam(3)
         nb     = 1
c
c        %--------------------------------------------%
c        | Revision 2 performs only implicit restart. |
c        %--------------------------------------------%
c
         iupd   = 1
         mode   = iparam(7)
c
c        %----------------%
c        | Error checking |
c        %----------------%
c
         if (n .le. 0) then
            ierr = -1
         else if (nev .le. 0) then
            ierr = -2
         else if (ncv .le. nev+1 .or. ncv .gt. n) then
            ierr = -3
         end if
c
c        %----------------------------------------------%
c        | NP is the number of additional steps to      |
c        | extend the length NEV Arnoldi factorization. |
c        %----------------------------------------------%
c
         np     = ncv - nev
c
         if (mxiter .le. 0)                     ierr = -4
         if (which .ne. 'LM' .and.
     &       which .ne. 'SM' .and.
     &       which .ne. 'LR' .and.
     &       which .ne. 'SR' .and.
     &       which .ne. 'LI' .and.
     &       which .ne. 'SI')                    ierr = -5
         if (bmat .ne. 'I' .and. bmat .ne. 'G') ierr = -6
c
         if (lworkl .lt. 3*ncv**2 + 6*ncv)       ierr = -7
         if (mode .lt. 1 .or. mode .gt. 4) then
                                                ierr = -10
         else if (mode .eq. 1 .and. bmat .eq. 'G') then
                                                ierr = -11
         else if (ishift .lt. 0 .or. ishift .gt. 1) then
                                                ierr = -12
         end if
c
c        %------------%
c        | Error Exit |
c        %------------%
c
         if (ierr .ne. 0) then
            info = ierr
            ido  = 99
            go to 9000
         end if
c
c        %------------------------%
c        | Set default parameters |
c        %------------------------%
c
         if (nb .le. 0)                         nb = 1
         if (tol .le. zero)                     tol = dlamch ('EpsMach')
c
         np     = ncv - nev
         nev0   = nev
c
c        %-----------------------------%
c        | Zero out internal workspace |
c        %-----------------------------%
c
         do 10 j = 1, 3*ncv**2 + 6*ncv
            workl(j) = zero
 10      continue
c
c        %-------------------------------------------------------%
c        | Pointer into WORKL for address of H, RITZR, RITZI,    |
c        | BOUNDS, Q, etc... and the remaining workspace.        |
c        | Also update pointer to be used on output.             |
c        | Memory is laid out as follows:                        |
c        | workl(1:ncv*ncv) := generated upper Hessenberg matrix |
c        | workl(ncv*ncv+1:ncv*ncv+ncv) := real part of Ritz     |
c        | workl(ncv*ncv+ncv+1:ncv*ncv+2*ncv) := imag part       |
c        | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds  |
c        | workl(ncv*ncv+3*ncv+1:2*ncv*ncv+3*ncv) := Q workspace |
c        | remaining ncv*ncv+3*ncv locations := mydnaup2's own   |
c        |    WORKL workspace (used by dneigh, dngets, dnconv)   |
c        %-------------------------------------------------------%
c
         ldh    = ncv
         ldq    = ncv
         ih     = 1
         ritzr  = ih     + ldh*ncv
         ritzi  = ritzr  + ncv
         bounds = ritzi  + ncv
         iq     = bounds + ncv
         iw     = iq     + ncv*ncv
         next   = iw     + ncv*(ncv+3)
c
         ipntr(4) = next
         ipntr(5) = ih
         ipntr(6) = ritzr
         ipntr(7) = ritzi
         ipntr(8) = bounds
         ipntr(14) = iw
      end if
c
c     %-------------------------------------------------------%
c     | Carry out the Implicitly restarted Arnoldi Iteration. |
c     %-------------------------------------------------------%
c
      call mydnaup2
     &   ( ido, bmat, n, which, nev0, np, tol, resid, mode, iupd,
     &     ishift, mxiter, v, ldv, workl(ih), ldh, workl(ritzr),
     &     workl(ritzi), workl(bounds), workl(iq), ldq, workl(iw),
     &     ipntr, workd, info )
c
c     %--------------------------------------------------%
c     | ido .ne. 99 implies use of reverse communication |
c     | to compute operations involving OP or shifts.    |
c     %--------------------------------------------------%
c
      if (ido .eq. 3) iparam(8) = np
      if (ido .ne. 99) go to 9000
c
      iparam(3) = mxiter
      iparam(5) = np
      iparam(9) = nopx
      iparam(10) = nbx
      iparam(11) = nrorth
c
c     %------------------------------------%
c     | Exit if there was an informational |
c     | error within mydnaup2.              |
c     %------------------------------------%
c
      if (info .lt. 0) go to 9000
      if (info .eq. 2) info = 3
c
      if (msglvl .gt. 0) then
         call ivout (logfil, 1, [mxiter], ndigit,
     &               '_naupd: number of update iterations taken')
         call ivout (logfil, 1, [np], ndigit,
     &               '_naupd: number of "converged" Ritz values')
         call dvout  (logfil, np, workl(ritzr), ndigit,
     &               '_naupd: real part of final Ritz values')
         call dvout  (logfil, np, workl(ritzi), ndigit,
     &               '_naupd: imaginary part of final Ritz values')
         call dvout  (logfil, np, workl(bounds), ndigit,
     &               '_naupd: corresponding error bounds')
      end if
c
      call arscnd (t1)
      tnaupd = t1 - t0
c
      if (msglvl .gt. 0) then
c
c        %--------------------------------------------------------%
c        | Version Number & Version Date are defined in version.h |
c        %--------------------------------------------------------%
c
         write (6,1000)
         write (6,1100) mxiter, nopx, nbx, nrorth, nitref, nrstrt,
     &                  tmvopx, tmvbx, tnaupd, tnaup2, tnaitr, titref,
     &                  tgetv0, tneigh, tngets, tnapps, tnconv
 1000    format (//,
     &      5x, '=============================================',/
     &      5x, '= Nonsymmetric implicit Arnoldi update code =',/
     &      5x, '= (nonsymmetric arrowhead restart, Givens)  =',/
     &      5x, '= Version Number:', ' 2.4' , 20x, ' =',/
     &      5x, '= Version Date:  ', ' 07/31/96' , 15x, ' =',/
     &      5x, '=============================================',/
     &      5x, '= Summary of timing statistics             =',/
     &      5x, '=============================================',//)
 1100    format (
     &      5x, 'Total number update iterations             = ', i5,/
     &      5x, 'Total number of OP*x operations            = ', i5,/
     &      5x, 'Total number of B*x operations             = ', i5,/
     &      5x, 'Total number of reorthogonalization steps  = ', i5,/
     &      5x, 'Total number of iterative refinement steps = ', i5,/
     &      5x, 'Total number of restart steps              = ', i5,/
     &      5x, 'Total time in user OP*x operation          = ', f12.6,/
     &      5x, 'Total time in user B*x operation           = ', f12.6,/
     &      5x, 'Total time in Arnoldi update routine       = ', f12.6,/
     &      5x, 'Total time in naup2 routine                = ', f12.6,/
     &      5x, 'Total time in basic Arnoldi iteration loop = ', f12.6,/
     &      5x, 'Total time in reorthogonalization phase    = ', f12.6,/
     &      5x, 'Total time in (re)start vector generation  = ', f12.6,/
     &      5x, 'Total time in Hessenberg eig. subproblem   = ', f12.6,/
     &      5x, 'Total time in getting the shifts           = ', f12.6,/
     &      5x, 'Total time in applying the shifts          = ', f12.6,/
     &      5x, 'Total time in convergence testing          = ', f12.6)
      end if
c
 9000 continue
c
      return
c
c     %-----------------%
c     | End of mydnaupd |
c     %-----------------%
c
      end
