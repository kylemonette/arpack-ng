c\BeginDoc
c
c\Name: dhaupd
c
c\Description:
c  Reverse communication interface for the Implicitly Restarted Arnoldi
c  iteration, using the nonsymmetric HARMONIC thick restart (dhaup2/
c  dhapps) in place of dnaupd's classical implicit-shift bulge-chasing
c  or dqaupd's exact-Ritz arrowhead restart. This is a local fork of
c  dqaupd with the HOUSE argument removed: the harmonic thick restart
c  only ever reduces via Householder reflectors.
c
c  Mode 1 (A*x = lambda*x, OP = A, B = I) is the only mode currently
c  exercised by this fork.
c
c\Usage:
c  call dhaupd
c     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
c       IPNTR, WORKD, WORKL, LWORKL, INFO )
c
c\Arguments
c  Identical to dqaupd -- see dqaupd.f/dnaupd.f for the full
c  description of every argument. There is no HOUSE argument.
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
c  2. A Unified View of Arrowhead Matrix Transformations and Lanczos
c     Restarts (and nonsymmetric companion), James Baglama,
c     Kyle Monette, Vasilije Perovic, (2026).
c
c\Routines called:
c     dhaup2  Local fork of dqaup2 that computes the harmonic thick
c             restart each cycle (see dhaup2.f for details) instead of
c             an exact-Ritz arrowhead restart.
c     ivout   ARPACK utility routine that prints integers.
c     arscnd  ARPACK utility routine for timing.
c     dvout   ARPACK utility routine that prints vectors.
c     dlamch  LAPACK routine that determines machine constants.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dhaupd
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
     &           ldh, ldq, levec, mode, msglvl, mxiter, nb,
     &           nev0, next, np, ritzi, ritzr, j
      save       bounds, ih, iq, ishift, iupd, iw, ldh, ldq,
     &           levec, mode, msglvl, mxiter, nb, nev0, next,
     &           np, ritzi, ritzr
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dhaup2 , dvout , ivout, arscnd, dstatn
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
         call dstatn
         call arscnd (t0)
         msglvl = mnaupd
c
         ierr   = 0
         ishift = iparam(1)
         mxiter = iparam(3)
         nb     = 1
c
         iupd   = 1
         mode   = iparam(7)
c
         if (n .le. 0) then
             ierr = -1
         else if (nev .le. 0) then
             ierr = -2
         else if (ncv .le. nev+1 .or.  ncv .gt. n) then
             ierr = -3
         else if (mxiter .le. 0) then
             ierr = -4
         else if (which .ne. 'LM' .and.
     &       which .ne. 'SM' .and.
     &       which .ne. 'LR' .and.
     &       which .ne. 'SR' .and.
     &       which .ne. 'LI' .and.
     &       which .ne. 'SI') then
            ierr = -5
         else if (bmat .ne. 'I' .and. bmat .ne. 'G') then
            ierr = -6
         else if (lworkl .lt. 3*ncv**2 + 6*ncv) then
            ierr = -7
         else if (mode .lt. 1 .or. mode .gt. 4) then
                                                ierr = -10
         else if (mode .eq. 1 .and. bmat .eq. 'G') then
                                                ierr = -11
         else if (ishift .lt. 0 .or. ishift .gt. 1) then
                                                ierr = -12
         end if
c
         if (ierr .ne. 0) then
            info = ierr
            ido  = 99
            go to 9000
         end if
c
         if (nb .le. 0)				nb = 1
         if (tol .le. zero)			tol = dlamch ('EpsMach')
c
         np     = ncv - nev
         nev0   = nev
c
         do 10 j = 1, 3*ncv**2 + 6*ncv
            workl(j) = zero
  10     continue
c
         ldh    = ncv
         ldq    = ncv
         ih     = 1
         ritzr  = ih     + ldh*ncv
         ritzi  = ritzr  + ncv
         bounds = ritzi  + ncv
         iq     = bounds + ncv
         iw     = iq     + ldq*ncv
         next   = iw     + ncv**2 + 3*ncv
c
         ipntr(4) = next
         ipntr(5) = ih
         ipntr(6) = ritzr
         ipntr(7) = ritzi
         ipntr(8) = bounds
         ipntr(14) = iw
c
      end if
c
      call dhaup2
     &   ( ido, bmat, n, which, nev0, np, tol, resid, mode, iupd,
     &     ishift, mxiter, v, ldv, workl(ih), ldh, workl(ritzr),
     &     workl(ritzi), workl(bounds), workl(iq), ldq, workl(iw),
     &     ipntr, workd, info )
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
      if (info .lt. 0) go to 9000
      if (info .eq. 2) info = 3
c
      if (msglvl .gt. 0) then
         call ivout (logfil, 1, [mxiter], ndigit,
     &               '_haupd: Number of update iterations taken')
         call ivout (logfil, 1, [np], ndigit,
     &               '_haupd: Number of wanted "converged" Ritz values')
         call dvout  (logfil, np, workl(ritzr), ndigit,
     &               '_haupd: Real part of the final Ritz values')
         call dvout  (logfil, np, workl(ritzi), ndigit,
     &               '_haupd: Imaginary part of the final Ritz values')
      end if
c
      call arscnd (t1)
      tnaupd = t1 - t0
c
 9000 continue
c
      return
c
c     %-----------------%
c     | End of dhaupd   |
c     %-----------------%
c
      end
