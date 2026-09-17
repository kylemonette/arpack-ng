c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dhaup2
c
c\Description:
c  Intermediate level interface called by dhaupd.
c
c  This is a local fork of dqaup2 that replaces the exact-Ritz
c  arrowhead restart with a HARMONIC thick restart. Each cycle, after
c  the Arnoldi expansion, this routine forms the harmonic matrix
c
c      G = H + RNORM**2 * (H'\E_{kplusp}) * E_{kplusp}^T
c
c  (RNORM = the current residual norm, E_{kplusp} the last unit
c  vector), computes its real Schur decomposition via LAPACK's DGEES
c  (G is not Hessenberg, so DLAHQR does not apply), and uses the
c  resulting harmonic Ritz values -- instead of dqaup2's exact Ritz
c  values of H itself -- both to decide how many Ritz pairs are
c  converged/kept each cycle and to select which are passed into
c  dhapps (the harmonic-restart fork of dqapps) for the thick restart.
c
c  Deciding how many to keep follows the two-stage scheme of the
c  reference MATLAB implementation (see \Remarks): a cheap pre-reorder
c  residual estimate (via eigenvectors of the un-reordered harmonic
c  Schur form) sizes the working count K from the base desired count
c  K0, and a separate, trustworthy post-reorder residual check (via
c  Schur vectors of a K0-only reordering) certifies convergence.
c  Because LAPACK's DTRSEN only supports a boolean (two-way) partition,
c  unlike the reference implementation's fine-grained ORDSCHUR
c  clustering, this routine calls DTRSEN twice per cycle: once on a
c  disposable copy of the harmonic Schur form to certify convergence
c  against the K0 best harmonic Ritz pairs, and once on the working
c  copy to select the K best for the actual restart.
c
c\Usage:
c  call dhaup2
c     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
c       ISHIFT, MXITER, V, LDV, H, LDH, RITZR, RITZI, BOUNDS,
c       Q, LDQ, WORKL, IPNTR, WORKD, INFO )
c
c\Arguments
c  Identical to dqaup2/dnaup2 -- see dnaup2.f for the full description
c  of every argument. There is no HOUSE argument: the harmonic thick
c  restart only ever reduces via Householder reflectors.
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
c     dgetv0   ARPACK initial vector generation routine.
c     dnaitr   ARPACK Arnoldi factorization routine.
c     dhapps   Local harmonic-restart fork of dqapps.
c     dnconv   ARPACK convergence of Ritz values routine.
c     dgetrf   LAPACK LU factorization.
c     dgetrs   LAPACK routine that solves a system using an LU
c              factorization from DGETRF (used here, with TRANS='T',
c              to solve H^T * x = e_{kplusp} without ever forming H^T
c              explicitly).
c     dgees    LAPACK routine that computes the real Schur form (and
c              Schur vectors) of a general (non-Hessenberg) matrix.
c              Used for the harmonic matrix G, which DLAHQR cannot
c              handle directly.
c     dtrevc   LAPACK routine that computes eigenvectors of a quasi-
c              triangular (real Schur) matrix.
c     dtrsen   LAPACK routine that reorders a real Schur factorization
c              so selected eigenvalues appear in the leading block.
c     dsortc   ARPACK sorting routine.
c     ivout    ARPACK utility routine that prints integers.
c     arscnd   ARPACK utility routine for timing.
c     dlacpy   LAPACK matrix copy routine.
c     dlaset   LAPACK matrix initialization routine.
c     dcopy    Level 1 BLAS that copies one vector to another.
c     ddot     Level 1 BLAS that computes the scalar product of two
c              vectors.
c     dnrm2    Level 1 BLAS that computes the norm of a vector.
c
c\Remarks
c  1. The reference MATLAB implementation (harmonics/thick_harm_arnoldi.m
c     in the nonsymmetric-arrowhead repository) computes a cheap
c     residual estimate from the eigenvectors of the un-reordered
c     harmonic Schur form to decide how many Ritz pairs to keep each
c     cycle, then a separate, trustworthy residual check (from Schur
c     vectors of the reordered form) to certify convergence -- noting
c     the former "is not trustworthy enough to certify convergence
c     because the last entry of the eigenvectors can get tiny". This
c     routine follows the same two-stage scheme (see \Description).
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dhaup2
     &   ( ido, bmat, n, which, nev, np, tol, resid, mode, iupd,
     &     ishift, mxiter, v, ldv, h, ldh, ritzr, ritzi, bounds,
     &     q, ldq, workl, ipntr, workd, info )
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
      integer    ido, info, ishift, iupd, mode, ldh, ldq, ldv, mxiter,
     &           n, nev, np
      Double precision
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(13)
      Double precision
     &           bounds(nev+np), h(ldh,nev+np), q(ldq,nev+np), resid(n),
     &           ritzi(nev+np), ritzr(nev+np), v(ldv,nev+np),
     &           workd(3*n), workl( (nev+np)*(nev+np+3) )
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
      logical    cnorm , getv0, initv, update, ushift
      integer    ierr  , iter , j    , kplusp, msglvl, nconv,
     &           nev0 , np0  , numcnv
      Double precision
     &           rnorm , eps23, harmanrm
      save       cnorm , getv0, initv, update, ushift,
     &           rnorm , iter , eps23, kplusp, msglvl, nconv ,
     &           nev0 , np0  , numcnv, harmanrm
c
c     %-----------------------%
c     | Local array arguments |
c     %-----------------------%
c
      integer    kp(4)
c
c     %--------------------------------------------------------------%
c     | Local workspace for the harmonic matrix G, its real Schur    |
c     | decomposition, and the block-selection bookkeeping used to   |
c     | decide K0/K each cycle (see \Description). Sized by LDH, the |
c     | same bound H itself uses (the actual size in use at any      |
c     | point is KPLUSP <= LDH). SAVEd since they must survive the   |
c     | possible reverse-communication exit for B*x (ushift is not   |
c     | used here -- ISHIFT=0/user shifts is not offered by the      |
c     | harmonic thick restart -- but the B*x exit still applies).   |
c     %--------------------------------------------------------------%
c
      integer    harmierr, harmnblk, harmm, harmm0, i, j2, iconj,
     &           lwork, numconv_est, k0, kwork, ipiv(ldh), stdierr,
     &           stdm
      integer    harmblockstart(ldh), harmblocklen(ldh),
     &           harmconsumed(ldh)
      logical    harmselect(ldh), harmbwork(ldh), stdselect(ldh)
      Double precision
     &           harmwork(7*ldh),
     &           harmsep1, harmsep2, hinvem(ldh), qwork(1),
     &           harmeig1r(ldh), harmeig1i(ldh), harmbidx(ldh),
     &           temp
c
c     %--------------------------------------------------------------%
c     | HARMG/HARMSCHUR/HARMSCHURVEC/HARMWR/HARMWI carry the         |
c     | harmonic matrix G and its Schur decomposition from their     |
c     | computation site (right after the Arnoldi expansion, where   |
c     | the residual estimates used to size K0/K are read off it) to |
c     | the restart-selection site below, so they are SAVEd          |
c     | allocatables (re-sized if LDH ever changes between problems).|
c     | HARMRESEST holds the cheap pre-reorder residual estimate;    |
c     | HARMCOPY0/HARMCOPYVEC0 are disposable copies used for the    |
c     | K0-only convergence check, which must not disturb the        |
c     | working decomposition used for the actual restart.           |
c     %--------------------------------------------------------------%
c
      Double precision, allocatable, save ::
     &           harmg(:,:), harmschur(:,:), harmschurvec(:,:),
     &           harmqorig(:,:), harmwr(:), harmwi(:), harmresest(:),
     &           harmcopy0(:,:), harmcopyvec0(:,:), harmres0(:),
     &           stdschur(:,:), stdschurvec(:,:), stdwr(:), stdwi(:),
     &           stdbnd(:)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy  , dgetv0 , dnaitr , dnconv ,
     &           dhapps , dvout  , ivout  , arscnd,
     &           dgees  , dgetrf , dgetrs , dtrevc , dtrsen ,
     &           dsortc , dlacpy , dlaset , dlahqr , dscal
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           ddot , dnrm2 , dlamch
      external   ddot , dnrm2 , dlamch, dlapy2
      logical    dhaupdselectdummy
      external   dhaupdselectdummy
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    min, max, abs, sqrt
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (ido .eq. 0) then
c
         call arscnd (t0)
c
         msglvl = mnaup2
c
         eps23 = dlamch ('Epsilon-Machine')
         eps23 = eps23**(2.0D+0  / 3.0D+0 )
c
         nev0   = nev
         np0    = np
c
         kplusp = nev + np
         nconv  = 0
         iter   = 0
         harmanrm = zero
c
         getv0    = .true.
         update   = .false.
         ushift   = .false.
         cnorm    = .false.
c
         if (info .ne. 0) then
            initv = .true.
            info  = 0
         else
            initv = .false.
         end if
      end if
c
   10 continue
c
      if (getv0) then
         call dgetv0  (ido, bmat, 1, initv, n, 1, v, ldv, resid, rnorm,
     &                ipntr, workd, info)
c
         if (ido .ne. 99) go to 9000
c
         if (rnorm .eq. zero) then
            info = -9
            go to 1100
         end if
         getv0 = .false.
         ido  = 0
      end if
c
      if (update) go to 20
c
      if (cnorm)  go to 100
c
      call dnaitr  (ido, bmat, n, 0, nev, mode, resid, rnorm, v, ldv,
     &             h, ldh, ipntr, workd, info)
c
      if (ido .ne. 99) go to 9000
c
      if (info .gt. 0) then
         np   = info
         mxiter = iter
         info = -9999
         go to 1200
      end if
c
c     %--------------------------------------------------------------%
c     |                                                              |
c     |           M A I N  ARNOLDI  I T E R A T I O N  L O O P       |
c     |           Each iteration implicitly restarts the Arnoldi     |
c     |           factorization in place.                            |
c     |                                                              |
c     %--------------------------------------------------------------%
c
 1000 continue
c
         iter = iter + 1
c
         if (msglvl .gt. 0) then
            call ivout (logfil, 1, [iter], ndigit,
     &           '_haup2: **** Start of major iteration number ****')
         end if
c
         np  = kplusp - nev
c
         ido = 0
   20    continue
         update = .true.
c
         call dnaitr  (ido  , bmat, n  , nev, np , mode , resid,
     &                rnorm, v   , ldv, h  , ldh, ipntr, workd,
     &                info)
c
         if (ido .ne. 99) go to 9000
c
         if (info .gt. 0) then
            np = info
            mxiter = iter
            info = -9999
            go to 1200
         end if
         update = .false.
c
c        %--------------------------------------------------------%
c        | Standard real Schur decomposition of H (via DLAHQR),   |
c        | computed every cycle purely so genuine (non-harmonic)  |
c        | eigenvalues/bounds of the CURRENT, unrestarted H are   |
c        | available to hand dneupd at the exit path below --     |
c        | dneupd re-runs DLAHQR on this same H internally at     |
c        | extraction time and matches it positionally against    |
c        | whatever RITZR/RITZI/BOUNDS it is given, so those must  |
c        | be genuine eigenvalues of H, not harmonic ones (see     |
c        | dqaup2's own exit-path comment for the same reason).    |
c        %--------------------------------------------------------%
c
         if (allocated(stdschur)) then
            if (size(stdschur,1) .ne. ldh) then
               deallocate (stdschur, stdschurvec, stdwr, stdwi, stdbnd)
            end if
         end if
         if (.not. allocated(stdschur)) then
            allocate (stdschur(ldh,ldh), stdschurvec(ldh,ldh),
     &                stdwr(ldh), stdwi(ldh), stdbnd(ldh))
         end if
c
         call dlacpy ('All', kplusp, kplusp, h, ldh, stdschur, ldh)
         call dlaset ('All', kplusp, kplusp, zero, one, stdschurvec,
     &                ldh)
         call dlahqr (.true., .true., kplusp, 1, kplusp, stdschur, ldh,
     &                stdwr, stdwi, 1, kplusp, stdschurvec, ldh,
     &                stdierr)
         if (stdierr .ne. 0) then
            info = -8
            go to 1200
         end if
c
         call dtrevc ('R', 'A', stdselect, kplusp, stdschur, ldh,
     &                harmwork, 1, q, ldq, kplusp, stdm,
     &                workl(kplusp+1), stdierr)
         if (stdierr .ne. 0) then
            info = -8
            go to 1200
         end if
c
         iconj = 0
         do 26 j2 = 1, kplusp
            if ( abs( stdwi(j2) ) .le. zero ) then
               temp = dnrm2( kplusp, q(1,j2), 1 )
               call dscal ( kplusp, one / temp, q(1,j2), 1 )
            else
               if (iconj .eq. 0) then
                  temp = dlapy2( dnrm2( kplusp, q(1,j2), 1 ),
     &                           dnrm2( kplusp, q(1,j2+1), 1 ) )
                  call dscal ( kplusp, one / temp, q(1,j2), 1 )
                  call dscal ( kplusp, one / temp, q(1,j2+1), 1 )
                  iconj = 1
               else
                  iconj = 0
               end if
            end if
   26    continue
c
         call dcopy (kplusp, stdschurvec(kplusp,1), ldh, stdbnd, 1)
         call dgemv ('T', kplusp, kplusp, one, q, ldq, stdbnd, 1,
     &               zero, workl, 1)
c
         iconj = 0
         do 33 j2 = 1, kplusp
            if ( abs( stdwi(j2) ) .le. zero ) then
               stdbnd(j2) = rnorm * abs( workl(j2) )
            else
               if (iconj .eq. 0) then
                  stdbnd(j2) = rnorm * dlapy2( workl(j2), workl(j2+1) )
                  stdbnd(j2+1) = stdbnd(j2)
                  iconj = 1
               else
                  iconj = 0
               end if
            end if
   33    continue
c
c        %--------------------------------------------------------%
c        | Allocate/resize the harmonic-decomposition workspace.  |
c        %--------------------------------------------------------%
c
         if (allocated(harmg)) then
            if (size(harmg,1) .ne. ldh) then
               deallocate (harmg, harmschur, harmschurvec, harmqorig,
     &                     harmwr, harmwi, harmresest, harmcopy0,
     &                     harmcopyvec0, harmres0)
            end if
         end if
         if (.not. allocated(harmg)) then
            allocate (harmg(ldh,ldh), harmschur(ldh,ldh),
     &                harmschurvec(ldh,ldh), harmqorig(ldh,ldh),
     &                harmwr(ldh), harmwi(ldh),
     &                harmresest(ldh), harmcopy0(ldh,ldh),
     &                harmcopyvec0(ldh,ldh), harmres0(ldh))
         end if
c
c        %----------------------------------------------------------%
c        | HINVEM = H(1:kplusp,1:kplusp)' \ E_{kplusp}, solved via  |
c        | an LU factorization of H itself and DGETRS with          |
c        | TRANS='T' (solves the transposed system without ever     |
c        | forming H^T explicitly). H is destroyed by DGETRF, so a  |
c        | copy (HARMSCHUR, reused purely as scratch here) is       |
c        | factored instead.                                        |
c        %----------------------------------------------------------%
c
         call dlacpy ('All', kplusp, kplusp, h, ldh, harmschur, ldh)
         call dlaset ('All', kplusp, 1, zero, zero, hinvem, ldh)
         hinvem(kplusp) = one
         call dgetrf (kplusp, kplusp, harmschur, ldh, ipiv, harmierr)
         if (harmierr .ne. 0) then
            info = -8
            go to 1200
         end if
         call dgetrs ('T', kplusp, 1, harmschur, ldh, ipiv, hinvem,
     &                ldh, harmierr)
         if (harmierr .ne. 0) then
            info = -8
            go to 1200
         end if
c
c        %----------------------------------------------------------%
c        | G = H + RNORM**2 * HINVEM * E_{kplusp}^T -- only the     |
c        | last column of G differs from H.                         |
c        %----------------------------------------------------------%
c
         call dlacpy ('All', kplusp, kplusp, h, ldh, harmg, ldh)
         do 25 i = 1, kplusp
            harmg(i,kplusp) = harmg(i,kplusp) + rnorm**2 * hinvem(i)
   25    continue
c
c        %----------------------------------------------------------%
c        | Real Schur decomposition of G (not Hessenberg, so DGEES  |
c        | is used instead of DLAHQR). Workspace query first.       |
c        %----------------------------------------------------------%
c
         call dgees ('V', 'N', dhaupdselectdummy, kplusp, harmg, ldh,
     &               harmm0, harmwr, harmwi, harmschurvec, ldh,
     &               qwork, -1, harmbwork, harmierr)
         lwork = max(int(qwork(1)), 8*kplusp)
         call dgees ('V', 'N', dhaupdselectdummy, kplusp, harmg, ldh,
     &               harmm0, harmwr, harmwi, harmschurvec, ldh,
     &               harmwork, lwork, harmbwork, harmierr)
         if (harmierr .ne. 0) then
            info = -8
            go to 1200
         end if
c
c        %----------------------------------------------------------%
c        | HARMG now holds T (the Schur form of G); HARMSCHURVEC    |
c        | holds Q (its orthogonal Schur vectors). Save both into   |
c        | persistent copies (HARMSCHUR/HARMQORIG) before DTREVC    |
c        | overwrites HARMSCHURVEC with eigenvectors of T below.    |
c        %----------------------------------------------------------%
c
         call dlacpy ('All', kplusp, kplusp, harmg, ldh, harmschur,
     &                ldh)
         call dlacpy ('All', kplusp, kplusp, harmschurvec, ldh,
     &                harmqorig, ldh)
c
c        %----------------------------------------------------------%
c        | Cheap pre-reorder residual estimate, from the            |
c        | eigenvectors of the (un-reordered) harmonic Schur form   |
c        | -- used only to size K0/K below, per \Remarks 1.          |
c        %----------------------------------------------------------%
c
         call dtrevc ('R', 'A', harmselect, kplusp, harmg, ldh,
     &                harmwork, 1, harmschurvec, ldh, kplusp, harmm,
     &                workl, harmierr)
         if (harmierr .ne. 0) then
            info = -8
            go to 1200
         end if
c
c        %----------------------------------------------------------%
c        | DTREVC overwrote HARMSCHURVEC with eigenvectors of T     |
c        | itself (in the Schur VECTOR basis) -- multiply by the    |
c        | ORIGINAL Schur vectors HARMQORIG to get eigenvectors of  |
c        | G, matching the reference implementation's PFULL = P*PE.|
c        %----------------------------------------------------------%
c
         call dgemmwrap (kplusp, harmqorig, ldh, harmschurvec, ldh,
     &                   workl, kplusp)
         do 27 j = 1, kplusp
            harmresest(j) = rnorm * workl(kplusp*(j-1)+kplusp)
   27    continue
c
         do 28 j = 1, kplusp
            harmanrm = max(harmanrm, abs(harmwr(j)), abs(harmwi(j)))
   28    continue
c
c        %----------------------------------------------------------%
c        | Restore HARMSCHURVEC to the genuine (orthogonal) Schur   |
c        | vectors of G -- DTREVC overwrote it above.               |
c        %----------------------------------------------------------%
c
         call dlacpy ('All', kplusp, kplusp, harmqorig, ldh,
     &                harmschurvec, ldh)
c
c        %----------------------------------------------------------%
c        | K0 = NEV0, adjusted (via the block-selection routine     |
c        | below) so a complex-conjugate pair is never split.       |
c        %----------------------------------------------------------%
c
         k0 = nev0
         call dhaupdselectblocks (kplusp, harmwr, harmwi, harmschur,
     &        ldh, which, k0, harmblockstart, harmblocklen,
     &        harmconsumed, harmselect, harmeig1r, harmeig1i,
     &        harmbidx, harmnblk)
c
         numconv_est = 0
         do 29 j = 1, kplusp
            if (harmselect(j) .and.
     &          abs(harmresest(j)) .lt. tol*max(harmanrm,eps23))
     &          numconv_est = numconv_est + 1
   29    continue
c
         kwork = k0 + min(numconv_est, (kplusp-k0)/2)
         kwork = min(kwork, kplusp-1)
         if (kwork .lt. k0) kwork = k0
         call dhaupdselectblocks (kplusp, harmwr, harmwi, harmschur,
     &        ldh, which, kwork, harmblockstart, harmblocklen,
     &        harmconsumed, harmselect, harmeig1r, harmeig1i,
     &        harmbidx, harmnblk)
c
c        %----------------------------------------------------------%
c        | Trustworthy post-reorder residual check, certifying      |
c        | convergence against the K0 best harmonic Ritz pairs.     |
c        | Uses a DISPOSABLE copy so the working decomposition      |
c        | (used for the actual restart below) is undisturbed.      |
c        %----------------------------------------------------------%
c
         call dlacpy ('All', kplusp, kplusp, harmschur, ldh,
     &                harmcopy0, ldh)
         call dlacpy ('All', kplusp, kplusp, harmschurvec, ldh,
     &                harmcopyvec0, ldh)
         call dcopy (kplusp, harmwr, 1, workl, 1)
         call dcopy (kplusp, harmwi, 1, workl(kplusp+1), 1)
         call dhaupdselectblocks (kplusp, workl, workl(kplusp+1),
     &        harmcopy0, ldh, which, k0, harmblockstart,
     &        harmblocklen, harmconsumed, harmselect, harmeig1r,
     &        harmeig1i, harmbidx, harmnblk)
         call dtrsen ('N', 'V', harmselect, kplusp, harmcopy0, ldh,
     &                harmcopyvec0, ldh, workl, workl(kplusp+1),
     &                harmm, harmsep1, harmsep2, harmwork, 7*ldh,
     &                ipiv, 1, harmierr)
         if (harmierr .ne. 0) then
            info = -8
            go to 1200
         end if
         do 31 j = 1, k0
            harmres0(j) = rnorm * harmcopyvec0(kplusp,j)
   31    continue
c
         nconv = 0
         do 32 j = 1, k0
            if (abs(harmres0(j)) .lt. tol*max(harmanrm,eps23))
     &          nconv = nconv + 1
   32    continue
         numcnv = k0
c
         if (msglvl .gt. 2) then
            kp(1) = kwork
            kp(2) = kplusp - kwork
            kp(3) = numcnv
            kp(4) = nconv
            call ivout (logfil, 4, kp, ndigit,
     &                  '_haup2: K, NP, NUMCNV, NCONV are')
         end if
c
         if ( (nconv .ge. numcnv) .or. (iter .gt. mxiter) ) then
c
            h(3,1) = rnorm
            nev = numcnv
            np = nconv
            if (iter .gt. mxiter .and. nconv .lt. numcnv) info = 1
c
c           %------------------------------------------------------%
c           | Hand dneupd genuine (non-harmonic) eigenvalues/bounds|
c           | of the current, unrestarted H -- see the standard    |
c           | DLAHQR computation above for why.                    |
c           %------------------------------------------------------%
c
            call dcopy (kplusp, stdwr, 1, ritzr, 1)
            call dcopy (kplusp, stdwi, 1, ritzi, 1)
            call dcopy (kplusp, stdbnd, 1, bounds, 1)
            go to 1100
c
         end if
c
         nev = kwork
         np = kplusp - kwork
c
c        %----------------------------------------------------------%
c        | Reorder the working Schur form so the KWORK selected     |
c        | ("desired") harmonic Ritz values occupy the leading      |
c        | block, updating the Schur vectors to match.               |
c        %----------------------------------------------------------%
c
         call dtrsen ('N', 'V', harmselect, kplusp, harmschur, ldh,
     &                harmschurvec, ldh, harmwr, harmwi, harmm,
     &                harmsep1, harmsep2, harmwork, 7*ldh, ipiv, 1,
     &                harmierr)
         if (harmierr .ne. 0 .or. harmm .ne. nev) then
            info = -8
            go to 1200
         end if
c
c        %----------------------------------------------------------%
c        | Apply the harmonic thick restart: build the extended     |
c        | basis from the NEV desired harmonic Schur pairs, reduce  |
c        | to upper Hessenberg form, and read the updated H, V, and |
c        | RESID off that result.                                   |
c        %----------------------------------------------------------%
c
         call dhapps (n, nev, np, v, ldv, h, ldh, resid,
     &        harmschur, ldh, harmschurvec, ldh, hinvem, workd)
c
         cnorm = .true.
         call arscnd (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy  (n, resid, 1, workd(n+1), 1)
            ipntr(1) = n + 1
            ipntr(2) = 1
            ido = 2
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy  (n, resid, 1, workd, 1)
         end if
c
  100    continue
c
         if (bmat .eq. 'G') then
            call arscnd (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
c
         if (bmat .eq. 'G') then
            rnorm = ddot  (n, resid, 1, workd, 1)
            rnorm = sqrt(abs(rnorm))
         else if (bmat .eq. 'I') then
            rnorm = dnrm2 (n, resid, 1)
         end if
         cnorm = .false.
c
      go to 1000
c
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
 1100 continue
c
      mxiter = iter
c
 1200 continue
      ido = 99
c
      call arscnd (t1)
      tnaup2 = t1 - t0
c
 9000 continue
c
      return
c
c     %-----------------%
c     | End of dhaup2   |
c     %-----------------%
c
      end
c
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dhaupdselectdummy
c
c\Description:
c  Dummy SELECT function argument for DGEES, required by its interface
c  but never invoked since SORT='N' is always passed.
c
c\EndDoc
c-----------------------------------------------------------------------
c
      logical function dhaupdselectdummy (wr, wi)
      Double precision wr, wi
      dhaupdselectdummy = .false.
      return
      end
c
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dhaupdselectblocks
c
c\Description:
c  Classifies the diagonal blocks (1x1 or 2x2) of a real Schur form,
c  sorts the block-representative eigenvalues by the WHICH preference,
c  and greedily marks whole blocks (from the most-wanted end) as
c  selected until their lengths sum to at least NEV -- growing NEV by
c  one first if the last block taken would otherwise split a complex-
c  conjugate pair. This is the same block-consumption scheme dqaup2
c  uses for the exact-Ritz arrowhead restart, factored out here since
c  dhaup2 needs it twice per cycle (once to size K0/K, once to select
c  the actual restart block).
c
c\Usage:
c  call dhaupdselectblocks
c     ( KPLUSP, WR, WI, SCHUR, LDSCHUR, WHICH, NEV, BLOCKSTART,
c       BLOCKLEN, CONSUMED, SELECT, EIG1R, EIG1I, BIDX, NBLK )
c
c\Arguments
c  KPLUSP  Integer.  (INPUT) Size of the Schur form.
c  WR, WI  Double precision arrays of length KPLUSP.  (INPUT) Real and
c          imaginary parts of the Schur form's eigenvalues.
c  SCHUR   Double precision KPLUSP by KPLUSP array.  (INPUT) The real
c          Schur form (only its subdiagonal is examined, to classify
c          1x1 vs 2x2 blocks).
c  LDSCHUR Integer.  (INPUT) Leading dimension of SCHUR.
c  WHICH   Character*2.  (INPUT) Shift/selection criteria, as in
c          dnaupd.
c  NEV     Integer.  (INPUT/OUTPUT) On input, the desired number of
c          selected eigenvalues; on output, possibly increased by one
c          to keep a complex-conjugate pair together.
c  BLOCKSTART, BLOCKLEN, CONSUMED  Integer arrays of length KPLUSP.
c          (WORKSPACE)
c  SELECT  Logical array of length KPLUSP.  (OUTPUT) .TRUE. at
c          positions belonging to a selected block.
c  EIG1R, EIG1I, BIDX  Double precision arrays of length KPLUSP.
c          (WORKSPACE)
c  NBLK    Integer.  (OUTPUT) Number of diagonal blocks found.
c
c\EndDoc
c-----------------------------------------------------------------------
c
      subroutine dhaupdselectblocks
     &   ( kplusp, wr, wi, schur, ldschur, which, nev, blockstart,
     &     blocklen, consumed, select, eig1r, eig1i, bidx, nblk )
c
      character  which*2
      integer    kplusp, ldschur, nev, nblk
      integer    blockstart(kplusp), blocklen(kplusp),
     &           consumed(kplusp)
      logical    select(kplusp)
      Double precision
     &           wr(kplusp), wi(kplusp), schur(ldschur,kplusp),
     &           eig1r(kplusp), eig1i(kplusp), bidx(kplusp)
c
      Double precision zero
      parameter (zero = 0.0D+0)
      integer    i, j, pick, m
      external   dsortc
c
      nblk = 0
      i = 1
   10 continue
      if (i .le. kplusp) then
         nblk = nblk + 1
         blockstart(nblk) = i
         if (i .lt. kplusp .and. abs(schur(i+1,i)) .gt. zero) then
            blocklen(nblk) = 2
            i = i + 2
         else
            blocklen(nblk) = 1
            i = i + 1
         end if
         eig1r(nblk) = wr(blockstart(nblk))
         eig1i(nblk) = abs(wi(blockstart(nblk)))
         bidx(nblk) = dble(nblk)
         go to 10
      end if
c
      call dsortc (which, .true., nblk, eig1r, eig1i, bidx)
c
      do 20 j = 1, nblk
         consumed(j) = 0
   20 continue
c
      m = 0
      j = nblk
   30 continue
      if (j .ge. 1 .and. m .lt. nev) then
         pick = nint(bidx(j))
         consumed(pick) = 1
         m = m + blocklen(pick)
         j = j - 1
         go to 30
      end if
c
      if (m .gt. nev) nev = m
c
      do 50 j = 1, nblk
         do 40 i = 0, blocklen(j)-1
            select(blockstart(j)+i) = (consumed(j) .eq. 1)
   40    continue
   50 continue
c
      return
      end
c
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dgemmwrap
c
c\Description:
c  C = A * B for square N by N A, B, C -- a thin DGEMM wrapper used
c  where an explicit temporary result buffer is clearer than an
c  in-place LAPACK call.
c
c\EndDoc
c-----------------------------------------------------------------------
c
      subroutine dgemmwrap (n, a, lda, b, ldb, c, ldc)
      integer n, lda, ldb, ldc
      Double precision a(lda,n), b(ldb,n), c(ldc,n)
      Double precision one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
      external dgemm
      call dgemm ('N', 'N', n, n, n, one, a, lda, b, ldb, zero, c, ldc)
      return
      end
