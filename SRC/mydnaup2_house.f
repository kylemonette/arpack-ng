c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: mydnaup2_house
c
c\Description:
c  Intermediate level interface called by mydnaupd.
c
c  This is a local fork of dnaup2 (arguments unchanged). Internally, it
c  computes the full real Schur decomposition of the current KEV+NP
c  upper Hessenberg matrix H (via LAPACK's DHSEQR) ONCE per iteration,
c  right after the Arnoldi expansion: the Ritz values and error bounds
c  used for convergence testing are read directly off it (replacing
c  dnaup2's dneigh call, whose internal QR pass this subsumes). At
c  shift-application time the SAME decomposition is reordered (via
c  LAPACK's DTRSEN) so the NEV desired eigenvalues occupy the leading
c  block, and that reordered Schur form/vectors are passed into
c  mydnapps_house (the House/DGEHRD-DORGHR variant of the nonsymmetric
c  arrowhead-restart fork of dnapps) instead of calling dnapps
c  directly with a shift list.
c
c\Usage:
c  call mydnaup2_house
c     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
c       ISHIFT, MXITER, V, LDV, H, LDH, RITZR, RITZI, BOUNDS,
c       Q, LDQ, WORKL, IPNTR, WORKD, INFO )
c
c\Arguments
c  Identical to dnaup2 -- see dnaup2.f for the full description of every
c  argument.
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
c     dgetv0   ARPACK initial vector generation routine.
c     dnaitr   ARPACK Arnoldi factorization routine.
c     mydnapps_house Local fork of dnapps: applies the nonsymmetric
c              arrowhead restart in place of bulge-chasing.
c     dnconv   ARPACK convergence of Ritz values routine.
c     dtrevc   LAPACK routine that computes eigenvectors of a quasi-
c              triangular (real Schur) matrix; used with the last row
c              of the Schur vectors to form the Ritz error bounds,
c              exactly as dneigh (now subsumed) did internally.
c     dngets   ARPACK reorder Ritz values and error bounds routine.
c     dsortc   ARPACK sorting routine.
c     ivout    ARPACK utility routine that prints integers.
c     arscnd   ARPACK utility routine for timing.
c     dmout    ARPACK utility routine that prints matrices
c     dvout    ARPACK utility routine that prints vectors.
c     dlamch   LAPACK routine that determines machine constants.
c     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dlahqr   LAPACK routine that computes the real Schur form (and
c              Schur vectors) of an upper Hessenberg matrix. Used
c              instead of the faster DHSEQR so the eigenvalue ordering
c              matches the DLAHQR call dneupd performs internally at
c              extraction time (a hard requirement of dneupd's
c              positional SELECT mechanism).
c     dtrsen   LAPACK routine that reorders a real Schur factorization
c              so selected eigenvalues appear in the leading block.
c     dsortc   ARPACK sorting routine (also used here to rank the fresh
c              Schur blocks by the user's WHICH preference).
c     dlacpy   LAPACK matrix copy routine.
c     dlaset   LAPACK matrix initialization routine.
c     dcopy    Level 1 BLAS that copies one vector to another.
c     ddot     Level 1 BLAS that computes the scalar product of two
c              vectors.
c     dnrm2    Level 1 BLAS that computes the norm of a vector.
c     dswap    Level 1 BLAS that swaps two vectors.
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
c FILE: mynaup2.F   SID: 2.8   DATE OF SID: 10/17/00   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine mydnaup2_house
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
      character  wprime*2
      logical    cnorm , getv0, initv, update, ushift
      integer    ierr  , iter , j    , kplusp, msglvl, nconv,
     &           nevbef, nev0 , np0  , nptemp, numcnv
      Double precision
     &           rnorm , temp , eps23
      save       cnorm , getv0, initv, update, ushift,
     &           rnorm , iter , eps23, kplusp, msglvl, nconv ,
     &           nevbef, nev0 , np0  , numcnv
c
c     %-----------------------%
c     | Local array arguments |
c     %-----------------------%
c
      integer    kp(4)
c
c     %--------------------------------------------------------------%
c     | Local workspace for the fresh real Schur decomposition of the|
c     | current KEV+NP upper Hessenberg H, computed just before      |
c     | applying shifts and passed into mydnapps. Sized by LDH, the  |
c     | same bound H itself uses (the actual size in use at any      |
c     | point is KPLUSP <= LDH).                                     |
c     %--------------------------------------------------------------%
c
      integer    arrowierr, arrownblk, arrowpick, arrowm, iconj
      integer    arrowblockstart(ldh), arrowblocklen(ldh),
     &           arrowconsumed(ldh), arrowiwork(1)
      logical    arrowselect(ldh)
      Double precision
     &           arrowwork(7*ldh),
     &           arrowsep1, arrowsep2,
     &           arroweig1r(ldh), arroweig1i(ldh), arrowbidx(ldh)
c
c     %--------------------------------------------------------------%
c     | ARROWSCHUR/ARROWSCHURVEC/ARROWWR/ARROWWI carry the fresh     |
c     | Schur decomposition of H from its computation site (right    |
c     | after the Arnoldi expansion, where the Ritz values/bounds    |
c     | are read off it) to the shift-application site, across the   |
c     | possible ISHIFT=0 reverse-communication exit in between, so  |
c     | they are SAVEd allocatables (re-sized if LDH ever changes    |
c     | between problems).                                           |
c     %--------------------------------------------------------------%
c
      Double precision, allocatable, save ::
     &           arrowschur(:,:), arrowschurvec(:,:),
     &           arrowwr(:), arrowwi(:), arrowbnd(:)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy  , dgetv0 , dnaitr , dnconv , dtrevc ,
     &           dngets , mydnapps_house , dvout  , ivout , arscnd,
     &           dlahqr , dtrsen , dsortc , dlacpy , dlaset,
     &           dgemv  , dscal
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           ddot , dnrm2 , dlapy2 , dlamch
      external   ddot , dnrm2 , dlapy2 , dlamch
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
c        %-------------------------------------%
c        | Get the machine dependent constant. |
c        %-------------------------------------%
c
         eps23 = dlamch ('Epsilon-Machine')
         eps23 = eps23**(2.0D+0  / 3.0D+0 )
c
         nev0   = nev
         np0    = np
c
c        %-------------------------------------%
c        | kplusp is the bound on the largest  |
c        |        Lanczos factorization built. |
c        | nconv is the current number of      |
c        |        "converged" eigenvlues.      |
c        | iter is the counter on the current  |
c        |      iteration step.                |
c        %-------------------------------------%
c
         kplusp = nev + np
         nconv  = 0
         iter   = 0
c
c        %---------------------------------------%
c        | Set flags for computing the first NEV |
c        | steps of the Arnoldi factorization.   |
c        %---------------------------------------%
c
         getv0    = .true.
         update   = .false.
         ushift   = .false.
         cnorm    = .false.
c
         if (info .ne. 0) then
c
c           %--------------------------------------------%
c           | User provides the initial residual vector. |
c           %--------------------------------------------%
c
            initv = .true.
            info  = 0
         else
            initv = .false.
         end if
      end if
c
c     %---------------------------------------------%
c     | Get a possibly random starting vector and   |
c     | force it into the range of the operator OP. |
c     %---------------------------------------------%
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
c
c           %-----------------------------------------%
c           | The initial vector is zero. Error exit. |
c           %-----------------------------------------%
c
            info = -9
            go to 1100
         end if
         getv0 = .false.
         ido  = 0
      end if
c
c     %-----------------------------------%
c     | Back from reverse communication : |
c     | continue with update step         |
c     %-----------------------------------%
c
      if (update) go to 20
c
c     %-------------------------------------------%
c     | Back from computing user specified shifts |
c     %-------------------------------------------%
c
      if (ushift) go to 50
c
c     %-------------------------------------%
c     | Back from computing residual norm   |
c     | at the end of the current iteration |
c     %-------------------------------------%
c
      if (cnorm)  go to 100
c
c     %----------------------------------------------------------%
c     | Compute the first NEV steps of the Arnoldi factorization |
c     %----------------------------------------------------------%
c
      call dnaitr  (ido, bmat, n, 0, nev, mode, resid, rnorm, v, ldv,
     &             h, ldh, ipntr, workd, info)
c
c     %---------------------------------------------------%
c     | ido .ne. 99 implies use of reverse communication  |
c     | to compute operations involving OP and possibly B |
c     %---------------------------------------------------%
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
     &           '_naup2: **** Start of major iteration number ****')
         end if
c
c        %-----------------------------------------------------------%
c        | Compute NP additional steps of the Arnoldi factorization. |
c        | Adjust NP since NEV might have been updated by last call  |
c        | to the shift application routine mydnapps_house.           |
c        %-----------------------------------------------------------%
c
         np  = kplusp - nev
c
         if (msglvl .gt. 1) then
            call ivout (logfil, 1, [nev], ndigit,
     &     '_naup2: The length of the current Arnoldi factorization')
            call ivout (logfil, 1, [np], ndigit,
     &           '_naup2: Extend the Arnoldi factorization by')
         end if
c
c        %-----------------------------------------------------------%
c        | Compute NP additional steps of the Arnoldi factorization. |
c        %-----------------------------------------------------------%
c
         ido = 0
   20    continue
         update = .true.
c
         call dnaitr  (ido  , bmat, n  , nev, np , mode , resid,
     &                rnorm, v   , ldv, h  , ldh, ipntr, workd,
     &                info)
c
c        %---------------------------------------------------%
c        | ido .ne. 99 implies use of reverse communication  |
c        | to compute operations involving OP and possibly B |
c        %---------------------------------------------------%
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
         if (msglvl .gt. 1) then
            call dvout  (logfil, 1, [rnorm], ndigit,
     &           '_naup2: Corresponding B-norm of the residual')
         end if
c
c        %--------------------------------------------------------%
c        | Compute the eigenvalues and corresponding error bounds |
c        | of the current upper Hessenberg matrix.                |
c        %--------------------------------------------------------%
c
         if (allocated(arrowschur)) then
            if (size(arrowschur,1) .ne. ldh) then
               deallocate (arrowschur, arrowschurvec, arrowwr,
     &                     arrowwi, arrowbnd)
            end if
         end if
         if (.not. allocated(arrowschur)) then
            allocate (arrowschur(ldh,ldh), arrowschurvec(ldh,ldh),
     &                arrowwr(ldh), arrowwi(ldh), arrowbnd(ldh))
         end if
c
c        %--------------------------------------------------------%
c        | DLAHQR (not the faster DHSEQR) is used deliberately:   |
c        | dneupd re-computes the Schur form of this same H with  |
c        | DLAHQR at extraction time and applies its positional   |
c        | SELECT/DTRSEN machinery assuming the Ritz arrays it    |
c        | was handed are in that SAME eigenvalue ordering --     |
c        | DHSEQR's multishift path can produce a different       |
c        | ordering, which silently breaks the extraction.        |
c        %--------------------------------------------------------%
c
         call dlacpy ('All', kplusp, kplusp, h, ldh, arrowschur, ldh)
         call dlaset ('All', kplusp, kplusp, zero, one,
     &                arrowschurvec, ldh)
         call dlahqr (.true., .true., kplusp, 1, kplusp,
     &                arrowschur, ldh, arrowwr, arrowwi,
     &                1, kplusp, arrowschurvec, ldh, arrowierr)
         if (arrowierr .ne. 0) then
            info = -8
            go to 1200
         end if
c
         call dcopy (kplusp, arrowwr, 1, ritzr, 1)
         call dcopy (kplusp, arrowwi, 1, ritzi, 1)
c
c        %--------------------------------------------------------%
c        | Eigenvectors of T (right, all) -> Q workspace. VL is   |
c        | not referenced for SIDE = 'R'; WORKL(KPLUSP+1:...) is  |
c        | free scratch at this point (its persistent copies are  |
c        | written just below).                                   |
c        %--------------------------------------------------------%
c
         call dtrevc ('R', 'A', arrowselect, kplusp, arrowschur, ldh,
     &                arrowwork, 1, q, ldq, kplusp, arrowm,
     &                workl(kplusp+1), arrowierr)
         if (arrowierr .ne. 0) then
            info = -8
            go to 1200
         end if
c
c        %--------------------------------------------------------%
c        | Scale the eigenvectors of T so their euclidean norms   |
c        | are one (DTREVC normalizes by largest component, with  |
c        | |(x,y)| = |x|+|y| for complex pairs) -- verbatim from  |
c        | dneigh.                                                |
c        %--------------------------------------------------------%
c
         iconj = 0
         do 37 j = 1, kplusp
            if ( abs( ritzi(j) ) .le. zero ) then
               temp = dnrm2( kplusp, q(1,j), 1 )
               call dscal ( kplusp, one / temp, q(1,j), 1 )
            else
               if (iconj .eq. 0) then
                  temp = dlapy2( dnrm2( kplusp, q(1,j), 1 ),
     &                           dnrm2( kplusp, q(1,j+1), 1 ) )
                  call dscal ( kplusp, one / temp, q(1,j), 1 )
                  call dscal ( kplusp, one / temp, q(1,j+1), 1 )
                  iconj = 1
               else
                  iconj = 0
               end if
            end if
   37    continue
c
c        %--------------------------------------------------------%
c        | Last components of H's eigenvectors = (last row of Z)  |
c        | times the eigenvectors of T; then the Ritz estimates   |
c        | are RNORM times their magnitudes -- verbatim from      |
c        | dneigh (complex pairs share one bound).                |
c        %--------------------------------------------------------%
c
         call dcopy (kplusp, arrowschurvec(kplusp,1), ldh, bounds, 1)
         call dgemv ('T', kplusp, kplusp, one, q, ldq, bounds, 1,
     &               zero, workl, 1)
c
         iconj = 0
         do 38 j = 1, kplusp
            if ( abs( ritzi(j) ) .le. zero ) then
               bounds(j) = rnorm * abs( workl(j) )
            else
               if (iconj .eq. 0) then
                  bounds(j) = rnorm * dlapy2( workl(j), workl(j+1) )
                  bounds(j+1) = bounds(j)
                  iconj = 1
               else
                  iconj = 0
               end if
            end if
   38    continue
c
c        %--------------------------------------------------------%
c        | Keep a copy of the bounds in DLAHQR (Schur)            |
c        | order -- dngets sorts BOUNDS together with RITZR/RITZI |
c        | below, and the exit path must hand dneupd the arrays   |
c        | in the original DLAHQR ordering (see the exit block).  |
c        %--------------------------------------------------------%
c
         call dcopy (kplusp, bounds, 1, arrowbnd, 1)
c
c        %----------------------------------------------------%
c        | Make a copy of eigenvalues and corresponding error |
c        | bounds obtained above.                             |
c        %----------------------------------------------------%
c
         call dcopy (kplusp, ritzr, 1, workl(kplusp**2+1), 1)
         call dcopy (kplusp, ritzi, 1, workl(kplusp**2+kplusp+1), 1)
         call dcopy (kplusp, bounds, 1, workl(kplusp**2+2*kplusp+1), 1)
c
c        %---------------------------------------------------%
c        | Select the wanted Ritz values and their bounds    |
c        | to be used in the convergence test.               |
c        | The wanted part of the spectrum and corresponding |
c        | error bounds are in the last NEV loc. of RITZR,   |
c        | RITZI and BOUNDS respectively. The variables NEV  |
c        | and NP may be updated if the NEV-th wanted Ritz   |
c        | value has a non zero imaginary part. In this case |
c        | NEV is increased by one and NP decreased by one.  |
c        %---------------------------------------------------%
c
         nev = nev0
         np = np0
         numcnv = nev
         call dngets  (ishift, which, nev, np, ritzr, ritzi,
     &                bounds, workl, workl(np+1))
         if (nev .eq. nev0+1) numcnv = nev0+1
c
c        %-------------------%
c        | Convergence test. |
c        %-------------------%
c
         call dcopy  (nev, bounds(np+1), 1, workl(2*np+1), 1)
         call dnconv  (nev, ritzr(np+1), ritzi(np+1), workl(2*np+1),
     &        tol, nconv)
c
         if (msglvl .gt. 2) then
            kp(1) = nev
            kp(2) = np
            kp(3) = numcnv
            kp(4) = nconv
            call ivout (logfil, 4, kp, ndigit,
     &                  '_naup2: NEV, NP, NUMCNV, NCONV are')
            call dvout  (logfil, kplusp, ritzr, ndigit,
     &           '_naup2: Real part of the eigenvalues of H')
            call dvout  (logfil, kplusp, ritzi, ndigit,
     &           '_naup2: Imaginary part of the eigenvalues of H')
            call dvout  (logfil, kplusp, bounds, ndigit,
     &          '_naup2: Ritz estimates of the current NCV Ritz values')
         end if
c
c        %---------------------------------------------------------%
c        | Count the number of unwanted Ritz values that have zero |
c        | Ritz estimates. If any Ritz estimates are equal to zero |
c        | then a leading block of H of order equal to at least    |
c        | the number of Ritz values with zero Ritz estimates has  |
c        | split off. None of these Ritz values may be removed by  |
c        | shifting. Decrease NP the number of shifts to apply. If |
c        | no shifts may be applied, then prepare to exit          |
c        %---------------------------------------------------------%
c
         nptemp = np
         do 30 j=1, nptemp
            if (bounds(j) .eq. zero) then
               np = np - 1
               nev = nev + 1
            end if
 30      continue
c
         if ( (nconv .ge. numcnv) .or.
     &        (iter .gt. mxiter) .or.
     &        (np .eq. 0) ) then
c
c           %------------------------------------------------%
c           | Prepare to exit. Put the converged Ritz values |
c           | and corresponding bounds in RITZ(1:NCONV) and  |
c           | BOUNDS(1:NCONV) respectively. Then sort. Be    |
c           | careful when NCONV > NP                        |
c           %------------------------------------------------%
c
c           %------------------------------------------%
c           |  Use h( 3,1 ) as storage to communicate  |
c           |  rnorm to _neupd if needed               |
c           %------------------------------------------%

            h(3,1) = rnorm
c
c           %----------------------------------------------------%
c           | Hand dneupd the Ritz values and bounds in the      |
c           | ORIGINAL DLAHQR (Schur) ordering of the final H,   |
c           | NOT sorted by WHICH. dneupd internally re-runs     |
c           | DLAHQR on this same H, marks its SELECT array at   |
c           | the POSITIONS jj of the converged wanted values in |
c           | the arrays we pass, and applies that SELECT to its |
c           | own Schur form via DTRSEN -- so the passed order   |
c           | must coincide with DLAHQR's, or the extraction     |
c           | reorders the wrong Schur blocks (returning         |
c           | unwanted eigenvalues and inconsistent vectors).    |
c           | The arrays were computed in exactly that order at  |
c           | the top of this iteration (ARROWWR/ARROWWI, plus   |
c           | the pristine bound copy ARROWBND) -- dngets only   |
c           | sorted the RITZR/RITZI/BOUNDS copies -- so simply  |
c           | restore them. dneupd does all wanted/converged     |
c           | selection itself from this data.                   |
c           %----------------------------------------------------%
c
            call dcopy (kplusp, arrowwr, 1, ritzr, 1)
            call dcopy (kplusp, arrowwi, 1, ritzi, 1)
            call dcopy (kplusp, arrowbnd, 1, bounds, 1)
c
            if (msglvl .gt. 1) then
               call dvout  (logfil, kplusp, ritzr, ndigit,
     &            '_naup2: Sorted real part of the eigenvalues')
               call dvout  (logfil, kplusp, ritzi, ndigit,
     &            '_naup2: Sorted imaginary part of the eigenvalues')
               call dvout  (logfil, kplusp, bounds, ndigit,
     &            '_naup2: Sorted ritz estimates.')
            end if
c
            if (iter .gt. mxiter .and. nconv .lt. numcnv) info = 1
c
            if (np .eq. 0 .and. nconv .lt. numcnv) info = 2
c
            np = nconv
            go to 1100
c
         else if ( (nconv .lt. numcnv) .and. (ishift .eq. 1) ) then
c
c           %-------------------------------------------------%
c           | Do not have all the requested eigenvalues yet.  |
c           | To prevent possible stagnation, adjust the size |
c           | of NEV.                                         |
c           %-------------------------------------------------%
c
            nevbef = nev
            nev = nev + min(nconv, np/2)
            if (nev .eq. 1 .and. kplusp .ge. 6) then
               nev = kplusp / 2
            else if (nev .eq. 1 .and. kplusp .gt. 3) then
               nev = 2
            end if
            if (nev .gt. kplusp - 2) then
               nev = kplusp - 2
            end if
c
            np = kplusp - nev
c
            if (nevbef .lt. nev)
     &         call dngets  (ishift, which, nev, np, ritzr, ritzi,
     &              bounds, workl, workl(np+1))
c
         end if
c
         if (msglvl .gt. 0) then
            call ivout (logfil, 1, [nconv], ndigit,
     &           '_naup2: no. of "converged" Ritz values at this iter.')
            if (msglvl .gt. 1) then
               kp(1) = nev
               kp(2) = np
               call ivout (logfil, 2, kp, ndigit,
     &              '_naup2: NEV and NP are')
               call dvout  (logfil, nev, ritzr(np+1), ndigit,
     &              '_naup2: "wanted" Ritz values -- real part')
               call dvout  (logfil, nev, ritzi(np+1), ndigit,
     &              '_naup2: "wanted" Ritz values -- imag part')
               call dvout  (logfil, nev, bounds(np+1), ndigit,
     &              '_naup2: Ritz estimates of the "wanted" values ')
            end if
         end if
c
         if (ishift .eq. 0) then
c
            ushift = .true.
            ido = 3
            go to 9000
         end if
c
   50    continue
c
         ushift = .false.
c
         if ( ishift .eq. 0 ) then
             call dcopy  (np, workl,       1, ritzr, 1)
             call dcopy  (np, workl(np+1), 1, ritzi, 1)
         end if
c
         if (msglvl .gt. 2) then
            call ivout (logfil, 1, [np], ndigit,
     &                  '_naup2: The number of shifts to apply ')
            call dvout  (logfil, np, ritzr, ndigit,
     &                  '_naup2: Real part of the shifts')
            call dvout  (logfil, np, ritzi, ndigit,
     &                  '_naup2: Imaginary part of the shifts')
            if ( ishift .eq. 1 )
     &          call dvout  (logfil, np, bounds, ndigit,
     &                  '_naup2: Ritz estimates of the shifts')
         end if
c
c        %----------------------------------------------------------%
c        | The full real Schur decomposition of the current KEV+NP  |
c        | upper Hessenberg H (T in ARROWSCHUR, Z in ARROWSCHURVEC, |
c        | eigenvalues in ARROWWR/ARROWWI, all still aligned with   |
c        | each other -- dngets sorted RITZR/RITZI, not these) was  |
c        | already computed by the single DHSEQR call at the top of |
c        | this iteration (see the OPTIMIZATION note there). H has  |
c        | not changed since, and DTREVC did not modify T or Z, so  |
c        | it is simply reused here -- no second decomposition is   |
c        | needed. The SAVE attribute on these arrays guarantees    |
c        | they survive the ISHIFT=0 reverse-communication exit     |
c        | above. (DTRSEN below reorders T and Z in place; that is  |
c        | fine, since they are recomputed fresh next iteration.)   |
c        %----------------------------------------------------------%
c
c        %----------------------------------------------------------%
c        | Classify each diagonal block (1x1 or 2x2) of the fresh   |
c        | Schur form. DHSEQR already returns each position's       |
c        | eigenvalue directly in ARROWWR/ARROWWI (both rows of a   |
c        | 2x2 block carry its complex-conjugate pair), so no extra |
c        | eigenvalue computation (e.g. DLANV2) is needed here.     |
c        %----------------------------------------------------------%
c
         arrownblk = 0
         i = 1
   55    continue
         if (i .le. kplusp) then
            arrownblk = arrownblk + 1
            arrowblockstart(arrownblk) = i
            if (i .lt. kplusp .and.
     &          abs(arrowschur(i+1,i)) .gt. zero) then
               arrowblocklen(arrownblk) = 2
               i = i + 2
            else
               arrowblocklen(arrownblk) = 1
               i = i + 1
            end if
            arroweig1r(arrownblk) = arrowwr(arrowblockstart(arrownblk))
            arroweig1i(arrownblk) =
     &         abs(arrowwi(arrowblockstart(arrownblk)))
            arrowbidx(arrownblk) = dble(arrownblk)
            go to 55
         end if
c
c        %----------------------------------------------------------%
c        | Select the desired blocks directly: sort the ARROWNBLK   |
c        | block-representative eigenvalues by the user's WHICH     |
c        | preference (largest real part, largest magnitude, etc.), |
c        | carrying each block's original index along in ARROWBIDX. |
c        | DSORTC sorts so the most-wanted blocks land at the END   |
c        | of the list; accumulate whole blocks from that end until |
c        | their lengths sum to at least NEV and mark them          |
c        | "consumed" (desired/kept).                               |
c        %----------------------------------------------------------%
c
         call dsortc (which, .true., arrownblk, arroweig1r, arroweig1i,
     &                arrowbidx)
c
         do 1210 j = 1, arrownblk
            arrowconsumed(j) = 0
 1210    continue
c
         arrowm = 0
         j = arrownblk
 1220    continue
         if (j .ge. 1 .and. arrowm .lt. nev) then
            arrowpick = nint(arrowbidx(j))
            arrowconsumed(arrowpick) = 1
            arrowm = arrowm + arrowblocklen(arrowpick)
            j = j - 1
            go to 1220
         end if
c
c        %----------------------------------------------------------%
c        | If the last block taken was a 2x2 block straddling the   |
c        | NEV boundary (a complex-conjugate pair that must not be  |
c        | split), ARROWM = NEV+1. Adopt the same policy dngets     |
c        | itself uses in that situation: grow NEV by one and       |
c        | shrink NP, rather than treating this as an error. The    |
c        | main loop re-derives NP from NEV at the top of the next  |
c        | iteration, so updating both here keeps them consistent.  |
c        %----------------------------------------------------------%
c
         if (arrowm .gt. nev) then
            nev = arrowm
            np  = kplusp - nev
         end if
c
c        %----------------------------------------------------------%
c        | Every index 1..KPLUSP belongs to exactly one block (the  |
c        | scan above partitions 1..KPLUSP with no gaps/overlaps),  |
c        | so ARROWSELECT can be set directly from ARROWCONSUMED in |
c        | one pass.                                                 |
c        %----------------------------------------------------------%
c
         do 1250 j = 1, arrownblk
            do 1245 i = 0, arrowblocklen(j)-1
               arrowselect(arrowblockstart(j)+i) =
     &            (arrowconsumed(j) .eq. 1)
 1245       continue
 1250    continue
c
c        %----------------------------------------------------------%
c        | Reorder the Schur form so the NEV selected ("desired")   |
c        | eigenvalues occupy the leading block, updating the Schur |
c        | vectors to match.                                        |
c        %----------------------------------------------------------%
c
         call dtrsen ('N', 'V', arrowselect, kplusp, arrowschur, ldh,
     &                arrowschurvec, ldh, arrowwr, arrowwi, arrowm,
     &                arrowsep1, arrowsep2, arrowwork, 7*ldh,
     &                arrowiwork, 1, arrowierr)
         if (arrowierr .ne. 0 .or. arrowm .ne. nev) then
            info = -8
            go to 1200
         end if
c
c        %-----------------------------------------------------------%
c        | Apply the nonsymmetric arrowhead restart: build the       |
c        | (NEV+1)x(NEV+1) arrowhead matrix from the NEV desired     |
c        | Schur eigenpairs, reduce it to upper Hessenberg form, and |
c        | read the updated H, V, and RESID off that result.         |
c        %-----------------------------------------------------------%
c
         call mydnapps_house (n, nev, np, v, ldv, h, ldh, resid, q, ldq,
     &        arrowschur, ldh, arrowschurvec, ldh, workd)
c
c        %---------------------------------------------%
c        | Compute the B-norm of the updated residual. |
c        | Keep B*RESID in WORKD(1:N) to be used in    |
c        | the first step of the next call to dnaitr . |
c        %---------------------------------------------%
c
         cnorm = .true.
         call arscnd (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy  (n, resid, 1, workd(n+1), 1)
            ipntr(1) = n + 1
            ipntr(2) = 1
            ido = 2
c
c           %----------------------------------%
c           | Exit in order to compute B*RESID |
c           %----------------------------------%
c
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy  (n, resid, 1, workd, 1)
         end if
c
  100    continue
c
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(1:N) := B*RESID            |
c        %----------------------------------%
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
         if (msglvl .gt. 2) then
            call dvout  (logfil, 1, [rnorm], ndigit,
     &      '_naup2: B-norm of residual for compressed factorization')
            call dmout  (logfil, nev, nev, h, ldh, ndigit,
     &        '_naup2: Compressed upper Hessenberg matrix H')
         end if
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
      nev = numcnv
c
 1200 continue
      ido = 99
c
c     %------------%
c     | Error Exit |
c     %------------%
c
      call arscnd (t1)
      tnaup2 = t1 - t0
c
 9000 continue
c
c     %-----------------------%
c     | End of mydnaup2_house |
c     %-----------------------%
c
      return
      end
