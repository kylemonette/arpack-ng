c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dhaup2
c
c\Description:
c  Intermediate level interface called by dhaupd.
c
c  Implements the thick-restarted Arnoldi iteration with harmonic Ritz
c  vectors (target zero). After each Arnoldi expansion, the harmonic
c  Ritz values, vectors, and residual norms of the current KEV+NP
c  factorization are computed from the harmonic pencil (see \Remarks).
c  These are used for the convergence test and for selecting the KEV
c  values to keep, exactly as dnaup2 uses Ritz values. The restart is
c  then applied by dhapps.
c
c  If the pencil becomes too ill-conditioned (see \Remarks), Ritz
c  values are used for the remainder of the run, and the restart is
c  applied by dqapps with Householder reflectors.
c
c\Usage:
c  call dhaup2
c     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
c       ISHIFT, MXITER, V, LDV, H, LDH, RITZR, RITZI, BOUNDS,
c       Q, LDQ, WORKL, IPNTR, WORKD, INFO, HARMFB )
c
c\Arguments
c  Identical to dnaup2 (see dnaup2.f), except that RITZR, RITZI and
c  BOUNDS hold harmonic Ritz values and residual norms, and ISHIFT must
c  be 1. Additionally:
c
c  HARMFB  Logical.  (OUTPUT)
c          .TRUE. if the iteration switched to Ritz values because the
c          harmonic pencil became too ill-conditioned.
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
c     dhapps   ARPACK harmonic thick restart routine.
c     dqapps   ARPACK arrowhead thick restart routine (Ritz values).
c     dnconv   ARPACK convergence of Ritz values routine.
c     dngets   ARPACK reorder Ritz values and error bounds routine.
c     dsortc   ARPACK sorting routine.
c     ivout    ARPACK utility routine that prints integers.
c     arscnd   ARPACK utility routine for timing.
c     dmout    ARPACK utility routine that prints matrices
c     dvout    ARPACK utility routine that prints vectors.
c     dlartg   LAPACK routine that generates a plane rotation.
c     dgghrd   LAPACK routine that reduces a matrix pencil to
c              Hessenberg-triangular form.
c     dhgeqz   LAPACK routine that computes the generalized real Schur
c              form of a Hessenberg-triangular pencil.
c     dtgevc   LAPACK routine that computes eigenvectors of a pencil in
c              generalized real Schur form.
c     dtgsen   LAPACK routine that reorders a generalized real Schur
c              form.
c     dlahqr   LAPACK routine that computes the real Schur form of an
c              upper Hessenberg matrix.
c     dtrevc   LAPACK routine that computes eigenvectors of a matrix in
c              real Schur form.
c     dtrsen   LAPACK routine that reorders a real Schur form.
c     dlamch   LAPACK routine that determines machine constants.
c     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dlacpy   LAPACK matrix copy routine.
c     dlaset   LAPACK matrix initialization routine.
c     dtrmm    Level 3 BLAS routine for triangular matrix times matrix.
c     drot     Level 1 BLAS that applies a plane rotation.
c     daxpy    Level 1 BLAS that computes a vector triad.
c     dcopy    Level 1 BLAS that copies one vector to another.
c     ddot     Level 1 BLAS that computes the scalar product of two
c              vectors.
c     dnrm2    Level 1 BLAS that computes the norm of a vector.
c     dscal    Level 1 BLAS that scales a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Remarks
c  1. Let m = KEV+NP and HBAR = [H; RNORM*e_m^T] = QH*[R; 0].
c     A harmonic Ritz pair (theta,y) satisfies
c         HBAR^T*(HBAR - theta*[I; 0])*y = 0,
c     which is equivalent to R*y = theta*B^T*y with B = QH(1:m,1:m).
c     The pencil (B^T, R) is reduced to generalized real Schur form,
c     with eigenvalues mu = 1/theta, and right Schur vectors Z. The
c     last column U of QH spans the null space of HBAR^T.
c  2. The residual norm of the harmonic Ritz pair (theta, V*y), with
c     ||y|| = 1, is ||HBAR*y - theta*[y; 0]||, computed directly.
c  3. The smallest singular value of B is |U(m+1)|. If it falls below
c     sqrt(eps), H is numerically singular and the harmonic values are
c     unreliable, so Ritz values are used from then on.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dhaup2
     &   ( ido, bmat, n, which, nev, np, tol, resid, mode, iupd,
     &     ishift, mxiter, v, ldv, h, ldh, ritzr, ritzi, bounds,
     &     q, ldq, workl, ipntr, workd, info, harmfb )
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
      logical    harmfb
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
      logical    cnorm , getv0, initv, update
      integer    ierr  , iter , j    , kplusp, msglvl, nconv,
     &           nevbef, nev0 , np0  , nptemp, numcnv
      Double precision
     &           rnorm , temp , eps23, harmtol
      save       cnorm , getv0, initv, update,
     &           rnorm , iter , eps23, kplusp, msglvl, nconv ,
     &           nevbef, nev0 , np0  , numcnv, harmtol
c
c     %-----------------------%
c     | Local array arguments |
c     %-----------------------%
c
      integer    kp(4)
c
c     %--------------------------------------------------------------%
c     | Workspace for the Schur forms and block selection, sized by  |
c     | LDH (the size in use is KPLUSP <= LDH).                      |
c     %--------------------------------------------------------------%
c
      integer    i, iconj, nblk, pick, msel, lwork
      integer    blockstart(ldh), blocklen(ldh), consumed(ldh),
     &           iwork(1)
      logical    select(ldh)
      Double precision
     &           work(8*ldh+16), eig1r(ldh), eig1i(ldh), bidx(ldh),
     &           dif(2), qdum(1), sep1, sep2, pl, pr, cs, sn, rtmp, den
c
c     %---------------------------------------------------------------%
c     | HS/HP:   Schur form (or the pencil's generalized Schur form). |
c     | HZ:      right Schur vectors.                                 |
c     | HQH:     orthogonal factor QH of HBAR (last column is U).     |
c     | HW:      HBAR, then HBAR times the eigenvectors.              |
c     | HWR/HWI: eigenvalues in Schur order.                          |
c     | HALR/HALI/HBET: generalized eigenvalues of the pencil.        |
c     | SAVEd, since they are used across the B*x reverse             |
c     | communication exit at the end of an iteration.                |
c     %---------------------------------------------------------------%
c
      Double precision, allocatable, save ::
     &           hs(:,:), hp(:,:), hz(:,:), hqh(:,:), hw(:,:),
     &           hwr(:), hwi(:), halr(:), hali(:), hbet(:)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy  , dgetv0 , dnaitr , dnconv , dngets ,
     &           dhapps , dqapps , dvout  , ivout  , arscnd ,
     &           dmout  , dsortc , dlartg , dgghrd , dhgeqz ,
     &           dtgevc , dtgsen , dlahqr , dtrevc , dtrsen ,
     &           dlacpy , dlaset , dtrmm  , drot   ,
     &           daxpy  , dscal
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
         harmtol = sqrt(eps23)
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
         harmfb = .false.
c
c        %---------------------------------------%
c        | Set flags for computing the first NEV |
c        | steps of the Arnoldi factorization.   |
c        %---------------------------------------%
c
         getv0    = .true.
         update   = .false.
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
c
         if (allocated(hs)) then
            if (size(hs,1) .ne. ldh) then
               deallocate (hs, hp, hz, hqh, hw, hwr, hwi, halr, hali,
     &                     hbet)
            end if
         end if
         if (.not. allocated(hs)) then
            allocate (hs(ldh,ldh), hp(ldh,ldh), hz(ldh,ldh),
     &                hqh(ldh+1,ldh+1), hw(ldh+1,ldh), hwr(ldh),
     &                hwi(ldh), halr(ldh), hali(ldh), hbet(ldh))
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
c     |           Each iteration thick restarts the Arnoldi          |
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
c        %-----------------------------------------------------------%
c        | Compute NP additional steps of the Arnoldi factorization. |
c        | Adjust NP since NEV might have been updated by the last   |
c        | restart.                                                  |
c        %-----------------------------------------------------------%
c
         np  = kplusp - nev
c
         if (msglvl .gt. 1) then
            call ivout (logfil, 1, [nev], ndigit,
     &     '_haup2: The length of the current Arnoldi factorization')
            call ivout (logfil, 1, [np], ndigit,
     &           '_haup2: Extend the Arnoldi factorization by')
         end if
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
     &           '_haup2: Corresponding B-norm of the residual')
         end if
c
         lwork = 8*ldh + 16
c
         if (.not. harmfb) then
c
c           %-------------------------------------------------%
c           | Givens QR of HBAR = [H; RNORM*e_m^T] in HW,     |
c           | accumulating QH explicitly in HQH.              |
c           %-------------------------------------------------%
c
            call dlacpy ('All', kplusp, kplusp, h, ldh, hw, ldh+1)
            do 21 j = 1, kplusp
               hw(kplusp+1,j) = zero
   21       continue
            hw(kplusp+1,kplusp) = rnorm
            call dlaset ('All', kplusp+1, kplusp+1, zero, one, hqh,
     &                   ldh+1)
            do 22 j = 1, kplusp
               call dlartg (hw(j,j), hw(j+1,j), cs, sn, rtmp)
               hw(j,j) = rtmp
               hw(j+1,j) = zero
               if (j .lt. kplusp) then
                  call drot (kplusp-j, hw(j,j+1), ldh+1, hw(j+1,j+1),
     &                       ldh+1, cs, sn)
               end if
               call drot (j+1, hqh(1,j), 1, hqh(1,j+1), 1, cs, sn)
   22       continue
c
c           %------------------------------------------------%
c           | Switch to Ritz values if |U(m+1)| is too small |
c           | (\Remarks 3).                                  |
c           %------------------------------------------------%
c
            if (abs(hqh(kplusp+1,kplusp+1)) .lt. harmtol) then
               harmfb = .true.
               if (msglvl .gt. 0) then
                  call dvout (logfil, 1, [hqh(kplusp+1,kplusp+1)],
     &                 ndigit, '_haup2: |U(m+1)| too small, using Ritz')
               end if
            end if
         end if
c
         if (.not. harmfb) then
c
c           %-------------------------------------------------%
c           | Generalized real Schur form (HS,HP) of the      |
c           | pencil (B^T, R), with right Schur vectors HZ.   |
c           %-------------------------------------------------%
c
            do 24 j = 1, kplusp
               do 23 i = 1, kplusp
                  hs(i,j) = hqh(j,i)
   23          continue
   24       continue
            call dlaset ('All', kplusp, kplusp, zero, zero, hp, ldh)
            call dlacpy ('Upper', kplusp, kplusp, hw, ldh+1, hp, ldh)
c
            call dgghrd ('N', 'I', kplusp, 1, kplusp, hs, ldh, hp, ldh,
     &                   qdum, 1, hz, ldh, ierr)
            call dhgeqz ('S', 'N', 'V', kplusp, 1, kplusp, hs, ldh,
     &                   hp, ldh, halr, hali, hbet, qdum, 1, hz, ldh,
     &                   work, lwork, ierr)
            if (ierr .ne. 0) then
               info = -8
               go to 1200
            end if
c
c           %------------------------------------------------%
c           | theta = 1/mu = BETA/(ALPHAR + i*ALPHAI). The    |
c           | denominator is nonzero since |U(m+1)| > 0.      |
c           %------------------------------------------------%
c
            do 25 j = 1, kplusp
               den = halr(j)**2 + hali(j)**2
               ritzr(j) = hbet(j) * halr(j) / den
               ritzi(j) = -hbet(j) * hali(j) / den
   25       continue
c
c           %------------------------------------------------%
c           | Eigenvectors of the pencil, back transformed by |
c           | HZ, stored in Q.                                |
c           %------------------------------------------------%
c
            call dlacpy ('All', kplusp, kplusp, hz, ldh, q, ldq)
            call dtgevc ('R', 'B', select, kplusp, hs, ldh, hp, ldh,
     &                   qdum, 1, q, ldq, kplusp, msel, work, ierr)
            if (ierr .ne. 0) then
               info = -8
               go to 1200
            end if
         else
c
c           %-------------------------------------------------%
c           | Real Schur form HS = HZ^T*H*HZ, and the         |
c           | eigenvectors of H stored in Q.                  |
c           %-------------------------------------------------%
c
            call dlacpy ('All', kplusp, kplusp, h, ldh, hs, ldh)
            call dlaset ('All', kplusp, kplusp, zero, one, hz, ldh)
            call dlahqr (.true., .true., kplusp, 1, kplusp, hs, ldh,
     &                   ritzr, ritzi, 1, kplusp, hz, ldh, ierr)
            if (ierr .ne. 0) then
               info = -8
               go to 1200
            end if
c
            call dlacpy ('All', kplusp, kplusp, hz, ldh, q, ldq)
            call dtrevc ('R', 'B', select, kplusp, hs, ldh, qdum, 1,
     &                   q, ldq, kplusp, msel, work, ierr)
            if (ierr .ne. 0) then
               info = -8
               go to 1200
            end if
         end if
c
c        %--------------------------------------------------%
c        | Scale the eigenvectors to unit euclidean norm,   |
c        | where a complex pair is stored as (real, imag)   |
c        | in consecutive columns.                          |
c        %--------------------------------------------------%
c
         iconj = 0
         do 26 j = 1, kplusp
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
   26    continue
c
         if (.not. harmfb) then
c
c           %------------------------------------------------%
c           | Residual norms ||HBAR*y - theta*[y; 0]||        |
c           | (\Remarks 2), with HW = HBAR*Q formed from the  |
c           | upper triangle of H and then its subdiagonal.   |
c           %------------------------------------------------%
c
            call dlacpy ('All', kplusp, kplusp, q, ldq, hw, ldh+1)
            call dtrmm ('Left', 'Upper', 'No transpose', 'Non-unit',
     &                  kplusp, kplusp, one, h, ldh, hw, ldh+1)
            do 27 i = 1, kplusp-1
               call daxpy (kplusp, h(i+1,i), q(i,1), ldq, hw(i+1,1),
     &                     ldh+1)
   27       continue
            call dcopy (kplusp, q(kplusp,1), ldq, hw(kplusp+1,1), ldh+1)
            call dscal (kplusp, rnorm, hw(kplusp+1,1), ldh+1)
c
            iconj = 0
            do 28 j = 1, kplusp
               if ( abs( ritzi(j) ) .le. zero ) then
                  call daxpy (kplusp, -ritzr(j), q(1,j), 1, hw(1,j), 1)
                  bounds(j) = dnrm2 (kplusp+1, hw(1,j), 1)
               else
                  if (iconj .eq. 0) then
                     call daxpy (kplusp, -ritzr(j), q(1,j), 1,
     &                           hw(1,j), 1)
                     call daxpy (kplusp, ritzi(j), q(1,j+1), 1,
     &                           hw(1,j), 1)
                     call daxpy (kplusp, -ritzr(j), q(1,j+1), 1,
     &                           hw(1,j+1), 1)
                     call daxpy (kplusp, -ritzi(j), q(1,j), 1,
     &                           hw(1,j+1), 1)
                     bounds(j) = dlapy2( dnrm2(kplusp+1, hw(1,j), 1),
     &                                   dnrm2(kplusp+1, hw(1,j+1), 1) )
                     bounds(j+1) = bounds(j)
                     iconj = 1
                  else
                     iconj = 0
                  end if
               end if
   28       continue
         else
c
c           %------------------------------------------------%
c           | Ritz estimates RNORM*|e_m^T*y|.                |
c           %------------------------------------------------%
c
            call dcopy (kplusp, q(kplusp,1), ldq, workl, 1)
            iconj = 0
            do 29 j = 1, kplusp
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
   29       continue
         end if
c
c        %------------------------------------------------%
c        | Eigenvalues in Schur order, used to select the |
c        | blocks kept by the restart.                    |
c        %------------------------------------------------%
c
         call dcopy (kplusp, ritzr, 1, hwr, 1)
         call dcopy (kplusp, ritzi, 1, hwi, 1)
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
     &                  '_haup2: NEV, NP, NUMCNV, NCONV are')
            call dvout  (logfil, kplusp, ritzr, ndigit,
     &           '_haup2: Real part of the harmonic Ritz values')
            call dvout  (logfil, kplusp, ritzi, ndigit,
     &           '_haup2: Imaginary part of the harmonic Ritz values')
            call dvout  (logfil, kplusp, bounds, ndigit,
     &           '_haup2: Residual norms of the harmonic Ritz pairs')
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
c
            h(3,1) = rnorm
c
c           %----------------------------------------------%
c           | To be consistent with dngets, we first do a  |
c           | pre-processing sort in order to keep complex |
c           | conjugate pairs together.  This is similar   |
c           | to the pre-processing sort used in dngets    |
c           | except that the sort is done in the opposite |
c           | order.                                       |
c           %----------------------------------------------%
c
            if (which .eq. 'LM') wprime = 'SR'
            if (which .eq. 'SM') wprime = 'LR'
            if (which .eq. 'LR') wprime = 'SM'
            if (which .eq. 'SR') wprime = 'LM'
            if (which .eq. 'LI') wprime = 'SM'
            if (which .eq. 'SI') wprime = 'LM'
c
            call dsortc  (wprime, .true., kplusp, ritzr, ritzi, bounds)
c
c           %----------------------------------------------%
c           | Now sort Ritz values so that converged Ritz  |
c           | values appear within the first NEV locations |
c           | of ritzr, ritzi and bounds, and the most     |
c           | desired one appears at the front.            |
c           %----------------------------------------------%
c
            if (which .eq. 'LM') wprime = 'SM'
            if (which .eq. 'SM') wprime = 'LM'
            if (which .eq. 'LR') wprime = 'SR'
            if (which .eq. 'SR') wprime = 'LR'
            if (which .eq. 'LI') wprime = 'SI'
            if (which .eq. 'SI') wprime = 'LI'
c
            call dsortc (wprime, .true., kplusp, ritzr, ritzi, bounds)
c
c           %--------------------------------------------------%
c           | Scale the Ritz estimate of each Ritz value       |
c           | by 1 / max(eps23,magnitude of the Ritz value).   |
c           %--------------------------------------------------%
c
            do 35 j = 1, numcnv
                temp = max(eps23,dlapy2 (ritzr(j),
     &                                   ritzi(j)))
                bounds(j) = bounds(j)/temp
 35         continue
c
c           %----------------------------------------------------%
c           | Sort the Ritz values according to the scaled Ritz  |
c           | estimates.  This will push all the converged ones  |
c           | towards the front of ritzr, ritzi, bounds          |
c           | (in the case when NCONV < NEV.)                    |
c           %----------------------------------------------------%
c
            wprime = 'LR'
            call dsortc (wprime, .true., numcnv, bounds, ritzr, ritzi)
c
c           %----------------------------------------------%
c           | Scale the Ritz estimate back to its original |
c           | value.                                       |
c           %----------------------------------------------%
c
            do 40 j = 1, numcnv
                temp = max(eps23, dlapy2 (ritzr(j),
     &                                   ritzi(j)))
                bounds(j) = bounds(j)*temp
 40         continue
c
c           %------------------------------------------------%
c           | Sort the converged Ritz values again so that   |
c           | the "threshold" value appears at the front of  |
c           | ritzr, ritzi and bound.                        |
c           %------------------------------------------------%
c
            call dsortc (which, .true., nconv, ritzr, ritzi, bounds)
c
            if (msglvl .gt. 1) then
               call dvout  (logfil, kplusp, ritzr, ndigit,
     &            '_haup2: Sorted real part of the eigenvalues')
               call dvout  (logfil, kplusp, ritzi, ndigit,
     &            '_haup2: Sorted imaginary part of the eigenvalues')
               call dvout  (logfil, kplusp, bounds, ndigit,
     &            '_haup2: Sorted ritz estimates.')
            end if
c
c           %------------------------------------%
c           | Max iterations have been exceeded. |
c           %------------------------------------%
c
            if (iter .gt. mxiter .and. nconv .lt. numcnv) info = 1
c
c           %---------------------%
c           | No shifts to apply. |
c           %---------------------%
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
c           %---------------------------------------%
c           | If the size of NEV was just increased |
c           | resort the eigenvalues.               |
c           %---------------------------------------%
c
            if (nevbef .lt. nev)
     &         call dngets  (ishift, which, nev, np, ritzr, ritzi,
     &              bounds, workl, workl(np+1))
c
         end if
c
         if (msglvl .gt. 0) then
            call ivout (logfil, 1, [nconv], ndigit,
     &           '_haup2: no. of "converged" Ritz values at this iter.')
            if (msglvl .gt. 1) then
               kp(1) = nev
               kp(2) = np
               call ivout (logfil, 2, kp, ndigit,
     &              '_haup2: NEV and NP are')
               call dvout  (logfil, nev, ritzr(np+1), ndigit,
     &              '_haup2: "wanted" Ritz values -- real part')
               call dvout  (logfil, nev, ritzi(np+1), ndigit,
     &              '_haup2: "wanted" Ritz values -- imag part')
               call dvout  (logfil, nev, bounds(np+1), ndigit,
     &              '_haup2: Ritz estimates of the "wanted" values ')
            end if
         end if
c
c        %----------------------------------------------------------%
c        | Partition the Schur form into its 1x1 and 2x2 diagonal   |
c        | blocks, each represented by its eigenvalue.              |
c        %----------------------------------------------------------%
c
         nblk = 0
         i = 1
   55    continue
         if (i .le. kplusp) then
            nblk = nblk + 1
            blockstart(nblk) = i
            if (i .lt. kplusp .and. abs(hs(i+1,i)) .gt. zero) then
               blocklen(nblk) = 2
               i = i + 2
            else
               blocklen(nblk) = 1
               i = i + 1
            end if
            eig1r(nblk) = hwr(blockstart(nblk))
            eig1i(nblk) = abs(hwi(blockstart(nblk)))
            bidx(nblk) = dble(nblk)
            go to 55
         end if
c
c        %----------------------------------------------------------%
c        | Sort the blocks by WHICH (most wanted at the end) and    |
c        | keep whole blocks from the end until at least NEV        |
c        | values are kept.                                         |
c        %----------------------------------------------------------%
c
         call dsortc (which, .true., nblk, eig1r, eig1i, bidx)
c
         do 60 j = 1, nblk
            consumed(j) = 0
   60    continue
c
         msel = 0
         j = nblk
   65    continue
         if (j .ge. 1 .and. msel .lt. nev) then
            pick = nint(bidx(j))
            consumed(pick) = 1
            msel = msel + blocklen(pick)
            j = j - 1
            go to 65
         end if
c
c        %----------------------------------------------------------%
c        | Keep a complex conjugate pair together at the NEV        |
c        | boundary by increasing NEV by one.                       |
c        %----------------------------------------------------------%
c
         if (msel .gt. nev) then
            nev = msel
            np  = kplusp - nev
         end if
c
         do 75 j = 1, nblk
            do 70 i = 0, blocklen(j)-1
               select(blockstart(j)+i) = (consumed(j) .eq. 1)
   70       continue
   75    continue
c
         if (.not. harmfb) then
c
c           %------------------------------------------------%
c           | Reorder the generalized Schur form so the NEV   |
c           | selected values lead, and apply the harmonic    |
c           | thick restart.                                  |
c           %------------------------------------------------%
c
            call dtgsen (0, .false., .true., select, kplusp, hs, ldh,
     &                   hp, ldh, halr, hali, hbet, qdum, 1, hz, ldh,
     &                   msel, pl, pr, dif, work, lwork, iwork, 1,
     &                   ierr)
            if (ierr .ne. 0 .or. msel .ne. nev) then
               info = -8
               go to 1200
            end if
c
            call dhapps (n, nev, np, v, ldv, h, ldh, resid, rnorm,
     &                   hz, ldh, hqh(1,kplusp+1))
         else
c
c           %------------------------------------------------%
c           | Reorder the Schur form so the NEV selected      |
c           | values lead, and apply the arrowhead restart.   |
c           %------------------------------------------------%
c
            call dtrsen ('N', 'V', select, kplusp, hs, ldh, hz, ldh,
     &                   hwr, hwi, msel, sep1, sep2, work, lwork,
     &                   iwork, 1, ierr)
            if (ierr .ne. 0 .or. msel .ne. nev) then
               info = -8
               go to 1200
            end if
c
            call dqapps (n, nev, np, v, ldv, h, ldh, resid, q, ldq,
     &                   hs, ldh, hz, ldh, workd, .true.)
         end if
c
c        %---------------------------------------------%
c        | Compute the B-norm of the updated residual. |
c        | Keep B*RESID in WORKD(1:N) to be used in    |
c        | the first step of the next call to dnaitr.  |
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
     &      '_haup2: B-norm of residual for compressed factorization')
            call dmout  (logfil, nev, nev, h, ldh, ndigit,
     &        '_haup2: Compressed upper Hessenberg matrix H')
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
c     %---------------%
c     | End of dhaup2 |
c     %---------------%
c
      return
      end
