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
c  Ritz values and residual norms of the current KEV+NP factorization
c  are computed from the real Schur form of an (KEV+NP+1) by
c  (KEV+NP+1) upper Hessenberg matrix HP (see \Remarks). These are used
c  for the convergence test and for selecting the KEV values to keep,
c  exactly as dnaup2 uses Ritz values. The restart is then applied by
c  dhapps.
c
c  If H becomes numerically singular (see \Remarks), Ritz values are
c  used for the remainder of the run, and the restart is applied by
c  dqapps with Householder reflectors.
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
c          .TRUE. if the iteration switched to Ritz values because H
c          became numerically singular, or the Schur form of HP could
c          not be computed or reordered (see \Remarks 3).
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
c  3. Krylov-Schur and Implicit Restarting: A Bordered Hessenberg
c     Perspective, James Baglama, Kyle Monette, Vasilije Perovic,
c     (2026).
c  4. A Unified View of Arrowhead Matrix Transformations and Lanczos
c     Restarts, James Baglama, Kyle Monette, Vasilije Perovic, (2026).
c
c\Routines called:
c     dgetv0   ARPACK initial vector generation routine.
c     dnaitr   ARPACK Arnoldi factorization routine.
c     dhapps   ARPACK harmonic thick restart routine.
c     dqapps   ARPACK bordered thick restart routine (Ritz values).
c     dnconv   ARPACK convergence of Ritz values routine.
c     dngets   ARPACK reorder Ritz values and error bounds routine.
c     dsortc   ARPACK sorting routine.
c     ivout    ARPACK utility routine that prints integers.
c     arscnd   ARPACK utility routine for timing.
c     dmout    ARPACK utility routine that prints matrices
c     dvout    ARPACK utility routine that prints vectors.
c     dlartg   LAPACK routine that generates a plane rotation.
c     dlahqr   LAPACK routine that computes the real Schur form of an
c              upper Hessenberg matrix.
c     dtrevc   LAPACK routine that computes eigenvectors of a matrix in
c              real Schur form.
c     dtrsen   LAPACK routine that reorders a real Schur form.
c     dlamch   LAPACK routine that determines machine constants.
c     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dlacpy   LAPACK matrix copy routine.
c     dlaset   LAPACK matrix initialization routine.
c     dtrmv    Level 2 BLAS routine for triangular matrix times vector.
c     dgemv    Level 2 BLAS routine for matrix vector multiplication.
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
c  1. Let m = KEV+NP and HBAR = [H; RNORM*e_m^T] = QH*[R; 0], where
c     QH = G_1*...*G_m is a product of plane rotations, G_j acting on
c     coordinates j and j+1. The last column U of QH spans the null
c     space of HBAR^T. With E = -U(1:m)/U(m+1), the (m+1) by (m+1)
c     upper Hessenberg matrix
c         HP = [HBAR, HBAR*E]
c     has as eigenvalues the harmonic Ritz values, i.e., the eigenvalues
c     of G = H + RNORM**2*H^{-T}*e_m*e_m^T, together with zero, whose
c     eigenvector is U. The restart keeps the NEV wanted harmonic Ritz
c     values and zero (see dhapps).
c  2. A unit eigenvector v of HP with eigenvalue theta /= 0 is
c     orthogonal to U and gives the harmonic Ritz vector y = [I, E]*v.
c     With r = |v(m+1)|/|U(m+1)|, ||y|| = sqrt(1+r^2) and the residual
c     norm of the harmonic Ritz pair (theta, V*y/||y||) is
c     |theta|*r/sqrt(1+r^2). Since U is both a left and a right
c     eigenvector of HP for zero, it is a Schur vector of HP. The zero
c     eigenvalue is therefore at the position of the largest entry of
c     HZ^T*U, which is identified this way rather than by magnitude,
c     since a harmonic Ritz value may also be close to zero.
c  3. Since 1/U(m+1)**2 = 1 + RNORM**2*||H^{-T}*e_m||**2, if |U(m+1)|
c     falls below sqrt(eps), H is numerically singular and the harmonic
c     values are unreliable, so Ritz values are used from then on.
c     Ritz values are also used from then on if DLAHQR, DTREVC, or
c     DTRSEN fails on HP; the current iteration is then redone with
c     Ritz values.
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
c     | LDH+1 (the size in use is at most KPLUSP+1 <= LDH+1). The    |
c     | arrays are SAVEd allocatables, allocated with HS below, so   |
c     | that they are not allocated and freed on every reverse-      |
c     | communication call.                                          |
c     %--------------------------------------------------------------%
c
      integer    i, iconj, nblk, pick, msel, lwork, mp, ntot, izero
      integer    iwork(1)
      integer, allocatable, save ::
     &           blockstart(:), blocklen(:), consumed(:)
      logical, allocatable, save ::
     &           select(:)
      Double precision
     &           qdum(1), sep1, sep2, rtmp, uinv
      Double precision, allocatable, save ::
     &           work(:), eig1r(:), eig1i(:), bidx(:),
     &           rc(:), rs(:)
c
c     %------------------------------------------------------------%
c     | HS:      Schur form of HP (or of H for Ritz values).       |
c     | HZ:      Schur vectors.                                    |
c     | HV:      eigenvectors of HS (harmonic Ritz values).        |
c     | HU:      last column U of the orthogonal factor QH of HBAR. |
c     | HWR/HWI: eigenvalues in Schur order.                       |
c     | SAVEd, since they are used across the B*x reverse          |
c     | communication exit at the end of an iteration.             |
c     %------------------------------------------------------------%
c
      Double precision, allocatable, save ::
     &           hs(:,:), hz(:,:), hv(:,:), hu(:), hwr(:), hwi(:)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy  , dgetv0 , dnaitr , dnconv , dngets ,
     &           dhapps , dqapps , dvout  , ivout  , arscnd ,
     &           dmout  , dsortc , dlartg , dlahqr , dtrevc ,
     &           dtrsen , dlacpy , dlaset , dtrmv  , daxpy  ,
     &           dscal  , dgemv
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      integer    idamax
      Double precision
     &           ddot , dnrm2 , dlapy2 , dlamch
      external   ddot , dnrm2 , dlapy2 , dlamch , idamax
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
            if (size(hs,1) .ne. ldh+1) then
               deallocate (hs, hz, hv, hu, hwr, hwi,
     &                     blockstart, blocklen, consumed, select,
     &                     work, eig1r, eig1i, bidx, rc, rs)
            end if
         end if
         if (.not. allocated(hs)) then
            allocate (hs(ldh+1,ldh+1), hz(ldh+1,ldh+1),
     &                hv(ldh+1,ldh+1), hu(ldh+1), hwr(ldh+1),
     &                hwi(ldh+1), blockstart(ldh+1), blocklen(ldh+1),
     &                consumed(ldh+1), select(ldh+1), work(8*ldh+16),
     &                eig1r(ldh+1), eig1i(ldh+1), bidx(ldh+1),
     &                rc(ldh), rs(ldh))
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
c        %--------------------------------------------------%
c        | The computation below restarts here with Ritz    |
c        | values if a harmonic computation fails.          |
c        %--------------------------------------------------%
c
   31    continue
c
         if (.not. harmfb) then
c
c           %--------------------------------------------------%
c           | Givens QR of HBAR = [H; RNORM*e_m^T], with the    |
c           | rotations stored in RC and RS. Since R is not     |
c           | needed, only the row of HBAR being reduced is     |
c           | kept, in WORK: rotation j combines WORK(j:m) with |
c           | row j+1 of H and keeps the second row.            |
c           %--------------------------------------------------%
c
            call dcopy (kplusp, h, ldh, work, 1)
            do 22 j = 1, kplusp-1
               call dlartg (work(j), h(j+1,j), rc(j), rs(j), rtmp)
               call dscal (kplusp-j, -rs(j), work(j+1), 1)
               call daxpy (kplusp-j, rc(j), h(j+1,j+1), ldh,
     &                     work(j+1), 1)
   22       continue
            call dlartg (work(kplusp), rnorm, rc(kplusp), rs(kplusp),
     &                   rtmp)
c
c           %------------------------------------------------%
c           | Switch to Ritz values if |U(m+1)| = |RC(m)| is |
c           | too small (\Remarks 3).                        |
c           %------------------------------------------------%
c
            if (abs(rc(kplusp)) .lt. harmtol) then
               harmfb = .true.
               if (msglvl .gt. 0) then
                  call dvout (logfil, 1, [rc(kplusp)],
     &                 ndigit, '_haup2: |U(m+1)| too small, using Ritz')
               end if
            end if
         end if
c
         if (.not. harmfb) then
c
c           %--------------------------------------------------%
c           | U = QH*e_{m+1} from the rotations (\Remarks 1).  |
c           %--------------------------------------------------%
c
            mp = kplusp + 1
            hu(mp) = rc(kplusp)
            temp = -rs(kplusp)
            do 23 j = kplusp-1, 1, -1
               hu(j+1) = rc(j) * temp
               temp = -rs(j) * temp
   23       continue
            hu(1) = temp
            uinv = one / hu(mp)
c
c           %--------------------------------------------------%
c           | HP = [HBAR, HBAR*E] in HS, where                 |
c           | HBAR*E = -HBAR*U(1:m)/U(m+1) is formed from the  |
c           | upper triangle of H and then its subdiagonal.    |
c           | Entries below the subdiagonal are zero.          |
c           %--------------------------------------------------%
c
            call dlaset ('Lower', mp-2, mp-2, zero, zero, hs(3,1),
     &                   ldh+1)
            call dlacpy ('Upper', kplusp, kplusp, h, ldh, hs, ldh+1)
            call dcopy (kplusp, hu, 1, hs(1,mp), 1)
            call dtrmv ('Upper', 'No transpose', 'Non-unit', kplusp,
     &                  h, ldh, hs(1,mp), 1)
            do 25 j = 1, kplusp-1
               hs(j+1,j) = h(j+1,j)
               hs(j+1,mp) = hs(j+1,mp) + h(j+1,j) * hu(j)
   25       continue
            hs(mp,kplusp) = rnorm
            hs(mp,mp) = rnorm * hu(kplusp)
            call dscal (mp, -uinv, hs(1,mp), 1)
c
c           %--------------------------------------------------%
c           | Real Schur form HS = HZ^T*HP*HZ, eigenvalues in  |
c           | HWR and HWI, and eigenvectors of HS in HV. If    |
c           | DLAHQR or DTREVC fails, switch to Ritz values.   |
c           %--------------------------------------------------%
c
            call dlaset ('All', mp, mp, zero, one, hz, ldh+1)
            call dlahqr (.true., .true., mp, 1, mp, hs, ldh+1, hwr, hwi,
     &                   1, mp, hz, ldh+1, ierr)
            if (ierr .eq. 0) then
               call dtrevc ('R', 'A', select, mp, hs, ldh+1, qdum, 1,
     &                      hv, ldh+1, mp, msel, work, ierr)
            end if
            if (ierr .ne. 0) then
               harmfb = .true.
               if (msglvl .gt. 0) then
                  call ivout (logfil, 1, [ierr], ndigit,
     &                 '_haup2: Schur form of HP failed, using Ritz')
               end if
               go to 31
            end if
c
c           %--------------------------------------------------%
c           | IZERO is the position of the zero eigenvalue of  |
c           | HP, where HZ^T*U = +-e_{IZERO} (\Remarks 2).     |
c           %--------------------------------------------------%
c
            call dgemv ('T', mp, mp, one, hz, ldh+1, hu, 1, zero,
     &                  work, 1)
            izero = idamax (mp, work, 1)
            if (hwi(izero) .ne. zero) then
               harmfb = .true.
               if (msglvl .gt. 0) then
                  call ivout (logfil, 1, [izero], ndigit,
     &                 '_haup2: zero eigenvalue of HP not found')
               end if
               go to 31
            end if
c
c           %--------------------------------------------------%
c           | Harmonic Ritz values in Schur order, skipping    |
c           | IZERO, with residual norms (\Remarks 2). WORK(j) |
c           | = HZ(m+1,:)*HV(:,j) is the last component of the |
c           | j-th eigenvector of HP, which is divided by the  |
c           | norm of the eigenvector, where a complex pair is |
c           | stored as (real, imag) in consecutive columns of |
c           | HV. Column j of HV is zero below row j (or j+1   |
c           | for the first column of a complex pair).         |
c           %--------------------------------------------------%
c
            call dgemv ('T', mp, mp, one, hv, ldh+1, hz(mp,1), ldh+1,
     &                  zero, work, 1)
c
            i = 0
            iconj = 0
            do 26 j = 1, mp
               if (j .ne. izero) then
                  i = i + 1
                  ritzr(i) = hwr(j)
                  ritzi(i) = hwi(j)
                  if (hwi(j) .eq. zero) then
                     rtmp = abs(work(j) * uinv) / dnrm2 (j, hv(1,j), 1)
                  else if (iconj .eq. 0) then
                     rtmp = dlapy2 (work(j), work(j+1)) * abs(uinv)
     &                      / dlapy2 (dnrm2 (j+1, hv(1,j), 1),
     &                                dnrm2 (j+1, hv(1,j+1), 1))
                     iconj = 1
                  else
                     iconj = 0
                  end if
                  bounds(i) = dlapy2 (hwr(j), hwi(j)) * rtmp
     &                        / dlapy2 (one, rtmp)
               end if
   26       continue
         else
c
c           %-------------------------------------------------%
c           | Real Schur form HS = HZ^T*H*HZ, and the         |
c           | eigenvectors of HS stored in Q.                 |
c           %-------------------------------------------------%
c
            call dlacpy ('All', kplusp, kplusp, h, ldh, hs, ldh+1)
            call dlaset ('All', kplusp, kplusp, zero, one, hz, ldh+1)
            call dlahqr (.true., .true., kplusp, 1, kplusp, hs, ldh+1,
     &                   ritzr, ritzi, 1, kplusp, hz, ldh+1, ierr)
            if (ierr .ne. 0) then
               info = -8
               go to 1200
            end if
c
            call dtrevc ('R', 'A', select, kplusp, hs, ldh+1, qdum, 1,
     &                   q, ldq, kplusp, msel, work, ierr)
            if (ierr .ne. 0) then
               info = -8
               go to 1200
            end if
c
c           %--------------------------------------------------%
c           | Ritz estimates RNORM*|e_m^T*HZ*y|, where y is a  |
c           | unit eigenvector of HS. WORKL(j) = HZ(m,:)*Q(:,j) |
c           | is divided by the norm of Q(:,j), where a complex |
c           | pair is stored as (real, imag) in consecutive    |
c           | columns of Q. Column j of Q is zero below row j  |
c           | (or j+1 for the first column of a complex pair). |
c           %--------------------------------------------------%
c
            call dgemv ('T', kplusp, kplusp, one, q, ldq, hz(kplusp,1),
     &                  ldh+1, zero, workl, 1)
            iconj = 0
            do 29 j = 1, kplusp
               if (ritzi(j) .eq. zero) then
                  bounds(j) = rnorm * abs(workl(j))
     &                        / dnrm2 (j, q(1,j), 1)
               else if (iconj .eq. 0) then
                  bounds(j) = rnorm * dlapy2 (workl(j), workl(j+1))
     &                        / dlapy2 (dnrm2 (j+1, q(1,j), 1),
     &                                  dnrm2 (j+1, q(1,j+1), 1))
                  bounds(j+1) = bounds(j)
                  iconj = 1
               else
                  iconj = 0
               end if
   29       continue
c
c           %------------------------------------------------%
c           | Eigenvalues in Schur order, used to select the |
c           | blocks kept by the restart.                    |
c           %------------------------------------------------%
c
            call dcopy (kplusp, ritzr, 1, hwr, 1)
            call dcopy (kplusp, ritzi, 1, hwi, 1)
         end if
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
c        | blocks, each represented by its eigenvalue. For harmonic |
c        | Ritz values, the Schur form is of order KPLUSP+1 and the |
c        | zero eigenvalue at IZERO is not a block, since it is     |
c        | always kept.                                             |
c        %----------------------------------------------------------%
c
         ntot = kplusp
         if (.not. harmfb) ntot = kplusp + 1
         nblk = 0
         i = 1
   55    continue
         if (i .le. ntot) then
            if (.not. harmfb .and. i .eq. izero) then
               i = i + 1
               go to 55
            end if
            nblk = nblk + 1
            blockstart(nblk) = i
            if (i .lt. ntot .and. abs(hs(i+1,i)) .gt. zero) then
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
         pick = 0
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
c        | boundary by increasing NEV by one. If that would leave   |
c        | NP = 0, drop the pair instead so that NP stays positive. |
c        %----------------------------------------------------------%
c
         if (msel .gt. nev) then
            if (msel .ge. kplusp) then
               consumed(pick) = 0
               msel = msel - blocklen(pick)
            end if
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
c           | Reorder the Schur form of HP so the NEV        |
c           | selected values and zero lead, and apply the   |
c           | harmonic thick restart. If DTRSEN cannot       |
c           | reorder, redo this iteration with Ritz values. |
c           %------------------------------------------------%
c
            select(izero) = .true.
            call dtrsen ('N', 'V', select, ntot, hs, ldh+1, hz, ldh+1,
     &                   hwr, hwi, msel, sep1, sep2, work, lwork,
     &                   iwork, 1, ierr)
            if (ierr .ne. 0 .or. msel .ne. nev+1) then
               harmfb = .true.
               if (msglvl .gt. 0) then
                  call ivout (logfil, 1, [ierr], ndigit,
     &                 '_haup2: DTRSEN failed on HP, using Ritz')
               end if
               go to 31
            end if
c
            call dhapps (n, nev, np, v, ldv, h, ldh, resid, rnorm,
     &                   hs, ldh+1, hz, ldh+1)
         else
c
c           %------------------------------------------------%
c           | Reorder the Schur form so the NEV selected     |
c           | values lead, and apply the bordered restart.   |
c           %------------------------------------------------%
c
            call dtrsen ('N', 'V', select, kplusp, hs, ldh+1, hz, ldh+1,
     &                   hwr, hwi, msel, sep1, sep2, work, lwork,
     &                   iwork, 1, ierr)
            if (ierr .ne. 0 .or. msel .ne. nev) then
               info = -8
               go to 1200
            end if
c
            call dqapps (n, nev, np, v, ldv, h, ldh, resid, rnorm,
     &                   q, ldq, hs, ldh+1, hz, ldh+1, workd, .true.)
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
