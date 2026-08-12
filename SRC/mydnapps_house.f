c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: mydnapps_house
c
c\Description:
c  House (Householder / DGEHRD-DORGHR) variant of mydnapps: given the
c  Arnoldi factorization
c
c     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
c
c  apply NP shifts implicitly via the NONSYMMETRIC ARROWHEAD restart
c  (see \Remarks), reducing the arrowhead matrix to upper Hessenberg
c  form with LAPACK's general Householder Hessenberg reduction
c  (DGEHRD/DORGHR) instead of the Givens chasing scheme used by
c  mydnapps/narrowgivens. The updated Arnoldi factorization becomes:
c
c     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
c
c\Usage:
c  call mydnapps_house
c     ( N, KEV, NP, V, LDV, H, LDH, RESID, Q, LDQ,
c       TSCHUR, LDTSCHUR, QSCHUR, LDQSCHUR, WORKD )
c
c\Arguments
c  Identical to mydnapps -- see mydnapps.f for the full description of
c  every argument.
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
c     dgehrd  LAPACK routine that reduces a general matrix to upper
c             Hessenberg form via Householder reflectors.
c     dorghr  LAPACK routine that generates the explicit orthogonal
c             matrix from DGEHRD's packed reflectors.
c     dlaset  LAPACK matrix initialization routine.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     dgemm   Level 3 BLAS routine for matrix-matrix multiplication.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     dscal   Level 1 BLAS that scales a vector.
c     dnrm2   Level 1 BLAS that computes the Euclidean norm of a vector.
c     arscnd  ARPACK utility routine for timing.
c
c\Remarks
c  1. UNLIKE narrowgivens (which chases the spike directly out of the
c     natural "downward-pointing" arrowhead D, no rotation needed),
c     LAPACK's general-purpose Householder Hessenberg reduction
c     (DGEHRD, exactly like MATLAB's builtin HESS) does *not* preserve
c     the required structure if applied directly to D. Following the
c     discussion of Reference 2, D must first be rotated via the
c     permutation Pi and transposed:
c
c         Pi*D^T*Pi = | delta            c*e_m^T*Qk   |
c                     | c*Qk^T*e_m       Pi*Sk^T*Pi   |
c
c     which moves the hub to position (1,1) and reverses the interior
c     block. DGEHRD/DORGHR are applied to THIS rotated matrix (called
c     DROT below); the result is then implicitly un-rotated by simply
c     reading its entries back out with the same index reversal.
c  2. Just as in mydnapps, KEV is never adjusted inside this routine
c     to avoid splitting a complex-conjugate pair; the caller must
c     handle that before calling.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine mydnapps_house
     &   ( n, kev, np, v, ldv, h, ldh, resid, q, ldq,
     &     tschur, ldtschur, qschur, ldqschur, workd )
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
      integer    kev, ldh, ldq, ldqschur, ldtschur, ldv, n, np
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           h(ldh,kev+np), q(ldq,kev+np), qschur(ldqschur,kev+np),
     &           resid(n), tschur(ldtschur,kev+np), v(ldv,kev+np),
     &           workd(2*n)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i, j, p, kplusp, msglvl, houseinfo, lwork
      Double precision
     &           rnorm, betak, qwork(1)
      Double precision
     &           drot(kev+1,kev+1), hbuf(kev+1,kev+1),
     &           tau(kev),
     &           qfinal(kev+np,kev), vnew(n,kev), qk1(kev,kev)
      Double precision, allocatable :: housework(:)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, dscal, dlaset, dgemv, dgemm, dlacpy,
     &           dgehrd, dorghr, arscnd
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dnrm2
      external   dnrm2
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      call arscnd (t0)
      msglvl = mnapps
c
      kplusp = kev + np
c
      if (np .eq. 0) go to 9000
c
c     %--------------------------------------------------------------%
c     | Workspace query: ask DGEHRD and DORGHR for their optimal     |
c     | (blocked) LWORK instead of hardcoding the unblocked minimum, |
c     | so LAPACK can use its blocked kernels for larger KEV. The    |
c     | query calls (LWORK = -1) only inspect dimensions -- no array |
c     | contents are referenced. 8*(KEV+1) is kept as a floor.       |
c     %--------------------------------------------------------------%
c
      call dgehrd (kev+1, 1, kev+1, drot, kev+1, tau, qwork, -1,
     &             houseinfo)
      lwork = int(qwork(1))
      call dorghr (kev+1, 1, kev+1, drot, kev+1, tau, qwork, -1,
     &             houseinfo)
      lwork = max(lwork, int(qwork(1)), 8*(kev+1))
      allocate (housework(lwork))
c
      rnorm = dnrm2(n, resid, 1)
c
c     %---------------------------------------------------------------%
c     | Build the ALREADY 180-degree-rotated (and transposed) (KEV+1) |
c     | by (KEV+1) arrowhead matrix DROT (hub at (1,1), border along  |
c     | row/column 1, interior block = Pi*Sk^T*Pi) -- mathematically  |
c     | identical to building the natural downward-pointing arrowhead |
c     | and then applying the rotation of \Remarks (1), but without   |
c     | ever forming that un-rotated matrix.                          |
c     %---------------------------------------------------------------%
c
      call dlaset ('All', kev+1, kev+1, zero, zero, drot, kev+1)
c
      do 50 j = 2, kev+1
         drot(1,j) = rnorm * qschur(kplusp, kev+2-j)
         drot(j,1) = drot(1,j)
   50 continue
c
      do 70 j = 2, kev+1
         do 60 i = 2, kev+1
            drot(i,j) = tschur(kev+2-j, kev+2-i)
   60    continue
   70 continue
c
c     %--------------------------------------------------------------%
c     | Reduce DROT to upper Hessenberg form via LAPACK's general    |
c     | Householder reduction. DGEHRD packs the Hessenberg values    |
c     | into the upper triangle + first subdiagonal of DROT, with    |
c     | Householder vector data below that -- copy the valid         |
c     | Hessenberg entries out to HBUF before DORGHR overwrites DROT |
c     | with the explicit orthogonal matrix.                         |
c     %--------------------------------------------------------------%
c
      call dgehrd (kev+1, 1, kev+1, drot, kev+1, tau,
     &             housework, lwork, houseinfo)
c
      call dlaset ('All', kev+1, kev+1, zero, zero, hbuf, kev+1)
      do 90 j = 1, kev+1
         do 80 i = 1, min(j+1,kev+1)
            hbuf(i,j) = drot(i,j)
   80    continue
   90 continue
c
      call dorghr (kev+1, 1, kev+1, drot, kev+1, tau,
     &             housework, lwork, houseinfo)
c
      deallocate (housework)
c
c     %--------------------------------------------------------------%
c     | Read the new leading KEV by KEV upper Hessenberg H, and      |
c     | BETAK (the new residual coupling), directly out of HBUF by   |
c     | un-rotating (the SAME index reversal used to build DROT,     |
c     | since Pi is an involution) -- entries strictly below H's own |
c     | subdiagonal are explicitly zeroed rather than read out of    |
c     | HBUF, since the corresponding HBUF locations there hold      |
c     | leftover Householder-vector data, not Hessenberg zeros.      |
c     %--------------------------------------------------------------%
c
      call dlaset ('All', kev, kev, zero, zero, h, ldh)
      do 110 j = 1, kev
         do 100 i = 1, min(j+1,kev)
            h(i,j) = hbuf(kev+2-j, kev+2-i)
  100    continue
  110 continue
c
      betak = hbuf(2,1)
c
c     %--------------------------------------------------------------%
c     | Build QFINAL = QSCHUR(:,1:kev) * QK1(1:kev,1:kev), the       |
c     | KPLUSP by KEV transformation carrying V's current KPLUSP-    |
c     | column basis directly onto the new KEV-column basis. QK1 is  |
c     | read out of the now-orthogonal DROT the same way H was read  |
c     | out of HBUF above (both need the same index reversal)--once  |
c     | QK1 is materialized explicitly (it's only KEV by KEV), the   |
c     | actual KPLUSP by KEV product is handed off to DGEM           |
c     %--------------------------------------------------------------%
c
      do 130 j = 1, kev
         do 120 p = 1, kev
            qk1(p,j) = drot(kev+2-p,kev+2-j)
  120    continue
  130 continue
c
      call dgemm ('N', 'N', kplusp, kev, kev, one, qschur, ldqschur,
     &            qk1, kev, zero, qfinal, kplusp)
c
      call dgemm ('N', 'N', n, kev, kplusp, one, v, ldv, qfinal,
     &            kplusp, zero, vnew, n)
      call dlacpy ('All', n, kev, vnew, n, v, ldv)
c
      if (rnorm .gt. zero) then
         call dscal (n, betak/rnorm, resid, 1)
      end if
c
      call arscnd (t1)
      tnapps = tnapps + (t1 - t0)
c
 9000 continue
      return
c
c     %-----------------------%
c     | End of mydnapps_house |
c     %-----------------------%
c
      end
