c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dqapps
c
c\Description:
c  Given the Arnoldi factorization
c
c     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
c
c  apply NP shifts implicitly (via the NONSYMMETRIC ARROWHEAD restart,
c  see \Remarks below) resulting in
c
c     A*(V_{k}*X) - (V_{k}*X)*(Hk) = rnew_{k}*e_{k}^T
c
c  where X = Qshr*Qk deflates the NP undesired eigenvalues of H (see
c  \References). The updated Arnoldi factorization becomes:
c
c     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
c
c  HOUSE selects how the (KEV+1) by (KEV+1) nonsymmetric arrowhead
c  matrix is reduced to upper Hessenberg form.
c
c\Usage:
c  call dqapps
c     ( N, KEV, NP, V, LDV, H, LDH, RESID, Q, LDQ,
c       TSCHUR, LDTSCHUR, QSCHUR, LDQSCHUR, WORKD, HOUSE )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Problem size, i.e. dimension of matrix A.
c
c  KEV     Integer.  (INPUT)
c          KEV+NP is the size of the input matrix H. KEV is the size
c          of the updated matrix HNEW. Unlike DNAPPS, this is a pure
c          INPUT here: the caller (dqaup2) is responsible for
c          incrementing KEV by one BEFORE calling this routine if
c          doing otherwise would split a complex-conjugate pair
c          across the KEV/NP boundary in TSCHUR (see \Remarks).
c
c  NP      Integer.  (INPUT)
c          Number of undesired eigenvalues being deflated.
c
c  V       Double precision N by (KEV+NP) array.  (INPUT/OUTPUT)
c          INPUT: V contains the current KEV+NP Arnoldi vectors.
c          OUTPUT: VNEW = V(1:n,1:KEV); the updated Arnoldi vectors
c          are in the first KEV columns of V.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  H       Double precision (KEV+NP) by (KEV+NP) array.  (INPUT/OUTPUT)
c          INPUT: not referenced (TSCHUR already carries everything
c          that is needed about the current KEV+NP upper Hessenberg
c          matrix -- kept as an argument only for interface parity
c          with DNAPPS).
c          OUTPUT: H contains the updated KEV by KEV upper Hessenberg
c          matrix in the leading KEV by KEV submatrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          INPUT: RESID contains the residual vector r_{k+p}.
c          OUTPUT: RESID is the updated residual vector rnew_{k}.
c
c  Q       Double precision KEV+NP by KEV+NP work array.  (WORKSPACE,
c          unused -- kept for interface parity with DNAPPS.)
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  TSCHUR  Double precision (KEV+NP) by (KEV+NP) array.  (INPUT)
c          The real Schur form of the current KEV+NP upper Hessenberg
c          matrix, S = QSCHUR^T * H * QSCHUR, computed externally by
c          the caller and REORDERED so that the KEV desired
c          eigenvalues occupy the leading KEV by KEV quasi-upper-
c          triangular block. Only TSCHUR(1:KEV,1:KEV) and
c          TSCHUR(KEV+NP,KEV+NP) are referenced (the latter only to
c          seed the arbitrary arrowhead hub value "delta"; see
c          \Remarks). If KEV falls in the middle of a 2x2 complex-
c          conjugate block, the caller must have already incremented
c          KEV by one so this never happens.
c
c  LDTSCHUR Integer.  (INPUT)
c          Leading dimension of TSCHUR exactly as declared in the
c          calling program.
c
c  QSCHUR  Double precision (KEV+NP) by (KEV+NP) array.  (INPUT)
c          The orthogonal Schur vectors corresponding to TSCHUR,
c          reordered to match (i.e. QSCHUR(:,1:KEV) spans the same
c          invariant subspace as TSCHUR(1:KEV,1:KEV)). Only row
c          KEV+NP and the first KEV columns are referenced.
c
c  LDQSCHUR Integer.  (INPUT)
c          Leading dimension of QSCHUR exactly as declared in the
c          calling program.
c
c  WORKD   Double precision work array of length 2*N.  (WORKSPACE,
c          unused -- kept for interface parity with DNAPPS.)
c
c  HOUSE   Logical.  (INPUT)
c          Selects which routine reduces the (KEV+1) by (KEV+1)
c          nonsymmetric arrowhead matrix to upper Hessenberg form:
c          .TRUE.  uses LAPACK's Householder reduction (DGEHRD/DORGHR).
c          .FALSE. uses one-way Givens chasing (dahhgv).
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
c     dahhgv  Local routine that reduces the nonsymmetric
c             arrowhead matrix directly to upper Hessenberg form via
c             a one-way Givens chasing scheme (no rotation of the
c             arrowhead is required -- see dahhgv.f). Used when
c             HOUSE=.FALSE.
c     dgehrd  LAPACK routine that reduces a general matrix to upper
c             Hessenberg form via Householder reflectors. Used when
c             HOUSE=.TRUE.
c     dorghr  LAPACK routine that generates the explicit orthogonal
c             matrix from DGEHRD's packed reflectors. Used when
c             HOUSE=.TRUE.
c     arscnd  ARPACK utility routine for timing.
c     dlaset  LAPACK matrix initialization routine.
c     dlacpy  LAPACK matrix copy routine.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     dgemm   Level 3 BLAS routine for matrix-matrix multiplication.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     dscal   Level 1 BLAS that scales a vector.
c     dnrm2   Level 1 BLAS that computes the Euclidean norm of a vector.
c
c\Remarks
c  1. This routine replaces DNAPPS's classical implicit-shift bulge-
c     chasing with the nonsymmetric arrowhead restart of Reference 2.
c     Given the caller's Schur decomposition of the current KEV+NP
c     upper Hessenberg matrix, reordered so the KEV desired
c     eigenvalues occupy TSCHUR's leading block, form the (KEV+1) by
c     (KEV+1) nonsymmetric arrowhead matrix
c
c         D = | TSCHUR(1:KEV,1:KEV)      c*QSCHUR(m,1:KEV)' |
c             | c*QSCHUR(m,1:KEV)              delta        |
c
c     where m = KEV+NP and c = RNORM = ||RESID||, reduce it to upper
c     Hessenberg form, and read the updated H, V, and RESID off that
c     result. "delta" is mathematically arbitrary (it never affects
c     the resulting HNEW, VNEW, RESID); TSCHUR(m,m) is used here
c     purely as a convenient, well-scaled placeholder.
c  2. When HOUSE=.FALSE., D is reduced directly (as built above) via
c     dahhgv. When HOUSE=.TRUE., LAPACK's general-purpose
c     Householder Hessenberg reduction (DGEHRD, exactly like MATLAB's
c     builtin HESS) does *not* preserve the required structure if
c     applied directly to D -- following the discussion of Reference
c     2, D must first be rotated via the permutation Pi and
c     transposed:
c
c         Pi*D^T*Pi = | delta            c*e_m^T*Qk   |
c                     | c*Qk^T*e_m       Pi*Sk^T*Pi   |
c
c     which moves the hub to position (1,1) and reverses the interior
c     block. DGEHRD/DORGHR are applied to THIS rotated matrix (called
c     DROT below); the result is then implicitly un-rotated by simply
c     reading its entries back out with the same index reversal.
c  3. Unlike DNAPPS, KEV is never adjusted inside this routine to
c     avoid splitting a complex-conjugate pair -- the caller must
c     arrange TSCHUR/QSCHUR (and increment KEV) BEFORE calling.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dqapps
     &   ( n, kev, np, v, ldv, h, ldh, resid, q, ldq,
     &     tschur, ldtschur, qschur, ldqschur, workd, house )
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
      logical    house
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
      logical    initd
      save       initd
      data       initd /.false./
      Double precision
     &           rnorm, betak, qwork(1)
      Double precision
     &           d(kev+1,kev+1), qarrow(kev+1,kev+1),
     &           drot(kev+1,kev+1), hbuf(kev+1,kev+1), tau(kev),
     &           qk1(kev,kev), qfinal(kev+np,kev)
      Double precision, allocatable :: housework(:), vnew(:,:)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, dscal, dlaset, dlacpy, dgemv, dgemm,
     &           dahhgv, dgehrd, dorghr, arscnd
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
c     %----------------------------------------------------------%
c     | debug.h/stat.h COMMON block variables aren't set up by   |
c     | any ARPACK driver program here (this routine is called   |
c     | directly from a MEX gateway) -- initialize the ones this |
c     | routine reads/writes (LOGFIL/NDIGIT/MNAPPS/TNAPPS) once, |
c     | the first time this routine is ever called.              |
c     %----------------------------------------------------------%
c
      if (.not. initd) then
         logfil = 6
         ndigit = -3
         mnapps = 0
         tnapps = 0.0D+0
         initd = .true.
      end if
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c
      call arscnd (t0)
      msglvl = mnapps
c
      kplusp = kev + np
c
c     %----------------------------------------------%
c     | Quick return if there are no shifts to apply |
c     %----------------------------------------------%
c
      if (np .eq. 0) go to 9000
c
      if (house) then
c
c        %--------------------------------------------------------------%
c        | Workspace query: ask DGEHRD and DORGHR for their optimal     |
c        | (blocked) LWORK instead of hardcoding the unblocked minimum, |
c        | so LAPACK can use its blocked kernels for larger KEV. The    |
c        | query calls (LWORK = -1) only inspect dimensions -- no array |
c        | contents are referenced. 8*(KEV+1) is kept as a floor.       |
c        %--------------------------------------------------------------%
c
         call dgehrd (kev+1, 1, kev+1, drot, kev+1, tau, qwork, -1,
     &                houseinfo)
         lwork = int(qwork(1))
         call dorghr (kev+1, 1, kev+1, drot, kev+1, tau, qwork, -1,
     &                houseinfo)
         lwork = max(lwork, int(qwork(1)), 8*(kev+1))
         allocate (housework(lwork))
      end if
c
c     %------------------------------------------------------%
c     | RNORM: the current residual norm, needed both as the |
c     | arrowhead spike scale and to normalize the residual  |
c     | direction when forming the updated residual below.   |
c     %------------------------------------------------------%
c
      rnorm = dnrm2(n, resid, 1)
c
      if (house) then
c
c        %---------------------------------------------------------------%
c        | Build the ALREADY 180-degree-rotated (and transposed) (KEV+1) |
c        | by (KEV+1) arrowhead matrix DROT (hub at (1,1), border along  |
c        | row/column 1, interior block = Pi*Sk^T*Pi) -- mathematically  |
c        | identical to building the natural downward-pointing arrowhead |
c        | and then applying the rotation of \Remarks (2), but without   |
c        | ever forming that un-rotated matrix.                          |
c        %---------------------------------------------------------------%
c
         call dlaset ('All', kev+1, kev+1, zero, zero, drot, kev+1)
c
         do 50 j = 2, kev+1
            drot(1,j) = rnorm * qschur(kplusp, kev+2-j)
            drot(j,1) = drot(1,j)
   50    continue
c
         do 70 j = 2, kev+1
            do 60 i = 2, kev+1
               drot(i,j) = tschur(kev+2-j, kev+2-i)
   60       continue
   70    continue
c
c        %--------------------------------------------------------------%
c        | Reduce DROT to upper Hessenberg form via LAPACK's general    |
c        | Householder reduction. DGEHRD packs the Hessenberg values    |
c        | into the upper triangle + first subdiagonal of DROT, with    |
c        | Householder vector data below that -- copy the valid         |
c        | Hessenberg entries out to HBUF before DORGHR overwrites DROT |
c        | with the explicit orthogonal matrix.                         |
c        %--------------------------------------------------------------%
c
         call dgehrd (kev+1, 1, kev+1, drot, kev+1, tau,
     &                housework, lwork, houseinfo)
c
         call dlaset ('All', kev+1, kev+1, zero, zero, hbuf, kev+1)
         do 90 j = 1, kev+1
            do 80 i = 1, min(j+1,kev+1)
               hbuf(i,j) = drot(i,j)
   80       continue
   90    continue
c
         call dorghr (kev+1, 1, kev+1, drot, kev+1, tau,
     &                housework, lwork, houseinfo)
c
         deallocate (housework)
c
c        %--------------------------------------------------------------%
c        | Read the new leading KEV by KEV upper Hessenberg H, and      |
c        | BETAK (the new residual coupling), directly out of HBUF by   |
c        | un-rotating (the SAME index reversal used to build DROT,     |
c        | since Pi is an involution) -- entries strictly below H's own |
c        | subdiagonal are explicitly zeroed rather than read out of    |
c        | HBUF, since the corresponding HBUF locations there hold      |
c        | leftover Householder-vector data, not Hessenberg zeros.      |
c        %--------------------------------------------------------------%
c
         call dlaset ('All', kev, kev, zero, zero, h, ldh)
         do 110 j = 1, kev
            do 100 i = 1, min(j+1,kev)
               h(i,j) = hbuf(kev+2-j, kev+2-i)
  100       continue
  110    continue
c
         betak = hbuf(2,1)
c
c        %--------------------------------------------------------------%
c        | Build QFINAL = QSCHUR(:,1:kev) * QK1(1:kev,1:kev), the       |
c        | KPLUSP by KEV transformation carrying V's current KPLUSP-    |
c        | column basis directly onto the new KEV-column basis. QK1 is  |
c        | read out of the now-orthogonal DROT the same way H was read  |
c        | out of HBUF above (both need the same index reversal).       |
c        %--------------------------------------------------------------%
c
         do 130 j = 1, kev
            do 120 p = 1, kev
               qk1(p,j) = drot(kev+2-p,kev+2-j)
  120       continue
  130    continue
c
         call dgemm ('N', 'N', kplusp, kev, kev, one, qschur, ldqschur,
     &               qk1, kev, zero, qfinal, kplusp)
c
      else
c
c        %--------------------------------------------------------%
c        | Build the natural (downward-pointing) nonsymmetric     |
c        | arrowhead matrix D directly -- no rotation is needed   |
c        | since dahhgv operates on this form as-is.        |
c        %--------------------------------------------------------%
c
         call dlaset ('All', kev+1, kev+1, zero, zero, d, kev+1)
c
         call dlacpy ('All', kev, kev, tschur, ldtschur, d, kev+1)
c
         do 700 i = 1, kev
            d(i,kev+1) = rnorm * qschur(kplusp,i)
            d(kev+1,i) = rnorm * qschur(kplusp,i)
  700    continue
c
         d(kev+1,kev+1) = tschur(kplusp,kplusp)
c
c        %--------------------------------------------------------%
c        | Reduce D directly to upper Hessenberg form.            |
c        %--------------------------------------------------------%
c
         call dahhgv (kev+1, d, kev+1, qarrow, kev+1)
c
c        %--------------------------------------------------------%
c        | Read the new leading KEV by KEV upper Hessenberg H,    |
c        | and BETAK (the new residual coupling), directly out of |
c        | the reduced result.                                    |
c        %--------------------------------------------------------%
c
         call dlacpy ('All', kev, kev, d, kev+1, h, ldh)
c
         betak = d(kev+1,kev)
c
c        %--------------------------------------------------------%
c        | Build QFINAL = QSCHUR(:,1:kev) * QARROW(1:kev,1:kev),  |
c        | the KPLUSP by KEV transformation carrying V's current  |
c        | KPLUSP-column basis directly onto the new KEV-column   |
c        | basis. Unlike the House variant (and the symmetric     |
c        | routines), QARROW is already in the natural, unreversed|
c        | order, so this is a plain matrix product               |
c        %--------------------------------------------------------%
c
         call dgemm ('N', 'N', kplusp, kev, kev, one, qschur, ldqschur,
     &               qarrow, kev+1, zero, qfinal, kplusp)
c
      end if
c
c     %--------------------------------------------------------%
c     | Update V: V(:,1:kev) <- V(:,1:kplusp) * QFINAL.        |
c     %--------------------------------------------------------%
c
      allocate (vnew(n,kev))
      call dgemm ('N', 'N', n, kev, kplusp, one, v, ldv, qfinal,
     &            kplusp, zero, vnew, n)
      call dlacpy ('All', n, kev, vnew, n, v, ldv)
      deallocate (vnew)
c
c     %--------------------------------------------------------%
c     | Update the residual vector. As shown in Reference 2,   |
c     | RESID_new is exactly BETAK times the ORIGINAL residual |
c     | direction -- i.e. simply rescaled, never re-oriented:  |
c     |    resid <- (betak / rnorm) * resid                    |
c     %--------------------------------------------------------%
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
c     %-----------------%
c     | End of dqapps |
c     %-----------------%
c
      end
