c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dhapps
c
c\Description:
c  Given the Arnoldi factorization
c
c     A*V_{k+p} - V_{k+p}*H_{k+p} = r_{k+p}*e_{k+p}^T,
c
c  apply the nonsymmetric HARMONIC thick restart (see \Remarks below),
c  keeping the KEV harmonic Ritz pairs selected by the caller, resulting
c  in the updated Arnoldi factorization
c
c     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
c
c\Usage:
c  call dhapps
c     ( N, KEV, NP, V, LDV, H, LDH, RESID, TSCHUR, LDTSCHUR,
c       QSCHUR, LDQSCHUR, HINVEM, WORKD )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Problem size, i.e. dimension of matrix A.
c
c  KEV     Integer.  (INPUT)
c          KEV+NP is the size of the input matrix H. KEV is the size
c          of the updated matrix HNEW. The caller is responsible for
c          incrementing KEV by one before calling this routine if doing
c          otherwise would split a complex-conjugate pair across the
c          KEV/NP boundary in TSCHUR (see \Remarks).
c
c  NP      Integer.  (INPUT)
c          Number of undesired harmonic Ritz values being dropped.
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
c          INPUT: the current KEV+NP upper Hessenberg matrix.
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
c  TSCHUR  Double precision (KEV+NP) by (KEV+NP) array.  (INPUT)
c          The real Schur form of the current KEV+NP harmonic matrix
c          G = H + RNORM**2 * (H'\E_{KEV+NP}) * E_{KEV+NP}^T, computed
c          externally by the caller and REORDERED so the KEV desired
c          harmonic Ritz values occupy the leading KEV by KEV quasi-
c          upper-triangular block. If KEV falls in the middle of a 2x2
c          complex-conjugate block, the caller must have already
c          incremented KEV by one so this never happens.
c
c  LDTSCHUR Integer.  (INPUT)
c          Leading dimension of TSCHUR exactly as declared in the
c          calling program.
c
c  QSCHUR  Double precision (KEV+NP) by (KEV+NP) array.  (INPUT)
c          The orthogonal Schur vectors corresponding to TSCHUR,
c          reordered to match (i.e. QSCHUR(:,1:KEV) spans the same
c          invariant subspace, under G, as TSCHUR(1:KEV,1:KEV)).
c
c  LDQSCHUR Integer.  (INPUT)
c          Leading dimension of QSCHUR exactly as declared in the
c          calling program.
c
c  HINVEM  Double precision array of length KEV+NP.  (INPUT)
c          H(1:KEV+NP,1:KEV+NP)' \ E_{KEV+NP}, i.e. the solution of
c          H^T * x = e_{KEV+NP}, precomputed by the caller (it is
c          needed there anyway, to form G before its Schur
c          decomposition).
c
c  WORKD   Double precision work array of length N.  (WORKSPACE)
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\References:
c  1. A Unified View of Arrowhead Matrix Transformations and Lanczos
c     Restarts (and nonsymmetric companion), James Baglama,
c     Kyle Monette, Vasilije Perovic, (2026).
c
c\Routines called:
c     dgehrd  LAPACK routine that reduces a general matrix to upper
c             Hessenberg form via Householder reflectors.
c     dorghr  LAPACK routine that generates the explicit orthogonal
c             matrix from DGEHRD's packed reflectors.
c     arscnd  ARPACK utility routine for timing.
c     dlaset  LAPACK matrix initialization routine.
c     dlacpy  LAPACK matrix copy routine.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     dger    Level 2 BLAS rank-one update.
c     dgemm   Level 3 BLAS routine for matrix-matrix multiplication.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     daxpy   Level 1 BLAS scalar times a vector plus a vector.
c     dscal   Level 1 BLAS that scales a vector.
c     ddot    Level 1 BLAS dot product.
c     dnrm2   Level 1 BLAS that computes the Euclidean norm of a vector.
c
c\Remarks
c  1. This routine replaces DNAPPS's classical implicit-shift bulge-
c     chasing with the nonsymmetric HARMONIC thick restart of
c     Reference 1. Given the caller's Schur decomposition of the
c     harmonic matrix G (built from the current KEV+NP upper
c     Hessenberg H), reordered so the KEV desired harmonic Ritz values
c     occupy TSCHUR's leading block, this routine:
c       a) builds the (KEV+NP+1) by (KEV+1) extended basis PEXT whose
c          first KEV columns are QSCHUR(:,1:KEV) padded with a zero
c          row, and whose last column is the arrow vector
c              Q = [ -RNORM**2 * HINVEM ; RNORM ]
c          Gram-Schmidt orthogonalized against those KEV columns and
c          normalized (RNORM = ||RESID||);
c       b) forms the (KEV+1) by KEV matrix
c              HNEW = PEXT' * H0 * QSCHUR(1:KEV+NP,1:KEV)
c          where H0 is the (KEV+NP+1) by (KEV+NP) matrix
c          [H; RNORM*E_{KEV+NP}^T] -- computed directly via BLAS calls
c          rather than a closed form, to mirror the reference MATLAB
c          implementation as closely as possible;
c       c) reduces HNEW (padded to a square (KEV+1) matrix with a
c          trailing zero column) to upper Hessenberg form via the SAME
c          180-degree-rotate-and-transpose trick DQAPPS uses for its
c          arrowhead matrix (see DQAPPS \Remarks 2) -- that trick only
c          depends on DGEHRD/DORGHR being applied to the natural
c          orientation of a matrix whose last row carries the residual
c          coupling, not on the matrix actually being arrowhead-shaped,
c          so it applies unchanged here.
c  2. Unlike DNAPPS, KEV is never adjusted inside this routine to
c     avoid splitting a complex-conjugate pair -- the caller must
c     arrange TSCHUR/QSCHUR (and increment KEV) BEFORE calling.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dhapps
     &   ( n, kev, np, v, ldv, h, ldh, resid,
     &     tschur, ldtschur, qschur, ldqschur, hinvem, workd )
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
      integer    kev, ldh, ldqschur, ldtschur, ldv, n, np
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           h(ldh,kev+np), hinvem(kev+np), qschur(ldqschur,kev+np),
     &           resid(n), tschur(ldtschur,kev+np), v(ldv,kev+np),
     &           workd(n)
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
      integer    i, j, kplusp, msglvl, houseinfo, lwork
      logical    initd
      save       initd
      data       initd /.false./
      Double precision
     &           rnorm, betak, znorm, qwork(1)
c
c     %-----------------------%
c     | Local Array Arguments |
c     %-----------------------%
c
      Double precision
     &           qorig(kev+np+1), mcoef(kev),
     &           temp1(kev+np,kev), hnew(kev+1,kev), hkpad(kev+1,kev+1),
     &           drot(kev+1,kev+1), hbuf(kev+1,kev+1), tau(kev),
     &           qk1(kev,kev), qfinal(kev+np,kev)
      Double precision, allocatable, save :: pext(:,:), housework(:),
     &           vnew(:,:)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, dscal, daxpy, dlaset, dlacpy, dgemv, dger,
     &           dgemm, dgehrd, dorghr, arscnd
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           ddot, dnrm2
      external   ddot, dnrm2
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %----------------------------------------------------------%
c     | debug.h/stat.h COMMON block variables aren't set up by   |
c     | any ARPACK driver program here. Initialize the ones this |
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
      call arscnd (t0)
      msglvl = mnapps
c
      kplusp = kev + np
c
      if (np .eq. 0) go to 9000
c
      if (allocated(pext)) then
         if (size(pext,1) .ne. kplusp+1) deallocate (pext)
      end if
      if (.not. allocated(pext)) allocate (pext(kplusp+1,kev+1))
c
      rnorm = dnrm2 (n, resid, 1)
c
c     %--------------------------------------------------------%
c     | Build the (KEV+NP+1) by (KEV+1) extended basis PEXT:   |
c     | first KEV columns are QSCHUR(:,1:KEV) padded with a    |
c     | zero row, last column is the Gram-Schmidt orthogonal-  |
c     | ized and normalized arrow vector                       |
c     |     QORIG = [ -RNORM**2 * HINVEM ; RNORM ]             |
c     %--------------------------------------------------------%
c
      call dlaset ('All', kplusp+1, kev+1, zero, zero, pext, kplusp+1)
      call dlacpy ('All', kplusp, kev, qschur, ldqschur, pext, kplusp+1)
c
      do 10 i = 1, kplusp
         qorig(i) = -rnorm**2 * hinvem(i)
   10 continue
      qorig(kplusp+1) = rnorm
c
      do 20 i = 1, kev
         mcoef(i) = ddot (kplusp+1, qorig, 1, pext(1,i), 1)
         call daxpy (kplusp+1, -mcoef(i), pext(1,i), 1, qorig, 1)
   20 continue
      znorm = dnrm2 (kplusp+1, qorig, 1)
      call dcopy (kplusp+1, qorig, 1, pext(1,kev+1), 1)
      call dscal (kplusp+1, one/znorm, pext(1,kev+1), 1)
c
c     %--------------------------------------------------------%
c     | HNEW = PEXT' * H0 * QSCHUR(1:kplusp,1:kev), where H0 is |
c     | [H; RNORM*E_{kplusp}^T]. Computed directly via BLAS,    |
c     | rather than a closed form, to mirror the reference      |
c     | MATLAB implementation as closely as possible. TEMP1 =   |
c     | H * QSCHUR(1:kplusp,1:kev) first, then HNEW = PEXT(1:   |
c     | kplusp,:)' * TEMP1, plus the rank-one contribution from |
c     | H0's last row (RNORM*E_{kplusp}^T), which only affects  |
c     | row KEV+1 of HNEW since PEXT(kplusp+1,1:kev) = 0.       |
c     %--------------------------------------------------------%
c
      call dgemm ('N', 'N', kplusp, kev, kplusp, one, h, ldh,
     &            qschur, ldqschur, zero, temp1, kplusp)
      call dgemm ('T', 'N', kev+1, kev, kplusp, one, pext, kplusp+1,
     &            temp1, kplusp, zero, hnew, kev+1)
      call daxpy (kev, pext(kplusp+1,kev+1)*rnorm, qschur(kplusp,1),
     &            ldqschur, hnew(kev+1,1), kev+1)
c
c     %--------------------------------------------------------%
c     | Pad HNEW with a trailing zero column, and build the    |
c     | ALREADY 180-degree-rotated (and transposed) (KEV+1) by |
c     | (KEV+1) matrix DROT -- see DQAPPS \Remarks 2 for why    |
c     | DGEHRD cannot be applied to HNEW directly. The identity |
c     | DROT(i,j) = HKPAD(kev+2-j,kev+2-i) is the same 180-deg  |
c     | rotation plus transpose used there, and is agnostic to  |
c     | what HKPAD actually represents.                         |
c     %--------------------------------------------------------%
c
      call dlaset ('All', kev+1, kev+1, zero, zero, hkpad, kev+1)
      call dlacpy ('All', kev+1, kev, hnew, kev+1, hkpad, kev+1)
c
      do 40 j = 1, kev+1
         do 30 i = 1, kev+1
            drot(i,j) = hkpad(kev+2-j, kev+2-i)
   30    continue
   40 continue
c
c     %--------------------------------------------------------------%
c     | Workspace query: ask DGEHRD and DORGHR for their optimal     |
c     | (blocked) LWORK instead of hardcoding the unblocked minimum. |
c     %--------------------------------------------------------------%
c
      call dgehrd (kev+1, 1, kev+1, drot, kev+1, tau, qwork, -1,
     &             houseinfo)
      lwork = int(qwork(1))
      call dorghr (kev+1, 1, kev+1, drot, kev+1, tau, qwork, -1,
     &             houseinfo)
      lwork = max(lwork, int(qwork(1)), 8*(kev+1))
      if (allocated(housework)) then
         if (size(housework) .lt. lwork) deallocate (housework)
      end if
      if (.not. allocated(housework)) allocate (housework(lwork))
c
      call dgehrd (kev+1, 1, kev+1, drot, kev+1, tau,
     &             housework, lwork, houseinfo)
c
      call dlaset ('All', kev+1, kev+1, zero, zero, hbuf, kev+1)
      do 60 j = 1, kev+1
         do 50 i = 1, min(j+1,kev+1)
            hbuf(i,j) = drot(i,j)
   50    continue
   60 continue
c
      call dorghr (kev+1, 1, kev+1, drot, kev+1, tau,
     &             housework, lwork, houseinfo)
c
c     %--------------------------------------------------------------%
c     | Read the new leading KEV by KEV upper Hessenberg H, and      |
c     | BETAK (the new residual coupling), directly out of HBUF by   |
c     | un-rotating -- same index reversal as DQAPPS.                |
c     %--------------------------------------------------------------%
c
      call dlaset ('All', kev, kev, zero, zero, h, ldh)
      do 80 j = 1, kev
         do 70 i = 1, min(j+1,kev)
            h(i,j) = hbuf(kev+2-j, kev+2-i)
   70    continue
   80 continue
c
      betak = hbuf(2,1)
c
c     %--------------------------------------------------------------%
c     | Build QFINAL = QSCHUR(:,1:kev) * QK1(1:kev,1:kev), the       |
c     | KPLUSP by KEV transformation carrying V's current KPLUSP-    |
c     | column basis directly onto the new KEV-column basis.         |
c     %--------------------------------------------------------------%
c
      do 100 j = 1, kev
         call dcopy (kev, drot(2,kev+2-j), -1, qk1(1,j), 1)
  100 continue
c
      call dgemm ('N', 'N', kplusp, kev, kev, one, qschur, ldqschur,
     &            qk1, kev, zero, qfinal, kplusp)
c
c     %--------------------------------------------------------%
c     | Compute the new residual direction BEFORE overwriting   |
c     | V. The extended basis has KPLUSP+1 "rows" (the current  |
c     | KPLUSP Arnoldi vectors, plus the not-yet-built (KPLUSP+ |
c     | 1)-th direction RESID/RNORM), so this is               |
c     |     WORKD = V(1:n,1:kplusp)*PEXT(1:kplusp,kev+1)        |
c     |           + PEXT(kplusp+1,kev+1) * (RESID/RNORM)        |
c     | -- using the ORIGINAL RESID (still un-touched at this   |
c     | point), which is exactly RNORM times that (kplusp+1)-th |
c     | direction. The result is (up to roundoff) already a     |
c     | unit vector orthogonal to the new V(1:n,1:kev), since    |
c     | PEXT's columns are orthonormal and V's current columns  |
c     | (extended by RESID/RNORM) are orthonormal. One          |
c     | defensive re-orthogonalization pass against the new V   |
c     | is still applied, matching the reference MATLAB         |
c     | implementation's own such pass.                         |
c     %--------------------------------------------------------%
c
      call dgemv ('N', n, kplusp, one, v, ldv, pext(1,kev+1), 1,
     &            zero, workd, 1)
      call daxpy (n, pext(kplusp+1,kev+1)/rnorm, resid, 1, workd, 1)
c
c     %--------------------------------------------------------%
c     | Update V: V(:,1:kev) <- V(:,1:kplusp) * QFINAL.        |
c     %--------------------------------------------------------%
c
      if (allocated(vnew)) then
         if (size(vnew,1) .ne. n .or. size(vnew,2) .lt. kev)
     &      deallocate (vnew)
      end if
      if (.not. allocated(vnew)) allocate (vnew(n,kev))
      call dgemm ('N', 'N', n, kev, kplusp, one, v, ldv, qfinal,
     &            kplusp, zero, vnew, n)
      call dlacpy ('All', n, kev, vnew, n, v, ldv)
c
      do 110 i = 1, kev
         znorm = ddot (n, workd, 1, v(1,i), 1)
         call daxpy (n, -znorm, v(1,i), 1, workd, 1)
  110 continue
      znorm = dnrm2 (n, workd, 1)
      call dcopy (n, workd, 1, resid, 1)
      if (znorm .gt. zero) call dscal (n, betak/znorm, resid, 1)
c
      call arscnd (t1)
      tnapps = tnapps + (t1 - t0)
c
 9000 continue
      return
c
c     %-----------------%
c     | End of dhapps   |
c     %-----------------%
c
      end
