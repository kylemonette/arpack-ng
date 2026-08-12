c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: mydsapps
c
c\Description:
c  Given the Arnoldi factorization
c
c     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
c
c  apply NP shifts implicitly resulting in
c
c     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
c
c  where Q is an orthogonal matrix of order KEV+NP. Q is the product of
c  rotations resulting from the NP bulge chasing sweeps.  The updated Arnoldi
c  factorization becomes:
c
c     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
c
c\Usage:
c  call mydsapps
c     ( N, KEV, NP, V, LDV, H, LDH, RESID, Q, LDQ,
c       EIGVAL, EIGVEC, LDEIGVEC, WORKD )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Problem size, i.e. dimension of matrix A.
c
c  KEV     Integer.  (INPUT)
c          INPUT: KEV+NP is the size of the input matrix H.
c          OUTPUT: KEV is the size of the updated matrix HNEW.
c
c  NP      Integer.  (INPUT)
c          Number of implicit shifts to be applied.
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
c  H       Double precision (KEV+NP) by 2 array.  (INPUT/OUTPUT)
c          INPUT: H contains the symmetric tridiagonal matrix of the
c          Arnoldi factorization with the subdiagonal in the 1st column
c          starting at H(2,1) and the main diagonal in the 2nd column.
c          OUTPUT: H contains the updated tridiagonal matrix in the
c          KEV leading submatrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RESID   Double precision array of length (N).  (INPUT/OUTPUT)
c          INPUT: RESID contains the the residual vector r_{k+p}.
c          OUTPUT: RESID is the updated residual vector rnew_{k}.
c
c  Q       Double precision KEV+NP by KEV+NP work array.  (WORKSPACE)
c          Work array used to accumulate the rotations during the bulge
c          chase sweep.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  EIGVAL  Double precision array of length KEV+NP.  (INPUT)
c          The KEV+NP eigenvalues of the input tridiagonal matrix H,
c          computed externally by the caller. MUST be arranged so
c          that the first KEV entries, EIGVAL(1:KEV), are the KEV
c          "desired" (to be kept) eigenvalues -- the remaining NP
c          entries (if present) are never referenced and may hold
c          the undesired eigenvalues in any order, or anything else.
c
c  EIGVEC  Double precision KEV+NP by KEV+NP array.  (INPUT)
c          The corresponding KEV+NP eigenvectors of H (as columns),
c          computed externally by the caller, arranged to match
c          EIGVAL: EIGVEC(:,1:KEV) must be the KEV desired
c          eigenvectors, in the same order as EIGVAL(1:KEV). Columns
c          KEV+1:KEV+NP (if present) are never referenced.
c
c  LDEIGVEC Integer.  (INPUT)
c          Leading dimension of EIGVEC exactly as declared in the
c          calling program.
c
c  WORKD   Double precision work array of length 2*N.  (WORKSPACE)
c          Distributed array used in the application of the accumulated
c          orthogonal matrix Q.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c
c\Routines called:
c     ivout   ARPACK utility routine that prints integers.
c     arscnd  ARPACK utility routine for timing.
c     dvout   ARPACK utility routine that prints vectors.
c     dlamch  LAPACK routine that determines machine constants.
c     dlartg  LAPACK Givens rotation construction routine.
c     dlacpy  LAPACK matrix copy routine.
c     dlaset  LAPACK matrix initialization routine.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     daxpy   Level 1 BLAS that computes a vector triad.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     dscal   Level 1 BLAS that scales a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Revision history:
c     12/16/93: Version ' 2.4'
c
c\SCCS Information: @(#)
c FILE: sapps.F   SID: 2.6   DATE OF SID: 3/28/97   RELEASE: 2
c
c\Remarks
c  1. In this version, each shift is applied to all the subblocks of
c     the tridiagonal matrix H and not just to the submatrix that it
c     comes from. This routine assumes that the subdiagonal elements
c     of H that are stored in h(1:kev+np,1) are nonegative upon input
c     and enforce this condition upon output. This version incorporates
c     deflation. See code for documentation.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine mydsapps
     &   ( n, kev, np, v, ldv, h, ldh, resid, q, ldq,
     &     eigval, eigvec, ldeigvec, workd )
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
      integer    kev, ldeigvec, ldh, ldq, ldv, n, np
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           eigval(kev+np), eigvec(ldeigvec,kev+np),
     &           h(ldh,2), q(ldq,kev+np), resid(n),
     &           v(ldv,kev+np), workd(2*n)
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
      integer    i, j, p, kplusp, msglvl
      Double precision
     &           rnorm, betak
      Double precision
     &           drot(kev+1,kev+1), qrot(kev+1,kev+1), qk1(kev,kev),
     &           qfinal(kev+np,kev), vnew(n,kev)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, dscal, dlaset, dlacpy, dgemv, dgemm,
     &           arrowgivens, arscnd
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
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c
      call arscnd (t0)
      msglvl = msapps
c
      kplusp = kev + np
c
c     %----------------------------------------------%
c     | Quick return if there are no shifts to apply |
c     %----------------------------------------------%
c
      if (np .eq. 0) go to 9000
c
c     %--------------------------------------------------------%
c     | ARROWHEAD RESTART (replaces the classical implicit-    |
c     | shift bulge-chasing entirely): given the caller's full |
c     | eigendecomposition of the current KEV+NP tridiagonal   |
c     | (EIGVAL, EIGVEC), already arranged so the first KEV    |
c     | entries/columns are the desired eigenpairs, form the   |
c     | (KEV+1) by (KEV+1) arrowhead matrix built from those   |
c     | KEV desired eigenpairs plus the residual attachment,   |
c     | reduce it directly to tridiagonal form via one-way     |
c     | Givens chasing (ARROWGIVENS, following Zha (1992)), and|
c     | read the new H, V, and RESID off that result. See "A   |
c     | Unified View of Arrowhead Matrix Transformations and   |
c     | Lanczos Restarts", Baglama, Monette, Perovic (2026).   |
c     %--------------------------------------------------------%
c
c     %------------------------------------------------------%
c     | RNORM: the current residual norm, needed both as the |
c     | arrowhead spike scale and to normalize the residual  |
c     | direction when forming the updated residual below.   |
c     %------------------------------------------------------%
c
      rnorm = dnrm2(n, resid, 1)
c
c     %-----------------------------------------------------------%
c     | Directly build the ALREADY 180-degree-rotated (KEV+1) by  |
c     | (KEV+1) arrowhead matrix DROT -- mathematically identical |
c     | to building the natural "downward-pointing" arrowhead     |
c     | and then rotating it 180 degrees, without forming the     |
c     | un-rotated matrix                                         |
c     %-----------------------------------------------------------%
c
      call dlaset ('All', kev+1, kev+1, zero, zero, drot, kev+1)
c
      do 50 j = 2, kev+1
         drot(1,j) = rnorm * eigvec(kplusp, kev+2-j)
         drot(j,1) = drot(1,j)
   50 continue
      do 60 i = 2, kev+1
         drot(i,i) = eigval(kev+2-i)
   60 continue
c
c     %------------------------------------------------------------%
c     | Tridiagonalize DROT via the one-way Givens chasing scheme; |
c     | DROT is overwritten in place with the tridiagonal result.  |
c     %------------------------------------------------------------%
c
      call arrowgivens (kev+1, drot, kev+1, qrot, kev+1)
c
c     %--------------------------------------------------------%
c     | Read the new leading KEV by KEV tridiagonal H, and     |
c     | BETAK, directly out of the rotated result              |
c     %--------------------------------------------------------%
c
      do 70 i = 1, kev
         h(i,2) = drot(kev+2-i,kev+2-i)
   70 continue
      do 80 i = 2, kev
         h(i,1) = drot(kev+2-i,kev+3-i)
   80 continue
c
      betak = drot(1,2)
c
c     %--------------------------------------------------------%
c     | Build QFINAL = EIGVEC(:,1:kev) * QK1(1:kev,1:kev), the |
c     | KPLUSP by KEV transformation carrying V's current      |
c     | KPLUSP-column basis directly onto the new KEV-column   |
c     | basis (QK1 read out of QROT the same way T was read    |
c     | out of DROT above).                                    |
c     %--------------------------------------------------------%
c
      do 120 j = 1, kev
         do 110 p = 1, kev
            qk1(p,j) = qrot(kev+2-p,kev+2-j)
  110    continue
  120 continue
c
      call dgemm ('N', 'N', kplusp, kev, kev, one, eigvec, ldeigvec,
     &            qk1, kev, zero, qfinal, kplusp)
c
c     %-------------------------------------------------%
c     | Update V: V(:,1:kev) <- V(:,1:kplusp) * QFINAL. |
c     %-------------------------------------------------%
c
      call dgemm ('N', 'N', n, kev, kplusp, one, v, ldv, qfinal,
     &            kplusp, zero, vnew, n)
      call dlacpy ('All', n, kev, vnew, n, v, ldv)
c
c     %--------------------------------------------------------%
c     | Update the residual vector. RESID_new is BETAK times.  |
c     | the ORIGINAL residual direction                        |
c     %--------------------------------------------------------%
c
      if (rnorm .gt. zero) then
         call dscal (n, betak/rnorm, resid, 1)
      end if
c
      call arscnd (t1)
      tsapps = tsapps + (t1 - t0)
c
 9000 continue
      return
c
c     %---------------%
c     | End of dsapps |
c     %---------------%
c
      end
