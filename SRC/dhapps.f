c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dhapps
c
c\Description:
c  Given the Arnoldi factorization
c
c     A*V_{m} - V_{m}*H_{m} = r_{m}*e_{m}^T,   m = KEV+NP,
c
c  apply a thick restart that retains the KEV desired harmonic Ritz
c  vectors (see \Remarks), resulting in the Arnoldi factorization
c
c     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T,   k = KEV,
c
c  where HNEW_{k} is upper Hessenberg.
c
c\Usage:
c  call dhapps
c     ( N, KEV, NP, V, LDV, H, LDH, RESID, RNORM, Z, LDZ, U )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Problem size, i.e. dimension of matrix A.
c
c  KEV     Integer.  (INPUT)
c          Size of the updated factorization. The caller ensures that
c          Z(:,1:KEV) does not split a complex conjugate pair.
c
c  NP      Integer.  (INPUT)
c          Number of harmonic Ritz values being discarded.
c
c  V       Double precision N by (KEV+NP) array.  (INPUT/OUTPUT)
c          INPUT: V contains the current KEV+NP Arnoldi vectors.
c          OUTPUT: the updated Arnoldi vectors are in the first KEV
c          columns of V.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  H       Double precision (KEV+NP) by (KEV+NP) array.  (INPUT/OUTPUT)
c          INPUT: the current KEV+NP upper Hessenberg matrix.
c          OUTPUT: the updated KEV by KEV upper Hessenberg matrix is in
c          the leading KEV by KEV submatrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          INPUT: RESID contains the residual vector r_{m}.
c          OUTPUT: RESID is the updated residual vector rnew_{k}.
c
c  RNORM   Double precision scalar.  (INPUT)
c          B-norm of the input RESID.
c
c  Z       Double precision (KEV+NP) by (KEV+NP) array.  (INPUT)
c          Right Schur vectors of the harmonic pencil, ordered so the
c          first KEV columns span the desired harmonic Ritz vectors.
c          Only the first KEV columns are referenced.
c
c  LDZ     Integer.  (INPUT)
c          Leading dimension of Z exactly as declared in the calling
c          program.
c
c  U       Double precision array of length KEV+NP+1.  (INPUT)
c          Unit vector spanning the null space of HBAR^T, where
c          HBAR = [H; RNORM*e_{m}^T].
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
c     dorghr  LAPACK routine that generates the orthogonal matrix from
c             the reflectors of DGEHRD.
c     arscnd  ARPACK utility routine for timing.
c     dlaset  LAPACK matrix initialization routine.
c     dlacpy  LAPACK matrix copy routine.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     dgemm   Level 3 BLAS routine for matrix-matrix multiplication.
c     dtrmm   Level 3 BLAS routine for triangular matrix times matrix.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     daxpy   Level 1 BLAS that computes a vector triad.
c     dscal   Level 1 BLAS that scales a vector.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     dvout   ARPACK utility routine that prints vectors.
c
c\Remarks
c  1. Let HBAR = [H; RNORM*e_{m}^T] and
c     Zk = Z(:,1:KEV). Since Zk spans an invariant subspace of the
c     harmonic matrix, HBAR*Zk lies in the span of the KEV+1
c     orthonormal columns of PK = [[Zk; 0] p], where p is U
c     orthogonalized against [Zk; 0]. Hence
c         A*V*Zk = [V RESID/RNORM]*PK*HK,   HK = PK^T*HBAR*Zk,
c     with HK a dense (KEV+1) by KEV matrix.
c  2. The bordered matrix [HK 0] is rotated by 180 degrees and
c     transposed, then reduced to upper Hessenberg form by DGEHRD.
c     After rotating back, the orthogonal factor is diag(Qk,1) and
c     the last row of the reduced matrix is BETAK*e_{k}^T. Then
c         VNEW = V*Zk*Qk,  HNEW = Qk^T*HK(1:KEV,:)*Qk,
c         rnew = BETAK*[V RESID/RNORM]*p.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dhapps
     &   ( n, kev, np, v, ldv, h, ldh, resid, rnorm, z, ldz, u )
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
      integer    kev, ldh, ldv, ldz, n, np
      Double precision
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           h(ldh,kev+np), resid(n), u(kev+np+1), v(ldv,kev+np),
     &           z(ldz,kev+np)
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
     &           betak, pnorm, qwork(1)
c
c     %-----------------------%
c     | Local Array Arguments |
c     %-----------------------%
c
      Double precision
     &           pvec(kev+np+1), coef(kev), hz(kev+np+1,kev),
     &           hk(kev+1,kev), drot(kev+1,kev+1), hbuf(kev+1,kev+1),
     &           tau(kev), qk(kev,kev), qfinal(kev+np,kev+1)
      Double precision, allocatable, save :: housework(:), vnew(:,:)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, daxpy, dscal, dlaset, dlacpy, dgemv, dgemm,
     &           dtrmm, dgehrd, dorghr, arscnd, dvout
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
c     | Initialize the debug.h/stat.h COMMON variables used here |
c     | the first time this routine is called.                   |
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
c     %-----------------------------------------------%
c     | Quick return if there are no values to remove |
c     %-----------------------------------------------%
c
      if (np .eq. 0) go to 9000
c
c     %-------------------------------------------------%
c     | Workspace query for DGEHRD and DORGHR, with     |
c     | 8*(KEV+1) as a floor.                           |
c     %-------------------------------------------------%
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
c     %--------------------------------------------------%
c     | PVEC = U orthogonalized against [Zk; 0], with    |
c     | two passes of classical Gram-Schmidt.            |
c     %--------------------------------------------------%
c
      call dcopy (kplusp+1, u, 1, pvec, 1)
      do 10 i = 1, 2
         call dgemv ('T', kplusp, kev, one, z, ldz, pvec, 1,
     &               zero, coef, 1)
         call dgemv ('N', kplusp, kev, -one, z, ldz, coef, 1,
     &               one, pvec, 1)
   10 continue
      pnorm = dnrm2 (kplusp+1, pvec, 1)
      call dscal (kplusp+1, one/pnorm, pvec, 1)
c
c     %--------------------------------------------------%
c     | HZ = HBAR*Zk, HBAR = [H; RNORM*e_m^T], using the |
c     | upper triangle of H and then its subdiagonal.    |
c     %--------------------------------------------------%
c
      call dlacpy ('All', kplusp, kev, z, ldz, hz, kplusp+1)
      call dtrmm ('Left', 'Upper', 'No transpose', 'Non-unit', kplusp,
     &            kev, one, h, ldh, hz, kplusp+1)
      do 20 i = 1, kplusp-1
         call daxpy (kev, h(i+1,i), z(i,1), ldz, hz(i+1,1), kplusp+1)
   20 continue
      call dcopy (kev, z(kplusp,1), ldz, hz(kplusp+1,1), kplusp+1)
      call dscal (kev, rnorm, hz(kplusp+1,1), kplusp+1)
c
c     %-----------------------------------------%
c     | HK = PK^T*HZ, PK = [[Zk; 0] PVEC]        |
c     %-----------------------------------------%
c
      call dgemm ('T', 'N', kev, kev, kplusp, one, z, ldz,
     &            hz, kplusp+1, zero, hk, kev+1)
      call dgemv ('T', kplusp+1, kev, one, hz, kplusp+1, pvec, 1,
     &            zero, hk(kev+1,1), kev+1)
c
c     %--------------------------------------------------------%
c     | DROT is [HK 0] rotated by 180 degrees and transposed:  |
c     | DROT(i,j) = [HK 0](kev+2-j,kev+2-i). Its first row is  |
c     | zero.                                                  |
c     %--------------------------------------------------------%
c
      call dlaset ('All', kev+1, kev+1, zero, zero, drot, kev+1)
      do 40 j = 1, kev+1
         do 30 i = 2, kev+1
            drot(i,j) = hk(kev+2-j, kev+2-i)
   30    continue
   40 continue
c
c     %--------------------------------------------------------%
c     | Reduce DROT to upper Hessenberg form. Copy the         |
c     | Hessenberg entries to HBUF before DORGHR overwrites    |
c     | DROT with the orthogonal factor.                       |
c     %--------------------------------------------------------%
c
      call dgehrd (kev+1, 1, kev+1, drot, kev+1, tau,
     &             housework, lwork, houseinfo)
      call dlaset ('All', kev+1, kev+1, zero, zero, hbuf, kev+1)
      do 60 j = 1, kev+1
         do 50 i = 1, min(j+1,kev+1)
            hbuf(i,j) = drot(i,j)
   50    continue
   60 continue
      call dorghr (kev+1, 1, kev+1, drot, kev+1, tau,
     &             housework, lwork, houseinfo)
c
c     %--------------------------------------------------------%
c     | Rotate back: the new KEV by KEV upper Hessenberg H,    |
c     | its subdiagonal coupling BETAK, and Qk.                |
c     %--------------------------------------------------------%
c
      call dlaset ('All', kev, kev, zero, zero, h, ldh)
      do 80 j = 1, kev
         do 70 i = 1, min(j+1,kev)
            h(i,j) = hbuf(kev+2-j, kev+2-i)
   70    continue
   80 continue
      betak = hbuf(2,1)
c
      do 90 j = 1, kev
         call dcopy (kev, drot(2,kev+2-j), -1, qk(1,j), 1)
   90 continue
c
c     %--------------------------------------------------------%
c     | QFINAL = [Zk*Qk PVEC(1:m)], so that V*QFINAL holds the |
c     | new V and, in its last column, V*PVEC(1:m).            |
c     %--------------------------------------------------------%
c
      call dgemm ('N', 'N', kplusp, kev, kev, one, z, ldz,
     &            qk, kev, zero, qfinal, kplusp)
      call dcopy (kplusp, pvec, 1, qfinal(1,kev+1), 1)
c
      if (allocated(vnew)) then
         if (size(vnew,1) .ne. n .or. size(vnew,2) .lt. kev+1)
     &      deallocate (vnew)
      end if
      if (.not. allocated(vnew)) allocate (vnew(n,kev+1))
      call dgemm ('N', 'N', n, kev+1, kplusp, one, v, ldv, qfinal,
     &            kplusp, zero, vnew, n)
c
c     %--------------------------------------------------------%
c     | RESID = BETAK*w, where w = [V RESID/RNORM]*PVEC is a   |
c     | unit vector orthogonal to the new V.                   |
c     %--------------------------------------------------------%
c
      call dscal (n, betak*pvec(kplusp+1)/rnorm, resid, 1)
      call daxpy (n, betak, vnew(1,kev+1), 1, resid, 1)
c
      call dlacpy ('All', n, kev, vnew, n, v, ldv)
c
      if (msglvl .gt. 1) then
         call dvout (logfil, 1, [betak], ndigit,
     &        '_happs: updated residual norm')
      end if
c
      call arscnd (t1)
      tnapps = tnapps + (t1 - t0)
c
 9000 continue
      return
c
c     %---------------%
c     | End of dhapps |
c     %---------------%
c
      end
