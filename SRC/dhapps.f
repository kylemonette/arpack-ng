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
c     ( N, KEV, NP, V, LDV, H, LDH, RESID, RNORM, TSCHUR, LDTSCHUR,
c       QSCHUR, LDQSCHUR )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Problem size, i.e. dimension of matrix A.
c
c  KEV     Integer.  (INPUT)
c          Size of the updated factorization. The caller ensures that
c          the leading KEV+1 by KEV+1 block of TSCHUR does not split a
c          complex conjugate pair.
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
c  H       Double precision (KEV+NP) by (KEV+NP) array.  (OUTPUT)
c          The updated KEV by KEV upper Hessenberg matrix is in the
c          leading KEV by KEV submatrix.
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
c  TSCHUR  Double precision (KEV+NP+1) by (KEV+NP+1) array.  (INPUT)
c          Real Schur form of the matrix HP of \Remarks 1, reordered so
c          that its leading KEV+1 by KEV+1 block holds the KEV desired
c          harmonic Ritz values and the zero eigenvalue. Only this
c          block is referenced.
c
c  LDTSCHUR Integer.  (INPUT)
c          Leading dimension of TSCHUR exactly as declared in the
c          calling program.
c
c  QSCHUR  Double precision (KEV+NP+1) by (KEV+NP+1) array.  (INPUT)
c          Schur vectors corresponding to TSCHUR. Only the first KEV+1
c          columns are referenced.
c
c  LDQSCHUR Integer.  (INPUT)
c          Leading dimension of QSCHUR exactly as declared in the
c          calling program.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\References:
c  1. Krylov-Schur and Implicit Restarting: A Bordered Hessenberg
c     Perspective, James Baglama, Kyle Monette, Vasilije Perovic,
c     (2026).
c  2. A Unified View of Arrowhead Matrix Transformations and Lanczos
c     Restarts, James Baglama, Kyle Monette, Vasilije Perovic, (2026).
c
c\Routines called:
c     dgehrd  LAPACK routine that reduces a general matrix to upper
c             Hessenberg form via Householder reflectors.
c     dorghr  LAPACK routine that generates the orthogonal matrix from
c             the reflectors of DGEHRD.
c     arscnd  ARPACK utility routine for timing.
c     dlaset  LAPACK matrix initialization routine.
c     dlacpy  LAPACK matrix copy routine.
c     dgemm   Level 3 BLAS routine for matrix-matrix multiplication.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     daxpy   Level 1 BLAS that computes a vector triad.
c     dscal   Level 1 BLAS that scales a vector.
c     dvout   ARPACK utility routine that prints vectors.
c
c\Remarks
c  1. Let HBAR = [H; RNORM*e_m^T] and let U be the unit vector spanning
c     the null space of HBAR^T. The (m+1) by (m+1) upper Hessenberg
c     matrix
c         HP = HBAR*[I, E],   E = -U(1:m)/U(m+1),
c     has the harmonic Ritz values and zero as its eigenvalues, and U
c     is its eigenvector for zero. Keeping the KEV desired harmonic
c     Ritz values and zero, with Schur vectors Zk = QSCHUR(:,1:KEV+1),
c     is an exact-shift restart of HP.
c  2. The (KEV+2) by (KEV+2) bordered matrix
c         D = | TSCHUR(1:KEV+1,1:KEV+1)  0 |
c             | QSCHUR(m+1,1:KEV+1)      0 |
c     is rotated by 180 degrees and transposed, then reduced to upper
c     Hessenberg form by DGEHRD. After rotating back, D = Qd*Hd*Qd^T
c     with Qd(KEV+2,KEV+2) = 1. Then X = Zk*Qd(1:KEV+1,1:KEV+1)
c     satisfies X(m+1,1:KEV) = 0, and
c         VNEW = V*X(1:m,1:KEV),  HNEW = Hd(1:KEV,1:KEV),
c         rnew = Hd(KEV+1,KEV)*[V RESID/RNORM]*X(:,KEV+1).
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dhapps
     &   ( n, kev, np, v, ldv, h, ldh, resid, rnorm, tschur, ldtschur,
     &     qschur, ldqschur )
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
      Double precision
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           h(ldh,kev+np), qschur(ldqschur,kev+np+1), resid(n),
     &           tschur(ldtschur,kev+np+1), v(ldv,kev+np)
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
      integer    i, j, kplusp, kp, nd, msglvl, houseinfo, lwork
      Double precision
     &           betak, qwork(1)
c
c     %-----------------------%
c     | Local Array Arguments |
c     %-----------------------%
c
      Double precision
     &           dflip(kev+2,kev+2), tau(kev+1), qk(kev+1,kev+1),
     &           xk(kev+np+1,kev+1)
      Double precision, allocatable, save :: housework(:), vnew(:,:)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, daxpy, dscal, dlaset, dlacpy, dgemm,
     &           dgehrd, dorghr, arscnd, dvout
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
      msglvl = mnapps
c
      kplusp = kev + np
      kp = kev + 1
      nd = kev + 2
c
c     %-----------------------------------------------%
c     | Quick return if there are no values to remove |
c     %-----------------------------------------------%
c
      if (np .eq. 0) go to 9000
c
c     %-------------------------------------------------%
c     | Workspace query for DGEHRD and DORGHR, with     |
c     | 8*ND as a floor.                                |
c     %-------------------------------------------------%
c
      call dgehrd (nd, 1, nd, dflip, nd, tau, qwork, -1, houseinfo)
      lwork = int(qwork(1))
      call dorghr (nd, 1, nd, dflip, nd, tau, qwork, -1, houseinfo)
      lwork = max(lwork, int(qwork(1)), 8*nd)
      if (allocated(housework)) then
         if (size(housework) .lt. lwork) deallocate (housework)
      end if
      if (.not. allocated(housework)) allocate (housework(lwork))
c
c     %--------------------------------------------------------%
c     | DFLIP is D of \Remarks 2 rotated by 180 degrees and    |
c     | transposed: DFLIP(i,j) = D(nd+1-j,nd+1-i). Its first   |
c     | row is zero and its first column holds the last row of |
c     | QSCHUR(:,1:KP) in reverse order.                       |
c     %--------------------------------------------------------%
c
      call dlaset ('All', 1, nd, zero, zero, dflip, nd)
      do 10 i = 2, nd
         dflip(i,1) = qschur(kplusp+1, nd+1-i)
   10 continue
      do 30 j = 2, nd
         do 20 i = 2, nd
            dflip(i,j) = tschur(nd+1-j, nd+1-i)
   20    continue
   30 continue
c
c     %---------------------------------------------------------%
c     | Reduce DFLIP to upper Hessenberg form, and rotate back  |
c     | to obtain the new KEV by KEV upper Hessenberg H and     |
c     | BETAK = Hd(KEV+1,KEV). These are read from entries on   |
c     | or above the subdiagonal of DFLIP, before DORGHR        |
c     | overwrites DFLIP with the orthogonal factor, from which |
c     | Qd(1:KP,1:KP) is then read.                             |
c     %---------------------------------------------------------%
c
      call dgehrd (nd, 1, nd, dflip, nd, tau, housework, lwork,
     &             houseinfo)
c
      if (kev .gt. 2) then
         call dlaset ('Lower', kev-2, kev-2, zero, zero, h(3,1), ldh)
      end if
      do 50 j = 1, kev
         do 40 i = 1, min(j+1,kev)
            h(i,j) = dflip(nd+1-j, nd+1-i)
   40    continue
   50 continue
      betak = dflip(3,2)
c
      call dorghr (nd, 1, nd, dflip, nd, tau, housework, lwork,
     &             houseinfo)
c
      do 60 j = 1, kp
         call dcopy (kp, dflip(2,nd+1-j), -1, qk(1,j), 1)
   60 continue
c
c     %------------------------------------------------------%
c     | XK = QSCHUR(:,1:KP)*Qd(1:KP,1:KP), and V*XK(1:m,:)   |
c     | holds the new V and, in its last column, the part of |
c     | the next basis vector in the range of V.             |
c     %------------------------------------------------------%
c
      call dgemm ('N', 'N', kplusp+1, kp, kp, one, qschur, ldqschur,
     &            qk, kp, zero, xk, kplusp+1)
c
      if (allocated(vnew)) then
         if (size(vnew,1) .ne. n .or. size(vnew,2) .lt. kp)
     &      deallocate (vnew)
      end if
      if (.not. allocated(vnew)) allocate (vnew(n,kp))
      call dgemm ('N', 'N', n, kp, kplusp, one, v, ldv, xk, kplusp+1,
     &            zero, vnew, n)
c
c     %------------------------------------------------------%
c     | RESID = BETAK*w, where w = [V RESID/RNORM]*XK(:,KP)  |
c     | is a unit vector orthogonal to the new V.            |
c     %------------------------------------------------------%
c
      if (rnorm .gt. zero) then
         call dscal (n, betak*xk(kplusp+1,kp)/rnorm, resid, 1)
      else
         call dscal (n, zero, resid, 1)
      end if
      call daxpy (n, betak, vnew(1,kp), 1, resid, 1)
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
