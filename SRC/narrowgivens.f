c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: narrowgivens
c
c\Description:
c  Reduces an M by M nonsymmetric arrowhead matrix D (pointing
c  southeast: hub D(M,M), spike along the last row/column, quasi-
c  upper-triangular leading block D(1:M-1,1:M-1)) to upper Hessenberg
c  form using Givens rotations in a one-way chasing scheme in the
c  style of Zha (see reference 2 below).
c
c  Unlike ARROWGIVENS (the symmetric routine, whose input must be
c  pre-rotated to put the hub at (1,1)), no such rotation of D is
c  needed here: the spike is chased directly out of its natural
c  "downward-pointing" form.
c
c\Usage:
c  call narrowgivens
c     ( M, D, LDD, Q, LDQ )
c
c\Arguments
c  M       Integer.  (INPUT)
c          Order of the arrowhead matrix D.  Must be >= 3.
c
c  D       Double precision M by M array.  (INPUT/OUTPUT)
c          INPUT: D contains the nonsymmetric arrowhead matrix
c              D = | S(1:k,1:k)      c*Q(m,1:k)' |,    k = M-1,
c                  | c*Q(m,1:k)         delta    |
c          where S is (real-Schur) quasi-upper-triangular. Only the
c          leading K by K block, the spike D(M,1:K)/D(1:K,M), and the
c          hub D(M,M) are read; all other entries are assumed zero and
c          are never referenced.
c          OUTPUT: D contains the resulting upper Hessenberg matrix,
c          i.e. Hk = Q'*D_old*Q, with leading K by K block Hk upper
c          Hessenberg, D(K,K+1) = D(K+1,K) = a scalar, and
c          D(K+1,K+1) = delta (unchanged).
c
c  LDD     Integer.  (INPUT)
c          Leading dimension of D exactly as declared in the calling
c          program.
c
c  Q       Double precision M by M array.  (OUTPUT)
c          On output, Q contains the accumulated orthogonal
c          transformation such that D_new = Q' * D_old * Q. Q is
c          initialized to the identity internally; the caller does
c          not need to (and should not) pre-initialize it. Because no
c          rotation of D is needed, Q's leading K by K block is
c          exactly Qk (the transform used to obtain Hk = Qk'*S*Qk),
c          and Q's last row/column reduce to e_{K+1}.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
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
c  2. H. Zha, A two-way chasing scheme for reducing an arrowhead matrix
c     to tridiagonal form, J. Numer. Lin. Alg. Appl., 1 (1992), 49-57.
c
c\Routines called:
c     dlaset  LAPACK matrix initialization routine.
c     dlartg  LAPACK Givens rotation construction routine.
c
c\Author
c     James Baglama
c     Kyle Monette
c     Vasilije Perovic
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine narrowgivens
     &   ( m, d, ldd, q, ldq )
c
      integer    ldd, ldq, m
      Double precision
     &           d(ldd,m), q(ldq,m)
c
      Double precision
     &           zero, one, tol
      parameter (zero = 0.0D+0, one = 1.0D+0, tol = 1.0D-20)
c
      integer    c, cp1, i, r
      Double precision
     &           cs, sn, rr, t1, t2
c
      external   dlaset, dlartg
      intrinsic  abs
c
c     %--------------------------------%
c     | Initialize Q to the identity.  |
c     %--------------------------------%
c
      call dlaset ('All', m, m, zero, one, q, ldq)
c
      do 100 r = m, 3, -1
c
c        %-------------------------------------------------------------%
c        | Chase the spike out of row r, one column at a time,         |
c        | starting from the border column 1 and working right until   |
c        | the entry just before the Hessenberg subdiagonal (column    |
c        | r-1) is cleared.                                            |
c        %-------------------------------------------------------------%
c
         do 90 c = 1, r-2
c
            if (abs(d(r,c)) .lt. tol) go to 90
c
            cp1 = c + 1
c
c           %--------------------------------------------------------%
c           | Construct the rotation (cs,sn) that zeros D(r,c) using |
c           | D(r,c+1) as pivot; this is a right rotation applied to |
c           | columns c:c+1, cf. dlartg(a,b,cs,sn,rr):               |
c           |    cs*a + sn*b = rr,   -sn*a + cs*b = 0                |
c           %--------------------------------------------------------%
c
            call dlartg (d(r,c), d(r,cp1), cs, sn, rr)
c
c           %--------------------------------------------------------%
c           | Column update (rows 1:r only): D(1:r,c:c+1) <- D * G'  |
c           | where G = [[sn,-cs],[cs,sn]] is the rotation embedded  |
c           | in columns/rows c,c+1.  This drives D(r,c) to zero and |
c           | D(r,c+1) to rr.                                        |
c           %--------------------------------------------------------%
c
            do 10 i = 1, r
               t1 =  sn*d(i,c) - cs*d(i,cp1)
               t2 =  cs*d(i,c) + sn*d(i,cp1)
               d(i,c) = t1
               d(i,cp1) = t2
   10       continue
c
c           %--------------------------------------------------------%
c           | Row update (full column range 1:m): D(c:c+1,1:m) <- G  |
c           | * D(c:c+1,1:m).  The full range is required since the  |
c           | border row/column may still have the spike             |
c           %--------------------------------------------------------%
c
            do 20 i = 1, m
               t1 =  sn*d(c,i) - cs*d(cp1,i)
               t2 =  cs*d(c,i) + sn*d(cp1,i)
               d(c,i) = t1
               d(cp1,i) = t2
   20       continue
c
c           %--------------------------------------------------------%
c           | Accumulate: Q(:,c:c+1) <- Q(:,c:c+1) * G'              |
c           %--------------------------------------------------------%
c
            do 30 i = 1, m
               t1 =  sn*q(i,c) - cs*q(i,cp1)
               t2 =  cs*q(i,c) + sn*q(i,cp1)
               q(i,c) = t1
               q(i,cp1) = t2
   30       continue
c
   90    continue
  100 continue
c
      return
c
c     %-----------------------%
c     | End of narrowgivens   |
c     %-----------------------%
c
      end
