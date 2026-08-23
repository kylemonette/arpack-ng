c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dahtgv
c
c\Description:
c  Tridiagonalizes an M by M symmetric "upward-pointing" arrowhead
c  matrix A using the one-way Givens chasing scheme of
c  Zha (see reference 2 below).
c
c\Usage:
c  call dahtgv
c     ( M, A, LDA, QMAT, LDQ )
c
c\Arguments
c  M       Integer.  (INPUT)
c          Order of the arrowhead matrix A.
c
c  A       Double precision M by M array.  (INPUT/OUTPUT)
c          INPUT: A contains the upward-pointing arrowhead matrix,
c          e.g. for M=4:
c              |a1 b2 b3 b4|
c              |b2 a2      |
c              |b3    a3   |
c              |b4       a4|
c          OUTPUT: A contains the resulting symmetric tridiagonal
c          matrix - subdiagonal in column 1 starting at A(2,1), main
c          diagonal in column 2.
c
c  LDA     Integer.  (INPUT)
c          Leading dimension of A exactly as declared in the calling
c          program.
c
c  QMAT    Double precision M by M array.  (OUTPUT)
c          On output, QMAT contains the accumulated orthogonal
c          transformation such that A_new = QMAT' * A_old * QMAT.
c          QMAT is initialized to the identity internally.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of QMAT exactly as declared in the
c          calling program.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\References:
c  1. A Unified View of Arrowhead Matrix Transformations and Lanczos
c     Restarts, James Baglama, Kyle Monette, Vasilije Perovic, (2026).
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
      subroutine dahtgv
     &   ( m, a, lda, qmat, ldq )
c
      integer    lda, ldq, m
      Double precision
     &           a(lda,m), qmat(ldq,m)
c
      Double precision
     &           zero, one
      parameter (zero = 0.0D+0, one = 1.0D+0)
c
      integer    i, ii, j, jj, k, r, c, uppercol, ncols, cols(4), jhi
      Double precision
     &           cs, sn, rr, t1, t2
c
      external   dlaset, dlartg, drot
c
c     %---------------------------------%
c     | Initialize QMAT to the identity |
c     %---------------------------------%
c
      call dlaset ('All', m, m, zero, one, qmat, ldq)
c
      do 100 k = 1, m-2
c
c        %-------------------------------------------------------%
c        | Eliminate A(i+1,1) using A(i,1) as pivot, rows i,i+1. |
c        %-------------------------------------------------------%
c
         i = m - k
         call dlartg (a(i,1), a(i+1,1), cs, sn, rr)
c
c        %-------------------------------------------------------%
c        | The set of columns this rotation touches: column 1,   |
c        | plus the contiguous range i:uppercol.                 |
c        %-------------------------------------------------------%
c
         uppercol = min(m, i + min(k,2))
         ncols = 1
         cols(1) = 1
         do 10 jj = i, uppercol
            if (jj .ne. 1) then
               ncols = ncols + 1
               cols(ncols) = jj
            end if
   10    continue
c
c        %-----------------------------------------------------%
c        | Apply the rotation to rows i,i+1 at these columns.  |
c        %-----------------------------------------------------%
c
         do 20 jj = 1, ncols
            c = cols(jj)
            t1       =  cs*a(i,c) + sn*a(i+1,c)
            t2       = -sn*a(i,c) + cs*a(i+1,c)
            a(i,c)   = t1
            a(i+1,c) = t2
   20    continue
c
c        %--------------------------------------------------------%
c        | Apply the rotation (from the right) to columns i,i+1   |
c        | at these same rows.                                    |
c        %--------------------------------------------------------%
c
         do 30 jj = 1, ncols
            r = cols(jj)
            t1       = a(r,i)*cs + a(r,i+1)*sn
            t2       = -a(r,i)*sn + a(r,i+1)*cs
            a(r,i)   = t1
            a(r,i+1) = t2
   30    continue
c
c        %----------------------%
c        | Accumulate into Q.  |
c        %----------------------%
c
         call drot (m, qmat(1,i), 1, qmat(1,i+1), 1, cs, sn)
c
c        %--------------------------------%
c        | Chase the bulge created above. |
c        %--------------------------------%
c
         do 90 j = 1, k-1
            r = i + j
            c = i + j - 1
            call dlartg (a(r,c), a(r+1,c), cs, sn, rr)
c
c           %-------------------------------------%
c           | Apply to rows r,r+1, columns c:jhi. |
c           %-------------------------------------%
c
            jhi = min(m, c+3)
            do 50 jj = c, jhi
               t1       =  cs*a(r,jj) + sn*a(r+1,jj)
               t2       = -sn*a(r,jj) + cs*a(r+1,jj)
               a(r,jj)   = t1
               a(r+1,jj) = t2
   50       continue
c
c           %------------------------------------------%
c           | Apply to columns r,r+1, rows c:jhi.      |
c           %------------------------------------------%
c
            do 60 ii = c, jhi
               t1        = a(ii,r)*cs + a(ii,r+1)*sn
               t2        = -a(ii,r)*sn + a(ii,r+1)*cs
               a(ii,r)   = t1
               a(ii,r+1) = t2
   60       continue
c
            call drot (m, qmat(1,r), 1, qmat(1,r+1), 1, cs, sn)
c
   90    continue
c
  100 continue
c
      return
c
c     %-----------------%
c     | End of dahtgv   |
c     %-----------------%
c
      end
