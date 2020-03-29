C Output from Public domain Ratfor, version 1.0
      subroutine rqfnb(n,p,a,y,rhs,d,u,beta,eps,wn,wp,nit,info)
      integer n,p,info,nit(3)
      double precision a(p,n),y(n),rhs(p),d(n),u(n),wn(n,9),wp(p,p+3)
      double precision one,beta,eps
      parameter( one = 1.0d0)
      call lpfnb(n,p,a,y,rhs,d,u,beta,eps,wn(1,1),wn(1,2), wp(1,1),wn(1,
     *3),wn(1,4),wn(1,5),wn(1,6), wp(1,2),wn(1,7),wn(1,8),wn(1,9),wp(1,3
     *),wp(1,4),nit,info)
      return
      end
      subroutine lpfnb(n,p,a,c,b,d,u,beta,eps,x,s,y,z,w, dx,ds,dy,dz,dw,
     *dr,rhs,ada,nit,info)
      integer n,p,pp,i,info,nit(3),maxit
      double precision a(p,n),c(n),b(p)
      double precision zero,one,mone,big,ddot,dmax1,dmin1,dxdz,dsdw
      double precision deltap,deltad,beta,eps,mu,gap,g
      double precision x(n),u(n),s(n),y(p),z(n),w(n),d(n),rhs(p),ada(p,p
     *)
      double precision dx(n),ds(n),dy(p),dz(n),dw(n),dr(n)
      parameter( zero = 0.0d0)
      parameter( one = 1.0d0)
      parameter( mone = -1.0d0)
      parameter( big = 1.0d+20)
      parameter( maxit = 500)
      nit(1)=0
      nit(2)=0
      nit(3)=n
      pp=p*p
      call dgemv('N',p,n,one,a,p,c,1,zero,y,1)
      do23000 i=1,n
      d(i)=one
23000 continue
23001 continue
      call stepy(n,p,a,d,y,ada,info)
      if(info .ne. 0)then
      return
      endif
      call dcopy(n,c,1,s,1)
      call dgemv('T',p,n,mone,a,p,y,1,one,s,1)
      do23004 i=1,n
      if(dabs(s(i)).lt.eps)then
      z(i)=dmax1(s(i), zero) + eps
      w(i)=dmax1(-s(i),zero) + eps
      else
      z(i)=dmax1(s(i), zero)
      w(i)=dmax1(-s(i),zero)
      endif
      s(i)=u(i)-x(i)
23004 continue
23005 continue
      gap = ddot(n,z,1,x,1)+ddot(n,w,1,s,1)
23008 if(gap .gt. eps .and. nit(1).lt.maxit)then
      nit(1)=nit(1)+1
      do23010 i = 1,n
      d(i) = one/(z(i)/x(i) + w(i)/s(i))
      ds(i)=z(i)-w(i)
      dz(i)=d(i)*ds(i)
23010 continue
23011 continue
      call dcopy(p,b,1,dy,1)
      call dgemv('N',p,n,mone,a,p,x,1,one,dy,1)
      call dgemv('N',p,n,one,a,p,dz,1,one,dy,1)
      call dcopy(p,dy,1,rhs,1)
      call stepy(n,p,a,d,dy,ada,info)
      if(info .ne. 0)then
      return
      endif
      call dgemv('T',p,n,one,a,p,dy,1,mone,ds,1)
      deltap=big
      deltad=big
      do23014 i=1,n
      dx(i)=d(i)*ds(i)
      ds(i)=-dx(i)
      dz(i)=-z(i)*(dx(i)/x(i) + one)
      dw(i)=-w(i)*(ds(i)/s(i) + one)
      if(dx(i).lt.0)then
      deltap=dmin1(deltap,-x(i)/dx(i))
      endif
      if(ds(i).lt.0)then
      deltap=dmin1(deltap,-s(i)/ds(i))
      endif
      if(dz(i).lt.0)then
      deltad=dmin1(deltad,-z(i)/dz(i))
      endif
      if(dw(i).lt.0)then
      deltad=dmin1(deltad,-w(i)/dw(i))
      endif
23014 continue
23015 continue
      deltap=dmin1(beta*deltap,one)
      deltad=dmin1(beta*deltad,one)
      if(min(deltap,deltad) .lt. one)then
      nit(2)=nit(2)+1
      mu = ddot(n,x,1,z,1)+ddot(n,s,1,w,1)
      g = mu + deltap*ddot(n,dx,1,z,1)+ deltad*ddot(n,dz,1,x,1) + deltap
     **deltad*ddot(n,dz,1,dx,1)+ deltap*ddot(n,ds,1,w,1)+ deltad*ddot(n,
     *dw,1,s,1) + deltap*deltad*ddot(n,ds,1,dw,1)
      mu = mu * ((g/mu)**3) /dfloat(2*n)
      do23026 i=1,n
      dr(i)=d(i)*(mu*(1/s(i)-1/x(i))+ dx(i)*dz(i)/x(i)-ds(i)*dw(i)/s(i))
23026 continue
23027 continue
      call dswap(p,rhs,1,dy,1)
      call dgemv('N',p,n,one,a,p,dr,1,one,dy,1)
      call dpotrs('U',p,1,ada,p,dy,p,info)
      call dgemv('T',p,n,one,a,p,dy,1,zero,u,1)
      deltap=big
      deltad=big
      do23028 i=1,n
      dxdz = dx(i)*dz(i)
      dsdw = ds(i)*dw(i)
      dx(i)= d(i)*(u(i)-z(i)+w(i))-dr(i)
      ds(i)= -dx(i)
      dz(i)= -z(i)+(mu - z(i)*dx(i) - dxdz)/x(i)
      dw(i)= -w(i)+(mu - w(i)*ds(i) - dsdw)/s(i)
      if(dx(i).lt.0)then
      deltap=dmin1(deltap,-x(i)/dx(i))
      endif
      if(ds(i).lt.0)then
      deltap=dmin1(deltap,-s(i)/ds(i))
      endif
      if(dz(i).lt.0)then
      deltad=dmin1(deltad,-z(i)/dz(i))
      endif
      if(dw(i).lt.0)then
      deltad=dmin1(deltad,-w(i)/dw(i))
      endif
23028 continue
23029 continue
      deltap=dmin1(beta*deltap,one)
      deltad=dmin1(beta*deltad,one)
      endif
      call daxpy(n,deltap,dx,1,x,1)
      call daxpy(n,deltap,ds,1,s,1)
      call daxpy(p,deltad,dy,1,y,1)
      call daxpy(n,deltad,dz,1,z,1)
      call daxpy(n,deltad,dw,1,w,1)
      gap = ddot(n,z,1,x,1)+ddot(n,w,1,s,1)
      goto 23008
      endif
23009 continue
      call daxpy(n,mone,w,1,z,1)
      call dswap(n,z,1,x,1)
      return
      end
      subroutine stepy(n,p,a,d,b,ada,info)
      integer n,p,pp,i,info
      double precision a(p,n),b(p),d(n),ada(p,p),zero
      parameter( zero = 0.0d0)
      pp=p*p
      do23038 j=1,p
      do23040 k=1,p
      ada(j,k)=zero
23040 continue
23041 continue
23038 continue
23039 continue
      do23042 i=1,n
      call dsyr('U',p,d(i),a(1,i),1,ada,p)
23042 continue
23043 continue
      call dposv('U',p,1,ada,p,b,p,info)
      return
      end

C Adding required subroutines from LAPACK and BLAS (netlib.org/LAPACK)


C DAPXY
      SUBROUTINE daxpy(N,DA,DX,INCX,DY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC mod
*     ..
      IF (n.LE.0) RETURN
      IF (da.EQ.0.0d0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
         m = mod(n,4)
         IF (m.NE.0) THEN
            DO i = 1,m
               dy(i) = dy(i) + da*dx(i)
            END DO
         END IF
         IF (n.LT.4) RETURN
         mp1 = m + 1
         DO i = mp1,n,4
            dy(i) = dy(i) + da*dx(i)
            dy(i+1) = dy(i+1) + da*dx(i+1)
            dy(i+2) = dy(i+2) + da*dx(i+2)
            dy(i+3) = dy(i+3) + da*dx(i+3)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
          dy(iy) = dy(iy) + da*dx(ix)
          ix = ix + incx
          iy = iy + incy
         END DO
      END IF
      RETURN
      END

C DCOPY

      SUBROUTINE dcopy(N,DX,INCX,DY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC mod
*     ..
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
         m = mod(n,7)
         IF (m.NE.0) THEN
            DO i = 1,m
               dy(i) = dx(i)
            END DO
            IF (n.LT.7) RETURN
         END IF
         mp1 = m + 1
         DO i = mp1,n,7
            dy(i) = dx(i)
            dy(i+1) = dx(i+1)
            dy(i+2) = dx(i+2)
            dy(i+3) = dx(i+3)
            dy(i+4) = dx(i+4)
            dy(i+5) = dx(i+5)
            dy(i+6) = dx(i+6)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            dy(iy) = dx(ix)
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      RETURN
      END

C DDOT

      DOUBLE PRECISION FUNCTION ddot(N,DX,INCX,DY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER incx,incy,n
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION dx(*),dy(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION dtemp
      INTEGER i,ix,iy,m,mp1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC mod
*     ..
      ddot = 0.0d0
      dtemp = 0.0d0
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
         m = mod(n,5)
         IF (m.NE.0) THEN
            DO i = 1,m
               dtemp = dtemp + dx(i)*dy(i)
            END DO
            IF (n.LT.5) THEN
               ddot=dtemp
            RETURN
            END IF
         END IF
         mp1 = m + 1
         DO i = mp1,n,5
          dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) +
     $            dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            dtemp = dtemp + dx(ix)*dy(iy)
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      ddot = dtemp
      RETURN
      END

C DGMEMV plus dependencies 

      SUBROUTINE dgemv(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*
*  -- Reference BLAS level2 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      parameter(one=1.0d+0,zero=0.0d+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC max
*     ..
*
*     Test the input parameters.
*
      info = 0
      IF (.NOT.lsame(trans,'N') .AND. .NOT.lsame(trans,'T') .AND.
     +    .NOT.lsame(trans,'C')) THEN
          info = 1
      ELSE IF (m.LT.0) THEN
          info = 2
      ELSE IF (n.LT.0) THEN
          info = 3
      ELSE IF (lda.LT.max(1,m)) THEN
          info = 6
      ELSE IF (incx.EQ.0) THEN
          info = 8
      ELSE IF (incy.EQ.0) THEN
          info = 11
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DGEMV ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((m.EQ.0) .OR. (n.EQ.0) .OR.
     +    ((alpha.EQ.zero).AND. (beta.EQ.one))) RETURN
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF (lsame(trans,'N')) THEN
          lenx = n
          leny = m
      ELSE
          lenx = m
          leny = n
      END IF
      IF (incx.GT.0) THEN
          kx = 1
      ELSE
          kx = 1 - (lenx-1)*incx
      END IF
      IF (incy.GT.0) THEN
          ky = 1
      ELSE
          ky = 1 - (leny-1)*incy
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF (beta.NE.one) THEN
          IF (incy.EQ.1) THEN
              IF (beta.EQ.zero) THEN
                  DO 10 i = 1,leny
                      y(i) = zero
   10             CONTINUE
              ELSE
                  DO 20 i = 1,leny
                      y(i) = beta*y(i)
   20             CONTINUE
              END IF
          ELSE
              iy = ky
              IF (beta.EQ.zero) THEN
                  DO 30 i = 1,leny
                      y(iy) = zero
                      iy = iy + incy
   30             CONTINUE
              ELSE
                  DO 40 i = 1,leny
                      y(iy) = beta*y(iy)
                      iy = iy + incy
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (alpha.EQ.zero) RETURN
      IF (lsame(trans,'N')) THEN
*
*        Form  y := alpha*A*x + y.
*
          jx = kx
          IF (incy.EQ.1) THEN
              DO 60 j = 1,n
                  temp = alpha*x(jx)
                  DO 50 i = 1,m
                      y(i) = y(i) + temp*a(i,j)
   50             CONTINUE
                  jx = jx + incx
   60         CONTINUE
          ELSE
              DO 80 j = 1,n
                  temp = alpha*x(jx)
                  iy = ky
                  DO 70 i = 1,m
                      y(iy) = y(iy) + temp*a(i,j)
                      iy = iy + incy
   70             CONTINUE
                  jx = jx + incx
   80         CONTINUE
          END IF
      ELSE
*
*        Form  y := alpha*A**T*x + y.
*
          jy = ky
          IF (incx.EQ.1) THEN
              DO 100 j = 1,n
                  temp = zero
                  DO 90 i = 1,m
                      temp = temp + a(i,j)*x(i)
   90             CONTINUE
                  y(jy) = y(jy) + alpha*temp
                  jy = jy + incy
  100         CONTINUE
          ELSE
              DO 120 j = 1,n
                  temp = zero
                  ix = kx
                  DO 110 i = 1,m
                      temp = temp + a(i,j)*x(ix)
                      ix = ix + incx
  110             CONTINUE
                  y(jy) = y(jy) + alpha*temp
                  jy = jy + incy
  120         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DGEMV .
*
      END

      LOGICAL FUNCTION lsame(CA,CB)
*
*  -- Reference BLAS level1 routine (version 3.1) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER ca,cb
*     ..
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC ichar
*     ..
*     .. Local Scalars ..
      INTEGER inta,intb,zcode
*     ..
*
*     Test if the characters are equal
*
      lsame = ca .EQ. cb
      IF (lsame) RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      zcode = ichar('Z')
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      inta = ichar(ca)
      intb = ichar(cb)
*
      IF (zcode.EQ.90 .OR. zcode.EQ.122) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
          IF (inta.GE.97 .AND. inta.LE.122) inta = inta - 32
          IF (intb.GE.97 .AND. intb.LE.122) intb = intb - 32
*
      ELSE IF (zcode.EQ.233 .OR. zcode.EQ.169) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
          IF (inta.GE.129 .AND. inta.LE.137 .OR.
     +        inta.GE.145 .AND. inta.LE.153 .OR.
     +        inta.GE.162 .AND. inta.LE.169) inta = inta + 64
          IF (intb.GE.129 .AND. intb.LE.137 .OR.
     +        intb.GE.145 .AND. intb.LE.153 .OR.
     +        intb.GE.162 .AND. intb.LE.169) intb = intb + 64
*
      ELSE IF (zcode.EQ.218 .OR. zcode.EQ.250) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
          IF (inta.GE.225 .AND. inta.LE.250) inta = inta - 32
          IF (intb.GE.225 .AND. intb.LE.250) intb = intb - 32
      END IF
      lsame = inta .EQ. intb
*
*     RETURN
*
*     End of LSAME
*
      END

      SUBROUTINE xerbla( SRNAME, INFO )
*
*  -- Reference BLAS level1 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER*(*)      SRNAME
      INTEGER            INFO
*     ..
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          len_trim
*     ..
*     .. Executable Statements ..
*
      WRITE( *, fmt = 9999 )srname( 1:len_trim( srname ) ), info
*
      stop
*
 9999 FORMAT( ' ** On entry to ', a, ' parameter number ', i2, ' had ',
     $      'an illegal value' )
*
*     End of XERBLA
*
      END

C DPOSV plus dependencies

      SUBROUTINE dposv( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
*
*  -- LAPACK driver routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           dpotrf, dpotrs, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      info = 0
      IF( .NOT.lsame( uplo, 'U' ) .AND. .NOT.lsame( uplo, 'L' ) ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( nrhs.LT.0 ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      ELSE IF( ldb.LT.max( 1, n ) ) THEN
         info = -7
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DPOSV ', -info )
         RETURN
      END IF
*
*     Compute the Cholesky factorization A = U**T*U or A = L*L**T.
*
      CALL dpotrf( uplo, n, a, lda, info )
      IF( info.EQ.0 ) THEN
*
*        Solve the system A*X = B, overwriting B with X.
*
         CALL dpotrs( uplo, n, nrhs, a, lda, b, ldb, info )
*
      END IF
      RETURN
*
*     End of DPOSV
*
      END

      SUBROUTINE dpotrf ( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK computational routine (version 3.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      parameter( one = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, JB, NB
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           lsame, ilaenv
*     ..
*     .. External Subroutines ..
      EXTERNAL           dgemm, dpotf2, dsyrk, dtrsm, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      info = 0
      upper = lsame( uplo, 'U' )
      IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -4
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DPOTRF', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      nb = ilaenv( 1, 'DPOTRF', uplo, n, -1, -1, -1 )
      IF( nb.LE.1 .OR. nb.GE.n ) THEN
*
*        Use unblocked code.
*
         CALL dpotf2( uplo, n, a, lda, info )
      ELSE
*
*        Use blocked code.
*
         IF( upper ) THEN
*
*           Compute the Cholesky factorization A = U'*U.
*
            DO 10 j = 1, n, nb
*
*              Update and factorize the current diagonal block and test
*              for non-positive-definiteness.
*
               jb = min( nb, n-j+1 )
 
               CALL dpotf2( 'Upper', jb, a( j, j ), lda, info )
 
               IF( info.NE.0 )
     $            GO TO 30
 
               IF( j+jb.LE.n ) THEN
*
*                 Updating the trailing submatrix.
*
                  CALL dtrsm( 'Left', 'Upper', 'Transpose', 'Non-unit',
     $                        jb, n-j-jb+1, one, a( j, j ), lda,
     $                        a( j, j+jb ), lda )
                  CALL dsyrk( 'Upper', 'Transpose', n-j-jb+1, jb, -one,
     $                        a( j, j+jb ), lda,
     $                        one, a( j+jb, j+jb ), lda )
               END IF
   10       CONTINUE
*
         ELSE
*
*           Compute the Cholesky factorization A = L*L'.
*
            DO 20 j = 1, n, nb
*
*              Update and factorize the current diagonal block and test
*              for non-positive-definiteness.
*
               jb = min( nb, n-j+1 )
 
               CALL dpotf2( 'Lower', jb, a( j, j ), lda, info )
 
               IF( info.NE.0 )
     $            GO TO 30
 
               IF( j+jb.LE.n ) THEN
*
*                Updating the trailing submatrix.
*
                 CALL dtrsm( 'Right', 'Lower', 'Transpose', 'Non-unit',
     $                       n-j-jb+1, jb, one, a( j, j ), lda,
     $                       a( j+jb, j ), lda )
 
                 CALL dsyrk( 'Lower', 'No Transpose', n-j-jb+1, jb,
     $                       -one, a( j+jb, j ), lda,
     $                       one, a( j+jb, j+jb ), lda )
               END IF
   20       CONTINUE
         END IF
      END IF
      GO TO 40
*
   30 CONTINUE
      info = info + j - 1
*
   40 CONTINUE
      RETURN
*
*     End of DPOTRF
*
      END

      SUBROUTINE dpotf2( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      parameter( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J
      DOUBLE PRECISION   AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      DOUBLE PRECISION   DDOT
      EXTERNAL           lsame, ddot, disnan
*     ..
*     .. External Subroutines ..
      EXTERNAL           dgemv, dscal, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, sqrt
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      info = 0
      upper = lsame( uplo, 'U' )
      IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -4
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DPOTF2', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
*
      IF( upper ) THEN
*
*        Compute the Cholesky factorization A = U**T *U.
*
         DO 10 j = 1, n
*
*           Compute U(J,J) and test for non-positive-definiteness.
*
            ajj = a( j, j ) - ddot( j-1, a( 1, j ), 1, a( 1, j ), 1 )
            IF( ajj.LE.zero.OR.disnan( ajj ) ) THEN
               a( j, j ) = ajj
               GO TO 30
            END IF
            ajj = sqrt( ajj )
            a( j, j ) = ajj
*
*           Compute elements J+1:N of row J.
*
            IF( j.LT.n ) THEN
               CALL dgemv( 'Transpose', j-1, n-j, -one, a( 1, j+1 ),
     $                     lda, a( 1, j ), 1, one, a( j, j+1 ), lda )
               CALL dscal( n-j, one / ajj, a( j, j+1 ), lda )
            END IF
   10    CONTINUE
      ELSE
*
*        Compute the Cholesky factorization A = L*L**T.
*
         DO 20 j = 1, n
*
*           Compute L(J,J) and test for non-positive-definiteness.
*
            ajj = a( j, j ) - ddot( j-1, a( j, 1 ), lda, a( j, 1 ),
     $            lda )
            IF( ajj.LE.zero.OR.disnan( ajj ) ) THEN
               a( j, j ) = ajj
               GO TO 30
            END IF
            ajj = sqrt( ajj )
            a( j, j ) = ajj
*
*           Compute elements J+1:N of column J.
*
            IF( j.LT.n ) THEN
               CALL dgemv( 'No transpose', n-j, j-1, -one, a( j+1, 1 ),
     $                     lda, a( j, 1 ), lda, one, a( j+1, j ), 1 )
               CALL dscal( n-j, one / ajj, a( j+1, j ), 1 )
            END IF
   20    CONTINUE
      END IF
      GO TO 40
*
   30 CONTINUE
      info = j
*
   40 CONTINUE
      RETURN
*
*     End of DPOTF2
*
      END

      LOGICAL FUNCTION disnan( DIN )
*
*  -- LAPACK auxiliary routine (version 3.7.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2017
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION, INTENT(IN) :: din
*     ..
*
*  =====================================================================
*
*  .. External Functions ..
      LOGICAL dlaisnan
      EXTERNAL dlaisnan
*  ..
*  .. Executable Statements ..
      disnan = dlaisnan(din,din)
      RETURN
      END

      LOGICAL FUNCTION dlaisnan( DIN1, DIN2 )
*
*  -- LAPACK auxiliary routine (version 3.7.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2017
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION, INTENT(IN) :: din1, din2
*     ..
*
*  =====================================================================
*
*  .. Executable Statements ..
      dlaisnan = (din1.NE.din2)
      RETURN
      END
       
      SUBROUTINE dscal(N,DA,DX,INCX)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,M,MP1,NINCX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC mod
*     ..
      IF (n.LE.0 .OR. incx.LE.0) RETURN
      IF (incx.EQ.1) THEN
*
*        code for increment equal to 1
*
*
*        clean-up loop
*
         m = mod(n,5)
         IF (m.NE.0) THEN
            DO i = 1,m
               dx(i) = da*dx(i)
            END DO
            IF (n.LT.5) RETURN
         END IF
         mp1 = m + 1
         DO i = mp1,n,5
            dx(i) = da*dx(i)
            dx(i+1) = da*dx(i+1)
            dx(i+2) = da*dx(i+2)
            dx(i+3) = da*dx(i+3)
            dx(i+4) = da*dx(i+4)
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         nincx = n*incx
         DO i = 1,nincx,incx
            dx(i) = da*dx(i)
         END DO
      END IF
      RETURN
      END

      SUBROUTINE dsyrk(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
*
*  -- Reference BLAS level3 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDC,N
      CHARACTER TRANS,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),C(LDC,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC max
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,L,NROWA
      LOGICAL UPPER
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      parameter(one=1.0d+0,zero=0.0d+0)
*     ..
*
*     Test the input parameters.
*
      IF (lsame(trans,'N')) THEN
          nrowa = n
      ELSE
          nrowa = k
      END IF
      upper = lsame(uplo,'U')
*
      info = 0
      IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
          info = 1
      ELSE IF ((.NOT.lsame(trans,'N')) .AND.
     +         (.NOT.lsame(trans,'T')) .AND.
     +         (.NOT.lsame(trans,'C'))) THEN
          info = 2
      ELSE IF (n.LT.0) THEN
          info = 3
      ELSE IF (k.LT.0) THEN
          info = 4
      ELSE IF (lda.LT.max(1,nrowa)) THEN
          info = 7
      ELSE IF (ldc.LT.max(1,n)) THEN
          info = 10
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DSYRK ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((n.EQ.0) .OR. (((alpha.EQ.zero).OR.
     +    (k.EQ.0)).AND. (beta.EQ.one))) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (alpha.EQ.zero) THEN
          IF (upper) THEN
              IF (beta.EQ.zero) THEN
                  DO 20 j = 1,n
                      DO 10 i = 1,j
                          c(i,j) = zero
   10                 CONTINUE
   20             CONTINUE
              ELSE
                  DO 40 j = 1,n
                      DO 30 i = 1,j
                          c(i,j) = beta*c(i,j)
   30                 CONTINUE
   40             CONTINUE
              END IF
          ELSE
              IF (beta.EQ.zero) THEN
                  DO 60 j = 1,n
                      DO 50 i = j,n
                          c(i,j) = zero
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 80 j = 1,n
                      DO 70 i = j,n
                          c(i,j) = beta*c(i,j)
   70                 CONTINUE
   80             CONTINUE
              END IF
          END IF
          RETURN
      END IF
*
*     Start the operations.
*
      IF (lsame(trans,'N')) THEN
*
*        Form  C := alpha*A*A**T + beta*C.
*
          IF (upper) THEN
              DO 130 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 90 i = 1,j
                          c(i,j) = zero
   90                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 100 i = 1,j
                          c(i,j) = beta*c(i,j)
  100                 CONTINUE
                  END IF
                  DO 120 l = 1,k
                      IF (a(j,l).NE.zero) THEN
                          temp = alpha*a(j,l)
                          DO 110 i = 1,j
                              c(i,j) = c(i,j) + temp*a(i,l)
  110                     CONTINUE
                      END IF
  120             CONTINUE
  130         CONTINUE
          ELSE
              DO 180 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 140 i = j,n
                          c(i,j) = zero
  140                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 150 i = j,n
                          c(i,j) = beta*c(i,j)
  150                 CONTINUE
                  END IF
                  DO 170 l = 1,k
                      IF (a(j,l).NE.zero) THEN
                          temp = alpha*a(j,l)
                          DO 160 i = j,n
                              c(i,j) = c(i,j) + temp*a(i,l)
  160                     CONTINUE
                      END IF
  170             CONTINUE
  180         CONTINUE
          END IF
      ELSE
*
*        Form  C := alpha*A**T*A + beta*C.
*
          IF (upper) THEN
              DO 210 j = 1,n
                  DO 200 i = 1,j
                      temp = zero
                      DO 190 l = 1,k
                          temp = temp + a(l,i)*a(l,j)
  190                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  200             CONTINUE
  210         CONTINUE
          ELSE
              DO 240 j = 1,n
                  DO 230 i = j,n
                      temp = zero
                      DO 220 l = 1,k
                          temp = temp + a(l,i)*a(l,j)
  220                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  230             CONTINUE
  240         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DSYRK .
*
      END

      SUBROUTINE dtrsm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*
*  -- Reference BLAS level3 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC max
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOUNIT,UPPER
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      parameter(one=1.0d+0,zero=0.0d+0)
*     ..
*
*     Test the input parameters.
*
      lside = lsame(side,'L')
      IF (lside) THEN
          nrowa = m
      ELSE
          nrowa = n
      END IF
      nounit = lsame(diag,'N')
      upper = lsame(uplo,'U')
*
      info = 0
      IF ((.NOT.lside) .AND. (.NOT.lsame(side,'R'))) THEN
          info = 1
      ELSE IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
          info = 2
      ELSE IF ((.NOT.lsame(transa,'N')) .AND.
     +         (.NOT.lsame(transa,'T')) .AND.
     +         (.NOT.lsame(transa,'C'))) THEN
          info = 3
      ELSE IF ((.NOT.lsame(diag,'U')) .AND. (.NOT.lsame(diag,'N'))) THEN
          info = 4
      ELSE IF (m.LT.0) THEN
          info = 5
      ELSE IF (n.LT.0) THEN
          info = 6
      ELSE IF (lda.LT.max(1,nrowa)) THEN
          info = 9
      ELSE IF (ldb.LT.max(1,m)) THEN
          info = 11
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DTRSM ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (m.EQ.0 .OR. n.EQ.0) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (alpha.EQ.zero) THEN
          DO 20 j = 1,n
              DO 10 i = 1,m
                  b(i,j) = zero
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
*
*     Start the operations.
*
      IF (lside) THEN
          IF (lsame(transa,'N')) THEN
*
*           Form  B := alpha*inv( A )*B.
*
              IF (upper) THEN
                  DO 60 j = 1,n
                      IF (alpha.NE.one) THEN
                          DO 30 i = 1,m
                              b(i,j) = alpha*b(i,j)
   30                     CONTINUE
                      END IF
                      DO 50 k = m,1,-1
                          IF (b(k,j).NE.zero) THEN
                              IF (nounit) b(k,j) = b(k,j)/a(k,k)
                              DO 40 i = 1,k - 1
                                  b(i,j) = b(i,j) - b(k,j)*a(i,k)
   40                         CONTINUE
                          END IF
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 100 j = 1,n
                      IF (alpha.NE.one) THEN
                          DO 70 i = 1,m
                              b(i,j) = alpha*b(i,j)
   70                     CONTINUE
                      END IF
                      DO 90 k = 1,m
                          IF (b(k,j).NE.zero) THEN
                              IF (nounit) b(k,j) = b(k,j)/a(k,k)
                              DO 80 i = k + 1,m
                                  b(i,j) = b(i,j) - b(k,j)*a(i,k)
   80                         CONTINUE
                          END IF
   90                 CONTINUE
  100             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*inv( A**T )*B.
*
              IF (upper) THEN
                  DO 130 j = 1,n
                      DO 120 i = 1,m
                          temp = alpha*b(i,j)
                          DO 110 k = 1,i - 1
                              temp = temp - a(k,i)*b(k,j)
  110                     CONTINUE
                          IF (nounit) temp = temp/a(i,i)
                          b(i,j) = temp
  120                 CONTINUE
  130             CONTINUE
              ELSE
                  DO 160 j = 1,n
                      DO 150 i = m,1,-1
                          temp = alpha*b(i,j)
                          DO 140 k = i + 1,m
                              temp = temp - a(k,i)*b(k,j)
  140                     CONTINUE
                          IF (nounit) temp = temp/a(i,i)
                          b(i,j) = temp
  150                 CONTINUE
  160             CONTINUE
              END IF
          END IF
      ELSE
          IF (lsame(transa,'N')) THEN
*
*           Form  B := alpha*B*inv( A ).
*
              IF (upper) THEN
                  DO 210 j = 1,n
                      IF (alpha.NE.one) THEN
                          DO 170 i = 1,m
                              b(i,j) = alpha*b(i,j)
  170                     CONTINUE
                      END IF
                      DO 190 k = 1,j - 1
                          IF (a(k,j).NE.zero) THEN
                              DO 180 i = 1,m
                                  b(i,j) = b(i,j) - a(k,j)*b(i,k)
  180                         CONTINUE
                          END IF
  190                 CONTINUE
                      IF (nounit) THEN
                          temp = one/a(j,j)
                          DO 200 i = 1,m
                              b(i,j) = temp*b(i,j)
  200                     CONTINUE
                      END IF
  210             CONTINUE
              ELSE
                  DO 260 j = n,1,-1
                      IF (alpha.NE.one) THEN
                          DO 220 i = 1,m
                              b(i,j) = alpha*b(i,j)
  220                     CONTINUE
                      END IF
                      DO 240 k = j + 1,n
                          IF (a(k,j).NE.zero) THEN
                              DO 230 i = 1,m
                                  b(i,j) = b(i,j) - a(k,j)*b(i,k)
  230                         CONTINUE
                          END IF
  240                 CONTINUE
                      IF (nounit) THEN
                          temp = one/a(j,j)
                          DO 250 i = 1,m
                              b(i,j) = temp*b(i,j)
  250                     CONTINUE
                      END IF
  260             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*B*inv( A**T ).
*
              IF (upper) THEN
                  DO 310 k = n,1,-1
                      IF (nounit) THEN
                          temp = one/a(k,k)
                          DO 270 i = 1,m
                              b(i,k) = temp*b(i,k)
  270                     CONTINUE
                      END IF
                      DO 290 j = 1,k - 1
                          IF (a(j,k).NE.zero) THEN
                              temp = a(j,k)
                              DO 280 i = 1,m
                                  b(i,j) = b(i,j) - temp*b(i,k)
  280                         CONTINUE
                          END IF
  290                 CONTINUE
                      IF (alpha.NE.one) THEN
                          DO 300 i = 1,m
                              b(i,k) = alpha*b(i,k)
  300                     CONTINUE
                      END IF
  310             CONTINUE
              ELSE
                  DO 360 k = 1,n
                      IF (nounit) THEN
                          temp = one/a(k,k)
                          DO 320 i = 1,m
                              b(i,k) = temp*b(i,k)
  320                     CONTINUE
                      END IF
                      DO 340 j = k + 1,n
                          IF (a(j,k).NE.zero) THEN
                              temp = a(j,k)
                              DO 330 i = 1,m
                                  b(i,j) = b(i,j) - temp*b(i,k)
  330                         CONTINUE
                          END IF
  340                 CONTINUE
                      IF (alpha.NE.one) THEN
                          DO 350 i = 1,m
                              b(i,k) = alpha*b(i,k)
  350                     CONTINUE
                      END IF
  360             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of DTRSM .
*
      END

      INTEGER FUNCTION ilaenv( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
*
*  -- LAPACK auxiliary routine (version 3.9.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2019
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    name, opts
      INTEGER            ispec, n1, n2, n3, n4
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            i, ic, iz, nb, nbmin, nx
      LOGICAL            cname, sname, twostage
      CHARACTER          c1*1, c2*2, c4*2, c3*3, subnam*16
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          char, ichar, int, min, real
*     ..
*     .. External Functions ..
      INTEGER            ieeeck, iparmq, iparam2stage
      EXTERNAL           ieeeck, iparmq, iparam2stage
*     ..
*     .. Executable Statements ..
*
      GO TO ( 10, 10, 10, 80, 90, 100, 110, 120,
     $        130, 140, 150, 160, 160, 160, 160, 160)ispec
*
*     Invalid value for ISPEC
*
      ilaenv = -1
      RETURN
*
   10 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      ilaenv = 1
      subnam = name
      ic = ichar( subnam( 1: 1 ) )
      iz = ichar( 'Z' )
      IF( iz.EQ.90 .OR. iz.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( ic.GE.97 .AND. ic.LE.122 ) THEN
            subnam( 1: 1 ) = char( ic-32 )
            DO 20 i = 2, 6
               ic = ichar( subnam( i: i ) )
               IF( ic.GE.97 .AND. ic.LE.122 )
     $            subnam( i: i ) = char( ic-32 )
   20       CONTINUE
         END IF
*
      ELSE IF( iz.EQ.233 .OR. iz.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR.
     $       ( ic.GE.145 .AND. ic.LE.153 ) .OR.
     $       ( ic.GE.162 .AND. ic.LE.169 ) ) THEN
            subnam( 1: 1 ) = char( ic+64 )
            DO 30 i = 2, 6
               ic = ichar( subnam( i: i ) )
               IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR.
     $             ( ic.GE.145 .AND. ic.LE.153 ) .OR.
     $             ( ic.GE.162 .AND. ic.LE.169 ) )subnam( i:
     $             i ) = char( ic+64 )
   30       CONTINUE
         END IF
*
      ELSE IF( iz.EQ.218 .OR. iz.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( ic.GE.225 .AND. ic.LE.250 ) THEN
            subnam( 1: 1 ) = char( ic-32 )
            DO 40 i = 2, 6
               ic = ichar( subnam( i: i ) )
               IF( ic.GE.225 .AND. ic.LE.250 )
     $            subnam( i: i ) = char( ic-32 )
   40       CONTINUE
         END IF
      END IF
*
      c1 = subnam( 1: 1 )
      sname = c1.EQ.'S' .OR. c1.EQ.'D'
      cname = c1.EQ.'C' .OR. c1.EQ.'Z'
      IF( .NOT.( cname .OR. sname ) )
     $   RETURN
      c2 = subnam( 2: 3 )
      c3 = subnam( 4: 6 )
      c4 = c3( 2: 3 )
      twostage = len( subnam ).GE.11
     $           .AND. subnam( 11: 11 ).EQ.'2'
*
      GO TO ( 50, 60, 70 )ispec
*
   50 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      nb = 1
*
      IF( subnam(2:6).EQ.'LAORH' ) THEN
*
*        This is for *LAORHR_GETRFNP routine
*
         IF( sname ) THEN
             nb = 32
         ELSE
             nb = 32
         END IF
      ELSE IF( c2.EQ.'GE' ) THEN
         IF( c3.EQ.'TRF' ) THEN
            IF( sname ) THEN
               nb = 64
            ELSE
               nb = 64
            END IF
         ELSE IF( c3.EQ.'QRF' .OR. c3.EQ.'RQF' .OR. c3.EQ.'LQF' .OR.
     $            c3.EQ.'QLF' ) THEN
            IF( sname ) THEN
               nb = 32
            ELSE
               nb = 32
            END IF
         ELSE IF( c3.EQ.'QR ') THEN
            IF( n3 .EQ. 1) THEN
               IF( sname ) THEN
*     M*N
                  IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                     nb = n1
                  ELSE
                     nb = 32768/n2
                  END IF
               ELSE
                  IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                     nb = n1
                  ELSE
                     nb = 32768/n2
                  END IF
               END IF
            ELSE
               IF( sname ) THEN
                  nb = 1
               ELSE
                  nb = 1
               END IF
            END IF
         ELSE IF( c3.EQ.'LQ ') THEN
            IF( n3 .EQ. 2) THEN
               IF( sname ) THEN
*     M*N
                  IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                     nb = n1
                  ELSE
                     nb = 32768/n2
                  END IF
               ELSE
                  IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                     nb = n1
                  ELSE
                     nb = 32768/n2
                  END IF
               END IF
            ELSE
               IF( sname ) THEN
                  nb = 1
               ELSE
                  nb = 1
               END IF
            END IF
         ELSE IF( c3.EQ.'HRD' ) THEN
            IF( sname ) THEN
               nb = 32
            ELSE
               nb = 32
            END IF
         ELSE IF( c3.EQ.'BRD' ) THEN
            IF( sname ) THEN
               nb = 32
            ELSE
               nb = 32
            END IF
         ELSE IF( c3.EQ.'TRI' ) THEN
            IF( sname ) THEN
               nb = 64
            ELSE
               nb = 64
            END IF
         END IF
      ELSE IF( c2.EQ.'PO' ) THEN
         IF( c3.EQ.'TRF' ) THEN
            IF( sname ) THEN
               nb = 64
            ELSE
               nb = 64
            END IF
         END IF
      ELSE IF( c2.EQ.'SY' ) THEN
         IF( c3.EQ.'TRF' ) THEN
            IF( sname ) THEN
               IF( twostage ) THEN
                  nb = 192
               ELSE
                  nb = 64
               END IF
            ELSE
               IF( twostage ) THEN
                  nb = 192
               ELSE
                  nb = 64
               END IF
            END IF
         ELSE IF( sname .AND. c3.EQ.'TRD' ) THEN
            nb = 32
         ELSE IF( sname .AND. c3.EQ.'GST' ) THEN
            nb = 64
         END IF
      ELSE IF( cname .AND. c2.EQ.'HE' ) THEN
         IF( c3.EQ.'TRF' ) THEN
            IF( twostage ) THEN
               nb = 192
            ELSE
               nb = 64
            END IF
         ELSE IF( c3.EQ.'TRD' ) THEN
            nb = 32
         ELSE IF( c3.EQ.'GST' ) THEN
            nb = 64
         END IF
      ELSE IF( sname .AND. c2.EQ.'OR' ) THEN
         IF( c3( 1: 1 ).EQ.'G' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.
     $          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )
     $           THEN
               nb = 32
            END IF
         ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.
     $          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )
     $           THEN
               nb = 32
            END IF
         END IF
      ELSE IF( cname .AND. c2.EQ.'UN' ) THEN
         IF( c3( 1: 1 ).EQ.'G' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.
     $          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )
     $           THEN
               nb = 32
            END IF
         ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.
     $          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )
     $           THEN
               nb = 32
            END IF
         END IF
      ELSE IF( c2.EQ.'GB' ) THEN
         IF( c3.EQ.'TRF' ) THEN
            IF( sname ) THEN
               IF( n4.LE.64 ) THEN
                  nb = 1
               ELSE
                  nb = 32
               END IF
            ELSE
               IF( n4.LE.64 ) THEN
                  nb = 1
               ELSE
                  nb = 32
               END IF
            END IF
         END IF
      ELSE IF( c2.EQ.'PB' ) THEN
         IF( c3.EQ.'TRF' ) THEN
            IF( sname ) THEN
               IF( n2.LE.64 ) THEN
                  nb = 1
               ELSE
                  nb = 32
               END IF
            ELSE
               IF( n2.LE.64 ) THEN
                  nb = 1
               ELSE
                  nb = 32
               END IF
            END IF
         END IF
      ELSE IF( c2.EQ.'TR' ) THEN
         IF( c3.EQ.'TRI' ) THEN
            IF( sname ) THEN
               nb = 64
            ELSE
               nb = 64
            END IF
         ELSE IF ( c3.EQ.'EVC' ) THEN
            IF( sname ) THEN
               nb = 64
            ELSE
               nb = 64
            END IF
         END IF
      ELSE IF( c2.EQ.'LA' ) THEN
         IF( c3.EQ.'UUM' ) THEN
            IF( sname ) THEN
               nb = 64
            ELSE
               nb = 64
            END IF
         END IF
      ELSE IF( sname .AND. c2.EQ.'ST' ) THEN
         IF( c3.EQ.'EBZ' ) THEN
            nb = 1
         END IF
      ELSE IF( c2.EQ.'GG' ) THEN
         nb = 32
         IF( c3.EQ.'HD3' ) THEN
            IF( sname ) THEN
               nb = 32
            ELSE
               nb = 32
            END IF
         END IF
      END IF
      ilaenv = nb
      RETURN
*
   60 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      nbmin = 2
      IF( c2.EQ.'GE' ) THEN
         IF( c3.EQ.'QRF' .OR. c3.EQ.'RQF' .OR. c3.EQ.'LQF' .OR. c3.EQ.
     $       'QLF' ) THEN
            IF( sname ) THEN
               nbmin = 2
            ELSE
               nbmin = 2
            END IF
         ELSE IF( c3.EQ.'HRD' ) THEN
            IF( sname ) THEN
               nbmin = 2
            ELSE
               nbmin = 2
            END IF
         ELSE IF( c3.EQ.'BRD' ) THEN
            IF( sname ) THEN
               nbmin = 2
            ELSE
               nbmin = 2
            END IF
         ELSE IF( c3.EQ.'TRI' ) THEN
            IF( sname ) THEN
               nbmin = 2
            ELSE
               nbmin = 2
            END IF
         END IF
      ELSE IF( c2.EQ.'SY' ) THEN
         IF( c3.EQ.'TRF' ) THEN
            IF( sname ) THEN
               nbmin = 8
            ELSE
               nbmin = 8
            END IF
         ELSE IF( sname .AND. c3.EQ.'TRD' ) THEN
            nbmin = 2
         END IF
      ELSE IF( cname .AND. c2.EQ.'HE' ) THEN
         IF( c3.EQ.'TRD' ) THEN
            nbmin = 2
         END IF
      ELSE IF( sname .AND. c2.EQ.'OR' ) THEN
         IF( c3( 1: 1 ).EQ.'G' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.
     $          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )
     $           THEN
               nbmin = 2
            END IF
         ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.
     $          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )
     $           THEN
               nbmin = 2
            END IF
         END IF
      ELSE IF( cname .AND. c2.EQ.'UN' ) THEN
         IF( c3( 1: 1 ).EQ.'G' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.
     $          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )
     $           THEN
               nbmin = 2
            END IF
         ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.
     $          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )
     $           THEN
               nbmin = 2
            END IF
         END IF
      ELSE IF( c2.EQ.'GG' ) THEN
         nbmin = 2
         IF( c3.EQ.'HD3' ) THEN
            nbmin = 2
         END IF
      END IF
      ilaenv = nbmin
      RETURN
*
   70 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      nx = 0
      IF( c2.EQ.'GE' ) THEN
         IF( c3.EQ.'QRF' .OR. c3.EQ.'RQF' .OR. c3.EQ.'LQF' .OR. c3.EQ.
     $       'QLF' ) THEN
            IF( sname ) THEN
               nx = 128
            ELSE
               nx = 128
            END IF
         ELSE IF( c3.EQ.'HRD' ) THEN
            IF( sname ) THEN
               nx = 128
            ELSE
               nx = 128
            END IF
         ELSE IF( c3.EQ.'BRD' ) THEN
            IF( sname ) THEN
               nx = 128
            ELSE
               nx = 128
            END IF
         END IF
      ELSE IF( c2.EQ.'SY' ) THEN
         IF( sname .AND. c3.EQ.'TRD' ) THEN
            nx = 32
         END IF
      ELSE IF( cname .AND. c2.EQ.'HE' ) THEN
         IF( c3.EQ.'TRD' ) THEN
            nx = 32
         END IF
      ELSE IF( sname .AND. c2.EQ.'OR' ) THEN
         IF( c3( 1: 1 ).EQ.'G' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.
     $          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )
     $           THEN
               nx = 128
            END IF
         END IF
      ELSE IF( cname .AND. c2.EQ.'UN' ) THEN
         IF( c3( 1: 1 ).EQ.'G' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.
     $          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )
     $           THEN
               nx = 128
            END IF
         END IF
      ELSE IF( c2.EQ.'GG' ) THEN
         nx = 128
         IF( c3.EQ.'HD3' ) THEN
            nx = 128
         END IF
      END IF
      ilaenv = nx
      RETURN
*
   80 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ilaenv = 6
      RETURN
*
   90 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ilaenv = 2
      RETURN
*
  100 CONTINUE
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ilaenv = int( real( min( n1, n2 ) )*1.6e0 )
      RETURN
*
  110 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ilaenv = 1
      RETURN
*
  120 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ilaenv = 50
      RETURN
*
  130 CONTINUE
*
*     ISPEC = 9:  maximum size of the subproblems at the bottom of the
*                 computation tree in the divide-and-conquer algorithm
*                 (used by xGELSD and xGESDD)
*
      ilaenv = 25
      RETURN
*
  140 CONTINUE
*
*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
*
*     ILAENV = 0
      ilaenv = 1
      IF( ilaenv.EQ.1 ) THEN
         ilaenv = ieeeck( 1, 0.0, 1.0 )
      END IF
      RETURN
*
  150 CONTINUE
*
*     ISPEC = 11: infinity arithmetic can be trusted not to trap
*
*     ILAENV = 0
      ilaenv = 1
      IF( ilaenv.EQ.1 ) THEN
         ilaenv = ieeeck( 0, 0.0, 1.0 )
      END IF
      RETURN
*
  160 CONTINUE
*
*     12 <= ISPEC <= 16: xHSEQR or related subroutines.
*
      ilaenv = iparmq( ispec, name, opts, n1, n2, n3, n4 )
      RETURN
*
*     End of ILAENV
*
      END

      INTEGER          FUNCTION ieeeck( ISPEC, ZERO, ONE )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            ispec
      REAL               one, zero
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      REAL               nan1, nan2, nan3, nan4, nan5, nan6, neginf,
     $                   negzro, newzro, posinf
*     ..
*     .. Executable Statements ..
      ieeeck = 1
*
      posinf = one / zero
      IF( posinf.LE.one ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      neginf = -one / zero
      IF( neginf.GE.zero ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      negzro = one / ( neginf+one )
      IF( negzro.NE.zero ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      neginf = one / negzro
      IF( neginf.GE.zero ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      newzro = negzro + zero
      IF( newzro.NE.zero ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      posinf = one / newzro
      IF( posinf.LE.one ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      neginf = neginf*posinf
      IF( neginf.GE.zero ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      posinf = posinf*posinf
      IF( posinf.LE.one ) THEN
         ieeeck = 0
         RETURN
      END IF
*
*
*
*
*     Return if we were only asked to check infinity arithmetic
*
      IF( ispec.EQ.0 )
     $   RETURN
*
      nan1 = posinf + neginf
*
      nan2 = posinf / neginf
*
      nan3 = posinf / posinf
*
      nan4 = posinf*zero
*
      nan5 = neginf*negzro
*
      nan6 = nan5*zero
*
      IF( nan1.EQ.nan1 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      IF( nan2.EQ.nan2 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      IF( nan3.EQ.nan3 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      IF( nan4.EQ.nan4 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      IF( nan5.EQ.nan5 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      IF( nan6.EQ.nan6 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      RETURN
      END

      INTEGER FUNCTION iparmq( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
*
*  -- LAPACK auxiliary routine (version 3.7.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2017
*
*     .. Scalar Arguments ..
      INTEGER            ihi, ilo, ispec, lwork, n
      CHARACTER          name*( * ), opts*( * )
*
*  ================================================================
*     .. Parameters ..
      INTEGER            inmin, inwin, inibl, ishfts, iacc22
      parameter( inmin = 12, inwin = 13, inibl = 14,
     $                   ishfts = 15, iacc22 = 16 )
      INTEGER            nmin, k22min, kacmin, nibble, knwswp
      parameter( nmin = 75, k22min = 14, kacmin = 14,
     $                   nibble = 14, knwswp = 500 )
      REAL               two
      parameter( two = 2.0 )
*     ..
*     .. Local Scalars ..
      INTEGER            nh, ns
      INTEGER            i, ic, iz
      CHARACTER          subnam*6
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          log, max, mod, nint, real
*     ..
*     .. Executable Statements ..
      IF( ( ispec.EQ.ishfts ) .OR. ( ispec.EQ.inwin ) .OR.
     $    ( ispec.EQ.iacc22 ) ) THEN
*
*        ==== Set the number simultaneous shifts ====
*
         nh = ihi - ilo + 1
         ns = 2
         IF( nh.GE.30 )
     $      ns = 4
         IF( nh.GE.60 )
     $      ns = 10
         IF( nh.GE.150 )
     $      ns = max( 10, nh / nint( log( real( nh ) ) / log( two ) ) )
         IF( nh.GE.590 )
     $      ns = 64
         IF( nh.GE.3000 )
     $      ns = 128
         IF( nh.GE.6000 )
     $      ns = 256
         ns = max( 2, ns-mod( ns, 2 ) )
      END IF
*
      IF( ispec.EQ.inmin ) THEN
*
*
*        ===== Matrices of order smaller than NMIN get sent
*        .     to xLAHQR, the classic double shift algorithm.
*        .     This must be at least 11. ====
*
         iparmq = nmin
*
      ELSE IF( ispec.EQ.inibl ) THEN
*
*        ==== INIBL: skip a multi-shift qr iteration and
*        .    whenever aggressive early deflation finds
*        .    at least (NIBBLE*(window size)/100) deflations. ====
*
         iparmq = nibble
*
      ELSE IF( ispec.EQ.ishfts ) THEN
*
*        ==== NSHFTS: The number of simultaneous shifts =====
*
         iparmq = ns
*
      ELSE IF( ispec.EQ.inwin ) THEN
*
*        ==== NW: deflation window size.  ====
*
         IF( nh.LE.knwswp ) THEN
            iparmq = ns
         ELSE
            iparmq = 3*ns / 2
         END IF
*
      ELSE IF( ispec.EQ.iacc22 ) THEN
*
*        ==== IACC22: Whether to accumulate reflections
*        .     before updating the far-from-diagonal elements
*        .     and whether to use 2-by-2 block structure while
*        .     doing it.  A small amount of work could be saved
*        .     by making this choice dependent also upon the
*        .     NH=IHI-ILO+1.
*
*
*        Convert NAME to upper case if the first character is lower case.
*
         iparmq = 0
         subnam = name
         ic = ichar( subnam( 1: 1 ) )
         iz = ichar( 'Z' )
         IF( iz.EQ.90 .OR. iz.EQ.122 ) THEN
*
*           ASCII character set
*
            IF( ic.GE.97 .AND. ic.LE.122 ) THEN
               subnam( 1: 1 ) = char( ic-32 )
               DO i = 2, 6
                  ic = ichar( subnam( i: i ) )
                  IF( ic.GE.97 .AND. ic.LE.122 )
     $               subnam( i: i ) = char( ic-32 )
               END DO
            END IF
*
         ELSE IF( iz.EQ.233 .OR. iz.EQ.169 ) THEN
*
*           EBCDIC character set
*
            IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR.
     $          ( ic.GE.145 .AND. ic.LE.153 ) .OR.
     $          ( ic.GE.162 .AND. ic.LE.169 ) ) THEN
               subnam( 1: 1 ) = char( ic+64 )
               DO i = 2, 6
                  ic = ichar( subnam( i: i ) )
                  IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR.
     $                ( ic.GE.145 .AND. ic.LE.153 ) .OR.
     $                ( ic.GE.162 .AND. ic.LE.169 ) )subnam( i:
     $                i ) = char( ic+64 )
               END DO
            END IF
*
         ELSE IF( iz.EQ.218 .OR. iz.EQ.250 ) THEN
*
*           Prime machines:  ASCII+128
*
            IF( ic.GE.225 .AND. ic.LE.250 ) THEN
               subnam( 1: 1 ) = char( ic-32 )
               DO i = 2, 6
                  ic = ichar( subnam( i: i ) )
                  IF( ic.GE.225 .AND. ic.LE.250 )
     $               subnam( i: i ) = char( ic-32 )
               END DO
            END IF
         END IF
*
         IF( subnam( 2:6 ).EQ.'GGHRD' .OR.
     $       subnam( 2:6 ).EQ.'GGHD3' ) THEN
            iparmq = 1
            IF( nh.GE.k22min )
     $         iparmq = 2
         ELSE IF ( subnam( 4:6 ).EQ.'EXC' ) THEN
            IF( nh.GE.kacmin )
     $         iparmq = 1
            IF( nh.GE.k22min )
     $         iparmq = 2
         ELSE IF ( subnam( 2:6 ).EQ.'HSEQR' .OR.
     $             subnam( 2:5 ).EQ.'LAQR' ) THEN
            IF( ns.GE.kacmin )
     $         iparmq = 1
            IF( ns.GE.k22min )
     $         iparmq = 2
         END IF
*
      ELSE
*        ===== invalid value of ispec =====
         iparmq = -1
*
      END IF
*
*     ==== End of IPARMQ ====
*
      END

C DPOTRS

      SUBROUTINE dpotrs( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      parameter( one = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           dtrsm, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      info = 0
      upper = lsame( uplo, 'U' )
      IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( nrhs.LT.0 ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      ELSE IF( ldb.LT.max( 1, n ) ) THEN
         info = -7
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DPOTRS', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 .OR. nrhs.EQ.0 )
     $   RETURN
*
      IF( upper ) THEN
*
*        Solve A*X = B where A = U**T *U.
*
*        Solve U**T *X = B, overwriting B with X.
*
         CALL dtrsm( 'Left', 'Upper', 'Transpose', 'Non-unit', n, nrhs,
     $               one, a, lda, b, ldb )
*
*        Solve U*X = B, overwriting B with X.
*
         CALL dtrsm( 'Left', 'Upper', 'No transpose', 'Non-unit', n,
     $               nrhs, one, a, lda, b, ldb )
      ELSE
*
*        Solve A*X = B where A = L*L**T.
*
*        Solve L*X = B, overwriting B with X.
*
         CALL dtrsm( 'Left', 'Lower', 'No transpose', 'Non-unit', n,
     $               nrhs, one, a, lda, b, ldb )
*
*        Solve L**T *X = B, overwriting B with X.
*
         CALL dtrsm( 'Left', 'Lower', 'Transpose', 'Non-unit', n, nrhs,
     $               one, a, lda, b, ldb )
      END IF
*
      RETURN
*
*     End of DPOTRS
*
      END

C DSWAP

      SUBROUTINE dswap(N,DX,INCX,DY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC mod
*     ..
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*       code for both increments equal to 1
*
*
*       clean-up loop
*
         m = mod(n,3)
         IF (m.NE.0) THEN
            DO i = 1,m
               dtemp = dx(i)
               dx(i) = dy(i)
               dy(i) = dtemp
            END DO
            IF (n.LT.3) RETURN
         END IF
         mp1 = m + 1
         DO i = mp1,n,3
            dtemp = dx(i)
            dx(i) = dy(i)
            dy(i) = dtemp
            dtemp = dx(i+1)
            dx(i+1) = dy(i+1)
            dy(i+1) = dtemp
            dtemp = dx(i+2)
            dx(i+2) = dy(i+2)
            dy(i+2) = dtemp
         END DO
      ELSE
*
*       code for unequal increments or equal increments not equal
*         to 1
*
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            dtemp = dx(ix)
            dx(ix) = dy(iy)
            dy(iy) = dtemp
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      RETURN
      END

C DSYR

      SUBROUTINE dsyr(UPLO,N,ALPHA,X,INCX,A,LDA)
*
*  -- Reference BLAS level2 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER INCX,LDA,N
      CHARACTER UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      parameter(zero=0.0d+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JX,KX
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC max
*     ..
*
*     Test the input parameters.
*
      info = 0
      IF (.NOT.lsame(uplo,'U') .AND. .NOT.lsame(uplo,'L')) THEN
          info = 1
      ELSE IF (n.LT.0) THEN
          info = 2
      ELSE IF (incx.EQ.0) THEN
          info = 5
      ELSE IF (lda.LT.max(1,n)) THEN
          info = 7
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DSYR  ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((n.EQ.0) .OR. (alpha.EQ.zero)) RETURN
*
*     Set the start point in X if the increment is not unity.
*
      IF (incx.LE.0) THEN
          kx = 1 - (n-1)*incx
      ELSE IF (incx.NE.1) THEN
          kx = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
      IF (lsame(uplo,'U')) THEN
*
*        Form  A  when A is stored in upper triangle.
*
          IF (incx.EQ.1) THEN
              DO 20 j = 1,n
                  IF (x(j).NE.zero) THEN
                      temp = alpha*x(j)
                      DO 10 i = 1,j
                          a(i,j) = a(i,j) + x(i)*temp
   10                 CONTINUE
                  END IF
   20         CONTINUE
          ELSE
              jx = kx
              DO 40 j = 1,n
                  IF (x(jx).NE.zero) THEN
                      temp = alpha*x(jx)
                      ix = kx
                      DO 30 i = 1,j
                          a(i,j) = a(i,j) + x(ix)*temp
                          ix = ix + incx
   30                 CONTINUE
                  END IF
                  jx = jx + incx
   40         CONTINUE
          END IF
      ELSE
*
*        Form  A  when A is stored in lower triangle.
*
          IF (incx.EQ.1) THEN
              DO 60 j = 1,n
                  IF (x(j).NE.zero) THEN
                      temp = alpha*x(j)
                      DO 50 i = j,n
                          a(i,j) = a(i,j) + x(i)*temp
   50                 CONTINUE
                  END IF
   60         CONTINUE
          ELSE
              jx = kx
              DO 80 j = 1,n
                  IF (x(jx).NE.zero) THEN
                      temp = alpha*x(jx)
                      ix = jx
                      DO 70 i = j,n
                          a(i,j) = a(i,j) + x(ix)*temp
                          ix = ix + incx
   70                 CONTINUE
                  END IF
                  jx = jx + incx
   80         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DSYR  .
*
      END
