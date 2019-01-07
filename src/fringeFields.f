      subroutine qqfringe(x,y,z,dy,bx,by,bz,fact1,fact2,qrad1,qrad2)
c*************************************************************************
c  qqfringe -                                                            *
c      calculates combined fringe fields of two adjacent quads, q1 & q2. *
c      the exit of q1 is located at (x, y, z) and the entrance of q2 is  *
c      dy downstream from the exit of q1.                                *
c      facti describes qi - bquad, bhex, boct, bdec, bddec               *
c      qradi radius of qi                                                *
c*************************************************************************
Cvardan      real fact1(5),fact2(5),b(3)
      real fact1,fact2,b(3)
      bx=0.
      by=0.
      bz=0.
c  exit of q1 field
      xdum=x
      ydum=y
      zdum=z
c      write(6,*)'fortran out qqfr fact=', fact1,'  grad=',qrad1
      call mpoles1(3,xdum,zdum,ydum,fact1,qrad1,b(1),b(3),b(2))
      bx=b(1)
      by=b(2)
      bz=b(3)
c      write(9,100)xdum,ydum,zdum,bx,by,bz
c 100  format(1x,'x,y,z=',3f10.3,' bi=',3f11.6)
c 101  format(1x,'      ',3f10.3,'    ',3f11.6)
c  entrance of q2 field
      ydum=dy-ydum
c      write(9,102)(fact2(i),i=1,5),qrad2
c 102  format(1x,'calling mpoles1  ',6f8.2)
      call mpoles1(1,xdum,zdum,ydum,fact2,qrad2,b(1),b(3),b(2))
      bx=bx+b(1)
      by=by-b(2)
      bz=bz+b(3)
c      write(9,101)xdum,ydum,zdum,bx,by,bz
      return
      end

      
      subroutine mpoles1(in,nx,ny,nz,ngradn,nrad,nbx,nby,nbz)
c stolen from raytrace jjl 10/8/86
c****                                                                   
c**** calculation of multipole(poles) field components                  
c****                                                                   
c****                                                                   
c****                                                                   
c**** 2 - quadrupole  (grad1)                                           
c**** 3 - hexapole    (grad2)                                           
c**** 4 - octapole    (grad3)                                           
c**** 5 - decapole    (grad4)                                           
c**** 6 - dodecapole  (grad5)                                           
c****                                                                   
c****                                                                   
      implicit real*8(a-h,o-z)                                          
cvardan      real nx,ny,nz,ngrad(5),nrad,nbx,nby,nbz
      real nx,ny,nz,nrad,ngrad(5),nbx,nby,nbz,ngradn
      common  /blckk91/  c0, c1, c2, c3, c4, c5                          
c      data c0, c1, c2, c3, c4, c5/0.1122d0,8.5d0,-1.4982d0,3.5882d0,
c     &     -2.1209d0,1.7230d0/
c  above are not quite right -JJL 9/4/98
      data c0, c1, c2, c3, c4, c5/0.1039d0,6.27108d0,-1.51247d0,
     &     3.59946d0,-2.1323d0,1.7230d0/
c      data ngrad / ngradn, 0, 0, 0, 0 /

c      write(6,*)'fortran out =', ngradn,'  rad=',nrad
c      write(6,*)in,',',nx,',',ny,',',nz,',',ngradn,',',nrad
      if (nrad.eq.0.)then
      write(6,*)' error in mpoles1,  nrad= 0.'
      call exit(0)
      endif
      ngrad(1) = ngradn
      ngrad(2) = 0
      ngrad(3) = 0
      ngrad(4) = 0
      ngrad(5) = 0

c     write(6,*)'mpole     nfrad=',ngrad,' nrad=',nrad
      grad1 = -ngrad(1)/nrad
      grad2 =  ngrad(2)/nrad**2
      grad3 = -ngrad(3)/nrad**3
      grad4 =  ngrad(4)/nrad**4
      grad5 = -ngrad(5)/nrad**5

      rad=nrad
      d = 2. * rad                                                      
      frh  = 1.d0
      fro  = 1.d0
      frd  = 1.d0
      frdd = 1.d0
      dh  = frh *d
      do  = fro *d
      dd  = frd *d
      ddd = frdd*d
      x = nx
      y = ny
      z = nz
      x2 = x*x                                                          
      x3 = x2*x                                                         
      x4 = x3*x                                                         
      x5 = x4*x                                                         
      x6 = x5*x
      x7 = x6*x
      y2 = y*y                                                          
      y3 = y2*y                                                         
      y4 = y3*y                                                         
      y5 = y4*y                                                         
      y6 = y5*y
      y7 = y6*y
c****
      s = z/d                                                           
      call bpls1( 2, d, s, re, g1, g2, g3, g4, g5, g6 )
c      write(6,*) d, s ,re, g1, g2, g3, g4, g5, g6
      b2x = grad1*( re*y - (g2/12.)*(3.*x2*y + y3) +                    
     1   (g4/384.)*(5.*x4*y + 6.*x2*y3 + y5 ) -                         
     2   (g6/23040.)*(7.*x6*y + 15.*x4*y3 + 9.*x2*y5 + y7)  )
      b2y = grad1*( re*x - (g2/12.)*(x3 + 3.*x*y2) +                    
     1   (g4/384.)*(x5 + 6.*x3*y2 + 5.*x*y4 ) -                         
     2   (g6/23040.)*(x7 + 9.*x5*y2 + 15.*x3*y4 + 7.*x*y6) )
      b2z = grad1*( g1*x*y - (g3/12.)*(x3*y + x*y3 ) +                  
     1   (g5/384.)*(x5*y +2.*x3*y3 + x*y5)  )
c****
c****
      ss = z/dh  + dsh
      call bpls1( 3, dh, ss, re, g1, g2, g3, g4, g5, g6 )
      b3x = grad2*( re*2.*x*y - (g2/48.)*(12.*x3*y + 4.*x*y3 ) )        
      b3y = grad2*( re*(x2-y2) - (g2/48.)*(3.*x4 + 6.*x2*y2 - 5.*y4 ) ) 
      b3z = grad2*( g1*(x2*y - y3/3.) - (g3/48.)*(3.*x4*y+2.*x2*y3-y5)) 
c****
c****
      ss = z/do  + dso
      call bpls1( 4, do, ss, re, g1, g2, g3, g4, g5, g6 )
      b4x = grad3*( re*(3.*x2*y - y3) - (g4/80.)*(20.*x4*y - 4.*y5 ) )  
      b4y = grad3*( re*(x3 - 3.*x*y2) - (g4/80.)*(4.*x5-20.*x*y4 ) )    
      b4z = grad3*g1*(x3*y - x*y3 )                                     
c****
c****
      ss = z/dd  + dsd
      call bpls1( 5, dd, ss, re, g1, g2, g3, g4, g5, g6 )
      b5x = grad4*re*(4.*x3*y - 4.*x*y3)                                
      b5y = grad4*re*(x4 - 6.*x2*y2 + y4 )                              
      b5z = grad4*g1*(x4*y - 2.*x2*y3 + y5/5. )                         
c****
c****
      ss = z/ddd + dsdd
      call bpls1( 6, ddd,ss, re, g1, g2, g3, g4, g5, g6 )
      b6x = grad5*re*(5.*x4*y - 10.*x2*y3 + y5 )                        
      b6y = grad5*re*(x5 - 10.*x3*y2 + 5.*x*y4 )                        
      b6z = 0.                                                          
c****
c****
      bx = b2x + b3x + b4x + b5x + b6x                                  
      by = b2y + b3y + b4y + b5y + b6y                                  
      bz = b2z + b3z + b4z + b5z + b6z                                  
      bt =   sqrt( bx*bx + by*by + bz*bz )                             
      nbx=bx
      nby=by
      nbz=bz
      return                                                            
      end                                                               


      subroutine bpls1 ( igp, d, s, re, g1, g2, g3, g4, g5, g6 )
c****
c****
c****
      implicit real*8 (a-h,o-z)
c****
c****
      common  /blckk91/  c0, c1, c2, c3, c4, c5                          
c****
c****
      s2 = s*s
      s3 = s2*s
      s4 = s2*s2
      s5 = s4*s
      cs = c0 + c1*s + c2*s2 + c3*s3 + c4*s4 + c5*s5                    
      cp1 =(c1 + 2.*c2*s + 3.*c3*s2 + 4.*c4*s3 + 5.*c5*s4) / d          
      cp2 = (2.*c2 + 6.*c3*s + 12.*c4*s2 + 20.*c5*s3  ) / (d*d)         
      cp3 = ( 6.*c3 + 24.*c4*s + 60.*c5*s2 ) / (d**3)                   
      cp4 = ( 24.*c4 + 120.*c5*s ) / (d**4)                             
c****
      cp5 = 120.*c5/(d**5)
c****
c****
c****
      if( abs(cs) .gt. 70. )  cs = sign(70.d0, cs )                   
      e = exp(cs)                                                      
c      write(6,*) 'cs=',cs
      re = 1./(1. + e)                                                  
      ere = e*re                                                        
      ere1= ere*re
      ere2= ere*ere1                                                    
      ere3= ere*ere2                                                    
      ere4= ere*ere3                                                    
c****
      ere5= ere*ere4
      ere6= ere*ere5
c****
c****
      cp12 = cp1*cp1                                                    
      cp13 = cp1*cp12                                                   
      cp14 = cp12*cp12                                                  
      cp22 = cp2*cp2                                                    
c****
      cp15 = cp12*cp13
      cp16 = cp13*cp13
      cp23 = cp2*cp22
      cp32 = cp3*cp3
c****
c****
      if( igp .eq. 6 ) return
      g1 = -cp1*ere1                                                    
c****
c****
      if( igp .eq. 5 ) return
      if( igp .eq. 4 ) go to 1
      g2 =-( cp2+cp12   )*ere1    + 2.*cp12 * ere2                      
      g3 =-(cp3 + 3.*cp1*cp2 + cp13  ) * ere1      +                    
     1   6.*(cp1*cp2 + cp13)*ere2 - 6.*cp13*ere3                        
c****
c****
      if( igp .eq. 3 ) return
1     g4 = -(cp4 + 4.*cp1*cp3 + 3.*cp22 + 6.*cp12*cp2 + cp14)*ere1  +   
     1   (8.*cp1*cp3 + 36.*cp12*cp2 + 6.*cp22 + 14.*cp14)*ere2    -     
     2   36.*(cp12*cp2 + cp14)*ere3       + 24.*cp14*ere4               
c****
c****
      if( igp .ne. 2 ) return
      g5 = (-cp5 - 5.*cp1*cp14 - 10.*cp2*cp3 - 10.*cp12*cp3 -
     1     15.*cp1*cp22 - 10.*cp13*cp2 - cp15)*ere1 +
     2     (10.*cp1*cp4 +20.*cp2*cp3 +60.*cp12*cp3 + 90.*cp1*cp22 +
     3     140.*cp13*cp2 +30.*cp15)*ere2 + (-60.*cp12*cp3 -
     4     90.*cp1*cp22 - 360.*cp13*cp2 - 150.*cp15)*ere3 +
     5     (240.*cp13*cp2 +240.*cp15)*ere4 + (-120.*cp15)*ere5
      g6 = (-6.*cp1*cp5 - 15.*cp2*cp4 - 15.*cp12*cp4 - 10.*cp32 -
     1     60.*cp1*cp2*cp3 - 20.*cp13*cp3 - 15.*cp23 - 45.*cp12*cp22 -
     2     15.*cp14*cp2 - cp16)*ere1 + (12.*cp1*cp5 + 30.*cp2*cp4 +
     3     90.*cp12*cp4 +20.*cp32 + 360.*cp1*cp2*cp3 +280.*cp13*cp3 +
     4     90.*cp23 + 630.*cp12*cp22 + 450.*cp14*cp2 + 62.*cp16)*ere2 +
     5     (-90.*cp12*cp4 - 360.*cp1*cp2*cp3 -720.*cp13*cp3 -90.*cp23 -
     6     1620.*cp12*cp22 -2250.*cp14*cp2 - 540.*cp16)*ere3 +
     7     (480.*cp13*cp3 + 1080.*cp12*cp22 + 3600.*cp14*cp2 +
     8     1560.*cp16)*ere4 + (-1800.*cp14*cp2 - 1800.*cp16)*ere5 +
     9     720.*cp16*ere6
c****
      return
      end
