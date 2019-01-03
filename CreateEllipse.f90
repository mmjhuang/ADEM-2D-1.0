     integer iseed, i, j, n
     real    rnd, r, rho, pi
     real data1(30,25000) 

     iseed = 425001

     r=5.e-4
     rho=2524
     din= r *3.0 * 2.
     n=1
     pi=2.*acos(0.)


do j=1,40
do i=1,55
     data1(1,n)=r/2.
     data1(2,n)=rho
     data1(4,n)= r*2.
     data1(3,n)=4*pi* data1(1,n)**2 * data1(4,n) *rho/3
     data1(5,n)=i*4.4*r / 2.  + R*3.
     data1(6,n)=din  + R*3. 
     data1(7,n)=0.
     data1(8,n)=0.

     data1(9,n)=0.

     call random_number(rnd)
     data1(10,n)=-1. + 2.*rnd
     call random_number(rnd)
     data1(11,n)=0. 
     data1(12,n)=0.
     data1(13,n)=0.

       do k=14,27
       data1(7,k)=0.
       enddo

     data1(28,n)=3. 
     n=n+1
end do
      din=din+3.0*r*2.
end do

! add baby particles
goto 1030   ! baby NOT added
do i=1,400*5
   if  ( abs( data1(1,i)/data1(4,i))>1.05 ) then
     data1(1,n)=r
     data1(2,n)=rho
     data1(3,n)=4*pi*r**3*rho/3
     data1(4,n)=0.
     data1(5,n)=data1(5,i)+ r* cos (data1(19,i)) 
     data1(6,n)=data1(6,i)+ r* sin (data1(19,i)) 
     data1(7,n)=0.
     data1(8,n)=0.
     data1(9,n)=0.
     data1(10,n)=data1(10,i)
     call random_number(rnd)
     data1(11,n)=0. 
     data1(12,n)=0.
     data1(13,n)=0.
       do k=14,27
       data1(7,k)=0.
       enddo
     data1(18,n)=i
     data1(28,n)=3. 
     data1(29,n)=0.
     data1(30,n)=0.
     n=n+1
   endif
end do
1030 continue
! end of adding babies 

print *, "REAL" , N-1

goto 1040 
! create the ghost particle in the left side
do i=1,1000
     data1(1,n)=r
     data1(2,n)=rho
     data1(3,n)=4*pi*r**3*rho/3
     data1(4,n)=0.
     data1(5,n)=0.5*R
     data1(6,n)=i*2.4*r
       do k=7,27
       data1(7,k)=0.
       end do
     data1(28,n)=-3. 
     data1(29,n)=0.
     data1(30,n)=0.
     n=n+1
end do
!-------

!create the ghost particle in the right side
do i=1,1000
     data1(1,n)=r
     data1(2,n)=rho
     data1(3,n)=4.*pi*r**3.*rho/3.
     data1(4,n)=0.
     data1(5,n)=7.5E-3 + 0.5*R
     data1(6,n)=i*2.4*r
       do k=7,27
       data1(7,k)=0.
       enddo   
     data1(28,n)=-3. 
     n=n+1
end do
!-------
     top= data1(6,n-1)


! -----create the ghost particle in the bottom
goto 1033
do i=1,6
     data1(1,n)=r
     data1(2,n)=rho
     data1(3,n)=4.*pi*r**3.*rho/3.
     data1(4,n)=0.
     data1(5,n)=i*2.3*r 
     data1(6,n)=0.5*R
       do k=7,27
       data1(7,k)=0.
       enddo        
     data1(28,n)=-3. 
     data1(29,n)=0.
     data1(30,n)=0.
     n=n+1
end do
1033 continue
!--------
   
!goto 1040

!-----create the ghost particle in the top
do i=1,6
     data1(1,n)=r
     data1(2,n)=rho
     data1(3,n)=4*pi*r**3*rho/3.
     data1(4,n)=0.
     data1(5,n)=i*2.3*r 
     data1(6,n)=top
       do k=7,27
       data1(7,k)=0.
       enddo     
     data1(28,n)=-3. 
     data1(29,n)=0.
     data1(30,n)=0.
     n=n+1
end do
!--------
1040 continue
! Codes to delelte meaningless info. can be added to here


open(20, file='file000000.dat')
 write (20, 110) data1
 110 format (30E16.8)
close(20)

PRINT *, "total", n-1	 
end
