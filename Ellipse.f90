!These works are writtten in Hong Kong Polytechnic University from Oct.2008-Oct.2009
!In order to make the flow chart clearly and add the fluid phase, it is rewrittn in NTNU
!Author: Yrjo^:) Jun HUANG	
!version Swendish Meatball (V.3.22)    
!11st Nov 2013. After coming back from Chengdu
!30 June, 2014, Canada
!10 Sept.2014, Montreal, Canada, V4.0 for ellipse particle
!20 Nov.2014, Shanghai, "timepara" is added. simulations can be restared from the previous data saved
!   subroutine delepart is divied into two subroutines, delepart and numfull
!01 Jan 2015, Shanghai, added subroutine "ghostwall", one wall only. The rest walls would be added by DONG SU
!10 Feb 2015, Change the nb from table index to chain index, 40% faster than before
!20 Sept.2015, Y.J aaded the rest four walls using subroutine ghostwall. 
!26 Sept.2015, Y.J added subroutine coat to avoid the particle overlap with neighboring particles. 
!in the previous version 20150210beta, there is a bug:
!if the number of particles created by particlecreate is fewer or more, sometime shows an error. 
!this version solved this problem, in subrotine: "for i =1 to 3000" 3000 was changed to numfull.

 
	module face		 
	integer, allocatable:: nb(:,:), nbx(:) , nby(:)	
	integer, allocatable:: mv(:,:,:)  	   
        end module	

!	=================== 
	program ApenDEM
!	=================== 
!The main code	
!input:   nump--maximum number of particles allowed; 
!         numcol--number of columns in array stma;
!         stma--main array;  xlower,xupper, ylower, yupper--compuational domain;  
!         tp0--initial time; tpfinal--finial time; rhop density;
!         tpout--output time gap;  fricp--surface roughness;
!         young--Young's modulus;  pois-- poision's ratio;
!         g0--shear modulus;  tau-- reaxlatiion time;
!         e--coefficient of resitution;  ginf--long term shear modulus;
!         youngr--relative young's modulus; 
!internal:istep--index for time steps between two screen outputs
!         iip--index for time steps between two files saved to harddisk 
!         kedm--time step of EDM / time step of TDM
!         iiv--time step for free particle / time step of TDM
!         numfull--number of rows in arrayment stma is not zero
!         timepara--parameters about time, including tp0, tpout and dtp0
!         nw -- particles close to each wall
!               (eg. nw(5000,2) means for 2 walls & each wall has 5000 neighboring particles)

	use face		
	implicit double precision (a-h, o-z)
	parameter (nump=25000)                  
	parameter (numcol=30)                 
	dimension stma(numcol, nump)   		
	dimension lsg(2,nump)	  
        dimension pai(numcol), paj(numcol), pak(numcol), pal(numcol)	
	dimension stma2(numcol, nump)
        dimension stma3(numcol, nump)
        dimension timepara(5)
        dimension nw(15000,4)
        dimension walld(7,4)
	allocate (nb(nump*50, 2)) ! allocate (nb(100,nump))
	allocate (nbx(nump))
	allocate (nby(nump))
   
       
	xlower=0.d0			
	xupper=1.e-1
	ylower=0.d0
	yupper=1.5      
	tpfinal=2.e-1	 
	tpout=1.e-3                                              	
	rhop= 2.5d3      
	fricp=2.d-1   		
	young=6.d10	                            
	pois=0.24d0		 
	g0=2.5d10        
	tau=1d-5		 
	e=0.01d0	 	
	ginf=0.25*g0      	
	pi=3.1415927    	
	youngr=young/(1.d0-pois*pois)/2	   	
	istep=1    ! index for screen 
	 
	open (10, file='file000030.dat')
	read (10,110) stma2
 110    format(30e16.8)
	close(10)
        stma=stma2

	open (13, file='time000030.dat')
  	read (13,113) timepara
 113    format(5e18.10)
	close(13) 
         
        tp0=timepara(3)   
        print *,'numfull 001=', timepara(5)                                                                                          
                                                             
        if (tp0<1.d-20)  then
         iip=1
          call delepart (stma, numcol, nump, numfull)              
          call addbaby (stma, numcol, nump, numfull)
          call addcoat (stma, numcol, nump, numfull, ratio)
          call getnumfull (stma, numcol, nump, numfull)   
          call reset (stma, numcol, nump, tp0)                          
          call intstat(nump, stma, numcol)   
          call field (stma, numcol, nump, tp0) 
          call dtpint (stma, youngr, vmax, numcol, nump, radmax, dtp0, pmassminr, radmin, iedm, dtfree)	                           
       	call searchsize (xlower, xupper, ylower, yupper, radmax,  mvx, mvy, cellp)                                    
	call search (stma, numcol, nump, xlower, xupper, ylower, yupper, cellp, lsg)     
        call BoundaryCellFix (stma, xlower, xupper, ylower, yupper, radmax,  mvx, mvy, cellp, nw, walld, numcol, nump, tp0)       
          iip=1
          iiv=1
        endif 
          
        if (tp0>=1.d-20) then 
        stma2=stma  
        call resetI (stma, numcol, nump, tp0) 
        call addcoat (stma, numcol, nump, numfull, ratio)
        iip=floor(timepara(4)+0.49)+1
        tpout=timepara(2) 
        print *, 'iip=', iip
        numfull = floor (timepara(5)+0.49)
        call dtpint (stma, youngr, vmax, numcol, nump, radmax, dtp0, pmassminr, radmin, iedm, dtfree)	      
   	call searchsize (xlower, xupper, ylower, yupper, radmax,  mvx, mvy, cellp)   
        call search (stma, numcol, nump, xlower, xupper, ylower, yupper, cellp, lsg)
        CALL BoundaryCellFix (stma, xlower, xupper, ylower, yupper, radmax,  mvx, mvy, cellp, nw, walld, numcol, nump, tp0)
        iiv=10000000  
        endif     

100	do while (tp0<tpfinal)  
	 if (iiv>=400) then	                        ! change from 500 to 400 to 20
	 call search (stma, numcol, nump, xlower, xupper, ylower, yupper, cellp, lsg)	                     
         call BoundaryCellFix (stma, xlower, xupper, ylower, yupper, radmax,  mvx, mvy, cellp, nw, walld, numcol, nump, tp0)                    
	 call dtp1 (stma, youngr, vmax, numcol, nump, dtp0, pmassminr, radmin, iedm, kedm, dtedm, tp0, tpminus1, dtfree)           
	   iiv=0
	   else
	   iiv=iiv+1
        endif
              
	call reset (stma, numcol, nump, tp0)      
        call field (stma, numcol, nump, tp0)
        call binary (stma, numcol, nump, pai, paj, pak, pal, alpha, e, youngr, dist, pi, fricp, dtp0, kedm, dtedm, deltamax, ratio)                     
	call ghostwall (stma, numcol, nump, pai, paj, pak, pal,  alpha,  e, youngr, dist, pi, fricp, dtp0, kedm, dtedm, nw, walld)                  	   
        call velocityPE(stma, numcol, nump, dtp0, kedm, dtedm, dtfree, iiv)    
	
            tp0=tp0+dtp0
    
 	call outputp (stma, nump, numcol, dtp0, tpout, istep, tp0, vmax, iip, dudt, pressF, partF, FricF, numfull, deltamax, ratio, stma2)  
           if (kedm==50) then
           kedm=1
           else
           kedm=kedm+1
           endif	 	
	enddo    
 
        print *, 'finish' 
	  	    
	end	   		  
!
!
!!=====================================
    subroutine delepart (stma, numcol, nump, numfull)
!=====================================  
!     make the mass of all ghost partile infinite
!     move all empty rows and rows for baby particles to the end of array stma
!     input:stma; 
!     output:stma

      use face
      implicit double precision (a-h, o-z)        
      dimension stma(numcol, nump)

      do i=1, nump
       if (int(stma(28,i)+0.43)==0  .or. int(stma(18,i)+0.43)>=0.98) then 
         do j=i+1, nump
	   do k=1, numcol
             stma(k,i)=stma(k,j)
           enddo
         enddo 
            do k=1, numcol
             stma(k,nump)=0.
            enddo   
       endif
       if (int(stma(28,i)+0.43)==-3)  stma(3,i)=1.d6
      enddo 

      goto 5030

      do i=1, nump        
        if   (int(stma(28,i)+0.43)==0) then
        numfull = i-1 
        goto 5030 
        endif
      enddo 
5030  continue        
          
      end
!
!
!!=====================================
    subroutine GetNumfull (stma, numcol, nump, numfull)
!=====================================  
!     If the code is not start from the inital state (baby particles were added)
!     using this subroutine to instead delpartII
!     Empty rows should be removed from the array stma
!     in the subroutine delepart already.

      use face
      implicit double precision (a-h, o-z)        
      dimension stma(numcol, nump)

      do i=1, nump        
        if   (int(stma(28,i)+0.43)==0) then
        numfull = i-1 
        goto 5030 
        endif
      enddo 
5030  continue    
        
      end


!!=====================================
    subroutine addbaby (stma, numcol, nump, numfull)
!=====================================  
!     make the mass of all ghost partile infinite
!     move all empty rows to the end of array stma
!     input:stma; 
!     output:stma

      use face
      implicit double precision (a-h, o-z)        
      dimension stma(numcol, nump)

      real :: a(3)=(/0.9,0.75,0.6/)
      real :: b(3)=(/0.98,1.88,2.60/)

      do i=1, nump        
        if   (int(stma(18,i)+0.43)/=0 .or. int(stma(28,i)+0.43)==0 ) then
        numfull = i-1 
        goto 5040 
        endif
      enddo 
5040  continue  

      numfull2=numfull

        do j=1, numfull  
       if (int(stma(28,j)+0.43)==3  .and.  abs( stma(4,j)/stma(1,j))>1.05) then 
       do 5050 m=1,3,1
        do i=-1,1,2 
     	stma(1,numfull2+1) = a(m)*stma(1,j) 
     	stma(2,numfull2+1) = stma(2,j)     	
	stma(3,numfull2+1) = stma(3,j)             
     	stma(4,numfull2+1) = stma(4,j)
        stma(5,numfull2+1) = stma(5,j)  +  i* b(m)* stma(1,j)   * cos(stma(19,j)) 
        stma(6,numfull2+1) = stma(6,j)  +  i* b(m)* stma(1,j)   * sin(stma(19,j)) 
     	stma(7,numfull2+1) = stma(7,j)
     	stma(8,numfull2+1) = stma(8,j)
     	stma(9,numfull2+1) = stma(9,j)
     	stma(10,numfull2+1) = stma(10,j)
     	stma(11,numfull2+1) = stma(11,j)
     	stma(12,numfull2+1) = stma(12,j)         
     	stma(13,numfull2+1) = i* b(m) * stma(1,j)
     	 do k=14,27
       	 stma(k,numfull2+1) = stma(k,j)
       	 enddo
     	stma(18,numfull2+1)=j
     	stma(28,numfull2+1)=3     
    	numfull2=numfull2+1
       enddo    	
5050       enddo       
!goto 12312   
12312 continue
       endif       
      enddo 
        
      end


!!=====================================
    subroutine addcoat (stma, numcol, nump, numfull, ratio)
!=====================================  
!     the particle diameters are times by a ratio, ratio
!     ratio -- is this ratio

      use face
      implicit double precision (a-h, o-z)        
      dimension stma(numcol, nump)

      ratio =1.005
      do i=1, nump        
       if (int(stma(28,i)+0.43)==0) goto 1488
       stma(1,i)=stma(1,i) * ratio
      enddo    
1488  continue    
      end


!=====================================
      subroutine reset (stma, numcol, nump, tp0)
!=====================================
!     before each time step, make the accleration zero
!     input:stma
!     output:stma

      use face
      implicit double precision (a-h, o-z)  
      dimension stma(numcol, nump)
           
	do 1405 i=1,nump	 
          if (int(stma(28,i)+0.43)==0) goto 1407
	stma(15,i)=0. 	   
        stma(16,i)=0.
	stma(17,i)=0.
	stma(25,i)=0. 	   
        stma(26,i)=0. 
	stma(27,i)=0.
                           
          if (int(stma(18,i)+0.43)>0) then
                  j=int(stma(18,i)+0.43)
        	  stma(5,i) = stma(5,j)  +   stma(13,i) * cos(stma(19,j))
                  stma(6,i) = stma(6,j)  +   stma(13,i) * sin(stma(19,j)) 
        	  stma(10,i) = stma(10,j) +  stma(13,i) * (-sin(stma(19,j))) * stma(22,j) ! * 0.00001
                  stma(11,i) = stma(11,j) +  stma(13,i) * cos(stma(19,j)) * stma(22,j)   ! * 0.00001
                  stma(19,i) = stma(19,j) 
                  stma(22,i) = stma(22,j)               
          endif

        if ((5007<=i) .and. (i<=5012) .and. (tp0<0.2) .and. stma(6,i)>0.32)    stma(6,i)=1.2-10. * 0.8* tp0        
   
1405	enddo
1407    continue
    
	end	
!
!

!=====================================
      subroutine resetI (stma, numcol, nump, tp0)
!=====================================
!     before each time step, make the accleration zero
!     input:stma
!     output:stma

      use face
      implicit double precision (a-h, o-z)  
      dimension stma(numcol, nump)
           
	do  i=1,nump 
        stma(1,i)=stma(1,i) - 1.d-6
        enddo

    
	end	
!
!=====================================
	subroutine field (stma, numcol, nump, tp0)
!=====================================
!        add accelration due to body force to each particle. 
!        Here, only gravitation is added
!input:  stma
!output: stma

       use face
       implicit double precision (a-h, o-z)  
       dimension stma(numcol, nump)
          
!        if (tp0<0.011) then
!        gravity=-(0.11-tp0)*10.*9.81 -9.81*0.001    
!        else
!        gravity=-9.81*0.001
!        endif
        gravity=-9.81 *12.

	do 1605 i=1,nump	 
        if (int(stma(28,i)+0.43)==0)  goto 1607	
        if (int(stma(28,i)+0.43)==3)  stma(16,i)=stma(16,i) + gravity
                                                       if (stma(6,i)>0.026) gravity=-1000.
                                                       !if (stma(6,i)>0.026)  stma(11,i)=-20.
1605	enddo
1607    continue   

	end	
!
!
!=====================================
	subroutine search(stma, numcol, nump, xlower, xupper, ylower, yupper, cellp, lsg)
!=====================================
!       search potential neighbors for each particle
! input:   stma;
!          cellp--size of background grid;
!          xlower,xupper, ylower, yupper--compuational domain;  
! output:  mv--list of all particles in each background grid;
!          nbx--number of neighbor for each particle; updated after search
!          nby--same as nbx, be updated in each time step
!          nb--list of number to each particle;
! internal:lsg--the location of all particles in the background
!          mvxv--the x-location for a particle
!          mvyv--the y-location for a particle    
!          iib--the first element in a row of nb     
!          ix, iy -- index for the neighboring cells for HA cell search, this idea is from SU Dong
!                    !! the orginial codes for HA-cell is from Dong Su, Y.J polished it.  

	use face		
	implicit double precision (a-h, o-z)   
	dimension lsg(2,nump) 
	dimension stma(numcol, nump)  
	
        nb=0  
	mv=0	
        nbx=0
        inb = size(nb,1)   ! size of nb
                                   
        do 1305 i=1,nump	
        if (int(stma(28,i)+0.43)==0) goto 1303		
                                 			 
	mvxv=int(stma(5,i)/cellp)+3
	mvyv=int(stma(6,i)/cellp)+3	 
			                          	                        
	lsg(1,i)=mvxv
	lsg(2,i)=mvyv   					
	mv(mvxv,mvyv,1)=mv(mvxv,mvyv,1)+1 	                	
	mv(mvxv,mvyv,(mv(mvxv,mvyv,1)+1))=i
                                 
1305	enddo	
1303    continue 
                          
	do 1320 i=1,nump
        if (int(stma(28,i)+0.43)==0) goto 1323		   	                         
   	  
        iy=0
        do 1325 ix=-1,0	
        if (stma(28, i)==0) goto 1330
        if  (mv((lsg(1,i)+ix), (lsg(2,i)+iy),1)==0) goto  1330
	 do 1328 jjp=2, mv((lsg(1,i)+ix),(lsg(2,i)+iy),1)+1	
        
          j=mv(lsg(1,i)+ix, lsg(2,i)+iy, jjp)
      
         if (int(stma(18,i)+0.4)==j .or.  int(stma(18,j)+0.4)==i .or. i==j)  goto 1332  
         if (int(stma(18,i)+0.4)==int(stma(18,j)+0.4) )  goto 1332  
         if (int(stma(28,i)+0.4)<0  .and. int(stma(28,j)+0.4)<0)   goto 1332                   
                dx=stma(5,i)-stma(5,j)
                dy=stma(6,i)-stma(6,j)
                dist2=dx**2+dy**2 
         if (dist2<cellp**2) then
              
               if(ix==0) then 
                  if (j<i) then    
                   nb(inb,1)=nb(inb,1)+1
                       if  (nb(inb,1) > inb-10)  then
                       print *, 'size of nb is small', nb(inb,1)
                       stop
                       endif
                   nb(nb(inb,1),1)=i
                   nb(nb(inb,1),2)=j  
                   nbx(j)=nbx(j) +1
                   nbx(i)=nbx(i) +1    
                  endif
               endif

               if(ix==-1) then
                  if (j/=i) then                  
                  nb(inb,1)=nb(inb,1)+1
                     if  (nb(inb,1) > inb-10)  then
                       print *, 'size of nb is small',  nb(inb,1)
                       stop
                     endif
                   nb(nb(inb,1),1)=i
                   nb(nb(inb,1),2)=j
                   nbx(j)= nbx(j) + 1
                   nbx(i)= nbx(i) + 1 
                  endif
                endif
          endif
1332   continue
1328   enddo
1330   continue 
1325   enddo 

        iy=1
        do 1326 ix=-1,1	
        if (stma(28, i)==0) goto 1331
        if  (mv((lsg(1,i)+ix), (lsg(2,i)+iy),1)==0) goto  1331
	 do 1329 jjp=2, mv((lsg(1,i)+ix),(lsg(2,i)+iy),1) + 1	
         j=mv(lsg(1,i)+ix, lsg(2,i)+iy, jjp)
           
         if (int(stma(18,i)+0.4)==j .or.  int(stma(18,j)+0.4)==i .or. i==j)  goto 1333 
         if (int(stma(18,i)+0.4)==int(stma(18,j)+0.4))  goto 1333 
         if (int(stma(28,i)+0.4)<0  .and. int(stma(28,j)+0.4)<0)   goto 1333                       
                dx=stma(5,i)-stma(5,j)
                dy=stma(6,i)-stma(6,j)
                dist2=dx**2+dy**2 
         if (dist2<cellp**2) then
                  if (j/=i) then                  
                  nb(inb,1)=nb(inb,1)+1
                     if  (nb(inb,1) > inb-10)  then
                       print *, 'size of nb is small', nb(inb,1)
                       stop
                     endif
                   nb(nb(inb,1),1)=i
                   nb(nb(inb,1),2)=j
                   nbx(j)=nbx(j) + 1
                   nbx(i)=nbx(i) + 1 
                  endif
          endif
1333   continue
1329   enddo
1331   continue 
1326   enddo                     	  
1323   continue  
               
1320   enddo 
      nby=nbx            
       end
!
!
!=====================================
	subroutine dtpint (stma, youngr, vmax, numcol, nump, radmax, dtp0, pmassminr, radmin, iedm, dtfree)		 
!=====================================
!       calculate the time step
! input:   stma
! output:  dtp0--time step for TDM algorithm
!          dtedm--time step for EDM algorithm
!          dtfree--time step for free flight particle
!          pmassmin--the minmum mass of all particle 
!          radmin--the radius of the smallest particle 
!          radmax--the radius of the biggest particle!         
!          pmassminr--ralative mass of two particles with mass of pmassmin
! internal:vmax--the velocity of the fastest particle !          
!          vsp--speed of each particle

        use face
	implicit double precision (a-h, o-z)   
	dimension stma(numcol, nump)  

	vmax=0.d0
	pmassmin=1.d2	   ! a big enough number
	radmin=stma(1,1) *100.		   ! a big enough number
	radmax=stma(1,1) * 0.001		   ! a small enough number
	do 1110 i=1, nump
          if (int(stma(28, i)+0.4)<=0) goto 1112
	  vsp=dsqrt(stma(10,i)*stma(10,i)+stma(11,i)*stma(11,i)) 
                                          !vsp: velocity of single particle
	  if (vsp>vmax) vmax=vsp 
                                            
	  if (stma(3,i)<pmassmin) pmassmin=stma(3,i)	 	                                     
	  if (stma(1,i)<radmin) radmin=stma(1,i)  
	  if (stma(1,i)>radmax) radmax=stma(1,i)

1110    enddo	
1112    continue		
		
	if (vmax<1.d0) vmax=1.d0	
	pmassminr=pmassmin/2.d0    
	
        dtp0=0.01 *2.87d0*(pmassminr*pmassminr/4/(0.5d0*radmin*youngr *youngr*vmax))**0.2        /2.
	if (dtp0>5.d-3*radmin/vmax) dtp0=5.d-3*radmin/vmax
        dtedm=50*dtp0
        dtfree=500.*dtp0
	end
!
!
!
!=====================================
	subroutine intstat(nump, stma, numcol)
!=====================================
! define/change the initial state of each particle in array stma
! input: stma
! putput: stma

	use face
        implicit double precision (a-h, o-z)   
	dimension stma(numcol, nump) 
			     				   
	do 1010 i=1, nump	
          if (int(stma(28,i)+0.43)==0)  goto 1020	
	stma(10,i)=0.0d0
        stma(19,i)=3.1415926/2.
	  if (int(stma(28,i)+0.43)==3)     stma(11,i)=stma(11,i) -10. ! initial velocity
          if (int(stma(28,i)+0.43)==-3)    stma(11,i)= 0. 
   
1010  enddo 
1020  continue        
	end
!
!
!
!=====================================
	subroutine dtp1 (stma, youngr, vmax, numcol, nump, dtp0, pmassminr, radmin, iedm, kedm, dtedm, tp0, tpminus1, dtfree)
!=====================================
! This subroutine is the same as subroutine dtpint, but simpiler
! cause  radmin, radmax, pmassminr and pmassmin are calcualted in dtpinf
! input:   stma
! output:  dtp0--time step for TDM algorithm
!          dtedm--time step for EDM algorithm
!          dtfree--time step for free flight particle        
! internal:vmax--the velocity of the fastest particle          
!          vsp--speed of each particle

        use face
	implicit double precision (a-h, o-z)   
	dimension stma(numcol, nump)  
                      
	vmax=1.d-1	  
	do 1910 i=1, nump
        if (int(stma(28,i)+0.43)==0)  goto 1918                                                        	
	vsp=dsqrt(stma(10,i)*stma(10,i)+stma(11,i)*stma(11,i)) 
	if (vsp>vmax) vmax=vsp 
1910    enddo	
1918    continue
	 	
        dtp0=(2.d-2)*2.87d0*(pmassminr*pmassminr/4/(0.5d0*radmin*youngr*youngr*vmax))**0.2   !contact time
	if (dtp0>.1d-2*radmin/vmax) dtp0=.1d-2*radmin/vmax  
        dtedm=dtp0*50. 
        dtfree=500.*dtp0

	end
!
!
!=========================
	subroutine binary (stma, numcol, nump, pai, paj, pak, pal,  alpha, e, youngr, dist, pi, fricp, dtp0, kedm, dtedm, deltamax, ratio)
!========================
!  This subroutine calculates the accelration (call TDM)
!  or velocity (call EDM) due to contact.
!  input: nb--list of neighbor
!         nbx--number of neigbors for each particle  
!         e, youngr, fricp, kedm, pi -- defined in the main program
!         dtedm-- time step for EDM 
!  outpur: stma
! internal: pai,paj--information for two colliding particle, to save memory        
!           alpha--angle of two particles to the cooridnates
!           dist--distance between two particles
!           coat-- thick of imagine coat
!    contactX & contactY -- cooridantes of the contact point

        use face
	implicit double precision (a-h, o-z)   
	dimension stma(numcol, nump)  	
	dimension pai(numcol), paj(numcol), pak(numcol), pal(numcol)
                
        nbx=nby     
        deltamax=0.d-4     
        inb = size(nb,1) 
  	   	                          
	
        do 1530 k0p=1, nb(inb,1)   
        i=nb(k0p,1)
        j=nb(k0p,2)  
     
    	  disqr=((stma(5,i)-stma(5,j))**2 + (stma(6,i)-stma(6,j))**2 )
	  dist=dsqrt(disqr) 
    
           if  ((nbx(i)>=1)  .or. (nbx(j)>=1) )  then              ! use tdm,   
           if (dist<(stma(1,i)+stma(1,j)) ) then
                kjp=int(stma(18,i) + 0.43)
                ljp=int(stma(18,j) + 0.43)   
           logi=0
           logj=0
           if (kjp/=0) logi=1
           if (ljp/=0) logj=1  
	   do 1511 jp=1,14    
	   pai(jp)=stma(jp,i)
	   paj(jp)=stma(jp,j)
              if (logi==1)  pak(jp)=stma(jp, kjp) 
              if (logj==1)  pal(jp)=stma(jp, ljp)               
1511	   end do

	   do 1512 jp=18,24
	   pai(jp)=stma(jp,i)
	   paj(jp)=stma(jp,j)
              if (logi==1)  pak(jp)=stma(jp, kjp) 
              if (logj==1)  pal(jp)=stma(jp, ljp)
1512	   enddo
	 	
	   do 1514 jp=15,17
	   pai(jp)=0.
	   paj(jp)=0.
              pak(jp)=0.
              pal(jp)=0.
1514	   end do
	   do 1515 jp=25,27
	   pai(jp)=0.
	   paj(jp)=0.
               pak(jp)=0.
               pal(jp)=0.
1515	   end do

	   if (abs(pai(5)-paj(5))<4.e-10) pai(5)=pai(5) - 1.d-10 		  
           alpha=datan((paj(6)-pai(6))/(paj(5)-pai(5))) 
	   if (pai(5)>paj(5))   alpha=alpha + pi           
           if (alpha<0) alpha=alpha+2*pi
           ! contact point
           contactX=stma(5,i)+ (stma(5,j)-stma(5,i)) * stma(1,i) / (stma(1,i)+stma(1,j))     
           contactY=stma(6,i)+ (stma(6,j)-stma(6,i)) * stma(1,i) / (stma(1,i)+stma(1,j))    

           if ((pai(1)+paj(1)-dist)>deltamax) deltamax=pai(1)+paj(1)-dist   
           if ((pai(1)+paj(1)-dist)>deltamax)  print *, deltamax
     
           call tdm_model (pai, paj, pak, pal, alpha, e, youngr, dist, pi, fricp, numcol, contactX, contactY, logi, logj) 
      
           nbx(i)=2        ! both two particles using tdm time step
           nbx(j)=2
			                              		   	
	   do 1542 jp=15,17
	   if (logi==0)   stma(jp,i)=pai(jp)+stma(jp, i)
	   if (logj==0)   stma(jp,j)=paj(jp)+ stma(jp,j)
           if (logi==1)   stma(jp,kjp)=stma(jp,kjp) + pak(jp)
           if (logj==1)   stma(jp,ljp)=stma(jp,ljp) + pal(jp)
1542      enddo 
	   do 1543 jp=25,27
	   if (logi==0)  stma(jp,i)=stma(jp,i)+pai(jp)
	   if (logj==0)  stma(jp,j)=stma(jp,j)+paj(jp)
           if (logi==1)  stma(jp,kjp)=stma(jp,kjp)+ pak(jp)
           if (logj==1)  stma(jp,ljp)=stma(jp,ljp)+ pal(jp)
1543      enddo 

          endif
         endif                     ! end tdm

1532     continue 
1530     enddo
    
       end	

!
!
!================================
	subroutine  ghostwall (stma, numcol, nump, pai, paj, pak, pal, alpha, e, youngr, dist, pi, fricp, dtp0, kedm, dtedm, nw, walld) 
!!================================!
!  judge wheather the particles contact with the ghostwall, if it is ture 
        use face
	implicit double precision (a-h, o-z)   
	dimension stma(numcol, nump) 
        dimension pai(numcol), paj(numcol), pak(numcol), pal(numcol)
        dimension nw(15000,4) 
        dimension walld(7,4)

        numpar=size(nw,1)
        numwall=size(nw,2) 
                 
        ftan=0.
        fnor=0.

        do ni=1, numwall

        if  (nw(numpar,ni)==0)  goto 1900
        do 1905 ii=1, nw(numpar,ni)

        i=nw(ii, ni)

        if (floor(STMA(28,I)+0.4)==3)  then      !it is a real solid particle       

          logi=0
          logj=0  
    
           do 1920 j=1,14
	   pai(j)=stma(j,i)  
           paj(j)=pai(j)
1920	   enddo 
    
           do 1922 j=15,17  !accelartion 
           pai(j)=0. 
1922       enddo
 
           paj(18)=0           

           do 1925 j=25,27   !augular accelartion
           pai(j)=0.
1925      enddo

           do 1924 j=10,12    !velocity
           paj(j)=0.
1924       enddo
      
           if (abs(walld(5,ni)/walld(6,ni))<1.e-18)  then 
           paj(5)=pai(5)+1.e-15
           paj(6)=-2.*walld(7,ni)-pai(6)
           endif 

           if (abs(walld(5,ni)/walld(6,ni))>1.e18)   then
           paj(5)= -2.*walld(7,ni)-pai(5)
           paj(6)=pai(6)+1.e-10   
           endif 

           if (abs(walld(5,ni)/walld(6,ni)) >=1.e-18 .and. abs(walld(5,ni)/walld(6,ni)) <=1.e18) then            
           paj(5)= pai(5)- 2*walld(5,ni)*(walld(5,ni)*pai(5)+walld(6,ni)*pai(6)+walld(7,ni))/  &
                    (walld(5,ni)**2 + walld(6,ni)**2)
           paj(6)= pai(6)- 2*walld(6,ni)*(walld(5,ni)*pai(5)+walld(6,ni)*pai(6)+walld(7,ni))/  &
                    (walld(5,ni)**2 + walld(6,ni)**2)
           endif
     
           dist=( (pai(6)-paj(6))**2 + (pai(5)-paj(5))**2)**0.5  
                       
            if (dist<2.*pai(1)) then
                
                  if ((pai(1)+paj(1)-dist)>deltamax) deltamax=pai(1)+paj(1)-dist   
                  if ((pai(1)+paj(1)-dist)>deltamax)  print *, deltamax

               kjp=int(stma(18,i) + 0.43)
               if (kjp/=0) then 
               logi=1
               do 1927 j=1,14
               pak(j)=stma(j,kjp)
1927           enddo 
               do 1929 j=15,17
               pak(j)=0.
1929           enddo

               do 1921 j=18,24
               pak(j)=stma(j,kjp)
1921           enddo

               do 1923 j=25,27
               pak(j)=0.
1923           enddo
               endif 
          
                contactX=0.5* (pai(5)+ paj(5)) 
                contactY=0.5* (pai(6)+ paj(6)) 
                  
           alpha=datan((paj(6)-pai(6))/(paj(5)-pai(5))) 
	   if (pai(5)>paj(5))   alpha=alpha + pi           
           if (alpha<0.) alpha=alpha+2*pi       
                     
          call tdm_model (pai, paj, pak, pal, alpha, e, youngr, dist, pi, fricp, numcol, contactX, contactY, logi, logj)	

           nbx(i)=2        ! this articles using tdm time step

           if (logi==0) then                   
	   do 1926 j=15,17
	   stma(j,i)=pai(j)+stma(j,i)	
        
1926       enddo 
	   do 1928 j=25,27
           stma(j,i)=pai(j)+stma(j,i)	
     
1928       enddo      
else
           stma(15,i)=0.
           stma(16,i)=0.        
           stma(25,i)=0.
           endif 

           if (logi==1) then      
                     
	   do 1937 j=15,17
           stma(j,i)=0.
	   stma(j,kjp)=pak(j)+stma(j,kjp)	  
1937       enddo 
	   do 1939 j=25,27
           stma(j,i)=0. 
           stma(j,kjp)=pak(j)+stma(j,kjp)	
1939       enddo 
           endif 
    
        endif  
        
        endif

1905   enddo

1900   continue

       enddo
                     
       end
!
!
!================================
      subroutine tdm_model (pai, paj, pak, pal, alpha, e, youngr, dist, pi, fricp, numcol, contactX, contactY, logi, logj)	   	
!===============================
! TDM collison model 
! input: pai,paj--information of two colliding particles
!        dist--distance between two particles
!
	implicit double precision (a-h, o-z)  
	dimension pai(numcol), paj(numcol), pak(numcol), pal(numcol)
	ui=dcos(alpha)*pai(10) + dsin(alpha)*pai(11)
	vi=dcos(alpha)*pai(11) - dsin(alpha)*pai(10)
	uj=dcos(alpha)*paj(10) + dsin(alpha)*paj(11)
 	vj=dcos(alpha)*paj(11) - dsin(alpha)*paj(10)	
			
        ci=vi+pai(22)*pai(1)      ! according to clockwise
	cj=vj-paj(22)*paj(1)           
   
	pmassr=(pai(3)*paj(3))/ (pai(3)+paj(3)) 
	radr=(pai(1)*paj(1))/ (pai(1)+paj(1))  
	coe1=4/3*youngr*dsqrt(radr)
	coe2=-2*0.9129*(dlog(e)/dsqrt(dlog(e)*dlog(e)+pi*pi)) *dsqrt(pmassr*2*youngr*radr**0.5)
		
	fnor=coe1*(pai(1)+paj(1)-dist)**1.5  + coe2*(ui-uj)*(pai(1)+paj(1)-dist)**0.5	
 
	if (fnor<0) fnor=0.d0								
	ftan=fnor*fricp 
	if (ci==cj) ftan=0.d0                                                                                          
        
        if (logi==0  )  then  
  	  if (ci<cj)  	  ftani=ftan 
	  if (ci>cj) 	  ftani=-ftan	 
	pai(15)=(dcos(alpha)*(-fnor)-dsin(alpha)*ftani)/pai(3)	 
	pai(16)=(dcos(alpha)*ftani +dsin(alpha)*(-fnor))/pai(3)
	pai(25)=ftani/(pai(1) * pai(3)*0.4)   
        endif
                                   
        if ( logj==0  )  then  
  	  if (ci<cj)  	  ftanj=-ftan	  
	  if (ci>cj) 	  ftanj=ftan
	paj(15)=(dcos(alpha)*(fnor)-dsin(alpha)*ftanj)/paj(3)  
	paj(16)=(dcos(alpha)*ftanj +dsin(alpha)*(fnor))/paj(3)
	paj(25)=  -ftanj/(paj(1) * paj(3)*0.4)  
        endif    
    
        if (logi==1  ) then
        ftotal= sqrt(fnor*fnor + ftan*ftan)                                                                                   
        if (ftan/=0.) then
        theta= atan (ftan/fnor)  
        else
        theta=1.e-9
        endif  
                                                                                 
        c2c= sqrt ( (contactX-pak(5))**2 + (contactY-pak(6))**2 )     ! center to conctat point
        beta = acos ( (c2c*c2c + pai(1)*pai(1) - pai(13)*pai(13)) / (2*pai(1)*c2c) )
        zita=beta + theta
        fnorp = ftotal * dcos(zita)        ! p for prime
        ftanp = ftotal * dsin(zita)   
                 
 	  if (ci<=cj)  	  then 
           ftani=ftan 
           xi=alpha-theta
          endif
	  if (ci>cj) 	 then 
          ftani=-ftan
           xi=alpha+theta
          endif
           if (xi<0.) xi=xi+2.*pi
           if (xi>2.*Pi) xi=xi-2.*pi

           chi=datan((contactY-pak(6))/(contactX-pak(5)))
           if (pak(5)>contactX) chi=chi+pi      
           if (chi<0.) chi=chi+2.*pi
     
          if  ( (((xi-chi)<pi) .and. ((xi-chi)>0.))  .or.  (((xi-chi)<-pi) .and. ((xi-chi)>-2*pi)) )       then
          ftanpi=  -ftanp
          else
          ftanpi=  ftanp
          endif

	pak(15)=(dcos(alpha)*(-fnor)-dsin(alpha)*ftani)/pak(3)
	pak(16)=(dcos(alpha)*ftani +dsin(alpha)*(-fnor))/pak(3)    
        pak(25) =  ftanpi * c2c / (pak(3)*( pak(1)*pak(1) +  pak(1)*pak(1)*4. ) *0.2)

        endif
      
       if (logj==1 ) then
        ftotal= sqrt(fnor*fnor + ftan*ftan)
                                                                       
        if (ftan/=0.) then
        theta= atan (ftan/fnor)  
        else
        theta=1.e-9
        endif  
                                                                      
        c2c= sqrt ( (contactX-pal(5))**2 + (contactY-pal(6))**2 )     ! center to conctat point
        beta = acos ( (c2c*c2c + paj(1)*paj(1) - paj(13)*paj(13)) / (2*paj(1)*c2c) )
        zita=beta  + theta
        fnorp = ftotal * cos(zita)        ! p for prime
        ftanp = ftotal * sin(zita)    
       
 	  if (ci<=cj)  	  then 
          ftanj=-ftan 
          xi=alpha-theta
          else 
          ftanj=ftan
          xi=alpha+theta
          endif

          if (xi<0.) xi=xi+2*pi
          if (xi>2*pi) xi=xi-2*pi
          
           chi=datan((contactY-pal(6))/(contactX-pal(5)))
           if (pal(5)>contactX) chi=chi+pi
           if (chi<0.) chi=chi+2.*pi
           if (chi>2.*pi) chi=chi-2.*pi

          if  ( (((xi-chi)<pi) .and. ((xi-chi)>0.))  .or.  (((xi-chi)<-pi) .and. ((xi-chi)>-2.*pi)) )       then
          ftanpj=  ftanp
          else
          ftanpj= -ftanp
          endif
       
	pal(15)=(dcos(alpha)*(fnor)-dsin(alpha)*ftanj)/pal(3)   
	pal(16)=(dcos(alpha)*ftanj +dsin(alpha)*(fnor))/pal(3) 
        pal(25) =  ftanpj *c2c / (pal(3)*(pal(1)*pal(1) + pal(1)*pal(1)*4. ) *0.2) 
        endif

	end	

!
!=====================================
subroutine outputp (stma, nump, numcol, dtp0, tpout, istep, tp0, vmax, iip, dudt, pressF, &
                    partF, FricF, numfull, deltamax, ratio, stma2)	   	
!===============================
	use face
	implicit double precision (a-h, o-z)   
	dimension stma(numcol, nump) 
	dimension stma2(numcol, nump) 
			 	 
	character * 17 fileout
	character * 17 file2out
        character * 17 time
	character * 6 nout0	   

	! save to computer 	
	if ((iip*tpout>=tp0-dtp0) .and. (iip*tpout<=tp0)) then
          if (tp0<3.* dtp0)  print *, 'finished one roop'
         stma2=stma  
         do i=1, nump        
         if (int(stma(28,i)+0.43)==0) goto 1477
         stma2(1,i)=stma(1,i)/ratio  
         enddo    
1477     continue 

	write (nout0, '(i6.6)') iip
	fileout='file'//nout0//'.dat'
	time='time'//nout0//'.dat'	      

	open(70, file=fileout)
	write (70, 170) stma2
	 170 format(30e16.8)
	close(70)

	open(73, file=time)
	write (73, 173) dtp0, tpout, tp0, real(iip), real(numfull)
	173 format(5e18.10)
        close(73)

	iip=iip+1
	endif 

	!! this is for screen
	nsrn=10000    	! very nsrn steps, show the result in sreeen	 
	if (istep==nsrn) then
	istep=0
	tke=0.d0
	omax=0.d0
	do i=1, nump	 
	tke=tke+0.5*stma(3,i)*(stma(10,i)*stma(10,i)+stma(11,i)*stma(11,i))
	if (dabs(stma(22,i))>omax) omax=dabs(stma(22,i))	 	
	enddo
	print*,'time=', tp0, 's    dt=', dtp0,'s'
	print*,'vmax=', vmax, '  kinetic energy=', tke, '     omax=', omax,  'deltamax=', deltamax	   
        Print *, '    '        
	endif
	istep=istep+1 
	 
	end
!
!
!=====================================
      subroutine VelocityPE (stma, numcol, nump,  dtp0, kedm, dtedm, dtfree, iiv)
!===============================
        use face
	implicit double precision (a-h, o-z)   
	dimension stma(numcol, nump) 	 
      
	do i=1, nump	              
        if (int(stma(28,i)+0.43)<=0   )  goto 1418	 
              
        if (nbx(i)>=2) then
	stma(10,i) = stma(10,i)+stma(15,i)*dtp0
	stma(11,i) = stma(11,i)+stma(16,i)*dtp0
	stma(22,i) = stma(22,i)+stma(25,i)*dtp0  
	stma(5,i) = stma(5,i)+stma(10,i)*dtp0 
        stma(6,i) = stma(6,i)+stma(11,i)*dtp0
	stma(19,i) = stma(19,i)+stma(22,i)*dtp0           
        endif

        if (nbx(i)==-1  .and. kedm==50 ) then
	stma(10,i) = stma(10,i)+stma(15,i)*dtedm
	stma(11,i) = stma(11,i)+stma(16,i)*dtedm
	stma(22,i) = stma(22,i)+stma(25,i)*dtedm   
	stma(5,i) = stma(5,i)+stma(10,i)*dtedm
	stma(6,i) = stma(6,i)+stma(11,i)*dtedm	
        stma(19,i) = stma(19,i)+stma(22,i)*dtedm	
        endif

        if (nbx(i)==0  .and. iiv==500 ) then
	stma(10,i) = stma(10,i)+stma(15,i)*dtfree
	stma(11,i) = stma(11,i)+stma(16,i)*dtfree
	stma(22,i) = stma(22,i)+stma(25,i)*dtfree  
	stma(5,i) = stma(5,i)+stma(10,i)*dtfree
	stma(6,i) = stma(6,i)+stma(11,i)*dtfree
        stma(19,i) = stma(19,i)+stma(22,i)*dtfree	
        endif
                	   
       enddo
1418   continue
		   
	end
!
!
!
!=====================================
	subroutine searchsize(xlower, xupper, ylower, yupper, radmax,  mvx, mvy, cellp)
!=====================================
	use face
	! integer, allocatable:: mv (:,:,:)	
	! mv is the grid for locate the particle for find their neighbours
	! this subroutine give the size of mv
	implicit double precision (a-h, o-z)  	 
		 
        cellp=2.5d0*radmax    *2.0
	mvx=int((xupper-xlower)/cellp)+4
	mvy=int((yupper-ylower)/cellp)+4 
	    	     
	allocate(mv(mvx,mvy,100))
	end

!
!=====================================
	subroutine BoundaryCellFix (stma, xlower, xupper, ylower, yupper, radmax,  mvx, mvy, cellp, nw, walld, numcol, nump, tp0)
!=====================================
!       give the neighobring particles to any fixed boundary (such as wall) 
	! integer, allocatable:: mv (:,:,:)	
	! mv is the grid for locate the particle for find their neighbours
	! this subroutine give the size of mv
        ! k is the id of each 2-D wall
        ! numwell-- number of walls in the codes
        ! all lines can be written in form of Ax+By+C=0   
        ! there are 6 elements in each row of array walld (wall data), which means
        ! x_min, x_max, y_min, y_max, A, B and C
        ! for lines in form of x=c, the 5th element equals 1.e-30 and 6th element equals c
        ! for lines in like y=c, the 5th element equals 1.e30 and 6th element equals c

	use face
	implicit double precision (a-h, o-z)  

 !       real :: walld(6,2)=reshape ( (/ -0.005, 0.01, -0.005, 0.005, 1.e-30,  0.0,     & 
 !                                      -0.005, 0.01, -0.005, 0.005, 1.e-30,  1.2   /),(/6,2/))  
        dimension stma(numcol, nump)
        dimension dis(4)
        dimension nw(15000,4)
        dimension walld(7,4)
        numpar=size(nw,1)
        numwall=size(nw,2) 

            !print *, 'BEFROE bfc', 'MV', mv(7, 723, 1), mv(7, 723, 7)
    
            do ii=1, nump
            stma(30,ii)=0.e0
            enddo
                                                                         !        print *, 'enter BCF'    
        walld(:,1)=(/-1.005,  0.085,   -0.005,  0.005,  1.e-20,  1.,      0.0/)
        walld(:,2)=(/-0.005,  0.085,    0.995,  1.005,  1.e-20,  1.,     -1.0 /)
        walld(:,3)=(/-0.003,  0.003,   -0.005,  1.2,    1.e0,    1.e-20,  -0.e0  /)
        walld(:,4)=(/ 0.062,  0.067,   -0.005,  1.2,    1.e0,    1.e-20,  -0.065  /)
    
        walld(3,2)=0.995- 20.* tp0
        walld(4,2)=1.005- 20.* tp0
        walld(7,2)= -(1.d0 - 20.d0 * tp0)    ! 20. is the velocity dropping down

        do ni=1,numwall
        nw(numpar, ni)=0
        ! if (abs(walld(5,ni))<1.e-15)   walld(5,ni)=1.e-20 
        ! if (walld(5,ni)>1.e15)    walld(5,ni)=1.e20 
        ! if (walld(5,ni)<-1.e15)   walld(5,ni)=-1.e20 
        ! if (abs(walld(6,ni))<1.e-15)   walld(6,ni)=1.e-20 
        ! if (walld(6,ni)>1.e15)    walld(6,ni)= 1.e20 
        ! if (walld(6,ni)<-1.e15)   walld(6,ni)=-1.e20 
        enddo		 

        do i=1, mvx
        do j=1, mvy

        do ni=1,numwall 

        x= (i+0.5 -3)* cellp
        y= (j+0.5 -3)* cellp

       if (x<walld(1,ni) .or. x>walld(2,ni) .or. y<walld(3,ni) .or. y>walld(4,ni))  goto 2500   
                                                                     
        node=1
           do ii=-1,1,2
           do jj=-1,1,2     
            xcorner=x+ii*cellp *1.5 
            ycorner=y+jj*cellp *1.5
            
             if (abs(walld(5,ni)/walld(6,ni))<1.e-18)  d=ycorner + walld(7,ni)     ! for curve y=const. 
             if (abs(walld(5,ni)/walld(6,ni))>1.e18)   d=xcorner + walld(7,ni)    ! for curve x=const.

             if (abs(walld(5,ni)/walld(6,ni)) >=1.e-18 .and. abs(walld(5,ni)/walld(6,ni)) <=1.e18) then
                 d= ((walld(5,ni)*xcorner + walld(6,ni)*ycorner + walld(7,ni)) )/ ((walld(6,ni)**2 +walld(5,ni)**2)**0.5)   
             endif  
             dis(node)=d     
            node=node+1     
           enddo
           enddo 
 
        dtd = minval(dis)*maxval(dis)
                                         
        if (dtd<=0.d2 ) then    
                    
         if (mv(i,j,1) ==0)  goto 2510                           
         do ia=2, (mv(i,j,1)+1)      
                     
         if (mv(i,j,ia)==0) goto 1520     
         nw(numpar,ni)=nw(numpar,ni)+1 
         nw(nw(numpar,ni),ni)=mv(i,j,ia)                   
           ! stma(30, mv(i,j,ia)) = 1.00*ni    !add color (given value) to each particles if they closing to a wall
           if (nw(numpar, ni)>(numpar-20)) then
           print *, 'arrangement NW is full. please change the size of NW or check the codes!',  nw(numpar,ni), numpar  
           stop
           endif
1520    continue
         enddo
2510    continue
        endif  

2500    continue

        enddo

        enddo
        enddo

	end

!=====================================
      subroutine VelocityMM (stma, numcol, nump,  dtp0, kedm, dtedm, stma3, dtfree)
!===============================
        use face
	implicit double precision (a-h, o-z)   
	dimension stma(numcol, nump), stma3(numcol, nump)
       
	do i=1, nump	              
        if (int(stma(28,i)+0.43)<=0   )  goto 1718	               
        if (nbx(i)>=2) then
	stma(10,i) = stma(10,i)+stma(15,i)*dtp0
	stma(11,i) = stma(11,i)+stma(16,i)*dtp0
	stma(22,i) = stma(22,i)+ 0.5* (stma(25,i)+stma3(25,i)) *dtp0  
	stma(5,i) = stma(5,i)+ 0.5*(stma(10,i)+stma3(10,i)) *dtp0
	stma(6,i) = stma(6,i)+ 0.5*(stma(11,i)+stma3(11,i))*dtp0	
	stma(19,i) = stma(19,i)+ 0.5*(stma(22,i)+stma3(22,i))*dtp0         
        endif

        if (nbx(i)==-1  .and. kedm==50 ) then
	stma(10,i) = stma(10,i)+stma(15,i)*dtedm
	stma(11,i) = stma(11,i)+stma(16,i)*dtedm
	stma(22,i) = stma(22,i)+ 0.5* (stma(25,i)+stma3(25,i)) *dtedm
	stma(5,i) = stma(5,i)+ 0.5*(stma(10,i)+stma3(10,i)) *dtedm
	stma(6,i) = stma(6,i)+ 0.5*(stma(11,i)+stma3(11,i))*dtedm
	stma(19,i) = stma(19,i)+ 0.5*(stma(22,i)+stma3(22,i))*dtedm
        endif

        if (nbx(i)==0  .and. iiv==500 ) then
	stma(10,i) = stma(10,i)+stma(15,i)*dtfree
	stma(11,i) = stma(11,i)+stma(16,i)*dtfree
	stma(22,i) = stma(22,i)+ 0.5* (stma(25,i)+stma3(25,i)) *dtfree
	stma(5,i) = stma(5,i)+ 0.5*(stma(10,i)+stma3(10,i)) *dtfree
	stma(6,i) = stma(6,i)+ 0.5*(stma(11,i)+stma3(11,i))*dtfree
	stma(19,i) = stma(19,i)+ 0.5* (stma(22,i)+stma3(22,i)) *dtfree
        endif                	   
	   if ((stma(10, i)>100) .or. (stma(11, i)>100)) then
	     print *, '~~~~~~~~~~~ERROR'		
	   endif
       enddo
1718   continue			
       stma3=stma	   
       end

!=====================================
      subroutine VelocityVV (stma, numcol, nump,  dtp0, kedm, dtedm, stma3, dtfree)
!===============================
        use face
	implicit double precision (a-h, o-z)   
	dimension stma(numcol, nump), stma3(numcol, nump)
      
	do i=1, nump	              
        if (int(stma(28,i)+0.43)<=0)  goto 1818	               
        if (nbx(i)>=2) then
	stma(5,i)=stma(5,i)+ stma3(10,i) *dtp0 + 0.5* stma3(15,i) * dtp0 * dtp0 
	stma(6,i)=stma(6,i)+ stma3(11,i) *dtp0 + 0.5* stma3(16,i) * dtp0 * dtp0 
	stma(19,i)=stma(19,i)+ stma3(22,i) *dtp0 + 0.5* stma3(25,i) * dtp0 * dtp0 
	stma(10,i)=stma(10,i)+ 0.5*(stma(15,i)+stma3(15,i)) *dtp0
	stma(11,i)=stma(11,i)+ 0.5*(stma(16,i)+stma3(16,i))*dtp0
	stma(22,i)=stma(22,i)+ 0.5* (stma(25,i)+stma3(25,i)) *dtp0           
        endif

        if (nbx(i)==-1  .and. kedm==50 ) then
	stma(5,i)=stma(5,i)+ stma3(10,i) *dtedm + 0.5* stma3(15,i) * dtedm * dtedm
	stma(6,i)=stma(6,i)+ stma3(11,i) *dtedm + 0.5* stma3(16,i) * dtedm * dtedm 
	stma(19,i)=stma(19,i)+ stma3(22,i) *dtedm + 0.5* stma3(25,i) * dtedm * dtedm 
	stma(10,i)=stma(10,i)+ 0.5*(stma(15,i)+stma3(15,i)) *dtedm
	stma(11,i)=stma(11,i)+ 0.5*(stma(16,i)+stma3(16,i))*dtedm
	stma(22,i)=stma(22,i)+ 0.5* (stma(25,i)+stma3(25,i)) *dtedm      
        endif

        if (nbx(i)==0  .and. iiv==500 ) then
	stma(5,i)=stma(5,i)+ stma3(10,i) *dtfree + 0.5* stma3(15,i) * dtfree * dtfree
	stma(6,i)=stma(6,i)+ stma3(11,i) *dtfree + 0.5* stma3(16,i) * dtfree * dtfree 
	stma(19,i)=stma(19,i)+ stma3(22,i) *dtfree + 0.5* stma3(25,i) * dtfree * dtfree 
	stma(10,i)=stma(10,i)+ 0.5*(stma(15,i)+stma3(15,i)) * dtfree
	stma(11,i)=stma(11,i)+ 0.5*(stma(16,i)+stma3(16,i)) * dtfree
	stma(22,i)=stma(22,i)+ 0.5*(stma(25,i)+stma3(25,i)) * dtfree    
        endif
                	   
	if ((stma(10,i)>100) .or. (stma(11,i)>100)) then
	  print *, '~~~~~~~~~~~ERROR'		
	endif
1818   continue
       enddo	
       stma3=stma			   
       end




!=====================================
	subroutine vmmax (stma, vmax, numcol, nump)		 
!=====================================
        use face
	implicit double precision (a-h, o-z)   
	dimension stma(numcol, nump)  

	vmax=0.d0
	do 11103 i=1, nump
          if (int(stma(28, i)+0.4)<=0) goto 11123
	  vsp=dsqrt(stma(10,i)*stma(10,i)+stma(11,i)*stma(11,i))                                           !vsp: velocity of single particle!
	  if (vsp>vmax) vmax=vsp          
11103   enddo		
11123   continue 
	print *, 'vmax is', vmax
         end 

