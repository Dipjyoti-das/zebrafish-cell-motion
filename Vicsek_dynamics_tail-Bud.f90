!!This code is written by DIPJYOTI DAS (dipjyoti.das@yale.edu, dipjyoti.das85@gmail.com)
!!The code performs Vicsek dynamics of soft particles (cells) in 2D.
!!! The cells are inside a 2D horse-shoe geometry resembeling zebrafish tailbud,
!!! and the boundaries are rigid and reflective.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


MODULE numz
  IMPLICIT NONE
  INTEGER,PARAMETER:: DP=KIND(1.0d0)
  REAL(DP),PARAMETER:: pi=3.14159265358979_DP
END MODULE numz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE constants
 USE numz
 IMPLICIT NONE
REAL(DP),PARAMETER::mean_noise=0.0_DP,var_noise=1.0_DP !// Noise parameters

REAL(DP),PARAMETER:: L=20.0_DP  !! initial length
REAL(DP),PARAMETER:: R=11.0_DP  !! Radius of PZ (=(L+2)/2)
REAL(DP),PARAMETER:: mu=1.0_DP  !! mobility
REAL(DP),PARAMETER:: Vo=1.00_DP  !! magnitude of self-propulsion velocity
REAL(DP),PARAMETER:: tau=1.0_DP  !! relaxation time for the alingment of self-propulsion vectors
REAL(DP),PARAMETER:: Frep=30.0_DP 
REAL(DP),PARAMETER:: Fadh=0.75000_DP !! Maximums of the Replusive and adhesive forces
REAL(DP),PARAMETER:: ao=(5.0_DP/6.0_DP)*0.5_DP !! average radius of cells
REAL(DP),PARAMETER:: pack_frac=0.95_DP  !! packing fraction of cells (initially)
REAL(DP),PARAMETER:: So=pi*ao**2 !!average cell area

!!!!%%%%%% initial numbers of cells in DM, PZ, PSM (all integers) %%%%%%%%%%%%
!!!! This determines how many cells are inside each region (initially)

INTEGER,PARAMETER:: NPtop=  NINT( 0.5_DP*(pi*R**2)*(pack_frac/So) )  
INTEGER,PARAMETER:: NP=   NINT( (1/5.0_DP)*L*(L+2.0)*(pack_frac/So) )
INTEGER,PARAMETER:: ND=   NINT( (1.0_DP/4.0_DP)*(L**2)*(pack_frac/So) )
INTEGER,PARAMETER:: NPSM= NINT( (1.0_DP/4.0_DP)*(L**2)*(pack_frac/So) )       
INTEGER,PARAMETER:: Ntot_initial = NP+ND+NPSM+NPtop
INTEGER,PARAMETER:: Narray=Ntot_initial*1000

!!!!! time & sampling %%%%%%%%%%%%%%%%%%%%%%%%
INTEGER,PARAMETER:: Maxstep=11000  !3600 !!total simulation time (in integral steps)
REAL(DP),PARAMETER:: h=0.005     !!time-increment
INTEGER,PARAMETER:: sample_step=36 !36 !! after how many step-intervals, we are going to sample the system

!!!!! CELL INPUT VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
INTEGER,PARAMETER:: N_in=NINT(L/(4.0_DP*ao)) !!! #cells are introduced in every (gamma/Vo) sec.
REAL(DP),PARAMETER:: gamma=1.0_DP
INTEGER,PARAMETER:: step_input = NINT(gamma/(1.0_DP*h*N_in)) !! Step #, after which 1 cell is introduced
END MODULE constants










!!!!!!!!!!!!!!!!!! MAIN PROGRAM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM horseshoe
 USE numz
 USE constants 

 IMPLICIT NONE 

!!!!!!! Difine variables !!!!!!!!!!!

!!!!%%%%%%%%%%%%%%%%%% position and velocity arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%
REAL(DP),DIMENSION(1:Narray)::x,y,Vx,Vy !!position & velocity of each particle
REAL(DP),DIMENSION(1:Narray)::phi !! self-propulsion direction of each particle
REAL(DP),DIMENSION(1:Narray)::radi !! radius of each particle
REAL(DP),DIMENSION(1:Narray)::dx,dy,dphi !!! Increments
REAL(DP),DIMENSION(1:Narray):: x_prev, y_prev
REAL(DP),DIMENSION(1:Narray):: eta !! noise-strength array for each cell
!!!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!!! sampling %%%%%%%%%%%%%%%%%%%%%%%%
INTEGER:: sample_counter=1
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!!%%%%% variables for the horse-shoe geometry %%%%%%%%
REAL(DP):: gap,xDMmin,xDMmax,xPSMmin,xPSMmax,xmin,xmax
REAL(DP):: ymax,yPZ,yDM,ymin,ylow,ytop,PZcenter
!!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!

!!!%%%%%%%% LOCAL variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INTEGER::idum,k,Ntot,step,i,j,n,Ntemp,dN
REAL(DP):: ran2,rc,fxsum,fysum,fx,fy,t
REAL(DP):: xtemp,ytemp,V_cross_n,V_mode,noise,Vel_mod
REAL(DP):: volume, volume_temp,dvol,dl
REAL(DP)::vx0,vy0,ux0,uy0
!REAL(DP):: Vo
!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!! Initialize all arrays to zero
x=0.0_DP
y=0.0_DP
Vx=0.0_DP
Vy=0.0_DP
phi=0.0_DP
radi=0.0_DP
dx=0.0_DP
dy=0.0_DP
dphi=0.0_DP
x_prev=0.0_DP
y_prev=0.0_DP

eta=0.900000_DP !!! Initialize noise array 

!!!!! Build the initial horse-shoe geometey!!!!!!!!!!!!!!
        gap=1.0_DP
	xDMmin=-L/4.0_DP
	xDMmax= L/4.0_DP
	xPSMmin=-L/4.0_DP -gap
	xPSMmax=+L/4.0_DP +gap
	xmin=xPSMmin-L/4.0_DP
	xmax=xPSMmax+L/4.0_DP   
    ymax= L/2.0+L/5.0 !(3.0_DP*L)/4.0_DP !(3.0_DP*L)/4.0_DP 
    ytop= ymax+R
    yPZ=   L/2.0_DP !(4.0_DP*L)/5.0_DP ! 
    ymin=0.0
    ylow=0.0
!!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OPEN(unit=21,file='Cell_Flow_noise-0.9.dat',status='unknown')

CALL SYSTEM_CLOCK(COUNT=idum)
!!!!!****** INITIALLY CELLS TAKEN RANDOMLY IN DIFFERENT REGIONS*************

!!!! DM region : put ND number of cells
     do k=1,ND
	   radi(k)=ao+0.1_DP*(ran2(idum)-0.5_DP)
	   rc=radi(k)
	   x(k)=xDMmin+rc+(xDMmax-xDMmin-2.0_DP*rc)*ran2(idum)
	   y(k)=2.0_DP*rc+(yPZ-rc)*ran2(idum)
	   phi(k)=2*pi*ran2(idum)  !pi/2.0_DP 
	   x_prev(k)=x(k)
	   y_prev(k)=y(k)
	enddo
	
!!! PZ Region: put NP number of cells

	do k=1+ND,(ND+NP)
	   radi(k)=ao+0.1_DP*(ran2(idum)-0.5_DP)
	   rc=radi(k)
	   x(k)=xmin+rc+(xmax-xmin-2.0_DP*rc)*ran2(idum)
	   y(k)=yPZ+rc+(ymax-yPZ-2.0_DP*rc)*ran2(idum)
	   phi(k)=2*pi*ran2(idum)
	   x_prev(k)=x(k)
	   y_prev(k)=y(k)
	enddo
	
!!! PZtop Region: put NPtop number of cells

	do k=(ND+NP+1),(ND+NP+NPtop)
	   radi(k)=ao+0.1_DP*(ran2(idum)-0.5_DP)
	   rc=radi(k)
	   x(k)= R*dsqrt(ran2(idum))*dcos(pi*ran2(idum))
	   y(k)= ymax + R*dsqrt(ran2(idum))*dsin(pi*ran2(idum))                         
	   phi(k)=2*pi*ran2(idum)
	   x_prev(k)=x(k)
	   y_prev(k)=y(k)
	enddo	
	
	
!!! PSM-left Region: put NPSM/2 number of cells

	do k=(ND+NP+NPtop+1),NINT(ND+NP+NPtop+NPSM/2.0)
	   radi(k)=ao+0.1_DP*(ran2(idum)-0.5_DP)
	   rc=radi(k)
	   x(k)=xmin+rc+(xPSMmin-xmin-2.0_DP*rc)*ran2(idum)
	   y(k)=2*rc+(yPZ-rc)*ran2(idum)
	   phi(k)= 2*pi*ran2(idum)  !3.0_DP*pi/2.0_DP
	   x_prev(k)=x(k)
	   y_prev(k)=y(k)
	enddo
	
!!! PSM-right Region: put NPSM/2 number of cells

	do k=1+NINT(ND+NP+NPtop+NPSM/2.0),Ntot_initial
           radi(k)=ao+0.1_DP*(ran2(idum)-0.5_DP)
	   rc=radi(k)
	   x(k)=xPSMmax+rc+(xmax-xPSMmax-2.0_DP*rc)*ran2(idum)
	   y(k)=2*rc+(yPZ-rc)*ran2(idum)
	   phi(k)=2*pi*ran2(idum)  !3.0_DP*pi/2.0_DP
	   x_prev(k)=x(k)
	   y_prev(k)=y(k)
	enddo
!!!!!***********************************************************************


!!! INITIALIZATION !!!!!!!!!!!!!

volume=((xmax-xmin)*ymax-2.0_DP*gap*yPZ)+ 0.500_DP*pi*R**2
 Ntot=Ntot_initial
 t=0.0_DP
sample_counter=1


DO step=1,Maxstep !!!! time-loop starts

t=t+h

   Do i=1,Ntot   !!!! BEGIN Cell Index Loop
	
        CALL SYSTEM_CLOCK(COUNT=idum)

         xtemp=x(i)
	     ytemp=y(i)
	     rc=radi(i)

   !!Calculate the force on ith cell due to all others **********************
       fxsum=0.0_DP
       fysum=0.0_DP
	 
        do j=1,Ntot   !!! BEGIN: Sum the forces Loop	         

	   call cell_force(x(i),y(i),x(j),y(j),radi(i),radi(j),fx,fy)	
	         fxsum=fxsum+fx
	         fysum=fysum+fy
	
        end do        !!!!END: Sum the forces Loop
        
        
   !!********************************************************************
              dx(i)=Vo*dcos(phi(i))*h+ mu*h*fxsum
	    	  dy(i)=Vo*dsin(phi(i))*h+ mu*h*fysum

               x(i)=x(i)+dx(i)
               y(i)=y(i)+dy(i)

       CALL boundary(ymax,yPZ,rc,xtemp,ytemp,x(i),y(i)) !!! reflecting boundary

 

     !!	******************UPGRADE THE ANGLES************************
	      Vx(i)=(x(i)-xtemp)/h
	      Vy(i)=(y(i)-ytemp)/h
	      V_cross_n=Vy(i)*dcos(phi(i))-Vx(i)*dsin(phi(i))
	      V_mode=dsqrt(Vx(i)**2+Vy(i)**2)
              
         !CALL  gasdev(noise,idum,mean_noise,var_noise)
              
              noise= pi*( 2.0_DP*ran2(idum)-1.0_DP)

	      dphi(i)=(dasin(V_cross_n/V_mode)*h)/tau + eta(i)*noise*dsqrt(h)
	      phi(i)=phi(i)+dphi(i)
     !!	******************UPGRADE THE ANGLES************************

     !!************* PRINT full config. after each sample_step******
     IF(mod(step,sample_step)==0)THEN
          vx0=(x(i)-x_prev(i))/(DBLE(sample_step*h))
	      vy0=(y(i)-y_prev(i))/(DBLE(sample_step*h))
	      ux0=vx0/dsqrt(vx0**2+vy0**2)
	      uy0=vy0/dsqrt(vx0**2+vy0**2) 
	      write(21,*) sample_counter,x(i),y(i),ux0,uy0,i
	      x_prev(i)=x(i)
	      y_prev(i)=y(i)
     END IF
     !!************************************************************

        END DO !!!! END Cell Index Loop
        
  IF(mod(step,sample_step)==0)THEN   
      sample_counter=sample_counter+1
  END IF
               
        
  Ntemp=Ntot
  Volume_temp=Volume
!!!%%%%%%% INPUT OF NEW CELLS AT ADM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 IF(mod(step,step_input) == 0)THEN
	
        do n=Ntot+1,Ntot+1
           radi(n)=ao+0.1_DP*(ran2(idum)-0.5_DP)
	   rc=radi(n)
	   x(n)=xDMmin+rc+(xDMmax-xDMmin-2.0_DP*rc)*ran2(idum)
	   y(n)=2.0_DP*rc !!!!!!  we put cells at lower layers, fixed y
	   phi(n)=pi/2.0_DP
	    x_prev(n)=x(n)
	    y_prev(n)=y(n)    
	          
	  end do
	  
    Ntot=Ntot+1
  
!!!%%%%%%%%%%%% ADJUST VOLUME TO INCORPORATE NEW CELLS %%%%%%%%%%%%
        dN=(Ntot-Ntemp)
	if(Ntot<Ntemp) dN=0
	dvol=DBLE(dN)*volume_temp/DBLE(Ntemp)
	
	dl=dvol/(xmax-xmin-2.0_DP*gap)
	
	ymax=ymax+dl
	ytop=ytop+dl
	yPZ=yPZ+dl
	ylow=ylow+dl !!!!These defines modified boundaries -- use this
!	yDM=yDM+dl 
volume=((xmax-xmin)*ymax-2.0_DP*gap*yPZ)+ 0.500_DP*pi*R**2
    

END IF	
!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

END DO !!!! time-loop ends

CLOSE(21)

END PROGRAM horseshoe

!!%%%%%%%%%%%%%%%%%%%%% MAIN PROGRAM ENDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

















!!!!/////// CELL-CELL force subroutine////////////////////////////////////

!!! This subroutine Computes the force on a cell by other cells

subroutine cell_force(x1,y1,x2,y2,r1,r2,fx0,fy0)
        USE numz
        USE constants
        IMPLICIT NONE
	
        REAL(DP),PARAMETER::del=1.2_DP
	REAL(DP),INTENT(IN)::x1,y1,x2,y2,r1,r2
        REAL(DP),INTENT(OUT):: fx0,fy0
        REAL(DP)::x12,y12,d12,req,r0

	x12=x1-x2
	y12=y1-y2
	d12=dsqrt(x12**2+y12**2)
	
	req=(r1+r2)
	r0=del*req
	
	if(d12<req .and. d12 > 0.0_DP)then
	  fx0=Frep*(x12/d12)*(req-d12)/req
	  fy0=Frep*(y12/d12)*(req-d12)/req
	endif
	
	if(d12<= r0 .and. d12 >= req)then
	  fx0=-Fadh*(x12/d12)*(d12-req)/(r0-req)
	  fy0=-Fadh*(y12/d12)*(d12-req)/(r0-req)
	endif
	
	if(d12 > r0 .or. d12==0.0_DP)then
	  fx0=0.0_DP
	  fy0=0.0_DP
	endif
	
	return
end subroutine cell_force

!!!!/////////////////////////////////////////////////////////////////////////



!!!!/////// Reflecting boundary condition subroutine////////////////////////////////////

!!! This subroutine imposes the reflecting boundary condition on each cell

subroutine boundary(ymax,yPZ,rc,xi,yi,xf,yf)
        USE numz
        USE constants
        IMPLICIT NONE
	
	REAL(DP),INTENT(IN)::ymax,yPZ,rc,xi,yi
        REAL(DP),INTENT(INOUT):: xf,yf
        REAL(DP):: gap,xDMmin,xDMmax,xPSMmin,xPSMmax,xmin,xmax,ymin
        REAL(DP):: Rad,Rf,alpha
        
        gap=1.0_DP
	xDMmin=-L/4.0_DP
	xDMmax= L/4.0_DP
	xPSMmin=-L/4.0_DP -gap
	xPSMmax=+L/4.0_DP +gap
	xmin=xPSMmin-L/4.0_DP
	xmax=xPSMmax+L/4.0_DP
        ymin=0.0_DP
        Rad=11.00_DP

         !! For DM & PSM regions :::::::::::::::::::::::::::
	  if( yi <= yPZ)then  

        if(xi >= xDMmin .and. xi <= xDMmax)then !!! DM region
        
	    if(xf <= (xDMmin+rc)) xf=2.0_DP*(xDMmin+rc)-xf
	    if(xf >= (xDMmax-rc)) xf=2.0_DP*(xDMmax-rc)-xf
        if(yf <= rc) yf=2.0_DP*rc-yf
        
         elseif(xi >= xPSMmax)then !! Right-PSM
            if(xf <= (xPSMmax+rc)) xf=2.0_DP*(xPSMmax+rc)-xf
	        if(xf >= (xmax-rc)) xf=2.0_DP*(xmax-rc)-xf
	        if(yf <= rc) yf=2.0_DP*rc-yf

         elseif( xi <= xPSMmin)then !! Left-PSM
            if(xf >= (xPSMmin-rc)) xf=2.0_DP*(xPSMmin-rc)-xf
            if(xf <= (xmin+rc)) xf=2.0_DP*(xmin+rc)-xf
            if(yf <= rc) yf=2.0_DP*rc-yf
             
          end if
           
        end if
	
       !!PZ region:::::::::::::::::::::::::::
       if(yi > yPZ .and. yi<= ymax)then

        if(xf <= (xmin+rc)) xf=2.0_DP*(xmin+rc)-xf
	    if(xf >= (xmax-rc)) xf=2.0_DP*(xmax-rc)-xf
	    
	    
	    !if(yf >= (ymax-rc)) yf=2.0_DP*(ymax-rc)-yf 
	    
	    !LEFT
	    if(xf >= (xPSMmin-rc) .and. xf <= (xDMmin+rc))then
	      if(yf <= (yPZ+rc)) yf=2.0_DP*(yPZ+rc)-yf
	    endif
	
	    !RIGHT
	    if(xf >= (xDMmax-rc) .and. xf <= (xPSMmax+rc))then
	      if(yf <= (yPZ+rc)) yf=2.0_DP*(yPZ+rc)-yf
	    endif

       end if
       
 if( yi> ymax)then      
!!!!! top rounded boundary !!!!!!!!
	    Rf=dsqrt(xf**2+ (yf-ymax)**2)
	    
	if(Rf >= (Rad-rc))then
	    
	    alpha= atan((yf-ymax)/xf)
	    
	    if(alpha>0)then
	    xf= (2*(Rad-rc) - Rf)*dcos(alpha)
	    yf= (2*(Rad-rc) - Rf)*dsin(alpha)+ ymax 
	    elseif(alpha<0)then
	    xf= (2*(Rad-rc) - Rf)*dcos(pi-abs(alpha))
	    yf= (2*(Rad-rc) - Rf)*dsin(pi-abs(alpha))+ ymax
	    end if
	   
	end if
 end if

	return
end subroutine boundary

!!!!/////////////////////////////////////////////////////////////////////////










!!!!/////// Uniform Random number generators////////////////////////////////////

FUNCTION ran2(idum)
  USE numz
  IMPLICIT NONE
  REAL(DP):: ran2
  !INTEGER,INTENT(inout),OPTIONAL::idum
  INTEGER,INTENT(inout)::idum
  INTEGER,PARAMETER::IM1=2147483563,IM2=2147483399,IMM1=IM1-1
  INTEGER,PARAMETER::IA1=40014,IA2=40692,IQ1=53668
  INTEGER,PARAMETER::IQ2=52774,IR1=12211,IR2=3791   
  INTEGER,PARAMETER::NTAB=32,NDIV=1+IMM1/NTAB
  REAL(DP),PARAMETER::AM=1.0_DP/IM1,EPS=1.2e-7,RNMX=1.0_DP-EPS
  INTEGER::idum2,j,k,iv(NTAB),iy
  SAVE iv,iy,idum2
  DATA idum2/123456789/, iv/NTAB*0/, iy/0/
  IF (idum<0) THEN
     idum=MAX(-idum,1)
     idum2=idum
      DO j=NTAB+8,1,-1
         k=idum/IQ1
         idum=IA1*(idum-k*IQ1)-k*IR1
         IF (idum<0) idum=idum+IM1
         IF (j.LE.NTAB) iv(j)=idum
      ENDDO
      iy=iv(1)
   ENDIF
   k=idum/IQ1
   idum=IA1*(idum-k*IQ1)-k*IR1
   IF (idum<0) idum=idum+IM1
   k=idum2/IQ2
   idum2=IA2*(idum2-k*IQ2)-k*IR2
   IF (idum2<0) idum2=idum2+IM2
   j=1+iy/NDIV
   iy=iv(j)-idum2
   iv(j)=idum
   IF(iy.LT.1)iy=iy+IMM1
   ran2=MIN(AM*iy,RNMX)
   RETURN
 END FUNCTION ran2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!/////// Gaussian random no. generator \\\\\\\\\\\\\\\\\\\\\\\\\\\
 Subroutine gasdev(g2,idum,mean,variance) 
  use numz
  Implicit none
      INTEGER,INTENT(INOUT)::idum
      REAL(DP),INTENT(IN)::variance,mean
      REAL(DP),INTENT(OUT)::g2
      REAL(DP)::ran2,g1
      INTEGER:: iset
      REAL(DP):: fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
    !  if (iset.eq.0) then
        DO
        v1=2.0_DP*ran2(idum)-1._DP
        v2=2._DP*ran2(idum)-1._DP
        rsq=v1**2+v2**2
        if((rsq<1.0).AND.(rsq/=0.0))EXIT
        ENDDO
        fac=variance*DSQRT(-2._DP*log(rsq)/(rsq))
        g1=v1*fac+mean
        g2=v2*fac+mean      
      END subroutine gasdev

!!!!!!!!!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
