INCLUDE "./pce-subs-2.f90"

!ifort pce-code.f90 -o pce.x
!./pce.x dirgap indgap thickmaxm < abs_coef > saida

program main

	use constants
	implicit none
	
	integer :: i,erro,iaux

	integer,parameter :: nsolar=2002
	integer :: noptics
	integer,parameter :: nvec=6001
	integer,parameter :: nth= 500

	double precision,dimension(nsolar,2) :: solarinc
	
	double precision,allocatable,dimension(:,:) :: abscoef

	double precision,dimension(nvec,2) :: isolar,isolarm,sepflux
	double precision,dimension(nvec,2) :: iabscoef,iabscoefsq
	double precision,dimension(nvec,2) :: bbir,bbpflux
	double precision,dimension(nvec,2) :: absc
	double precision,dimension(nvec-1,2) :: absc2		

	double precision,dimension(nvec,2) :: auxvec !vetor utilizado para integracoes

	double precision :: e0, ef
	double precision :: res,ephoton,eaux
	double precision :: dirgap,indgap

	double precision :: j0r,j0,jsc1,voc,j01,j1,j2,ff

	double precision :: pin,pm,voct,vmaxt,vmax

	double precision :: deltag,fr

	double precision :: pce,pcemax,thickness,negpower

	double precision,dimension(7) :: aflag 

	character(len=1) :: flag

	double precision :: thickmaxm !maximum thicknes in  m

	integer :: imax
	
	double precision :: pceaux	


	double precision :: faux1,faux2,faux3,faux4,faux5,faux6,faux7,faux8,faux9
	
	CHARACTER(LEN=30) :: Format


	!variaveis getarg

  	character*30 arg
  	integer ios


      CALL getarg(1, arg)
		read(arg,*,iostat=ios) dirgap

      CALL getarg(2, arg)
		read(arg,*,iostat=ios) indgap

      CALL getarg(3, arg)
		read(arg,*,iostat=ios) thickmaxm
		
      CALL getarg(4, arg)
		read(arg,*,iostat=ios) noptics		


	OPEN(UNIT=100, FILE= "am1.5G.dat",STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de entrada solar"

	OPEN(UNIT=200, FILE= "PCE-Limit.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de saida PCE-Limit" 
   	    	   	
	

	deltag= dirgap-indgap
	fr = dexp(-deltag/(keV*tkv))


	do i=1,nsolar

		read(100,*) solarinc(i,1),solarinc(i,2)

	end do

		!read(*,*) flag
		
	allocate(abscoef(noptics,2))


	do i=1,noptics

		read(*,*) abscoef(i,1),aflag(2),aflag(3),aflag(4),aflag(5)
		abscoef(i,2) = (aflag(2)+aflag(3)+aflag(4))*1E04		
		!abscoef(i,2) = aflag(2)+aflag(3)+aflag(4)+2.0*(aflag(5)+aflag(6)+aflag(7))

		aflag = 0.0

		!write(203,*) abscoef(i,1),abscoef(i,2)
	end do

	!homogeneizando os dados de input mediante interpolacao	
	e0 = solarinc(1,1)
	ef = solarinc(nsolar,1)

	do i=1,nvec

		ephoton = e0+(ef-e0)*dble((i-1)/(nvec-1.0))

		eaux= ((c*hev)/(ephoton+0.00000001))*10**9

		call interp1D(nsolar,solarinc,ephoton,res)
		
		isolar(i,1) = ephoton
		isolar(i,2) = res
		
		call interp1D(noptics,abscoef,eaux,res)
		

		iabscoefsq(i,1) = ephoton
		
		if (eaux .le. dirgap) then
		
		iabscoef(i,2) = 0.0
		
		else
		
		iabscoef(i,2) = res*100.00
		
		end if
		
		iabscoef(i,1) = ephoton
		
		if (eaux .ge. dirgap ) then
		iabscoefsq(i,2) = 1.0
		else
		iabscoefsq(i,2) = 0.0
		end if	
		
	

	end do

	do i=1,nvec
	isolarm(i,1)=isolar(i,1)*1E-9
	isolarm(i,2)=isolar(i,2)
	end do

	do i=1,nvec
		sepflux(i,1) = isolarm(i,1)
		sepflux(i,2) = isolarm(i,2)*(isolarm(i,1)/(h*c))
	end do

	call blackbody(nvec,isolarm,bbir,bbpflux)

	!calculo de pin, potencia incidente
	call integral1D(nvec,isolar,pin)

	!write(*,*) pin

	!call slme(nvec,bbpflux,sepflux,iabscoef,isolar,thick,pin,fr,pm,jsc1,voc,pcemax)
	
	write(200,*) "SQ-Limit"
	call sq(nvec,bbpflux,sepflux,iabscoefsq,isolar,pin,1d0,pm,jsc1,j01,j1,j2,vmaxt,voct,ff,pce)
	write(200,*) '#Jsc:',jsc1,'W/Vm^{2}'
	write(200,*) '#Vmax:',vmaxt,'V'
	write(200,*) '#Voc:',voct,'V'
	write(200,*) "#SQ-PCE",pce*100,"%"
	write(200,*) "#FF",ff*100,"%"	
	
	write(200,*) 	
	write(200,*) "SLME-Limit"	
	call sq(nvec,bbpflux,sepflux,iabscoefsq,isolar,pin,fr,pm,jsc1,j01,j1,j2,vmaxt,voct,ff,pce)
	write(200,*) '#Jsc:',jsc1,'W/Vm^{2}'
	write(200,*) '#Vmax:',vmaxt,'V'
	write(200,*) '#Voc:',voct,'V'
	write(200,*) '#fr:',fr	
	write(200,*) "#SLME_max-PCE",pce*100,"%"
	write(200,*) "#FF",ff*100,"%"			
	
	pceaux = 0.0

	write(*,*)

	
	!thickmaxm=thickmaxm*1E-6

	!call absorbance(nvec,iabscoef,thickmaxm,absc)
	!absc2=absc(2:nvec,:)
	
	!do i=1,nvec-1
	
	!	write(*,*) absc2(i,1),absc2(i,2)
	
	!end do

	write(*,*) "# thickness,pce,jmax,j0,jsc,vmax,voc,ff"
	
	do i=1,nth
		
		thickness = thickmin+(thickmaxm-thickmin)*dble((i-1.0)/(nth-1.0))
		
		call absorbance(nvec,iabscoef,thickness,absc)
		absc2=absc(2:nvec,:)

		call slme(nvec-1,voct,bbpflux,sepflux,absc2,isolar,pin,fr,pm,jsc1,j01,j1,j2,vmax,voc,ff,pce)		

	Format = "(2E15.3,6E15.3)"
	!write(*,"(8F10.3)") real(thickness*1E6),real(pce*100.0),real(j1),real(j01),real(jsc1),real(voc),real(j2),real(pm)
	write(*,Format) real(thickness*1E6),real(pce*100.0),real(j1),real(j01),real(jsc1),real(vmax),&
			 real(voc),real(ff*100)

	!absc= 0.0


		if (pce .gt. pceaux) then
			pceaux = pce

			faux1 = thickness
			faux2 = fr
			faux3 = jsc1
			faux4 = j01
			faux5 = j1
			faux6 = voc
			faux7 = pce
			faux8 = vmax
			faux9 = ff			

		else
			continue
		end if
		
	!write(*,*) thickness,pce
	end do


	write(*,*) '#Jsc:',faux3,'W/Vm^{2}'
	write(*,*) '#Vmax:',faux8,'V'
	write(*,*) '#Voc:',faux6,'V'	
	write(*,*) '#fr:',faux2
	write(*,*) "#SLME-max",faux7*100.00,"%"
	write(*,*) '#Thickness:',faux1*1E6,'micro m'
	write(*,*) "#FF",faux9*100.00,"%"	


	deallocate(abscoef)

	close(100)
	close(200)




end program main




