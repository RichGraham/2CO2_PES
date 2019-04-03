! Test program
use PES_2CO2_details
implicit none
 
double precision rab(9), PES_2CO2, xStar(9)

xStar(1:9)=(/0.20971571,  0.22605512,  0.22877215,  0.20158676,  0.2298365 , &
     0.2503203 ,  0.18510657,  0.2182286 ,  0.25385964/)
!xStar(1:9)=(/ 0.3018619 ,  0.33131717,  0.31699997,  0.30609471,  0.34708315, &
!     0.34030303,  0.27652502,  0.3131383 ,  0.31583372 /)

!xStar(1:9)=(/ 0.16918005,  0.20963065,  0.27441895,  0.15582496,  0.18869334, &
!        0.23773219,  0.14123042,  0.16638506,  0.20074532 /)

!xStar(1:9)=(/ 0.37918623,  0.43618246,  0.39957494,  0.32684715,  0.42808256, &
!     0.48455317,  0.26286132,  0.34577986,  0.44897025 /)

!xStar(1:9)=(/0.3791862300000000129962530, 0.4361824599999999940713735, &
!     0.3995749399999999895705116, 0.3268471499999999752006374, 0.4280825600000000008549250,&
!     0.4845531699999999775130277, 0.2628613200000000094114228, 0.3457798599999999944465401, &
!     0.4489702500000000151558766/)

!xStar(1:9)=(/ 0.24927635,  0.29087296,  0.3123213 ,  0.277394  ,  0.30196993, &
!        0.29306223,  0.28166435,  0.27929086,  0.25204695 /)

call load_GP_2CO2_Data
call fixedAngleSlice

!e=PES_2CO2_GP( xStar)
!e=PES_2CO2( rab)
!write(6,*)e

end
!

subroutine fixedAngleSlice()
  use PES_2CO2_details
#ifdef DEBUG
    use testing_2CO2
#endif

  implicit none
  double precision rab(9)
  integer i, itot
  double precision  r, beta1,beta2, alpha2, e, PES_2CO2, PI

  PI=4.D0*DATAN(1.D0)
    
  itot=500

  !T shape
  !beta1 = acos(0.0)
  !beta2 = acos(1.0)
  !alpha2 = 0.0

  !I shape
  !beta1 = acos(1.0)
  !beta2 = acos(1.0)
  !alpha2 = 0.0
  
  !Para shape
  !beta1 = acos(0.0)
  !beta2 = acos(0.0)
  !alpha2 = 0.0
  
  !X shape
  beta1 = acos(0.0)
  beta2 = acos(0.0)
  alpha2 = 0.5 * PI

  
  open (unit=15, file="PES_2CO2_Out.dat ", status='replace')
  
  do i=0, itot

     ! specify centre-to-centre separation
     r = (  1.5 + 20.0*i/(1.0*itot) ) 

     call computeDistances(r,alpha2,beta1,beta2,rab)

#ifdef DEBUG
     call testPerms( rab )
#endif
     
     e=PES_2CO2( rab)
     !e_GP = PES_2CO2_GP( xStar)
     write(15,*) r , e 
     
  enddo
  write(6,*)'Written to file: PES_2CO2_Out.dat '
  close(15)

end subroutine fixedAngleSlice
  
  
  
subroutine computeDistances(r,alpha2,beta1, beta2, rDist)
  use PES_2CO2_details
  implicit none
  double precision  rDist(9), r, beta1, beta2, alpha2
  double precision O1(3),O2(3),O3(3),O4(3),C1(3),C2(3)
  double precision sqLength

  O1(:)=0
  O2(:)=0
  O3(:)=0
  O4(:)=0

  C1(:)=0
  C2(:)=0
  
  !CO2 1
  !C is already at the origin
  O1(1)=lCO * sin( beta1)
  O1(3)=lCO *cos(beta1)
  
  O2(1)=-lCO * sin( beta1)
  O2(3)=-lCO * cos(beta1)
  
 
  !CO2 2
  C2(3)= r
  
  O3(1) =  lCO * sin(beta2) * cos(alpha2)
  O3(2) =  lCO * sin(beta2) * sin(alpha2)
  O3(3) =  r  +  lCO * cos(beta2)
  
  O4(1) =  -lCO * sin(beta2) * cos(alpha2)
  O4(2) =  -lCO * sin(beta2) * sin(alpha2)
  O4(3) =  r  -  lCO * cos(beta2)
  
  
  rDist(1)= SQRT( sqLength(O1,O3) )
  rDist(2)= SQRT( sqLength(O1,C2) )
  rDist(3)= SQRT( sqLength(O1,O4) )


  rDist(4)= SQRT( sqLength(C1,O3) )
  rDist(5)= SQRT( sqLength(C1,C2) )
  rDist(6)= SQRT( sqLength(C1,O4) )

  rDist(7)= SQRT( sqLength(O2,O3) )
  rDist(8)= SQRT( sqLength(O2,C2) )
  rDist(9)= SQRT( sqLength(O2,O4) )

  
   
end subroutine computeDistances


double precision function sqLength( u,v)
  !returns the square of the distance between two vectors
  implicit none
  double precision sum
  double precision u(3), v(3)
  integer i

  sum=0.0
  do i = 1,3 
    sum = sum  +  (  v(i)-u(i)  )**2
  enddo
  
  sqLength= sum
end function sqLength
