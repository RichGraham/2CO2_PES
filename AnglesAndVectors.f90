module AnglesAndVectors

contains

  
subroutine computeDistances(r12,r13,r23,alpha1,beta1, alpha2, beta2,rab)
  use PES_details
  use AnglesAndVectors
  implicit none
  
  double precision, dimension(:)::  rab
  double precision  r12,r13,r23, beta1, alpha1, beta2, alpha2
  !double precision  O1(3),O2(3),C1(3),O3(3),C2(3),O4(3),Ar(3), cosTheta13
  double precision :: O1(3)=0.0,O2(3)=0.0,C1(3)=0.0,O3(3)=0.0,C2(3)=0.0,O4(3)=0.0,Ar(3)=0.0, cosTheta13

  !==== CO2 1 ====
  !C is already at the origin
  O1(1)=lCO * sin( beta1) * cos(alpha1);
  O1(2)=lCO * sin( beta1) * sin(alpha1);
  O1(3)=lCO *cos(beta1);

  O2(1)=-lCO * sin( beta1)* cos(alpha1);
  O2(2)=-lCO * sin( beta1) * sin(alpha1);
  O2(3)=-lCO * cos(beta1);

  !print *,O1
  !print *,O2
  
  !==== CO2 2 ====
  C2(3)= r12;

  O3(1)=lCO * sin( beta2) * cos(alpha2);
  O3(2)=lCO * sin( beta2) * sin(alpha2);
  O3(3)=lCO *cos(beta2)  + r12;
  
  O4(1)=-lCO * sin( beta2)* cos(alpha2);
  O4(2)=-lCO * sin( beta2) * sin(alpha2);
  O4(3)=-lCO * cos(beta2) + r12;  

  !==== Ar ====
  cosTheta13 = (r12*r12 + r13*r13 - r23*r23) / 2.0/r12/r13;
  Ar(1) = r13 * sqrt( abs(1.0 - cosTheta13*cosTheta13));
  Ar(3) = r13 * cosTheta13;

  !compute interatomic distances
  rab(1) = sqrt( sqLength(O1,O3)  );
  rab(2) = sqrt( sqLength(O1,C2)  );
  rab(3) = sqrt( sqLength(O1,O4)  );
  !=
  rab(4) = sqrt( sqLength(C1,O3)  );
  rab(5) = sqrt( sqLength(C1,C2)  );
  rab(6) = sqrt( sqLength(C1,O4)  );
  !=
  rab(7) = sqrt( sqLength(O2,O3)  );
  rab(8) = sqrt( sqLength(O2,C2)  );
  rab(9) = sqrt( sqLength(O2,O4)  );
  
  rab(10) = sqrt(  sqLength(Ar,O1)  );
  rab(11) = sqrt(  sqLength(Ar,C1)  );
  rab(12) = sqrt(  sqLength(Ar,O2)  );
  
  rab(13) = sqrt(  sqLength(Ar,O3)  );
  rab(14) = sqrt(  sqLength(Ar,C2)  );
  rab(15) = sqrt(  sqLength(Ar,O4)  );

  !print *, rab(
  
end subroutine computeDistances



subroutine computeDistancesWithSlide(r12,r13,r23,alpha1,beta1,cosThetaTilde,deltaX, rab)
  use PES_details
  use AnglesAndVectors
  implicit none
  double precision, dimension(:)::  rab
  double precision  r12,r13,r23, beta1, alpha1,cosThetaTilde, deltaX
  double precision O1(3),O2(3),C(3),Ar1(3),Ar2(3), cosTheta13

  C(:)=0.0
  O1(:)=0.0
  O2(:)=0.0
  Ar1(:)=0.0
  Ar2(:)=0.0

  !C begin at the origin but slides along ThetaTilde a distance deltaX
  C(1)=deltaX*Sqrt(1.0- cosThetaTilde**2)
  C(3)=deltaX*cosThetaTilde
  
  
  O1(1)=C(1) + lCO * sin( beta1) * cos(alpha1);
  O1(2)=lCO * sin( beta1) * sin(alpha1);
  O1(3)=C(3) + lCO *cos(beta1);
  
  O2(1)=C(1) -lCO * sin( beta1)* cos(alpha1);
  O2(2)=-lCO * sin( beta1) * sin(alpha1);
  O2(3)=C(3) -lCO * cos(beta1);
  
  
  !Ar1
  Ar1(3)= r12;


  !Ar2
  cosTheta13 = (r12*r12 + r13*r13 - r23*r23) / 2.0/r12/r13;
  Ar2(1) = r13 * sqrt( abs(1.0 - cosTheta13*cosTheta13));
  Ar2(3) = r13 * cosTheta13;
    
  !print*,'*********************',Ar2, cosTheta13

  !compute interatomic distances
  rab(1) = sqrt( sqLength(Ar1,O1)  );
  rab(2) = sqrt( sqLength(Ar1,C)  );
  rab(3) = sqrt(  sqLength(Ar1,O2)  );
  
  rab(4) = sqrt(  sqLength(Ar2,O1)  );
  rab(5) = sqrt(  sqLength(Ar2,C)  );
  rab(6) = sqrt(  sqLength(Ar2,O2)  );
  
  rab(7) = sqrt(  sqLength(Ar2,Ar1)  );
  
end subroutine computeDistancesWithSlide

