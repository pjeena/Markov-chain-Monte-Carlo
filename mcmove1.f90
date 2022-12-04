Subroutine montecarlo(config,coordinates,matrix,neigh,Jij_matrix,L,Ns,iT)
integer,intent(in):: L,Ns
double precision,intent(in):: iT
double precision,intent(inout),dimension(0:Ns-1,0:2)::config
integer,intent(in),dimension(0:Ns-1,0:2)::coordinates
double precision,intent(in),dimension(0:4*L-1,0:4*L-1,0:4*L-1)::matrix
double precision,intent(in),dimension(0:Ns-1,0:5)::neigh
double precision,intent(in),dimension(0:Ns-1,0:5,0:2,0:2):: Jij_matrix
!double precision,dimension(0:Ns-1,0:2)::config
!double precision,intent(out),dimension(0:3)::output
double precision Ene,Mag,energy,r1,r2,j1,j2,j3,j4
double precision,dimension(0:2):: nb,n_spin,magn
double precision u,v,cost,v1,v2,x,E1,M1,E2,M2,R,phi,Pi
double precision,dimension(0:2) :: s_rotate,s,a,w,p
integer n,s1,eqSteps,mcSteps,t
integer,dimension(0:2) :: a1,a2,a3
double precision,dimension(0:2,0:2) :: j01,j02,j03,j12,j13,j23,Rot
double precision,dimension(0:2,0:2) :: j0_1, j0_2,j0_3,j1_2,j1_3,j2_3




Pi = 3.1415927 




     do j=0,Ns-1
          call random_number(u)                           
          n = int(FLOOR(Ns*u))        
 

!          s = config(n,:)                        ! choosing spin

          s1 = mod(n,4)                     
          n_spin =matmul(Jij_matrix(n,0,:,:),config(int(neigh(n,0)),:))+&
                        matmul(Jij_matrix(n,1,:,:),config(int(neigh(n,1)),:))+&
                        matmul(Jij_matrix(n,2,:,:),config(int(neigh(n,2)),:))+&
                        matmul(Jij_matrix(n,3,:,:),config(int(neigh(n,3)),:))+&   
                        matmul(Jij_matrix(n,4,:,:),config(int(neigh(n,4)),:))+&
                        matmul(Jij_matrix(n,5,:,:),config(int(neigh(n,5)),:))

          

          nb = n_spin                  !nearest neighbour effective spin
          call random_number(R)
          x = norm2(n_spin)*iT/2.0

          theta = DACOS(-LOG(R*exp(x)+ (1.0-R)*exp(-x))/x)      
          
          call random_number(w)
          w = w/norm2(w)
          p = w - (dot_product(n_spin,w)/norm2(n_spin)**2)*n_spin
          p = p/norm2(p)
          
          n_spin = n_spin/norm2(n_spin)
          call random_number(phi)
          phi = 2.0*Pi*phi
          Rot(0,0) = COS(phi) + (1-COS(phi))*n_spin(0)**2
          Rot(0,1) = n_spin(0)*n_spin(1)*(1-COS(phi)) - n_spin(2)*SIN(phi)
          Rot(0,2) = n_spin(0)*n_spin(2) *(1-COS(phi)) + n_spin(1)*SIN(phi)
          Rot(1,0) = n_spin(1)*n_spin(0)*(1-COS(phi)) + n_spin(2)*SIN(phi)
          Rot(1,1) = COS(phi) + (1-COS(phi))*n_spin(1)**2
          Rot(1,2) = n_spin(1)*n_spin(2) *(1-COS(phi)) - n_spin(0)*SIN(phi)
          Rot(2,0) = n_spin(2)*n_spin(0) *(1-COS(phi)) - n_spin(1)*SIN(phi)
          Rot(2,1) = n_spin(2)*n_spin(1) *(1-COS(phi)) + n_spin(0)*SIN(phi)
          Rot(2,2) = COS(phi) + (1-COS(phi))*n_spin(2)**2

          p = matmul(Rot,p)
          s = (COS(theta)*n_spin/norm2(n_spin)) + (SIN(theta)*p/norm2(p))
          s = s/(2.0*norm2(s))

          config(n,:) = s                                !modified configuration
     end do

end Subroutine montecarlo
