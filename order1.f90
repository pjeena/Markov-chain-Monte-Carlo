Subroutine calcenergy(config,coordinates,matrix,neigh,Jij_matrix,L,Ns,output)

integer,intent(in):: L,Ns

double precision,intent(in),dimension(0:Ns-1,0:2)::config
integer,intent(in),dimension(0:Ns-1,0:2)::coordinates
double precision,intent(in),dimension(0:4*L-1,0:4*L-1,0:4*L-1)::matrix
double precision,intent(in),dimension(0:Ns-1,0:5)::neigh
double precision,intent(in),dimension(0:Ns-1,0:5,0:2,0:2):: Jij_matrix
double precision,intent(out),dimension(0:6) ::output
double precision Ene,Mag,energy,j1,j2,j3,j4
double precision,dimension(0:2):: nb,n_spin,magn,s_rotate,s,m_t1b,m_t2
double precision,dimension(0:1):: m_e
double precision E1,M1,E2,M2,ma2,me,mt1a,mt1b,mt2,metheta
integer n,s1,t
integer,dimension(0:2) :: a1,a2,a3
double precision,dimension(0:2,0:2) :: j01,j02,j03,j12,j13,j23
double precision,dimension(0:2,0:2) :: j0_1, j0_2,j0_3,j1_2,j1_3,j2_3



E1=0.0
M1=0.0
E2=0.0
M2=0.0


!*************************************************CALCULATING
!ENERGY***************************************************************** 
     energy = 0.0
     do n1 =0,Ns-1
          s = config(n1,:)
          s1 = mod(n1,4)

          n_spin =matmul(Jij_matrix(n1,0,:,:),config(int(neigh(n1,0)),:))+&
                        matmul(Jij_matrix(n1,1,:,:),config(int(neigh(n1,1)),:))+&
                        matmul(Jij_matrix(n1,2,:,:),config(int(neigh(n1,2)),:))+&
                        matmul(Jij_matrix(n1,3,:,:),config(int(neigh(n1,3)),:))+&   
                        matmul(Jij_matrix(n1,4,:,:),config(int(neigh(n1,4)),:))+&
                        matmul(Jij_matrix(n1,5,:,:),config(int(neigh(n1,5)),:))



          energy = energy + dot_product(n_spin,s)
     end do
     Ene = energy/2


     magn = (/0.0,0.0,0.0/)
     do n2=0,Ns-1
        magn = magn + config(n2,:)

     end do
     Mag = 0.5*norm2(magn)

!***************************************m_A2*****************************************
     ma2 = 0.0
     do i=0,Ns-1,4
           ma2 = ma2 + config(i,0) + config(i,1)+ config(i,2) &
                         +config(i+2,0) - config(i+2,1)- config(i+2,2) &
                         -config(i+3,0) + config(i+3,1)- config(i+3,2) &
                         -config(i+1,0) - config(i+1,1)+ config(i+1,2) 

     end do
     ma2 = ABS(ma2)/(2.0*sqrt(3.0))
     

!***************************************m_e*****************************************
     m_e = (/0.0,0.0/)
     do i=0,Ns-1,4
           m_e(0) = m_e(0) -2.0*config(i,0) + config(i,1)+ config(i,2) &
                        -2.0*config(i+2,0) - config(i+2,1)- config(i+2,2) &
                        +2.0*config(i+3,0) + config(i+3,1)- config(i+3,2) &
                        +2.0*config(i+1,0) - config(i+1,1)+ config(i+1,2)
                      
           m_e(1) = m_e(1) -config(i,1) + config(i,2)+ config(i+2,1)-config(i+2,2) &
                         -config(i+3,1)- config(i+3,2) +config(i+1,1)+config(i+1,2) 
     end do

     m_e(0) =  m_e(0)/(2.0*sqrt(6.0))
     m_e(1) = m_e(1)/(2.0*sqrt(2.0))
     me = norm2(m_e)


!***************************************m_t1b*****************************************

     m_t1b = (/0.0,0.0,0.0/)
     do i=0,Ns-1,4
           m_t1b(0) = m_t1b(0) + config(i,1) + config(i,2)- config(i+2,1)-config(i+2,2) &
                      - config(i+3,1)+ config(i+3,2) + config(i+1,1)-config(i+1,2)

           m_t1b(1) = m_t1b(1) + config(i,0) + config(i,2)- config(i+2,0)+config(i+2,2) &
                      - config(i+3,0)- config(i+3,2) + config(i+1,0)-config(i+1,2)

           m_t1b(2) = m_t1b(2) + config(i,0) + config(i,1)- config(i+2,0)+config(i+2,1) &
                      + config(i+3,0)- config(i+3,1) - config(i+1,0)-config(i+1,1)

     end do
     mt1b = norm2(m_t1b) /(2.0*sqrt(2.0))                


!***************************************m_t2*****************************************
     m_t2 = (/0.0,0.0,0.0/)             
     do i=0,Ns-1,4
           m_t2(0) = m_t2(0) - config(i,1) + config(i,2)+ config(i+2,1)-config(i+2,2) &
                      + config(i+3,1)+ config(i+3,2) - config(i+1,1)-config(i+1,2)

           m_t2(1) = m_t2(1) + config(i,0) - config(i,2)- config(i+2,0)-config(i+2,2) &
                      - config(i+3,0)+ config(i+3,2) + config(i+1,0)+config(i+1,2)

           m_t2(2) = m_t2(2) - config(i,0) + config(i,1)+ config(i+2,0)+config(i+2,1) &
                      - config(i+3,0)- config(i+3,1) + config(i+1,0)-config(i+1,1)

     end do
     mt2 = norm2(m_t2) /(2.0*sqrt(2.0))            
     

     metheta = me * COS(6.0*DATAN(m_e(1)/m_e(0)))







output(0) = Ene
output(1) = Mag
output(2) = ma2
output(3) = me
output(4) = mt1b
output(5) = mt2
output(6) = metheta

end Subroutine calcenergy
