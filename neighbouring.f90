Subroutine neighbour(coordinates,matrix,L,Ns,output)

integer,intent(in):: L,Ns


integer,intent(in),dimension(0:Ns-1,0:2)::coordinates
double precision,intent(in),dimension(0:4*L-1,0:4*L-1,0:4*L-1)::matrix
double precision,intent(out),dimension(0:Ns-1,0:5) ::output
double precision Ene,Mag,energy,j1,j2,j3,j4
double precision,dimension(0:2):: nb,n_spin,magn,s_rotate,s,m_t1b,m_t2
double precision,dimension(0:1):: m_e
double precision E1,M1,E2,M2,ma2,me,mt1a,mt1b,mt2,metheta
integer n,s1,t
integer,dimension(0:2) :: a1,a2,a3
double precision,dimension(0:2,0:2) :: j01,j02,j03,j12,j13,j23




!*************************************************CALCULATING
!ENERGY***************************************************************** 
do n1 =0,Ns-1
               s1 = mod(n1,4)


               if (s1==0) then 
                    a1 =coordinates(n1,:)-(coordinates(n1+1,:)-coordinates(n1,:))
                    a2 =coordinates(n1,:)-(coordinates(n1+2,:)-coordinates(n1,:))
                    a3 =coordinates(n1,:)-(coordinates(n1+3,:)-coordinates(n1,:))
                    do i2=0,2
                         a1(i2) = modulo(a1(i2) + 4*L, 4*L)
                         a2(i2) = modulo(a2(i2) + 4*L, 4*L)
                         a3(i2) = modulo(a3(i2) + 4*L, 4*L)
                    end do
					
                    output(n1,0) = n1+1
                    output(n1,1) = n1+2
				            output(n1,2) = n1+3
					          output(n1,3) = int(matrix(a1(0),a1(1),a1(2)))
					          output(n1,4) = int(matrix(a2(0),a2(1),a2(2)))
					          output(n1,5) = int(matrix(a3(0),a3(1),a3(2)))
					

               else if (s1==1) then 
                    a1 =coordinates(n1,:)-(coordinates(n1-1,:)-coordinates(n1,:))
                    a2 =coordinates(n1,:)-(coordinates(n1+1,:)-coordinates(n1,:))
                    a3 =coordinates(n1,:)-(coordinates(n1+2,:)-coordinates(n1,:))
                    do i2=0,2
                         a1(i2) = modulo(a1(i2) + 4*L, 4*L)
                         a2(i2) = modulo(a2(i2) + 4*L, 4*L)
                         a3(i2) = modulo(a3(i2) + 4*L, 4*L)
                    end do
					          output(n1,0) = n1-1
					          output(n1,1) = n1+1
					          output(n1,2) = n1+2
					          output(n1,3) = int(matrix(a1(0),a1(1),a1(2)))
					          output(n1,4) = int(matrix(a2(0),a2(1),a2(2)))
					          output(n1,5) = int(matrix(a3(0),a3(1),a3(2)))


               else if (s1==2) then 
                    a1 =coordinates(n1,:)-(coordinates(n1-2,:)-coordinates(n1,:))
                    a2 =coordinates(n1,:)-(coordinates(n1-1,:)-coordinates(n1,:))
                    a3 =coordinates(n1,:)-(coordinates(n1+1,:)-coordinates(n1,:))
                    do i2=0,2
                         a1(i2) = modulo(a1(i2) + 4*L, 4*L)
                         a2(i2) = modulo(a2(i2) + 4*L, 4*L)
                         a3(i2) = modulo(a3(i2) + 4*L, 4*L)
                    end do
					          output(n1,0) = n1-2
				            output(n1,1) = n1-1
					          output(n1,2) = n1+1
				            output(n1,3) = int(matrix(a1(0),a1(1),a1(2)))
					          output(n1,4) = int(matrix(a2(0),a2(1),a2(2)))
				            output(n1,5) = int(matrix(a3(0),a3(1),a3(2)))



               else 
                    a1 =coordinates(n1,:)-(coordinates(n1-3,:)-coordinates(n1,:))
                    a2 =coordinates(n1,:)-(coordinates(n1-2,:)-coordinates(n1,:))
                    a3 =coordinates(n1,:)-(coordinates(n1-1,:)-coordinates(n1,:))
                    do i2=0,2
                         a1(i2) = modulo(a1(i2) + 4*L, 4*L)
                         a2(i2) = modulo(a2(i2) + 4*L, 4*L)
                         a3(i2) = modulo(a3(i2) + 4*L, 4*L)
                    end do
					          output(n1,0) = n1-3
					          output(n1,1) = n1-2
					          output(n1,2) = n1-1
					          output(n1,3) = int(matrix(a1(0),a1(1),a1(2)))
					          output(n1,4) = int(matrix(a2(0),a2(1),a2(2)))
					          output(n1,5) = int(matrix(a3(0),a3(1),a3(2)))
               end if
               
end do



end Subroutine neighbour
