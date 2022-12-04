Subroutine interaction(neigh,J_combined,L,Ns,output)

integer,intent(in):: L,Ns

double precision,intent(in),dimension(0:Ns-1,0:5)::neigh
double precision,intent(in),dimension(0:Ns-1,0:5,0:3):: J_combined
double precision,intent(out),dimension(0:Ns-1,0:5,0:2,0:2) ::output
double precision j1,j2,j3,j4
double precision,dimension(0:2,0:2) :: j01,j02,j03,j12,j13,j23
double precision,dimension(0:2,0:2) :: j0_1, j0_2,j0_3,j1_2,j1_3,j2_3

     do n1 =0,Ns-1
          s1 = mod(n1,4)

          if (s1==0) then 

               j1 = J_combined(n1,0,0)
               j2 = J_combined(n1,0,1)
               j3 = J_combined(n1,0,2)
               j4 = J_combined(n1,0,3)               
               j01 =transpose(reshape((/j1,j3,-j4,j3,j1,-j4,j4,j4,j2/),shape(j01)))
          
          
               j1 = J_combined(n1,1,0)
               j2 = J_combined(n1,1,1)
               j3 = J_combined(n1,1,2)
               j4 = J_combined(n1,1,3)               
               j02 = transpose(reshape((/j2,j4,j4,-j4,j1,j3,-j4,j3,j1/),shape(j02)))
          
          
               j1 = J_combined(n1,2,0)
               j2 = J_combined(n1,2,1)
               j3 = J_combined(n1,2,2)
               j4 = J_combined(n1,2,3)               
               j03 = transpose(reshape((/j1,-j4,j3,j4,j2,j4,j3,-j4,j1/),shape(j03)))
               
               
               j1 = J_combined(n1,3,0)
               j2 = J_combined(n1,3,1)
               j3 = J_combined(n1,3,2)
               j4 = J_combined(n1,3,3)               
               j0_1 =transpose(reshape((/j1,j3,-j4,j3,j1,-j4,j4,j4,j2/),shape(j01)))
          
          
               j1 = J_combined(n1,4,0)
               j2 = J_combined(n1,4,1)
               j3 = J_combined(n1,4,2)
               j4 = J_combined(n1,4,3)               
               j0_2 = transpose(reshape((/j2,j4,j4,-j4,j1,j3,-j4,j3,j1/),shape(j02)))
          
          
               j1 = J_combined(n1,5,0)
               j2 = J_combined(n1,5,1)
               j3 = J_combined(n1,5,2)
               j4 = J_combined(n1,5,3)               
               j0_3 = transpose(reshape((/j1,-j4,j3,j4,j2,j4,j3,-j4,j1/),shape(j03)))                        
    
               output(n1,0,:,:) = j01
               output(n1,1,:,:) = j02
               output(n1,2,:,:) = j03
               output(n1,3,:,:) = j0_1
               output(n1,4,:,:) = j0_2
               output(n1,5,:,:) = j0_3
                              

						
						
          else if (s1==1) then 


			   
               j1 = J_combined(n1,0,0)
               j2 = J_combined(n1,0,1)
               j3 = J_combined(n1,0,2)
               j4 = J_combined(n1,0,3)               
               j01 =transpose(reshape((/j1,j3,-j4,j3,j1,-j4,j4,j4,j2/),shape(j01)))
          
          
               j1 = J_combined(n1,1,0)
               j2 = J_combined(n1,1,1)
               j3 = J_combined(n1,1,2)
               j4 = J_combined(n1,1,3)               
               j12 = transpose(reshape((/j1,-j4,-j3,j4,j2,-j4,-j3,j4,j1/),shape(j12)))
          
          
               j1 = J_combined(n1,2,0)
               j2 = J_combined(n1,2,1)
               j3 = J_combined(n1,2,2)
               j4 = J_combined(n1,2,3)               
               j13 = transpose(reshape((/j2,j4,-j4,-j4,j1,-j3,j4,-j3,j1/),shape(j13)))    
               
               j1 = J_combined(n1,3,0)
               j2 = J_combined(n1,3,1)
               j3 = J_combined(n1,3,2)
               j4 = J_combined(n1,3,3)               
               j0_1 =transpose(reshape((/j1,j3,-j4,j3,j1,-j4,j4,j4,j2/),shape(j01)))
          
          
               j1 = J_combined(n1,4,0)
               j2 = J_combined(n1,4,1)
               j3 = J_combined(n1,4,2)
               j4 = J_combined(n1,4,3)               
               j1_2 = transpose(reshape((/j1,-j4,-j3,j4,j2,-j4,-j3,j4,j1/),shape(j12)))
          
          
               j1 = J_combined(n1,5,0)
               j2 = J_combined(n1,5,1)
               j3 = J_combined(n1,5,2)
               j4 = J_combined(n1,5,3)               
               j1_3 = transpose(reshape((/j2,j4,-j4,-j4,j1,-j3,j4,-j3,j1/),shape(j13)))                         			   
               output(n1,0,:,:) = transpose(j01)
               output(n1,1,:,:) = j12
               output(n1,2,:,:) = j13
               output(n1,3,:,:) = transpose(j0_1)
               output(n1,4,:,:) = j1_2
               output(n1,5,:,:) = j1_3
                


          else if (s1==2) then 
 
			   
               j1 = J_combined(n1,0,0)
               j2 = J_combined(n1,0,1)
               j3 = J_combined(n1,0,2)
               j4 = J_combined(n1,0,3)               
               j02 = transpose(reshape((/j2,j4,j4,-j4,j1,j3,-j4,j3,j1/),shape(j02)))
          
          
               j1 = J_combined(n1,1,0)
               j2 = J_combined(n1,1,1)
               j3 = J_combined(n1,1,2)
               j4 = J_combined(n1,1,3)               
               j12 = transpose(reshape((/j1,-j4,-j3,j4,j2,-j4,-j3,j4,j1/),shape(j12)))
          
          
               j1 = J_combined(n1,2,0)
               j2 = J_combined(n1,2,1)
               j3 = J_combined(n1,2,2)
               j4 = J_combined(n1,2,3)               
               j23 = transpose(reshape((/j1,-j3,j4,-j3,j1,-j4,-j4,j4,j2/),shape(j23)))               			   
			   
			   
               j1 = J_combined(n1,3,0)
               j2 = J_combined(n1,3,1)
               j3 = J_combined(n1,3,2)
               j4 = J_combined(n1,3,3)               
               j0_2 = transpose(reshape((/j2,j4,j4,-j4,j1,j3,-j4,j3,j1/),shape(j02)))
          
          
               j1 = J_combined(n1,4,0)
               j2 = J_combined(n1,4,1)
               j3 = J_combined(n1,4,2)
               j4 = J_combined(n1,4,3)               
               j1_2 = transpose(reshape((/j1,-j4,-j3,j4,j2,-j4,-j3,j4,j1/),shape(j12)))
          
          
               j1 = J_combined(n1,5,0)
               j2 = J_combined(n1,5,1)
               j3 = J_combined(n1,5,2)
               j4 = J_combined(n1,5,3)               
               j2_3 = transpose(reshape((/j1,-j3,j4,-j3,j1,-j4,-j4,j4,j2/),shape(j23)))  
    

               output(n1,0,:,:) = transpose(j02)
               output(n1,1,:,:) = transpose(j12)
               output(n1,2,:,:) = j23
               output(n1,3,:,:) = transpose(j0_2)
               output(n1,4,:,:) = transpose(j1_2)
               output(n1,5,:,:) = j2_3


          else 
		   
               j1 = J_combined(n1,0,0)
               j2 = J_combined(n1,0,1)
               j3 = J_combined(n1,0,2)
               j4 = J_combined(n1,0,3)               
               j03 = transpose(reshape((/j1,-j4,j3,j4,j2,j4,j3,-j4,j1/),shape(j03)))
          
          
               j1 = J_combined(n1,1,0)
               j2 = J_combined(n1,1,1)
               j3 = J_combined(n1,1,2)
               j4 = J_combined(n1,1,3)               
               j13 = transpose(reshape((/j2,j4,-j4,-j4,j1,-j3,j4,-j3,j1/),shape(j13)))
          
          
               j1 = J_combined(n1,2,0)
               j2 = J_combined(n1,2,1)
               j3 = J_combined(n1,2,2)
               j4 = J_combined(n1,2,3)               
               j23 = transpose(reshape((/j1,-j3,j4,-j3,j1,-j4,-j4,j4,j2/),shape(j23))) 
               
               
               j1 = J_combined(n1,3,0)
               j2 = J_combined(n1,3,1)
               j3 = J_combined(n1,3,2)
               j4 = J_combined(n1,3,3)               
               j0_3 = transpose(reshape((/j1,-j4,j3,j4,j2,j4,j3,-j4,j1/),shape(j03)))
          
          
               j1 = J_combined(n1,4,0)
               j2 = J_combined(n1,4,1)
               j3 = J_combined(n1,4,2)
               j4 = J_combined(n1,4,3)               
               j1_3 = transpose(reshape((/j2,j4,-j4,-j4,j1,-j3,j4,-j3,j1/),shape(j13)))
          
          
               j1 = J_combined(n1,5,0)
               j2 = J_combined(n1,5,1)
               j3 = J_combined(n1,5,2)
               j4 = J_combined(n1,5,3)               
               j2_3 = transpose(reshape((/j1,-j3,j4,-j3,j1,-j4,-j4,j4,j2/),shape(j23))) 


               output(n1,0,:,:) = transpose(j03)
               output(n1,1,:,:) = transpose(j13)
               output(n1,2,:,:) = transpose(j23)
               output(n1,3,:,:) = transpose(j0_3)
               output(n1,4,:,:) = transpose(j1_3)
               output(n1,5,:,:) = transpose(j2_3)

          end if

     end do



end Subroutine interaction
