Subroutine structurefactor(config,coordinates,Ns,q,output)

integer,intent(in):: Ns

!double precision,intent(in),dimension(0:k_size-1):: h,k
double precision,intent(in),dimension(0:Ns-1,0:2)::config
integer,intent(in),dimension(0:Ns-1,0:2)::coordinates

        
double precision,intent(in),dimension(0:2):: q
double precision,dimension(0:2):: spin_pl

double precision ::sq
double precision,intent(out) ::output

double precision,dimension(0:2) :: x,y

Pi = 3.141592653589793




        x = (/0.0,0.0,0.0 /)
        y = (/0.0,0.0,0.0 /)
!        q = ((2.0*Pi)/(4.0*L))*(/h(i),h(i),k(j)/)
        do ii=0,Ns-1
             
            spin_pl = config(ii,:) - (DOT_PRODUCT(q,config(ii,:))/norm2(q)**2)*q
            spin_pl = spin_pl/(2.0*norm2(spin_pl))
            x = x + spin_pl*COS(DOT_PRODUCT(q,coordinates(ii,:)))
            y =  y + spin_pl*SIN(DOT_PRODUCT(q,coordinates(ii,:)))

        end do
        sq =  (norm2(x)**2 + norm2(y)**2)/Ns
  

output = sq
end Subroutine structurefactor
