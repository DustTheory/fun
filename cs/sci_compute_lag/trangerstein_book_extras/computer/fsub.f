        subroutine sub(m,n,vector,matrix)
        integer m,n
        double precision matrix(m,n),vector(n)
        integer i,j
        
        print *, "in sub"
        print *, "vector = ",(vector(j),j=1,n)
        do i=1,m
          print *, "matrix = ",(matrix(i,j),j=1,n)
        enddo

        return
        end
