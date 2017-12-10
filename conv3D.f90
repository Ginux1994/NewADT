
    
    
    
    
    
    subroutine dNdxConv(r, s, numNd, nod5, N, area, coords)
    ! subroutine param declare
    real, intent(in):: r, s, area, coords(3,4)
    integer, intent(in):: numNd, nod5
    real, intent(out)::N(4)
    
    ! param in the subroutine
    real:: dNdxi(2, 4)
    
    end subroutine 