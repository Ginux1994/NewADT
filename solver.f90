	module solver
	implicit none
	integer:: N
	real, allocatable:: Bmat(:,:), g(:)
	
contains
	subroutine JacobiSolver(A, b, res, N)
	implicit none 
	real:: A(N,N), b(N), res(N), Xpre(N), Xnow(N), tole, temp, normX
    integer, intent(in)::N
	integer:: i, j, k, l, iteMax, iteLoop
	allocate(Bmat(N,N), g(N))
	do i=1,N
		do j=1,N
            if(i==j) then
                Bmat(i,j) = 0
            else 
			    Bmat(i,j) = -A(i,j)/A(i,i)
            end if
		enddo
		g(i) = b(i)/A(i,i)
	enddo
	
	iteMax = 100
	normX = 0.0
    tole = 1e-7
	do iteLoop=1,iteMax	
        normX = 0.0
		do i=1,N
			temp = 0.0
			do j=1,N
				temp = temp + Bmat(i,j)*Xpre(j)
			enddo
			Xnow(i) = temp + g(i)
			normX = normX + sqrt((Xnow(i) - Xpre(i))*(Xnow(i) - Xpre(i)))
		enddo
		
		if(normX<tole) then
			res = Xnow
			exit
		else 
			Xpre = Xnow
		end if
	enddo
	
	
	
	end subroutine
	
	
	
	end module