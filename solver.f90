    module choSolver
    !----------------------------------------module coment
    !  Version     :  V1.0
    !  Coded by    :  syz
    !  Date        :
    !-----------------------------------------------------
    !  Description :   Cholesky分解计算对称正定方程模块
    !
    !-----------------------------------------------------

    contains

    subroutine solve(A,b,x,N)
    implicit real(a-z)
    integer::N
    real::A(N,N),b(N),x(N)
    real::L(N,N),y(N),LT(N,N)
    !LT 为L的转置矩阵
    integer::i,j

    call  chol(A,L,N)

    call  downtri(L,b,y,N)

    do i=1,N
        do j=1,N
            LT(i,j)=L(j,i)
        end do
    end do

    call uptri(LT,y,x,N)  !这一步已经算出了x

    end subroutine solve

    subroutine  chol(A,L,N)
    !---------------------------------subroutine  comment
    !  Version   :  V1.0
    !  Coded by  :  syz
    !  Date      :
    !-----------------------------------------------------
    !  Purpose   :  Cholesky分解子程序
    !-----------------------------------------------------
    integer::N
    real::A(N,N),L(N,N)
    integer::i,j,k

    L=0

    L(1,1)=sqrt(a(1,1))


    L(2:,1)=a(2:,1)/L(1,1)


    do j=2,N

        s=0
        do k=1,j-1
            s=s+L(j,k)**2
        end do

        L(j,j)=sqrt(a(j,j)-s)

        !注意i范围
        do i=j+1,N

            s=0
            do k=1,j-1
                s=s+L(i,k)*L(j,k)
            end do

            L(i,j)=(a(i,j)-s)/L(j,j)

        end do

    end do

    end subroutine


    subroutine uptri(A,b,x,N)
    !---------------------------------subroutine  comment
    !  Version   :  V1.0
    !  Coded by  :  syz
    !  Date      :  2010-4-8
    !-----------------------------------------------------
    !  Purpose   :  上三角方程组的回带方法
    !                 Ax=b
    !-----------------------------------------------------

    implicit real(a-z)

    integer::i,j,k,N

    real::A(N,N),b(N),x(N)

    x(N)=b(N)/A(N,N)

    !回带部分
    do i=n-1,1,-1

        x(i)=b(i)
        do j=i+1,N
            x(i)=x(i)-a(i,j)*x(j)
        end do
        x(i)=x(i)/A(i,i)

    end do

    end subroutine uptri


    subroutine downtri(A,b,x,N)
    !---------------------------------subroutine  comment
    !  Version   :  V1.0
    !  Coded by  :  syz
    !  Date      :  2010-4-9
    !-----------------------------------------------------
    !  Purpose   :  下三角方程组的回带方法
    !                 Ax=b
    !-----------------------------------------------------

    implicit real(a-z)
    integer::i,j,N
    real::A(N,N),b(N),x(N)

    x(1)=b(1)/a(1,1)

    do k=2,N
        x(k)=b(k)
        do i=1,k-1
            x(k)=x(k)-a(k,i)*x(i)
        end do
        x(k)=x(k)/a(k,k)

    end do

    end subroutine downtri

    end module choSolver

    
    
    module solver
	implicit none

contains
	subroutine JacobiSolver(A, b, res, N)
	implicit none 
	real:: Bmat(N,N), g(N), A(N,N), b(N), res(N), Xpre(N), Xnow(N), tole, temp, normX
    integer, intent(in)::N
	integer:: i, j, k, l, iteMax, iteLoop
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
	
	
    
    subroutine choleskySolver(A, b, res, N)
    implicit none 
	real:: A(N,N), L_(N,N), b(N), res(N), d, s
    integer, intent(in)::N
	integer:: i, j, k, l, iteMax, iteLoop
    logical:: isspd
    
    do j=1,N
        d=0.0
        do k=1,j-1
            s = 0.0
            do i=1,k
                s = s + L_(k,i)*L_(j,i)
            enddo
            s = (A(j,k) - s)/L_(k,k)
            L_(j,k) = s
            d = d + s*s
            isspd = isspd .and. (A(k,j) == A(j,k))
        enddo
        d = A(j,j) - d
        isspd = isspd .and. (d > 0.0);
        if(d>0.0) then
            L_(j,j) = sqrt(d)
        else 
            L_(j,j) = 0.0
        end if
        
        do k=j+1,N
            L_(j,k) = 0.0
        enddo
    enddo
    
    res = b
    do k=1,N
        do i=1,k-1
            res(k) = res(k) - res(k)*L_(k,i)
        enddo
        res(k) = res(k)/L_(k,k)
    enddo
    
    do k=N,1,-1
        do i=k+1,N
            res(k) = res(k) - res(k)*L_(i,k)
        enddo
        res(k) = res(k)/L_(k,k)
    enddo
    
    
    
    
    end subroutine
    
	end module