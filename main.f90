program test     
    use solCtrlInf
    use mainCtrlInf
    use timeStepInf
    use nodeInf
    use ndSHInf
    use initInf
    use heatFlowCtrlInf
    use timeFuncInf
    use heatFlowInf
    use element

    use solver
    implicit none
    
    real:: res(100), mat(100,100), r(100)
    integer:: i,j,k
    k=0
    do i=1,100
        do j=1,100
            k = k+1
            mat(i,j) =k
        enddo
        r(i) = k
    enddo
    call JacobiSolver(mat,r,res,100)
! read in the information, no calculation(ind=0)    
    iper = 1; ind = 0
    call ADINI
    
	ind = 1; isref = 0; index=0; nReformStep = 0;
    call assem_element
    
	call execute_heatFlow


    

    
    
    
    
    
    

    
end program
 