   
    
    subroutine mapLocalDOF(numNd, ndStatus)
    use mainCtrlInf, only: numEquation
    use solCtrlInf, only:DOFMap, numMasterDOF, numSlaveDOF
    implicit none ! 变量类型声明和implicit none必须一起出现
    ! subr variable declartion
    integer, intent(in):: numNd
    integer, intent(in):: ndStatus(numNd,1)


    ! inner variable
    integer:: nodeLoop, DOFloop, nodDOFLoop = 1, ndDOF = 1, indMaster = 0, indSlave = 0

!    allocate (ndStatus(numNd, 1)) ! 动态数组大小声明必须在 变量声明之后
    allocate(DOFMap(numNd))
    do nodeLoop=1,numNd
        do DOFloop=1, ndDOF
            if(ndStatus(nodeLoop, DOFloop) == 0) then ! ignore this node dof
                DOFMap(nodDOFLoop) = 0
            else if(ndStatus(nodeLoop, DOFloop) == 1) then ! master node dof
                indMaster = indMaster + 1
                DOFMap(nodDOFLoop) = indMaster
            else if(ndStatus(nodeLoop, DOFloop) == 2) then ! slave node dof
                indSlave = indSlave + 1
                DOFMap(nodDOFLoop) = -indSlave
            end if
            nodDOFLoop = nodDOFLoop + 1
        enddo
    enddo 
    numMasterDOF = indMaster; numSlaveDOF = indSlave
    numEquation = numMasterDOF + numSlaveDOF
    end subroutine mapLocalDOF
   
 !************************************************************************************************************************       
    subroutine assemStaMatK_k(elmMat, ndID)
	
    use solCtrlInf
    use nodeInf, only: dofID
    implicit none
    real:: elmMat(8,8)
    integer:: ndID(8)
    integer:: outterDOFID, innerDOFID, i, j   
    
    do i =1,8
        outterDOFID = dofID(1, ndID(i))
        do j=1,8
            innerDOFID = dofID(1, ndID(j))
            staMatK_k(outterDOFID,innerDOFID) = staMatK_k(outterDOFID,innerDOFID) + elmMat(i,j)
        enddo
    enddo   
    end subroutine assemStaMatK_k
    
 !************************************************************************************************************************        
   subroutine assemStaMatK_c(elmMat, ndID)
	
    use solCtrlInf
    use nodeInf, only: dofID
    implicit none
    real:: elmMat(4,4)
    integer:: ndID(4)
    integer:: outterDOFID, innerDOFID, i, j   
    
    do i =1,4
        outterDOFID = dofID(1, ndID(i))
        do j=1,4
            innerDOFID = dofID(1, ndID(j))
            staMatK_c(outterDOFID,innerDOFID) = staMatK_c(outterDOFID,innerDOFID) + elmMat(i,j)
        enddo
    enddo   
    end subroutine assemStaMatK_c
 !************************************************************************************************************************        
       
    
    subroutine assemDynMatK_k(elmMat, ndID)
    use solCtrlInf
    use nodeInf, only: dofID
    implicit none
    real:: elmMat(8,8)
    integer:: ndID(8)   
    integer:: outterDOFID, innerDOFID, i, j   

    do i =1,8
        outterDOFID = dofID(1, ndID(i))
        do j=1,8
            innerDOFID = dofID(1, ndID(j))
            dynMat(outterDOFID,innerDOFID) = dynMat(outterDOFID,innerDOFID) + elmMat(i,j)
        enddo
    enddo   
    end subroutine assemDynMatK_k
	
 !************************************************************************************************************************        
        
    
    subroutine assemDynMatK_c_r(elmMat, ndID)
    use solCtrlInf
    use nodeInf, only: dofID
    implicit none
    real:: elmMat(4,4)
    integer:: ndID(4)   
    integer:: outterDOFID, innerDOFID, i, j    

    do i =1,4
        outterDOFID = dofID(1, ndID(i))
        do j=1,4
            innerDOFID = dofID(1, ndID(j))
            dynMat(outterDOFID,innerDOFID) = dynMat(outterDOFID,innerDOFID) + elmMat(i,j)
        enddo
    enddo   
    end subroutine assemDynMatK_c_r

 
!************************************************************************************************************************       
    
    subroutine assemElmC0Mat(elmMat, ndID)
    use solCtrlInf
    use nodeInf, only: dofID
    implicit none
    real:: elmMat(8,8)
    integer:: ndID(8)   
    integer:: outterDOFID, innerDOFID, i, j   
    
    do i =1,8
        outterDOFID = dofID(1, ndID(i))
        do j=1,8
            innerDOFID = dofID(1, ndID(j))
            C0mMat(outterDOFID,innerDOFID) = C0mMat(outterDOFID,innerDOFID) + elmMat(i,j)

        enddo
    enddo   
    end subroutine assemElmC0Mat

!************************************************************************************************************************    
    subroutine assemStaVec_k(elmVec, ndID)
    use mainCtrlInf, only: nNd
    use nodeInf, only: dofID
    use solCtrlInf
    implicit none
    real:: elmVec(8)
    integer:: ndID(nNd)
   
    integer:: outterDOFID, i    
    do i =1,8
        outterDOFID = dofID(1, ndID(i))
        staVec(outterDOFID) = staVec(outterDOFID) + elmVec(i)
    enddo
       
    end subroutine assemStaVec_k
	
!************************************************************************************************************************    
    subroutine assemStaVec_c(elmVec, ndID)
    use mainCtrlInf, only: nNd
    use nodeInf, only: dofID
    use solCtrlInf
    implicit none
    real:: elmVec(4)
    integer:: ndID(nNd)   
    integer:: outterDOFID, i    
    do i =1,4
        outterDOFID = dofID(1, ndID(i))
        staVec(outterDOFID) = staVec(outterDOFID) + elmVec(i)
    enddo
       
    end subroutine assemStaVec_c
!************************************************************************************************************************    
    
    subroutine assemDynVec_k(elmVec, ndID)
	use mainCtrlInf, only: nNd
    use nodeInf, only: dofID
    use solCtrlInf
    implicit none
    real:: elmVec(8)
    integer:: ndID(nNd)
   
    integer:: outterDOFID, i    
    do i =1,8
        outterDOFID = dofID(1, ndID(i))
        dynVec(outterDOFID) = dynVec(outterDOFID) + elmVec(i)
    enddo
	end subroutine assemDynVec_k
	
	subroutine assemDynVec_c_r(elmVec, ndID)
    use nodeInf, only: dofID
    use solCtrlInf
    implicit none
    real:: elmVec(4)
    integer:: ndID(4)
   
    integer:: outterDOFID, i    
    do i =1,4
        outterDOFID = dofID(1, ndID(i))
        dynVec(outterDOFID) = dynVec(outterDOFID) + elmVec(i)
    enddo
	end subroutine assemDynVec_c_r
	
