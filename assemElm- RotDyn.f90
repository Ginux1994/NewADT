   
    
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
    
    subroutine assemElmMat(elmMat, lm)
    use mainCtrlInf, only: nNd
    use solCtrlInf, only: numTotalEqua, DOFMap, numMasterDOF, numSlaveDOF, KmMat, KsMat, KsmMat    
    implicit none
    real:: elmMat(8,8)
    integer:: lm(1)   
    integer:: outterDOFID, innerDOFID, i, j   
    

    do i =1,8
        outterDOFID = DOFMap(lm(i))
        do j=1,8
            innerDOFID = DOFMap(lm(j))
            if (outterDOFID>0) then
                if (innerDOFID>0) then 
                    KmMat(outterDOFID,innerDOFID) = KmMat(outterDOFID,innerDOFID) + elmMat(i,j)
                    if (outterDOFID /= innerDOFID)  KmMat(innerDOFID,outterDOFID) = KmMat(outterDOFID,innerDOFID)
                else if (innerDOFID<0) then
                    KsmMat(outterDOFID,-innerDOFID) = KsmMat(outterDOFID,-innerDOFID) + elmMat(i,j)
                end if
            else if (outterDOFID<0) then
                if (innerDOFID>0) then 
                    KsmMat(innerDOFID,-outterDOFID) = KsMat(innerDOFID,-outterDOFID) + elmMat(i,j)
                else if (innerDOFID<0) then
                    KsMat(-outterDOFID,-innerDOFID) = KsMat(-outterDOFID,-innerDOFID) + elmMat(i,j)
                    if (outterDOFID /= innerDOFID) KsMat(-innerDOFID,-outterDOFID) = KsMat(-outterDOFID,-innerDOFID)
                end if
            end if
        enddo
    enddo   
    end subroutine assemElmMat
!************************************************************************************************************************    
    subroutine assemElmVec(elmVec, lm)
    use mainCtrlInf, only: nNd
    use solCtrlInf, only: numTotalEqua, DOFMap, numMasterDOF, numSlaveDOF, Fs, Fm
    implicit none
    real:: elmVec(8), lm(nNd)
   
    integer:: outterDOFID, i    
    do i =1,8
        outterDOFID = DOFMap(lm(i))
        if (outterDOFID>0) then
            Fm(outterDOFID) = Fm(outterDOFID) + elmVec(i)
        else if (outterDOFID<0) then
            Fs(-outterDOFID) = Fs(-outterDOFID) + elmVec(i)
        end if
    enddo
       
    end subroutine assemElmVec