 
    
!*****************************************************************************  
    subroutine equitIte(effKmat, outHeatF, tempInc, tempPre) 
    use mainCtrlInf, only: numEquation, nEquitTStep, nLatht, CmatType, haveSHnode, nLrElmGrp, nNonlrElmGrp
    use solCtrlInf !, only: C1Mat, C0Mat, icount
    use heatFlowCtrlInf, only: haveLoad, haveLnConvOrFixTempLoad
    use ndSHInf, only: dql, qlv
    use timeStepInf, only: a0, convergenFlag
    
    implicit none
    integer:: iteLoop, i, j
    
    real:: workVec(numEquation), outHeatF(numEquation), tempInc(numEquation), tempPre(numEquation), tempNow(numEquation), &
            effKmat(numEquation, numEquation)
    real:: norm_tempInc = 0.0,norm_tempIncNow = 0.0, norm_tempNow = 0.0, norm_outHeatF = 0.0
    real:: norm_temp=0.0, tol_input = 0.0, tol_now, rnorm = 0.0
    
    workVec = 0.0 !call initRvec(workVec, numEquation)   
    icount = 3
    if (nLatht>0) then
        ! call pahse, calculate workVec
    end if
    tempInc = tempInc + workVec
    tempNow = tempPre + tempInc
    
    do while(norm_tempIncNow<tol_now)
        
! 1. get the 外部作用热载荷        
        if(haveLoad>0) then
            ! read 3
        else 
            do i=1,numEquation
               outHeatF(i) = 0.0
            end do
        end if
! 2. get the 外部对流及固定温度载荷        
        if (haveLnConvOrFixTempLoad>0) then
            outHeatF = outHeatF - matmul(effKmat, tempNow)
        end if
! 3. consider the 热容的作用        
        if (haveSHnode>0) then
            workVec = a0*tempInc
            if (CmatType == 1) then ! lump capacity matrix
                do i=1,numEquation
                    outHeatF(i) = outHeatF(i) - C1Mat(1)*workVec(i)
                enddo
            else if (CmatType == 2) then ! consistent capacity matrix  
                call mult(outHeatF, C0mMat, workVec, numEquation, numEquation)
                !outHeatF = outHeatF - C0mat*workVec
            end if
            
            tempNow = tempPre + tempInc
        end if
        
! 4. consider the 潜热的作用
        if (nLatht>0) then
            do i=1, nLatht
                do j=1, numEquation
                    outHeatF(j) = outHeatF(j) - dql(j,i)*a0
                enddo
            enddo
        end if
! 5. deal with the 非线性单元矩阵，如辐射的作用
        if (nNonlrElmGrp>0) then 
            do i = 1, nNonlrElmGrp
                ! read the control file 
                 ! call element
            end do 
        end if
! 6. calculate the increment of the temperatures
        do i=1, numEquation
            norm_outHeatF = norm_outHeatF + outHeatF(i)*outHeatF(i)
        end do
        if (iteLoop>(nEquitTStep/2+2).or.(norm_outHeatF<rnorm)) then 
                ! call colsol
        else 
            write(*,*) 'out of balance heat flows larger than incremental heat flows after iteration,i5/4x,20hstop if number eq. 0'
            convergenFlag = 2
            return 
        end if  
        
! 7. check convergence or not
        if (nLatht>0) then
            ! call phase
        end if
        ! 现在 outHeatF 存储着 新计算出来的t时刻下，第i个迭代步下的 温度增量
        ! 计算当前迭代步下 的 温度、温度增量及 温度范数、增量范数
        do i=1,numEquation
            tempInc(i) = tempInc(i) + outHeatF(i)
            tempNow(i) = tempNow(i) + outHeatF(i)
            norm_tempNow = tempNow(i)*tempNow(i) !dnorm
            norm_tempIncNow = outHeatF(i)*outHeatF(i)
        enddo
        ! 更新收敛准则
        if (norm_temp<norm_tempNow) norm_temp=norm_tempNow
        tol_now = norm_temp*tol_input
        ! 判断是否达到最大迭代次数
        iteLoop = iteLoop + 1      
        if (iteLoop>nEquitTStep) then
            write(*,*) 'iteration limit reached with no convergence stop   of solution  '
            return 
        end if
        
        
        icount = 2
        if (nLatht>0) then
            ! 计入潜热增量qlv
            qlv = qlv - dql
            dql = 0.0
        end if
        
       
    end do
  
    end subroutine equitIte
    
    
!*****************************************************************************       
    subroutine updateTemp(tempPre, tempInc)
    use solCtrlInf, only: klin, nLrElm, nNonLrElm
    use mainCtrlInf, only: analyType, timeIntType, alpha, numEquation
    implicit none
!timeIntType 
        ! eq.1, euler backward method
        ! eq.2, euler forward method 
        ! eq.3, trapezoidal rule     
        ! eq.4, alpha, family method 
    real:: tempPre(numEquation), tempInc(numEquation), beta, fact
    integer:: i
    if(klin==0) then ! linear analysis
        if(timeIntType<3) then  ! 欧拉前推，后插
            do i=1,numEquation
                tempPre(i) = tempInc(i)
            enddo
        else 
            fact = alpha
            if (alpha==0.0) fact = 1.0
            beta = (1.0-fact)/fact
            fact = 1.0/fact
            do i=1,numEquation
                tempPre = tempInc(i)*fact - tempPre(i)*beta
            enddo
        end if
    else ! non linear analysis
        if(timeIntType==2) then ! 欧拉前推
            do i=1,numEquation
                tempPre(i) = tempInc(i)
            enddo
        else if(timeIntType>=3) then
            fact = alpha
            if (alpha==0.0) fact = 1.0
            fact = 1.0/fact
            do i=1,numEquation
                tempPre = tempInc(i)*fact + tempPre(i)
            enddo
        else 
            do i=1,numEquation
                tempPre(i) = tempPre(i) + tempInc(i)
            enddo
        end if
    end if  
    
    end subroutine updateTemp
    
!*****************************************************************************      
    
    subroutine mult(vecA, matB, vecC, row, column)
    implicit none
    integer:: row,column
    real:: vecA(column), matB(row,column), vecC(column)
    integer::i, j   
    vecA = vecA - matmul(matB, vecC)
    end subroutine mult

 !*****************************************************************************    
    subroutine addla(ndHeatFlow, lm)
    use mainCtrlInf, only: nNd
    use solCtrlInf, only: numTotalEqua, DOFMap, numMasterDOF, numSlaveDOF, Fs, Fm
    implicit none
    real:: ndHeatFlow(8), lm(nNd)
   
    integer:: outterDOFID, i    
    do i =1,8
        outterDOFID = DOFMap(lm(i))
        if (outterDOFID>0) then
            Fm(outterDOFID) = Fm(outterDOFID) - ndHeatFlow(i)
        else if (outterDOFID<0) then
            Fs(-outterDOFID) = Fs(-outterDOFID) - ndHeatFlow(i)
        end if
    enddo      
    end subroutine addla
 !*****************************************************************************    
    subroutine addma(ndHeatFlow, lm)
    use mainCtrlInf, only: nNd
    use solCtrlInf, only: numTotalEqua, DOFMap, numMasterDOF, numSlaveDOF, Fs, Fm
    implicit none
    real:: ndHeatFlow(8), lm(nNd)
   
    integer:: outterDOFID, i    
    do i =1,8
        outterDOFID = DOFMap(lm(i))
        if (outterDOFID>0) then
            Fm(outterDOFID) = Fm(outterDOFID) + ndHeatFlow(i)
        else if (outterDOFID<0) then
            Fs(-outterDOFID) = Fs(-outterDOFID) + ndHeatFlow(i)
        end if
    enddo      
    end subroutine addma 
    

