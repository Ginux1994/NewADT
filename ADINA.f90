 
    
!*****************************************************************************  
    subroutine equitIte!(Phi_incTol)
    use mainCtrlInf, only: numEquation, nMaxEqStep, nLatht, CmatType, haveSHnode, nLrElmGrp, nNonlrElmGrp, tol_input
    use solCtrlInf !, only: C1Mat, C0Mat, icount
    use heatFlowCtrlInf, only: haveLoad, haveLnConvOrFixTempLoad
    use heatFlowInf, only: nFixTempNd, execute_FixTempNd, fixTempNdInf, fixTempNdValue!, loadStepLoop
    use ndSHInf, only: dql, qlv
    use timeStepInf, only: a0, convergenFlag
    use LU
    use element, only: assem_element
    implicit none
    integer:: iteLoop, i, j, fixTempNdID
    
    real:: workVec(numEquation), outHeatF(numEquation), tempPre(numEquation), tempNow(numEquation), &
            effKmat(numEquation, numEquation)
    real:: norm_PhiIncTol = 0.0,norm_PhiIncNow = 0.0, norm_outHeatF = 0.0
    real:: norm_Now=0.0, tol_now 
    real:: dummy
    workVec = 0.0 !call initRvec(workVec, numEquation)   
    icount = 3
    if (nLatht>0) then
        ! call pahse, 调用潜热会计算一次温度增量，即温度增量因潜热有一个初值，但暂时先赋零值
		
    end if
    
    Phi_incTol = 0.0; Phi_incNow = 0.0
    
    do iteLoop=1, nMaxEqStep 
        
! 1. 收敛指标 及 中间矩阵变量 归零:
        norm_PhiIncTol = 0.0; norm_PhiIncNow = 0.0; norm_outHeatF = 0.0
        norm_Now=0.0;

		dynMat = 0.0; dynVec = 0.0;	
        
! 2. get the 外部对流及固定温度载荷        
        if (haveLnConvOrFixTempLoad>0) then

            dynVec = dynVec - matmul(staMatK_k, Phi_now)

        end if
! 3. consider the 热容的作用        
        if (haveSHnode>0) then
            workVec = a0*Phi_incTol
            if (CmatType == 1) then ! lump capacity matrix
                do i=1,numEquation
                    dynVec(i) = dynVec(i) - C1Mat(i)*workVec(i)
                enddo
            else if (CmatType == 2) then ! consistent capacity matrix  
                do i=1,numEquation  ! dynVec = dynVec - matmul(C0Mat, Phi_now)
                    dummy=0.0
                    do j=1,numEquation
                        dummy = dummy + C0mMat(i,j)*workVec(j)
                    enddo
                    dynVec(i) = dynVec(i) - dummy
                enddo
                !outHeatF = outHeatF - C0mat*workVec
            end if

        end if
        
! 4. consider the 潜热的作用
        if (nLatht>0) then
            do i=1, nLatht
                do j=1, numEquation
                    dynVec(j) = dynVec(j) - dql(j,i)*a0
                enddo
            enddo
        end if
! 5. deal with the 非线性单元矩阵，如辐射的作用
		ind=4
        do i = 1, nNonlrElmGrp
                ! read the control file 
            call assem_element !计算并累加 非线性单元产生的 动态mat和动态vec
        end do              
        totMat = staMat + dynMat; totVec = staVec + dynVec
        norm_outHeatF = norm2(totVec)
		if(nFixTempNd>0) call execute_FixTempNd   
! 6. calculate the increment of the temperatures
        

        if (iteLoop<(nMaxEqStep).or.(norm_outHeatF<rnorm)) then             
                call LUsolve(totMat,totVec,Phi_incNow,numEquation)            
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
        
        !
        do i=1,nFixTempNd
            fixTempNdID = fixTempNdInf(1,i)
            Phi_incNow(fixTempNdID) = 0.0   
        enddo
            
        Phi_incTol = Phi_incTol + Phi_incNow
        Phi_now = Phi_now + Phi_incNow
        norm_PhiIncTol = norm2(Phi_incTol)
        norm_PhiIncNow = norm2(Phi_incNow)

        ! 更新收敛准则
        if (norm_Now<norm_PhiIncNow) norm_Now=norm_PhiIncTol
        tol_now = norm_Now*tol_input
        ! 判断是否达到最大迭代次数
        if(norm_PhiIncNow<0.000001) exit  !if(norm_PhiIncNow<tol_now) exit
        ! iteLoop = iteLoop + 1      
        if (iteLoop>nMaxEqStep) then
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
    
    do i=1,nFixTempNd
        fixTempNdID = fixTempNdInf(1,i)
        Phi_now(fixTempNdID) = fixTempNdValue(loadLoop,fixTempNdID)
    enddo
    
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
    
    subroutine mult1(vecA, matB, vecC, row, column)
    implicit none
    integer:: row,column
    real:: vecA(column), matB(row,column), vecC(column)
    integer::i, j   
    vecA = vecA - matmul(matB, vecC)
    end subroutine mult1

