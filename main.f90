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
    use LU

    implicit none
    integer:: i,j,k, convergenceFlag, iit, ksref, nMatReform, kequit = 0, ite = 0
    integer:: loopNum=0
! read in the information, no calculation(ind=0)    
    character(len=20) outputFile
    outputFile = "OUT.DAT"
    open(unit = 6, file = outputFile)
    
    
    iper = 1;
    ind = 0
    call ADINI
    
do iper=1,nFixTimeStep
    
    totMat=0.0; totVec = 0.0;
    staMat=0.0; staVec = 0.0; dynMat = 0.0; dynVec = 0.0;
    staMatK_c = 0.0; staMatK_k = 0.0; staMatC0 = 0.0;
    C0mMat = 0.0; C1Mat = 0.0; Fm_step=0.0; Phi_incTol=0.0; Phi_pre = 0.0; Phi_incNow=0.0

    ind = 1; isref = 0; index=0; nReformStep = 0;
    ! 只计算刚度阵，不计算对流右端项
    call assem_element ! 计算sta传导刚度矩阵 staMatK_k; 累加sta对流刚度矩阵到staMat
        
    
    
    ind=2
    if(iper>1) call ADINI


    call execute_heatFlow  ! 累加第iper总步下各个loadStepLoop时间步下的集中热流+对流右端项，都是sta静态的，存到Fm_step中
                            ! 而对流放在dyn中，并每loadStepLoop都会更新，而staVec直接从Fm_step中取就行了
    
    do loadLoop=1,stepNumAt(iper)

        staVec = Fm_step(:,loadLoop)
        
    ! 开始进行时间步迭代循环
    
        timeNow = timeStart
        ind = 4
    
        ksref = 0; kequit = 0; ite = 0; icount = 2;

        staMat = staMatK_c + staMatK_k + staMatC0
        
        
! 时间积分循环        

        timeNow = timeNow + dt
        
        ksref = ksref + 1
        iref = isref - ksref    ! if(klin.ne.0) isref = nMaxStep + 1
        if(iref==0) ksref = 0
        if(loadLoop==1) iref = 0   ! iref为0时，不重新形成传导矩阵
        
        kequit = kequit + 1
        iit = iequit - kequit   !  if(klin==0) iequit = nEquitTStep; 否则iequit=1
                                ! iit为1，则需要迭代，为0不需要迭代，一般第一个时间步不需要
        if(iit==0) kequit = 0
        
        nMatReform = 0  ! 记录非线性传到矩阵重新形成的次数
        convergenceFlag = 0     ! 记录平衡步迭代是否收敛，=0收敛
        
! now we begin to iteration
        
! 1. 计算此刻的热流向量
        Phi_pre = Phi_now
        ! 热容右端项
        if(CmatType==0) then
            ! staVec = staVec - C0mMat*Phi_incTol
            staVec = staVec + matmul(staMatC0, Phi_pre)
        else if(CmatType==1) then
            do i=1,numEquation
                staVec(i) = staVec(i) + staMatC0(i,i)*Phi_pre(i)
            enddo
        end if
! 2. 计算等效非线性 传到矩阵
        ind = 4
        call assem_element ! 总装非线性矩阵到 dynMat 中
        totMat = staMat + dynMat
        totVec = staVec + dynVec
        rnorm = norm2(totVec)
        
        if((klin==0).and.(timeIntType.ne.2)) then
            if(nFixTempNd>0) call execute_FixTempNd   
            call LUsolve(totMat,totVec,Phi_now,numEquation) ! 线性问题直接计算就ok，不用迭代
        end if
!        call updateTemp(Phi_now, Phi_incTol)
        loopNum = loopNum+1
        if((klin>0).and.(iit==0).and.(alpha/=0.0)) then
            call equitIte
        end if
! 更新解向量，重置中间矩阵    
        ! 温度更新在迭代时已经做好了，这里不需要再做Phi_Now = Phi_Now + Phi_incTol
        
        ! 直接给值覆盖，没有累加，所以不需要归零了 totMat = 0.0; totVec = 0.0;
        staVec = 0.0; staMat = 0.0! staMat里已经有对流刚度阵且不变，不能归零
		dynMat = 0.0; dynVec = 0.0;
        write(6,*) loopNum, (Phi_now(i), i=1,numEquation)
    enddo !do loadLoop=1,stepNumAt(iper)
enddo !do iper=1,nFixTimeStep   

    
    end program
 
    subroutine mult(vecA, matB, vecC)
    use mainCtrlInf, only:numEquation
    integer::i,j,k
    real::vecA(numEquation), vecC(numEquation), matB(numEquation,numEquation), dummy
    do i=1,numEquation
        dummy=0.0
        do j=1,numEquation
            dummy = dummy + matB(i,j)*vecC(j)
        enddo
        vecA(i) = vecA(i) - dummy
    enddo
    end subroutine mult