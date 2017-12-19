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
    ! ֻ����ն��󣬲���������Ҷ���
    call assem_element ! ����sta�����նȾ��� staMatK_k; �ۼ�sta�����նȾ���staMat
        
    
    
    ind=2
    if(iper>1) call ADINI


    call execute_heatFlow  ! �ۼӵ�iper�ܲ��¸���loadStepLoopʱ�䲽�µļ�������+�����Ҷ������sta��̬�ģ��浽Fm_step��
                            ! ����������dyn�У���ÿloadStepLoop������£���staVecֱ�Ӵ�Fm_step��ȡ������
    
    do loadLoop=1,stepNumAt(iper)

        staVec = Fm_step(:,loadLoop)
        
    ! ��ʼ����ʱ�䲽����ѭ��
    
        timeNow = timeStart
        ind = 4
    
        ksref = 0; kequit = 0; ite = 0; icount = 2;

        staMat = staMatK_c + staMatK_k + staMatC0
        
        
! ʱ�����ѭ��        

        timeNow = timeNow + dt
        
        ksref = ksref + 1
        iref = isref - ksref    ! if(klin.ne.0) isref = nMaxStep + 1
        if(iref==0) ksref = 0
        if(loadLoop==1) iref = 0   ! irefΪ0ʱ���������γɴ�������
        
        kequit = kequit + 1
        iit = iequit - kequit   !  if(klin==0) iequit = nEquitTStep; ����iequit=1
                                ! iitΪ1������Ҫ������Ϊ0����Ҫ������һ���һ��ʱ�䲽����Ҫ
        if(iit==0) kequit = 0
        
        nMatReform = 0  ! ��¼�����Դ������������γɵĴ���
        convergenceFlag = 0     ! ��¼ƽ�ⲽ�����Ƿ�������=0����
        
! now we begin to iteration
        
! 1. ����˿̵���������
        Phi_pre = Phi_now
        ! �����Ҷ���
        if(CmatType==0) then
            ! staVec = staVec - C0mMat*Phi_incTol
            staVec = staVec + matmul(staMatC0, Phi_pre)
        else if(CmatType==1) then
            do i=1,numEquation
                staVec(i) = staVec(i) + staMatC0(i,i)*Phi_pre(i)
            enddo
        end if
! 2. �����Ч������ ��������
        ind = 4
        call assem_element ! ��װ�����Ծ��� dynMat ��
        totMat = staMat + dynMat
        totVec = staVec + dynVec
        rnorm = norm2(totVec)
        
        if((klin==0).and.(timeIntType.ne.2)) then
            if(nFixTempNd>0) call execute_FixTempNd   
            call LUsolve(totMat,totVec,Phi_now,numEquation) ! ��������ֱ�Ӽ����ok�����õ���
        end if
!        call updateTemp(Phi_now, Phi_incTol)
        loopNum = loopNum+1
        if((klin>0).and.(iit==0).and.(alpha/=0.0)) then
            call equitIte
        end if
! ���½������������м����    
        ! �¶ȸ����ڵ���ʱ�Ѿ������ˣ����ﲻ��Ҫ����Phi_Now = Phi_Now + Phi_incTol
        
        ! ֱ�Ӹ�ֵ���ǣ�û���ۼӣ����Բ���Ҫ������ totMat = 0.0; totVec = 0.0;
        staVec = 0.0; staMat = 0.0! staMat���Ѿ��ж����ն����Ҳ��䣬���ܹ���
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