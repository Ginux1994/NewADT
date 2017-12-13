  
!********************************************************************       
    module mainCtrlInf ! 2.主控信息
    implicit none
    character(len = 100) title
    ! card 2.1
    ! 节点总数，线性单元组数，非线性单元组数， 求解类型，开始时间，
    integer:: nNd, nLrElmGrp, nNonlrElmGrp, calType, nFixTimeStep, startT
    integer:: nElmGrp, numEquation
    ! card 2.2
    ! 分析类型，集中节点比热数，相变界面数
    integer:: analyType, nSHnd, nLatht
    integer:: CmatType
    ! card 2.3
    ! solve the eigen or not
    integer:: eigenIndex
    ! card 2.4
    ! 重新形成矩阵之间的时间步数，平衡迭代之间的时间步数，平衡迭代最大次数，收敛误差
    integer:: nReformStep, nEquitTStep,nMaxEqStep
    real:: tolerance
    ! card 2.5
    ! 时间积分类型，α
    integer:: timeIntType
    real:: alpha
    ! card 2.6   
    integer:: NPB, IPRLH
    integer:: haveSHnode
    integer, allocatable:: npar(:,:)
contains      
        subroutine read_mainCtrlInf
        implicit none
            read(5,*) title
            read(5,*) nNd, nLrElmGrp, nNonlrElmGrp, calType, nFixTimeStep, startT
            nElmGrp = nLrElmGrp + nNonlrElmGrp
            allocate(npar(20,nElmGrp))
            read(5,*) analyType, nSHnd, nLatht
            if ((nFixTimeStep==0).or.(analyType==0)) nFixTimeStep = 1
            read(5,*) eigenIndex
            read(5,*) nReformStep, nEquitTStep, nMaxEqStep, tolerance
			if(nReformStep==0) nReformStep=1
			if(nEquitTStep) nEquitTStep = 1
			if(nMaxEqStep==0) nMaxEqStep=15
            read(5,*) timeIntType, alpha
            read(5,*) NPB, IPRLH
			
            if (timeIntType<=0) timeIntType=1
            if (timeIntType==1) then
                alpha = 1.0
            else if (timeIntType==2) then
                alpha = 0.0
            else if (timeIntType==3) then
                alpha = 0.5
            end if
			
            if (nSHnd>0) haveSHnode = 1  
			
            if (analyType==1) then 
                CmatType = 1
            else if(analyType==2) then
                CmatType = 0
            end if
        end subroutine read_mainCtrlInf
    end module mainCtrlInf
   
!********************************************************************           
    module solCtrlInf
    ! iper为当前迭代步ID，一共nFixTimeStep步（时间步长不变的时段数），每个固定时间步下的时间步长在timeStepInf中定义
    integer:: numLinearEqua, numNonLnEqua, numTotalEqua
    integer:: ind, index, iper, icount, nste, isref, iequit, nLrElm, nNonLrElm, klin

    real:: timePre, timeStart, timeNow
    integer:: numMasterDOF, numSlaveDOF
    real, allocatable:: KsMat(:, :), KmMat(:, :), KsmMat(:, :), C0sMat(:, :), C0mMat(:, :), C0smMat(:, :),C1Mat(:)
    real, allocatable:: Fs(:), Fm(:)  
    real, allocatable:: Phi(:)
    integer, allocatable:: DOFMap(:)!, ndStatus(:,1)

    ! character(len=5):: ind
contains    
    subroutine init
    use mainCtrlInf
        allocate(KmMat(numMasterDOF, numMasterDOF), Fm(numMasterDOF))
		allocate(C0mMat(numMasterDOF, numMasterDOF))
        allocate(C1Mat(numMasterDOF), Phi(numMasterDOF))
        KmMat = 0.0; Fm = 0.0; C1Mat = 0.0; Phi = 0.0
        numTotalEqua = numLinearEqua + numNonLnEqua
        isref = nReformStep; iequit = nEquitTStep;
        if(nLatht>0) klin=1
    end subroutine init
    
    end module    
    
    
    
!********************************************************************       
    module timeStepInf ! 3.时间步长信息
    implicit none
    ! 第i个求解时间段的时间步步数，及步长增量
    integer, allocatable:: stepNumAt(:)
    integer:: convergenFlag, nMaxStep
    real, allocatable:: dtAt(:)
    real:: a0, dt, dt_alpha, startTime
contains     
        subroutine read_timeStepInf
            use mainCtrlInf, only: timeIntType, alpha, nFixTimeStep
			use solCtrlInf
            integer:: i
            allocate(stepNumAt(nFixTimeStep), dtAt(nFixTimeStep))
            read(5,*)
            read(5,*) (stepNumAt(i), i=1, nFixTimeStep)
            read(5,*) (dtAt(i), i=1, nFixTimeStep)
			nste = stepNumAt(1)
			dt = dtAt(1)
            nMaxStep = 0			
            do i=1, nFixTimeStep
                if(stepNumAt(i)>nMaxStep) nMaxStep = stepNumAt(i)
            enddo
			if(klin.ne.0) isref = nMaxStep + 1
			
			
            call resetTimeIntgConst
        end subroutine read_timeStepInf
        
        subroutine resetTimeIntgConst  
            use mainCtrlInf, only: alpha, timeIntType
            use solCtrlInf, only: iper
            if(timeIntType<=4) then
                dt_alpha = dtAt(iper)*alpha
                if(alpha==0.0) a0 = 1.0/dtAt(iper)
                if(alpha>0.0) a0 = 1.0/dt_alpha
            end if     
        end subroutine resetTimeIntgConst
        
    end module timeStepInf
!********************************************************************    
    module nodeInf ! 4.节点信息
    implicit none
    integer:: coordSysType, ndDOF , indMaster = 0, indSlave = 0
    integer, allocatable:: dofID(:,:)
    real, allocatable:: x(:), y(:), z(:)    
contains
        subroutine read_nodeInf
            use mainCtrlInf, only: nNd, numEquation
            use solCtrlInf, only: numMasterDOF, numSlaveDOF
            integer:: i, nodeID, equaLoop=0
            ndDOF = 1
            allocate(dofID(ndDOF,nNd), x(nNd), y(nNd), z(nNd))
            do i=1,nNd
                read(5,*) nodeID, dofID(1,nodeID), x(nodeID), y(nodeID), z(nodeID)
                if(dofID(1,nodeID)==0) then
                    indMaster = indMaster + 1
                    dofID(1,nodeID) = indMaster
                else if(dofID(1,nodeID)==1) then
                    indSlave = indSlave + 1
                    dofID(1,nodeID) = -indSlave
                else if(dofID(1,nodeID)==-1) then
                    dofID(1,nodeID) = 0
                end if
            enddo
            numEquation = indMaster
            numMasterDOF = indMaster; numSlaveDOF = indSlave
        end subroutine read_nodeInf
    end module nodeInf
    
!********************************************************************    
    module ndSHInf ! 5.集中节点热容信息
    implicit none
    ! n号节点的集中节点热容
    real, allocatable:: ndSH(:), dql(:,:), qlv(:,:) !ndSH(nSHnd), dql(numEquation,nSHnd)
	integer, allocatable:: SHNdID(:)
contains
    subroutine read_ndSHInf
        use mainCtrlInf, only: nSHnd, nLatht, numEquation
        use nodeInf, only: dofID, ndDOF
        integer:: i, j, nodeID, id
        real:: nodeSH

        if (nSHnd>0) then
            allocate(ndSH(nSHnd), SHNdID(nSHnd), dql(numEquation,nLatht), qlv(numEquation,nLatht))
            ndSH = 0.0
            do i=1,nSHnd
                read(5,*) nodeID, nodeSH				
                do j=1,ndDOF
                    id = dofID(j, nodeID)
					SHNdID(i) = id
                    if(id.ne.0) ndSH(i) = nodeSH
                enddo
            enddo
        end if               
    end subroutine read_ndSHInf
    end module ndSHInf
!********************************************************************    
    
    module initInf ! 6.初始信息
    implicit none
    ! 0-自动形成0初始条件，初始温度条件
    integer:: autoZeroInit
    real, allocatable:: tempV(:)
contains
        subroutine read_initInf
            use mainCtrlInf, only: numEquation
            integer:: i
            allocate(tempV(numEquation))
            read(5,*) autoZeroInit
            if(autoZeroInit==0) then
                tempV = 0.0
            else 
                read(5,*) (tempV(i), i=1,numEquation)     
            end if   
        end subroutine read_initInf
    end module initInf
    
    
   
!********************************************************************    
    module heatFlowCtrlInf ! 7.热流控制卡 
    implicit none
    ! 7.1 control card
    ! 曲线数目，曲线点最大数目，对流类型，固定节点温度数目，对流节点数，辐射节点数目，集中热流数目，面分布热流数目，体分布热流数目
    integer:: nCurve, maxCurvePoint, convType, nFixTempNd, nConvNd, nRadiaNd, nLoad1, nLoad2, nLoad3, nHeatGenElmGrp, nHeatGenElm
    ! 7.2 phase transform temp相变界面温度，及增量、初始相变标志
    real, allocatable:: phaseTransfTemp(:), phaseTransfDt(:)
    integer, allocatable::initPhaseFlag(:)
    integer:: haveLoad=0, haveLnConvOrFixTempLoad=0, nscr   
contains
        subroutine read_heatFlowCtrlInf
            use mainCtrlInf, only: nLatht
            integer:: i
            allocate(phaseTransfTemp(nLatht), phaseTransfDt(nLatht), initPhaseFlag(nLatht))
            read(5,*) nCurve, maxCurvePoint, convType, nFixTempNd, nConvNd, nRadiaNd, nLoad1, nLoad2, nLoad3, nHeatGenElmGrp, nHeatGenElm
            if(nLatht>0) then
                read(5,*) (phaseTransfTemp(i), i=1, nLatht)
                read(5,*) (phaseTransfDt(i), i=1, nLatht)
                read(5,*) (initPhaseFlag(i), i=1, nLatht) 
            end if
            if((nConvNd + nRadiaNd + nLoad1 + nLoad2 + nLoad3)>0) haveLoad=1
            if((convType==1).or.(nFixTempNd>0)) haveLnConvOrFixTempLoad=1 !convType==1表示线性对流
            nscr = nFixTempNd + nConvNd + nRadiaNd
        end subroutine read_heatFlowCtrlInf
    end module
!*****************************************************************************    
   
    module timeFuncInf ! 8. 时间函数曲线信息，包括外部热流的以及材料的
    implicit none
    real, allocatable:: timeAt(:,:), timeFuncAt(:,:)  
contains
        subroutine read_timeFuncInf
            use heatFlowCtrlInf, only: nCurve, maxCurvePoint
            integer:: curveID, nCurvePoints, i, j
        
            allocate(timeAt(maxCurvePoint, nCurve), timeFuncAt(maxCurvePoint, nCurve))
            do i=1,nCurve   
                read(5,*) curveID, nCurvePoints                 
                read(5,*) (timeAt(j, curveID), timeFuncAt(j, curveID), j=1,nCurvePoints)
            end do
        end subroutine read_timeFuncInf
    end module