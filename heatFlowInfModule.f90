    
!*****************************************************************************    
    module heatFlowInf
    use heatFlowCtrlInf, only: nscr, nFixTempNd, nConvNd, nRadiaNd, nLoad1, nLoad2, nLoad3, &
                                nHeatGenElmGrp, nHeatGenElm, maxCurvePoint, haveLrConv
	use solCtrlInf, only: iper
	use mainCtrlInf, only: startT_input, nFixTimeStep
	use timeStepInf
	use timeFuncInf
    implicit none
    integer, allocatable:: fixTempNdInf(:,:), convNdInf(:,:),radiaNdInf(:,:), load1Inf(:,:), heatGenInf(:,:), heatGenElmID(:,:)
    integer:: loadStepLoop, loadStepNow
	real, allocatable::  fixTempNdValue(:,:), convNdValue(:,:), radiaNdValue(:,:), load1Value(:,:), heatGenValue(:,:)
contains
        subroutine read_heatFlowInf
        use timeStepInf, only: stepNumAt, dtAt, nMaxStep
        use mainCtrlInf, only: startT_input
        use solCtrlInf, only: timeStart
        if(iper==1) then
            timeStart = startT_input
            !allocate(fixTempNdInf(4, nFixTempNd), fixTempNdValue(stepNumAt(iper), nFixTempNd))! , fixTempNdValue()
            !allocate(convNdInf(4,nConvNd), convNdValue(stepNumAt(iper), nConvNd))
            !allocate(radiaNdInf(4,nRadiaNd), radiaNdValue(stepNumAt(iper), nRadiaNd))
            !allocate(load1Inf(4,nLoad1), load1Value(stepNumAt(iper), nLoad1))
            !allocate(heatGenInf(4,nHeatGenElmGrp), heatGenElmID(nHeatGenElm, nHeatGenElmGrp), heatGenValue(stepNumAt(iper), nHeatGenElm))               
            allocate(fixTempNdInf(4, nFixTempNd), fixTempNdValue(nMaxStep, nFixTempNd))! , fixTempNdValue()
            allocate(convNdInf(4,nConvNd), convNdValue(nMaxStep, nConvNd))
            allocate(radiaNdInf(4,nRadiaNd), radiaNdValue(nMaxStep, nRadiaNd))
            allocate(load1Inf(4,nLoad1), load1Value(nMaxStep, nLoad1))
            allocate(heatGenInf(4,nHeatGenElmGrp), heatGenElmID(nHeatGenElm, nHeatGenElmGrp), heatGenValue(nMaxStep, nHeatGenElm))   
        else 
            timeStart = timeStart + stepNumAt(iper-1)*dtAt(iper-1)
        end if
            if(nFixTempNd>0) call read_fixTempNdFlow
            if(nConvNd>0) call read_ConvNdFlow
            if(nRadiaNd>0) call read_RadiaNdFlow
            if(nLoad1>0) call read_Load1
			if(nHeatGenElm>0) call read_HeatGen
        end subroutine read_heatFlowInf

        subroutine read_fixTempNdFlow
        use heatFlowCtrlInf, only: nCurve, maxCurvePoint
        use solCtrlInf, only: timeStart
        use timeFuncInf, only: timeAt
        use timeStepInf, only: stepNumAt
            integer:: i, j, k, l, m, ndID, curveID, fac, kn, pointLoop
			real:: time_now, time_pre, timeTmp, startTime, sloop, arg
            


            do i=1, nFixTempNd
                read(5,*) ndID, curveID, fac, startTime, kn
                if (fac==0.0) fac = 1.0
                fixTempNdInf(1,i) = ndID; fixTempNdInf(2,i) = curveID; fixTempNdInf(3,i) = fac; fixTempNdInf(4,i) = startTime;
                time_now = timeStart

                do j=1,stepNumAt(iper)
                    time_pre = time_now + dt_alpha
                    time_now = time_now + dtAt(iper)
                    timeTmp = time_pre - startTime
                    if(timeTmp<=0.0) return
                    pointLoop = 0

                    do k=2,maxCurvePoint
                        pointLoop = pointLoop + 1
                        if((timeTmp>=timeAt(k-1, curveID)).and.(timeTmp<=timeAt(k, curveID))) exit
                    enddo
                    sloop = (timeFuncAt(pointLoop+1, curveID) - timeFuncAt(pointLoop, curveID))/&
                        (timeAt(pointLoop+1, curveID) - timeAt(pointLoop, curveID))
                    arg = timeFuncAt(pointLoop, curveID) + sloop*(timeTmp - timeAt(pointLoop, curveID))
                    fixTempNdValue(j, i) = arg*fac
                enddo

            enddo

                   
        end subroutine read_fixTempNdFlow
    

	
        subroutine read_ConvNdFlow
        use solCtrlInf, only: timeStart
            integer:: i, j, k, pointLoop, ndID, curveID, fac, startTime, kn
            real:: arg, sloop, time_now, time_pre, timeTmp
            
            

            do i=1,nConvNd
                read(5,*) ndID, curveID, fac, startTime, kn
                if (fac==0.0) fac = 1.0
                convNdInf(1,i) = ndID; convNdInf(2,i) = curveID; convNdInf(3,i) = fac; convNdInf(4,i) = startTime;

                
                
                
                time_now = timeStart
                do j=1,stepNumAt(iper)
                    time_pre = time_now + dt_alpha
                    time_now = time_now + dtAt(iper)
                    timeTmp = time_pre - startTime
                    if(timeTmp<=0.0) return
                    pointLoop = 0
                    do k=2,maxCurvePoint
                        pointLoop = pointLoop + 1
                        if((timeTmp>=timeAt(k-1, curveID)).and.(timeTmp<=timeAt(k, curveID))) exit
                    enddo
                    sloop = (timeFuncAt(pointLoop+1, curveID) - timeFuncAt(pointLoop, curveID))/&
                        (timeAt(pointLoop+1, curveID) - timeAt(pointLoop, curveID))
                    arg = timeFuncAt(pointLoop, curveID) + sloop*(timeTmp - timeAt(pointLoop, curveID))
                    convNdValue(j, i) = arg*fac
                enddo

            enddo

                      
        end subroutine read_ConvNdFlow
    
        subroutine read_RadiaNdFlow
        use solCtrlInf, only: timeStart
            integer:: i, k, j, pointLoop
            integer:: ndID, curveID, fac, startTime, kn
            real:: arg, sloop, time_now, time_pre, timeTmp

            do i=1,nRadiaNd
                read(5,*) ndID, curveID, fac, startTime, kn
                if (fac==0.0) fac = 1.0
                radiaNdInf(1,i) = ndID; radiaNdInf(2,i) = curveID; radiaNdInf(3,i) = fac; radiaNdInf(4,i) = startTime;

                time_now = timeStart
                do j=1,stepNumAt(iper)
                    time_pre = time_now + dt_alpha
                    time_now = time_now + dtAt(iper)
                    timeTmp = time_pre - startTime
                    if(timeTmp<=0.0) return
                    pointLoop = 0
                    do k=2,maxCurvePoint
                        pointLoop = pointLoop + 1
                        if((timeTmp>=timeAt(k-1, curveID)).and.(timeTmp<=timeAt(k, curveID))) exit
                    enddo
                    sloop = (timeFuncAt(pointLoop+1, curveID) - timeFuncAt(pointLoop, curveID))/&
                        (timeAt(pointLoop+1, curveID) - timeAt(pointLoop, curveID))
                    arg = timeFuncAt(pointLoop, curveID) + sloop*(timeTmp - timeAt(pointLoop, curveID))
                    radiaNdValue(j, i) = arg*fac
                enddo

            enddo

                    
        end subroutine read_RadiaNdFlow
    
        subroutine read_Load1
        use solCtrlInf, only: timeStart
            integer:: i, j, k, pointLoop, ndID, curveID, fac, startTime, kn
            real:: arg, sloop, time_now, time_pre, timeTmp
            load1Value = 0.0
            do i=1, nLoad1
                read(5,*) ndID, curveID, fac, startTime, kn
                if (fac==0.0) fac = 1.0
                load1Inf(1,i) = ndID; load1Inf(2,i) = curveID; load1Inf(3,i) = fac; load1Inf(4,i) = startTime;

                time_now = timeStart
                do j=1,stepNumAt(iper)
                    time_pre = time_now + dt_alpha
                    time_now = time_now + dtAt(iper)
                    timeTmp = time_pre - startTime
                    if(timeTmp<=0.0) return
                    pointLoop = 0
                    do k=2,maxCurvePoint
                        pointLoop = pointLoop + 1
                        if((timeTmp>=timeAt(k-1, curveID)).and.(timeTmp<=timeAt(k, curveID))) exit
                    enddo
                    sloop = (timeFuncAt(pointLoop+1, curveID) - timeFuncAt(pointLoop, curveID))/&
                        (timeAt(pointLoop+1, curveID) - timeAt(pointLoop, curveID))
                    arg = timeFuncAt(pointLoop, curveID) + sloop*(timeTmp - timeAt(pointLoop, curveID))
                    load1Value(j, i) = load1Value(j, i) + arg*fac
                enddo
            enddo
              
        end subroutine read_Load1
        
		subroutine read_HeatGen		
        use timeStepInf
        use solCtrlInf, only: timeStart
        use timeFuncInf, only: timeAt, timeFuncAt
		    integer:: i, j, k, l, pointLoop,&
                      linearType, grpID, nElm,&
                      elmID, curveID, fac, startTime, kn
            integer:: stepNumNow
            real:: arg, sloop, time_now, time_pre, timeTmp
            stepNumNow = stepNumAt(iper)
                      
            do i=1, nHeatGenElmGrp
				
                read(5,*) linearType, grpID, nElm
				heatGenInf(1,i) = linearType; heatGenInf(2,i) = grpID; heatGenInf(3,i) = nElm;
				
				do j=1, nElm
					read(5,*) elmID, curveID, fac, startTime, kn
                	heatGenElmID(j, i) = elmID	
					if (fac==0.0) fac = 1.0
					time_now = timeStart
					do l=1,stepNumAt(iper)
						time_pre = time_now + dt_alpha
						time_now = time_now + dtAt(iper)
						timeTmp = time_pre - startTime
						if(timeTmp<=0.0) return
						pointLoop = 0
						do k=2,maxCurvePoint
							pointLoop = pointLoop + 1
							if((timeTmp>=timeAt(k-1, curveID)).and.(timeTmp<=timeAt(k, curveID))) exit
						enddo
						sloop = (timeFuncAt(pointLoop+1, curveID) - timeFuncAt(pointLoop, curveID))/&
							(timeAt(pointLoop+1, curveID) - timeAt(pointLoop, curveID))
						arg = timeFuncAt(pointLoop, curveID) + sloop*(timeTmp - timeAt(pointLoop, curveID))
						heatGenValue(l, j) = arg*fac ! 热生成数据按照单元顺序存储，其相应的单元编号通过heatGenElmID(j, i) = elmID	获取
                    enddo
                enddo
            enddo

		end subroutine read_HeatGen
		

		
		
        subroutine read_Load2
        use solCtrlInf, only: timeStart
        
        
        end subroutine read_Load2
                
        subroutine read_Load3
        use solCtrlInf, only: timeStart
        
        
        end subroutine read_Load3
             
!*********************************************************************************************    
        
        subroutine execute_heatFlow
        use element, only: execute_HeatGen
        use convElm, only: assem_convElm
        use timeStepInf, only: stepNumAt, nMaxStep, a0
        use solCtrlInf, only: iper, Fm_step, numMasterDOF, &
                              C0mMat, C1Mat, staMatC0, staMat, staVec
        use mainCtrlInf, only: numEquation, CmatType, timeIntType
        integer:: i, j, k
        real:: a0_
        
        
        do loadStepLoop=1, stepNumAt(iper)
            if(nFixTempNd>0) call execute_FixTempNd
            if(haveLrConv>0) call assem_convElm  ! 线性对流矩阵及产生的热流在 组装刚度阵的时候就已经完成了，不再在这个地方重新计算
            ! if(nRadiaNd>0) ! call execute_RadiaNdFlow
            if(nLoad1>0) call execute_Load1
			if(nHeatGenElm>0) call execute_HeatGen
            Fm_step(:, loadStepLoop) = staVec
            staVec = 0.0
        enddo
        
        a0_ = a0
        if(timeIntType==2) a0_ = -a0
        if(CmatType==0) then
            do i=1, numEquation
                do j=1, numEquation
                    staMatC0(i,j) = a0_*C0mMat(i,j)
                enddo
            enddo
        else if(CmatType==1) then
            do i=1, numEquation               
                staMatC0(i,i) = a0_*C1Mat(i)
            enddo
        end if
        end subroutine execute_heatFlow
		
        subroutine execute_EffectHtFlow
        use mainCtrlInf, only: timeIntType, analyType, nSHnd
        use solCtrlInf, only: klin
if(timeIntType==2) then
    ! 显式积分过程

    
else 
    if((klin/=0).and.(analyType>0)) then
       ! call mult(vecA, matB, vecC, row, column)
    else
        if(nSHnd>0) then
        
        end if
    end if
    
    
    
endif
        
        end subroutine execute_EffectHtFlow
        
        
        subroutine execute_FixTempNd
		use mainCtrlInf, only: numEquation
		use solCtrlInf, only: staVec, staMatK_k, Phi_incNow
		real:: a0_, bigNum
		integer:: i, j, k, fixTempNdID, load1NdID

        do i=1,nFixTempNd
            fixTempNdID = fixTempNdInf(1,i)
            if(loadStepLoop==1) then
                bigNum=1.e10 
            else 
                bigNum=1.0
            end if
            staVec(fixTempNdID) = staMatK_k(fixTempNdID, fixTempNdID)*fixTempNdValue(loadStepLoop,i)*bigNum   
            staMatK_k(fixTempNdID, fixTempNdID) = staMatK_k(fixTempNdID, fixTempNdID)*bigNum                
        enddo
    
        end subroutine execute_FixTempNd
		
		subroutine execute_Load1
		use solCtrlInf, only: C0mMat, C1Mat, staVec
		real:: a0_
		integer:: j, load1NdID
			if(nLoad1>0) then
                do j=1,nLoad1
                    load1NdID = load1Inf(1, j)
                    staVec(load1NdID) = staVec(load1NdID) + load1Value(loadStepLoop,j)
                enddo
            end if
		end subroutine execute_Load1
		
		
		
		subroutine executeConvNdFlow
		
    		
		end subroutine executeConvNdFlow
        
        
        
        
        
    end module heatFlowInf

 