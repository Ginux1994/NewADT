    
!*****************************************************************************    
    module heatFlowInf
    use heatFlowCtrlInf, only: nscr, nFixTempNd, nConvNd, nRadiaNd, nLoad1, nLoad2, nLoad3, nHeatGenElmGrp, nHeatGenElm, maxCurvePoint
	use solCtrlInf, only: iper
	use mainCtrlInf, only: startT, nFixTimeStep
	use timeStepInf
	use timeFuncInf
    implicit none
    integer, allocatable:: fixTempNdInf(:,:), convNdInf(:,:),radiaNdInf(:,:), load1Inf(:,:), heatGenInf(:,:), heatGenElmID(:,:)
	real, allocatable::  fixTempNdValue(:,:), convNdValue(:,:), radiaNdValue(:,:), load1Value(:,:), heatGenValue(:,:)
contains
        subroutine read_heatFlowInf
            if(nFixTempNd>0) call read_fixTempNdFlow
            if(nConvNd>0) call read_ConvNdFlow
            if(nRadiaNd>0) call read_RadiaNdFlow
            if(nLoad1>0) call read_Load1
			if(nHeatGenElm>0) call read_HeatGen
        end subroutine read_heatFlowInf

        subroutine read_fixTempNdFlow
        use heatFlowCtrlInf, only: nCurve, maxCurvePoint
            integer:: i, j, k, l, m, ndID, curveID, fac, kn, pointLoop
			real:: time_now, time_pre, timeTmp, startTime, sloop, arg
            allocate(fixTempNdInf(4, nFixTempNd), fixTempNdValue(stepNumAt(iper), nFixTempNd))! , fixTempNdValue()


            do i=1, nFixTempNd
                read(5,*) ndID, curveID, fac, startTime, kn
                if (fac==0.0) fac = 1.0
                fixTempNdInf(1,i) = ndID; fixTempNdInf(2,i) = curveID; fixTempNdInf(3,i) = fac; fixTempNdInf(4,i) = startTime;
                time_now = startT

                do j=1,stepNumAt(iper)
                    time_pre = time_now + dt_alpha
                    time_now = time_now + dtAt(iper)
                    timeTmp = time_pre - startTime
                    if(timeTmp<=0.0) return
                    pointLoop = 0
                    write(*,*) timeAt(1,1)
                    do l=1,nCurve
                        do m=1,maxCurvePoint
                            write(*,*) timeAt(m,l)
                        enddo
                    enddo
                    do k=2,maxCurvePoint
                        pointLoop = pointLoop + 1
                        if((timeTmp>=timeAt(k-1, curveID)).and.(timeTmp<=timeAt(k, curveID))) exit
                        pause
                    enddo
                    sloop = (timeFuncAt(pointLoop+1, curveID) - timeFuncAt(pointLoop, curveID))/&
                        (timeAt(pointLoop+1, curveID) - timeAt(pointLoop, curveID))
                    arg = timeFuncAt(pointLoop, curveID) + sloop*(timeTmp - timeAt(pointLoop, curveID))
                    fixTempNdValue(j, i) = arg*fac
                enddo

            enddo

                   
        end subroutine read_fixTempNdFlow
    

	
        subroutine read_ConvNdFlow
            integer:: i, j, k, pointLoop, ndID, curveID, fac, startTime, kn
            real:: arg, sloop, time_now, time_pre, timeTmp
            allocate(convNdInf(4,nConvNd), convNdValue(stepNumAt(iper), nConvNd))
            

            do i=1,nConvNd
                read(5,*) ndID, curveID, fac, startTime, kn
                if (fac==0.0) fac = 1.0
                convNdInf(1,i) = ndID; convNdInf(2,i) = curveID; convNdInf(3,i) = fac; convNdInf(4,i) = startTime;

                time_now = startT
                do j=1,stepNumAt(iper)
                    time_pre = time_now + dt_alpha
                    time_now = time_now + dtAt(iper)
                    timeTmp = time_pre - startTime
                    if(timeTmp<=0.0) return
                    pointLoop = 0
                    do k=2,maxCurvePoint
                        pointLoop = pointLoop + 1
                        if((timeTmp>=timeAt(k-1, curveID)).and.(timeTmp<=timeAt(k, curveID))) exit
                        pause
                    enddo
                    sloop = (timeFuncAt(pointLoop+1, curveID) - timeFuncAt(pointLoop, curveID))/&
                        (timeAt(pointLoop+1, curveID) - timeAt(pointLoop, curveID))
                    arg = timeFuncAt(pointLoop, curveID) + sloop*(timeTmp - timeAt(pointLoop, curveID))
                    convNdValue(j, i) = arg*fac
                enddo

            enddo

                      
        end subroutine read_ConvNdFlow
    
        subroutine read_RadiaNdFlow
            integer:: i, k, j, pointLoop
            integer:: ndID, curveID, fac, startTime, kn
            real:: arg, sloop, time_now, time_pre, timeTmp
            allocate(radiaNdInf(4,nRadiaNd), radiaNdValue(stepNumAt(iper), nRadiaNd))


            do i=1,nRadiaNd
                read(5,*) ndID, curveID, fac, startTime, kn
                if (fac==0.0) fac = 1.0
                radiaNdInf(1,i) = ndID; radiaNdInf(2,i) = curveID; radiaNdInf(3,i) = fac; radiaNdInf(4,i) = startTime;

                time_now = startT
                do j=1,stepNumAt(iper)
                    time_pre = time_now + dt_alpha
                    time_now = time_now + dtAt(iper)
                    timeTmp = time_pre - startTime
                    if(timeTmp<=0.0) return
                    pointLoop = 0
                    do k=2,maxCurvePoint
                        pointLoop = pointLoop + 1
                        if((timeTmp>=timeAt(k-1, curveID)).and.(timeTmp<=timeAt(k, curveID))) exit
                        pause
                    enddo
                    sloop = (timeFuncAt(pointLoop+1, curveID) - timeFuncAt(pointLoop, curveID))/&
                        (timeAt(pointLoop+1, curveID) - timeAt(pointLoop, curveID))
                    arg = timeFuncAt(pointLoop, curveID) + sloop*(timeTmp - timeAt(pointLoop, curveID))
                    radiaNdValue(j, i) = arg*fac
                enddo

            enddo

                    
        end subroutine read_RadiaNdFlow
    
        subroutine read_Load1
            integer:: i, j, k, pointLoop, ndID, curveID, fac, startTime, kn
            real:: arg, sloop, time_now, time_pre, timeTmp
            allocate(load1Inf(4,nLoad1), load1Value(stepNumAt(iper), nLoad1))
            
                if (nLoad1>0) then 
                    do i=1, nLoad1
                        read(5,*) ndID, curveID, fac, startTime, kn
                        if (fac==0.0) fac = 1.0
                        load1Inf(1,i) = ndID; load1Inf(2,i) = curveID; load1Inf(3,i) = fac; load1Inf(4,i) = startTime;
						
						time_now = startT						
						do j=1,stepNumAt(iper)
							time_pre = time_now + dt_alpha
							time_now = time_now + dtAt(iper)
							timeTmp = time_pre - startTime
							if(timeTmp<=0.0) return
							pointLoop = 0
							do k=2,maxCurvePoint
								pointLoop = pointLoop + 1
								if((timeTmp>=timeAt(k-1, curveID)).and.(timeTmp<=timeAt(k, curveID))) exit
								pause
							enddo						
							sloop = (timeFuncAt(pointLoop+1, curveID) - timeFuncAt(pointLoop, curveID))/&
                                    (timeAt(pointLoop+1, curveID) - timeAt(pointLoop, curveID))
							arg = timeFuncAt(pointLoop, curveID) + sloop*(timeTmp - timeAt(pointLoop, curveID))
							load1Value(j, i) = arg*fac
						enddo
						
						
                    enddo
                end if
               
                
        end subroutine read_Load1
        
		subroutine read_HeatGen		
        use timeStepInf
		    integer:: i, j, k, l, pointLoop,&
                      linearType, grpID, nElm,&
                      elmID, curveID, fac, startTime, kn
            integer:: stepNumNow
            real:: arg, sloop, time_now, time_pre, timeTmp
            stepNumNow = stepNumAt(iper)
            allocate(heatGenInf(4,nHeatGenElmGrp), heatGenElmID(nHeatGenElm, nHeatGenElmGrp), heatGenValue(stepNumNow, nHeatGenElm))             
            do i=1, nHeatGenElmGrp
				
                read(5,*) linearType, grpID, nElm
				heatGenInf(1,i) = linearType; heatGenInf(2,i) = grpID; heatGenInf(3,i) = nElm;
				
				do j=1, nElm
					read(5,*) elmID, curveID, fac, startTime, kn
                	heatGenElmID(j, i) = elmID	
					if (fac==0.0) fac = 1.0
					time_now = startT
					do l=1,stepNumAt(iper)
						time_pre = time_now + dt_alpha
						time_now = time_now + dtAt(iper)
						timeTmp = time_pre - startTime
						if(timeTmp<=0.0) return
						pointLoop = 0
						do k=2,maxCurvePoint
							pointLoop = pointLoop + 1
							if((timeTmp>=timeAt(k-1, curveID)).and.(timeTmp<=timeAt(k, curveID))) exit
							pause
						enddo
						sloop = (timeFuncAt(pointLoop+1, curveID) - timeFuncAt(pointLoop, curveID))/&
							(timeAt(pointLoop+1, curveID) - timeAt(pointLoop, curveID))
						arg = timeFuncAt(pointLoop, curveID) + sloop*(timeTmp - timeAt(pointLoop, curveID))
						heatGenValue(j, i) = arg*fac ! 热生成数据按照单元顺序存储，其相应的单元编号通过heatGenElmID(j, i) = elmID	获取
                    enddo
                enddo
            enddo

		end subroutine read_HeatGen
		

		
		
        subroutine read_Load2
        
        
        
        end subroutine read_Load2
                
        subroutine read_Load3
        
        
        
        end subroutine read_Load3
        
        
            
        subroutine execute_heatFlow
        use element, only: execute_HeatGen
            ! if(nConvNd>0) ! call execute_ConvNdFlow
            ! if(nRadiaNd>0) ! call execute_RadiaNdFlow
            if(nLoad1>0) call execute_Load1
			if(nHeatGenElm>0) call execute_HeatGen
		
        end subroutine execute_heatFlow
		
        
        
		
		subroutine execute_Load1
		use timeStepInf, only: a0
		use mainCtrlInf, only: numEquation, CmatType, timeIntType
		use solCtrlInf, only: Fm, KmMat, C0mMat, C1Mat
		real:: r(nFixTempNd), a0_
		integer:: i, j, k, fixTempNdID, load1NdID
			if(nFixTempNd>0) then
				do i=1,nFixTempNd
					fixTempNdID = fixTempNdInf(1,i)
					do j=1,nLoad1
						load1NdID = load1Inf(1, j)
						Fm(load1NdID) = Fm(load1NdID) + load1Value(iper,j)
						Fm(fixTempNdID) = Fm(fixTempNdID) + fixTempNdValue(iper, fixTempNdID)*1.e10
					enddo
				enddo 
            end if
            a0_ = a0
			if(timeIntType==2) a0_ = -a0
			if(CmatType==0) then 
				do i=1, numEquation
					do j=1, numEquation
						KmMat(i,j) = a0_*C0mMat(i,j)
					enddo
				enddo 
			else if(CmatType==1) then
				do i=1, numEquation
					KmMat(i,i) = a0_*C1Mat(i)
                    do j=1,numEquation
                        write(*,"(F9.3)") KmMat(i,j)
                    enddo
                    write(*,*)
				enddo
			end if
		end subroutine execute_Load1
		
		
		
		subroutine executeConvNdFlow
		
    		
		end subroutine executeConvNdFlow
		
		

		
        
        
    end module heatFlowInf
