   
!****************************************************************************    
    module Gauss
    real::gaussPoint(4,4), gaussWeight(4,4)
    data gaussPoint / 0.d0,   0.d0,   0.d0,   0.0,&
      -.5773502691896d0,  .5773502691896d0,0.0,  0.d0,&
      -.7745966692415d0,  .0000000000000d0, .7745966692415d0,0.0,&
      -.8611363115941d0,  -.3399810435849d0, .3399810435849d0,    .8611363115941d0/      
    data gaussWeight /2.0,0.0,0.0,0.0,&
                   1.0,1.0,0.0,0.0,&
       .5555555555556d0,.8888888888889d0,.5555555555556d0, 0.0d0,&
       .3478548451375d0, .6521451548625d0, .6521451548625d0,.3478548451375d0/
    end module
    
    module solCtrlInf
    integer:: index, iper, numLinearEqua, numNonLnEqua, numTotalEqua
    integer:: nLrElm, nNonLrElm, klin
    integer:: icount  
    integer:: numMasterDOF, numSlaveDOF
    real, allocatable:: coords(:,:)
    real, allocatable:: KsMat(:, :), KmMat(:, :), KsmMat(:, :), C0mat(:, :), C1mat(:)
    real, allocatable:: Fs(:), Fm(:)    
    integer, allocatable:: DOFMap(:)!, ndStatus(:,1)

    character(len=5):: ind
contains    
    subroutine init
        allocate(KmMat(numMasterDOF, numMasterDOF))
        numTotalEqua = numLinearEqua + numNonLnEqua
    end subroutine init
    
    end module
    

    
    module var
    integer:: iref
    end module

    
    
!********************************************************************       
    module mainCtrlInf ! 2.������Ϣ
    ! card 2.1
    ! �ڵ����������Ե�Ԫ�����������Ե�Ԫ������ ������ͣ���ʼʱ�䣬
    integer:: nNd, nLrElmGrp, nNonlrElmGrp, calType, nFixTimeStep, startT
    integer:: nElmGrp, numEquation
    ! card 2.2
    ! �������ͣ����нڵ����������������
    integer:: analyType, nSHnd, nLatht
    integer:: CmatType=0
    ! card 2.3
    ! solve the eigen or not
    integer:: eigenIndex
    ! card 2.4
    ! �����γɾ���֮���ʱ�䲽����ƽ�����֮���ʱ�䲽����ƽ��������������������
    integer:: nReformStep, nEquitTStep,nMaxEqStep
    real:: tolerance
    ! card 2.5
    ! ʱ��������ͣ���
    integer:: timeIntType, alpha
    ! card 2.6   
    integer:: NPB, IPRLH
    integer:: haveSHnode
    integer, allocatable:: npar(:,:)
contains      
        subroutine read_mainCtrlInf
            read(5,*) nNd, nLrElmGrp, nNonlrElmGrp, calType, nFixTimeStep, startT
            nElmGrp = nLrElmGrp + nNonlrElmGrp
            allocate(npar(20,nElmGrp))
            if ((nFixTimeStep==0).or.(calType==1)) nFixTimeStep = 1
            read(5,*) analyType, nSHnd, nLatht
            read(5,*) eigenIndex
            read(5,*) nReformStep, nEquitTStep, nMaxEqStep, tolerance
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
                CmatType = 2
            end if
        end subroutine read_mainCtrlInf
    end module mainCtrlInf
    
!********************************************************************       
    module timeStepInf ! 3.ʱ�䲽����Ϣ
    ! ��i�����ʱ��ε�ʱ�䲽����������������
    integer, allocatable:: stepNumAt(:)
    integer:: convergenFlag, nMaxStep
    real, allocatable:: dtAt(:), a0
contains     
        subroutine read_timeStepInf
            use mainCtrlInf, only: timeIntType, alpha, nFixTimeStep
            integer:: i
            allocate(stepNumAt(nFixTimeStep), dtAt(nFixTimeStep))
            read(5,*)
            read(5,*) (stepNumAt(i), i=1, nFixTimeStep)
            read(5,*) (dtAt(i), i=1, nFixTimeStep)
            nMaxStep = 0
            do i=1, nFixTimeStep
                if(stepNumAt(i)>nMaxStep) nMaxStep = stepNumAt(i)
            enddo
            
            if(timeIntType==4) then
                if(alpha==0.0) a0 = 1.0/dtAt(0)
                if(alpha>0.0) a0 = 1.0/(dtAt(0)*alpha)
            end if           
        end subroutine read_timeStepInf
    end module timeStepInf
!********************************************************************    
    module nodeInf ! 4.�ڵ���Ϣ
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
    module ndSHInf ! 5.���нڵ�������Ϣ
    ! n�Žڵ�ļ��нڵ�����
    real, allocatable:: ndSH(:), dql(:,:), qlv(:,:) !ndSH(nSHnd), dql(numEquation,nSHnd)
contains
    subroutine read_ndSHInf
        use mainCtrlInf, only: nSHnd, nLatht, numEquation
        use nodeInf, only: dofID, ndDOF
        integer:: i, j, nodeID, id
        real:: nodeSH

        if (nSHnd>0) then
            allocate(ndSH(numEquation), dql(numEquation,nLatht), qlv(numEquation,nLatht))
            ndSH = 0.0
            do i=1,nSHnd
                read(5,*) nodeID, nodeSH
                do j=1,ndDOF
                    id = dofID(j, nodeID)
                    if(id.ne.0) ndSH(id) = nodeSH
                enddo
            enddo
        end if               
    end subroutine read_ndSHInf
    end module ndSHInf
!********************************************************************    
    
    module initInf ! 6.��ʼ��Ϣ
    ! 0-�Զ��γ�0��ʼ��������ʼ�¶�����
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
    module heatFlowCtrlInf ! 7.�������ƿ� 
    ! 7.1 control card
    ! ������Ŀ�����ߵ������Ŀ���������ͣ��̶��ڵ��¶���Ŀ�������ڵ���������ڵ���Ŀ������������Ŀ����ֲ�������Ŀ����ֲ�������Ŀ
    integer:: nCurve, maxCurvePoint, convType, nFixTempNd, nConvNd, nRadiaNd, nLoad1, nLoad2, nLoad3, nHeatGenElmGrp, nHeatGenElm
    ! 7.2 phase transform temp�������¶ȣ�����������ʼ����־
    real, allocatable:: phaseTransfTemp(:), phaseTransfDt(:)
    integer, allocatable::initPhaseFlag(:)
    integer:: haveLoad=0, haveLnConvOrFixTempLoad=0, nscr   
contains
        subroutine read_heatFlowCtrlInf
            use mainCtrlInf, only: nLatht
            allocate(phaseTransfTemp(nLatht), phaseTransfDt(nLatht), initPhaseFlag(nLatht))
            read(5,*) nCurve, maxCurvePoint, convType, nFixTempNd, nConvNd, nRadiaNd, nLoad1, nLoad2, nLoad3, nHeatGenElmGrp, nHeatGenElm
            if(nLatht>0) then
                read(5,*) (phaseTransfTemp(i), i=1, nLatht)
                read(5,*) (phaseTransfDt(i), i=1, nLatht)
                read(5,*) (initPhaseFlag(i), i=1, nLatht) 
            end if
            if((nConvNd + nRadiaNd + nLoad1 + nLoad2 + nLoad3)>0) haveLoad=1
            if((convType==1).or.(nFixTempNd>0)) haveLnConvOrFixTempLoad=1 !convType==1��ʾ���Զ���
            nscr = nFixTempNd + nConvNd + nRadiaNd
        end subroutine read_heatFlowCtrlInf
    end module
!*****************************************************************************    
   
    module timeFuncInf ! 8. ʱ�亯��������Ϣ�������ⲿ�������Լ����ϵ�
    real, allocatable:: timeAt(:,:), timeFuncAt(:,:)  
contains
        subroutine read_timeFuncInf
            use heatFlowCtrlInf, only: nCurve, maxCurvePoint
            integer:: curveID, nCurvePoints, i, j
        
            allocate(timeAt(nCurve, maxCurvePoint), timeFuncAt(nCurve, maxCurvePoint))
            do i=1,nCurve   
                read(5,*) curveID, nCurvePoints                 
                read(5,*) (timeAt(curveID, j), timeFuncAt(curveID, j), j=1,nCurvePoints)
            end do
        end subroutine read_timeFuncInf
    end module
    
    
!*****************************************************************************    
    module heatFlowInf
    use heatFlowCtrlInf, only: nscr
    real, allocatable:: fixTempNdInf(:,:), convNdInf(:,:), radiaNdInf(:,:)
contains
        subroutine read_heatFlowInf
            call read_fixTempNdFlow
            call read_ConvNdFlow
            call read_RadiaNdFlow
        end subroutine read_heatFlowInf
    
        subroutine read_fixTempNdFlow
            use solCtrlInf, only: iper
            use heatFlowCtrlInf, only: nFixTempNd, nConvNd, nRadiaNd
            integer:: i, ndID, curveID, fac, startTime, kn
            allocate(fixTempNdInf(4,nFixTempNd))
            if (iper==1) then 
                if (nFixTempNd>0) then 
                    do i=1, nFixTempNd
                        read(5,*) ndID, curveID, fac, startTime, kn
                        if (fac==0.0) fac = 1.0
                        fixTempNdInf(1,i) = ndID; fixTempNdInf(2,i) = curveID; fixTempNdInf(3,i) = fac; fixTempNdInf(4,i) = startTime;
                    enddo
                end if
            end if           
        end subroutine read_fixTempNdFlow
    
        subroutine read_ConvNdFlow
            use solCtrlInf, only: iper
            use heatFlowCtrlInf, only: nConvNd
            integer:: i, ndID, curveID, fac, startTime, kn
            allocate(convNdInf(nConvNd,4))
            if (iper==1) then 
                if (nConvNd>0) then    
                    do i=1,nConvNd
                        read(5,*) ndID, curveID, fac, startTime, kn
                        if (fac==0.0) fac = 1.0
                        convNdInf(i,1) = ndID; convNdInf(i,2) = curveID; convNdInf(i,3) = fac; convNdInf(i,4) = startTime;
                    enddo
                end if
            end if           
        end subroutine read_ConvNdFlow
    
        subroutine read_RadiaNdFlow
            use solCtrlInf, only: iper
            use heatFlowCtrlInf, only: nRadiaNd
            integer:: ndID, curveID, fac, startTime, kn
            allocate(radiaNdInf(nRadiaNd,4))
            if (iper==1) then 
                if (nConvNd>0) then      
                    do i=1,nRadiaNd
                        read(5,*) ndID, curveID, fac, startTime, kn
                        if (fac==0.0) fac = 1.0
                        radiaNdInf(i,1) = ndID; radiaNdInf(i,2) = curveID; radiaNdInf(i,3) = fac; radiaNdInf(i,4) = startTime;
                    enddo
                end if
            end if           
        end subroutine read_RadiaNdFlow
    
        
    end module heatFlowInf

!*****************************************************************************    

    module thrdCondElm
! npar information    
    integer:: nthrdCondElm, nonLrMat, deathType, nMaxNd, nGauP_rs, nGauP_t, matType, nMatGrp, nMatK, nMatC
! 3d conduction element information
    integer, allocatable:: nElmNd(:), nGeomNd(:), elmMatID(:), elmNdID(:,:)
    integer, allocatable:: lm(:,:)
    real, allocatable:: propK(:,:), propC(:,:), rlh(:,:), specHeat, deaTime(:)
    real, allocatable:: xyz(:,:)
! variable for assem element
     
contains
        subroutine read_thrdCondElm(elmGrpID)
            use mainCtrlInf, only: nElmGrp, npar, nLatht, nNd
            use nodeInf, only: dofID, x, y, z
            integer:: matID, elmID, i, j, k, ndID, idLoop=1, ii, jj, kk
            integer:: elmGrpID
            nthrdCondElm = npar(2,elmGrpID); nonLrMat=npar(3,elmGrpID); deathType=npar(4,elmGrpID); 
            nMaxNd=npar(7,elmGrpID); nGauP_rs=npar(10,elmGrpID); nGauP_t=npar(11,elmGrpID); 
            matType=npar(15,elmGrpID); nMatGrp=npar(16,elmGrpID); nMatK=npar(17,elmGrpID); 
            nMatC=npar(18,elmGrpID);
        
            allocate(propK(nMatK,nMatGrp), propC(nMatC,nMatGrp), rlh(nLatht,nMatGrp))
            allocate(nElmNd(nthrdCondElm), nGeomNd(nthrdCondElm), elmMatID(nthrdCondElm), deaTime(nthrdCondElm))
            allocate(elmNdID(21,nthrdCondElm), xyz(3,nNd), lm(nMaxNd, nthrdCondElm))
    ! read material inf
            do i=1,nMatGrp
                read(5,*) matID
                if(nLatht>0) read(5,*) (rlh(j,matID), j=1,nLatht)
                read(5,*) (propK(j,matID), j=1,nMatK)
                read(5,*) (propC(j,matID), j=1,nMatC)
            enddo   
    ! read element inf
            do i=1,nthrdCondElm ! �����д�����Ԫѭ�������Լ������Ե�
                read(5,*) elmID, nElmNd(elmID), nGeomNd(elmID), dummy2, elmMatID(elmID), dummy2, kg, deaTime(elmID)
                if(nElmNd(elmID)==0) nElmNd(elmID) = nMaxNd
                if(nGeomNd(elmID)==0) nGeomNd(elmID) = nMaxNd
                read(5,*) (elmNdID(j,nthrdCondElm), j=1,8)
                read(5,*) (elmNdID(j,nthrdCondElm), j=9,21)
                idLoop = 1
                do k=1,nElmNd(elmID) ! ��i�ŵ�Ԫ���нڵ�ѭ��
                    ndID = elmNdID(k,nthrdCondElm)
                    xyz(1, ndID) = x(ndID); xyz(2, ndID) = Y(ndID); xyz(3, ndID) = Z(ndID); 
                    lm(idLoop, elmID) = dofID(1, ndID)
                    idLoop = idLoop + 1
                enddo   
            enddo
        
        end subroutine read_thrdCondElm
    
        subroutine assem_thrdCondElm
        use solCtrlInf, only: index, KmMat
        use nodeInf, only:dofID
            integer:: elmID, elmNdNum, ndID(8)
            real:: coords(3,8)
            real:: tempNow(8), tempInc(8), res(8), N(8), dNdxi(8), Bmat(3,8), Dmat(3,3), Cmat0(8,8), Cmat1(8), Kmat(8,8)
            if (index==1) then
                ! assem linear conductivity matrix
                do elmID=1, nthrdCondElm
                    elmNdNum = nElmNd(elmID)
                    nod9 = elmNdNum - 8
                    call getCoordsByElmID(elmID, elmNdNum, coords)
                    call cal3DCondMat(elmID, elmNdNum, nod9, coords, propK, tempNow, Kmat, res)
                    call getNdIDByElmID(elmID, 8, ndID)
                    call assemElmMat(Kmat, ndID)
                enddo
            else if((index == 2).or. (index==3))then 
                ! assem linear capacity matrix
                do elmID=1, nthrdCondElm
                    elmNdNum = nElmNd(elmID)
                    nod9 = elmNdNum - 8
                    call getCoordsByElmID(elmID, elmNdNum, coords)
                    call cal3DHeatCapMat(elmID, elmNdNum, nod9, coords, propC, tempNow, tempInc, Cmat0, CMat1, res)
                    call getNdIDByElmID(elmID, 8, ndID)
                    if(CmatType==0) then
                        call assemElmMat(Cmat0, ndID)
                    else if(CmatType==1) then
                        do i=1,8
                        j = dofID(1, ndID(i))
                        KmMat(j,j) = KmMat(j,j) + Cmat1(i)
                        end do
                    end if
                enddo
                !call
            end if       
         end subroutine assem_thrdCondElm
    
        subroutine getCoordsByElmID(elmID, nElmNd, coords)
            use nodeInf, only: x, y, z
            integer, intent(in):: elmID, nElmNd
            real, intent(out):: coords(3,8)
            integer:: ndID(8), i
            call getNdIDByElmID(elmID, nElmNd, ndID)
            do i=1,8
                coords(1,i) = x(ndID(i))
                coords(2,i) = y(ndID(i))
                coords(3,i) = z(ndID(i))
            enddo    
        end subroutine getCoordsByElmID
    
        subroutine getNdIDByElmID(elmID, nElmNd, ndID)
            integer, intent(out):: ndID(8)
            integer, intent(in):: elmID, nElmNd
            integer:: i
            do i=1,nElmNd
                ndID(i) = elmNdID(i,elmID)
            enddo       
        end subroutine getNdIDByElmID
    
    end module thrdCondElm
!*****************************************************************************     
    
    module element
    use thrdCondElm
contains
        subroutine read_element
            use mainCtrlInf, only: nElmGrp, npar
            use thrdCondElm
            integer:: elmGrpID, nparTemp(20)
            do elmGrpID=1,nElmGrp
                nparTemp = 0
                read(5,*) nparTemp
                select case(nparTemp(1))
                case(3)
                    npar(:, elmGrpID) = nparTemp
                    call read_thrdCondElm(elmGrpID)
                end select
            enddo
        end subroutine read_element
   
        subroutine exceuteAssem(elmType)
        integer:: elmType
            select case(elmType)
            case(1)
                    ! call assem_onedCondElm
            case(2)
                    ! call assem_twodCondElm
            case(3)
                    call assem_thrdCondElm
            case(4)
                    ! call assem_convElm
            case(5)
                    ! call assem_radiaElm
            end select
        end subroutine exceuteAssem
        
        subroutine assem_element
            use mainCtrlInf, only: nNd, nSHnd, nLatht, numEquation, npar, nLrElmGrp, nNonLrElmGrp, nElmGrp, &
                CmatType
            use solCtrlInf, only: index, DOFMap, numMasterDOF, &
            numSlaveDOF, KmMat, KsMat, KsmMat , numLinearEqua, numNonLnEqua, numTotalEqua  
            use heatFlowCtrlInf, only: haveLnConvOrFixTempLoad, nFixTempNd
            use heatFlowInf, only:fixTempNdInf
        
            implicit none
            integer:: i, j, k, elmGrpLoop, fileID
if(index==1) then           
            if(haveLnConvOrFixTempLoad==1) then ! ��������Զ�����̶��¶ȣ�����assem linear matrix
        
                do elmGrpLoop=1,nElmGrp

                    if (npar(1, elmGrpLoop)==5) cycle ! radiation element, ignore             
                    call exceuteAssem(npar(1, elmGrpLoop))

                end do
                ! deal withh the fix temperature node
                do i=1,nFixTempNd
                    j = fixTempNdInf(1, i)
                    KmMat(j,j) = KmMat(j,j) + 1.0e10
                end do
            
                index=2       
            else if((nSHnd>0).and.(nLatht.ne.0)) then ! ��װ���ݾ���
 
                do elmGrpLoop=1,nElmGrp
                
                    if (npar(1, elmGrpLoop)>3) cycle ! radiation element, ignore               
                    call exceuteAssem(npar(1, elmGrpLoop))
                        ! �������ݼ���󣬲�����ӵ��ն����ϣ�������Ϊ ���غɼ���
                end do
            end if
else if(index==3) then
            if(nLatht>0) then
                ! ��װǱ������
                    
            end if  
else if(index==4) then
            
            
end if            
            end subroutine assem_element
    end module element
!*****************************************************************************       
    
    
    
    module onedCondElm   
    end module onedCondElm
            
    module twodCondElm    
    end module twodCondElm
    
    module convElm
    end module convElm
    
    module radiaElm
    end module radiaElm
    
!*****************************************************************************    

     
!*****************************************************************************    
    subroutine cal3DHeatCapMat(elmID, ndNum, nod9, coords, propC, tempNow, tempInc, Cmat0, CMat1, res)
    use solCtrlInf, only: index, icount
    use mainCtrlInf, only:CmatType
    use Gauss
    use var
    !* param table of the function()
    integer, intent(in)::elmID, ndNum, nod9
    real, intent(in)::coords(3,8), propC(1), tempNow(8), tempInc(8)
    real, intent(out):: Cmat0(8,8), Cmat1(8), res(8)
    
    !* param table in the function
    ! 1. shap function:
    real:: N(8), dNdxi(3, 8), Jmat(3,3), Jinv(3,3), Jdet
    ! 2. gauss integration loop variable
    real::r, s, t, weightTotal, fac
    integer::  xLoop, yLoop, zLoop, i, j, k, numGaussP = 2
    ! 3. D K matrix, and the temporary variable
    real:: Dmat(3,3), Bmat(3,8), elmTemp = 0.0, averageTempInc = 0.0
    
if(index<4) then !linear analysis or nonlinear
    
    ! calculate linear C matrix
    do xLoop=1,numGaussP
         r = gaussPoint(xLoop, numGaussP)
        do yLoop = 1,numGaussP
            s = gaussPoint(yLoop, numGaussP)
            do zLoop = 1,numGaussP
                t = gaussPoint(zLoop, numGaussP)
                weightTotal = gaussWeight(xLoop, numGaussP)*gaussWeight(yLoop, numGaussP)*gaussWeight(zLoop, numGaussP)
                call caldNdx(r, s, t, nod9, elmID, coords, N, dNdxi, Jmat, Jinv, Jdet)
                fac = weightTotal*Jdet*specHeat
                
                if(CmatType>=2) then
                   ! call initRmat(Cmat0,8,8)
                    do i=1,8
                        do j=1,8
                            Cmat0(i,j) = Cmat0(i,j) + N(i)*N(j)*fac
                        enddo
                    enddo
                else 
                    fac = fac/8
                   ! call initRmat(Cmat1, 8 ,8)
                    do i=1,8
                        Cmat1(i) = fac
                    enddo
                end if              
            enddo   
        enddo
    enddo

else 
    
    ! calculate non linear C matrix
    do xLoop=1,numGaussP
         r = gaussPoint(xLoop, numGaussP)
        do yLoop = 1,numGaussP
            s = gaussPoint(yLoop, numGaussP)
            do zLoop = 1,numGaussP
                t = gaussPoint(zLoop, numGaussP)
                weightTotal = gaussWeight(xLoop, numGaussP)*gaussWeight(yLoop, numGaussP)*gaussWeight(zLoop, numGaussP)
                call caldNdx(r, s, t, nod9, elmID, coords, N, dNdxi, Jmat, Jinv, Jdet)
                
                do i=1,8
                    elmTemp = elmTemp + N(i)*tempNow(i)
                enddo
                ! get the nonlinear C-specHeat at the temperature of elmTemp
                fac = weightTotal*Jdet*specHeat
                
                if(CmatType>=2) then
                    ! calculate the consistent matrix
                    if(icount>=2) then
                        ! calculate the right hand
                        do i=1,8
                            averageTempInc = averageTempInc + N(i)*tempInc(i)
                        enddo
                        do i=1,8
                            res(i) = res(i) + averageTempInc*N(i)*fac
                        enddo
                        continue
                    end if 
                    
                    if(iref == 0) then
                        do i=1,8
                            do j=1,8
                            Cmat0(i,j) = Cmat0(i,j) + N(i)*N(j)*fac
                            enddo
                        enddo
                        continue
                    end if
                else 
                    ! calcualte the lumped matrix
                    fac = fac/8
                    if(iope==2) then
                        do i=1,8
                            Cmat1(i) = Cmat1(i) + fac
                            res(i) = res(i) - fac*tempNow(i)
                        enddo                       
                        continue
                    else 
                        if(icount<=2) then
                            do i=1,8
                                res(i) = res(i) + fac*tempInc(i)
                            enddo
                            continue
                        end if
                        if(iref==0) then
                            do i=1,8
                                Cmat0(i,i) = Cmat0(i,i) + fac
                            enddo
                            continue
                        end if 
                     
                    end if ! if(iope==2) then  
                end if ! if(CmatType>=2) then
                
            enddo ! do zLoop = 1,numGaussP   
        enddo !do yLoop = 1,numGaussP
    enddo !do xLoop=1,numGaussP
    
end if !(index<4) then !linear analysis or nonlinear
       
    end subroutine cal3DHeatCapMat
   
    
    
 !*****************************************************************************    
    subroutine cal3DCondMat(elmID, ndNum, nod9, coords, propK, tempNow, Kmat, res)
    use Gauss
    use solCtrlInf, only: index, icount
    use mainCtrlInf, only: timeIntType
    use var
    implicit none
    !* param table of the function()
    integer, intent(in)::elmID, ndNum, nod9
    real, intent(in)::coords(3,8), propK(1), tempNow(8)
    real, intent(out)::Kmat(8,8), res(8)
    
    !* param table in the function
    ! 1. shap function:
    real:: N(8), dNdxi(3, 8), Jmat(3,3), Jinv(3,3), Jdet
    ! 2. gauss integration loop variable
    real::r, s, t, weightTotal, fac
    integer::  xLoop, yLoop, zLoop, i, j, k, numGaussP = 2
    ! 3. D K matrix, and the temporary variable
    real:: Dmat(3,3), Bmat(3,8), Ktemp(8,8), elmTemp = 0.0
    
    Dmat = 0.0 ! call initRmat(Dmat, 3, 3)
    Kmat = 0.0 ! call initRmat(Kmat, 8, 8)
    Ktemp = 0.0 ! call initRmat(Ktemp, 8, 8)
    

if (index<4) then !calculate linear element's conductivity matrix
    ! initial D matrix

    Dmat(1,1) = propK(1);Dmat(2,2) = propK(1);Dmat(3,3) = propK(1); !conduction coefficient matrix
    
    ! calculate linear K matrix
    do xLoop=1,numGaussP
         r = gaussPoint(xLoop, numGaussP)
        do yLoop = 1,numGaussP
            s = gaussPoint(yLoop, numGaussP)
            do zLoop = 1,numGaussP
                t = gaussPoint(zLoop, numGaussP)
                weightTotal = gaussWeight(xLoop, numGaussP)*gaussWeight(yLoop, numGaussP)*gaussWeight(zLoop, numGaussP)
                call calBmat(r, s, t, nod9, elmID, coords, N, dNdxi, Jmat, Jinv, Jdet, Bmat)
                fac = weightTotal*Jdet
    
                Kmat=Kmat + matmul(matmul(transpose(Bmat),Dmat),Bmat)*fac ! B^T*D*B*(|J|*w)             
                
            enddo   
        enddo
    enddo
else ! calculate the nonlinear element matrix
    
    Dmat(1,1) = propK(1);Dmat(2,2) = propK(2);Dmat(3,3) = propK(3);
    do xLoop=1,numGaussP
        r = gaussPoint(xLoop, numGaussP)
        do yLoop = 1,numGaussP
            s = gaussPoint(yLoop, numGaussP)
            do zLoop = 1,numGaussP
                t = gaussPoint(zLoop, numGaussP)
                weightTotal = gaussWeight(xLoop, numGaussP)*gaussWeight(yLoop, numGaussP)*gaussWeight(zLoop, numGaussP)
                call calBmat(r, s, t, nod9, elmID, coords, N, dNdxi, Jmat, Jinv, Jdet, Bmat)
                fac = weightTotal*Jdet
    
                do i=1, ndNum
                    elmTemp = elmTemp + N(i)*tempNow(i)
                enddo
                ! calculate the nonlinear D matrix here
                
                ! calculate the heat flow here               
                Ktemp=Ktemp + matmul(matmul(transpose(Bmat),Dmat),Bmat)*fac ! B^T*D*B*(|J|*w)             
                res = res + matmul(Ktemp,tempNow)
                if(timeIntType==2) cycle
                if((icount <= 2).and.(iref==0)) then
                    Kmat = Kmat + Ktemp ! rebuild k matrix
                end if
                    
            enddo   
        enddo
    enddo
        
end if
    end subroutine cal3DCondMat
   
    
    
!*****************************************************************************    
    subroutine cal3DHeatGenVec(elmID, ndID, ndNum, nod9, coords, heatGen, res)
    use Gauss !, only:gaussPoint,gaussWeight
    implicit none
    
    real, intent(in)::coords(3,8), heatGen
    integer, intent(in)::elmID, ndID, nod9, ndNum
    real, intent(out):: res(:)
    
    ! loop variable
    integer:: xLoop, yLoop, zLoop, i, j, k
    ! gauss integration data
    integer:: numGaussP = 2
    real::r, s, t, weightTotal, fac
    ! shape function
    real:: N(8), dNdxi(3,8), Jmat(3,3), Jinv(3,3), Jdet
    
    res = 0.0 !call initRvec(res, ndNum)
    do xLoop=1,numGaussP
         r = gaussPoint(xLoop, numGaussP)
        do yLoop = 1,numGaussP
            s = gaussPoint(yLoop, numGaussP)
            do zLoop = 1,numGaussP
                t = gaussPoint(zLoop, numGaussP)
                weightTotal = gaussWeight(xLoop, numGaussP)*gaussWeight(yLoop, numGaussP)*gaussWeight(zLoop, numGaussP)
         
                call caldNdx(r, s, t, nod9, elmID, coords, N, dNdxi, Jmat, Jinv, Jdet)
                fac = weightTotal*Jdet*heatGen
                do k = 1,ndNum
                    res(k) = res(k) + N(k)*fac
                enddo
            enddo
        enddo
    enddo
    
    return
    end subroutine cal3DHeatGenVec
                
      
    
 !*****************************************************************************    
    subroutine calBmat(r, s, t, nod9, elmID, coords, N, dNdxi, Jmat, Jinv, Jdet, Bmat)

    implicit none
    integer:: i,j,k !loop param
    ! declare the type and dimension of the function parameter
    real, intent(out):: Bmat(3,8)
    real, intent(out):: N(8), dNdxi(3,8), Jmat(3,3), Jinv(3,3), Jdet
    real, intent(in)::r, s, t, coords(3,8)
    integer:: nod9, elmID
    real::tempDouble
    
    call caldNdx(r, s, t, nod9, elmID, coords, N, dNdxi, Jmat, Jinv, Jdet)
       
    Bmat = 0.0 !call initRmat(Bmat, 3, 8)
    Bmat = Bmat + matmul(Jinv, dNdxi)
    
    return
    end subroutine calBmat
    


!***************************************************************************** 
    subroutine caldNdx(r, s, t, nod9, elmID, coords, N, dNdxi, Jmat, Jinv, Jdet)
    ! calculate the shap function and derivatives.
    implicit none
    integer:: i,j,k !loop param
    ! declare the type and dimension of the function parameter
    real, intent(out):: N(8), dNdxi(3,8), Jmat(3,3), Jinv(3,3), Jdet
    real, intent(in)::r, s, t, coords(3,8)
    integer:: nod9, elmID
    ! declare the temporal variable of this fuction
    real::rPlus, sPlus, tPlus, rMinus, sMinus, tMinus, rSqure, sSqure, tSqure
    real::tempDouble
    rPlus=1.0 + r
    sPlus=1.0 + s
    tPlus=1.0 + t
    rMinus=1.0 - r
    sMinus=1.0 - s
    tMinus=1.0 - t
    rSqure=1.0 - r*r
    sSqure=1.0 - s*s
    tSqure=1.0 - t*t

    N(1)=0.125*rPlus*sPlus*tPlus
    N(2)=0.125*rMinus*sPlus*tPlus
    N(3)=0.125*rMinus*sMinus*tPlus
    N(4)=0.125*rPlus*sMinus*tPlus
    N(5)=0.125*rPlus*sPlus*tMinus
    N(6)=0.125*rMinus*sPlus*tMinus
    N(7)=0.125*rMinus*sMinus*tMinus
    N(8)=0.125*rPlus*sMinus*tMinus

    dNdxi(1,1)= 0.125*sPlus*tPlus
    dNdxi(1,2)=-dNdxi(1,1)
    dNdxi(1,3)=-0.125*sMinus*tPlus
    dNdxi(1,4)=-dNdxi(1,3)
    dNdxi(1,5)= 0.125*sPlus*tMinus
    dNdxi(1,6)=-dNdxi(1,5)
    dNdxi(1,7)=-0.125*sMinus*tMinus
    dNdxi(1,8)=-dNdxi(1,7)

    dNdxi(2,1)= 0.125*rPlus*tPlus
    dNdxi(2,2)= 0.125*rMinus*tPlus
    dNdxi(2,3)=-dNdxi(2,2)
    dNdxi(2,4)=-dNdxi(2,1)
    dNdxi(2,5)= 0.125*rPlus*tMinus
    dNdxi(2,6)= 0.125*rMinus*tMinus
    dNdxi(2,7)=-dNdxi(2,6)
    dNdxi(2,8)=-dNdxi(2,5)

    dNdxi(3,1)= 0.125*rPlus*sPlus
    dNdxi(3,2)= 0.125*rMinus*sPlus
    dNdxi(3,3)= 0.125*rMinus*sMinus
    dNdxi(3,4)= 0.125*rPlus*sMinus
    dNdxi(3,5)=-dNdxi(3,1)
    dNdxi(3,6)=-dNdxi(3,2)
    dNdxi(3,7)=-dNdxi(3,3)
    dNdxi(3,8)=-dNdxi(3,4)

    Jmat = matmul(dNdxi, transpose(coords))
    Jdet = Jmat(1,1)*Jmat(2,2)*Jmat(3,3)+ Jmat(1,2)*Jmat(2,3)*Jmat(3,1)+ Jmat(1,3)*Jmat(2,1)*Jmat(3,2)- Jmat(1,3)*Jmat(2,2)*Jmat(3,1)- Jmat(1,2)*Jmat(2,1)*Jmat(3,3)- Jmat(1,1)*Jmat(2,3)*Jmat(3,2)

    !  compute the inverse of Jmat
    Jinv(1,1)=( Jmat(2,2)*Jmat(3,3) - Jmat(2,3)*Jmat(3,2))/Jdet
    Jinv(2,1)=(-Jmat(2,1)*Jmat(3,3) + Jmat(2,3)*Jmat(3,1))/Jdet
    Jinv(3,1)=( Jmat(2,1)*Jmat(3,2) - Jmat(2,2)*Jmat(3,1))/Jdet
    Jinv(1,2)=(-Jmat(1,2)*Jmat(3,3) + Jmat(1,3)*Jmat(3,2))/Jdet
    Jinv(2,2)=( Jmat(1,1)*Jmat(3,3) - Jmat(1,3)*Jmat(3,1))/Jdet
    Jinv(3,2)=(-Jmat(1,1)*Jmat(3,2) + Jmat(1,2)*Jmat(3,1))/Jdet
    Jinv(1,3)=( Jmat(1,2)*Jmat(2,3) - Jmat(1,3)*Jmat(2,2))/Jdet
    Jinv(2,3)=(-Jmat(1,1)*Jmat(2,3) + Jmat(1,3)*Jmat(2,1))/Jdet
    Jinv(3,3)=( Jmat(1,1)*Jmat(2,2) - Jmat(1,2)*Jmat(2,1))/Jdet

    return
    end subroutine caldNdx

