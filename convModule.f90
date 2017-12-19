

!*****************************************************************************    

    module convElm
	
    implicit none
! npar information    
    integer:: nConvElm, nonLrMat, deathType, nConvNd, nConvLine, nConvSurf, convLineType,&
			  convCoefType, nConvCoefGrp, nConvCoef, curveType
! convtion element information
    integer, allocatable:: convNdID(:,:), convSFNdID(:,:), surfHID(:)
    real, allocatable:: prop_LrH(:), prop_NonLrH(:,:), propNonLrH_now(:), areaConv(:,:), rlh(:,:), specHeat, deaTime(:), &
	xyz_conv(:,:)
    real:: sigma
	
	real, allocatable::propNonH_now(:)
! variable for assem element
     
contains
        subroutine read_convElm(elmGrpID)
            use mainCtrlInf, only: nElmGrp, npar, nLatht, nNd
            use nodeInf, only: dofID, x, y, z
            integer:: elmGrpID, matGrpID, surfID, i, j, k, ndID, idLoop=1, ii, jj, kk
            integer:: ndIDTemp(8), dummy2, kg, nod9
			npar(3,elmGrpID) = 1; deathType=npar(4,elmGrpID);
			nConvNd = npar(5,elmGrpID); nConvLine = npar(6,elmGrpID); nConvSurf=npar(7,elmGrpID); 
			convLineType=npar(9,elmGrpID);!轴对称=0，平面=1
            convCoefType=npar(15,elmGrpID); nConvCoefGrp=npar(16,elmGrpID); nConvCoef=npar(17,elmGrpID); 

            allocate(prop_NonLrH(2*nConvCoef,nConvCoefGrp), prop_LrH(nConvCoefGrp))
            allocate(deaTime(nConvSurf))
            allocate(convSFNdID(8,nConvSurf), xyz_conv(3,nConvSurf*8))
            allocate(surfHID(nConvSurf))
    ! read material inf
            if(convCoefType==1) then ! 常 对流系数，不随温度变化
                do i=1,nConvCoefGrp
                    read(5,*) matGrpID, prop_LrH(matGrpID)
                enddo
            else  ! 变 对流系数，随温度变化
                allocate(propNonLrH_now(nConvCoef))
                do i=1,nConvCoefGrp              
                    read(5,*) matGrpID, curveType
                    read(5,*) (prop_NonLrH(j, matGrpID), j=1, 2*nConvCoef)
                enddo   
            end if
    ! read element inf
            do i=1,nConvSurf ! 对所有传导单元循环，线性及非线性的
                read(5,*) surfID, (ndIDTemp(j), j=1,8), matGrpID, kg, deaTime(surfID)  
                convSFNdID(:, surfID) = ndIDTemp
                surfHID(i) = matGrpID
            enddo	
            xyz_conv = 0.0
			do i=1,nConvSurf ! 对i号单元所有节点循环
				do j=1,8
					ndID = convSFNdID(j,i)
					if(ndID.ne.0) then
                        xyz_conv(1, ndID) = x(ndID); xyz_conv(2, ndID) = Y(ndID); xyz_conv(3, ndID) = Z(ndID);
                    end if
				enddo
            enddo
			
        end subroutine read_convElm
		
		subroutine assem_convElm
        use solCtrlInf
        use mainCtrlInf, only: analyType, CmatType, timeIntType
        use nodeInf, only:dofID
        use heatFlowCtrlInf, only: nConvNd
        use heatFlowInf, only:convNdValue, convNdInf
        implicit none
            integer:: sfID, elmNdNum, ndID(4), nodeIDtemp1, nodeIDtemp2, i, j, k, &
                        matGrpID
            real:: coords(3,4), h
            real:: tempNd(4), tempEnv(4), res(4), N(4), convK(4,4)

if (ind==1) then
                ! assem linear convection matrix
                do sfID=1, nConvSurf                    
                    matGrpID = surfHID(sfID)
                    h = prop_LrH(matGrpID)
                    call getCoordsBySFID(sfID, 4, coords)
                    call getNdIDBySFID(sfID, 8, ndID)                   
                    call calLrConvK(4, coords, tempNd, tempEnv, h, res, convK)            
                    call assemStaMatK_c(convK, ndID)
                    
                enddo
else if(ind==2) then
                ! assem linear convection 载荷
                do sfID=1, nConvSurf
                    matGrpID = surfHID(sfID)
                    h = prop_LrH(matGrpID)
                    call getCoordsBySFID(sfID, 4, coords)
                    call getNdIDBySFID(sfID, 8, ndID)         
                    
                    do i=1,4
                        nodeIDtemp1 = ndID(i)
                        do j=1, nConvNd
                            nodeIDtemp2 = convNdInf(1,j) ! 过滤出受对流节点的输入的环境温度
                            if(nodeIDtemp1==nodeIDtemp2) tempEnv(i) = convNdValue(1, j)
                        enddo
                        tempNd(i) = Phi_now(nodeIDtemp1)
                    enddo
                    
                    call calLrConvK(4, coords, tempNd, tempEnv, h, res, convK)
                    call assemDynVec_c_r(res, ndID)        ! 累加对流右端项到 动态vec中，并在每次时间迭代前归零           
                enddo   
                
else if(ind==4) then

            allocate(propNonH_now(2*nConvCoef))    
            ! assem nonlinear convection matrix
                do sfID=1, nConvSurf
                    
                    matGrpID = surfHID(sfID)
                    propNonH_now = prop_NonLrH(:,matGrpID)
                    call getCoordsBySFID(sfID, 4, coords)
                    call getNdIDBySFID(sfID, 4, ndID)
                    
                    do i=1,4
                        nodeIDtemp1 = ndID(i)
                        do j=1, nConvNd
                            nodeIDtemp2 = convNdInf(1,j) ! 过滤出受对流节点的输入的环境温度
                            if(nodeIDtemp1==nodeIDtemp2) tempEnv(i) = convNdValue(1, j)
                        enddo
                        tempNd(i) = Phi_now(nodeIDtemp1)
                    enddo
                                       
                    call calNonLrConvK(4, coords, tempNd, tempEnv, res, convK)
                    
                    if(timeIntType==2) then
                        call assemDynVec_c_r(res, ndID)
                    else 
                        call assemDynVec_c_r(res, ndID)
                        if(icount>2)  return
                        if(isref==0) call assemDynMatK_c_r(convK, ndID)
                    end if
                    
                enddo
end if

		end subroutine assem_convElm

        subroutine getCoordsBySFID(sfID, nSfNd, coords)
            use nodeInf, only: x, y, z
            integer, intent(in):: nSfNd, sfID
            real, intent(out):: coords(3,4)
            integer:: ndID(4), i
            call getNdIDBySFID(sfID, nSfNd, ndID)
            do i=1,nSfNd
                coords(1,i) = x(ndID(i))
                coords(2,i) = y(ndID(i))
                coords(3,i) = z(ndID(i))
            enddo    
        end subroutine getCoordsBySFID
    
        subroutine getNdIDBySFID(sfID, nSfNd, ndID)
            integer, intent(out):: ndID(4)
            integer, intent(in):: sfID, nSfNd
            integer:: i
            do i=1,4
                ndID(i) = convSFNdID(i,sfID)
            enddo       
        end subroutine getNdIDBySFID
		


		subroutine calNonLrConvK(ndNum, coords, tempNd, tempEnv, res, convK)
		use Gauss
        use solCtrlInf, only: ind, icount, isref
        use mainCtrlInf, only: timeIntType
        implicit none
        !* param table of the function()
        integer, intent(in)::ndNum
        real, intent(in)::coords(3,4), tempNd(4), tempEnv(4)
        real, intent(out)::convK(4,4), res(4)
        real:: dt, dc, propH
        !* param table in the function
        ! 1. shap function:
        real:: N(4), dNdxi(2, 4)
        ! 2. gauss integration loop variable
        real::r, s, t, weightTotal, fac
        integer::  xLoop, yLoop, i, j, k, numGaussP = 3, nod5=0
        ! 3. D K matrix, and the temporary variable
        real:: Bmat(2,4), tempNd_=0.0, tempEnv_=0.0, tempDelt=0.0, &
				tempOffset, convCoef, areaConv

        convK = 0.0 ! call initRmat(Kmat, 8, 8)
        res = 0.0
    
        
        do xLoop=1,numGaussP
             r = gaussPoint(xLoop, numGaussP)
            do yLoop = 1,numGaussP
                s = gaussPoint(yLoop, numGaussP)
                    weightTotal = gaussWeight(xLoop, numGaussP)*gaussWeight(yLoop, numGaussP)
                    call dNdxi_conv(r, s, ndNum, nod5, coords, N, areaConv)
                    fac = weightTotal*areaConv
					! 计算环境平均温度和节点平均温度
					do i=1,ndNum
						tempEnv_ = tempEnv_ + tempEnv(i)*N(i)
						tempNd_ = tempNd_ + tempNd(i)*N(i)
						tempDelt = tempDelt + (tempNd(i) - tempEnv(i))*N(i)
                    enddo
                    
                    do j=2,nConvCoef
                        if(tempDelt<=propNonH_now(j)) then
                            dc = propNonH_now(nConvCoef + j) - propNonH_now(nConvCoef + j - 1)
                            dt = propNonH_now(j) - propNonH_now(j - 1)
                            propH = propNonH_now(nConvCoef + j - 1) + dc*(tempDelt - propNonH_now(j - 1))/dt
                            exit
                        end if
                    enddo
					fac = fac*propH
					! 计算右端项
					do i=1,ndNum
						res(i) = res(i) + fac*N(i)*tempDelt
					enddo
					! 选择性的计算辐射刚度矩阵
					if((icount<=2).and.(isref==0)) then
						if(timeIntType.ne.2) then       
                            do i=1,ndNum
                                do j=1,ndNum
                                    convK(i,j) = convK(i,j) + N(i)*N(j)*fac
                                enddo
                            enddo
						end if
					end if
            enddo  ! do yLoop = 1,numGaussP
        enddo ! do xLoop=1,numGaussP
        end subroutine calNonLrConvK


		subroutine calLrConvK(ndNum, coords,  tempNd, tempEnv, propH, res, convK)
		use Gauss
        use solCtrlInf, only: ind, icount, isref
        use mainCtrlInf, only: timeIntType
        implicit none
        !* param table of the function()
        integer, intent(in)::ndNum
        real, intent(in)::coords(3,4), tempNd(4), tempEnv(4)
        real, intent(out)::convK(4,4), res(4)
        real:: propH
        !* param table in the function
        ! 1. shap function:
        real:: N(4), dNdxi(2, 4)
        ! 2. gauss integration loop variable
        real::r, s, t, weightTotal, fac
        integer::  xLoop, yLoop, i, j, k, numGaussP = 3, nod5=0
        ! 3. D K matrix, and the temporary variable
        real:: Bmat(2,4), tempNd_=0.0, tempEnv_=0.0, tempDelt=0.0, &
				tempOffset, convCoef, areaConv

        convK = 0.0 ! call initRmat(Kmat, 8, 8)
        res = 0.0
    
        
        do xLoop=1,numGaussP
             r = gaussPoint(xLoop, numGaussP)
            do yLoop = 1,numGaussP
                s = gaussPoint(yLoop, numGaussP)
                    weightTotal = gaussWeight(xLoop, numGaussP)*gaussWeight(yLoop, numGaussP)
                    call dNdxi_conv(r, s, ndNum, nod5, coords, N, areaConv)
                    fac = weightTotal*areaConv*propH
					! 计算环境平均温度和节点平均温度
					do i=1,ndNum
						tempEnv_ = tempEnv_ + tempEnv(i)*N(i)
						tempNd_ = tempNd_ + tempNd(i)*N(i)
						tempDelt = tempDelt + (tempNd(i) - tempEnv(i))*N(i)
                    enddo

					! 计算右端项
					do i=1,ndNum
						res(i) = res(i) + fac*N(i)*tempDelt
					enddo
					! 选择性的计算对流刚度矩阵      
                    do i=1,ndNum
                        do j=1,ndNum
                                    convK(i,j) = convK(i,j) + N(i)*N(j)*fac
                        enddo
                    enddo

            enddo  ! do yLoop = 1,numGaussP
        enddo ! do xLoop=1,numGaussP
        end subroutine calLrConvK
		
	
		subroutine dNdxi_conv(r, s, ndNum, nod5, coords, N, areaConv)
		
		       ! calculate the shap function and derivatives.
			implicit none
			integer:: i,j,k !loop param
			! declare the type and dimension of the function parameter
			real:: N(4), dNdxi(2,4), areaConv, Jmat(3,3)
			real::r, s, t, coords(3,4)
			integer:: nod5, elmID, ndNum
			! declare the temporal variable of this fuction
			real::rPlus, sPlus, rMinus, sMinus, rSqure, sSqure
			real:: temp(3)
			rPlus=1.0 + r
			sPlus=1.0 + s
			rMinus=1.0 - r
			sMinus=1.0 - s
			rSqure=1.0 - r*r
			sSqure=1.0 - s*s

			N(1)=0.25*rPlus*sPlus
			N(2)=0.25*rMinus*sPlus
			N(3)=0.25*rMinus*sMinus
			N(4)=0.25*rPlus*sMinus

			dNdxi(1,1)= 0.25*sPlus
			dNdxi(1,2)=-dNdxi(1,1)
			dNdxi(1,3)=-0.25*sMinus
			dNdxi(1,4)=-dNdxi(1,3)

			dNdxi(2,1)= 0.25*rPlus
			dNdxi(2,2)= 0.25*rMinus
			dNdxi(2,3)=-dNdxi(2,2)
			dNdxi(2,4)=-dNdxi(2,1)
			Jmat = 0.0
			do i=1,2
				do j=1,3
					do k=1,4
						Jmat(i,j) = Jmat(i,j) + dNdxi(i,k)*coords(j,k)
					enddo
				enddo
			enddo
			
			temp = 0.0
			do i=1,3
				temp(1) = temp(1) + Jmat(1,i)*Jmat(1,i)
				temp(2) = temp(2) + Jmat(1,i)*Jmat(2,i)
				temp(3) = temp(3) + Jmat(2,i)*Jmat(2,i)
			enddo
			areaConv = sqrt(temp(1)*temp(3)-temp(2)*temp(2))
			return
		end subroutine dNdxi_conv
    
    end module convElm
