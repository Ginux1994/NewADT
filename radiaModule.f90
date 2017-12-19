

!*****************************************************************************    

    module radiaElm
	
    implicit none
! npar information    
    integer:: nRadiaElm, nonLrMat, deathType, nRadiaNd, nRadiaLine, nRadiaSurf, radiaLineType,&
			  radiaEmissType, nRadiaEmissGrp, nRadiaEmiss, unitType
! radiation element information
    integer, allocatable:: radiaNdID(:,:), radiaSFNdID(:,:), shapFLine(:,:), thickRadia(:,:), shapFSF(:)
    real, allocatable:: propEmiss(:,:), shapFNd(:,:), areaRadia(:,:), rlh(:,:), specHeat, deaTime(:), &
	xyz_radia(:,:)
    real:: sigma
	
	
! variable for assem element
     
contains
        subroutine read_radiaElm(elmGrpID)
            use mainCtrlInf, only: nElmGrp, npar, nLatht, nNd
            use nodeInf, only: dofID, x, y, z
            integer:: elmGrpID, matID, surfID, i, j, k, ndID, idLoop=1, ii, jj, kk
            integer:: ndIDTemp(8), dummy2, kg, nod9
			npar(3,elmGrpID) = 1; deathType=npar(4,elmGrpID);
			nRadiaNd = npar(5,elmGrpID); nRadiaLine = npar(6,elmGrpID); nRadiaSurf=npar(7,elmGrpID); 
			radiaLineType=npar(9,elmGrpID);!轴对称=0，平面=1
            radiaEmissType=npar(15,elmGrpID); nRadiaEmissGrp=npar(16,elmGrpID); nRadiaEmiss=npar(17,elmGrpID); 

			allocate(shapFSF(nRadiaSurf))
            allocate(propEmiss(nRadiaEmiss,nRadiaEmissGrp))
            allocate(deaTime(nRadiaSurf))
            allocate(radiaSFNdID(8,nRadiaSurf), xyz_radia(3,nRadiaSurf*8))
    ! read material inf
			read(5,*) unitType, sigma
            do i=1,nRadiaEmissGrp              
				if(radiaEmissType==1) then
					read(5,*) matID, propEmiss(matID, i)
				else 
					read(5,*) matID
					! read(5,*) (propEmiss(matID, i), )
                end if
            enddo   
    ! read element inf
            do i=1,nRadiaSurf ! 对所有传导单元循环，线性及非线性的
                read(5,*) surfID, (ndIDTemp(j), j=1,8), dummy2, kg, shapFSF(surfID), deaTime(surfID)  
                radiaSFNdID(:, surfID) = ndIDTemp
            enddo	
            xyz_radia = 0.0
			do i=1,nRadiaSurf ! 对i号单元所有节点循环
				do j=1,8
					ndID = radiaSFNdID(j,i)
					if(ndID.ne.0) then
                        xyz_radia(1, ndID) = x(ndID); xyz_radia(2, ndID) = Y(ndID); xyz_radia(3, ndID) = Z(ndID);
                    end if
				enddo
            enddo
			
        end subroutine read_radiaElm
		
		subroutine assem_radiaElm
        use solCtrlInf
        use mainCtrlInf, only: analyType, CmatType, timeIntType
        use nodeInf, only:dofID
        use heatFlowCtrlInf, only: nRadiaNd
        use heatFlowInf, only:radiaNdValue, radiaNdInf
            integer:: sfID, elmNdNum, ndID(4), nodeIDtemp1, nodeIDtemp2, i, j, k
            real:: coords(3,4), shapF, emmision
            real:: tempNd(4), tempEnv(4), res(4), N(4), radiaK(4,4)
            if (ind>0) then
                do sfID=1, nRadiaSurf
                    
                    call getCoordsBySFID(sfID, 4, coords)
                    call getNdIDBySFID(sfID, 8, ndID)
                    shapF = shapFSF(sfID); emmision = propEmiss(sfID, 1)
                    
                    do i=1,4
                        nodeIDtemp1 = ndID(i)
                        do j=1, nRadiaNd
                            nodeIDtemp2 = radiaNdInf(1,j) ! 过滤出受辐射节点的输入的环境温度
                            if(nodeIDtemp1==nodeIDtemp2) tempEnv(i) = radiaNdValue(1, j)
                        enddo
                        tempNd(i) = Phi_now(nodeIDtemp1)
                    enddo
                    
                    call calRadiaK(4, emmision, shapF, coords, tempNd, tempEnv, res, radiaK)
                    
                    if(timeIntType==2) then
                        call assemDynVec_c_r(res, ndID)
                    else 
                        call assemDynVec_c_r(res, ndID)
                        if(icount>2)  return
                        if(isref==0) call assemDynMatK_c_r(radiaK, ndID)
                    end if
                    
                enddo
                
            end if

		end subroutine assem_radiaElm

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
                ndID(i) = radiaSFNdID(i,sfID)
            enddo       
        end subroutine getNdIDBySFID
		


		subroutine calRadiaK(ndNum, propEmiss, shapF, coords, tempNd, tempEnv, res, radiaK)
		use Gauss
        use solCtrlInf, only: ind, icount, isref
        use mainCtrlInf, only: timeIntType
        implicit none
        !* param table of the function()
        integer, intent(in)::ndNum
        real, intent(in)::coords(3,4), propEmiss, shapF, tempNd(4), tempEnv(4)
        real, intent(out)::radiaK(4,4), res(4)
    
        !* param table in the function
        ! 1. shap function:
        real:: N(4), dNdxi(2, 4)
        ! 2. gauss integration loop variable
        real::r, s, t, weightTotal, fac
        integer::  xLoop, yLoop, i, j, k, numGaussP = 3, nod5=0
        ! 3. D K matrix, and the temporary variable
        real:: Bmat(2,4), tempNd_=0.0, tempEnv_=0.0, tempDelt=0.0, &
				tempOffset, radiaCoef, radiaCoef_, emmision, areaRadia

        radiaK = 0.0 ! call initRmat(Kmat, 8, 8)
        res = 0.0
    
		tempOffset = 459.69
		if(unitType==1) then
			tempOffset = 273.16
		else 
			tempOffset = 0.0
		end if
        
        do xLoop=1,numGaussP
             r = gaussPoint(xLoop, numGaussP)
            do yLoop = 1,numGaussP
                s = gaussPoint(yLoop, numGaussP)
                    weightTotal = gaussWeight(xLoop, numGaussP)*gaussWeight(yLoop, numGaussP)
                    call dNdxi_radia(r, s, ndNum, nod5, coords, N, areaRadia)
                    tempEnv_ = 0.0; tempNd_ = 0.0; tempDelt = 0.0
                    fac = weightTotal*areaRadia
					! 计算环境平均温度和节点平均温度
					do i=1,ndNum
						tempEnv_ = tempEnv_ + tempEnv(i)*N(i)
						tempNd_ = tempNd_ + tempNd(i)*N(i)
						tempDelt = tempDelt + (-tempNd(i) + tempEnv(i))*N(i)
					enddo
					! 偏移原始值
					tempEnv_ = tempEnv_ + tempOffset; tempNd_ = tempNd_ + tempOffset;
					if(radiaEmissType==1) then
						emmision = propEmiss;
					else
						! call 
					end if		
					! 计算辐射系数
					radiaCoef = emmision*shapF*sigma*(tempEnv_*tempEnv_ + tempNd_*tempNd_)*(tempEnv_ + tempNd_)
					fac = fac*radiaCoef
					! 计算右端项
					do i=1,ndNum
						res(i) = res(i) + fac*tempDelt*N(i)
					enddo
					! 选择性的计算辐射刚度矩阵
					if((icount<=2).and.(isref==0)) then
						if(timeIntType.ne.2) then       
                            do i=1,ndNum
                                do j=1,ndNum
                                    radiaK(i,j) = radiaK(i,j) + N(i)*N(j)*fac
                                enddo
                            enddo
						end if
					end if
            enddo  ! do yLoop = 1,numGaussP
        enddo ! do xLoop=1,numGaussP
        end subroutine calRadiaK
		
		
		subroutine dNdxi_radia(r, s, ndNum, nod5, coords, N, areaRadia)
		
		       ! calculate the shap function and derivatives.
			implicit none
			integer:: i,j,k !loop param
			! declare the type and dimension of the function parameter
			real:: N(4), dNdxi(2,4), areaRadia, Jmat(3,3)
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
			areaRadia = sqrt(temp(1)*temp(3)-temp(2)*temp(2))
			return
		end subroutine dNdxi_radia
    
    end module radiaElm
