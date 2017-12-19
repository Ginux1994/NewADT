
    module thrdCondElm
	
    implicit none
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
            integer:: elmGrpID, dummy2, kg, nod9
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
            do i=1,nthrdCondElm ! 对所有传导单元循环，线性及非线性的
                read(5,*) elmID, nElmNd(elmID), nGeomNd(elmID), dummy2, elmMatID(elmID), dummy2, kg, deaTime(elmID)
                if(nElmNd(elmID)==0) nElmNd(elmID) = nMaxNd
                if(nGeomNd(elmID)==0) nGeomNd(elmID) = nMaxNd
                read(5,*) (elmNdID(j,i), j=1,8)
!                read(5,*) (elmNdID(j,i), j=9,21)
                idLoop = 1
                do k=1,nElmNd(elmID) ! 对i号单元所有节点循环
                    ndID = elmNdID(k,i)
                    xyz(1, ndID) = x(ndID); xyz(2, ndID) = Y(ndID); xyz(3, ndID) = Z(ndID); 
                    lm(idLoop, elmID) = dofID(1, ndID)
                    idLoop = idLoop + 1
                enddo   
            enddo
        
        end subroutine read_thrdCondElm
    
        subroutine assem_thrdCondElm
        use solCtrlInf
        use mainCtrlInf, only: analyType, CmatType, nLatht, nSHnd
        use nodeInf, only:dofID
		use ndSHInf
            integer:: elmID, elmNdNum, ndID(8), ndIDtemp, i, j, k, nod9
            real:: coords(3,8)
            real:: tempNow(8), tempInc(8), res(8), N(8), dNdxi(8), Bmat(3,8), Dmat(3,3), Cmat0(8,8), Cmat1(8), Kmat(8,8)
            if (ind==1) then
                ! assem linear conductivity matrix
                do elmID=1, nthrdCondElm
                    elmNdNum = nElmNd(elmID)
                    nod9 = elmNdNum - 8
                    call getCoordsByElmID(elmID, elmNdNum, coords)
                    call cal3DCondMat(elmID, elmNdNum, nod9, coords, propK, tempNow, Kmat, res)
                    call getNdIDByElmID(elmID, 8, ndID)
                    call assemStaMatK_k(Kmat, ndID)
                enddo
            end if
            if((analyType>0).and.(ind == 2)) then 
                ! assem linear capacity matrix
                do elmID=1, nthrdCondElm
                    elmNdNum = nElmNd(elmID)
                    nod9 = elmNdNum - 8
                    call getCoordsByElmID(elmID, elmNdNum, coords)
                    call cal3DHeatCapMat(elmID, elmNdNum, nod9, coords, propC, tempNow, tempInc, Cmat0, CMat1, res)
                    call getNdIDByElmID(elmID, 8, ndID)
                    if(CmatType==0) then
                        call assemElmC0Mat(Cmat0, ndID) 
                    else if(CmatType==1) then
                        do i=1,8
                        j = dofID(1, ndID(i))
                        C1Mat(j) = C1Mat(j) + CMat1(i)
                        end do 
                    end if
										
                enddo
				if(nSHnd>0) then	
					do i=1,nSHnd
						ndIDtemp = SHNdID(i)						
						if(CmatType==0) then
							C0mMat(ndIDtemp, ndIDtemp) = C0mMat(ndIDtemp, ndIDtemp) + ndSH(i)
						else 
							C1Mat(ndIDtemp) = C1Mat(ndIDtemp) + ndSH(i)
						end if
					enddo
				end if
                if(nLatht>0) then
					! 累加相变产生的等效热容
				else 
					return
				end if
            end if
			! 后续这些在时间积分中进行
            if(ind==4) then
                ! calculate the nonlinear final system conductivity and effictive heat flows
            end if
			
			if(ind==3) then
				! 频域分析
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
   
		
		subroutine execute_thrdHtGen(grpID)
		use heatFlowInf, only: heatGenInf, heatGenElmID, heatGenValue, loadStepLoop
		implicit none
		integer:: grpID, elmID, nElm, ndID(8)
		integer:: i, j, k
		real:: coords(3,8), heatGen, res(8)
        
		nElm = heatGenInf(3, grpID)
		do i=1,nElm
			elmID = heatGenElmID(i, grpID)
			call getCoordsByElmID(elmID, 8, coords)
			heatGen = heatGenValue(loadStepLoop, grpID)
            call cal3DHeatGenVec(elmID, 8, 0, coords, heatGen, res)
            call getNdIDByElmID(elmID, 8, ndID)
			call assemStaVec_k(res, ndID)
		enddo		
		end subroutine execute_thrdHtGen
		
    end module thrdCondElm

        
!*****************************************************************************    
    subroutine cal3DHeatCapMat(elmID, ndNum, nod9, coords, propC, tempNow, tempInc, Cmat0, CMat1, res)
        use solCtrlInf, only: ind, icount, isref
        use mainCtrlInf, only:CmatType, timeIntType
        use Gauss

        implicit none
        !* param table of the function()
        integer::elmID, ndNum, nod9
        real:: coords(3,8), propC(1), tempNow(8), tempInc(8)
        real:: Cmat0(8,8), Cmat1(8), res(8)
    
        !* param table in the function
        ! 1. shap function:
        real:: N(8), dNdxi(3, 8), Jmat(3,3), Jinv(3,3), Jdet
        ! 2. gauss integration loop variable
        real::r, s, t, weightTotal, fac
        integer::  xLoop, yLoop, zLoop, i, j, k, numGaussP = 3
        ! 3. D K matrix, and the temporary variable
        real:: Dmat(3,3), Bmat(3,8), elmTemp = 0.0, averageTempInc = 0.0, specHeat
    
        Cmat0=0.0; Cmat1=0.0; res=0.0;
if(ind<4) then !linear analysis or nonlinear
        specHeat = propC(1)
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
                
                    if(CmatType==0) then
                       ! call initRmat(Cmat0,8,8)
                        do i=1,8
                            do j=1,8
                                Cmat0(i,j) =  Cmat0(i,j) + N(i)*N(j)*fac
                            enddo
                        enddo
                    else if(CmatType==1) then
                        fac = fac/8
                       ! call initRmat(Cmat1, 8 ,8)
                        do i=1,8
                            Cmat1(i) = fac
                        enddo
                    else 
                        pause
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
                    elmTemp=0.0
                    do i=1,8
                        elmTemp = elmTemp + N(i)*tempNow(i)
                    enddo
                    ! get the nonlinear C-specHeat at the temperature of elmTemp
                    fac = weightTotal*Jdet*specHeat
                
                    if(CmatType>=2) then
                        ! calculate the consistent matrix
                        if(icount>=2) then
                            ! calculate the right hand
                            averageTempInc = 0.0
                            do i=1,8
                                averageTempInc = averageTempInc + N(i)*tempInc(i)
                            enddo
                            do i=1,8
                                res(i) = averageTempInc*N(i)*fac
                            enddo
                            continue
                        end if 
                    
                        if(isref == 0) then
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
                        if(timeIntType==2) then
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
                            if(isref==0) then
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
    
end if !(ind<4) then !linear analysis or nonlinear
       
        end subroutine cal3DHeatCapMat
   
    
    
!*****************************************************************************    
    subroutine cal3DCondMat(elmID, ndNum, nod9, coords, propK, tempNow, Kmat, res)
        use Gauss
        use solCtrlInf, only: ind, icount, isref
        use mainCtrlInf, only: timeIntType
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
    

if (ind<4) then !calculate linear element's conductivity matrix
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
                    if((icount <= 2).and.(isref==0)) then
                        Kmat = Kmat + Ktemp ! rebuild k matrix
                    end if
                    
                enddo   
            enddo
        enddo
        
end if
        end subroutine cal3DCondMat


!*****************************************************************************
    subroutine cal3DHeatGenVec(elmID, ndNum, nod9, coords, heatGen, res)
        use Gauss !, only:gaussPoint,gaussWeight
        implicit none

        real:: coords(3,8), heatGen, res(8)
        integer:: elmID, nod9, ndNum

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
                        res(k) = N(k)*fac
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
        if(Jdet<0.0) Jdet = -Jdet
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
