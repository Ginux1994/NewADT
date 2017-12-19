   
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
    

    

!*****************************************************************************    


	
!*****************************************************************************     
! 把elcal和elment整合到一个elment里面，并把信息读入和计算分开为两个子函数
    module element
    use thrdCondElm
    use radiaElm
    use convElm
    implicit none
contains
        subroutine read_element
            use mainCtrlInf, only: nElmGrp, npar
            integer:: elmGrpID, nparTemp(20)
            do elmGrpID=1,nElmGrp
                nparTemp = 0
                read(5,*) nparTemp
                select case(nparTemp(1))
                case(3)
                    npar(:, elmGrpID) = nparTemp
                    call read_thrdCondElm(elmGrpID)
                case(4)
                    npar(:, elmGrpID) = nparTemp
                    call read_convElm(elmGrpID)
				case(5)
					npar(:, elmGrpID) = nparTemp
					call read_radiaElm(elmGrpID)
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
                    call assem_convElm
            case(5)
                    call assem_radiaElm
            end select
        end subroutine exceuteAssem
        
           
        subroutine execute_HeatGen
		use thrdCondElm
        use mainCtrlInf, only: npar
        use heatFlowInf, only: heatGenInf, nHeatGenElmGrp
			integer:: i, j, k, elmGrpLoop, grpID
			
			do elmGrpLoop=1,nHeatGenElmGrp
                grpID = heatGenInf(2,elmGrpLoop)
				select case(npar(1, grpID))
				case(1)
                    ! call assem_onedCondElm
				case(2)
                    ! call assem_twodCondElm
				case(3)
                    call execute_thrdHtGen(grpID)
				end select
            enddo	
		end subroutine execute_HeatGen
        
        
        subroutine assem_element
            use mainCtrlInf, only: nNd, nSHnd, nLatht, numEquation, npar, nLrElmGrp, nNonLrElmGrp, nElmGrp, &
                CmatType, analyType
            use solCtrlInf, only: ind, DOFMap, numMasterDOF, &
            numSlaveDOF, staMatK_k, numLinearEqua, numNonLnEqua, numTotalEqua, &
            klin, rnorm
            use heatFlowCtrlInf, only: haveLnConvOrFixTempLoad, nFixTempNd
            use heatFlowInf, only:fixTempNdInf
        
            implicit none
            integer:: i, j, k, elmGrpLoop, fileID
! ind=1, 组装等效现行刚度矩阵
if(ind==1) then   
    
            if(haveLnConvOrFixTempLoad==1) then ! 如果有线性对流或固定温度，继续assem linear matrix
        
                do elmGrpLoop=1,nElmGrp

                    if (npar(1, elmGrpLoop)==5) cycle ! radiation element, ignore             
                    call exceuteAssem(npar(1, elmGrpLoop))

                end do
                ! deal withh the fix temperature node
                do i=1,nFixTempNd
                    j = fixTempNdInf(1, i)
                    ! staMatK_k(j,j) = staMatK_k(j,j) * 1.0e10
                end do
            end if
            
ind=2 ! 开始组装 热容矩阵    

            if((analyType>0).and.(nLatht==0)) then ! 如果analyType=1，瞬态分析，或者有相变分析，则 总装热容矩阵
 
                do elmGrpLoop=1,nElmGrp
                
                    if (npar(1, elmGrpLoop)>3) cycle ! 忽略对流辐射单元组           
                    call exceuteAssem(npar(1, elmGrpLoop)) ! 具体组装集中热容还是分布热容，在单元内部消化
                        ! 集中热容计算后，不立马加到刚度阵上，而是作为 外载荷加入
                end do
            else if((analyType>0).and.(nLatht>0)) then
                ! 总装热容矩阵+潜热产生的等效热容向量
            end if
            
else if(ind==3) then ! 进行频域分析
            if(klin>0) then
                ! 总装潜热向量
                    
            end if  
!******* ind = 4 时开始进行 非线性传到矩阵组装，及等效热流向量组装； 时间积分过程中的刚度矩阵组装 **********************************************************************                
else if(ind==4) then  
    
            if(klin==0) return
            do elmGrpLoop=1,nElmGrp
                ! 目前先考虑单元非线性，不考虑材料非线性，后续npar卡片应该按照线性类型分开
                if (npar(1, elmGrpLoop)==5) call exceuteAssem(npar(1, elmGrpLoop))
            enddo

            
end if            
            end subroutine assem_element
    end module element
!*****************************************************************************       
    
    
    
    module onedCondElm   
    end module onedCondElm
            
    module twodCondElm    
    end module twodCondElm


    
!*****************************************************************************    
     
