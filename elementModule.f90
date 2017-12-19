   
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
! ��elcal��elment���ϵ�һ��elment���棬������Ϣ����ͼ���ֿ�Ϊ�����Ӻ���
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
! ind=1, ��װ��Ч���иնȾ���
if(ind==1) then   
    
            if(haveLnConvOrFixTempLoad==1) then ! ��������Զ�����̶��¶ȣ�����assem linear matrix
        
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
            
ind=2 ! ��ʼ��װ ���ݾ���    

            if((analyType>0).and.(nLatht==0)) then ! ���analyType=1��˲̬���������������������� ��װ���ݾ���
 
                do elmGrpLoop=1,nElmGrp
                
                    if (npar(1, elmGrpLoop)>3) cycle ! ���Զ������䵥Ԫ��           
                    call exceuteAssem(npar(1, elmGrpLoop)) ! ������װ�������ݻ��Ƿֲ����ݣ��ڵ�Ԫ�ڲ�����
                        ! �������ݼ���󣬲�����ӵ��ն����ϣ�������Ϊ ���غɼ���
                end do
            else if((analyType>0).and.(nLatht>0)) then
                ! ��װ���ݾ���+Ǳ�Ȳ����ĵ�Ч��������
            end if
            
else if(ind==3) then ! ����Ƶ�����
            if(klin>0) then
                ! ��װǱ������
                    
            end if  
!******* ind = 4 ʱ��ʼ���� �����Դ���������װ������Ч����������װ�� ʱ����ֹ����еĸնȾ�����װ **********************************************************************                
else if(ind==4) then  
    
            if(klin==0) return
            do elmGrpLoop=1,nElmGrp
                ! Ŀǰ�ȿ��ǵ�Ԫ�����ԣ������ǲ��Ϸ����ԣ�����npar��ƬӦ�ð����������ͷֿ�
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
     
