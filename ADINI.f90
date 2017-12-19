 
    subroutine ADINI
    
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
!    use thrdCondElm
    implicit none
!    use condElmCtrlInf, only: coords
    character(len=20) fileName
    integer:: nJumpLine, i
    fileName = "ADINAT1.DAT"
    open(unit = 5, file = fileName)
    
    if(iper==1) then
		call read_mainCtrlInf
		call read_timeStepInf
		call read_nodeInf
		call read_ndSHInf
		call read_initInf
		call read_heatFlowCtrlInf
		call read_timeFuncInf
		call init
        call read_heatFlowInf
		call read_element

    else 
        rewind 5
        nJumpLine = 11+nNd  ! call read_initInf
        do i=1,nJumpLine
            read(5,*)
        enddo
        
        
		call read_heatFlowCtrlInf
		call read_timeFuncInf
		call init
        call resetTimeIntgConst        
        call read_heatFlowInf



	end if
 
    end subroutine ADINI