 
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
    use thrdCondElm
!    use condElmCtrlInf, only: coords
    character(len=20) fileName
    fileName = "ADINAT.DAT"
    open(unit = 5, file = fileName)
    
    if(iper>1) then
        nste = stepNumAt(iper)
        startTime = previousTime
        call resetTimeIntgConst
    end if
    
    call read_mainCtrlInf
    call read_timeStepInf
    call read_nodeInf
    call read_ndSHInf
    call read_initInf
    call read_heatFlowCtrlInf
    call read_timeFuncInf
    ! nste = stepNumAt(iper)
    
    call read_heatFlowInf
    call read_element
    
    call init
    end subroutine ADINI