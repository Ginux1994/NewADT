program test 
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
    
    
    
    call read_mainCtrlInf
    call read_timeStepInf
    call read_nodeInf
    call read_ndSHInf
    call read_initInf
    call read_heatFlowCtrlInf
    call read_timeFuncInf
    iper = 1
    call read_heatFlowInf
    call read_element
    index = 1
    call init
    call assem_element
    

    
end program