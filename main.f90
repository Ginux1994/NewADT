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
    
    
    
    iper = 1
    index = 1
    call ADINI
    
    call assem_element
    
	call exceuteLoad1

    
    
    
    
    
    
    

    
end program
    
   