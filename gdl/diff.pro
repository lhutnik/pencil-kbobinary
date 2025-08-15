function diff, arr
    compile_opt strictarr

    n = n_elements(arr)             ; Return number of element in array
    
    ; Special case for limited arrays
    if n le 1 then begin        ; If elements <=1,
        empty_arr = arr[0:0]*0  ; Provide a zero-length array
        return, empty_arr[0:-1]    
    endif
     
    ; Otherwise, return differences for 2+ elements
    return, arr[1:*] - arr[0:(n-2)]
end
