function uniq, arr
    compile_opt strictarr

    if n_elements(arr) eq 0 then return, -1  ; Check if array is empty

    if n_elements(arr) eq 1 then return, [0] ; Scalar output

    diffs = [1, diff(arr)]                   ; Find differences between consecutive elements

    idx = where(diffs ne 0, count)           ; Find positions with differences

    if count gt 0 then idx = idx - 1

    idx[0] = 0                               ; Include 0th element 

    return, idx
end
