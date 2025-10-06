function UNIQ, ARRAY, IDX

; Check the arguments.                                                                                                                                                                                                                                                  
  s = size(ARRAY)
  if (s[0] eq 0) then return, 0         ;A scalar                                                                                                                                                                                                                       
  if n_params() ge 2 then begin         ;IDX supplied?                                                                                                                                                                                                                  
     q = array[idx]
     indices = where(q ne shift(q,-1), count)
     if (count GT 0) then return, idx[indices] $
     else return, n_elements(q)-1
  endif else begin
     indices = where(array ne shift(array, -1), count)
     if (count GT 0) then return, indices $
     else return, n_elements(ARRAY)-1
  endelse
end
