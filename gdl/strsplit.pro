function strsplit, string, delimiter, extract=extract
    compile_opt strictarr

    ; If no delimiter is given, default to whitespace
    if n_elements(delimiter) eq 0 then delimiter = ' '

    ; Use stregex to split based on delimiter
    pattern = '[^' + delimiter + ']+'
    tokens = stregex(string, pattern, /extract)

    return, tokens
end
