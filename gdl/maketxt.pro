pro maketxt

pc_read_pstalk,obj=fs
; Tag names are written as the following:
; xp yp zp vpx vpy vpz rhopswarm aps ux uy uz rho potself                                                                                                                                                                                 
print, tag_names(fs) ; Print all tags of read in object 

tmp = size(fs.t)    
nt=tmp[1]            ; Number of times stalked
print, tmp[0]

tmp = size(fs.ipar) 
npar=tmp[1]          ; Number of particles stalked
print, tmp[0]

print, nt
print, npar

print, "The middle ipar column is:"
print, fs.ipar[*,0,0]

print, "The sizes of ipar, xp, and aps arrays are the following:"
print, size(fs.ipar)
print, size(fs.xp)
print, size(fs.aps)

;help, fs.xp
print, fs.xp[*,0,0] ; xp data for first particle?
;print, fs.xp[0,0,*]

for it=0,nt-1 do begin ; For all times stalked
   fname='file_'+string(it,format="(i04)")+'.txt' ; Save printf to this file
   filepath='~/pencil-kbobinary/gravitational-scattering/output/'+fname
   openw,lun,filepath,/get_lun,/append

   printf,lun,'# Timestep=',it," Time=",fs.t[it],format='%s %i %s %f'
   printf,lun,'# particle, x, y, rhopswarm, aps'
   for k=0,npar-1 do begin ; For all particles  
         ;if (fs.aps[k,it] ne 0) then begin
      printf,lun,fs.ipar[k],fs.xp[it,k],fs.yp[it,k],fs.rhopswarm[it,k],fs.aps[it,k],format='%i %f %f %f %f'
         ;endif
   endfor
   close,lun
   free_lun,lun
endfor

end
