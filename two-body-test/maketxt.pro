pro maketxt

pc_read_pstalk,obj=fs

;xp yp zp vpx vpy vpz rhopswarm aps ux uy uz rho potself                                                                                                                                                                                 
print, tag_names(fs) ; Print all tags of read in object 

tmp = size(fs.t)    
nt=tmp[1]            ; Number of times stalked
tmp = size(fs.ipar) 
npar=tmp[1]          ; Number of particles stalked

for it=1,nt-1 do begin ; For all timesteps
   fname='file_'+string(it,format="(i04)")+'.txt' ; Save printf to this file
   filepath='~/pencil-kbobinary/gravitational-scattering/output/'+fname
   openw,lun,filepath,/get_lun,/append

   printf,lun,'# Timestep=',it," Time=",fs.t[it],format='%s %i %s %f'
   printf,lun,'# particle, x, y, rhopswarm, aps'
   for k=0,npar-1 do begin ; For all particles  
      ;if (fs.aps[k,it] ne 0) then begin
      printf,lun,fs.ipar[k],fs.xp[k,it],fs.yp[k,it],fs.rhopswarm[k,it],fs.aps[k,it],format='%i %f %f %f %f'
      ;endif
   endfor
   close,lun
   free_lun,lun
endfor

end
