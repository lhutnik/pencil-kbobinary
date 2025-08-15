pro makecsv

pc_read_pstalk,obj=fs

;xp yp zp vpx vpy vpz rhopswarm aps ux uy uz rho potself                                                                                                                                                                                 
print, tag_names(fs)

tmp = size(fs.t)
nt=tmp[1]
tmp = size(fs.ipar)
npar=tmp[1]

for it=1,nt-1,2 do begin

   fname='file_'+string(it,format="(i04)")+'.txt'
   openw,lun,fname,/get_lun,/append

   printf,lun,'# Timestep=',it," Time=",fs.t[it],format='%s %i %s %f'
   printf,lun,'# particle, x, y, rhopswarm'
   for k=0,npar-1 do begin
      if (fs.aps[k,it] ne 0) then begin
         printf,lun,fs.ipar[k],fs.xp[k,it],fs.yp[k,it],fs.rhopswarm[k,it],format='%i %f %f %f'
      endif
   endfor
   close,lun
   free_lun,lun
endfor

end
