;This function is designed to allow for the review of the dewarped and
;smoothed fiber files run en masse through the pipeline, as errors
;present in the smoothed stacks are causing severe problems in GBP

;NOTE: This function requires that the FIRST hard drive is connected to
;the computer if /hd is selected



pro fiber_review, hd=hd
;/hd = this keyword will search the hard drive rather than the current
;data amassed on the computer

  if keyword_set(hd) then begin
     dir = '/var/run/media/bbordwell/FIRST/Lick_Oct2011/Dewarped_Data_&_Processed_Fibers/*/'
     tag = 'hd'
  endif else begin
     dir = '~/urap/scripts/idl/current_code/data/'
     tag = 'cc'
  endelse





  openw, 110, tag+'_fiber_review.txt'
  printf, 110, 'This file contains an assessment of processed fiber files run through the FIRST IDL pipeline as of 140212'
  printf,110,""
  printf, 110, "Fibers marked GOOD have no NaNs in the smoothed image, BAD do. In the name of this file, hd implies that files on the hard disk were scanned and cc implies that files in the current_code directory have been scanned."
  printf,110,""
  printf,110,""
  printf,110,""

  for a = 0, 1 do begin
     if a then files = file_search(dir+'*fibre*norm_stack*') $
     else files = file_search(dir+'*fibre*smoothed_stack*')
  
     for i = 0, n_elements(files)-1 do begin
        im = mrdfits(files[i])
        printf, 110, files[i], ':'
        for j = 0, 8 do begin
           if total(where(im[*,*,j] NE im[*,*,j])) NE -1 then rep = "BAD" else rep = "GOOD"
           printf, 110, j, '           ',rep
        endfor
        printf, 110, ""
     endfor
  endfor


  close, 110
end
