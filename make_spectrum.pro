pro make_spectrum, star_name, star_ID, HIPPARCOS=hipparcos, SCALE=scale, $
                   SAVE_PATH=save_path, STAR_TYPE=star_type
;+
;
;FUNCTION: MAKE_SPECTRUM, STAR_NAME, STAR_HD, [...SCALE=scale,
;                         /HIPPARCOS, SAVE_PATH='PATH']
;
;PURPOSE:
;  This procedure will take the star name used in tagging the data files
;  and the corresponding ID number and create a spectrum to be used in
;  DSP_FIT using Pickles spectra and the Hipparcos catalog.
;
;INPUT:
;  star_name = The string used as in the data file naming system to
;              identify the specific target. 
;  star_ID = The number by which the star is identified in stellar
;            catalogues. The default catalog is the Henry Draper
;            catalog. 
;  /scale = The scale factor by which the spectrum should be multiplied
;          to convert it to physical units. Spectra are normalized to
;          their value at 5556A by default. 
;  /hipparcos = If keyword is set then star_ID will be referenced
;               against the Hipparcos catalog.
;  /save_path = A string specifying where the new spectra should be
;               saved. It is expected that you will want to save your
;               spectra in the same location as where the Hipparcos
;               catalog file and the Pickles spectra are saved.
;  /star_type = A string specifying the star type if it is known. In the
;               case that this keyword is specified, star_ID can be set
;               to a null string...or anything really.
;
;OUTPUT:
;  (None) = Program will save the interpolated spectrum as a two column
;           text file with the wavelength values in angstroms (1st
;           column) and the normalized flux values (2nd column).
;
;NOTES:(None)
;-
  
  if ~keyword_set(save_path) then begin
     path = '~/urap/scripts/idl/current_code/reference_spectra/'
  endif else path = save_path
  



  ;Identifying the spectral type of the star...
  if ~keyword_set(star_type) then begin
     hip_cat = mrdfits(path+'Hip_cat.fit',1)
     star_ID = long(star_ID)
     if keyword_set(Hipparcos) then begin
        star_ID = long(repstr(strupcase(star_ID),'HIP',''))
        HIP = long(hip_cat.HIP)
        ind = where(HIP EQ star_ID)
     endif else begin
        star_ID = long(repstr(strupcase(star_ID),'HD',''))
        HD = long(hip_cat.HD)
        ind = where(HD EQ star_ID)
        if ind EQ -1 then stop
     endelse
     
     star_type = (strlowcase((hip_cat.sptype)[ind]))[0] ;scalar
  endif else star_type = strlowcase(star_type)




  ;Obtaining the raw stellar spectrum file...
  again: sp_file = strcompress(path+'pickles/uk'+star_type+'.dat',/remove_all)
  file = file_search(sp_file)
  ans = ""  &  check = 0
  while file EQ "" do begin
     ans = ""  &  check = 1
     print, 'Identified type, '+ star_type+', not present.'
     print, 'Which stellar type should be used?'
     print, '(Answer "closest" if you would like the program to select the '
     print, ' most similar spectral type. Do not give this answer for '
     print, ' spectral types involving complicated characters)'
     print
     read,'Use: ', ans

     if ans EQ 'closest' then break
     star_type = strlowcase(ans)
     sp_file = strcompress(path+'pickles/uk'+star_type+'.dat',/remove_all)
     file = file_search(sp_file)  
  endwhile

  if ans EQ 'closest' then begin
     files = file_search(strcompress(path+'pickles/uk*.dat',/remove_all))
     spawn, 'echo ~', tilde
     sp_file = repstr(sp_file,'~', tilde)
     files = [files,sp_file]  &  files = files[sort(files)]
     ind = where(files EQ sp_file)
     val = strlen(path+'pickles/uk')
     test = [ind-1,ind+1]  &  test = files[test]  &  score = [val,val]
     
     for i = 0,1 do begin
        while strcmp(test[i], sp_file, score[i]) do score[i]++
     endfor

     dummy = max(score,loc)  &  sp_file = test[loc]
  endif

  if check then begin
     print, "I will use: ",sp_file
     ans = ""  &  read, "Is that okay? (Y/N)", ans
     ans = strupcase(ans)  &  if ans EQ 'N' then goto, again
  endif



  ;Reading in and modifying the values in the spectrum file...
  readcol, sp_file, lambda, flux, comment = '#'
  inds = where(lambda GE 5500 AND lambda LE 1e4)
  lambda = lambda[inds]  &  flux = flux[inds]
  if keyword_set(scale) then flux *= scale


  ;Creating the new spectrum file...
  save_name = strcompress(path+'Pickles_'+star_name+'.txt',/remove_all)
  openw, unit, save_name, /get_lun
  for i = 0, n_elements(lambda)-1 do printf, unit, lambda[i], flux[i]
  close, unit

  print, "Spectrum complete!"
end
