function get_fibers, file, RIGHT=right, REBIN=rebin, SHOW=show
;+
;
;FUNCTION: GET_FIBERS, FILE, DIR, FILE2, [.../REBIN, /SHOW]
;
;PURPOSE:
;  This program generates smoothed fiber profiles to be used in the
;  creation of the P2VM matrix.
;
;METHOD:
;  0] Dewarped fiber profiles are loaded in and collapsed via a median
;     in the third dimension.
;  1] Each column of the image is median-smoothed with a width of about
;     25 pixels (this is the "fit" image).
;  2] Each column is divided by the total of the values in the column to
;     normalize its value (this is the "norm" image). 
;  3] Each row in the imag is fit with a 5th degree polynomial and
;     replaced with this fit (this is the "smooth" image). 
;
;INPUT:
;  file = The string for the image file, given as: '7160_20111019_Aldebaran_0'
;  /right = Perform fiber smoothing using the right-dewarped fiber image
;  /rebin = Rebin the fiber file, given as a string of the size of each bin. 
;  /show = Display the final smoothed fiber result and save a .ps file.
;
;OUTPUT:
;  finished_image = a stack of the three images (fit, normalized,
;                   smoothed with the polyfit)
;NOTES:
;  -Originally, normalized Gaussian profiles were fit to the fiber
;  images for use in GET_BEST_PARAMS, but when more files were run
;  through later parts of the program it was discovered that these
;  profiles were not very good representations of the raw fiber
;  files. The issue turned out to be due (as was an issue in GET_DEFORM)
;  to GAUSSFIT and its issues with smallscale variation. Instead,
;  columns were smoothed with a median boxcar and normalized to the
;  total. Continuity was imposed across rows using a fifth degree
;  polynomial fit to these smoothed columns.
;
;-
    ;Load in file...
    if keyword_set(right) then part = 'right' else part = 'left'
    dir = 'data/'
    dwfile = strcompress(dir+file+'_'+part+'_dewarped.fits.gz',/remove_all)
   
    ;SAVE MIN, MAX OF IMAGE.  STUPID VALUES PROBLEM?
    tries = 0
    tryagain: image = mrdfits(dwfile, 0, hdr, /fscale, /silent, status=check)
    if check LT 0 then begin
       ;Checking for the file...
       chk = file_search(dwfile)
       if chk EQ "" then begin
          print, "Dewarped file not found, skipping GET_FIBERS"
          goto, theend
       endif
          

       ;Avoiding transient issues with file reading...
       tries++  &  wait, 1
       if tries LT 3 then goto, tryagain else goto, theend
    endif


    stat = fxpar(hdr,'STATUS')
    if strpos(stat,'REJECT') NE -1 then begin
       i = strmid(file,strlen(file)-1,1)
       (scope_varfetch('fib_status_token',level=-1))[i]=0
       print, 'Cannot smooth the fiber profile: ', file
       goto, theend
    endif
    image = median(image, dimension=3) 
    siz = size(image)  &  nl = siz[1]  &  np = siz[2]
    original = image
 



    ;Smooth and impose continuity...
    fit = image  &  fit_smooth = fit

    for i = 0, nl-1 do fit[i,*] = msmooth(image[i,*],np*.05)
    for i = 0, nl-1 do begin
       col = image[i,*]  &  col -= min(col)
       pop = []  &  for n = 0, np-1 do pop = [pop,intarr(col[n]/10.+1)+n]
       coeff = [max(col)-min(col), avg(pop), stddev(pop), min(col)]
       
       mpp_col = mpfitpeak(findgen(np), col, mpp_coeff, nterms=4, estimates=coeff)
       fit[i,*] = mpp_col
     endfor

    fit_norm = fit
    for i = 0, nl-1 do fit_norm[i,*] /= total(fit[i,*])
    for i = 0, np-1 do begin
       dummy = poly_fit(findgen(nl), fit[*,i], 5, yfit = dummmy)
       fit_smooth[*,i] = dummmy
    endfor
    for i = 0, nl-1 do fit_smooth[i,*] /= total(fit_smooth[i,*])



    ;Remove any leftover negatives...
    if total(fit LT 0) NE 0 then remove_negs, fit
    if total(fit_norm LT 0) NE 0 then remove_negs, fit_norm
    if total(fit_smooth LT 0) NE 0 then remove_negs, fit_smooth




    ;Rebin and save...
    if keyword_set(rebin) then begin
       rb = fix(rebin)  &  xdim = nl/rb
       rebinned = fltarr(xdim, np)  &  rebinned2 = rebinned
       for i = 0, xdim-1 do begin
          rebinned[i,*] = total(fit_smooth[rb*i : rb*i+rb-1,*], 1)/rb
          rebinned2[i,*] = total(fit_norm[rb*i : rb*i+rb-1,*], 1)/rb
       endfor
       rb = strcompress(string(rb),/remove_all)

       
       shdr = hdr  &  nhdr = hdr 
       fxaddpar, shdr, 'STATUS', stat, ' BINNED, SMOOTHED FIBER' 
       fxaddpar, nhdr, 'STATUS', stat, ' BINNED, NORMALIZED FIBER' 
       inds = where(shdr NE strjoin(strarr(80)+' '))
       shdr = shdr[inds]  &  nhdr = nhdr[inds]


       writefits, dir+file+'_'+part+'_fitted_smoothed_rebinned'+ $
                  rb+'.fits', rebinned, shdr, /compress
       writefits, dir+file+'_'+part+'_fitted_rebinned'+$
                  rb+'.fits', rebinned2, nhdr, /compress
     endif




    ;Save files...
    shdr = hdr  &  nhdr = hdr  
    fxaddpar, shdr, 'STATUS', stat, ' SMOOTHED FIBER' 
    fxaddpar, nhdr, 'STATUS', stat, ' NORMALIZED FIBER' 
    inds = where(shdr NE strjoin(strarr(80)+' '))
    shdr = shdr[inds]  &  nhdr = nhdr[inds]
    writefits, dir+file+'_'+part+'_fitted_smoothed.fits', $
               fit_smooth, shdr, /compress
    writefits, dir+file+'_'+part+'_fitted.fits', $
               fit_norm, nhdr,/compress



    
    ;Make pretty pictures...
    ;if keyword_set(show) then lim = 1 else lim = 0
    ;for i = 0, lim do begin
    ;   if i EQ 0 then psopen, file+'_'+part+'_fiber_fitting.ps',$
    ;                          /inch, xsize=11, ysize=7, /color
    ;   !p.multi = [0, 4, 1]
    ;   display, image, title = "Cropped", charsize=2 
    ;   display, fit, title = "Fitted", charsize=2,$
    ;            min = min(fit), max=1.05*max(fit)
    ;   display, fit_norm, title = "Normalized Fit", charsize=2, $
    ;            min = min(fit_norm), max=1.05*max(fit_norm)
    ;   display, fit_smooth, title = "Polyfitted Fit", charsize=2,$
    ;            min = min(fit_norm), max=1.05*max(fit_norm)
    ;   if i EQ 0 then psclose 
    ;endfor

    


    ;Returning the final image product
    finished_image = {fs: fit_smooth, fn: fit_norm, sh: shdr, nh: nhdr}
    return, finished_image
    theend: return, 0
end








pro remove_negs, image
;+
;
;FUNCTION: REMOVE_NEGS, IMAGE
;
;PURPOSE:
;  To remove the occasional negative value left over after
;
;METHOD:
;  0] Individual slices of the image are taken and negative numbers
;     identified
;  1] The 9 closest neighbors are scanned for negative numbers, and
;     those that are not negative are medianed to replace the central
;     negative number (using a non-treated image box).
;  2] After all negative numbers are replaced, the treated image slice
;     is added back into the image cube.
;
;INPUT:
;  image = a 3D data cube
;
;OUTPUT:
;  (None) = the image array itself is updated
;
;NOTES: 
;  -Started off with a box of the nearest neighbors, but after testing
;  decided to expand outwards to a 5x5 box. 
;-

  siz = size(image)  &  nl = siz[1]  &  np = siz[2]  &  ni = siz[3]
  if siz[0] EQ 2 then ni = 0
  for a = 0, ni-1 do begin
     temp = image[*,*,a]  &  inds = where(temp LT 0)
     for b = 0, n_elements(inds)-1 do begin
        ind = array_indices(temp, inds[b])  &  indt = ind
        if indt[0] LE 1 then indt[0] = 2  
        if indt[1] LE 1 then indt[1] = 2
        if indt[0] GE nl-2 then indt[0] = nl-3
        if indt[1] GE np-2 then indt[1] = np-3

        box = image[indt[0]-2:indt[0]+2,indt[1]-2:indt[1]+2,a]
        check = where(box GE 0)  
        if total(check) EQ -1 then stop, "CLUMP!" 
        val = median(box[check])
        temp[ind[0],ind[1]] = val
     endfor

     image[*,*,a] = temp
  endfor
end
