pro correct_deform, file, num, dir, RIGHT=right, SHOW=show
;+
;
;FUNCTION: CORRECT_DEFORM, FILE, NUM, DIR, [.../MAKEDARK, /RIGHT, /SHOW]
;
;PURPOSE:
;  To correct the warp introduced into the image by the optics based on
;  the deformation profile measured in get_deform.
;
;METHOD:
;  0] Image is read in, the center band identified using BAND_FINDER and
;     its exposure time used to identify a relevant master dark to be
;     scaled for dark correction.
;  1] Scaled master dark (by the ratio of the medians of the band
;     region) is subtracted from each frame
;  2] Deformation profile is loaded in, and used to pad the array that
;     will be used to shift the rows of the image individually. Spline
;     interpolation is used to shift each row appropriately.
;  3] The leftmost and rightmost columns containing a 0 are identified
;     and used to isolate the image region of the array.
;  4] Check for oversubtraction of the dark, update the header and save
;     the new, dewarped fits file.
;
;INPUT:
;  file = The stub of the data sequence to be dewarped.
;         Example: '4019_20111016_Aldebaran'
;  num = The number of the particular datacube in the sequence to be dewarped.
;  dir = The path to the directory in which the raw data file can be found.
;  /right = A keyword denoting that the right half of the image is being
;           treated, rather than the left as is default.
;  /show = A keyword to plot results of image dewarping to the screen.  
;
;OUTPUT:
;  (None) - writes dewarped image cube to a file. 
;           Example: 7160_20111019_Aldebaran_0_left_dewarped.fits.gz
;
;NOTES:
;  -Removed makedark keyword in favor or having the program create a
;  dark if a relevant dark file was not located upon search. 
;  -Three methods were tried for the identification of bad pixels using
;  the dark frames. The first method (in use) looked for pixels that
;  were above a 3 sigma cut around the median in the majority of dark
;  frames available. The second method (not used) performed this same
;  check in a median-stacked image. The third method (not used) aimed to
;  identify bad pixels based on pixel-to-pixel variation in different
;  dimensions.
;  -Prior to identification of the issues with the dark files,
;  correct_deform was set up to search for any dark files that might
;  possibly be relevant to the image to be dark corrected, in absence of
;  an easy method of referencing the log files for each night. After
;  exploration of the rapid variation of the dark files, it was decided
;  that a better solution involved creating general master darks based
;  on exposure times using darks from the same night, and scaling a
;  normalized version based on the median of the band region of the image.
;  -Once long ago, an issue of having completely miscut images was
;  introduced by cutting the image (horizontally) prior to applying the
;  spline. This was remedied and the cutting moved to the EDGE_CORRECT
;  procedure in order to have uniform cuts amongst a single sequence.
;  -Changed fitting to be centered around 0 because it was silly to use
;  so many columns in dewarped.
;  -Used to use while loops to identify columns with padding columns
;  mixed in...but I finally found a better solution.
;  -Removing the band in this step, rather than edge correction, because
;  is actually an issue post-dewarping.
;
;-
    print, "Dewarping..."

    ;Reading in raw data file...
    file = repstr(file, '.fits.gz', '')
    date = strmid(file,5,8)
    if strpos(dir,'/',/reverse_search) NE strlen(dir)-1 then dir += '/'
    image=mrdfits(strcompress(dir+file+'_'+num+'.fits.gz',/remove_all),$
                  0, hdr, /fscale, /silent) + 2.^15
    siz = size(image)  &  nl = siz[1]  &  np = siz[2]  &  ni = siz[3]
    



    ;Checking the image status and working with /RIGHT...
    stat = fxpar(hdr, 'STATUS')   
    if logical_true(stat) EQ 0 then goto, theend


    chk = fxpar(hdr,'IMGTYPE')  &  stat = strcompress(stat,/remove_all)
    pos = strpos(stat,',')  
    if keyword_set(right) then begin
       part = 'right'  &  poss = [pos+1,strlen(stat)-pos-1]
    endif else begin
       part = 'left'  &  poss = [0,pos]
    endelse
    stat = strmid(stat,poss[0], poss[1])

    if stat EQ 'REJECT' then begin
       if strpos(chk,'FIBER') NE -1 then tokn = 'fib' else tokn = 'stack'
       tokn += '_status_token'    
       (scope_varfetch(tokn,level=-1))[fix(num)]=0


       print, "Image is reject, skipping correct_deform for this image half"
       goto, theend
    endif


    fxaddpar, hdr, 'SIDE', strupcase(part), ' SELECTED HALF OF IMAGE',$
              before = 'STATUS'





    ;Locating and scaling the dark frame...
    dt = string(fxpar(headfits(dir+file+'_0.fits.gz'), 'SHUTTER'))
    darkfile = strcompress('data/master_dark_'+date+'_'+dt+'.fits.gz',/remove_all)
    spectrum = median(median(image, dimension = 3), dimension =2)
    band = band_finder(spectrum)
    dark = mrdfits(darkfile, /silent)
    scale = median(median(image[band[0]:band[1],*,*], dimension =3))/ $
            median(dark[band[0]:band[1],*])
    dark *= scale
    ;REPORT MIN/MAX DARKS IN LOG
    


    ;Cutting the image and dark, then subtracting dark...
    ;if keyword_set(right) then image = image[nl/2.:*,*,*] $
    ;else image = image[0:nl/2.-1,*,*]
    ;if keyword_set(right) then dark = dark[nl/2.:*,*] $
    ;else dark = dark[0:nl/2.-1,*]
    for i = 0, ni-1 do image[*,*,i] -= dark    
    if keyword_set(right) then image = image[band[1]:*,*,*] $
    else image = image[0:band[0],*,*]

    ;Checking the dark correction for oversubtraction...
    check = avg(image LT -1*sqrt(median(dark)))
    check = fix(check*100)/100.
    chk_str = strmid(strcompress(string(check),/remove_all),0,3)
    fxaddpar, hdr, 'DARKCHK', chk_str, $
              ' CHECK FOR NEGATIVES AFTER SUBTRACTION',before='STATUS'
    case 1 of
       check GE .27: stat = 'REJECT'
       check GE .21: stat = 'DUBIOUS'
       check GT .16: if stat NE 'DUBIOUS' then stat = 'SUSPECT'          
       else: break
    endcase
    image[where(image EQ 0)] += .1



    ;Loading in the deformation profile...
    restore, 'data/'+file+'_'+part+'_deform.sav'
    ;x = findgen(np)
    ;output=fitting[0]-floor(fitting[0]) + $
    ;       fitting[1] * x^1 + $
    ;       fitting[2] * x^2 + $
    ;       fitting[3] * x^3 + $
    ;       fitting[4] * x^4 + $
    ;       fitting[5] * x^5
    ;output is now preloaded to deal with bad fit issues


    ;Padding the image with columns of 0s prior to dewarping...
    dewarped = [fltarr(max(abs(output)) + 1, np, ni),$
                image,$
                fltarr(max(abs(output)) + 1, np, ni)]


    ;Shifting row by row...
    for i = 0, n_elements(output)-1 do begin
       arr = findgen(n_elements(dewarped[*,i,*])) 
       dewarped[*,i,*] = interpol(dewarped[*,i,*], arr, $
                                  arr + 1*output[i], /spline) 
    endfor



    
    ;Identifying columns shared by all rows...
    slice = total(dewarped[*,*,0] EQ 0,2)  &  nm = n_elements(slice)
    inds = where(slice EQ 0)
    locs = intarr(n_elements(inds))
    for i = 0, n_elements(locs)-1 do begin
       j = inds[i]  & d = 0 
       while total(slice[j:j+d]) EQ 0 AND j+d LT nm do d++
       locs[i] = d  
    endfor
    rhs = max(locs, lhs)  &  lhs = inds[lhs]  &  rhs += lhs
    dewarped = dewarped[lhs:rhs,*,*]  &  image = dewarped  




    ;Saving...
    scale = strcompress(string(scale),/remove_all) 
    scale = strmid(scale,0,strpos(scale, '.')+4)
    if strlen(scale) GT 9 then stop, "WTF"
    fxaddpar, hdr, 'DSCALE', scale, ' SCALING FACTOR FOR DARK',$
              before = 'STATUS'
    band = strcompress(strjoin(string(fix(band)),':'),/remove_all)
    fxaddpar, hdr, 'DBAND', band, ' INDICES OF CENTER BAND', before='STATUS'
    fxaddpar, hdr, 'STATUS', stat,' DEWARPED', before='END'
    hdr = hdr[where(hdr NE strjoin(strarr(80)+' '))]

    if (size(dewarped))[1] LT 100 then stop,"HA"





    writefits, strcompress('data/'+file+'_'+num+'_'+part+'_dewarped1.fits',$
                           /remove_all), dewarped, hdr, /compress





    ;Display before and after...
    if keyword_set(show) then begin
       !p.multi = [0,2,1]
       image = median(image, dimension=3)
       display, image, min = min(image), max=max(image)
       display, median(dewarped, dimension=3), min = min(image), max=max(image)
       !p.multi = 0
    endif

    theend: print, "Done!"
end







function band_finder, spectrum
;+
;
;FUNCTION: BAND_FINDER, SPECTRUM
;
;PURPOSE:
;  This procedure identifies the edges of the central dark band on the
;  CCD to allow proper scaling of the manufactured dark.
;
;METHOD:
;  0] Obtain a characteristic spectrum from a data cube, smooth it
;  broadly to eliminate possible noise issues, and slice
;  horizontally to select the region where the band should be expected
;  to occur (100-300).
;  1] Find an approximate first derivative and identify points where the
;  slope is less than .1* the standard deviation of the slope to
;  identify band candidates.
;  2] Find the longest (almost) continuous stretch of points that
;  satisfy the last condition (note: the almost is to denote that one
;  point can be an outlier)
;  3] Identify the edges of that population and return the indices. 
;
;INPUT:
;  spectrum = A horizontal profile of the data cube collapsed such that
;             the structure in the bright parts of the image is
;             well-characterized.  
;
;OUTPUT:
;  bounds = a two element array containing the left and right indices of
;           the dark band of the image. 
;
;NOTES: 
;  -It is important to note that the spectrum is best created by median
;  collapsing through the third and second dimensions of the image,
;  otherwise images that are off center may have unexpected flat regions
;  that throw off this method.
;
;-  
  ;Getting the spectrum smooth and relevant...
  sm_spec = msmooth(spectrum, n_elements(spectrum)*.1)
  sm_spec = sm_spec[100:300] 




  ;Using the slope to find flat spots...
  ydif = abs(sm_spec[1:*]-(shift(sm_spec,1))[1:*])
  ydif[where(ydif LT stdev(ydif)*.1)] = 0

  


  ;Finding the biggest flat spot...
  nyd = n_elements(ydif)  &  locs = fltarr(nyd)  &  check = locs
  for i = 0, nyd-1 do begin
     d = 0  &  val = stdev(ydif)*.15
     while total(ydif[i:i+d]) LT val AND d LT nyd-i-1 do d++
     locs[i] = d  &  check[i] = avg(sm_spec[i+1:i+1+d])
  endfor




  ;Getting its indices...
  inds = where(check LT avg(sm_spec))
  locs = locs[inds]
  rhs = max(locs, lhs)  &  lhs = inds[lhs]+101  &  rhs += lhs
  bounds = [lhs,rhs]

  return, bounds
end








pro edge_correct, file, dir, RIGHT=right, FIBER=fiber, RESCUE=rescue
;+
;
;FUNCTION: EDGE_CORRECT, FILE, DIR, [.../RIGHT, /FIBER]
;
;PURPOSE:
;  This procedure loads in all of the dewarped images and ensures a
;  uniform number of pixels across the spectral axis in each data
;  cube. Cuts are made based on the minimum of either the far edge, or a
;  2.5sigma cut from the maximum so as to retain the maximum number of
;  pixels horizontally while removing columns too dark to be useful.
;
;METHOD:
;  0] Creates giant image stack out of all images.
;  1] Removes any columns from the center band using the most
;     conservative band indices from BAND_FINDER.
;  2] Removes any other dark regions left over from dewarping the image
;     by removing all edge values at the minimum.
;  3] Identifies the outside edge value and a 2.5sigma cut from the max
;     using a median spectrum and uses the minimum of the two as the
;     lower limit on allowed values for the image.
;  4] Uses a median-smoothed version of the spectrum to identify where
;     the image sits above this lower limit (avoids issues with
;     features) and identify the index limits of the largest region
;     meeting this criterion.
;  5] Utilizes these limits to cut the stack of images and write each
;     relevant datacube to its final dewarped .fits file.
;
;INPUT:
;  file = the file descriptor up to the cube number, i.e. 1040_20111015_Aldebaran
;  dir = the location of the file to be loaded, i.e. data/
;  /right = denotes that the right side of the images will be treated
;  /fiber = denotes that the files being treated are fiber files
;  /rescue = denotes that some issue occured with the correction of the
;            files, so rather than rerunning the whole program, they
;            will be corrected based off of the good files.
;
;OUTPUT:
;  (None)- Updates the dewarped fits files so that all observations with
;          the same initial four characters share the same dimensions
;
;NOTES:
;  -Originally had a final step in which negative numbers were removed
;  and replaced with the median of their neighbors using REMOVE_NEGS,
;  but found that it was too complicated and too time-consuming to be
;  worth doing for every frame in each data cube. Instead, REMOVE_NEGS
;  will be used on the fibers after the smoothing process in
;  GET_FIBERS. Since a median frame is used, and it is smoothed and
;  continuity enforced across columns, this should evade most of the
;  issues that were occuring with its general implementation.
;-
    ;Reading in the images...
    print, "Reading in files for slicing..."
    if keyword_set(right) then part = 'right' else part = 'left'
    if keyword_set(fiber) then index=8 else index=9
    if keyword_set(fiber) then len=50 else len=100




    ;Allowing easy repair if some tomfoolery happens during correction...
    if keyword_set(rescue) then begin
       fil1 = dir+file+'_*_'+part+'_dewarped1.fits.gz'
       fil = dir+file+'_*_'+part+'_dewarped.fits.gz'
       files1 = file_search(fil1)  &  nf1 = n_elements(files1)
       fil = (file_search(fil))[0]  &  hdr = headfits(fil)
       
       edg_lin = fxpar(hdr,'EDGES')
       edg = fix(strsplit(edg_lin,':',/extract))
       for i = 0, nf1-1 do begin
          img = mrdfits(files1[i],0,hdr)  &  img = img[edg[0]:edg[1],*,*]
          filename = repstr(files1[i],'d1.fits','d.fits')
          writefits, repstr(filename,'.gz',''), img, hdr,/compress
       endfor
       goto, theend
    endif
          
       


    ;Otherwise starting from the beginning...
    nl= 512  &  np = nl
    fil = dir+file+'_0_'+part+'_dewarped1.fits.gz'
    hdr = headfits(fil)
    n = fxpar(hdr,'NAXIS3')  &  if n NE len then len = n
    giantstack=fltarr(nl,np,(index+1)*len)
    bnd_info = strarr(index+1)
    for i = 0, index do begin
       fil = strcompress(dir+file+'_'+string(i)+'_'+part+ $
                         '_dewarped1.fits.gz',/remove_all)
       dummy = mrdfits(fil,0,hdr)  &  sz = size(dummy)
       if i EQ 0 then hdrs = strarr(index+1,n_elements(hdr))


       hdrs[i,*] = hdr
       bnd_info[i] = fxpar(hdr, 'DBAND ')
       giantstack[0:sz[1]-1,*,i*len:(i+1)*len-1] = dummy
    endfor

    ;Sorting out the band cuts...
    pos = strpos(bnd_info, ':') 
    if keyword_set(right) then begin
       rhs = fix(reform((strmid(bnd_info,pos+1,3))[0,*]))
       mx = max(rhs)  &  rh = rhs-mx  &  ind = where(rh NE 0)
       for i = 0, n_elements(ind)-1 do $ 
          giantstack[*,*,ind[i]*len: (ind[i]+1)*len-1] =$
          giantstack[shift(findgen(nl),rh[ind[i]]),*,ind[i]*len: (ind[i]+1)*len-1]
       ;if mx LE nl-1 then mx = 0 else mx += 1 ;keeping numbers straight
       giantstack = giantstack[0:nl-1-mx,*,*]
    endif else begin
       ;can't keep more than was cut from an earlier image
       lhs = min(fix(reform((strmid(bnd_info,0,pos))[0,*])))
       restore, 'data/'+file+'_'+part+'_deform.sav'
       x = findgen(512)
       output=fitting[0]-floor(fitting[0]) + $
              fitting[1] * x^1 + $
              fitting[2] * x^2 + $
              fitting[3] * x^3 + $
              fitting[4] * x^4 + $
              fitting[5] * x^5
       mx = max(abs(output),loc)  &  lhs -= output[loc]
       giantstack = giantstack[0:lhs,*,*] 
    endelse


    ospec = median(median(giantstack, dimension = 3), dimension = 2)    
    spec = msmooth(ospec, n_elements(ospec)*.1)
    mn = min(spec, loc)  &  ind = where(spec GT mn)
    if keyword_set(right) then ind = (where(ind GT loc))[0] $
    else ind = (where(ind LT loc))[-1] 




    ;Keeping only values GT max-2.5sig or the outside boundary...
    inds = where(ospec NE 0)  &  mx = max(ospec[inds])  &  std = stdev(ospec[inds])
    val0 = mx - 2.5*std  &  inds = where(spec NE 0)
    if keyword_set(right) then val = min([val0,spec[inds[-1]]]) $
    else val = min([val0, spec[inds[0]]])
    crit = spec GT val  &  nd = n_elements(spec)  &  score = intarr(nd)
    for i = 0, nd-1 do begin
       d= 0  &  while crit[i+d] NE 0 AND i+d LT nd-1 do d++  &  score[i] = d
    endfor
    rhs = max(score,lhs)  &  rhs += lhs
    giantstack = giantstack[lhs:rhs,*,*]
    band = strcompress(strjoin(string(fix([lhs,rhs])),':'),/remove_all)





    ;Updating the files with the newly cut images and removing negative numbers...
    for i = 0, index do begin
       dewarped = giantstack[*,*,i*len:(i+1)*len-1]
       num = strcompress(string(i),/remove_all)  
       hdr = reform(hdrs[i,*])
       fxaddpar, hdr, 'EDGES', band, ' INDICES OF EDGES OF FRINGES', before='STATUS' 
       fxaddpar, hdr, 'STATUS', fxpar(hdr,'STATUS'), ' DEWARPED+EC', before='END'
       hdr = hdr[where(hdr NE strjoin(strarr(80)+' '))]
       fil = 'data/'+file+'_'+num+'_'+part+'_dewarped.fits.gz'
       spawn, 'rm '+repstr(fil, 'd.fits', 'd1.fits')
       writefits, repstr(fil,'.gz',''), dewarped, hdr, /compress 
    endfor

    theend: print, "Edges corrected..."
end





