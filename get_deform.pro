pro get_deform, file, dir, RIGHT=right, FIBER=fiber, SHOW=show
;+
;
;FUNCTION: GET_DEFORM, FILE, DIR, [...\RIGHT, \FIBER, \SHOW]
;
;PURPOSE:
;  This procedure will identify the deformation profile (a.k.a. the warp
;  of the image, an instrumental effect).
;
;METHOD:
;  0] All images in a sequence are loaded into a giantstack which is
;     collapsed using a average through the 3rd dimension. Not all
;     frames are used in the case of a fiber file. Vertical cuts in
;     brightness are established to constrain rows considered in the
;     next step to those bright enough to contain the feature.
;  1] The minimum of the deepest telluric feature in each row is
;     identified and used as line center. The center is saved as part of
;     a vector of positions corresponding to the row. 
;  2] Positions greater than 50 pixels from the median value of this
;     vector are then thrown out.  
;  3] A fifth degree polynomial is fit to this series of positions, 
;     positions greater than 1 pixel from it thrown out, and a final
;     fifth  degree polynomial fit to this dataset. 
;  4] This final fit will be used to dewarp the raw image file in
;     correct_deform and is saved as a .sav file.
;
;INPUT:
;  file = The original raw image file. File should omit the _#.fits.gz. 
;         Example: "4129_20111016_Aldebaran"
;  dir = the directory in which the raw image file can be found
;  /right = keyword to run this protocol on the right half of the image
;  /fiber = keyword to run this protocol for a fiber file
;  /show = keyword to display the fitted profile over the raw image
;
;OUTPUT:
;  (None) = Protocol will save a .sav file into the data directory with the
;           coefficients of the final fitted fifth degree polynomial
;           (fitting) as well as the locations of the vertical value
;           cuts (lims). 
;           Example: "4129_20111016_Aldebaran_0_left_deform.sav"
;
;NOTES:
;  -Add 2.^15 when reading in raw images to make all values positive.
;  -The function gaussfit is incredibly sensitive to small scale
;  variations in the profile of the feature. To make the procedure more
;  robust, the original method of fitting a Gaussian profile to the
;  telluric feature was replaced by:
;    1] obtaining a median smoothed version of the row
;    2] subtracting out the smoothed row to isolate the line
;    3] identifying the position of the deepest feature in this
;    continuum-subtracted row
;  -Prior to the implementation of the spectral power density fit to
;  obtain spatial frequencies, there was an attempt to improve the speed
;  of the pipeline by cutting the image vertically based on a lower
;  bound 2.5-sigma cut in value to remove noise. After the SPD fit came
;  into use, it was more valuable to save these values to retain the
;  speed of the FFT. Any old notation referring to a "cut_check" option
;  refers to these attempts.
;  -These vertical cuts were retained during the process of fitting the
;  line locations so as to avoid fitting the telluric profile in regions
;  where it did not provide a significant contribution.
;  -One of the more difficult features of this part of the pipeline was
;  avoiding falsely identifying the drop off between the left and right
;  halves of the image as part of our telluric feature. To solve this
;  issue, the feature was identified by subtracting out a
;  median-smoothed spectrum and identifying the minimum.
;  -Sean's original code made horizontal cuts in each row from
;  pixel 135 to 160 so as to provide a smaller region for the Gaussian
;  fit. Given gaussfit's issues with oversensitivity to small
;  features, it seemed prudent to rely on the fifth degree polynomial
;  fit to keep things smooth, and to just use the location of the
;  minimum of the feature as the location of the line center. The idea
;  of using a centroid method was explored but ultimately rejected given
;  the typical narrowness of the feature (2-5 pixels, spanning the
;  trough below zero after continuum subtraction).
;  -Changed WRAPPER and GET_DEFORM such that GET_DEFORM is always run,
;  so that all image files have modified headers
;
;-

  
  ;Reading in the images...
  print, "Reading in files..."
  
  if keyword_set(right) then part = 'right' else part = 'left'
  if keyword_set(fiber) then n_img = 450 else n_img = 1000
  if keyword_set(fiber) then index = 8 else index = 9
  nl = 512  &  np = nl  &  bad_img_check = 0  &  count = 0
  
  giantstack=fltarr(nl,np,n_img)
  for i = 0, index do begin
     num = string(fix(i))  
     filename = strcompress(dir+file+"_"+num+".fits.gz", /remove_all)
     again: image = mrdfits(filename, 0, hdr, /fscale, status=chk)+2.^15
     ;THERE SHOULD BE A CHECK FOR NEGATIVE BRIGHTNESS VALUES HERE
     if chk EQ -2 then begin
        chk = 0  &  wait, 2
        goto, again
     endif else if chk LT 0 then stop

     stat = fxpar(hdr, 'STATUS')
     if ~logical_true(stat) then begin
        img_eval, image, hdr, filename    
        if keyword_set(fiber) then imtype = 'FIBER' $
        else imtype = 'SCIENCE IMAGE'
        fxaddpar, hdr, 'IMGTYPE', imtype, ' TYPE OF IMAGE RECORDED',$
                  after='GAIN'
        hdr = hdr[where(hdr NE strjoin(strarr(80)+' '))]
        writefits, repstr(filename,'.gz',''), image-2.^15, hdr, /compress
     endif else count++ 



     val = 'REJECT'
     if keyword_set(right) then stat = strpos(stat,','+val) $
     else stat = strpos(stat,val+',')
     if stat NE -1 then begin
        if keyword_set(fiber) then tokn = 'fib' else tokn = 'stack' 
        tokn += '_status_token'  &  bad_img_check = 1    
        (scope_varfetch(tokn,level=-1))[i]=0
        
        
        print, "Image is reject, skipping get_deform for this image half"
     endif
     

     siz = size(image)  &  set = siz[3]
     if keyword_set(fiber) then begin
        if i EQ 0 and set NE 50 then giantstack=fltarr(siz[1],siz[2],siz[3]*(index+1))
     endif
     giantstack[* , * , set*i : set*(i+1)-1] = image
  endfor
  if bad_img_check then goto, theend
  if count EQ index+1 AND $
     file_search('data/'+file+'_'+part+'_deform.sav') NE "" $
  then goto, theend
  image=giantstack
  

  print, "Done!"
  print




  ;Cropping the image...
  print, "Cropping..."
  if keyword_set(right) and n_elements(image[*,0,0]) GT nl/2.-1 then begin
     image=image[nl/2.:*, *, *]
  endif else image=image[0:nl/2.-1, *, *]
  print, "Done!"
  print
  

  ;THERE SHOULD BE A CHECK FOR CENTERING AND FLUX/CONTRAST HERE
  ;--->IMG_EVAL

  ;Computing the average of all the images..
  print, "Computing average..."
  if keyword_set(fiber) then begin
     frames = get_fiber_contrast(image)
     continuum = avg(image[*,*,frames],2)
  endif else continuum=avg(image, 2)
  print, "Done!"
  print
     



  ;Locating the location of the deep feature in each row...
  print, "Computing line locations..."
  lineloc = fltarr((size(continuum))[2])
  flux = avg(continuum, 0)




  ;Determining vertical value cuts for fitting of line profile...
  m = max(flux,loc)  &  sd = stdev(flux)  &  val = m - 2.5*sd
  spot_b = (where(abs(flux[0:loc]-val) EQ min(abs(flux[0:loc]-val))))[0]
  spot_t = (where(abs(flux[loc:*]-val) EQ min(abs(flux[loc:*]-val))))[-1] + loc
  if spot_b EQ -1 then spot_b = 0  &  if spot_t EQ -1 then spot_t = n_elements(flux)-1 

  for i = spot_b, spot_t do begin
     wave = continuum[*,i]      ;take a slice

     ;Subtract continuum from the data to isolate the line...
     x = findgen(n_elements(wave)) 
     form = msmooth(wave, 50)
     wave= wave-form         
     dummy = min(wave, loc)
     lineloc[i] = loc
  endfor 

  ;Setting all rows outside the vertical value cut to the nearest identified line location...
  lineloc[0:spot_b] = lineloc[spot_b] 
  lineloc[spot_t:*] = lineloc[spot_t] 
  x=findgen(n_elements(lineloc))

  ;Getting rid of huge outliers and only using measured line centers...
  ;(i.e. within vertical value cuts)
  m = median(lineloc[spot_b:spot_t]) 
  std = min([10, stddev(lineloc[spot_b:spot_t])])
  loc = where(abs(lineloc-m) LT std)
  loc = loc[where(loc GE spot_b AND loc LE spot_t)]
  lineloc_full = lineloc  &  lineloc = lineloc[loc] 
  x_full = x  &  x = x[loc]

                                ;THERE SHOULD BE A CHECK ON THE NUMBER
                                ;OF POINTS INVOLVED IN THE FIT HERE
  
  fitting = poly_fit(x, lineloc, 5)
  output=fitting[0]  + $
         fitting[1] * x_full^1 + $
         fitting[2] * x_full^2 + $
         fitting[3] * x_full^3 + $
         fitting[4] * x_full^4 + $
         fitting[5] * x_full^5
  



  ;Select pixels that are less than one pixel away from the fitted line...
  goodpixels=where(abs(output-lineloc_full) lt 1)
  lineloc = lineloc_full  &  lineloc = lineloc[goodpixels] 
  m = median(lineloc) 
  x = x_full  &  x = x[goodpixels]

                                ;THERE SHOULD BE A CHECK ON THE NUMBER
                                ;OF POINTS INVOLVED IN THE FIT HERE

  fitting = poly_fit(x, lineloc, 5)
  output=fitting[0]  + $
         fitting[1] * x_full^1 + $
         fitting[2] * x_full^2 + $
         fitting[3] * x_full^3 + $
         fitting[4] * x_full^4 + $
         fitting[5] * x_full^5
  
  if spot_b NE 0 then begin
     if avg((abs(output-lineloc_full))[0:spot_b]) GT stddev(lineloc) $
     then output[0:spot_b] = output[spot_b+1]
  endif
  
  if spot_t NE np-1 then begin
     if avg((abs(output-lineloc_full))[spot_t:*]) GT stddev(lineloc) $
     then output[spot_t:*] = output[spot_t-1]
  endif
  



  ;Saving a plot of the fitted deformation profile if SHOW keyword is set...
  ;lim = keyword_set(show)
  ;for i = 0, lim do begin
  ;   if i eq 0 then begin
  ;      psopen, strcompress('data/'+file+'_'+part+".ps",/remove_all), /inch,$
  ;              xsize=6, ysize=6, /color, /xstyle, /ystyle 
  ;      loadct, 39
  ;   endif
  ;   plot, x, lineloc, xtitle="Pixel", ytitle ="Location of Line", psym=3, $
  ;         title="Deformation Profile"
  ;   oplot, x, output, color=100
  ;   if i eq 0 then psclose
  ;endfor




  ;Saving the fitted deformation profile...
  print, "Saving deformation profile..."
  lims = [spot_b, spot_t]
  save, filename=strcompress('data/'+file+"_"+part+"_deform.sav", /remove_all), $
        fitting, output, lims, /verbose
  print, "Done!"
  theend: print
end








function get_fiber_contrast, data
;+
;FUNCTION: GET_FIBER_CONTRAST(DATA)
;
;PURPOSE:
;  This function checks attempts to measure the level of blurring that
;  occurs when averaging a 3D image through the third dimension, and
;  identifies the starting frame and number of subsequent frames that
;  produces an average with the best contrast.
; 
;METHOD:
;  0] Identifies brightest row from average frame and creates a
;     "characteristic" row for each frame out of the nearest fifth of
;     the image.
;  1] Sorts the frames on the basis of the depth of the telluric feature
;     being used to dewarp the images.
;  2] Adds frames together one by one in order of decreasing feature
;     depth until the ratio of the average feature depth to the average
;     dispersion of the spectrum in the brightest row of the summed frame
;     begins to decrease.
;
;INPUT:
;  data = a 3D array where the first two dimensions are the horizontal
;         and vertical axes of the CCD and the third dimension represents
;         unique frames.
;
;OUTPUT:
;  frames = a 1D array of positions (corresponding to the 3rd dimension
;           of the image) of frames that should be combined to maximize 
;           contrast 
;
;NOTES: 
;  -Originally this function was called KEEP_CONTRAST and was designed
;  for use in get_spatial_freq, and so a method that looked for a
;  contrast ratio between the closest local minimum and the peak of a
;  column made sense. As this function has been repurposed for
;  identifying the best combination of frames to create contrast in raw
;  fiber files (to locate the telluric feature necessary to dewarp
;  them), it has been renamed GET_CONTRAST and much of the method
;  fundamentally altered.
;  -After some experimentation with methods to best generate contrast,
;  it was decided that this function was best designed to be very
;  specialized to the problem. Rather than establishing a general
;  brightness and "average pixel" metric, a criterion based on the depth
;  of the telluric feature was used to establish the frames with the
;  best contrast.
;-

  siz = size(data)  &  ndim = siz[0]
  if ndim LT 3 then stop, "Image must have 3 dimensions."
  nl = siz[1]  &  np = siz[2]  &  nF = siz[3]
  



  ;Creating the characteristic row for each frame
  avg_column = avg(avg(data,2),0)
  dummy = max(avg_column, bright_row)
  if bright_row-.1*np LT 0 then bright_row = .1*np
  if bright_row+.1*np GE np then bright_row = np-.1*np-1
  row_stack = total(data[*,bright_row-.1*np : bright_row+.1*np,*],2)




  ;score and stack
  row_difs = row_stack
  for a = 0, nF-1 do row_difs[*,a] = row_stack[*,a]-msmooth(row_stack[*,a], .2*nl)
  frame_score = abs(min(row_difs, dimension = 1))
  frame_sort = reverse(sort(frame_score)) ;descending order




  ;while loop to combine frames and assess score
  ib = 50  &  i = ib  &  base = total(data[*,*,frame_sort[0:ib-1]],3)
  last_score = -1  &  new_score = 0
  while new_score GT last_score AND i LT nF-1 do begin
     last_score = new_score
     i++

     stack = base+total(data[*,*,frame_sort[ib:i]],3)
     dummy = max(avg(stack,0), bright_row)
     check_row = stack[*,bright_row]
     check_row -= msmooth(check_row, n_elements(check_row)*.2)
     avg_min = min(check_row)/i  &  avg_dispersion = total(abs(check_row))/i
     new_score = abs(avg_min)/avg_dispersion  
  endwhile


  ;NEED TO HAVE A CHECK TO WARN IF A VERY LOW NUMBER OF FRAMES IS USED


  ;Reporting
  frames = frame_sort[0:i]
  return, frames
end







function msmooth, arr, wid
;+
;FUNCTION: MSMOOTH(ARR, WID)
;
;PURPOSE:
;  To remove small-scale features and noise in a 1D vector using a running
;  median approach. 
;
;INPUT:
;  arr = the 1D array to be smoothed
;  wid = the width of the kernel to smooth the array by
;
;OUTPUT:
;  new_arr = an array of the same size as arr, filled with the running
;  median results.
;-

  new_arr = arr
  for a = wid/2., n_elements(arr)-wid/2.-1 do begin
     new_arr[a] = median(arr[ a-wid/2. : a+wid/2. ])
  endfor

  return, new_arr
end









pro img_eval, img, hdr, ofile
;+
;
;FUNCTION: IMG_EVAL(IMG)
;
;PURPOSE:
;  This function will perform the necessary checks to rate each image
;  for its suitability to run through the pipeline based on its
;  centering, brightness, and other features.
;
;INPUT:
;  img = The data cube to be evaluated.
;  hdr = The hdr corresponding to the data cube being evaluated.
;  ofile = The filename corresponding to the data cube being evaluated.
;
;OUTPUT:
;  lines = The line to be added to the header to mark the image as
;          sound, suspect, or dubious. The line with the vertical center
;          of the fringes
;
;NOTES:
;  ;also need to work on spectra for all of our stars/finding analogous
;stars/creating a reference log for each star-->see laptop
;-


  flat = median(img,dimension=3)  &  siz = size(flat)  
  nl = siz[1]  &  np = siz[2]  &  split = dblarr(nl/2.,np,2) 
  split[*,*,0] = flat[0:nl/2-1,*]  &  split[*,*,1] = flat[nl/2:*,*]
  value = fltarr(2,4,2) ;left/right, 4 criteria, value/score 
  stats = []


  for i = 0, 1 do begin
     ;criterion 1 = vertical alignment
     ;   centered = stdev fits on either side
     ;   suspect = stdev on one side only, enough to fit all baselines
     ;   bad = not all baselines fit (do we need both halves?)

     ;Checking the vertical centering of the fringes in the image
     hmed = median(split[*,*,i],dimension =1)   
     dummy = max(hmed,loc)  &  value[i,0,0] = loc
     case 1 of 
        abs(loc-256) LE 40: value[i,0,1] = 1
        abs(loc-256) GT 40 AND abs(loc-256) LE 150: value[i,0,1] = 0
        abs(loc-256) GT 150: value[i,0,1] = -1
     endcase
                                ;leaving the cuts, but convolving the
                                ;results of crit2 because else too strict
     

     ;1.5 sigma cut to give a bit more space
     ;criterion 2 = vertical spread
     ;  longest possible wavelength is about 250 pixels...
     ;  215 gives wiggle room for .85um
     pop=[]  &  for j = 0, np-1 do pop = [pop,intarr(hmed[j]/1000.)+j]
     spd = stdev(float(pop))  &  inds = loc+[-1,1]*spd*1.5  &  value[i,1,0] = spd
     if inds[0] LT 0 then inds[0] = 0  
     if inds[1] GT np-1 then inds[1] = np-1
     case 1 of 
        abs(total([-1,1]*inds)) GE 250: value[i,1,1]=1 
        abs(total([-1,1]*inds)) LT 215: value[i,1,1]=-1
        else: value[i,1,1]=0
     endcase
     ;hokay hokay, so we do this with the hist-pop assembly and the finding
     ;of the stdev which can be used to define a range around the maxloc
     ;position. Take the minimum of these left and right places (if it
     ;exists), and if they define a range of 250+ (how much leeway? I say
     ;more than right on the money, give it 240...? 1um = 253, .85um = 215)




     ;criterion 3 = feature depth-->dep
     ;   significant = >.15-.2
     ;   suspicious = .05-.15
     ;   trash = LT .05
     spec = split[*,loc,i]    
     sm = msmooth(spec,.1*n_elements(spec)) 
     dep = abs(min(spec-sm,loc)/sm[loc])  &  value[i,2,0] = dep
     case 1 of 
        dep GT .15: value[i,2,1] = 1
        dep GE .05 AND dep LE .15: value[i,2,1] = 0
        dep LT .05: value[i,2,1] = -1
     endcase
     ;depth of deep feature (as a measure of brightness of entire image)
     ;  sadly, this criterion is something best performed after get_deform
     ;...unless we just use the brightest row (because otherwise dual)




     ;criterion 4 = flatness of spectrum-->rat
     ;   all good = LT .5
     ;   suspect = .5-.75/.85
     ;   give up = GT .75
     ;   give up completely = GT .9
     rat = min(spec)/max(spec)  &  value[i,3,0] = rat
     case 1 of 
        rat LT .5: value[i,3,1] = 1
        rat GE .85: value[i,3,1] = -1
        else: value[i,3,1] = 0
     endcase


     ;Rating...
     ;Note: The first criterion can't set the stage for a rejected
     ;      file alone because the spread has to be considered
     rate = ['SOUND','SUSPECT','DUBIOUS','REJECT']
     test = reform(fix(value[i,*,1]))
     case 1 of
        total(test EQ 1) EQ 4: stats = [stats,rate[0]]

        total(test[1:*] EQ -1) NE 0: stats = [stats,rate[3]]
        
        total(test[1:*] EQ 0) EQ 3 and test[0] LE 0: stats= [stats,rate[2]]

        else: stats = [stats,rate[1]]
        ;total(test[0:1] EQ 1) EQ 2: stats=[stats,rate[1]]

        ;total(test[0:1] EQ [0,1]) EQ 2: begin
        ;   if total(test[2:3] EQ 0) EQ 2 then stats = [stats,rate[2]] $
        ;   else stats = [stats,rate[1]]
        ;end
        
        

        ;total(test[0:1] EQ [1,0]) EQ 2: begin
        ;   if total(test[2:3] EQ 1) EQ 2 then stats = [stats,rate[1]] $
        ;   else stats = [stats,rate[2]]
        ;end
        
        ;else: stats = [stats,rate[2]]
     endcase
  endfor




  check = total(where(total(total(img,2),1) EQ 512.^2 * 2.^15) NE -1)
  if check NE 0 then begin
     stats = ['REJECT','REJECT']
     openu, unit, 'zero_files.txt', /get_lun, /append
     printf, unit, ofile
     close, unit
  endif

;and need to note that these ratings will be different for the
;left/right side...so comment and program needs to deal with both
  vert_tag = strcompress(strjoin(string(fix(value[*,0,0])),','),/remove_all)
  fxaddpar, hdr, 'VERTCEN', vert_tag, $
            ' VERTICAL CENTERING OF FRINGES, (L,R)', before='END'
  
  spd_tag = strcompress(strjoin(string(fix(value[*,1,0])),','),/remove_all)
  fxaddpar, hdr, 'VERTSPD', spd_tag, ' VERTICAL SPREAD OF FRINGES, (L,R)'

  dep = string(reform(abs(value[*,2,0])))
  dep = reform((strmid(dep, strpos(dep,'.'), [4,4]))[0,*])
  dep_tag = strcompress(strjoin(dep,','),/remove_all)
  fxaddpar, hdr, 'FEATURE', dep_tag, ' DEPTH OF TELLURIC FEATURE, (L,R)'

  rat = string(reform(abs(value[*,3,0])))
  rat = reform((strmid(rat, strpos(rat,'.'), [4,4]))[0,*])
  rat_tag  = strcompress(strjoin(rat,','),/remove_all)
  fxaddpar, hdr, 'RANGE', rat_tag, ' TYPICAL MIN/MAX RATIO, (L,R)'

  
  stat_tag = strjoin(stats,',')  &  cmnt = ' READY FOR PROCESSING'
  if total(stats EQ 'REJECT') EQ 2 then cmnt = " CAN'T BE PROCESSED"
  fxaddpar, hdr, 'STATUS', stat_tag, cmnt
  hdr = hdr[where(hdr NE strjoin(strarr(80)+' '))]

end
