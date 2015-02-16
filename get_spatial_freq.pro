function get_spatial_freq, file, right=right, bin = bin, dir = dir
;This is the global function. It adjusts the wavelengths to the column
;location and then fits wavelengths to the fringes. 

;NOTE: ELSA'S UPDATED VERSION OF THIS FUNCTION IS LOCATED AT THE
;END OF THIS CODE

  common constants, p_det, foc, tau, xvg, p_vg, glass, surf, init,part
   p_det = 0.016; pixel size on the detector (mm)
   foc = 150.; focal length (mm)
   tau = .050; integration time (s)
   print, "The integration time is ", strcompress(string(round(tau*1e3)), /remove_all), " milliseconds. If this is not correct, please update the code."
   XVG = [2,3,7,14,27,29,37,43,46] ; fiber positions in the v-groove
   p_VG = 0.037 ; v-groove pitch (mm)
   glass = 2. ; glass is the glass of the prism (BK7 for version 9 or SF2 for version 18. For the 9-fiber instrument, glass=1)
   surf = 9.*!pi*3.^2/(4.*37) ;surface collecting photons
   init = strmid(file, 0, 4)
   if keyword_set(right) then part = 'right' else part = 'left'

   star = strmid(file, 14, strpos(file, '_', /reverse_search)-14) 

   file=repstr(file, '.fits.gz', '')

   if keyword_set(dir) EQ 1 then dir = dir else dir = 'data/'
   if keyword_set(right) then part = 'right' else part = 'left'
   if keyword_set(bin) then begin
         image=mrdfits(strcompress(dir+file+'_'+part+'_dewarped_bin'+string(fix(bin))+'.fits.gz',/remove_all))
   endif else begin
         image=mrdfits(dir+file+'_'+part+'_dewarped.fits.gz')
   endelse


   data = avg(image, 2);collapse the stack (vertical average)

   ;FIND THE WAVELENGTH THAT CORRESPONDS TO EACH COLUMN
   if ~keyword_set(bin) then bin=0
   L = calib_lambda(data, file, star,right=keyword_set(right), bin = bin) 
   
   ;FIT THE VERTICAL ENVELOPE
   yarr = findgen(n_elements(data[0,*]))
   env = fltarr(n_elements(data[*,0]), n_elements(data[0,*]))
   for i = ceil(n_elements(data[*,0])/2), n_elements(data[*,0])-1 do begin
       result = gaussfit(yarr, transpose(image(i,*)), nterms=3)
       env[i,*] = result
   endfor
   for i = floor(n_elements(data[*,0,0])/2),0,-1 do begin
       result = gaussfit(yarr, transpose(image(i,*)), nterms=3)
       env[i,*] = result
   endfor
   
   ;FIT THE FRINGES
   blurbit = fit_freq(image, L, env, file, right=keyword_set(right), bin = bin)
  
   line = 'rm '+init+'_*'+part+'.sav'
   spawn, line  ;this command will remove all files created to circumvent fitting madness

   if keyword_set(bin) then begin
      save, blurbit, filename = strcompress(dir+file+'_'+part+'_blur'+'_bin'+string(fix(bin))+'.sav',/remove_all)
   endif else begin
      save, blurbit, filename = strcompress(dir+file+'_'+part+'_blur.sav',/remove_all)
   endelse
   ;first element is frame to start with, second is the number to average through

   Sig = 0
   return, Sig
end



function calib_lambda, data, file, star, right = right, bin = bin, dir = dir
;Returns the wavelengths corresponding to the spectrum of the data.

  if keyword_set(dir) then begin
     dir = dir 
     dirr = strmid(dir, 0, strlen(dir) -strpos(dir, 'data')-2)+'reference_spectra/'
  endif else begin
     dir = 'data/'
     dirr = 'reference_spectra/'
  endelse

   spectr = avg(data[*,(size(data))[2] /4 : (size(data))[2] *3/4], 1) ;1d array that is the average spectrum of the data. The center half with the most flux was selected.

   if keyword_set(bin) then nbin = double(bin) else nbin = double(1)

   xdim = n_elements(spectr)
   xarr = findgen(n_elements(spectr))*nbin
   
   ;spectr2 = [spectr[0:130],spectr[150:n_elements(spectr)-1]]    ;remove the absorption line for fitting purposes (for fit_5th)  
   ;xarr2 = [xarr[0:130],xarr[150:n_elements(xarr)-1]]
   ;fit = poly_fit(xarr2, spectr2, 5, yfit = fit_5th)

   fit = poly_fit(xarr, spectr, 5, yfit = fit_5th)

   if keyword_set(bin) then begin
      spectrum_theo_g = spectrum_theo(star, xdim, .643, 68*!dtor, 1.e-6, fit_5th, 0, bin = bin) 
   endif else begin
      spectrum_theo_g = spectrum_theo(star, xdim, .643, 68*!dtor, 1.e-6, fit_5th, 0) 
   endelse
   ;all constants obtained from L212 in Elsa's 5.get_spatial_freq.i

   T = min(spectrum_theo_g)/min(continuum_detection(10,spectrum_theo_g))

   print, "Computing the X^2 map."

   cent_wav_arr = findgen(101)/1000. + 0.65 ;0.65 to 0.75 ;NEW: extended range .55 to .9 because not fitting incidence angle->changed back, .65 to .75
   inc_angl_arr = (findgen(101)/100.*20 + 50)*!dtor ;50 to 70
   x2 = fltarr(n_elements(cent_wav_arr), n_elements(inc_angl_arr))
   
   fname = strcompress(dirr+'Pickles_'+star+'.txt',/remove_all)
   readcol, fname, lambda_A1, flux1, /silent
   stel_spec = fltarr(n_elements(flux1), 2)
   stel_spec[*,0] = lambda_A1
   stel_spec[*,1] = flux1
   ;loading in spectral information now to avoid slowing the loop

   for i = 0, n_elements(cent_wav_arr)-1 do begin
      for j = 0, n_elements(inc_angl_arr)-1 do begin
         if keyword_set(bin) then begin
            spectrum = spectrum_theo(star, xdim, cent_wav_arr[i], inc_angl_arr[j], T, fit_5th, stel_spec, batch = 'boost', bin = bin)
         endif else begin
            spectrum = spectrum_theo(star, xdim, cent_wav_arr[i], inc_angl_arr[j], T, fit_5th, stel_spec, batch = 'boost')
         endelse
           x2[i,j] = total((spectr/continuum_detection(10,spectr) - spectrum/continuum_detection(10,spectrum))^2.)
          ;getting a "chi-square" based on a continuum normalized version of the spectra
      endfor
  endfor

   loc = max(where(x2 eq min(x2[where(x2 EQ x2)]))) ;avoiding NaN, don't ask me how they occasionally get in there
   loc = array_indices(x2, loc)

   cent_wav = cent_wav_arr[loc[0]] ;best guess of the central wavelength
   inc_angl = inc_angl_arr[loc[1]] ;best guess of the incidence angle

   ymap = fltarr(n_elements(x2[*,0]), n_elements(x2[0,*]))
   ymap[loc[0], loc[1]] = 10000
   print, "The best guess of the central wavelength is:", cent_wav
   print, "The best guess of the incidence angle is:", inc_angl

   if keyword_set(show) then begin
      loadct, 19
      display, x2, min = min(x2)*1.1, max = max(x2)*.9, xrange=[0.65, 0.75], yrange=[50, 70], background=!white, color=0, $
               title=textoidl('\chi^2')+' Map', $
               xtitle="Central wavelength ["+textoidl('\mu')+"m]", $
               ytitle="Incidence angle on the prism (deg)", charsize=2
      loadct, 0
      xyouts, cent_wav, inc_angl, "O", color = 0
   endif   

   if keyword_set(bin) then lambda_arr = wavelaw(xdim, cent_wav, inc_angl, bin = bin) else lambda_arr = wavelaw(xdim, cent_wav, inc_angl) 
   
   if keyword_set(bin) then begin
      if keyword_set(right) then begin
         save, cent_wav, inc_angl, filename =strcompress(dir+file+'_right_CLfit_bin'+string(fix(bin))+'.sav',/remove_all)
      endif else begin
         save, cent_wav, inc_angl, filename =strcompress(dir+file+'_left_CLfit_bin'+string(fix(bin))+'.sav',/remove_all)
      endelse 
   endif else begin
      if keyword_set(right) then begin
         save, cent_wav, inc_angl, filename =strcompress(dir+file+'_right_CLfit.sav')
      endif else begin
         save, cent_wav, inc_angl, filename =strcompress(dir+file+'_left_CLfit.sav')
      endelse 
   endelse
   ;saving the values of the fitted incidence angle and central wavelength

   return, lambda_arr
end




function spectrum_theo, star_type, xdim, cent_wav, inc_angl, T, fit_5th, batch_store, batch = batch, bin = bin, dirt = dirt
;Model function for the spectrum, comprising stellar spectrum and atmosphere transmission.

  common constants; pix_size, foc, tau, xvg, p_vg, glass, surf
 
  if keyword_set(dirt) EQ 1 then dirt = dirt else dirt = 'reference_spectra/'

  h=6.626068e-34
  c=299792458.

  if keyword_set(batch) EQ 1 then begin
  ;the batch keyword was added to run spectrum_theo more quickly in loops
     lambda_A = batch_store[*,0]
     Flux = batch_store[*,1]
  endif else begin
     readcol, dirt+'Pickles_'+star_type+'.txt', lambda_A, flux, /silent
  endelse
  ;lambda_A in A, Flux in erg.s-1.cm-2.A-1

  Flux *= 10. ;  Flux in W.m-2.µm-1
  lambda_um = lambda_A * 1.e-4 ;  wavelengths in µm
  ;nlambda_um = n_elements(lambda_um)
  
  ldif_um=lambda_um*0.
  for i = 0, n_elements(lambda_um)-2 do begin
     ldif_um[i] = abs(lambda_um[i+1] - lambda_um[i])
  endfor
  ldif_um[n_elements(ldif_um)-1] = ldif_um[n_elements(ldif_um)-2]
  ;creating an array of the delta lambdas

  Flux=Flux*lambda_um*ldif_um/(h*c) ;photons per m^2

  if keyword_set(bin) then L=wavelaw(xdim, cent_wav, inc_angl, bin = bin) else L = wavelaw(xdim, cent_wav, inc_angl) 
  ;Wavelengths corresponding to the given parameters
 
  xarr = findgen(n_elements(lambda_um))
  polyfit, xarr, Flux, 3, coeff, scoeff, yfit  ;fitting the global shape of the spectrum

  Flux = Flux/yfit*min(yfit) ;This is a change from Elsa's code that we believe will be more effective in flattening out the spectrum
  Flux=interpol(Flux, lambda_um, L) ; interpolate to get the flux at each pixel

  Tatm=T_atm(L, dir = dirt)     ; ;returns the percentage transmission of each wavelength

  Flux=1.e-6*Flux*tau*Surf*Tatm*fit_5th/max(fit_5th)*T ; //make flux curvy like the observed spectrum

  return, flux
end




function wavelaw, xdim, cent_wav, inc_angle, bin = bin
;generates the wavelength values corresponding to the pixels
  common constants; pix_size, foc, tau, xvg, p_vg, glass, surf

  A=60.*!dtor ;        // prism angle (rad)  
  ld = findgen(1001)/1000.*0.6 + 0.4 ;0.4 to 1

  ;Thorlabs prisms used during the 2011 runs are in F2 and not SF2
  ; Returns the refractive index of the F2 glass for wavelengths L (in um)
  ; It uses Sellmeier formula for F2 glass.
  ; Source: http://refractiveindex.info/?group=SCHOTT&material=F2
  n = sqrt( 1 $
            + 1.34533359*ld^2/(ld^2-0.00997743871) $
            + 0.209073176*ld^2/(ld^2-0.0470450767) $
            + 0.937357162*ld^2/(ld^2-111.886764) ) 
  
  ;inc_angle *= !dtor ;moving this step up because it's being looped
  
  ;Deviation as a function of lambda (n(lambda))
  D = inc_angle+asin(n*sin(A-asin(sin(inc_angle)/n)))-A ;
  ;if you are trying to work out the math here, there is a !pi term missing, but it won't affect the tan
 
 ;null deviation : centered on l0
  ind = min(where(abs(ld-cent_wav) eq min(abs(ld-cent_wav))))
  if ind ne -1 then begin
     D0 = D[ind]     
     D -= D0
  endif else begin
     print, "Problem! Problem!"
     stop
  endelse

  ;position on the CCD as a function of lambda
  x=foc*tan(D); (mm)
  x_px = x/p_det ; (mm)

  if keyword_set(bin) then nbin = double(bin) else nbin = double(1) ;making sure shortening due to binning is accounted for

  X = findgen(xdim)*nbin-floor((xdim*nbin)/2.)
  L=interpol(ld, x_px, X) 
  L = L[sort(L)]
  return, L
end



function t_atm, L, dirt = dirt 
;finds the transmission as a function of wavelength
  if keyword_set(dirt) EQ 1 then dirt = dirt else dirt = 'reference_spectra/'

  readcol, dirt + 'Tatm_short.txt', lambda, transmission, /silent  
  tatm = interpol(transmission, lambda, L)

return, tatm
end



function fit_freq, data, L, env, file, right = right, bin = bin, dir = dir
; Returns the wavenumbers resulting from the fitting of interferograms.

 ;This is the remaining section of code to be completed...look
 ; at the central part of Elsa's 5.Get_spatial_freq.i for her code.
  common constants ; pix_size, foc, tau, xvg, p_vg, glass, surf
  common constants2, nl, np, nb, ni
  sizdat = size(data)
  nl = sizdat[1]
  np = sizdat[2]
  ni = sizdat[3]

 
  if keyword_set(dir) then dir = dir else dir = 'data/'

  frames = keep_contrast(data)
  if frames[0]+frames[1] GE ni then frames[0] = ni-frames[1]-1
  fitnum = nl/10
 
  data = data[*,*, frames[0]:frames[0]+frames[1]]
  ;When averaging over many frames the fringes may become quite blurred due to
  ;poor seeing during observation. Keep contrast checks how much blurring would
  ;develop ranging over various numbers of frame and finds the number of frames
  ;that maximizes SNR while maintaining contrast

  datam = avg(data, 2)
  
  B = get_UV(xvg)*p_VG ;baselines in mm
  nB = n_elements(B)

  s = 1./L

  bgrnd = dblarr(nl) ;background flux
  mu0 = bgrnd+.1 ;uncorrelated flux
  Re = dblarr(nB, nl)+.1 ;real part
  Im = Re ;imaginary part
  Sigg = s ;wavenumber
  Sig = Sigg*0.
  stdev_fit = bgrnd ;standard deviation in fit of fringes
  pVGfit = bgrnd ;fit of pitch value
  XVGfit = dblarr(9,nl) ;fit of baseline locations in terms of pitch
  X2 = bgrnd                        ;chi square values of fit
  x = dindgen(np) - p_det*(np/2.-1) ;pixel location in mm
  fit = dblarr(nl, np)
  bgrnd = bgrnd +10
  param_err = dblarr(3+2*nB, nl) ;takes the standard deviations of the various parameters
  XVG_err = XVGfit
  ;setting up the arrays to take the parameters of the fit and various errors
  
  print, "Fitting the interferograms"
  mx = max(avg(datam, 1), nlm)
  ;fitting from the brightest pixel with respect to columns

  a=[Sigg[nlm],bgrnd[nlm], mu0[nlm], Re[*,nlm], Im[*,nlm]]
  ;setting initial parameters

  if keyword_set(bin) then fitnum /= double(bin) 
  ;making sure that our current fitting distribution is maintained 

  for i=nlm, fitnum*floor(nl/fitnum)-1, fitnum do begin
     ;fitnum was put in place due to the fact that it took far too long to fit every
     ;column, so an optimal number of columns to fit was found and used instead

     param=[Sigg[i],bgrnd[i], mu0[i], Re[*,i], Im[*,i]] 
     parinfon = [replicate({mpprint: 0}, n_elements(param))] ;avoids printing so many lines
     col_env = double(env[i,*])
     save, col_env, filename = init+'_'+strcompress(string(Sigg[i]),/remove_all)+'_env'+part+'.sav'
     
     ;Couldn't force mpfitfun to take the envelope values any other way

     a = mpfitfun('fringes',x, reform(datam[i,*]), err, param, weights = x*0+double(1),nprint = 1000, parinfo = parinfon, perror = param_stdev)
     
     Sig[i] = a[0]
     Re[*,i] = a[3:nB+2]
     Im[*,i] = a[nB+3:2*nB+2]
     bgrnd[i] = a[1]
     mu0[i] = a[2]
     param_err[*,i] = param_stdev
     stdev_fit[i] = sqrt((total((reform(datam[i,*])- fringes(x,a))^2))/(n_elements(x)-n_elements(a)))


     ;FITTING OF FIBER POSITIONS IN THE V-GROOVE (XVG)
     param2 = double([XVG, a[0]])
     col_envf = double([reform(env[i,*]),a])
     save, col_envf, filename = init+'_'+strcompress(string(a[0]),/remove_all)+'_envf'+part+'.sav'
     parinfof = [replicate({mpprint: 0}, n_elements(param2))]
   
     b = mpfitfun('fringes_fib',x, datam[i,*], err, param2, functargs = fcnaf, nprint = 1000, parinfo = parinfof, perror = param2_stdev)
    
     X2[i] = total((reform(datam[i,*])- fringes_fib(x,b))^2)
     XVGfit[*,i] = b
     XVGfit[*,i] += avg(XVG)-avg(XVGfit[*,i])
     fit[i,*]=fringes_fib(x,b)
     XVG_err[*,i] = param2_stdev
  endfor
  for i = nlm, round(fitnum*(nlm/fitnum-floor(nlm/fitnum))), -1.*fitnum do begin
    
     param=[Sigg[i],bgrnd[i], mu0[i], Re[*,i], Im[*,i]]
     parinfon = [replicate({mpprint: 0}, n_elements(param))]
     col_env = double(env[i,*])
     save, col_env, filename = init+'_'+strcompress(string(Sigg[i]),/remove_all)+'_env'+part+'.sav'
     a = mpfitfun('fringes',x, reform(datam[i,*]), err, param, nprint = 1000, parinfo = parinfon, perror = param_stdev)
          
     Sig[i] = a[0]
     Re[*,i] = a[3:nB+2]
     Im[*,i] = a[nB+3:2*nB+2]
     bgrnd[i] = a[1]
     mu0[i] = a[2]
     param_err[*,i] = param_stdev
     stdev_fit[i] = sqrt((total((reform(datam[i,*])- fringes(x,a))^2))/(n_elements(x)-n_elements(a)))

     param2 = double([XVG, a[0]])
     col_envf = double([reform(env[i,*]),a])
     save, col_envf, filename = init+'_'+strcompress(string(a[0]),/remove_all)+'_envf'+part+'.sav'
     parinfof = [replicate({mpprint: 0}, n_elements(param2))]
  
     b = mpfitfun('fringes_fib',x, datam[i,*], err, param2, nprint = 1000, parinfo = parinfof, perror = param2_stdev)
    
     X2[i] = total((reform(datam[i,*])- fringes_fib(x,b))^2)
     XVGfit[*,i] = b
     XVGfit[*,i] += avg(XVG)-avg(XVGfit[*,i])
     fit[i,*]=fringes_fib(x,b)
     XVG_err[*,i] = param2_stdev
  endfor
  
  interp= where(Sig NE 0)

  Sig = interpol(Sig[interp], interp, dindgen(nl),/spline)
  bgrnd = interpol(bgrnd[interp], interp, dindgen(nl),/spline)
  mu0 = interpol(mu0[interp], interp, dindgen(nl),/spline)
  for i = 0, nB -1 do begin
     Im[i,*] = interpol(reform(Im[i,interp]), interp, dindgen(nl),/spline)
     Re[i,*] = interpol(reform(Re[i,interp]), interp, dindgen(nl),/spline)
  endfor
  for i = 0, 8 do begin
     XVGfit[i,*] = interpol(reform(XVGfit[i,interp]), interp, dindgen(nl),/spline)
  endfor
  for i = 0, np-1 do begin
     fit[*,i] =  interpol(reform(fit[interp,i]), interp, dindgen(nl),/spline)
  endfor
  ;interpolating through all of the arrays containing fitted values to get information for each column

  cols_fit = where(total(param_err,1) NE 0)
  ;saving the numbers of the columns that were fit for reference later
 
  param_err = param_err[*, where(total(param_err, 1) NE 0)]
  XVG_err = XVG_err[*, where(total(XVG_err, 1) NE 0)]
  ;removing empty rows

  XVG2 = avg(XVGfit, 1)


  ;FINDING POLYNOMIAL LAW AND CHECKING GOODNESS OF FIT
  x = dindgen(nl)
  a = dblarr(3)+1.
  w = x
  w(where(avg(datam,1) GT median(avg(datam,1))*.9)) = 1 ;weighting the fit
  x = x + 1
  Sigfit = curvefit(x, Sig, w, a, Sig_stdev_fit, /noderivative, function_name = 'poly2pro')
  
  if total(abs(1/Sigfit-L)) GT nl*.01 then begin
     print, "Wavelength fit not good enough"
  endif
  ;checking fit of wavelength
  

;WORK IN YOUR PART VARIABLE HERE
  ;SAVING THE RESULTS:  
  file=repstr(file, '.fits.gz', '')
  if keyword_set(bin) then begin
     if keyword_set(right) then begin
        save_name=strcompress(dir+file+ '_right_bin'+string(fix(bin)), /remove_all)
     endif else begin
        save_name=strcompress(dir+file+ '_left_bin'+string(fix(bin)), /remove_all)
     endelse
  endif else begin
     if keyword_set(right) then begin
        save_name=strcompress(dir+file+ '_right', /remove_all)
     endif else begin
        save_name=strcompress(dir+file+ '_left', /remove_all)
     endelse
  endelse  
  ;creating the save name skeleton
   
  Pvgfit = Sigfit*(1./Sig)*p_vg ;getting p_vg as a function of wavelength (non-physical)

  ; The fiber positions (XVG)
  filex=strcompress(save_name+"_XVGfit.txt", /remove_all) 
  openw, 110, filex
  printf,110, XVGfit[*,0], XVGfit[*,1], XVGfit[*,2], XVGfit[*,3], XVGfit[*,4], XVGfit[*,5], XVGfit[*,6], XVGfit[*,7], XVGfit[*,8] ;this will allow for the usage of readcol for each baseline-->see method in xblock_comp.pro
  close, 110

  print, "XVGfit results saved in: "      
  print, "   "+save_name+"_XVGfit.txt"   
  
  ;The fitted wavelengths and polynomial law
  ;this is a log for human purposes, will next save a file for program purposes
  Lfit = 1./Sig

  filep=strcompress(save_name+"_Spatial_freq.txt", /remove_all)           
  openw, 110, filep  
  printf, 110, "p_VG (mm):";
  printf, 110, p_VG ;
  printf, 110, ""
  printf, 110, "Information from interferogram fit:"
  printf, 110, "Fitted columns: "
  printf, 110, cols_fit
  printf, 110, "Wavenumbers (um^-1):"
  printf, 110, Sig ;fitted on interferograms
  printf, 110, "Wavenumber Standard Deviation of fit (um^-1):"
  printf, 110, reform(param_err[0,*])
  printf, 110, "Wavelengths (um):";
  printf, 110, Lfit ;fitted on interferograms  
  printf, 110, "Background Flux:"
  printf, 110, bgrnd
  printf, 110, "Background Flux stdev:"
  printf, 110, reform(param_err[1,*])
  printf, 110, "Uncorrelated Flux:"
  printf, 110, mu0
  printf, 110, "Uncorrelated Flux stdev:"
  printf, 110, reform(param_err[2,*])
  printf, 110, "avg(XVGfit,1):" ;
  printf, 110, XVG2;
  printf, 110, "Further information (XVGfit, Re, Im, fit, errors on all, etc.) saved in:"
  printf, 110, strcompress(save_name+"_Spatial_freq_extr.sav", /remove_all)
  printf, 110, ""
  printf, 110, "Information from polynomial fit: "
  printf, 110, "Wavenumbers (um^-1):"
  printf, 110, Sigfit
  printf, 110, "Wavenumber stdev (um^-1):"
  printf, 110, Sig_stdev_fit ;fitted on interferograms
  printf, 110, "V-groove pitch (um) based on polynomial law" ;
  printf, 110, Pvgfit ;  // new 18.06 : save the Pvg depending on pixel, from polynomial law
  close, 110 

  otro_err = param_err[3:*, *]
  save, XVGfit, XVG_err, fit, stdev_fit, Re, Im, otro_err, filename = strcompress(save_name+"_Spatial_freq_extr.sav", /remove_all)
  ;saving everything else we got out of the fit in case we want it later

  print, "Spatial frequencies results saved in: " ;
  print, "   "+save_name+"_Spatial_freq.txt" ;
  
  Pvg = Pvgfit
  XVG = XVG2
  save, p_VG, XVG, Sig, Pvg, filename = strcompress(save_name+'_Spatial_freq.sav', /remove_all)
  ;The desired information for get_best_params is saved in this file

  return, frames
end


function get_UV, xvg
;finds the values of all of the baselines

  nx = n_elements(xvg)
  r = (nx-1)*nx ;in reality uv will be half this number
  uv = dblarr(r)
 
  for i = 0, nx-2 do begin
     for j = i+1, nx-1 do begin
     uv[j+i*(nx-1)] = abs(xvg[i] - xvg[j])
     endfor
  endfor
  ;finding all of the separations

  uv = uv(where(uv NE 0))
  return, uv
end



function fringes, x, a
;creates model function for fringes as the sum of cos and sin functions
;for a good description of the background behind this function read the
;data reduction section of the first on-sky results paper (Huby, 2012)
common constants
common constants2
  restore, init+'_'+strcompress(string(a[0]),/remove_all)+'_env'+part+'.sav'
  envl = col_env
  B = get_UV(xvg)*p_VG            ;baselines in mm
 
  nB = n_elements(B)
  sig=a[0]*1000
  bgrnd=a[1]
  mu0 = a[2]
  Re = a[3:nB+2]
  Im = a[nB+3:2*nB+2]

  ;setting up arrays of parameters to be used in the functions 
  ;reform functions are used to ensure proper multiplication of "matrices"

   fct_cos = cos(2*!pi*x*sig##reform(B, nB, 1)/foc)
   fct_sin = sin(2*!pi*x*sig##reform(B, nB, 1)/foc)

  ;interior of functions = 2*pi*(location on ccd)*(spatial frequency), where
  ;spatial frequency = wavenumber*baseline/(distance to image plane=foc)
 
  y = (fct_cos##Re+fct_sin##Im+mu0)*envl+bgrnd 
  y = reform(y)
 
  return, y
 
end



function fringes_fib, x, f
;models the fringes but instead makes the baselines a free parameter so
;that they can be fitted using the fitted parameter values from the
;first pass
common constants
common constants2
  restore, init+strcompress(string(f[-1]),/remove_all)+'_envf'+part+'.sav'
  envf = col_envf

  as = envf[np :n_elements(envf)-1]
  envf = envf[0:np-1]
  xvg1=f[0:-2]
  bs = get_UV(xvg1)*p_VG
  
  sig=as[0]*1000
  bgrnd=as[1]
  mu0 = as[2]
  Re = as[3:nB+2]
  Im = as[nB+3:2*nB+2]
  
  fct_cos = cos(2*!pi*x*sig##reform(bs, nB, 1)/foc)
  fct_sin = sin(2*!pi*x*sig##reform(bs, nB, 1)/foc)
  
  ;interior of functions = 2*pi*(location on ccd)*(spatial frequency), where
  ;spatial frequency = wavenumber*baseline/(distance to image plane=foc)
 
  y = (fct_cos##Re+fct_sin##Im+mu0)*envf+bgrnd 
  y = reform(y)
  return, y
 
end



pro poly2pro, x, a, y
;the second degree polynomial for curvefit for the fitting of the polynomial law
  y = a[2]*x^2 +a[1]*x+a[0]
end



function keep_contrast, data
;this function checks the level of blurring that occurs from averaging
;of an image through the frames and returns the optimal number of frames
;to average through
 
;METHOD 1:
;0*: Edited so that the program goes through and tries starting at each frame
;1: Goes through the columns and finds the ratio of the max value to the
;value of the closest loccal min to get an idea of contrast
;2: Finds ratio of the brightest pixel to the average pixel in each
;column, to get a rough idea of SNR
;3: Finds worst column (where difference between contrast and SNR is minimized)->blur score is that difference*SNR(column)
;4: Finds number of frames at which the blur score is maximized

  meta_score = fltarr((size(data))[3],2)
  quip = "search"
  for q = 1, (size(data))[3]-1 do begin
     score = fltarr((size(data))[3])
     for i = q, (size(data))[3]-1 do begin
        data_slice = data[*, *, q-1 :i] ;takes the slice of frames to work with
        data_av = avg(data_slice, 2) ;averages the slice to get our average frame we will be checking the blurring of
        diff = fltarr((size(data))[1],2)
        for k = 0, (size(data))[1]-1 do begin 
           bright = max(data_av[k,*], maxloc)
           contrast = fltarr((size(data))[2])
           for j = maxloc-1, 0, -1 do begin
              value = data_av[k, j]
              if value LE 0 then begin
                 value = 1.     ;minimum number of counts possible, stops explosion
              endif
              contrast[j] = bright/value ;column positions in contrast are the columns defining the range tested      
           endfor
           for j = maxloc+1, (size(data))[2]-1 do begin
              value = data_av[k, j]
              if value LE 0 then begin
                 value = 1.     ;minimum number of counts possible, stops explosion
              endif
              contrast[j] = bright/value ;column positions in contrast are the columns defining the range tested      
           endfor
           
           troughs = where(contrast[*] GT shift(contrast[*], 1) AND contrast[*] GT shift(contrast[*,0], -1))
           if (size(data))[2]-1 LT maxloc+10 then r = -1. else r = 1.
           if total(troughs) EQ -1 then troughs = maxloc+10*r
           trough = min(abs(maxloc - troughs)) 
           trough = troughs[where(troughs EQ maxloc + trough OR troughs EQ maxloc - trough)] ;locating the closest minimum to the central peak
           if n_elements(trough) GT 1 then  trough = trough[0] ;forcing IDL to pick the left side if left and right have the minima in the same relative location
           contrast = contrast[trough] ;finding the contrast at that peak
           
           if maxloc - trough GT 0  then begin
              normal = avg(data_av[k, trough : maxloc])
           endif else begin
              normal = avg(data_av[k, maxloc :trough])
           endelse
           ;figuring out whether the minimum is on the left or right hand side and finding the respective "normal" value

           SNR = bright/normal
           diff[k,0] = contrast-SNR
           diff[k,1] = 1/contrast
        endfor

        worst = min(diff[*,0], minloc)
        score[i] = diff[minloc,0]*diff[minloc,1] ;finding the worst column and getting its score
        if score[i] LT score[i-1] then goto, escape ;jumping out of the loop if best score has been achieved
     endfor
     
     escape: meta_score[q,0] = score[i-1]
     meta_score[q,1] = i - (q-1) ;getting the scores for each starting frame as well
  endfor

  best = max(meta_score[*,0], mloc)
  print, "The best number of frames to average through is " + string(meta_score[mloc,1]+1) + " with a score of " + string(best)
  frnum = [mloc,meta_score[mloc,1]+1]
  
  return, frnum
end


function continuum_detection, nmed, spectr
;name was retained from Elsa's code, but what this function
;really does is boxcar smooth the spectrum

   nl = n_elements(spectr)   
   spectr_med = spectr

   for i = round(nmed/2.), nl-round(nmed/2.)-1 do begin
       spectr_med[i] = median(spectr[i-round(nmed/2.):i+round(nmed/2.)])
   endfor

   return, spectr_med
end

function get_spatial_freq_upd, file_stack, dirsave, right = right, show = show
;this get_spatial_freq function has been updated to take into account
;changes that have been implemented in Elsa's code

                                ;FIRST READ IN ALL CORRECTED FILES AND
                                ;CREATE "AVERAGE" IMAGE, USE TO FIND ENVELOPE
  print, "Reading all corrected files"
  init = strmid(file_stack[0], 0, 4)
  date2 = strmid(file_stack[0], 5, 8)
  star = strmid(file_stack[0], 14, strpos(file_stack[0], '_', 14)-14)

  if keyword_set(right) EQ 1 then part = right else part = left
  file_avg = file_search(strcompress(dirsave+init+'all_data_avg_corr_'+part+'.fits.gz',/remove_all))
  if file_avg NE "" then begin
      data = mrdfits(file_avg[0])
      print, "Averaged file located"
  endif
  if file_avg EQ "" then begin
     img1 = avg(mrdfits(file_stack[0]), 2)
     siz_img1 = size(img1)
     data = fltarr(siz_img1[1], siz_img1[2], n_elements(file_stack))
     data[*,*,0] = img1
     for i = 1, n_elements(file_stack)-1 do begin
        data[*,*,i] = avg(mrdfits(file_stack[i]), 2)
     endfor
     data = avg(data,2)
     writefits, strcompress(dirsave+init+'all_data_avg_corr_'+part+'.fits.gz',/remove_all), data
     datam = data
  endif

  siz_files = size(file_stack)
  if siz_files EQ 1 then file_name = file_stack[round(siz_files[1]/2.)]
  if siz_files EQ 2 then file_name = file_stack[round(siz_files[1]/2.), round(siz_files[2]/2.)]
  if siz_files EQ 0 then file_name = file_stack

  print, "Reading in file: ", file_name
  data = mrdfits(strcompress(dirsave+file_name+'.fits.gz',/remove_all))
  
  siz_data = size(data)
  nl = siz_data[1]
  np = siz_data[2]
  ni = siz_data[3]

  if ni GE 100 then data= data[*,*,0:99]
  
  max_test = fltarr(np, ni)
  for i = 0, nl-1 do begin
     for j = 0, ni -1 do begin
        max_test[i,j] = max(data[i, *, j])
     endfor
  endfor
  max_test = avg(max_test)
  if max_test LE 5000 then data = data[*,*,0:14]
  ;testing the number of counts to work around the case of low SNR

  ;FITTING THE ENVELOPE
  env = datam
  xp = findgen(np)+1
  nlm = round(nl/2.)
  for i = nlm, nl-1 do begin
     env[i, *] = gaussfit(xp, datam[i,*], nterms = 3)
  endfor
  for i = nlm-1, 0, -1 do begin
     env[i, *] = gaussfit(xp, datam[i,*], nterms = 3)
  endfor

  if keyword_set(show) then display, env

  save_name = strcompress(dirsave+init+'_'+date2+'_'+star+'_'+part,/remove_all)

  ;FITTING THE SPECTRUM BASED ON RANDOM MIDDLE FILE
  L = calib_lambda(data, file_name)
  Sig = 1./L
  
  ;FITTING THE FRINGES BASED ON RANDOM MIDDLE FILE
  Sig = fit_freq(data, L, env, file_name, 18)

  return, Sig
end
