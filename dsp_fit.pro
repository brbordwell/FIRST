pro get_spatial_freq, file, RIGHT=right, BIN=bin, DIR=dir
;This is the global function. It adjusts the wavelengths to the column
;location and then fits wavelengths to the fringes. 

;NOTE: ELSA'S UPDATED VERSION OF THIS FUNCTION IS LOCATED AT THE
;END OF THIS CODE
;+
;
;FUNCTION: GET_SPATIAL_FREQ, FILE, [.../RIGHT,BIN='binsize',DIR='path']
;
;PURPOSE:
;  To set up the common constants required in subsequent functions, and
;  prep for the fitting of the spatial frequencies prior to running the
;  fitting procedure FIT_FREQ.
 
;INPUT:
;  file = The file tag of the dewarped image file.
;         Example: 1040_20111014_Aldebaran_0
;       
;OUTPUT:
;  (None) = See the files saved in FIT_FREQ
;
;NOTES:
;  -For individual methods, reference the individual functions.
;  -The common constants are:
;     p_det = pixel size on the detector (mm)
;     foc = focal length (mm)
;     XVG = fiber positions in the v-groove
;     p_VG = v-groove pitch (mm)
;     glass = glass is the glass of the prism (BK7 for version 9 or SF2
;             for version 18. For 9-fibers, glass=1; for 18, glass=2)
;     surf = surface collecting photons (m^2)
;            Formula: (Npupil_used)*pi*Dtel^2/(Npupil*4?)
;     tau = the integration time for the image (s)
;     tag = the composed filetag containing the sequence number, date,
;           star, file number and side of the image.
;           Example: 1040_20111014_Aldebaran_0_left  
;     ni,np,nl =  ;number of images, pixels and wavelengths in the image
;  -Changed the fitting of the vertical envelope to the method employed
;   in the updated GET_FIBERS code
;-


  ;Setting up a common block with required constants...
  common constants, p_det, foc, p_vg, xvg, glass, surf, tau, tag, rb, hdr, star, part, nl, np, ni
  

  if keyword_set(right) then part = 'right' else part = 'left'
  if keyword_set(dir) then dir = dir else dir = 'data/'
  if ~keyword_set(bin) then rb = '' $
  else rb = strcompress('_bin'+string(fix(bin)),/remove_all)
  tag = dir+file+'_'+part
  star = strmid(file, 14, strpos(file, '_', /reverse_search)-14) 

  p_det = 0.016  &  foc = 150.  &  p_VG = 0.037  
  XVG = [2,3,7,14,27,29,37,43,46] 
  glass = 'F2'  
  surf = 9.*!pi*3.^2/(37*4.) 
  

  




  ;Loading in and collapsing the image...
  filename = tag+'_dewarped'+rb+'.fits.gz'
  if file_search(filename) EQ "" then begin
     print, "File not found. Cannot perform wavelaw fit."
    goto, theend
  endif
  image = mrdfits(filename, 0 ,hdr)
  data = avg(image, 2)
  siz = size(image)  &  nl = siz[1]  &  np = siz[2]  &  ni = siz[3]

  stat = strcompress(fxpar(hdr, 'STATUS'),/remove_all)
  if stat EQ 'REJECT' then begin
     num = strpos(file,'_',/reverse_search)
     num = strcompress(strmid(file,num+1,5),/remove_all)
     (scope_varfetch('stack_status_token',level=-1))[fix(num)]=0
     print, "Image is reject, skipping DSP_fit"
     goto, theend
  endif




  tau = fxpar(hdr,'SHUTTER')/1000.
  print, "The integration time is ", $
         strcompress(fix(string(tau*1e3)), /remove_all), $
         " milliseconds. If this is not correct, please update the code."




  ;Find the wavelength corresponding to each column...
  if ~keyword_set(bin) then bin=0
  fil = file_search(tag+'_CLfit2.sav')
  if fil NE "" then restore, fil else L=calib_lambda(data, bin=bin) 
  if n_elements(L) EQ 1 then goto, theend



  ;Fit the vertical envelope...
  x_out = dindgen(np)  
  env_vars = dblarr(nl,4)

  
;  -As was observed in GET_FIBERS, envelope fitting is a delicate task.
;   Originally I began with trying to replace GAUSSFIT (GF) with a
;   smoothed and interpolated (at the ends) (SIP) version of the
;   envelopes, but upon considering what this replacement would require
;   further along in the pipeline (FFT of the envelope rather than using
;   fitted Gaussian variables) I decided that it was worth more to look
;   deeper at possible methods of fitting a Gausian. I found the
;   function MPFITPEAK (MPP), and compared:
;     MPP fit of the raw column
;     MPP fit of the SIP version
;     GF fit of the raw column
;   While no one function was completely ideal all of the time, both raw
;   column fits together provided sufficiently alike fits to the data to
;   be usable in the pipeline. The simplest way to choose the best
;   function was to sum the absolute value of the differentials for each
;   column. 
;   NEWS: Upon playing with things even further, I found that MPP and GF
;   give almost identical results when presented with the same estimates
;   of the coefficients for a Gaussian with a constant vertical
;   offset. This, luckily, obviates the need for a differential comparison.
;   I chose to use MPP because in case there was any deviation, MPP
;   tends to be more robust given (no) bad estimates.



  if file_search(tag+'_env.sav') EQ "" then begin
     for i = 0, nl-1 do begin
        col = data[i,*]  &  col -= min(col)
        pop = []  &  for j = 0, np-1 do pop = [pop,intarr(col[j]/10.+1)+j]
        coeff = [max(col)-min(col), avg(pop), stddev(pop), min(col)]
        mpp_coeff = coeff 
        
        mpp_col = mpfitpeak(x_out, col, mpp_coeff, nterms=4, estimates=coeff)
        env_vars[i,*] = mpp_coeff
     endfor
     save, env_vars, filename = tag+'_env.sav'
  endif


  ;Fit the fringes (2nd pass wavelaw fit)...
  if file_search(tag+'_Spatial_freq.sav') EQ "" then begin
     hedr = fit_freq2(image, L, env_vars) ;, file, dir)
     hedr = hedr[where(strpos(hedr,strjoin(strarr(80)+' ')) EQ -1)]
     filename = tag+'_dewarped'+rb+'.fits'
     writefits, filename, image, hedr, /compress
  endif
  theend: print
end






function DSP_fringes, x, c
;+
;
;FUNCTION: DSP_FRINGES(X,C)
;
;PURPOSE:
;  Returns a function model for the fringe pattern (which is composed of
;  the sum of cosine and sine functions)
;
;INPUT:
;  x = The positions (in mm) over which the function ranges on the CCD.
;  c = The various constants of the function to be fit, which are varied
;      by mpfitfun.
;
;OUTPUT:
;  y = The fringe pattern for the given values of x.
;
;NOTES:(None)
;
;-

  ;Setting up constants...
  common constants
  B = get_UV(XVG)*p_VG  &  nB = n_elements(B)
  sig = c[0]*1000  &  sf = sig/foc
  back = c[1]  &  mu0 = c[2]  &  re = abs(c[3:nB+2])
  Gc = 1./sqrt(4*!pi^2*(c[-1]*p_det)^2)


  ;Positive and negative peaks...
  fct = dblarr(np)
  a = [1., sf, Gc, back] 
  for k = 0, nB-1 do begin
     a_p = a  &  a_p[1] *= B[k]  &   a_m = a_p  &  a_p[1] *= -1
     fct = [[fct], [gaussian(x,a_p)^2+gaussian(x,a_m)^2]]
  endfor
  y = fct[*,1:*]

  
  ;0 frequency peak...
  a_0 = a  &  a[1] *= 0  &  fct_0 = gaussian(x,a_0)^2


  ;Composing the final function...
  for d = 0, np-1 do y[d,*] *= re ;QUESTIONS
  for d = 0, nB-1 do y[*,d] += mu0^2*fct_0 ;QUESTIONS
  y = total(y,2)

  return, y
end








function fit_freq2, data, L, env
;+
;
;FUNCTION: FIT_FREQ2(DATA,L,ENV)
;
;PURPOSE:
;  Returns the wavenumbers resulting from the fitting of the DSP
;
;INPUT:
;  data= the cube of data to be fitted (nl x np x ni dims)
;  L = the wavelaw obtained from the spectrum fitting
;  env = the fitted envelope variables for a Gaussian with a vertical
;        offset (nl x 4)
;
;OUTPUT:
;
;NOTES:
;  p_VG = the V-groove pitch
;  p_det = the pixel size on the detector
;  star_type = defining the stellar spectrum
;  save_name = the name the graphs will be saved under
;  version = "9" or "18", defining the version of the instrument in use
;  NOTE: tau (integration time), p_det (pixel size on the detector, mm)
;  and XVG (fiber positions in the V-groove) must be externalized
;  variables
;
;-
  common constants ;grabbing externalized variables

  XVG = double(XVG)
  datam = avg(data, 2) ;averaging through the stack of images
  B = get_UV(XVG)*p_VG  &  nB = n_elements(B)
  s = 1./L 


  ;Computing the mean DSP...
  file_dsp = (file_search(tag+'_all_data_DSPm.fits.gz'))[0]
  if file_dsp EQ "" then begin
     print, "Computing DSP..."
     DSParr = data*0. 
     for i = 0, ni-1 do begin
        dsp = abs(fft(reform(data[*,*,i]), dimension = 2))^2
        DSParr[*,*,i] = dsp
        if i/(ni/5.)-i/fix(ni/5.) EQ 0 then print, "i = ", i
     endfor
     DSPm = avg(DSParr,2) 

     hedr = hdr
     fxaddpar, hedr, 'STATUS', fxpar(hdr,'STATUS'),' AVG SPECTRAL POWER DENSITY DIST.'
     hedr = hedr[where(hedr NE strjoin(strarr(80)+' '))]
     writefits, tag+'_all_data_DSPm.fits', DSPm, hedr, /compress
  endif else begin
     print, "DSP already computed"
  
     print, "Loading DSP..."
     DSPm = mrdfits(file_dsp)
  endelse
stop

  ;Subtracting the photon noise component off the DSP...
  val = reform(median(DSPm[*,150:-150],dimension = 2))
  for i = 0, n_elements(nl)-1 do DSPm[*,i] -= val 
  ;INDICES WAT????????????


  ;Setting up the parameters for the fit...
  ;(the derivations behind these values will be explained in the notes)
  Back = dblarr(nl)
  amp = sqrt(2D*!pi)*env[*,0]*env[*,2]
  Re = dblarr(nB, nl)  &  for i = 0, nB-1 do Re[i,*] += amp
  Mu0 = amp + env[*,3]  


  Sig = double(s)  &  Pvg = Sig
  X2 = Back  &  stdev_fit = Back
  Re_err = Re
  param_err = dblarr(3+nB, nl)

  x = (dindgen(np)-np/2.+1)*p_det
  FIT = dblarr(nl,np)




  ;Fitting the DSP...
  print, "Fitting the DSP..."
  dummy = max(median(datam, dimension=2), nlm) ;Starting with the max SNR


  DL = 8. ;only fitting one spectral channel
  fe = 1./p_det ;sampling frequency

  if np/2. EQ fix(np)/2 then u = x*fe^2/np else u = (x*fe-1)*fe/(np-1) 
  u = shift(u, floor(np/2.)+1)
  Sigg = Sig*0.
  xf = [findgen(fix(nlm/DL)+1)*DL,findgen(fix(nl-1-nlm)/DL+1)*DL+nlm]+1
  xf = xf[uniq(xf)] 

  print, "Beginning from ", nlm, " and increasing..."
  for q = 0, 1 do begin
     if q EQ 0 then set = [nlm, nl-1, DL] else set = [nlm, 0, -DL]
     if q EQ 1 then print, "...and decreasing..."
     for j = set[0], set[1], set[2] do begin
        if j/10. EQ fix(j)/10 then print, "  Channel ", j, "/", nl




        col = reform(DSPm[j,*])
        param = [Sig[j], Back[j], Mu0[j], Re[*,j], env[j,2]]        
        parinfon = [replicate({mpprint: 0, fixed: 0}, nB+4)] 
        parinfon[3:*].fixed = 1 ;1st 3 parameters are free
        
        a = mpfitfun('DSP_fringes', u, col, err, param, nprint = 1e3, $
                     weights=u*0+1D, parinfo = parinfon, $
                     perror = param_stdev, status = issues)
        
        if issues GE 1 and issues NE 5 then begin 
           param_err[*,j] = param_stdev[0:-2]  
        endif
        dof = np-(nB+4)+1


        param = a  &  parinfon.fixed = ~parinfon.fixed
;        parinfon[0:2].fixed = 1
        parinfon[-1].fixed = 1  ;fixed->free (except env) 

        a = mpfitfun('DSP_fringes', u, col, err, param, nprint = 1e3, $
                     weights = u*0+double(1), parinfo = parinfon, $
                     perror = param_stdev2, status = issues)
        
        Sigg[j] = a[0]  &  Back[j] = a[1]  &  Mu0[j] = a[2]
        Re[*,j] = a[3:nB+2]  &  
        if issues GE 1 and issues NE 5 then begin 
           Re_err[*,j] = param_stdev2[3:nB+2]
        endif else xf = xf[where(xf NE j)]
        stdev_fit[j] = sqrt(total(col-DSP_fringes(u,a)^2)/dof) 
        FIT[j,*] = DSP_fringes(u,a)
     endfor
  endfor
  writefits, tag+"_all_data_DSPm_fit.fits", FIT, /compress
  



  ;FITTING OF THE WAVELENGTH FUNCTION LAMBDA = F(PIXEL)
   &  x = findgen(nl)+1  &  Siggg = 1./Sigg[xf-1]


  cut = where(L[xf-1]/Siggg LT 4. AND L[xf-1]/Siggg GT .25) ;making a cut of bad fits
  xfc = xf[cut]  &  Siggg = Siggg[cut]
  if n_elements(cut) LT 5 then begin
     print, "Too many bad points to fit"
     fxaddpar, hedr, 'STATUS', 'REJECT', ' TOO MANY BAD PTS FOR WAVELAW FIT'
     goto, the_end
  endif
  
  r = poly_fit(xfc,Siggg,3,chisq = Sig_stdev_fit)
  xr = fltarr(4,nl)+1  &  xr[1,*] = x  &  xr[2,*] = x^2  &  xr[3,*] = x^3 
  Sigfit = reform(xr##r)  &  Sigfit = 1./Sigfit


  test = total(abs(1./Sigfit-L))  &  hedr = hdr
  test = strcompress(string(test),/remove_all)
  fxaddpar, hedr, 'WLAW_CHK', strmid(test,0,5), $
            ' TOTAL DIFF. BTWN 1ST AND 2ND PASS', before = 'STATUS'

  if test GT .01*nl then begin
     fxaddpar, hedr, 'STATUS', 'REJECT', ' WAVELAW FITTED, BAD FIT'
     print, "Wavelength fit is not good enough."
  endif else begin
     stat = fxpar(hedr,'STATUS')
     fxaddpar, hedr, 'STATUS', stat, ' WAVELAW FITTED, GOOD FIT'
     print, "Good wavelength fit."
  endelse



  ;CORRECTION OF THE V-GROOVE PITCH TO GET VALUES DEPENDING ON THE PIXELS
  Pvgfit = Sigfit*L*p_VG

  if keyword_set(bin) then bit = '_bin'+string(fix(bin)) else bit = ''
  save_name = strcompress(tag+bit, /remove_all)
  ;creating the save name skeleton

  
  ;The fitted wavelengths and polynomial law
  ;this is a log for human purposes, will next save a file for program purposes
  Lfit = 1./Sigfit  &  Sig = Sigfit  &  cols_fit = xf


  filep=strcompress(save_name+"_Spatial_freq.txt", /remove_all)           
  openw, 110, filep  
  printf, 110, "p_VG (mm):";
  printf, 110, p_VG ;
  printf, 110, ""
  printf, 110, "Information from interferogram fit:"
  printf, 110, "Fitted columns: "
  printf, 110, cols_fit
  printf, 110
  printf, 110, "Wavenumbers (um^-1):"
  printf, 110, Sig ;fitted on interferograms
  printf, 110
  printf, 110, "Wavenumber Standard Deviation of fit (um^-1):"
  printf, 110, reform(param_err[0,*])
  printf, 110
  printf, 110, "Wavelengths (um):";
  printf, 110, Lfit ;fitted on interferograms  
  printf, 110
  printf, 110, "Background Flux:"
  printf, 110, Back
  printf, 110
  printf, 110, "Background Flux stdev:"
  printf, 110, reform(param_err[1,*])
  printf, 110
  printf, 110, "Uncorrelated Flux:"
  printf, 110, mu0
  printf, 110
  printf, 110, "Uncorrelated Flux stdev:"
  printf, 110, reform(param_err[2,*])
  printf, 110
  printf, 110, "Further information (Re, Re_err, fit, errors on all, etc.) saved in:"
  printf, 110, strcompress(save_name+"_Spatial_freq_extr.sav", /remove_all)
  printf, 110
  printf, 110
  printf, 110, "Information from polynomial fit: "
  printf, 110, "Wavenumbers (um^-1):"
  printf, 110, Sigfit
  printf, 110
  printf, 110, "Wavenumber stdev (um^-1):"
  printf, 110, Sig_stdev_fit ;fitted on interferograms
  printf, 110
  printf, 110, "V-groove pitch (um) based on polynomial law" ;
  printf, 110, Pvgfit ;  // new 18.06 : save the Pvg depending on pixel, from polynomial law
  close, 110 


  otro_err = param_err[3:*, *]
  save, Re_err, fit, stdev_fit, Re, otro_err, filename = strcompress(save_name+"_Spatial_freq_extr.sav", /remove_all)
  ;saving everything else we got out of the fit in case we want it later


  print, "Spatial frequencies results saved in: " ;
  print, "   "+save_name+"_Spatial_freq.txt" ;
  

  Pvg = Pvgfit
  save, p_VG, XVG, Sig, Pvg, filename = strcompress(save_name+'_Spatial_freq.sav', /remove_all)
  ;The desired information for get_best_params is saved in this file

  the_end: return, hedr
end

;NOTE TO SELF:
;On the IDL page for the FFT there is an extensive description of the
;runtime concerns for the FFT, which we could use to have our pipeline
;optimize for the number of pixels to cut to versus bad signal to fit
;    http://www.physics.nyu.edu/grierlab/idl_html_help/F4.html






function calib_lambda, data, BIN=bin, SHOW=show
;+
;
;FUNCTION: CALIB_LAMBDA(DATE, [...BIN = 'binsize']
;
;PURPOSE:
;  Returns a first pass guess of the wavelengths corresponding to the
;  spectral axis of the data.
;
;METHOD:
;
;INPUT:
;  data = a 2D collapsed version of the data cube being fit
;  /bin = if keyword is set with the size of the bin (preferably a
;         string, but the program converts whatever data type to
;         whatever it needs to avoid trouble)
;  /show = if keyword is set, the results of the chi square will be
;          displayed and the image saved.
;
;OUTPUT:
;  lambda_arr = the wavelengths corresponding to the spectral axis of the data
;
;NOTES:
;  - cent_wav_arr = findgen(201)/1000. + 0.65 ;0.65 to 0.75 
;    NEW: extended range .55 to .9 because not fitting incidence angle
;    ->changed back, .65 to .75; NEW: changed to .65 to .85 for the sake of the RHS
;  - ;trying to assess the center "half"
  ;1d array that is the average spectrum of the data. The center half with the most flux was selected.
  ;THERE'S NO GUARANTEE THAT THIS IS ALIGNED IN THE CENTER
;spectrum_theo_g--> ;all constants obtained from L212 in Elsa's 5.get_spatial_freq.i
;
;XDIM = NL, CORRECT AND REPLACE
;  ;is continuum detection essentially msmooth?Y
;getting a "chi-square" based on a continuum normalized version of the spectra
;
;avoiding NaN, don't ask me how they occasionally get in there
;
;-
  common constants
  common CL_constants, lambda_A, Flux0, dirr

  if keyword_set(dir) then begin
     dir = dir 
     dirr = strmid(dir, 0, strlen(dir) -strpos(dir, 'data')-2)+'reference_spectra/'
  endif else begin
     dir = 'data/'  &  dirr = 'reference_spectra/'
  endelse




  ;Getting the center "half" of the image...
  vc = fxpar(hdr,'VERTCEN')  
  if part EQ "right" then ind = [5,3] else ind = [0,3]
  vc = strmid(vc,ind[0],ind[1])  &  wid = np/4. 
  if vc-wid LT 0 then bs = 0 else bs = vc-wid
  if vc+wid GT np-1 then ts = np-1 else ts = vc+wid
  spectr = avg(data[*,bs:ts],1) 

 


  ;Loading in the spectrum info to fit...
  fname = strcompress(dirr+'Pickles_'+star+'.txt',/remove_all)
  readcol, fname, lambda_A, flux0, /silent




  ;Calculating the constants for the fit...
  if keyword_set(bin) then nbin = double(bin) else nbin = 1D
  nl = n_elements(spectr)  &  xarr = findgen(nl)*nbin
  fit = poly_fit(xarr, spectr, 5, yfit = fit_5th)
  spectrum_theo_g = spectrum_theo(.643, 68*!dtor, 1.e-6, fit_5th, bin=bin)
  T = min(spectrum_theo_g)/min(msmooth(spectrum_theo_g, nl*.05))*1.e-6 



  ;Setting up the arrays for the Chi-square fit...
  if part EQ 'right' then L0_arr = findgen(121)/1000. + 0.68 $
  else L0_arr = findgen(101)/1000. + 0.65
  I0_arr = (findgen(51)/100.*20 + 50)*!dtor
  x2 = fltarr(n_elements(L0_arr), n_elements(I0_arr))
   

  ;Using "Chi-square" minimization to fit for L0 and I0...
  print, "Computing the X^2 map."
  print


  for i = 0, n_elements(L0_arr)-1 do begin
     for j = 0, n_elements(I0_arr)-1 do begin
        spectrum = spectrum_theo(L0_arr[i], I0_arr[j], T, fit_5th, bin=bin)
        

        msr = msmooth(spectr,10)  &  dat = (spectr-msr)/min(msr)
        msm = msmooth(spectrum,10)  &  theor = (spectrum-msm)/min(msm)
 

        x2[i,j] = total((dat - theor)^2.)         
    endfor
  endfor


  ind = where(x2 EQ x2)  &  if total(ind) EQ -1 then goto, fail
  loc = max(where(x2 EQ min(x2[ind])))  &  loc = array_indices(x2, loc)

  cent_wav = L0_arr[loc[0]]     ;best guess of the central wavelength
  inc_angl = I0_arr[loc[1]]     ;best guess of the incidence angle




  ;Showing the results of this fit...
  print, "The best guess of the central wavelength is: ", cent_wav
  print, "The best guess of the incidence angle is: ", inc_angl


  if keyword_set(show) then begin
     loadct, 19
     display, x2, min = min(x2)*1.1, max = max(x2)*.9, xrange=[0.65, 0.75], yrange=[50, 70], background=!white, color=0, $
              title=textoidl('\chi^2')+' Map', $
              xtitle="Central wavelength ["+textoidl('\mu')+"m]", $
              ytitle="Incidence angle on the prism (deg)", charsize=2
     loadct, 0
     xyouts, cent_wav, inc_angl, "O", color = 0
  endif   
  



  ;Saving the results of the fit...
  lambda_arr = wavelaw(cent_wav, inc_angl, bin = bin)  &  L = lambda_arr
  if keyword_set(bin) then bit = strcompress('_bin'+string(fix(bin)),/remove_all) else bit = ''
  save, L, cent_wav, inc_angl, filename = tag+'_CLfit2'+bit+'.sav'

  return, lambda_arr
  fail: return, 0
end




function spectrum_theo, L0, I0, T, fit_5th, BIN=bin
;+
;
;FUNCTION: SPECTRUM_THEO(L0, I0, T, FIT_5TH, [...BIN='binsize'])
;
;PURPOSE:
;  Model function for the spectrum, comprising stellar spectrum and
;  atmosphere transmission.
;
;METHOD:
;
;INPUT:
;  L0 = central wavelength
;  I0 = incidence angle of the light
;  T = ratio of the deepest stellar feature to the continuum
;  fit_5th = the shape of the continuum
;  /bin = if keyword is set, columns have been binned by the given
;         binsize (should be a string, but code is robust to anything)
;
;OUTPUT:
;
;NOTES:
;  -As they are read in, lambda_A is in Amperes and flux is in
;   erg/(s*cm^2*A). They are converted to µm and W/(µm*m^2).
;  -This is a change from Elsa's code that we believe will be more 
;   effective in flattening out the spectrum (L606)
;  -interpolate to get the flux at each pixel
;  -returns the percentage transmission of each wavelength
;  -make flux curvy like the observed spectrum
;
;-
  common constants
  common CL_constants
  h = 6.626068e-34  &  c = 299792458.




  ;Ensuring that spectral info is available...
  if n_elements(lambda_A) EQ 0 OR n_elements(flux0) EQ 0 then $
     readcol, dirr+'Pickles_'+star+'.txt', lambda_A, flux0, /silent


  ;Changing units...
  Flux = Flux0 * 10. 
  lambda_um = lambda_A * 1.e-4 


  ;Creating an array of the delta lambdas...
  ldif_um=lambda_um*0.  
  for i = 0, n_elements(lambda_um)-2 do begin
     ldif_um[i] = abs(lambda_um[i+1] - lambda_um[i])
  endfor
  ldif_um[-1] = ldif_um[-2]
  

  ;Obtaining corresponding wavelengths...
  Flux = Flux*lambda_um*ldif_um/(h*c) ;photons per m^2
  L = wavelaw(L0, I0, bin=bin) 


  ;Getting the atmospheric transmission...
  readcol, dirr + 'Tatm_short.txt', lambda, transmission, /silent  
  t_atm = interpol(transmission, lambda, L)
 
  ;Fitting the global shape of the spectrum...
  xarr = findgen(n_elements(lambda_um))
  coeff = poly_fit(xarr, flux, 3, yfit=fity)  
  Flux = Flux/fity*min(fity)  &  Flux = interpol(Flux, lambda_um, L) 
  Flux *= tau * Surf * fit_5th/max(fit_5th) * T * t_atm

  return, Flux
end




function wavelaw, L0, I0, BIN=bin
;+
;
;FUNCTION: WAVELAW(L0,IO,[...BIN='binsize'])
;
;PURPOSE:
;  Generates the wavelength values corresponding to the pixels
;
;INPUT:
;  L0 = central wavelength ()
;  I0 = incidence angle of the light (rad)
;  /bin = if keyword is set, columns have been binned by the given
;         binsize (should be a string, but code is robust to anything)
;
;OUTPUT:
;
;NOTES:
;  -Thorlabs prisms used during the 2011 runs are in F2 and not SF2
;   Returns the refractive index of the F2 glass for wavelengths L (in um)
;   It uses Sellmeier formula for F2 glass.
;   Source: http://refractiveindex.info/?group=SCHOTT&material=F2
;  -In determining the deviation as a function of lambda, you may notice
;   that there is a missing !pi to be added in. This term, however, does
;   not affect the tangent and therefore has been left out.
;  -Usage of nbin (vs bin) is to avoid overwriting variable values in
;   other functions. Use of nbin in general is to account for
;   "shortening" of the horizontal axis due to binning.
 ;0.4 to 1
; no sorting step at last line because of how interpol works
;-
  
  common constants
  common CL_constants


  ;Finding the index of refraction as a function of wavelength...
  ld = findgen(1001)/1000.*0.6 + 0.4
  case 1 of
     glass EQ 'Si': n = sqrt( 1 $
                              + 0.6961663*ld^2/(ld^2-0.0684043^2) $
                              + 0.4079426*ld^2/(ld^2-0.1162414^2) $
                              + 0.8974794*ld^2/(ld^2-9.896161^2) )
     glass EQ 'BK7': n = sqrt( 1 $
                               + 1.03961212*ld^2/(ld^2-0.00600069867) $
                               + 0.231792344*ld^2/(ld^2-0.0200179144) $
                               + 1.01046945*ld^2/(ld^2-103.560653) )
     glass EQ 'SF2': n = sqrt( 1 $
                               + 1.40301821*ld^2/(ld^2-0.0105795466) $
                               + 0.231767504*ld^2/(ld^2-0.0493226978) $
                               + 0.939056586*ld^2/(ld^2-112.405955) )
     glass EQ 'F2': n = sqrt( 1 $
                              + 1.34533359*ld^2/(ld^2-0.00997743871) $
                              + 0.209073176*ld^2/(ld^2-0.0470450767) $
                              + 0.937357162*ld^2/(ld^2-111.886764) ) 
     else: stop, "Glass type given is not recognized."
  endcase



  
  ;Finding the deviation as a function of lambda
  A=60.*!dtor ;prism angle (rad)
  D = I0+asin(n*sin(A-asin(sin(I0)/n)))-A 

 
  ;Ensuring that the null deviation is centered on L0...
  ind = min(where(abs(ld-L0) EQ min(abs(ld-L0))))
  if ind NE -1 then begin
     D0 = D[ind]  &  D -= D0
  endif else stop, "Wat."




  ;Finding position on the CCD as a function of lambda...
  if bin then nbin = double(bin) else nbin = double(1) 
  x_px = foc*tan(D)/p_det
  X = findgen(nl)*nbin-fix(nl*nbin)/2
  L = interpol(ld, x_px, X)  &  L = L[sort(L)]


  return, L
end























