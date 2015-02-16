function fit_binary, file, num, S=S, B=b, UNQ_TAG=unq_tag, LOG=log, $
                     WAVE_LIM=wave_lim, $
                     RIGHT=right, BIN=bin, DISPLAY=display, SAVE=save
;+
;
;FUNCTION: FIT_BINARY(FILE, [.../S, /DISPLAY, /SAVE]
;
;PURPOSE:
;  To fit the separation and position angle of the binary as observed by
;  FIRST based on the closure phases.
;
;INPUT:
;  file = If init = the four character indentifier string for a given
;         file, then file = A string of the form: 
;         (init for target)_(init for calibrator)_(date)
;         Example: file would be the two number sequences beginning
;                   "7129_7123_20111017_S1_b100_(stuff)" 
;  num = The string of a number associated with the data cubes that were
;        used in DO_CALIB to perform the calibration. May be ''.
;  /S = The threshold used in GET_BEST_PARAMS. If keyword is not set,
;       then the default value of S = 1 will be used.
;  /b = The bundle size used in GET_BEST_PARAMS. If keyword is not set,
;       then the default value of b = 100 will be used.
;  /unq_tag = A unique tag added on to the filename (stuff in the
;             example above).
;  /log = If keyword is set, then the values of rlim used will be
;         assumed to be exponents of powers of 10. 
;  /wave_lim = Wavelength limits to be covered in the fit.
;  /right = If keyword is set, right side of image is considered. 
;  /bin = If keyword is set, then bin is equal to the size by which the
;         spectral axis of the image was binned.
;  /display = If keyword is set then program will plot the results of
;             the binary fit. 
;  /save = If keyword is set then 
;
;OUTPUT:
;
;NOTES:
;   ;don't need lg limits because already synced in do_calib
  ;l1g = lglim[0]  &  l2g = lglim[1]
  ;l1d = ldlim[0]  &  l2d = ldlim[1]
;-

  ;common year, lglim, ldlim, xlim, ylim, rlim, Npos, Nr, Nlc, Nl_pts
  ;lglim = the left hand limits on the wavelaw, ldlim = same on the right
  ;have to obtain these from somewhere








  ;Defining the limits of the fit...
  xlim = fltarr(2)
  xlim[0] = .010
  xlim[1] = .070
  rlim = [0.8,1.25]
  ;xlim[0] = .041                ;mas, diffraction limit at 600 micron
  ;xlim[1] = .400 ;mas, FOV of telescope at 800 micron
  ylim = xlim  &  
  Npos = (xlim[1]-xlim[0])/.001
;rlim = [.2,3.]  &  
  Nr = (rlim[1]-rlim[0])/.01
  if keyword_set(bin) then Nl_pts = bin else Nl_pts = 1.
  ;Elsa has externalized values for the wavelength limits, I made them options



  ;Loading in the calibrated closure phases...
  if keyword_set(right) then part = "right" else part = "left"  
  year = strmid(file,10+5,4)
  if keyword_set(S) then S = 'S'+string(S) else S = 'S1'
  if keyword_set(b) then b = 'b'+string(b) else b = 'b100'
  if keyword_set(unq_tag) then u = unq_tag+'_' else u = ''

  cp_name = file+"_"+S+"_"+b+"_"+u+part
  restore, cp_name+"_CPm"+num+".sav" 
  if keyword_set(wave_lim) then begin
      inds = where(L GE wave_lim[0] AND L LE wave_lim[1])
      L = L[inds]  &  CPm = CPm[inds,*,*]
   endif else inds = where(L GT 0.65)
  L = L[inds]  &  CPm = CPm[inds,*,*]

  nCP = (size(CPm))[2]  &  nl = n_elements(L)  
  if keyword_set(Nl_pts) then Nlc = 0 else Nlc = nl
  nF = 9  &  nB = fix(nf*(nf-1)/2.)
  
  L = quik_bin(L,3.)  &  CPt = fltarr(nl/3., nCP,2)
  for i = 0, nCP-1 do begin
     CPt[*,i,0] = quik_bin(CPm[*,i,0],3.)
     CPt[*,i,1] = quik_bin(CPm[*,i,1],3.)
  endfor
  CPm = CPt  &  nl /=3.  &  Nlc /= 3.
  
  ;If keyword is set to display, then displaying closure phases...
  if keyword_set(display) then begin
     window, 0
     display, reform(CPm[*,*,0]), $
              xtitle = "Spectral Channel (pixel)", $
              ytitle = "Baseline Triangle", $
              title = "Closure Phases: "+date2+", "+star+" ("+part+")"

     loadct, 39
     oplot, [0,nCP], l1*[1,1]-1, linestyle = 1, color = 250
     oplot, [0,nCP],n_elements(L)*[1,1]+l1, linestyle = 1, color = 100
     loadct, 0
  endif


  save_name = cp_name+"_bifit.sav"
  print,"Saving under: "+save_name
  print, "Npos = "+string(Npos)+",  Nr = "+string(Nr)




  ;Setting up the X, Y and r grids...
  Xp = findgen(Npos)/(Npos-1)*(xlim[1]-xlim[0])+xlim[0]
  Yp = findgen(Npos)/(Npos-1)*(ylim[1]-ylim[0])+ylim[0]
  if keyword_set(log) then begin
     r = 10^(findgen(Nr)/(Nr-1)*(rlim[1]-rlim[0])+rlim[0])
  endif else begin
     r = findgen(Nr)/(Nr-1)*(rlim[1]-rlim[0])+rlim[0]
  endelse 






  ;Wavelength binning...
  if Nl_pts NE 1 then begin
     nlb = fix(floor(float(nl)/float(Nl_pts)))
     nres = nl-nlb*Nl_pts

     Lb = L[0:nlb-1]*0
     CPmb = CPm[0:nlb-1,*,*]*0
     ind = findgen(Nl_pts)*Nl_pts
     for k = 0, Nl_pts-1 do begin
        inds = ind + k
        Lb += L[inds]
        CPmb[*,*,0] += CPm[inds,*,0]
        CPmb[*,*,1] += 1./CPm[inds,*,1]^2
     endfor

     Lb /= Nl_pts  &  CPmb[*,*,0] /= Nl_pts
     CPmb[*,*,1] = 1./sqrt(CPmb[*,*,1])

                                ;Elsa keeps the last point being
                                ;binned by binning it with
                                ;whatever's left...I chose to cut
                                ;it off instead

     nl = nlb  &  L = Lb  &  CPm = CPmb
  endif



  ;Getting the baselines in the pupil plane...
  UV2D = get_UV2d_by_year(year, L, Dpup=Dpup, right=keyword_set(right))
  Dtel = 7*Dpup
  t = findgen(1000)/999.*2*!pi  &  xt = cos(t)  &  yt = sin(t)




  ;Initializing...
  chi2map = dblarr(Npos,Npos,Nr,Nl)
  A = fltarr(3,Nl)  &  afit = dblarr(Npos,Npos,Nl)
  CPmod = reform(CPm[*,*,0]*0)
  arr = fltarr(3)
  if n_elements(selec) EQ 0 then selec = findgen(nCP)

  for i = 0, Npos-1 do begin
     arr[1] = Xp[i]  
     
     for j = 0, Npos-1 do begin
        arr[2] = Yp[j]
        
        for k = 0, Nr-1 do begin
           arr[0] = r[k]  &  CPsim = CP_binary(UV2D,arr)

           D = abs(transpose(CPm[*,selec,0])-CPsim[selec,*])
           ;plot, L, CPm[*,70,0], psym=2,/xstyle
           ;oplot, L, CPsim[70,*], color = (k+1)/(1.*nr)*150+100
           ;wait, .001
           D += (D GT 180.)*(360. - 2.*D)

           batch = (D/transpose(CPm[*,selec,1]))^2
           cut = fltarr(nCP*0.8,nl)
           for lm = 0,nl-1 do begin
              inds = sort(batch[*,lm])
              cut[*,lm] = batch[inds[0:nCP*0.8-1]]
           endfor

           Chi2map[i,j,k,*] = reform(total(cut,1))
        endfor
     endfor
  endfor
  
       
  if keyword_set(display) then begin
     print, "*** RESULTS ***[ "+string(l1)+" - "+string(l2)+" ]************" 
     print, "Chi2 min   :\t"+string(min(Chi2))                                   
     print, "Chi2 median:\t"+string(median(Chi2[*]))                                 
     print, "\n Optimal flux ratio:\t"+string(a[1])                                 
     print, "Optimal X position:\t"+string(a[2])                                     
     print, "Optimal Y position:\t"+string(a[3])                                     
     print, "Separation        :\t"+string(sqrt(a[2]^2+a[3]^2))                      
     print, "*********************************************"    
  endif
  

  Chi2 = total(reform(min(Chi2map,dimension = 3)), 3)


  dummy = min(min(Chi2, dimension = 2),indX)
  dummy = min(min(Chi2,dimension=1),indY)
  print, "x= "+string(Xp[indX])+" y="+string(Yp[indY])+" sep="+string(sqrt(Xp[indX]^2+Yp[indY]^2))
  dummy = min(Chi2map[indX,indY,*,*],mnx, dimension=3)
  Rfit = r[mnx]
  print, "Avg Flux ratio = ", avg(Rfit)


  print, "CP reconstruction..."
  for c = 0, nl-1 do begin
     ;l1 = l[c]  &  l2 = l[c+1]
     ind = c;where(L LT l2 AND L GE l1)
     a = [Rfit[c],Xp[indX],Yp[indY]]
     ;if Arr[4,c] NE 0 then 
     CPmod[ind,*] = CP_binary(UV2D[*,*,ind],a)
  endfor

  print, 'Saving...'
  ;if keyword_set(save) then begin


  save, filename=save_name, Chi2map, Xp, Yp, r, L, A, Dtel, CPmod 
  print, "SAVE NAME= "+save_name                                              
  
  print, "...DONE!"                                 
  return, 1

end


function quik_bin, val, box

  new_val = fltarr(n_elements(val)/box)
  for i = 0, n_elements(val)/box-1 do begin
     new_val[i] = avg(val[i*box:(i+1)*box-1])
  endfor
  return, new_val
end


function get_UV2D_by_year, year, L, DPUP=dpup, RIGHT=right, SEGMENT=segment
;+
;
;FUNCTION: GET_UV2D_BY_YEAR(YEAR,L,[.../RIGHT, SEGMENT=SEGMENT]
;
;PURPOSE:
;  To get the baselines of the pupil plane as a function of the year and
;  wavelaw. 
;
;INPUT:
;  year = The year in which the observations were taken.
;  L = The wavelaw associated with the calibrated closure phases.
;  /Dpup = The ratio of the pupils to the telescope diameter.
;  /right = If keyword is set, the right hand side of the v-groove will
;           be considered.
;  /segment = If the segment is already known 
;
;
;OUTPUT:
;  uv2d = The baselines in the pupil plane.
;
;NOTES: (NONE)
;
;-
  ;Values specific to FIRST...
  sides = ['left','right'] 
  if keyword_set(right) then part = sides[1] else part = sides[0]

  seg_info = {yr2011:[[19,23,16,14,10,32,27,17,36],$
                                [26,34,31,30,13,20,21,24,37]],$
              yr2012:[[21,18,22,8,35,17,36,31,26],$
                                [9,13,10,20,11,29,24,33,28]],$
              yr2013:[[21,18,22,8,35,17,36,31,26],$
                               [9,13,10,20,11,29,24,33,28]]}
  years = [2011,2012,2013]  &  dpup = [8.2/7, 3/7., .365]
  



  ;Matching these values to input...
  yind = where(years EQ year)  &  pind = where(sides EQ part)
  if ~keyword_set(segment) then segment = seg_info.(yind)  
  segment = segment[*,pind]  &  dpup = dpup[yind]




  ;Finding the 2 dimensional baselines...
  nf = n_elements(segment) & nB = round(nf*(nf-1)/2.)
  xy = get_fibpos(segment)
  uv = transpose(get_uv2d(xy[0,*], xy[1,*]))


  c_rad_to_mas = !radeg*3600.*1000.


  nl = n_elements(L)
  uv2d = fltarr(nB,2,nl)
  for i = 0,nl-1 do uv2d[*,*,i] = uv*(dpup/(L[i]*1e-6)/c_rad_to_mas)[0]
  
  return, uv2d
end








function get_fibpos, segment, INTEGER=integer
;+
;
;FUNCTION: GET_FIBPOS(SEGMENT, [.../INTEGER])
;
;PURPOSE:
;  This function will return the coordinates of the segments in the
;  pupil plane.
;
;INPUT:
;  segment = The subpupil numbers for the segments being considered. 
;  /integer= If keyword is set, function will return the coordinates in
;            integer numbers, otherwise function will return the
;            coordinate in their real proportions.
;
;OUTPUT:
;  result = The coordinates of the segments in the pupil plane.
;
;NOTES: (NONE)
;
;-
  
  x0 = [6,8,7,5,4,5,7,10,9,8,6,4,3,2,3,4,6,8,9,12,11,10,9,7,5,3,2,1,0,1,$
        2,3,5,7,9,10,11]/2.  
  if keyword_set(integer) then begin
     y0 = [3,3,4,4,3,2,2,3,4,5,5,5,4,3,2,1,1,1,2,3,4,5,6,6,6,6,5,4,3,2,$
           1,0,0,0,0,1,2]
  endif else begin
     y0 = [3,3,4,4,3,2,2,3,4,5,5,5,4,3,2,1,1,1,2,3,4,5,6,6,6,6,5,4,3,2,$
           1,0,0,0,0,1,2]*sin(!pi/3.)
  endelse


  x = x0[segment]  &  y = y0[segment]
  result = transpose([[x],[y]])

  return, result
end








function get_UV2D, x, y, POSITIVE=positive
;+
;
;FUNCTION: GET_UV2D(X,Y,[.../POSITIVE])
;
;PURPOSE:
;  This function will return the baselines corresponding to the given
;  coordinates in the pupil plane.
;
;INPUT:
;  x,y = coordinates of the subpupils in the pupil plane
;  /pos = if you want all the coordinates to have a positive x coordinate
;
;OUTPUT:
;  uv = The baselines corresponding to the coodinates in the pupil plane.
;
;NOTES:
;
;-

  nx = n_elements(x)  &  uv = []


  for i = 0, nx-2 do begin
     for j = i+1, nx-1 do begin
        if keyword_set(positive) then begin
           case 1 of
              x[j] GT x[i]:  uv = [[uv],[abs(x[j]-x[i]), y[j]-y[i]]]
              x[j] LT x[i]: uv = [[uv],[abs(x[j]-x[i]), y[i]-y[j]]]
              else: uv = [[uv],[abs(x[j]-x[i]), abs(y[i]-y[j])]]
           endcase
        endif else begin
           uv = [[uv],[x[j]-x[i],[y[j]-y[i]]]]
        endelse
     endfor
  endfor


  return, uv
end








function cp_binary, uv, a
;+
;
;FUNCTION: CP_BINARY(UV,A)
;
;INPUT:
;  uv = An array of 3 dimensions such that:
;       1 = frequency
;       2 = x / y
;       3 = wavelength
;  a = A 3-element vector such that:
;      a(1) = flux ratio                                                                          
;      a(2) = x position of the companion (central star being at (0,0)) 
;      a(3) = y position of the companion  
;
;OUTPUT:
;  CP = The closure phases in degrees corresponding to a binary
;       system with parameters uv and a.
;
;NOTES: (NONE)                   
;
;-

  nF = 9. ;number of fibers
  r = abs(a[0])  &  xp = a[1]  &  yp = a[2]

  Vc = reform(1./(1.+r)*(1.+r*exp(-complex(0,2)*!pi*1000.*(xp*UV[*,0,*]+yp*UV[*,1,*]))))

  indCLO = get_CLO(nF)                                       
  Bi = Vc[indCLO[0,*],*]*Vc[indCLO[1,*],*]*conj(Vc[indCLO[2,*],*]) 
  CP = atan(Bi,/phase)*!radeg                           

  return, CP
end  
