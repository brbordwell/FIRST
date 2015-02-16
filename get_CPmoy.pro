pro get_CPmoy, file, S, b, DIR=dir, RIGHT=right, BIN=bin, SHOW=show
;+
;
;FUNCTION: GET_CPMOY, FILE, S, B, [...DIR='path', BIN=#, /RIGHT, /SHOW]
;
;PURPOSE:
;
;INPUT:
;  file = The file for which the closure phases will be calculated.
;  S = The threshold:
;     -If S is LE 1, then that fraction of measurements will be used. 
;     -If S is GT 1, then bispectrum moduli above S will be kept.
;     -If S is EQ 0, then the program will default to using the average
;      measurement from a stack averaged through the stack of images in
;      a column half way through the spectral axis. 
;  b = The bundle size
;  /dir = The path to where files are to be loaded from, and saved.
;  /right = If keyword is set, the right side of the image is used.
;  /bin = The size of the bins into which the columns are combined.
;  /show = If keyword is set, the results will be plotted. 
;
;OUTPUT:
;  CP = array of dimensions: nl by nCP by nbmax by 4
;       4: -1 for the mean values per bundle
;          -2 for the error per bundle (simple rms and rms weighted with
;           bispectrum modulus)
;          -1 for the position among the images  
;  CPmoy1 = Phase of averaged bispectrum, array of dimensions: 
;           nl by nCP by 2,  
;           2: -1 for mean value of complex part of bispectrum
;              -1 for error (rms/sqrt(N))
;  CPmoy2 = array of dimensions: nl by nCP by 2
;           2: -1 mean value estimated from image bundles
;              -1 for error (rms/sqrt(N)) corrected by X2 calculation
;
;NOTES:
;  -nl = dimensions of spectral axis (in pixels)
;  -nCP = number of triangles used for closure phases
;  -nbmax = maximum number of bundles
;
;- 

  if bin EQ 0 or bin EQ "" then nbin = "" $
  else nbin = strcompress('_bin'+string(fix(nbin)),/remove_all)
  if ~keyword_set(dir) then dir='data/'
  if keyword_set(right)  then part = 'right' else part = 'left' 



  
  ;Loading bispectrum...
  print, "Loading bispectrum..."

  bifile = strcompress(dir+file+"_"+part+nbin+"_Bi.sav", /remove_all)
  restore, bifile ;Restoring Bi and XVG.


  nF = n_elements(XVG[*,0])     
  Bvg = get_UV(XVG)  &    nB = n_elements(Bvg)          
  indCLO = get_CLO(nf)  &  nCP = n_elements(indCLO[0,*]) 


  Barr = temporary(Bi)  &  sizB = size(Barr)  
  nl = sizB[1]  &  nCP = sizB[2]  &  ni = sizB[3]




  ;Establishing threshold...
  Bmod = abs(Barr)
  if S EQ 0 then begin
     S = avg(avg(Bmod[round(nl/2.),*,*], 2)) 
     print, "Default Threshold", string(S)
  endif


  ;Getting user input if the threshold wasn't specified...
  if keyword_set(show) then begin
      window, 1
      plot, avg(Bmod[round(nl/2.), *, *],2), psym = 2, $
            title = "Bispectrum Modulus", xtitles = "Baseline Triangle"
      oplot, [0, nCP+1], [1.,1.]*S, psym = 4, color = !red


      S = 0  &  ok = 0
      while ok NE 1 do begin
          print, "Threshold can be chosen to be either GT 1, "
          print, "and only measurements with bispectrum modulus"
          print, "greater than S (threshold) will be kept,"
          print, "or threshold can be chosen to be LE 1, "
          print, "and S will then correspond to"
          print, "the fraction of the amount of measurements"
          print, "that will be taken into account."
          
          read, "Threshold (S)? ", S
          if S LT 1 AND S GT 0 then begin
             print, "Selecting " + $
                    strcompress(string(S*100), /remove_all) + $
                    "% of the data?"
          endif else begin
             print, "Threshold at " + $
                    strcompress(string(S),/remove_all) + "?"
          endelse

          ok = 0  &  read, "yes (1), no (0) ", ok
      endwhile
  endif

  
  ;Getting the bundles worked in...
  if S LT 1 AND S GT 0 then begin
      print, "Selected " + string(S*100)+"% of the data"
      nbmax = fix(ni/double(b)*double(S))
  endif else begin
      print, "Threshold = "+ string(S)
      nbmax = fix(ni/double(b))
  endelse
  
  print, "Bundles of b = " + strcompress(string(b),/remove_all) + $
         " meas. -> nbmax = " + strcompress(string(nbmax),/remove_all)
  
  
  
  
  ;Capturing the closure phases...
  CParr = dblarr(nl, nCP, nbmax+1, 4)
  CPmoy1 = dblarr(nl, nCP, 2)  &  CPmoy2 = dblarr(nl, nCP, 4)
  SNR = dblarr(nl, nCP)
  nis = round((1-S)*ni)


  for l = 0, nl-1 do begin
      if l mod 10 EQ 0 then print, strcompress("l=   "+string(l)+"/"+ $
                                               string(nl),/remove_all)

      for k = 0, nCP-1 do begin
          if S LE 1 AND S GT 0 then ind = (sort(Bmod[l,k,*]))[nis:*] $
          else ind = where(Bmod[l,k,*] GT S)

    
          ;Summing the bundles...
          if total(ind) NE -1 then begin
             n = n_elements(ind)

             bbb = Barr[l,k,ind]
             
             i = imaginary(avg(bbb))  &  r = real_part(avg(bbb))
             CPmoy1[l,k,0] = atan(i,r)*!radeg
             
             i = imaginary(bbb)  &  r = real_part(bbb)
             cp = atan(i, r)*!radeg - CPmoy1[l,k,0]
             
             if n EQ 1 then CPmoy1[l,k,1]  = 0 $
             else CPmoy1[l,k,1] = stdev(cp)




             if n GE 2 then begin
                ;Added dummy variables to make the code more readable...
                i = imaginary(bbb)  &  r = real_part(bbb)
                cov = avg(i*r)-avg(i)*avg(r)
                ang = CPmoy1[l,k,0]*!dtor  
                sr = stdev(r)^2  &  si = stdev(i)^2

                btm = (si*cos(ang)^2+sr*sin(ang)^2-sin(2*ang)*cov)
                SNR[l,k] = sqrt(n)*abs(avg(bbb))/sqrt(btm)
             endif else SNR[l,k] = 0



                         
             if n GE b then begin
                n = n_elements(ind)  &  npos = floor(n/double(b))
                CParr[l,k,0,*] = npos 

                for j = 0, npos -1 do begin
                   ;Finding the phase of Bi avg'd over b pts...
                   pos = ind[j*b :(j+1)*b-1]  &  P = round(median(pos))
                   CParr[l,k,j+1,2] = P
                   
                   bbb = Barr[l,k,pos]
                   i = imaginary(avg(bbb))  &  r = real_part(avg(bbb))
                   cpm = atan(i, r)
                   CParr[l,k,j+1,0] = cpm*!radeg

                   i = imaginary(bbb)  &  r = real_part(bbb)
                   CParr[l,k,j+1,1] = stdev(atan(i,r)-cpm)*!radeg/sqrt(b)

                   cov = avg(i*r)-avg(i)*avg(r)
                   si = stdev(i)^2  &  sr = stdev(r)^2 
                   btm = sin(cpm)^2*si+cos(cpm)^2*sr - sin(2*cpm)*cov
                   CParr[l,k,j+1,3] = sqrt(b)*abs(avg(bbb))/sqrt(btm)
                       



                   if npos GT 1 then begin
                      sigma = reform(CParr[l,k,1:npos, 1])
                      m = reform(CParr[l,k,1:npos,0])
                      mW = warp_around_zero(m)
                      m0 = mW[0]  &  mW = mW[1:*] 
                      r = mean_estimator(mW, sigma)
                      
                      CPmoy2[l,k,0] = r[0]+m0
                      CPmoy2[l,k,1] = r[1] 
                      if r[2] GT 1 then CParr[l,k,*,1] *= r[2]

                      sigma = atan(1./CParr[l,k,1:npos,3])*!radeg
                      r = mean_estimator(mW, sigma)
                      CPmoy2[l,k,2] = r[0]+m0
                      CPmoy2[l,k,3] = r[1]
                   endif else begin
                      CPmoy2[l,k,0] = CParr[l,k,1,0]
                      CPmoy2[l,k,1] = CParr[l,k,1,2]
                      CPmoy2[l,k,2] = CParr[l,k,1,0]
                      CPmoy2[l,k,3] = CParr[l,k,1,2]
                   endelse
                endfor
             endif
          endif
       endfor
   endfor




  ;Logging in everything...
  file1 = strcompress(dir+file+'_'+part+nbin+"_S"+string(S)+ $
                      "_b"+string(b)+"_CP_log.txt", /remove_all)
  openw, 110, file1
  printf, 110, "S = "+string(S)
  printf, 110, "If S < 1, S represents the percentage of data that"
  printf, 110, "are kept. Bispectrum modulus is the selection parameter"
  printf, 110, "for each closure phase."
  printf, 110, " b = " + string(b)+ "  # of meas. averaged per point"
  printf, 110
  printf, 110, "CP: "
  printf, 110, "dim1: nl, spectral channels;"
  printf, 110, "      1. Mean values per bundle"
  printf, 110, "dim2: nCP, closure phase numbers; "
  printf, 110, "      2. Error per bundle (simple rms)"
  printf, 110, "dim3: nbmax, images; "
  printf, 110, "      3. Error per bundle (rms weighted by Bi)"
  printf, 110, "dim4: 4; "
  printf, 110, "      4. Position among the images"
  printf, 110
  printf, 110, "Stored variables: CParr, S, b"
  close, 110
  
  save, CParr, S, b, filename = repstr(file1,'_log.txt','.sav')




  file2 = strcompress(dir+file+'_'+part+nbin+"_S"+string(S) + $
                      "_CPmoy1_log.txt", /remove_all)
  openw, 110, file2
  printf, 110, "S = "+string(S)
  printf, 110, "If S < 1, S represents the percentage of data "
  printf, 110, "that are kept. Bispectrum modulus is the selection"
  printf, 110, "parameter for each closure phase"
  printf, 110
  printf, 110, "CPmoy1: "
  printf, 110, "dim1: nl, spectral channels;"
  printf, 110, "      1. Mean value (arg(avg(Bi))"
  printf, 110, "dim2: nCP, closure phase numbers;" 
  printf, 110, "      2. error (rms/sqrt(N))"
  printf, 110, "dim3: 2"
  printf, 110
  printf, 110, "Stored variables: CPmoy1, S, SNR"
  close, 110
  
  save, CPmoy1, S, SNR, filename = repstr(file2,'_log.txt','.sav')



  
  file3 = strcompress(dir+file+'_'+part+nbin+"_S"+string(S) + $
                      "_CPmoy2_log.txt", /remove_all)
  openw, 110, file3
  printf, 110, "S = "+string(S)
  printf, 110, "If S < 1, S represents the percentage of data"
  printf, 110, "that are kept. Bispectrum modulus is the selection"
  printf, 110, "parameter for each closure phase"
  printf, 110
  printf, 110, "CPmoy2: "
  printf, 110, "dim1: nl, spectral channels; "
  printf, 110, "      1. Mean estimated from the image bundles"
  printf, 110, "dim2: nCP, closure phase numbers;"
  printf, 110, "      2. Error corrected by X2 calculation(rms/sqrt(N))"
  printf, 110, "dim3: 2"
  printf, 110
  printf, 110, "Stored variables: CPmoy2, S, SNR, b"
  close, 110
  
  save, CPmoy2, S, SNR, b, filename = repstr(file3,'_log.txt','.sav')
end








function get_UV, XVG
;+
;
;FUNCTION: GET_UV(XVG)
;
;PURPOSE:
;  Finds the lengths of all the baselines
;
;INPUT:
;  XVG = The positions along the v-groove.
;
;OUTPUT:
;  UV = An array of all the baseline lengths.
;
;NOTES:
;  -The order of the baseline lengths is important, leave it.
;  -In reality there will be r/2 baselines.
;
;-

  nx = n_elements(XVG)
  r = (nx-1)*nx  &  uv = dblarr(r)
 

  for i = 0, nx-2 do begin
     for j = i+1, nx-1 do begin
        UV[j+i*(nx-1)] = abs(XVG[i] - XVG[j])
     endfor
  endfor


  UV = UV[where(UV NE 0)]
  return, uv
end








function get_ijk, cp, indCLO, XVGi, Bvgi
;+
;
;FUNCTION: GET_IJK(CP, INDCLO, XVG, BVG)
;
;PURPOSE:
;  Gets indices of the fibers forming the baseline triangle number cp.
;
;INPUT:
;  cp = Index of the baseline triangle. 
;  indCLO = indices of the baselines forming the baseline triangles.  
;  XVG = Rounded positions of the fibers.
;  Bvg = rounded value of the baselines formed by the fibers.
;
;OUTPUT:
;  IJK = The indices of the triangle used to find the closure phase.
;
;NOTES:
;  -The order of the baseline lengths is important, leave it.
;  -In reality there will be r/2 baselines.
;
;-


  indCLO = fix(abs(indCLO))

  ijk = fltarr(3)
  ijk[0:1] = get_ij(Bvgi[indCLO[0,cp]], XVGi)
  int = get_ij(Bvgi[indCLO[1, cp]], XVGi)
  ijk[2] = int[1]
  
  return, ijk
end








function get_ij, B, XVG
;+
;
;FUNCTION: GET_IJ(B, XVG)
;
;PURPOSE:
;  To identify the indices corresponding to a specific baseline.
;
;INPUT:
;  B = Baselines (in pitch)
;  XVG = Fiber positions in v-groove (mm)
;
;OUTPUT:
;  [A,C] = The indices involved in the baseline.
;
;NOTES:
;
;-
  nXVG = n_elements(XVG)


  if n_elements(B) GT 1 then begin
     for a = 0, nXVG -2 do begin
        for c = a+1, nXVG-1 do begin
           if where(B EQ XVG[c] - XVG[a]) NE -1 then begin    
              return, [a,c]
           endif 
        endfor
     endfor
  endif else begin
     for a = 0, nXVG -2 do begin
        for c = a+1, nXVG-1 do begin
           if XVG[c] - XVG[a] EQ B then begin    
              return, [a,c]
           endif 
        endfor
     endfor
  endelse
end








function get_CLO, nXVG
;+
;
;FUNCTION: GET_CLO(nXVG)
;
;PURPOSE:
;  This function will obtain the indices of baselines to be used for 
;  triangles of closure phase calculations.
;
;INPUT:
;  nXVG = The number of positions along the v-groove, # of fibers.
;
;OUTPUT:
;  CLO = An array of the all of the indices of the baselines to be used
;        to identify the triangles for the calculation of the closure
;        phases. 
;
;NOTES: (None)
;
;-

  CLO = [0,0,0]
  for i = 0, nXVG-3 do begin
     for j = i+1, nXVG-2 do begin
        for k = j+1, nXVG-1 do begin
           CLO = [[CLO], [i*nXVG-((i+1)*i)/2.+j-i, $
                          j*nXVG-((j+1)*j)/2.+k-j, $
                          i*nXVG-((i+1)*i)/2.+k-i]]
        endfor
     endfor
  endfor


  CLO = CLO[*,1:*]
  CLO -= 1  ;accounting for Yorick->IDL


  return, CLO
end








function mean_estimator, m, e
;+
;
;FUNCTION: MEAN_ESTIMATOR(M,E)
;
;PURPOSE:
;  To calculate an error-weighted mean and find the final error
;  involved. 
;
;INPUT:
;  M = Mean
;  E = Error
;
;OUTPUT:
;  Returns a vector r:
;    r[0] = mean of the series m weighted by the errors e
;    r[1] = final error on the mean value (corrected if chi2 GT 1)
;    r[2] = chi2 (correction applied if chi2 GT 1)
;
;NOTES: (NONE)
;
;-

    r = fltarr(3)
    indr = where(e NE 0)
    if total(indr) NE -1 then begin
       nlb = n_elements(indr)
       if nlb GT 1 then begin
          mean = total(m[indr]/e[indr]^2)/total(1./e[indr]^2)
          r[0] = mean
          
          chi2 = 1./(nlb-1)*total((mean-m[indr])^2/e[indr]^2)
          r[2] = chi2
          
          err = sqrt(1./total(1./e[indr]^2))
          if chi2 GE 1 then r[1] = sqrt(err^2*chi2) else r[1] = err
       endif
    endif


    return, r
end








function warp_around_zero, m
;+
;
;FUNCTION: WARP_AROUND_ZERO(M)
;
;PURPOSE:
;
;
;INPUT:
;  m = a vector of n phases in degrees between -+180 deg
;
;OUTPUT:
;  Mz = A vector of n+1 elements,
;       mz[0] = mean of the phases derived from complex vectors
;       mz[rest] = the same phases as m centered on zero mean value
;
;NOTES: (NONE)
;
;-

   p = fltarr(2, n_elements(m))
   p[0,*] = cos(m*!dtor)  &  p[1,*] = sin(m*!dtor)

   pm = avg(p,1)
   m0 = (atan(pm[1], pm[0]))[0]
   
   p[0,*] = cos(m*!dtor-m0)  &  p[1,*] = sin(m*!dtor-m0)
   md = atan(p[1,*], p[0,*])*!radeg

   Mz = fltarr(1+n_elements(md))
   Mz[0] = m0*!radeg  &  Mz[1:*] = md

   return, Mz
end








function CP_binary, a, UV, indCLO
;+
;
;FUNCTION: CP_BINARY(A, UV, INDCLO)
;
;PURPOSE:
;
;INPUT:
;  a = Coefficients provided by the main function,
;    a[0] = flux ratio
;    a[1] = x position of companion with respect to partner
;    a[2] = y position of companion with respect to partner
;  UV = array of baseline lengths, 
;    UV[0,*] = frequency
;    UV[1,*] = x/y
;    UV[2,*] = wavelength
;  indCLO = The indices of the closure phase triangles.
;
;OUTPUT:
;  Returns closure phases in degrees corresponding to a binary system
;  with the above input parameters.
;
;
;NOTES:
;  -Have not yet implemented this, just seems useful for performing
;   comparisons later.
;
;-
  r = abs(a[0])  &  xp = a[1]  &  yp = a[2]

  Vc = 1./(1.+r)*(1+r*exp(-2.*complex(0,1)*!pi*(xp*UV[*,0,*]+yp*UV[*,1,*])))

  Bi = Vc[indCLO[0,*],*]*Vc[indCLO[1,*],*]*conj(Vc[indCLO[2,*],*])
  CP = atan(imaginary(Bi), real_part(Bi))*!radeg

  return, CP
end
