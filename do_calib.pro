pro easy_call, seq_tgt, seq_cal, FNAME=fname, RIGHT=right, S=s, B=b, DIR=dir
;+
;
;FUNCTION: EASY_CALL, SEQ_TGT, SEQ_CAL, [...,FNAME='TAG',/RIGHT]
;
;PURPOSE:
;  This function makes calling do_calib easier by loading the requisite
;  values beforehand and matching up the wavelaws and bispectra.
;
;INPUT:
;  seq_tgt, seq_cal = The four character identifiers, date and star of
;  the target and calibration files. Example: 1000_20111019_Vega
;  /fname = A tag to identify the saved calibrated CPs with.
;  /right = If keyword is set, then working with right half of image.
;  /S = If keyword is set, then
;  /b = If keyword is set, then 
;  /dir = 
;
;OUTPUT: (NONE)
;
;NOTES: (NONE)
;
;-
  if keyword_set(right) then part = 'right' else part = 'left'
  if ~keyword_set(S) then S = 1  &  if ~keyword_set(b) then b = 100
  if ~keyword_set(dir) then dir = 'data/'

  ;Gathering the target info...
  info_tgt = align_columns(seq_tgt, part, dir, S=S, b=b)
  CPm_tgt = info_tgt.CP_avg
  L_tgt = info_tgt.wavelaw
  SNR_tgt = info_tgt.SNR
  dummy = size(temporary(info_tgt))


  ;Gathering the calibrator info...
  info_cal = align_columns(seq_cal, part, dir, S=S, b=b)
  CPm_cal = info_cal.CP_avg
  L_cal = info_cal.wavelaw
  SNR_cal = info_cal.SNR
  dummy = size(temporary(info_tgt))




  ;Lining up the wavelaws (even if not linear)...
  match_wavelengths, L_tgt, L_cal, SNR_tgt, SNR_cal, CPm_tgt, CPm_cal


  ;Generating the save_name and calibrating the closure phases...
  date = strmid(seq_tgt, 5, 8)
  if keyword_set(S) then S = 'S'+string(S) else S = 'S1'
  if keyword_set(b) then b = 'b'+string(b) else b = 'b100'
  if keyword_set(fname) then fname += '_' else fname = '' 
  fname = strcompress(date+'_'+S+'_'+b+'_'+fname,/remove_all)
  fname = strmid(seq_tgt,0,4)+'_'+strmid(seq_cal,0,4)+'_'+fname+part


  do_calib, CPm_tgt, CPm_cal, L_tgt, L_cal, SNR_tgt, SNR_cal, fname, 2
end








pro do_calib, CPm_tgt, CPm_cal, L_tgt, L_cal, SNR_tgt, SNR_cal, $
              savename, multi
;+
;
;FUNCTION: DO_CALIB, CPM_TGT, CPM_CAL, L_TGT, L_CAL, SNR_TGT, SNR_CAL,
;                    SAVENAME, MULTI 
;
;PURPOSE:
;  This function will calibrate the closure phases of a target to those
;  of the calibrator star
;
;INPUT:
;  CPm_tgt, CPm_cal = Array of dimensions nl by nCP by 4 by ns
;                     containing the closure phases to calibrate.
;                     (# of columns by # of CP triangles by 4 by # of
;                      sequences) 
;  L_tgt, L_cal = Array of dimensions nl by ns containing the fitted
;                 wavelaws aligned and sliced to share the same
;                 coverage.
;  SNR_tgt, SNR_cal = Array of dimensions nl by nCP by ns containing the
;                     SNRs from get_CPmoy.
;  savename = The calibrated closure phases will be saved in: 
;             savename +'_CPm.sav' 
;  multi = how often a calibrated CP is saved:
;           0 if only an averaged calibrated CP is saved
;           1 if every sequence of calibrated CP is saved
;           2 if only different sequences for the target are saved, with
;             the calibrator's sequences averaged
;
;
;OUTPUT:

;
;NOTES:
;  -As the alignment of the wavelaws still will contain minor
;   differences (although much better than before), I maintained the
;   average of the two wavelaws to achieve the final wavelaw applied.
;-


  siz = size(CPm_tgt)  &  nl = siz[1]  &  nCP = siz[2]  
  ns_tgt = siz[4]  &  ns_cal = ns_tgt




  ;Cleanly working out the avg where there is GT 1 sequence...
  CPm = fltarr(nl,nCP,2)
  CP_cal_moy = CPm  &  CP_tgt_moy = CP_cal_moy

  
  if ns_tgt GT 1 then begin
     for k = 0, nCP-1 do begin
        for l = 0, nl-1 do begin
           W = Warp_around_zero(CPm_tgt[l,k,0,*])
           moy = W[0]
           CPnew = W[1:*]
           CP_tgt_moy[l,k,*] = (mean_estimator(CPnew, CPm_tgt[l,k,1,*]))[0:1]
           CP_tgt_moy[l,k,0] += moy
        endfor
     endfor
  endif else CP_tgt_moy = reform(CPm_tgt[*,*,0:1,0])


  if ns_cal GT 1 then begin
     for k = 0, nCP-1 do begin
        for l = 0, nl-1 do begin
           W = (Warp_around_zero(CPm_cal[l,k,0,*]))
           moy = w[0]
           CPnew = w[1:*]
           CP_cal_moy[l,k,*] = (mean_estimator(CPnew, CPm_cal[l,k,1,*]))[0:1]
           CP_cal_moy[l,k,0] += moy
        endfor
     endfor
  endif else CP_cal_moy= reform(CPm_cal[*,*,0:1,0])




  ;Finding the average results across all sequences...
  CPm[*,*,0] = CP_tgt_moy[*,*,0]-CP_cal_moy[*,*,0]
  CPm[*,*,1] = sqrt(CP_tgt_moy[*,*,1]^2+CP_cal_moy[*,*,1]^2)
  L = (total(L_tgt,2)+total(L_cal,2))/(ns_tgt+ns_cal)
  SNR = []

  file = 'data/'+savename+'_CPm.sav'
  save, filename = file, CPm, L, SNR
  print, "CPm, SNR and L saved in ", file




  ;Finding the results sequence by sequence...
  if multi EQ 1 then begin
     if ns_tgt NE ns_cal then begin
        print, "ERROR: mode multi = 1 requires that the number of" + $
               " sequences in CPm_cal and CPm_tgt are " + $
               "greater than 1 and are equal."
        CPm = 0 
        goto, just_end
     endif

     for i = 0, ns_tgt-1 do begin
        CPm[*,*,0] = CPm_tgt[*,*,0,i] - CPm_cal[*,*,0,i]
        CPm[*,*,1] = sqrt(CPm_tgt[*,*,1,i]^2 + CPm_cal[*,*,1,i]^2)
        L = (L_tgt[*,i]+L_cal[*,i])/2.
        SNR = []

        num = strcompress(string(i),/remove_all)
        file = 'data/'+savename+'_CPm'+num+'.sav'
        save, filename = file, CPm, L, SNR
        print, "CPm, SNR and L saved in ", file
     endfor
  endif
  



  ;Finding the results with the calibrator sequences averaged...
  if multi EQ 2 then begin
     CPcalib = CPm

     for k = 0, nCP-1 do begin
        for l = 0, nl-1 do begin
           W = warp_around_zero(CPm_cal[l,k,0,*])
           moy = W[0]
           CPnew = W[1:*]

           CPcalib[l,k,*] = (mean_estimator(CPnew, CPm_cal[l,k,1,*]))[0:1]
           CPcalib[l,k,0] += moy
        endfor
     endfor

     CPcalib[*,*,1] = sqrt(reform(total(CPm_cal[*,*,1,*]^2,4)))/ns_cal
     Lc = avg(L_cal,1)

     for i = 0, ns_tgt-1 do begin
        CPm[*,*,0] = CPm_tgt[*,*,0,i] - CPcalib[*,*,0]
        CPm[*,*,1] = sqrt(CPm_tgt[*,*,1,i]^2 + CPcalib[*,*,1]^2)
        L = (L_tgt[*,i]+Lc)/2.
        SNR = []

        num = strcompress(string(i),/remove_all)
        file = 'data/'+savename+'_CPm'+num+'.sav'
        save, filename = file, CPm, L, SNR
        print, "CPm, SNR and L saved in ", file
     endfor
  endif

  just_end: print
end








function warp_around_zero, m
;+
;
;FUNCTION: WARP_AROUND_ZERO(M)
;
;PURPOSE:
;  To properly average the closure phases while accounting for the 0
;  values properly.
;
;INPUT:
;  m = A vector of n phases in degrees between -+180 degrees. 
;
;OUTPUT:
;  Mz = a vector of n+1 elements, where:
;       mz[0] = Mean of the phases derived from complex vectors.
;       mz[1:*] = The same phases as m centered on zero mean value.
;
;NOTES: (NONE)
;
;-
  nm = n_elements(m)
  p = fltarr(2, nm)  
  p[0,*] = cos(m*!dtor)  &  p[1,*] = sin(m*!dtor)


  pm = avg(p,1)  &  m0 = (atan(pm[1], pm[0]))[0]
  p[0,*] = cos(m*!dtor-m0)  &  p[1,*] = sin(m*!dtor-m0)
  md = atan(p[1,*], p[0,*])*!radeg


  Mz = fltarr(1+nm)  
  Mz[0] = m0*!radeg  &  Mz[1:*] = md

  return, Mz
end








function mean_estimator, m, e
;+
;
;FUNCTION: MEAN_ESTIMATOR(M,E)
;
;PURPOSE:
;  To perform an average taking into account the estimated error.
;
;INPUT:
;  m =
;  e = 
;
;OUTPUT: 
;  r = A vector such that:
;      r[0] = mean of the series m weighted by the errors e
;      r[1] = final error on the mean value (corrected if chi2 GT 1) 
;      r[2] = chi2 (correction applied if chi2 GT 1)
;
;NOTES: (NONE)
;
;-

    r = fltarr(3)
    indr = where(e NE 0)

    if total(indr) NE -1 then begin
        if n_elements(indr) GT 1 then begin
            mean = total(m[indr]/e[indr]^2)/total(1./e[indr]^2)
            r[0] = mean

            nlb = n_elements(indr)
            chi2 = 1./(nlb-1)*total((mean-m[indr])^2/e[indr]^2)
            r[2] = chi2

            err = sqrt(1./total(1./e[indr]^2))
            if chi2 GE 1 then r[1] = sqrt(err^2*chi2) else r[1] = err
        endif
    endif

    return, r
end








function align_columns, seq, part, dir, S=s, B=b
;+
;
;FUNCTION: ALIGN_COLUMNS(SEQ, PART)
;
;PURPOSE:
;  To perform the work of loading in and lining up the wavelaw, closure
;  phase and SNR arrays based on the columns that were tossed out during
;  the process of fitting the real and imaginary parts of the
;  interferogram. 
;
;INPUT:
;  seq = The four character identifier string pertaining to the file.
;  part = A string specifying the left or right part of the image is in
;         use. 
;  dir = 
;  /S = The threshold used in GET_BEST_PARAMS. If keyword is not set,
;       then the default value of S = 1 will be used.
;  /b = The bundle size used in GET_BEST_PARAMS. If keyword is not set,
;       then the default value of b = 100 will be used. 
;
;OUTPUT:
;  struct = A structure containing the stacked average closure phase
;           array, wavelaw array, and SNR array with the tags CP_avg,
;           wavelaw and SNR. 
;
;NOTES: (NONE)
;
;-
  if keyword_set(S) then S_str = strcompress(string(S),/remove_all) $
  else S_str = '1'
  if keyword_set(b) then b_str = strcompress(string(b),/remove_all) $
  else b_str = '100'
  



  ;Identifying where the different value arrays will be cut...
  lims = intarr(2,10)
  num = strcompress(string(intarr(10)),/remove_all)
  tag = dir+seq+'*_'+num+'_'+part

  for i = 0, 9 do begin
     if (file_search(tag[i]+'_Bi.sav'))[0] NE "" then begin
        restore, tag[i]+'_Bi.sav'
        hdr = img_hdr
     endif else begin
        if (file_search(tag[i]+'_dewarped.fits.gz'))[0] NE "" then begin
           hdr = headfits(tag[i]+'_dewarped.fits.gz')
        endif else stop, "ERROR: Missing header file"
     endelse

     dead_inds = fix(strsplit(fxpar(hdr, 'DEADINDS'),',',/extract))
     dead_inds[0] += 1  &  if dead_inds[1] NE -1 then dead_inds[1] -= 1
     lims[*,i] = dead_inds
  endfor

  cut = [max(lims[0,*]), min(lims[1,*])] ;For the wavelaw...
  lims[0,*] = cut[0]-lims[0,*]  &  lims[1,*] -= cut[1] ;For all else...




  ;Loading in desired values and slicing appropriately...
  for i = 0, 9 do begin
     restore, tag[i]+'_S'+S_str+'_b'+b_str+'_CPmoy1.sav'
     restore, tag[i]+'_S'+S_str+'_b'+b_str+'_CP.sav'  
     restore, tag[i]+'_Spatial_freq.sav'  

     siz = size(CParr)  &  n_arr = cut[1]-cut[0]+1
     if i EQ 0 then begin
        CPm = make_array([n_arr,siz[2],siz[4], 10], /double)
        L = make_array([n_arr, 10], /double)
        SNRarr = make_array([n_arr, (size(SNR))[2], 10], /double)
     endif
     
     CPm[*,*,*,i] = CParr[lims[0,i]:-lims[1,i]-1,*,1,*]
     L[*,i] = (1./Sig[cut[0]:cut[1]])
     ;SNRarr[*,*,i] = SNR[lims[0,i]:-lims[1,i]-1,*]
  endfor




  ;Creating and returning the structure...
  struct = {CP_avg: CPm, wavelaw: L, SNR: SNRarr}
  return, struct
end








pro match_wavelengths, La, Lb, SNRa, SNRb, CPam, CPbm
;+
;
;FUNCTION: MATCH_WAVELENGTHS, LA, LB, SNRA, SNRB, CPAM, CPBM
;
;PURPOSE:
;  Aligns the 3 arrays used in DO_CALIB based on the individual wavelaws
;  so that all non-linearities are accounted for.
;
;INPUT:
;  La, Lb = Target and calibrator wavelaws.
;  SNRa, SNRb = Target and calibrator SNRs.
;  CPma, CPmb = Target and calibrator averaged closure phases.
;
;OUTPUT:
;  Changes are passed back through input arrays.
;
;NOTES: (NONE)
;
;-

  rng = [(size(La))[1], (size(Lb))[1]]  &  stk = (size(La))[2]
  lis = []  &  indi = []
  for d = 0, stk-1 do begin
     ;Finding the differences between the wavelaws element by element...
     a = fltarr(2,n_elements(La[*,0]))+1  &  a[0,*] = La[*,d]
     b = fltarr(n_elements(Lb[*,0]),2)-1  &  b[*,1] = Lb[*,d]
     dif = abs(a##b)


     ;Finding the best matching pixels between the to wavelaws 
     inds1 = []  &  inds2 = []
     for f = 0, rng[0]-1 do $
        inds1 = [inds1,(where(dif[*,f] EQ min(dif[*,f])))[0]]
     for f = 0, rng[1]-1 do $
        inds2 = [inds2,(where(dif[f,*] EQ min(dif[f,*])))[0]]
     inds = {one: inds1, two: inds2}


     ;Making sure the only matches are the closest...
     dummy = min(rng, loc)
     include = intarr(rng[loc])  &  arr = (inds.(loc))
     for e = 0, rng[loc]-1 do begin
        poss = where(arr EQ e)  &  chk = poss EQ inds1[e]
        include[e] = poss[where(chk EQ 1)] 
     endfor
     include = include[where(include NE -1)]
     n_inc = n_elements(include)  &  inds = intarr(2,n_inc)
     inds[loc,*] = arr[include]  &  inds[(~loc),*] = include

     
     indi = [[indi], [inds]]  &  lis =[lis,n_inc]
  endfor




  ;Getting the ranges of indices to grab out of our index stacks...
  emp = fltarr(2, stk)
  for d = 1, stk-1 do emp[0,d] = total(lis[0:d-1])
  emp[1,*] = emp[0,*]+min(lis)-1



  
  ;Indexing all the arrays...
  nla = rng[0]  &  nlb = rng[1]  &  rng = min(lis)-1
  nCP = n_elements(SNRa[0,*,0])

  for d = 0, stk-1 do begin
     list_a = indi[0,emp[0,d]:emp[1,d]]
     list_b = indi[1,emp[0,d]:emp[1,d]]

     La[0:rng,d] = La[list_a,d]  &  Lb[0:rng,d] = Lb[list_b,d]

     for e = 0, nCP-1 do begin
        SNRa[0:rng,e,d] = SNRa[list_a,e,d] 
        SNRb[0:rng,e,d] = SNRb[list_b,e,d]

        for f = 0, 1 do begin
           CPam[0:rng,e,f,d] = CPam[list_a,e,f,d]
           CPbm[0:rng,e,f,d] = CPbm[list_b,e,f,d]
        endfor
     endfor
  endfor
  
  La = La[0:rng,*]  &  Lb =Lb[0:rng,*]
  CPam = CPam[0:rng,*,*,*]  &  CPbm = CPbm[0:rng,*,*,*]
  SNRa = SNRa[0:rng,*,*]  &  SNRb = SNRb[0:rng,*,*]
  ;BOOM! PIXEL BY PIXEL MATCHING!!
     




  ;Getting the wavelengths consistent in the sequence stacks...
  wavf = La[-1,*] & wav0 = La[0,*]
  fin = min(wavf[where(wavf NE 0)]) & beg = max(wav0)
     
  ns = n_elements(La[0,*])
  spots = fltarr(2, 2*ns)
  for d = 0, ns-1 do begin
     dum = min(abs(la[*,d]-beg),loc)  &  spots[0,d] = loc
     dum = min(abs(la[*,d]-fin),loc)  &  spots[1,d] = loc
     dum = min(abs(lb[*,d]-beg),loc)  &  spots[0,ns+d] = loc
     dum = min(abs(lb[*,d]-fin),loc)  &  spots[1,ns+d] = loc 
  endfor

  ;Finding the maximum overlap and starting position...
  rng = min(spots[1,*]-spots[0,*])  &  locs = findgen(rng+1)
  spots = reform(spots[0,*])  
  
  ;Indexing...
  for d = 0, ns-1 do begin
     list_a = locs+spots[d]  &  list_b = locs+spots[d+ns]
     La[0:rng,d] = La[list_a,d]  &  Lb[0:rng,d] = Lb[list_b,d]
     CPam[0:rng,*,*,d] = CPam[list_a,*,*,d]
     CPbm[0:rng,*,*,d] = CPbm[list_b,*,*,d]
     SNRa[0:rng,*,d] = SNRa[list_a,*,d]
     SNRb[0:rng,*,d] = SNRb[list_b,*,d]
  endfor
  
  ;Shortening...
  La = La[0:rng,*]  &  Lb =Lb[0:rng,*]
  CPam = CPam[0:rng,*,*,*]  &  CPbm = CPbm[0:rng,*,*,*]
  SNRa = SNRa[0:rng,*,*]  &  SNRb = SNRb[0:rng,*,*]
end


