pro get_best_params, files, fibres_init, DIR=dir, RIGHT=right, $
                          BIN=bin, DISPLAY=display
;+
;
;FUNCTION: GET_BEST_PARAMS, DIR, FILES, FIBRES_INIT, BIN, METHOD, 
;                          [.../RIGHT, /DISPLAY, /COMPARE]
;PURPOSE:
;  Returns the complex bispectrum of data based on calculation of the
;  V2PM and P2VM using the smoothed fiber profiles.
;
;INPUT:
;  files = An array of a stack of files to be run through
;          GET_BEST_PARAMS with the given fiber stack. 
;  fibres_init = The four character string that identifies the smoothed
;                fiber stack to be used in the construction of the P2VM
;                matrix used in the fit the carrying waves.
;  /dir = A path given to where all the relevant files are saved and
;         where the results will also be saved.
;  /right = Working with the right side of the image. 
;  /bin = The bin size used to combine columns in the image. 
;  /display = Plot results.
;
;OUTPUT:
;
;NOTES:
;files should be an array of all the subsequences of each sequence in
;the form: "1040_20111013_Aldebaran_0"
;-  ;save_name should be an array just like files
;-  ;file_freq should be an array just like files
;
;-
  ;Setting up important variables
  common M_constants, p_det, foc, p_VG, save_name
  p_det = .016  &  foc = 150.  &  p_VG = .037
  init = strmid(files[0], 0, 4) 
  if ~keyword_set(dir) then dir = 'data/'
  if keyword_set(right) then part = 'right' else part = 'left'
  save_name = strcompress(dir +files+'_'+ part, /remove_all) 


  print, "GET BEST PARAMS..."
  if bin EQ 0 then begin
     bin = ""
     print, "Data without binning"
  endif else begin
     bin = "_bin"+string(bin)
     print, "Binned data with nbin =", string(nbin)
  endelse




  ;Loading in the fiber files...  
  file_ffit = (file_search(dir+fibres_init+"*"+part+$
                          bin+"_fitted_smoothed_stack.fits.gz"))[0]


  if file_ffit EQ "" then begin
     print, "ERROR: No spatial fiber file"
     goto, theend
  endif
  if n_elements(file_ffit) GT 1 then begin
     print, "ERROR: Several matching fiber files exist:"
     print, transpose(file_ffit)
     stop
  endif else begin
     F = mrdfits(file_ffit, 0, fib_hdr)
     fibers = F  &  siz = size(F)  &  nlf = siz[1]  &  nF = siz[3]
  endelse




  ;Loading in the spatial frequencies...
  file_freq = file_search(dir+files+'_'+part+bin+"_Spatial_freq.sav")
  if n_elements(file_freq) NE 10 then goto, theend ;temp stopgap so can sleep
  if (file_freq)[0] EQ "" then begin
     print, "ERROR: No spatial frequencies file"
     goto, theend
  endif
  restore, file_freq[0] ;will restore: p_VG, XVG, Sig, Pvg




  ;Finding the triangle indices...
  nXVG = n_elements(XVG)  &  nB = (nXVG)*(nXVG-1)/2.
  indCLO = get_CLO(nXVG)  &  nCP = (size(indCLO))[2]


  IJK = fltarr(3, nCP)
  XVGint = round(XVG)
  Bvg = get_UV(XVG)  &  Bvgint = round(Bvg)
  for a = 0, nCP -1 do  IJK[*, a] = get_ijk(a, indCLO, XVGint, Bvgint)




  ;Opening up a log to keep track of the calculations...
  filetxt = strcompress(dir+strmid(files[0], 0, strlen(files[0])-2)+$
                        "_"+part+bin+"_GBPlog.txt", /remove_all)
  openw, 110, filetxt
  

  N = n_elements(files)
  for s = 0, N-1 do begin
     ;Loading in the data file...
     printf, 110, files[s]+"-------------------------------"
     file = dir+files[s]+'_'+part+'_dewarped'+bin+'.fits.gz'
     data = mrdfits(file, 0, img_hdr) 
     siz = size(data)  &  nl = siz[1]  &  np = siz[2]  &  ni = siz[3]


     ;Getting the indices to align everything by absolute wavelength...
     if s EQ 0 then cut = align_lambda(nl, nlf, img_hdr, fib_hdr, $
                                       right=keyword_set(right))


     ;Finding the matrices to obtain the carrying waves...
     fits_file = save_name[s]+'_P2VM-V2PM.fits' 
     if file_search(fits_file) EQ "" then begin
        ;Finding the P2VM...
        fibers = F
        M = get_M_fib(fibers, file_freq[s], cut)
        if M[0] EQ 0 then goto, theend

        ;Finding the V2PM...
        params = get_best_param(data, M, img_hdr, cut, matrix=V2PM) 
        if n_elements(params) EQ 1 then begin
           fxaddpar, img_hdr, 'STATUS', 'REJECT', $
                     ' COULD NOT CREATE V2PM'
           print, 'Failed to fit carrying waves.'


           printf, 110, 'Failed to fit carrying waves.'
           printf, 110
           goto, next_file
        endif else begin
           fxaddpar, img_hdr, 'STATUS', fxpar(img_hdr,'STATUS'),$
                     ' V2PM GENERATED'
           print, 'Fitted carrying waves!'
           
           
           matrices = [[[M]],[[V2PM]]]
           writefits, fits_file, matrices, img_hdr, /compress
        endelse
     endif else begin
        restore, save_name[s]+'_MU.sav'  &  goto, calc_BI
     endelse




     ;Finding and saving Mu...
     Re = params[*,0:nB-1,*]  &  Im = params[*,nB :2*nB-1,*]
     mu = dcomplex(Re, Im)
     dummy = size(temporary(Re))  &  dummy = size(temporary(Im)) ;Freeing memory

     file = strcompress(save_name[s]+"_MU.sav", /remove_all)
     save, mu, filename = file

     printf, 110, "MU_ij saved as:"                
     printf, 110, "   "+save_name[s]+"_MU.sav"    
     printf, 110, "   under the variable MU"      
     printf, 110




     ;Finding and saving the bispectrum...
     calc_BI: Bi = mu[*,indCLO[0,*],*] * $
          mu[*,indCLO[1,*],*] * $
          conj(mu[*, indCLO[2,*],*]) 
     file = strcompress(save_name[s]+"_Bi.sav", /remove_all)
     save, Bi, XVG, img_hdr, filename = file
     dummy = 0*temporary(mu)  &  dummy = 0*temporary(Bi)


     printf, 110, "Bispectrum saved as:"           
     printf, 110, "   "+save_name[s]+"_Bi.sav"    
     printf, 110, "   under the variable Bi"      
     next_file: printf, 110
  endfor



  
  close, 110
  theend: print
end








function get_ij, BL, XVG
;+
;
;FUNCTION: GET_IJ(B, XVG)
;
;PURPOSE:
;  To find the indices of the positions in the v-groove corresponding to
;  a specific baseline.
;
;INPUT:
;  B = baselines (in pitch)
;  XVG = fiber positions in v-groove (mm)
;
;OUTPUT:
;
;NOTES:(NONE)
;
;-

  nXVG = n_elements(XVG)  &  nBL = n_elements(BL)
  if nBL GT 1 then begin
     inds = intarr(2,nBL)
     

     for a = 0, nBL-1 do begin
        for b = 0, nXVG -2 do begin
           for c = b+1, nXVG-1 do begin
              test = where(BL EQ abs(XVG[c]-XVG[b]))
              if test NE -1 then inds[*,a] = [b,c]
           endfor
        endfor
     endfor
   

     return, inds
  endif else begin
     

     for a = 0, nXVG -2 do begin
        for b = a+1, nXVG-1 do begin
           if abs(XVG[b] - XVG[a]) EQ BL then return, [a,b]
        endfor
     endfor
  endelse
end








function get_M_fib, fiber, freq_file, lineup
;+
;
;FUNCTION: GET_M_FIB(FIBER, XVG, SIG, PVG, P_DET, FOC)
;
;PURPOSE:
;  This function will take the corrected fiber data and find the P2VM.  
;
;INPUT:
;  fiber = the deform-corrected fiber file from correct_deform
;  XVG = the (fitted?) baseline array (unitless)
;  Sig = the wavenumber array (inverse mm)
;  Pvg = the fitted vgroove pitch as a function of pixel (micron with
;        Elsa, I think we should use mm for consistency)
;  p_det = the dimensions of the pixel (mm)
;  foc = the focal length (mm)
;
;OUTPUT:
;  M = P2VM
;
;NOTES:
;  -Have had issues in the past with the calculation of the variable E,
;   due to taking the square root of fiber values that can be
;   negative. Due to this issue, a lot of time was put into ensuring
;   that negative values are not present in the fibers by really
;   exploring the dark correction and smoothing of the fibers.
;  -Bvg should by nl x nB 
;  -Last column of M is uncorrelated term, with no background
;  -Freq should also be nl x nB, no idea why Elsa's using four indices 
;
;-
  common M_constants
  restore, freq_file

  ; Loading in variables...
  B = get_UV(round(XVG))  &  nB = n_elements(B)
  Bvg = reform(B, 1, nB) ## reform(Pvg, n_elements(Pvg), 1) 

  ;This should not be necessary when the full pipeline is working properly
  if lineup[1,0] GT n_elements(Bvg)-1 then begin
     if lineup[1,0]-n_elements(Bvg)+1 GT 2 then begin
        print, "Wavelaw doesn't match, no Bi/CP"
        goto, theend
     endif else lineup[1,*] -= lineup[1,0]-n_elements(Bvg)+1
  endif

  Bvg = Bvg[lineup[0,0]:lineup[1,0],*]

  fiber = fiber[lineup[0,1]:lineup[1,1],*,*]
  siz = size(fiber)  &  nl = siz[1]  &   np = siz[2]
  

  ; Creating the matrix and adding the uncorrelated term...
  M = dblarr(2*nB+1, np, nl)
  M[2.*nB,*,*] = transpose(total(fiber, 3)) 
  x = (findgen(np) - floor(np/2))*p_det


  ; Obtaining the spatial frequencies...
  freq = Bvg*0.
  for i = 0, nB -1 do begin
     freq[*,i]  = Bvg[*,i]*1000./(1/sig*foc) 
  endfor


  ;Obtaining the pitch pairs...
  sorted_locs = intarr(2,nB)
  for i =0, nB-1 do sorted_locs[*,i] = get_ij(B[i], fix(XVG))

  
  ;Adding in the carrying waves
  for i = 0, nl-2 do begin     
     for j = 0, nB-1 do begin
        loc = sorted_locs[*,j]
        E = reform(sqrt(fiber[i,*,loc[0]]*fiber[i,*,loc[1]]))
        M[j,*,i] = cos(2*!pi*freq[i,j]*x)*E
        M[nB+j,*,i] = sin(2*!pi*freq[i,j]*x)*E
     endfor
  endfor


  return, M
  theend: return, 0
end








function get_M_Env, data_S, XVG
;+
;
;FUNCTION: GET_M_FIB(FIBER, XVG, SIG, PVG, P_DET, FOC)
;
;PURPOSE:
;  This function will take the fitted envelope information and find the
;  P2VM. (Useful in cases where the fiber profiles are insufficient)  
;
;INPUT:
;  fiber = the deform-corrected fiber file from correct_deform
;  XVG = the (fitted?) baseline array (unitless)
;  Sig = the wavenumber array (inverse mm)
;  Pvg = the fitted vgroove pitch as a function of pixel (micron with
;        Elsa, I think we should use mm for consistency)
;  p_det = the dimensions of the pixel (mm)
;  foc = the focal length (mm)
;
;OUTPUT:
;  M = P2VM
;
;NOTES:
;  -Have had issues in the past with the calculation of the variable E,
;   due to taking the square root of fiber values that can be
;   negative. Due to this issue, a lot of time was put into ensuring
;   that negative values are not present in the fibers by really
;   exploring the dark correction and smoothing of the fibers.
;  -Freq result should be nl x nB, no idea why Elsa's using four
;   indices.
;
;-
  ;Gathering up variables...
  common M_constants
  restore, freq_file
  restore, repstr(freq_file, 'Spatial_freq', 'env')  


  siz = size(data_S)  &  nl = siz[1]  &  np = siz[2]  &  ni = siz[3]
  nXVG = n_elements(XVG)
  B = get_UV(round(XVG))  &  nB = n_elements(B)
  Bvg = reform(B, 1, nB) ## reform(Pvg, n_elements(Pvg), 1)


  ;Recreating the envelope profiles...
  x = dindgen(np)  &  Env = dblarr(nl,np)
  for i = 0, nl-1 do begin
     ev = reform(env_vars[i,*])
     Env[i,*] = ev[3]+ev[0]*exp(-0.5*((x-ev[1])/ev[2])^2)
  endfor


  ;Getting the frequencies of the waves...
  freq = Bvg*0.
  for i = 0, nB -1 do begin
     freq[*,i]  = Bvg[*,i]/(1/sig*foc) 
  endfor

  
  ;Filling in the P2VM...
  M= dblarr(2*nB+1, np, nl)  &  x = (x - floor(np/2))*p_det
  for i = 0, nl-1 do begin
     for j = 0, nB -1 do begin
        IJ = get_ij(B[i], round(XVG))
        for k = 0, np-1 do begin
           M[j,k,i] = cos(2*!pi*freq[j,i]*x[k])*Env[i,k]
           M[nB+j,k,i] = sin(2*!pi*freq[j,i]*x[k])*Env[i,k]
        endfor
     endfor
     M[2.*nB,*,i] = Env[i,*]
  endfor


  return, M
end








function get_best_param, data_S, M, hdr, lineup, MATRIX=matrix
;+
;
;FUNCTION: GET_BEST_PARAM(DATA_S, M, [.../SAVEIT])
;
;PURPOSE:
;  To fit the real and imaginary parts of the interferogram using the
;  P2VM and a singular value decomposition. 
;
;INPUT:
;  data_S = The array of data to be from which the carrying waves will
;           be fitted.
;  M = The P2VM matrix created from the fiber profiles or carrying wave 
;      envelopes. Default = fiber profiles 
;  hdr = Set to the header, the updated header will be returned to
;         the variable this keyword is set to.
;  /matrix = If keyword is set, the V2PM matrix will be returned to the 
;            variable this keyword is set to.
;
;OUTPUT:
;  a = The fitted real and imaginary components of the interferogram
;
;NOTES:
;  -M = 2nB+1 x np x nl, s = 2nB+1, u = 2nB+1 x 2nB+1, vt = np x np
;  -Mi = np x 2nB+1
;
;-

  ;Setting major constants...
  data_S = data_S[lineup[0,0]:lineup[1,0],*,*]
  siz = size(data_S)  &  nl=siz[1]  &  np=siz[2]  &  ni=siz[3]


  ;Checking for the inevitable bad results...
  if total(data_S EQ data_S) NE n_elements(data_S) then $
     stop, "This data has NaNs"
  if total(M EQ M OR M EQ 0) NE n_elements(M) then $
     stop, "This fiber matrix has NaNs"


  ;Setting up arrays for the fit...
  ns = (size(M))[1]  &  a = dblarr(nl, ns, ni)
  fit = data_S*0  &  matrix = M*0

  
  for i = 0, nl-1 do begin
     ;Using SVDC to find the V2PM...
     svdc, reform(temporary(M[*,*,i])), s, u, vt 

     
     ;Checking for the inevitable NaNs...
     catch, NaN_bug
     if total([-142,-137,-134] EQ NaN_bug) NE 0 then begin
        stop, "SVDC has somehow introduced NaNs into the fit"
        catch, /cancel
     endif else if NaN_bug NE 0 then catch, /cancel
     
     
     ;Calculating the V2PM... 
     si = dblarr(ns, ns)  &  for q = 0, ns-1 do si[q,q] = 1./s[q]
     V2PM = vt##si##transpose(u) 
     matrix[*,*,i] = V2PM




     ;If NaNs exist in the V2PM... 
     lim = 10  &  chk = (i LE lim OR i GE nl -lim)
     if total(V2PM NE V2PM) NE 0 AND chk then begin
        print, "ERROR: We have NaNs!"
        print, "       Rectifying in the jankiest fashion imaginable!"
        dead_inds = [-1,-1]  
        
        ;...we'll just cut off the edges...
        if i LE lim then begin
           data_S = data_S[i+1:*,*,*]
           fit = fit[i+1:*,*,*]
           a = a[i+1:*,*,*]
           dead_inds[0] = i
        endif else begin
           data_S = data_S[0:i-1,*,*]
           fit = fit[0:i-1,*,*]
           a = a[0:i-1,*,*]
           dead_inds[1] = i
        endelse
        goto, outta_here_really                   
     endif else if total(V2PM NE V2PM) GT 0 then stop, "We've got issues." 
     ;...unless there are more than that
      




     ;Finding the carrying waves...
     for j = 0, ni-1 do begin
        a_best_fit = V2PM##data_S[i,*,j] 
        a[i,*,j] = a_best_fit
     endfor


     ;I still don't know why this is necessary
     outta_here: if i EQ nl-1 then goto, outta_here_really 
  endfor

  
  

  ;Putting lost indices into the header
  outta_here_really: print
  if n_elements(dead_inds) EQ 0 then dead_inds = [-1,-1]
  dead_inds = strcompress(strjoin(string(dead_inds),','),/remove_all)
  fxaddpar, hdr, 'DEADINDS', dead_inds, $
            ' EDGES OF IMAGE FITTED IN GBP', before = 'STATUS'


  return, a 
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








function align_lambda, i_nl, f_nl, i_hdr, f_hdr, RIGHT=right
;+
;
;FUNCTION: ALIGN_LAMBDA(I_NL, F_NL, I_HDR, F_HDR, [.../RIGHT])
;
;PURPOSE:
;  To find the indices over which the image and fiber files share the
;  save wavelengths (on the assumption that they initially obey the same
;  wavelaw solution).
;
;INPUT:
;  i_nl, f_nl = The number of columns in the images and fibers.
;  i_hdr, f_hdr = The image/fiber headers.
;  /right = If keyword is set, right side of image is in use.
;
;OUTPUT:
;  align = The indices for the image and fiber files for use in
;          calculating the P2VM in get_M_fib.
;
;NOTES:(None)
;
;-
  align = intarr(2,2)           ;top = image, bottom = fiber

  
  ;Obtaining the indices of the edge_correct cuts...
  edg = fxpar(i_hdr, 'EDGES')  &  edg = strsplit(edg,':',/extract)
  fedg = fxpar(f_hdr, 'EDGES')  &  fedg = strsplit(fedg,':',/extract)




  ;If image is right side, accounting for dark band...
  if keyword_set(right) then begin
     bnd = fxpar(i_hdr, 'DBAND')  &  bnd = strmid(bnd, 4, 3)
     fbnd = fxpar(f_hdr, 'DBAND')  &  fbnd = strmid(fbnd, 4, 3)
     
     band = fix([bnd,fbnd]) + 1
  endif else band = [0,0]
  
  
  
  
  ;Getting the left hand limits...
  edges = fix([edg[0], fedg[0]]) + band + 1  
  left = max(edges, loc)
  align[0, loc] = 0  &  align[0, (~loc)] = left-edges[(~loc)]

  

  ;Getting the right hand limits...
  edges2 = edges + fix([edg[1], fedg[1]]) - 1  
  right = min(edges2,loc)

  lens = [i_nl, f_nl] - reform(align[0,*]) - 1

  right = min([right - edges[loc] - 1, lens+align[0,loc]])
  shift = align[0,(~loc)]-align[0,loc]
  align[1, loc] = right  &  align[1, (~loc)] = right-shift
  if align[1,loc]-align[0,loc] NE align[1,~loc]-align[0,~loc] then begin
     bit = min([align[1,loc]-align[0,loc],align[1,~loc]-align[0,~loc]])
     align[1,*] = align[0,*]+bit
  endif



  return, align
end
