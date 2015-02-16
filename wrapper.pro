pro wrapper, file, LIST=list, RIGHT=right, REBIN=rebin, RERUN=rerun, RUN=run, SHOW=show, STORE=store
;+
;
;FUNCTION: WRAPPER, FILE, [.../RIGHT, /FIBER, REBIN='binsize', /PREP, RERUN = '000000',/SHOW,/STORE]
;
;PURPOSE:
;  This is the wrapper that makes it easy to process images, since it takes
;  care of calling the various procedures with the correct syntax. 
;
;METHOD:
;
;INPUT:
;  The input can be:
;  file = a string (e.g. 'myfile.fits.gz'), or 
;         a string array (e.g. ['myfile1.fits.gz', 'myfile2.fits.gz']
;  /list = if keyword is set, user will be prompted for a list
;          name. Lists are saved in the form: (list name)_filelist.txt 
;          Note: Each image name should be placed on a separate line
;                in that file. Note that all images contained in the
;                list must have the same keywords, exposure time, etc.
;  /right = if keyword is set, the right side of the image file will be
;           processed in all programs for which that is relevant. 
;  /rebin = if keyword is set, image columns will be binned by the
;           string given. Example: rebin = '4' 
;  /rerun = if keyword is set, all programs corresponding to the
;           relevant position in the string will be reprocessed. Input
;           is given as 'rerun='xxxxxx'', where x is either 0
;           (don't rerun) or 1 (rerun) and the six positions in the 
;           string indicate the sequential main functions in the
;           pipeline. Default is that no function is rerun if data
;           products for the function already exist.
;  /run = keyword functions exactly as rerun, but designates which
;         functions will be run or not run in the case where no data
;         products are existant. Takes over the function of the older
;         "prep" keyword 
;  /show = if keyword is set, then images of results will be displayed
;          to user. If using screen, do not set this keyword.
;  /store = if keyword is set, then it will be assumed that the images
;           are being read off a hard-disk, and pertinent results will
;           be moved there as soon as there is no use for them locally.
;
;OUTPUT:(None)
;
;NOTES:(None)
;
;-
    ;If we are working with a list...
    if ~n_params() and keyword_set(list) then begin 
        print, "We've got a list!"

        starlist = ""
        read, "Which list? ", starlist
        list = starlist
        starlist = strcompress(starlist+'_filelist.txt', /remove_all)
        readcol, starlist, images, format='a', /silent
        print, n_elements(images), " files will be processed."  
        
        for i = 0, n_elements(images)-1 do begin
           
           img = images[i]
           pos0 = strpos(img, '/',/reverse_search)+1
           filecrop = strmid(img, pos0,strlen(img)-pos0-10) 
           print, filecrop
          
           dojunk, img, right=keyword_set(right), rebin=keyword_set(rebin),$
                   rerun=rerun, run=run, no_show=~keyword_set(show)
        endfor


    ;If an array, or single file, was inputted...
    endif else if n_params() gt 0 and (size(file))[0] eq 0 then begin 
       if n_elements(file) EQ 1 then print, "We've got a string!" $
       else print, "We've got an array!"

       images = file
       for i = 0, n_elements(images)-1 do begin
          dojunk, images[i], right=keyword_set(right), rebin=keyword_set(rebin),$
                  rerun=rerun, run=run, no_show=~keyword_set(show)
       endfor


    ;If input was improperly given...
    endif else begin
        print, "Improper input."
        print, "The input can be a string (e.g. 'myfile.fits.gz'),"
        print, "  a string array (e.g. ['myfile1.fits.gz', 'myfile2.fits.gz']),"
        print, "  or use the /list keyword to obtain file names from "
        print, "  '(tag)_filelist.txt', where tag is provided by the user."
        print, "Keywords include 'rerun' and 'rebin'. Include 'rerun='xxxxxx'',"
        print, "  where x is either 0 (don't rerun) or 1 (rerun) and the six "
        print, "  positions in the string indicate the sequential main functions"
        print, "  in the pipeline, to rerun any of these functions, and include"
        print, "  'rebin=n' (n is integer) if you want to rebin the fitted fiber image."
    endelse
end








pro dojunk, file, RIGHT=right, REBIN=rebin, RERUN=rerun, RUN=run, no_show=no_show
;+
;
;FUNCTION: DOJUNK, FILE, [.../RIGHT, REBIN='BINSIZE',RERUN='XXXXXX',/PREP,/NO_SHOW]
;
;PURPOSE:
;  This procedure processes one sequence through the full pipeline
;  (barring the use of keywords),
;
;INPUT:
;  file = a string (e.g. 'myfile.fits.gz')
;  /right = If keyword is set, then the right side of the image file
;           will be processed in all programs for which that is
;           relevant.             
;  /rebin = If keyword is set, then the image columns will be binned by
;           the string given. Example: rebin = '4' 
;  /rerun = If keyword is set, then all programs corresponding to the
;           relevant position in the string will be reprocessed. Input
;           is given as 'rerun='xxxxxx'', where x is either 0
;           (don't rerun) or 1 (rerun) and the six positions in the 
;           string indicate the sequential main functions in the
;           pipeline. Default is that no function is rerun if data
;           products for the function already exist.
;  /run = Keyword functions exactly as rerun, but designates which
;         functions will be run or not run in the case where no data
;         products are existant. Takes over the function of the older
;         "prep" keyword. 
;  /no_show = If keyword is set, then all programs will be run in
;             "screen mode", i.e. no commands which require the use of
;             x11 forwarding will be run.
;
;OUTPUT:
;
;NOTES:
;  -This procedure processes one sequence at a time, and thus needs to
;   be called multiple times if there are multiple sequences to be run.
;  -Made a change to this procedure to do checking for prerun files
;   internally in GET_DEFORM so that all files have up-to-date header
;   files. 
;
;-
  STORE = 1
  
  ;This variable is modified via SCOPE_VARFETCH in each function when 
  ;the status of a file suddenly becomes "reject" (or dubious and 
  ;dysfunctional)
  ;Longname to avoid any matches through silliness
  ;10 places for 10 files, Boolean values
   stack_status_token = intarr(10)+1
   
   t0 = systime(/seconds)       ;setting the initial time
   if keyword_set(rerun) then begin
      ns = strlen(rerun)  &  if ns NE 6 then stop, "Improper input: RERUN."
      checked = ~fix(strmid(rerun,indgen(ns),intarr(ns)+1))
   endif else checked = intarr(6)+1
   if keyword_set(run) then begin
      ns = strlen(run)  &  if ns NE 6 then stop, "Improper input: RUN."
      runit = ~fix(strmid(run,indgen(ns),intarr(ns)+1))
   endif else runit = intarr(6)
   if keyword_set(rebin) then rebin = double(rebin) else rebin = 0
   if keyword_set(right) then side = 'right' else side = 'left'
   if keyword_set(no_show) then show = 0 else show = 1


   dirloc = strpos(file, '/', /reverse_search)  
   filedir = strmid(file, 0, dirloc+1) 
   filecrop = strmid(file, dirloc+1, strlen(file)-dirloc-11)  ;strlen(_0.fits.gz') = 10
   s = strlen(filecrop)
   if strmid(filecrop,s-1,1) EQ '_' then filecrop = strmid(filecrop,0,s-1) 
   if keyword_set(store) then begin
      temp = strpos(strmid(filedir,0,strlen(filedir)-1),'/',/reverse_search)
      cut = strmid(filedir,temp+1)
      stor_loc = repstr(filedir,cut,'Processed/'+cut)
   endif
   ;obtain the path information and clip it off




   ;Getting the deformation profile...
   print, ""  &   print, "GETTING THE DEFORMATION PROFILE"
    
   check = strcompress('data/'+filecrop+'_*'+side+'_deform.sav',/remove_all)
   check = file_search(check)  & check = logical_true(check)
   if check*checked[0] OR runit[0] then  print, "Skipping get_deform for ", filecrop $
   else get_deform, filecrop, filedir, right=keyword_set(right),show=show
                                ;want to edit the header in get_deform
                                ;so as to flag files with a warning or
                                ;no go, but must remember these changes
                                ;in correct_deform because of hardcoded
                                ;header stuff--->FIXED

   
   ;Processing relevant darks...
   if ~runit[1] then master_dark, filecrop +'_0.fits.gz', filedir




   ;Correcting the deformation...
   print, ""  &  print, "CORRECTING THE DEFORMATION"
   for i = 0, 9 do begin
      j =strcompress(string(fix(i)),/remove_all)
      checka = file_search('data/'+filecrop+'_'+j+'_'+side+'_dewarped1.fits.gz')
      checkb = file_search('data/'+filecrop+'_'+j+'_'+side+'_dewarped.fits.gz')
      check = (checka+checkb)[0]

      if (check NE "")*checked[1] OR runit[1] OR total(~stack_status_token) then begin
         print, "Skipping correct_deform for: ", filecrop+'_'+j
      endif else begin
         correct_deform, filecrop, j, filedir, right = keyword_set(right),show=show
      endelse      
   endfor

   check = file_search('data/'+filecrop+'_*_'+side+'_dewarped1.fits.gz')
   s = ~runit[1]
   case s*total((check NE "")*stack_status_token) of
      10: begin
         edge_correct, filecrop, 'data/', right = keyword_set(right)
         if keyword_set(store) then begin
            for i = 0, 9 do begin
               nem = filecrop+'_'+strcompress(string(fix(i)),/remove_all)+$
                     '_'+side+'_dewarped.fits.gz'
               im = median(mrdfits('data/'+nem,0,dwhdr),dimension=3)
               nem = repstr(nem,'dewarped','dw_stack')
               writefits, '~/../../'+stor_loc+repstr(nem,'.gz',''), $
                       temporary(im), dwhdr, /compress
            endfor
         endif
      end
      0: print, "Skipping edge_correction for: ", filecrop
      else: begin
         if total(stack_status_token) EQ 10 then begin 
            edge_correct, filecrop, 'data/', right = keyword_set(right), /rescue
            if keyword_set(store) then begin
               for i = 0, 9 do begin
                  nem = filecrop+'_'+strcompress(string(fix(i)),/remove_all)+$
                        '_'+side+'_dewarped.fits.gz'
                  im = median(mrdfits('data/'+nem,0,dwhdr),dimension=3)
                  nem = repstr(nem,'dewarped','dw_stack')
                  writefits, '~/../../'+stor_loc+repstr(nem,'.gz',''), $
                             temporary(im), dwhdr, /compress
               endfor
            endif
         endif else print, "Not all files dewarped, can't correct edges."
      end
   endcase


   ;Moving files around...
   if keyword_set(store) and ~runit[1] then $
      spawn, 'mv data/'+filecrop+'_'+side+'_deform.sav'+' '+stor_loc


   ;Binning if desired...
   if keyword_set(rebin) then begin
      print, ""  &  print, "BINNING DEWARPED FILES"
      rb = strcompress(string(fix(rebin)),/remove_all)

      for i = 0, 9 do begin
         j = strcompress(string(fix(i)),/remove_all)  &  bin_file = filecrop+'_'+j
         check = file_search('data/'+bin_file+'_'+side+'_dewarped_bin'+rb+'.fits.gz')
         if (check[0] NE "")*checked[2] OR runit[2] then print, "Skipping binning file: ", bin_file $
         else binning, filedir, bin_file, rebin
      endfor
   endif




   ;Getting spatial frequencies...
   print, ""  &  print, "OBTAINING SPATIAL FREQUENCIES"

   for i = 0, 9 do begin
      j = strcompress(string(fix(i)),/remove_all)
      gsf_file = filecrop+'_'+j
      if keyword_set(rebin) then bin_str = '_bin'+rb+'_' else bin_str = ''

      check = file_search('data/'+gsf_file+'_'+side+'_'+bin_str+'Spatial_freq.sav')

      if (check[0] NE "")*checked[3] OR runit[3] OR ~stack_status_token[i] then begin
         print, "Skipping get_spatial_freq for: ", gsf_file+bin_str
      endif else get_spatial_freq, gsf_file, right = keyword_set(right), bin = keyword_set(rebin)
      
      if keyword_set(store) and ~runit[3] then $
         spawn, 'mv data/'+gsf_file+'_'+side+'_all_data_DSPm.fits.gz '+stor_loc
   endfor




   ;BINNING CLEANLY IMPLEMENTED UP UNTIL THIS POINT*************************************************************************************************************
   ;Preparing fiber files... 
   print, ""  &  print, "PROCESSING RELEVANT FIBERS"
   fib_init = fib_check(filecrop, filedir, right=keyword_set(right), check=checked[0:1], rund=runit[0:1])

   ;until we have thought on this a bit more, we are going to have to stop processing
   ;after working with the fibers and get_spatial_freq
   if total(stack_status_token) NE 10 then goto, theend

   ;Finding bispectra...
   print, ""  &  print, "FINDING BISPECTRA"

   gbp_file = strarr(10)
   for i = 0, 9 do gbp_file[i] = strcompress(filecrop+'_'+string(fix(i)),/remove_all)

   check = file_search(strcompress('data/'+filecrop+'_*_'+side+'_Bi.sav',/remove_all))

   if (n_elements(check) EQ 10)*checked[4] OR runit[4] then begin
      print, "Skipping get_best_params for: ", strmid(gbp_file[0], 0, strlen(gbp_file[0])-2)
   endif else begin
      get_best_params, gbp_file, fib_init, bin=rebin, right=keyword_set(right)
      if keyword_set(store) then begin
         ;spawn, 'rm data/'+filecrop+'*_'+side+'_dewarped.fits.gz '
         spawn, 'mv data/'+fib_init+'*fitted_smoothed_stack.fits.gz '+stor_loc
      endif
   endelse    



   ;finding closure phases
   print, "FINDING CLOSURE PHASES"

   check1 = file_search(strcompress('data/'+filecrop+'_9_'+side+'_Bi.sav',/remove_all))
   bundle = 100
   if logical_true(check1) then begin 
      for i = 0, 9 do begin
         check = file_search(strcompress('data/'+filecrop+'_'+string(i)+'_'+side+'_*b'+string(bundle)+'*CPmoy2.sav',/remove_all))

         if logical_true(check[0])*checked[5] or runit[5] then begin
            print, "Skipping get_CPmoy for: ", gbp_file[i]
         endif else get_CPmoy, gbp_file[i], 1, bundle, right=keyword_set(right), bin=rebin
      endfor

      if keyword_set(store) and ~runit[5] then begin
         spawn, 'mv data/'+filecrop+'*_'+side+'_Bi.sav '+stor_loc
         spawn, 'mv data/'+filecrop+'*_'+side+'_Spatial_freq.sav '+stor_loc
         spawn, 'mv data/'+filecrop+'*_'+side+'*_CP.sav '+stor_loc
      endif
   endif


   theend: tf= systime(/seconds)
   print, "Finished ", file
   print, "Ran in ", (tf-t0)/60., " minutes."
end





pro master_dark, file, dir
;+
;
;FUNCTION: MASTER_DARK, FILE, DIR
;
;PURPOSE:
;  To generate a composite dark file out of all the dark files with the
;  same exposure time as the target from the same night
;
;METHOD:
;  0] The file name, and exposure time of the image, is used to locate
;     the appropriate dark files.
;  1] All similar darks are stacked and median combined. This median
;     combined image is written to a fits file for use in correct_deform.
;  2] If darks can not be found on the night of the observation, darks
;     from the next closest night are used.
;
;INPUT:
;  file = The file for which a dark is needed.
;  dir = The directory in which the raw data file is located.
;
;OUTPUT:
;  (None) = a dark file is saved as 'master_dark_(exposure time).fits.gz
;
;NOTES:
;  -Rather than dividing the median-combined dark image by the median of
;   the whole image in this procedure to normalize the dark, the median
;   of the band region measured in correct_deform is used during
;   dark_correction. 
;
;-

  ;Checking to see if a dark has already been made...
  if strmid(dir,strlen(dir)-1,1) NE '/' then dir+='/'
  if strpos(dir, '~') NE -1 then begin
     spawn, 'echo ~',tilde
     dir = repstr(dir, '~',tilde)
  endif
  
  hdr = headfits(dir+file)
  if n_elements(hdr) EQ 1 then begin
     ;files are moved away if preprocessing finished
     hdr = headfits('data/'+repstr(file,'.fits.gz','_left_dewarped.fits.gz'))
     if n_elements(hdr) EQ 1 then stop, "Major problem"
  endif
  star = strcompress(strmid(file,14,strpos(file,'_',14)-14),/remove_all)
  date = strmid(file,5,8)


  dt = strcompress(string(fxpar(hdr,'SHUTTER')),/remove_all)
  filenaem = 'data/master_dark_'+date+'_'+dt+'.fits.gz'
  check = file_search(filenaem)




  ;If dark does not exist, generating dark...
  if check EQ "" then begin
     pos = strpos(file, '_', /reverse_search)
     filetag = strmid(file,5,9)
     
     dfiles = file_search(dir+'*'+filetag+'*dark*.fits.gz')
     zfiles = dfiles[where(strpos(dfiles,'_0.') NE -1)]
     dfile = ''
     for d = 0, n_elements(zfiles)-1 do begin
        hedr = headfits(zfiles[d])
        t = strcompress(string(fxpar(hedr,'SHUTTER')),/remove_all)
        if t EQ dt then begin
           pos = strlen(dir)
           init = strmid(zfiles[d], pos, 4)
           dfile = [dfile, init] 
        endif
     endfor
     if n_elements(dfile) EQ 1 then stop, "No appropriate darks?" $
     else dfile = dfile[1:*]
     ind0 = (where(strpos(dfiles,star) NE -1))[0]
     ndf = n_elements(dfiles)  &  list = sort(abs(findgen(ndf)-ind0))
     if n_elements(list) GT 20 then begin
        list=list[0:19]  &  dfiles = dfiles[list]
        init = strmid(dfiles,strpos(dfiles,'/',/reverse_search)+1,4)
        initu = init[uniq(init,sort(init))]  &  dfile = ''
        for d = 0, n_elements(initu)-1 do begin
           if total(strpos(dfiles,initu[d]) NE -1) mod 2 EQ 0 then begin
              dfile = [dfile,initu[d]]
           endif
        endfor
     endif

    
     print, "Processing dark frames and starting bad pixel search..."
     dstack = fltarr(512,512,2)
     for d = 0, n_elements(dfile)-1 do begin
        files = dfiles[where(strpos(dfiles, dfile[d]) NE -1)]
        condition = n_elements(files) EQ 4
                
        dark0 = mrdfits(files[0], /fscale, /silent) + 2.^15
        dstack = [[[dstack]],[[dark0]]]
        dark1 = mrdfits(files[1], /fscale, /silent) + 2.^15
        dstack = [[[dstack]],[[dark1]]]
        if condition then begin
           dark2 = mrdfits(files[2], /fscale, /silent) + 2.^15
           dstack = [[[dstack]],[[dark2]]]
           dark3 = mrdfits(files[3], /fscale, /silent) + 2.^15
           dstack = [[[dstack]],[[dark3]]]
        endif
        
        dark0f = median(dark0, dimension=3, /even)
        dark1f = median(dark1, dimension=3, /even)
        if condition then begin
           dark2f = median(dark2, dimension=3, /even)
           dark3f = median(dark3, dimension=3, /even)
        endif
        
       
       
    
        loc = where( abs(dark0f - median(dark0f)) gt 3.*stddev(dark0f) )
        badpixels_0 = intarr((size(dark0f))[1], (size(dark0f))[2])
        if max(loc) ne -1 then badpixels_0[loc] = 1
        
        loc = where( abs(dark1f - median(dark1f)) gt 3.*stddev(dark1f) )
        badpixels_1 = intarr((size(dark1f))[1], (size(dark1f))[2])
        if max(loc) ne -1 then badpixels_1[loc] = 1
        
        if condition then begin
           loc = where( abs(dark2f - median(dark2f)) gt 3.*stddev(dark2f) )
           badpixels_2 = intarr((size(dark2f))[1], (size(dark2f))[2])
           if max(loc) ne -1 then badpixels_2[loc] = 1
           
           loc = where( abs(dark3f - median(dark3f)) gt 3.*stddev(dark3f) )
           badpixels_3 = intarr((size(dark3f))[1], (size(dark3f))[2])
           if max(loc) ne -1 then badpixels_3[loc] = 1
        endif
        badpixels = badpixels_0 + badpixels_1 
        if condition then badpixels = badpixels + badpixels_2 + badpixels_3
        
        if condition then cut = 3 else cut = 1
        loc = where(badpixels lt cut) ; where only (1) 2 or fewer of the frames found a bad pixel
        badpixels[*,*] = 1
        badpixels[loc] = 0
        print, strcompress(round(total(badpixels))), " bad pixels were found."
        save, filename='data/'+dfile[d]+'_badpixelmap.sav', badpixels ;, /verbose
        
        ;REPORT MIN/MAX DARKS IN LOG
     endfor
     dstack = dstack[*,*,2:*]
     dark = median(dstack, dimension = 3)
     print, "Done with dark frame and bad pixel shenanigans."
     
     ;clean the header up before writing the fits file
     fxaddpar, hdr, 'STATUS','SOUND',' MASTER DARK',before='END'

     ndf = n_elements(dfile)  &  n_tags = ndf/3+logical_true(ndf mod 3)
     tags = 'DFILES' + strcompress(string(indgen(n_tags)),/remove_all)
     
     if n_tags GT 1 then begin
        for i = 0, n_tags-2 do begin
           fxaddpar, hdr, tags[i], strjoin(dfile[i*3:(i+1)*3-1],','), $
                     ' COMPONENT DARK INITS', before='STATUS'
        endfor
     endif
     fxaddpar, hdr, tags[-1], strjoin(dfile[(n_tags-1)*3:*],','), $
               ' COMPONENT DARK INITS', before='STATUS'

     hdr = hdr[where(hdr NE strjoin(strarr(80)+' '))]
     writefits, repstr(filenaem,'.gz',''), dark, hdr, /compress
  endif else print, "Appropriate dark located..."
end








function fib_check, file, dir, RIGHT= right, CHECKD= checkd, RUND= rund
;+
;FUNCTION: FIB_CHECK(FILE,DIR,[.../RIGHT,CHECKD='XX',RUND='XX')
;
;PURPOSE:
;  Processes all fiber files, if not already processed, for a single
;  target designated by file.
;
;INPUT:
;  file = The file for which fiber files must be processed.
;  dir = 
;  /right = 
;  /checkd = 
;  /rund = 
;
;OUTPUT:
;
;NOTES:
;-obtaining info for GBP until observation log info is incorporated
;- ;10=_0.fits.gz
;
;-
  ;Getting relevant target info from file name...
  if strmid(dir,strlen(dir)-1,1) NE '/' then dir+='/'
  if keyword_set(right) then side = 'right' else side = 'left'
  und1 = strpos(file, '_')
  init = fix(strmid(file, 0, und1))
  filetag = strmid(file, und1+1, strlen(file)-(und1+1)) 




  ;Looking for processed fiber files and obtaining init value...
  check = file_search('data/*_'+filetag+'_fibre*_'+side+$
                      '_fitted_smoothed_stack.fits.gz')
  inits = strmid(check, strlen('data/'), 4)  &  inits_i=fix(inits)
  diff_arr = abs(inits_i-init)
  fibstart = inits[(where(diff_arr EQ min(diff_arr)))[0]]





  ;If fibers have not been processed, dewarping and smoothing...
  if logical_true(check[0]) then begin
     print, "Fiber profiles have been processed already"      
  endif else begin
     ;Finds all the relevant fiber files, and gathers information on 
     ;...the corresponding sequence number
     print, "Processing fibers profiles..."
     fib_status_token = intarr(9)+1  &  n_fst = 9

     ;dir='data/'
     ;files =  file_search('data/*_'+filetag+'_fibre_0_'+side+'_dewarped.fits.gz')  
     files = file_search(dir+'*_'+filetag+'_fibre_0.fits.gz')  
     if files[0] EQ "" then goto, theend
     dlen = strlen(dir)  &  init_fib_u = strmid(files, dlen, 4)  
     N = n_elements(files)
     ;Init_fib contains all of the unique seq #s for the fiber files





     for i = 0, N-1 do begin
        ;Obtaining deformation profile...
        ;dlen = 5  &  val = strlen('_?_'+side+'_dewarped.fits.gz')
        ;filecrop = strmid(files[i], dlen, strlen(files[i])-val-dlen)
        filecrop = strmid(files[i], dlen, strlen(files[i])-dlen-10)
        print, "Working on ", filecrop, "..."

        check = file_search('data/'+filecrop+'*'+side+'*deform*')
        if (check NE "")*checkd[0] OR rund[0] then begin
           print, 'Skipping get_deform for: ',filecrop
        endif else begin
           get_deform, filecrop, dir, /fiber, right=keyword_set(right)
        endelse



        ;Checking to make sure a general dark has been generated...
        if ~rund[1] then master_dark, filecrop +'_0.fits.gz', dir




        ;Dewarping fiber profiles...
        for j = 0, n_fst-1 do begin
           k = strcompress(string(fix(j)),/remove_all)
           check = (file_search('data/'+filecrop+'_'+k+ $
                               '_'+side+'_dewarped.fits.gz')) NE ""
           check2 = (file_search('data/'+filecrop+'_'+k+ $
                               '_'+side+'_dewarped1.fits.gz')) NE ""

           
           chk = total(~fib_status_token)
           if (check OR check2)*checkd[1] OR rund[1] OR chk then begin
              print, "Skipping correct_deform for ", filecrop+'_'+k
           endif else begin
              correct_deform, filecrop, k, dir, right=keyword_set(right)
           endelse
        endfor



        
        ;Correcting the edges of the dewarped profiles...
        check = file_search('data/'+filecrop+'_*_'+ $
                            side+'_dewarped1.fits.gz')
        s = ~rund[1]
        case s*total((check NE "")*fib_status_token) of
           n_fst: edge_correct, filecrop, 'data/', $
                                right = keyword_set(right), /fiber

           0: print, "Skipping edge_correction for: ", filecrop

           else: begin
              if total(fib_status_token) EQ n_fst then begin 
                 edge_correct, filecrop, 'data/', /rescue, /fiber, $
                               right = keyword_set(right)
              endif else begin
                 print, "Not all files dewarped, can't correct edges."
              endelse
           end
        endcase



        
        ;Smoothing the fiber profiles for get_best_params...
        if n_elements(img_sm) NE 0 then dum = size(temporary(img_sm))
        if n_elements(img_nrm) NE 0 then dum = size(temporary(img_nrm))
        if total(fib_status_token) EQ n_fst then begin
           for j = 0, n_fst-1 do begin
              k = strcompress(string(fix(j)),/remove_all)
              chk = (file_search('data/'+filecrop+'_'+k+'_'+ $
                                side+'_fitted.fits.gz'))[0]
              if chk EQ "" then begin
                 procd_fib = get_fibers(filecrop+'_'+k, $
                                        right = keyword_set(right))
              endif else begin
                 print, "Fiber already smoothed."
                 fn = mrdfits(chk, 0, nhdr)
                 sfile = repstr(chk,'ted.','ted_smoothed.')
                 fs = mrdfits(sfile, 0, shdr)

                 procd_fib = {fs: fs, fn: fn, sh: shdr, nh: nhdr}
              endelse


              if size(procd_fib,/type) EQ 8 then begin
                 if n_elements(img_sm) EQ 0 then begin
                    siz = size(procd_fib.fs)
                    img_sm = fltarr(siz[1], siz[2], n_fst)
                    img_nrm = img_sm
                 endif
              
                 img_sm[*,*,j] = procd_fib.fs
                 img_nrm[*,*,j] = procd_fib.fn
              endif
           endfor
        endif




        ;Save image stacks as .fits files
        if n_elements(img_sm) NE 0 AND size(procd_fib,/type) EQ 8 then begin
           chk = total(reform(total(total(img_sm,2),1)) EQ 0)
           if chk EQ 0 then begin
              if total(img_sm LT 0) GT 0 then img_sm -= min(img_sm)
              if total(img_nrm LT 0) GT 0 then img_nrm -= min(img_nrm)

              filecrop = repstr(filecrop, 'fibre_', 'fibre')
              beg = 'data/'+filecrop+'_'+side+'_fitted_'
              fin = '_stack.fits'  & tag = ['smoothed','norm']
              smooth = beg+tag[0]+fin  &  norm = beg+tag[1]+fin
              

              shdr = procd_fib.sh  &  nhdr = procd_fib.nh
              fxaddpar, shdr, 'STATUS', $
                        fxpar(shdr, 'STATUS',comment=dum), dum+' STACK'
              fxaddpar, nhdr, 'STATUS', $
                        fxpar(nhdr, 'STATUS',comment=dum), dum+' STACK'
              writefits, smooth, img_sm, shdr, /compress
              writefits, norm, img_nrm, nhdr, /compress
           endif else begin
              print,"Could not stack fiber profiles for: ",filecrop
           endelse
        endif else print, "No good fiber profiles for: ", filecrop
     endfor

     print, "Fiber profiles have now been processed"
  endelse
  
  return, fibstart
  theend: print, "No relevant fibers found?"
  return, 'bumpkiss'
end



