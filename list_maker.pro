pro cheat, name
  files = file_search('data/*'+name+'*_5_*dewarped*')
  files = repstr(files, 'left_dewarped','')
  files = repstr(files, 'right_dewarped','')
  files = files[uniq(files)]
  files = files[where(strpos(files,'fibre') EQ -1)]
  
  openw,110, strcompress(name+'_filelist.txt',/remove_all)
  for i = 0, n_elements(files)-1 do begin
     printf, 110, files[i]
  endfor
  close, 110
end


pro list_maker, star, date = date, all=all
;+
;PURPOSE:
;This procedure allows me to make a file list for use in wrapper based
;on locating all the files 
;
;INPUT:
;star: the star to locate the files for
;
;OUTPUT:
;none, but saves a list as a star+filelist.txt file
;-

  run = 'Lick_July2012'
  if keyword_set(date) then begin
     files = file_search(strcompress('~/urap/rawdata/'+run+'/*/*_'+date+'_'+star+'_5.fits.gz',/remove_all))
     if keyword_set(all) then star = 'all'
     name = date+'_'+star
  endif else begin
     files = file_search(strcompress('~/urap/rawdata/'+run+'/*/*_*_'+star+'_5.fits.gz',/remove_all))
     if keyword_set(all) then star = 'all'
     name = star
  endelse
  ;finding all the data involving that star (only grabbing those w/ 5 because
  ;that will acct for all non fib/dark and wrapper doesn't need the number)

  files = files[where(strpos(files,'fibre') EQ -1)]
  openw,110, strcompress(name+'_filelist.txt',/remove_all)
  for i = 0, n_elements(files)-1 do begin
     printf, 110, files[i]
  endfor
  close, 110

  print, n_elements(files), " files located and written to: ", strcompress(name+'_filelist.txt',/remove_all)

  ;writing a properly formatted .txt file of all of the file names to run the wrapper on
end

pro list_maker_rb, star
;+
;PURPOSE:
;This procedure acts like list_maker, but instead looks for processed
;files to bin for a star
;
;INPUT:
;star: the star to look for processed files for
;
;OUTPUT:
;none, but saves a list as star+binlist.txt
;-

  files = file_search(strcompress('data/*'+star+'_dewarped.fits.gz',/remove_all))  
  files = strmid(files, 0, strlen(files[0])-5)
  ;gathering appropriate files and cutting them down to easy to work with size

  openw,110, strcompress(star+'_binlist.txt',/remove_all)
  for i = 0, n_elements(files)-1 do begin
     printf, 110, files[i]
  endfor
  close, 110

  print, n_elements(files), " files located and written to: ", strcompress(star+'_binlist.txt',/remove_all)
  ;writing a properly formatted .txt file of all of the file names to run the rebin_all on
end

