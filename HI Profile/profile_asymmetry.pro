function profile_asymmetry,velocity,flux,v0,vlow,vhigh,rms=rms,plot=plot,inspect=inspect

  if n_elements(plot) eq 0 then plot=0

  ;get rms if it's not already defined
  if n_elements(rms) eq 0 then begin
     sel=(velocity lt vlow) or (velocity gt vhigh)
     rms = robust_sigma(flux[sel])
  endif


  sel1=where(velocity ge vlow and velocity le v0,count1)
  sel2=where(velocity ge v0 and velocity le vhigh,count2)

  if (count1 ge 2) and (count2 ge 2) then begin

     ;ignoring channel sizes on these sums. Ok to first order.
     flow = total(flux[sel1])
     eflow = rms*sqrt(count1)
     fhigh = total(flux[sel2])
     efhigh = rms*sqrt(count2)

     ratio = flow/fhigh
     eratio = sqrt((eflow/fhigh)^2 + (flow*efhigh/fhigh^2)^2)

     asym=alog10(flow/fhigh)
     easym=1/alog(10)*eratio/ratio

;show plot if keyword is set
     if n_elements(plot) gt 0 then begin
        plot,velocity,flux,xtitle='velocity',ytitle='flux density [Jy]',/xsty,/ysty,xrange=[v0-1000,v0+1000]
        oplot,velocity,flux,psym=4,symsize=0.5
        oplot,[v0,v0],[-1000,1000],color=cgcolor('green'),thick=3
        oplot,velocity[sel1],flux[sel1],color=cgcolor('red')
        oplot,velocity[sel2],flux[sel2],color=cgcolor('blue')
        oplot,[-9e9,9e9],[1.5*rms,1.5*rms],linestyle=2
     endif

  endif else begin
     
     asym = -999
     easym = -999

  endelse

  if keyword_set(inspect) then stop

  return, [asym,easym]

end


pro asymmetry_himanga,id,asym,easym,v0=v0,path=path,plot=plot

  ;id should be SDSS ID


  asym=-99
  if 1-keyword_set(path) then path = './'
  
                                ;define par file
  par = path + id+'_par.sav'
  
                                ;define spectrum file
  spec = path + 'mangaHI-'+id+'.fits'
  
  restore,par
  s = mrdfits(spec,1,hdr)
  if n_elements(v0) eq 0 then v0 = sxpar(hdr,'OBJ_VEL')
  rms = par.statinfo.rms
  
  if n_elements(s.vhi) gt 1000 then begin
     print,'error in spectrum size. Skipping'
     return
  endif

  awvinds = [par.awvinfo.xmin,par.awvinfo.xmax] - round(150/4)
  if awvinds[0] lt 0 or awvinds[1] lt 0 then return

  vmin = min([s.vhi[awvinds[0]],s.vhi[awvinds[1]]])
  vmax = max([s.vhi[awvinds[0]],s.vhi[awvinds[1]]])

  out = profile_asymmetry(s.vhi,s.fhi,v0,vmin,vmax,rms=rms,plot=1)
  asym = out[0]
  easym = out[1]

end

function find_edges,velocity,flux,v0,rms=rms,threshold=threshold,maxcount=maxcount,reverse=reverse,plot=plot,center_offset_start=center_offset_start


  if n_elements(maxcount) eq 0 then maxcount = 3
  if n_elements(threshold) eq 0 then threshold=1.5
  if n_elements(center_offset_start) eq 0 then center_offset_start=50

  ;find index corresponding to v0
  minsep = min(abs(velocity-v0),ind0)

  ;do quick meausre of rms
  if n_elements(rms) eq 0 then begin
     sel=where((velocity lt (v0 - 300)) or (velocity gt (v0 + 300)))
;     rms = robust_sigma(flux[sel])
     rms = median(abs(flux[sel] - median(flux)))/0.67
     print,'calculated rms: ',rms
  endif

  if keyword_set(plot) then begin
     plot,flux
     oplot,[-9e9,9e9],threshold*[rms,rms],color=cgcolor('green'),linestyle=2
  endif



  edges = [-1,-1]
  for step=-1,1,2 do begin
     
     if keyword_set(reverse) then dir = -1 $
                                        else dir = 1

     ind=ind0
     if keyword_set(reverse) then begin
        ind = ind0+center_offset_start*step
        if keyword_set(plot) then oplot,[ind,ind],[-1000,1000],color=cgcolor('red')
     endif

     done=0
     counter=0
     while done eq 0 do begin
        ind=ind+dir*step

        if 1-keyword_set(reverse) and (flux[ind] lt threshold*rms) then begin
           counter += 1
        endif else if keyword_set(reverse) and (flux[ind] gt threshold*rms) then begin
           counter += 1
        endif else counter = 0

        if counter eq maxcount then begin
           print,'edge found'
           ind = ind-maxcount*step*dir ;go back to the first instance where we were below the threshold
           if keyword_set(plot) then begin
              oplot,[ind,ind],[-1000,1000],color=cgcolor('turquoise')
           endif

           if step eq -1 then edges[0] = ind $
           else edges[1] = ind
           done=1

        endif else if (1-keyword_set(reverse)) and (counter ne maxcount) and ((ind eq 0) or (ind eq n_elements(flux)-1)) then begin
           print,'Failed to find edge'
           if step eq -1 then edges[0] = ind $
           else edges[1] = ind
           done=1
        endif else if (keyword_set(reverse)) and (counter ne maxcount) and (ind eq ind0) then begin
           print,'failed to find edge'
           if step eq -1 then edges[0] = ind0+center_offset_start*step $
           else edges[1] = ind0+center_offset_start*step
           done=1

        endif
     endwhile
  endfor
  print,'edges = ',edges
  return, edges

end



pro asymmetry_alfalfa,id,asym,easym,v0=v0,path=path,plot=plot,center_offset_start = center_offset_start
    
  if n_elements(v0) eq 0 then begin
     print,'must specify v0'
     return
  endif

  spec = path+id+'.fits'
  
  s=mrdfits(spec,1,hdr)
  
  ;need to find vmin/vmax from the data
  edges = find_edges(s.vhelio,s.flux,v0,rms=rms,threshold=1.5,/reverse,/plot,maxcount=4,center_offset_start = center_offset_start)
  vmin = min(s[edges].vhelio)
  vmax = max(s[edges].vhelio)

  out = profile_asymmetry(s.vhelio,s.flux,v0,vmin,vmax,rms=rms,plot=1)
  asym=out[0]
  easym=out[1]
  print,vmin,vmax
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

use_optv=0

;;;now calculate asym for everything in bbrd, gbt and alfalfa
cat = mrdfits('/users/dstark/breakbrd/master/catalogs/bbrd_HI_catalog_v2.fits',1) ; HI catalog
t=mrdfits('/users/dstark/breakbrd/bbrd_tables_additions_092322.fits',1) ; main catalog

;;;;GBT first;;;;
skip=0
if skip eq 0 then begin
gcat = cat[where(cat.session ne 'ALFALFA')]
id_arr = gcat.sdss_name
asymarr = fltarr(n_elements(gcat))
easymarr = fltarr(n_elements(gcat))
alfalfa=intarr(n_elements(gcat))

for i=0,n_elements(gcat)-1 do begin
   ;find the observing name
   sel=where(strtrim(t.sdss_name,2) eq strtrim(gcat[i].sdss_name,2),count)
   if count gt 0 and gcat[i].logmhi gt 0 then begin
      print,'finding asymmetry for '+gcat[i].sdss_name

      ;determine which central velocity to use
      if use_optv eq 1 then use_vel = t[sel].z_new*2.998e5 $
      else use_vel = gcat[i].vhi

      asymmetry_himanga,t[sel].shortname,asym,easym,v0=use_vel,path='/users/dstark/breakbrd/master/reduction_files/',/plot
      asymarr[i] = asym
      easymarr[i] = easym
      print,'asym,easym: ',asym,easym

   endif
endfor
endif


;;;now ALFALFA

acat = cat[where(cat.session eq 'ALFALFA')]
;asymarr_alfalfa = fltarr(n_elements(gcat))
;easymarr_alfalfa = fltarr(n_elements(gcat))

for i=0,n_elements(acat)-1 do begin
   ;find the observing name
   sel=where(strtrim(t.sdss_name,2) eq strtrim(acat[i].sdss_name,2),count)
   if count gt 0 and acat[i].logmhi gt 0 then begin
      print,'finding asymmetry for '+acat[i].sdss_name
;      asymmetry_alfalfa,'a100_'+strjoin(strsplit(t[sel].sdss_name,/extract)),asym,easym,v0=t[sel].z_new*2.998e5,path='/users/dstark/breakbrd/alfalfa_crossmatch/spectra_official/',/plot,center_offset_start = round(acat[i].wf50/2/5.25*2)

      if use_optv eq 1 then use_vel = t[sel].z_new*2.998e5 $
      else use_vel = acat[i].vhi
      print,'using velocity: ',use_vel

      asymmetry_alfalfa,'a100_'+strjoin(strsplit(t[sel].sdss_name,/extract)),asym,easym,v0=use_vel,path='/users/dstark/breakbrd/alfalfa_crossmatch/spectra_official/',/plot,center_offset_start = round(acat[i].wf50/2/5.25*2)
      ;asymarr_alfalfa[i] = asym
      ;easymarr_alfalfa[i] = easym
      print,'asym,easym: ',asym,easym
      asymarr = [asymarr,asym]
      easymarr = [easymarr,easym]
      alfalfa = [alfalfa,1]
      id_arr = [id_arr,acat[i].sdss_name]

   endif
endfor

;put this info into a structure and save as a fits file
table = {id:id_arr[0],alfalfa:alfalfa[0],asym:asymarr[0],easym:easymarr[0]}
table = replicate(table,n_elements(id_arr))
table.id = id_arr
table.alfalfa = alfalfa
table.asym = asymarr
table.easym = easymarr

if use_optv eq 1 then outfile = 'breakbrd_asyms_optv0_02Dec2022.fits' $
                                else outfile = 'breakbrd_asyms_hiv0_02Dec2022.fits'

mwrfits,table,outfile,/create


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;now do the same things for the HI MaNGA survey;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

cat=mrdfits('/users/dstark/17AGBT012/master/catalogs/mangahi_dr3_111321.fits',1)
gcat = cat[where(cat.session ne 'ALFALFA')]
acat = cat[where(cat.session eq 'ALFALFA')]
drp=mrdfits('/users/dstark/17AGBT012/idl_routines/drpall-v3_0_1.fits',1)

skip=0
if skip eq 0 then begin
;gbt
id_arr = gcat.plateifu
asymarr = fltarr(n_elements(gcat))
easymarr = fltarr(n_elements(gcat))
alfalfa=intarr(n_elements(gcat))

for i=0,n_elements(gcat)-1 do begin
   ;find the observing name
   sel=where(drp.plateifu eq strtrim(gcat[i].plateifu,2),count)
   if count gt 0 then begin
   print,sel
   ;vhel=drp[sel].z*2.998e5
   vhel=gcat[i].vhi

   if use_optv eq 1 then use_vel = drp[sel].z*2.998e5 $
   else use_vel = gcat[i].vhi

   if (gcat[i].logmhi gt 0) and file_test('/users/dstark/17AGBT012/master/reduction_files/'+strtrim(gcat[i].plateifu,2)+'_par.sav') then begin
      print,'finding asymmetry for '+gcat[i].plateifu
      asymmetry_himanga,gcat[i].plateifu,asym,easym,v0=use_vel,path='/users/dstark/17AGBT012/master/reduction_files/',/plot
      asymarr[i] = asym
      easymarr[i] = easym
      print,'asym,easym: ',asym,easym

   endif
endif

endfor
endif

if use_optv eq 1 then outfile = 'himanga_asyms_optv0_02Dec2022.fits' $
else outfile = 'himanga_asyms_hiv0_02Dec2022.fits'

;put this info into a structure and save as a fits file
table = {id:id_arr[0],alfalfa:alfalfa[0],asym:asymarr[0],easym:easymarr[0]}
table = replicate(table,n_elements(id_arr))
table.id = id_arr
table.alfalfa = alfalfa
table.asym = asymarr
table.easym = easymarr

mwrfits,table,outfile,/create

;alfalfa
;id_arr = gcat.plateifu
;asymarr = fltarr(n_elements(acat))
;easymarr = fltarr(n_elements(acat))
;alfalfa=intarr(n_elements(acat))

;for i=0,n_elements(acat)-1 do begin
   ;find the observing name
   ;sel=where(strtrim(t.sdss_name,2) eq strtrim(acat[i].sdss_name,2),count)
 ;  if  acat[i].logmhi gt 0 then begin
  ;    print,'finding asymmetry for '+acat[i].plateifu
;      asymmetry_alfalfa,'a100_'+strjoin(strsplit(t[sel].sdss_name,/extract)),asym,easym,v0=t[sel].z_new*2.998e5,path='/users/dstark/breakbrd/alfalfa_crossmatch/spectra_official/',/plot,center_offset_start = round(acat[i].wf50/2/5.25*2)
   ;   asymmetry_alfalfa,acat[i].plateifu,asym,easym,v0=acat[i].vhi,path='/users/dstark/17AGBT012/master/reduction_files/',/plot,center_offset_start = round(acat[i].wf50/2/5.25*2)
      ;asymarr_alfalfa[i] = asym
      ;easymarr_alfalfa[i] = easym
    ;  print,'asym,easym: ',asym,easym
      
   ;endif
;endfor

end
