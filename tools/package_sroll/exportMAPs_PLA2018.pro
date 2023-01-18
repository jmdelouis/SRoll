;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
; procedures to export various MAP products to EFDD FITS files.  Params that
; are likely to change from one release to the next, like group names or 
; source directory, are at the beginning of each script.  
; 
; Currently handled:         module name
; - [detector maps]
; - frequency map            freqmap
; - Gal and ptSrc masks      skymask
; - Zodi correction maps     zodicorr
; - cmb maps: nils/smica/sevem   cmb method
; - CO type 1/2/3            co1/2/3
; - Dust opacity map         dustops
; - Sky power spectra and covar matrices     powspec
; - Mask for sky power spectra               psmask
; - CMB power spectrum and cov mat           cmbspec
;
;   Loops
; - Loop for freq maps       freqmaps
; - Loop over detsets        pspectra
; 
; and finally a module to check the EFDD file headers and structure
; - check, 'filename' [, /quiet to not print header]
;
;-----------------------------------------------------------------------------

function bplsuffix, suff
   case suff of
       'full': tag = 'GHz_dustleak_ground_'
       'year-1': tag = 'GHz_dustleak_ground_year1_'
       'year-2': tag = 'GHz_dustleak_ground_year2_'
       'halfmission-1': tag = 'GHz_dustleak_ground_hm1_'
       'halfmission-2': tag = 'GHz_dustleak_ground_hm2_'
       'survey-1': tag = 'GHz_dustleak_ground_survey1_'
       'survey-2': tag = 'GHz_dustleak_ground_survey2_'
       'survey-3': tag = 'GHz_dustleak_ground_survey3_'
       'survey-4': tag = 'GHz_dustleak_ground_survey4_'
       'survey-5': tag = 'GHz_dustleak_ground_survey5_'
   endcase
   return, tag
end

function r2tag, suff
    case suff of   ; prepare tags for round-2 files
        'full': tag = '_full'
        'full-ringhalf-1': tag = '_full_hr1'
        'full-ringhalf-2': tag = '_full_hr2'
        'full-halfmission-1': tag = '_full_hm1'
        'full-halfmission-2': tag = '_full_hm2'
        'year-1': tag = '_yr1'
        'year-2': tag = '_yr2'
        'year-1-ringhalf-1': tag = '_yr1_hr1'
        'year-1-ringhalf-2': tag = '_yr1_hr2'
        'year-2-ringhalf-1': tag = '_yr2_hr1'
        'year-2-ringhalf-2': tag = '_yr2_hr2'
        'halfmission-1': tag = '_hm1'
        'halfmission-2': tag = '_hm2'
        'halfmission-1-ringhalf-1': tag = '_hm1_hr1'
        'halfmission-1-ringhalf-2': tag = '_hm1_hr2'
        'halfmission-2-ringhalf-1': tag = '_hm2_hr1'
        'halfmission-2-ringhalf-2': tag = '_hm2_hr2'
    endcase 
    return, tag
end 

function fsl_corr, freq
   case strmid(freq,0,3) of
       '100': corr = 1.00087d
       '143': corr = 1.00046d
       '217': corr = 1.00043d
       else:  corr = 1.00d
   endcase 
   return, corr
end 

; apply fslcorrection only to valid pixels 
function submap, name, bpcorr, fslcorr, com

   map = float(pioread(name , com=com))
   bad = where(map eq -1.63750e+30, compl=good, nbad, ncompl=ngood)
   if ngood gt 0 then begin
       if n_elements(bpcorr) gt 10 then map[good] = (map[good] - bpcorr[good]) * fslcorr else map[good] = map[good] * fslcorr
   endif else begin
       if n_elements(bpcorr) gt 10 then map = (map - bpcorr) * fslcorr else map = map * fslcorr
   endelse 
   
   return, map
end 

; compute difference of valid pixels only
function diffmap, m1, m2

   bpval = -1.63750e+30
   bad = where(map eq bpval or m2 eq bpval, compl=good, nbad, ncompl=ngood)
   diff = m1 - m2
   if ngood gt 0 then diff[bad]  = bpval
   
   return, map
end 


function convertcov,name
;====================================================================================
; COMPUTE COVARIANCE AS PLANCK 2015 FROM RD12_RC4
;====================================================================================
;cov=convertcov('/redtruck/delouis/PROD_RD12_RC4/100ds1_ful.all_ful.RD12_RC4.COV.fits')
;cov is a 6 column maps:II,IQ,IU,QQ,QU,UU

read_fits_map,name,m2

res=fltarr(12l*2048l*2048l,3,3)
for i=0l,12l*2048l*2048l-1l do begin
    if i mod (2048l*2048l) eq (2048l*2048l-1) then print,string((i+1)/(12l*2048.*2048.)*100,format='(F7.2)'),'% Done'
    if m2[i,0] gt 0 then res(i,*,*)=invert([[m2[i,0],m2[i,1],m2[i,2]], $
                                            [m2[i,1],m2[i,3],m2[i,5]], $
                                            [m2[i,2],m2[i,5],m2[i,4]]],/double)
end

p1=[0,1,2,1,2,2]
p2=[0,0,0,1,1,2]
; XTRACT 6 Columns
; 6 colums :II,IQ,IU,QQ,QU,UU]
ref=fltarr(12l*2048l*2048l,6)
for i=0,5 do ref(*,i)=res(*,p1[i],p2[i])

return,ref
end



;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; export RC4 (PR3) products
; /test to only print input map names and output file names (DO THIS FIRST)
;-----------------------------------------------------------------------------

pro pr3maps, test=test, stop=stop, verb=verb

   if defined(test) then test=1 else test=0
   if defined(verb) then verb=verb  else verb=0

   freqs = ['100', '143', '217', '353', '545', '857']
   swbs  = ['143-5','143-6','143-7',  '217-1','217-2','217-3','217-4',  '353-1','353-2','353-7','353-8', $
            '545-1','545-2','545-4',  '857-1','857-2','857-3','857-4']

   nt=0 
   ; this is used to name output maps
   cover = ['full', 'halfmission-1','halfmission-2', 'year-1','year-2', $
            'survey-1','survey-2','survey-3','survey-4','survey-5']
   ; for PR3:
   cover = ['full', 'halfmission-1','halfmission-2']

   ;---------------------------------
   ; full channel maps
   ;---------------------------------
   nn=0

;   for c=2,2 do for f=0,0 do begin   ;; for testing
   for c=0,2 do for f=0,3 do begin   ;; for PR3
       rc4map, freqs[f], "GHz", cover[c], test=test, verb=verb, stop=stop    &     nn=nn+1
   endfor 
   print, " >> ... exported ", strtrim(nn,2), " full channel maps"   &  nt=nt+nn

   ;---------------------------------
   ; some PSB only maps at 353 
   ;---------------------------------
   print, ""   &    nn=0
;   for c=2,2 do for f=3,3 do begin  ;; for testing
   for f=3,3 do for c=0,2 do begin  ;; for PR3
       rc4map, freqs[f], "psb", cover[c], fit=".psb_ful", test=test, verb=verb   &    nn=nn+1
   endfor
   print, " >> ... exported ",strtrim(nn,2), " selected PSB only maps at 353"   &  nt=nt+nn
   
   ;---------------------------------
   ; odd-even maps
   ;---------------------------------
   print, ""   &    nn=0
   types = ['GHz', 'ds1','ds2','ds3','ds4','ds5','ds6']
   for c=0,0 do for f=0,3 do for t=0,0 do begin               ; for PR3: no detsets
       rc4map, freqs[f], types[t], 'evenring', test=test, verb=verb    &     nn=nn+1
       rc4map, freqs[f], types[t], 'oddring',  test=test, verb=verb    &     nn=nn+1
   endfor 
   rc4map, '353', 'psb', 'evenring', fit=".psb_ful", test=test, verb=verb    &     nn=nn+1
   rc4map, '353', 'psb', 'oddring',  fit=".psb_ful", test=test, verb=verb    &     nn=nn+1
   print, " >> ... exported ", strtrim(nn,2), " odd-even ring maps"   &  nt=nt+nn

   ;---------------------------------
   ; finish up
   ;---------------------------------
   print, "   ==>> exported ", strtrim(nt,2), " maps total"

   return
end 

;-----------------------------------------------------------------------------
; export single RC4 frequency map
; freq: 3-digit freq. or  detname
; type: GHz or ds3 for detector selection (not used for single bolo)
; suff: full or survey-n ... is coverage
; fit:  selection used for destripe fit (def: .all_full, .all_hr1,2 for halfring)
; oe:   odd/eve/all
;-----------------------------------------------------------------------------

pro rc4map, freq, type, suff, fit=fit, oe=oe, $
            stop=stop, test=test, verb=verb

   if defined(test) then test=test  else test=0
   if defined(verb) then verb=verb  else verb=0
   if defined(oe)   then oe  =oe    else oe='full'
   if defined(stop) then stop=stop     else stop=0

   if not defined(fit) then fit = ".all_ful" 
   if (strmid(fit,5,2) eq 'hr') then hr=1 else hr=0

   relnum = "_R3.01"   &   procver = "RD12_RC4"
   Nside = 2048LL  &   Npix = 12*Nside*Nside
   indir = "/redtruck/delouis/PROD_RD12_RC4/"
   hitdi = "/mnt/gpfs3/opsman/Releases/HFI_PR3/"   ;/redtruck/SimuData/DEC16v1_fits/"
   outdir = "/redtruck/opsman/Releases/HFI_PR3_beta/"
   hpinv  = -1.63750e+30   ; Healpix invalid

   ; ttag is appended to freq to build input filename
   if type eq "GHz" then ttag = "GHz" else ttag=type

   case suff of
       "":         tag=ttag
       "full":     tag = ttag+"_ful"
       'year-1':   tag = ttag+'_yr1'
       'year-2':   tag = ttag+'_yr2'
       'survey-1': tag = ttag+'_s1'
       'survey-2': tag = ttag+'_s2'
       'survey-3': tag = ttag+'_s3'
       'survey-4': tag = ttag+'_s4'
       'survey-5': tag = ttag+'_s5'
       'halfmission-1': tag = ttag+'_hm1'
       'halfmission-2': tag = ttag+'_hm2'
       'oddring':     tag = ttag+'_odd'
       'evenring':    tag = ttag+'_even'
       else: begin
           print, " ERROR: suffix "+suff+" unknown ... quitting"
           return
       end 
   endcase
   osuff = suff
   if suff eq "oddring" then suff="full-oddring"
   if suff eq "evenring" then suff="full-evenring"

   case freq of
       '100': begin & bw =  '33' & rf = 100.89 & end
       '143': begin & bw =  '46' & rf = 142.88 & end
       '217': begin & bw =  '65' & rf = 221.16 & end
       '353': begin & bw = '102' & rf = 357.5  & end
       '545': begin & bw = '171' & rf = 555.2  & end
       '857': begin & bw = '245' & rf = 866.8  & end
   endcase

   if fix(freq) le 400 then pol=1 else pol=0
   if strmid(freq,3,1) eq "-" and strlen(freq) eq 5 then pol=0  ; single SWB bolo map

   ; inputs
   if verb then print, " >> Read data for freq = "+freq
   if pol then begin
       signame = freq+tag+fit+'.RD12_RC4.P.fits'
       covname = freq+tag+fit+'.RD12_RC4.COV.fits'
   endif else begin
       signame = freq+tag+fit+'.RD12_RC4.I.fits'
       covname = freq+tag+fit+'.RD12_RC4.II.fits'
   endelse 
   print, suff
   if suff eq "full" then begin
       ;hitname = "RD12RC4_HIT_"+freq+"ghz_full.fits"
       hitname = 'HFI_SkyMap_'+ freq +'_'+strtrim(Nside,2)+'_R3.00_'+suff+'.fits'
   endif else begin
       ;hitname = "RD12RC4_HIT_"+freq+"ghz_"+strmid(tag,4,9)+".fits"
       hitname = 'HFI_SkyMap_'+ freq +'_'+strtrim(Nside,2)+'_R3.00_'+suff+'.fits'
   endelse 

   if verb then begin
       print, "    - sigfile: "+signame
       print, "    - hitfile: "+hitname
       print, "    - covfile: "+covname
   endif 
   if not test then begin
       sig = mrdfits(indir+signame, 1, hdsig, /sil)
       hit = mrdfits(hitdi+hitname, 1, hdhit, /sil)
       if pol then $
         cov = convertcov(indir+covname) $
       else $
         cov = mrdfits(indir+covname, 1, hdcov, /sil)
   endif else begin
       spawn, "\ls "+indir+signame, outs
       spawn, "\ls "+hitdi+hitname, outh
       spawn, "\ls "+indir+covname, outc
;       print, outs+ ";  "+ outh+ ";  "+ outc
   endelse 

   ; outputs: modify ttag for output names
   if ttag eq "GHz" then ttag=""
   if ttag eq "psb" then ttag="-psb"
   name = 'HFI_SkyMap_'+ freq +ttag+'_'+strtrim(Nside,2)+relnum
;   name = 'HFI_SkyMap_'+ freq +'_'+strtrim(Nside,2)+relnum
   com = 'begin=0;end='+strtrim(Npix-1,2)
   if not hr then outfile = name +'_'+suff+'.fits' $
   else outfile = name +'_'+suff+'-ringhalf-'+strmid(fit,7,1)+'.fits' 
   
   if not test then begin
       if pol then $
         map = replicate({i_stokes:0.0, q_stokes: 0.0, u_stokes: 0.0, hits: 0L, $
           ii_cov: 0.0, iq_cov: 0.0, iu_cov: 0.0,  qq_cov: 0.0, qu_cov: 0.0, uu_cov: 0.0}, Npix) $
       else $
         map = replicate({ i_stokes: 0.0,  hits: 0L, ii_cov: 0.0 }, Npix) 
   
       if verb then print, " >> Write structure and fill it in"
       map.i_stokes = reorder(reform(sig.(0), npix,1),  /r2n)
       if pol then begin
           map.q_stokes = reorder(reform(sig.q_polarisation, npix,1),  /r2n)
           map.u_stokes = reorder(reform(sig.u_polarisation, npix,1),  /r2n)
      ;     map.hits     = reorder(reform(hit.i_stokes, npix,1),  /r2n)
           map.hits     = hit.hits
        ;   map.ii_cov   = reorder(reform(cov.II, npix,1), /r2n)  ;;   &  nz=where(map.ii_cov gt 0)  &  map[nz].ii_cov = 1./map[nz].ii_cov
        ;   map.iq_cov   = reorder(reform(cov.IQ, npix,1), /r2n)  ;;   &  nz=where(map.iq_cov gt 0)  &  map[nz].iq_cov = 1./map[nz].iq_cov
        ;   map.iu_cov   = reorder(reform(cov.IU, npix,1), /r2n)  ;;   &  nz=where(map.iu_cov gt 0)  &  map[nz].iu_cov = 1./map[nz].iu_cov
        ;   map.qq_cov   = reorder(reform(cov.QQ, npix,1), /r2n)  ;;   &  nz=where(map.qq_cov gt 0)  &  map[nz].qq_cov = 1./map[nz].qq_cov
        ;   map.qu_cov   = reorder(reform(cov.QU, npix,1), /r2n)  ;;   &  nz=where(map.qu_cov gt 0)  &  map[nz].qu_cov = 1./map[nz].qu_cov
        ;   map.uu_cov   = reorder(reform(cov.UU, npix,1), /r2n)  ;;   &  nz=where(map.uu_cov gt 0)  &  map[nz].uu_cov = 1./map[nz].uu_cov
           map.ii_cov   = reorder(reform(cov[*,0], npix,1), /r2n)  
           map.iq_cov   = reorder(reform(cov[*,1], npix,1), /r2n)  
           map.iu_cov   = reorder(reform(cov[*,2], npix,1), /r2n)  
           map.qq_cov   = reorder(reform(cov[*,3], npix,1), /r2n)  
           map.qu_cov   = reorder(reform(cov[*,4], npix,1), /r2n)  
           map.uu_cov   = reorder(reform(cov[*,5], npix,1), /r2n)  
       endif else begin                                    
      ;     map.hits     = reorder(reform(hit.i_stokes, npix,1),  /r2n)
           map.hits     = hit.hits
           map.ii_cov   = reorder(reform(cov.(0), npix,1),/r2n)  ;;&  nz=where(map.ii_cov gt 0)  &  map[nz].ii_cov = 1./map[nz].ii_cov
       endelse 
       ; replace inpainted values by HPinvalid
       zz = where(map.ii_cov eq 0, nsel)  
       if nsel gt 0 then begin
           print, " >> replace ", strtrim(nsel, 2), " inpainted values"
           map[zz].i_stokes = hpinv  
           if pol then begin 
               map[zz].q_stokes = hpinv  &  map[zz].u_stokes = hpinv
           endif
       endif 
   endif 

   if stop then begin
;       help, /str, map 
       print, form='(" >> Exporting ", a-36, a-39, a-31, "  ==>  ", a0)', signame, covname, hitname, outfile
       print, " - mean signal: ", mean(map.I_stokes), "; input signal units: ", sxpar(hdsig, "TUNIT1")
       print, " - mean covar:  ", mean(map.II_cov), ";  input covar units: ", sxpar(hdcov, "TUNIT1")
;       stop 
   endif

   if not test then begin
       if verb then print, " >> Write header and export data"
       fxbhmake, hdr, Npix,       'FREQ-MAP',' Extension name', /init, /date, /extver
       sxdelpar, hdr, 'EXTNAME'
       sxaddpar, hdr, 'COMMENT',  '  ' , after='TFORM3'
       sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
       sxaddpar, hdr, 'COMMENT',  '  '
       
       if fix(freq) lt 400 then begin
           unitI= "Kcmb"   &   unitII ="Kcmb^2"
       endif else begin
           unitI= "MJy/sr"   &   unitII ="(Mjy/sr)^2"
       endelse 
       
       sxaddpar, hdr, 'TUNIT1',   unitI
       if pol then begin
           sxaddpar, hdr, 'TUNIT2',   unitI 
           sxaddpar, hdr, 'TUNIT3',   unitI 
           sxaddpar, hdr, 'TUNIT4',   '   '
           sxaddpar, hdr, 'TUNIT5',   unitII
           sxaddpar, hdr, 'TUNIT6',   unitII
           sxaddpar, hdr, 'TUNIT7',   unitII
           sxaddpar, hdr, 'TUNIT8',   unitII
           sxaddpar, hdr, 'TUNIT9',   unitII
           sxaddpar, hdr, 'TUNIT10',  unitII
           sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT10'
       endif else begin
           sxaddpar, hdr, 'TUNIT2',   '   '
           sxaddpar, hdr, 'TUNIT3',   unitII 
           sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT3'
       endelse 
       
       sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
       sxaddpar, hdr, 'COMMENT',  '  '
       sxaddpar, hdr, 'EXTNAME',  'FREQ-MAP',    ' Extension name'
       sxaddpar, hdr, 'PIXTYPE',  'HEALPIX'
       sxaddpar, hdr, 'POLCCONV', 'COSMO',       ' Polarization convention'
       sxaddpar, hdr, 'COORDSYS', 'GALACTIC',    ' Coordinate system'
       sxaddpar, hdr, 'ORDERING', 'NESTED',      ' Healpix ordering'
       sxaddpar, hdr, 'NSIDE',     Nside,        ' Healpix Nside'
       sxaddpar, hdr, 'FIRSTPIX',  0,            ' First pixel # (0 based)'
       sxaddpar, hdr, 'LASTPIX',   Npix-1,       ' Last pixel # (0 based)'
       sxaddpar, hdr, 'FILENAME',  outfile,      ' FITS filename'
       sxaddpar, hdr, 'BAD_DATA',  HPinv,        ' HEALPIX bad pixel value'
       sxaddpar, hdr, 'FREQ',     freq,          ' reference frequency'
       sxaddpar, hdr, 'PROCVER',  procver,       ' Product version'
       ; ---- addt'l keywords for Aladin etc. requested 6.jun ------
       sxaddpar, hdr, 'UNITFREQ', 'GHz ',        ' frequency units'
       sxaddpar, hdr, 'BNDCTR',   freq,          ' band center, same as FREQ'           ; string - this is the band "label"
       sxaddpar, hdr, 'RESTFRQ',  rf,            ' effective frequency'                 ; number: v_cen from Locke's table
       sxaddpar, hdr, 'BNDWID',   bw,            ' effective bandwidth (approximate)'   ; number: BW from Locke's table
       ; ---- END addt'l keywords for Aladin etc. requested 6.jun ------
       sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
       sxaddpar, hdr, 'COMMENT',  'Further details in the Planck Legacy Archive and Explanatory Supplement '
       sxaddpar, hdr, 'COMMENT',  'http://www.cosmos.esa.int/web/planck/pla/'
       sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
   endif

;   print, ' >> Exporting '+outdir+outfile
   print, form='(" >> Exporting ", a-36, a-39, a-31, "  ==>  ", a0)', signame, covname, hitname, outfile
   print, "       - with " + covname + " and " + hitname
   if not test then begin 
       mwrfits, map, outdir+outfile, hdr, /create
   ;    print, " - mean signal: ", mean(map.I_stokes), "; input/output signal units: ", sxpar(hdsig, "TUNIT1"), " / ", unitI
   ;    print, " - mean covar:  ", mean(map.II_cov), ";  input/output covar units: ", sxpar(hdcov, "TUNIT1"), " / ", unitII
   ;    print, " - rest freq:  ", rf, "; bandwidth: ", bw
   endif 
   if verb then print, " ----------------------------------------------------------------------------- "
   map = 0    ; recover memory allocation
   
   if stop then stop
   return
end

pro viewall

  spawn, "ls -1 /mnt/gpfs3/opsman/Releases/HFI_PR3_beta/*.fits", out
  planck_parchment, 1
  pxs=1200  

   for k=0,n_elements(out)-1 do begin
       mm = mrdfits(out[k], 1, hd)
     ;  tit = sxpar(hd, "PROCVER")+" "+sxpar(hd,"FREQ")+" GHz"  
       pos = strpos(out[k], "HFI_Sky")  &  tit= strmid(out[k],  pos,99)
       pos = strpos(tit, ".fits")    &  tit= strmid(tit,  0, pos)
       name = tit+".ps"
       mollview, col=253, pxs=1200, /hist,/nested, mm.(0), tit=tit, png=tit+".png"
       mollview, col=253, pxs=1200, /hist,/nested, mm.(0), tit=tit, ps=tit+".ps"
       spawn, "ps2pdf "+tit
   end
   return
end

pro rc4bias, test=test, verb=verb

   if defined(test) then test=test  else test=0
   if defined(verb) then verb=verb  else verb=0

   freqs=['100', '217','353']
   types=['', '_vs_NOISE_LOG10', '_vs_SIGNAL_LOG10']
   tags =['-nominal_', '-noiseRatio_','-signalRatio_']
   indir='/redtruck/delouis/PROD_RD12_RC4/CHECK/'
   relnum = "_R3.00"   &   procver = "RD12_RC4"
   Nside = 2048LL  &   Npix = 12*Nside*Nside
   outdir = "/redtruck/opsman/Releases/HFI_PR3_beta/"
   hpinv  = -1.63750e+30        ; Healpix invalid
   
   for i=0,2 do begin
       freq= freqs[i]
       for j=0,2 do begin
           type=types[j] & tag=tags[j]
           infile=indir+freq+'_CO_SUBPIX'+type+'.fits'
           outfile='HFI_BiasMap-CO'+tag+freq+'_2048_R3.00-full.fits'
           
           spawn, "\ls "+infile, string
           print, "  ==>  "+ outfile
           
           if not test then begin
               sig = mrdfits(infile, 1, hd)
               map = replicate({i_stokes:0.0, q_stokes: 0.0, u_stokes: 0.0}, Npix)
            
               map.i_stokes = reorder(reform(sig.temperature, npix,1),  /r2n)
               map.q_stokes = reorder(reform(sig.q_polarisation, npix,1),  /r2n)
               map.u_stokes = reorder(reform(sig.u_polarisation, npix,1),  /r2n)
               
               fxbhmake, hdr, Npix,       'BIAS-MAP',' Extension name', /init, /date, /extver
               sxaddpar, hdr, 'COMMENT',  '  ' , after='TFORM3'
               sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
               sxaddpar, hdr, 'COMMENT',  '  '
               sxaddpar, hdr, 'TUNIT1',   ''
               sxaddpar, hdr, 'TUNIT2',   ''
               sxaddpar, hdr, 'TUNIT3',   ''
               sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
               sxaddpar, hdr, 'COMMENT',  '  '
               sxaddpar, hdr, 'EXTNAME',  'BIAS-MAP',    ' Extension name'
               sxaddpar, hdr, 'PIXTYPE',  'HEALPIX'
               sxaddpar, hdr, 'POLCCONV', 'COSMO',       ' Polarization convention'
               sxaddpar, hdr, 'COORDSYS', 'GALACTIC',    ' Coordinate system'
               sxaddpar, hdr, 'ORDERING', 'NESTED',      ' Healpix ordering'
               sxaddpar, hdr, 'NSIDE',     Nside,        ' Healpix Nside'
               sxaddpar, hdr, 'FIRSTPIX',  0,            ' First pixel # (0 based)'
               sxaddpar, hdr, 'LASTPIX',   Npix-1,       ' Last pixel # (0 based)'
               sxaddpar, hdr, 'FILENAME',  outfile,      ' FITS filename'
               sxaddpar, hdr, 'BAD_DATA',  HPinv,        ' HEALPIX bad pixel value'
               sxaddpar, hdr, 'FREQ',     freq,          ' reference frequency'
               sxaddpar, hdr, 'PROCVER',  procver,       ' Product version'
               ; ---- addt'l keywords for Aladin etc. requested 6.jun ------
               ; sxaddpar, hdr, 'UNITFREQ', 'GHz ',        ' frequency units'
               ; sxaddpar, hdr, 'BNDCTR',   freq,          ' band center, same as FREQ'           ; string - this is the band "label"
               ; sxaddpar, hdr, 'RESTFRQ',  rf,            ' effective frequency'                 ; number: v_cen from Locke's table
               ; sxaddpar, hdr, 'BNDWID',   bw,            ' effective bandwidth (approximate)'   ; number: BW from Locke's table
               ; ---- END addt'l keywords for Aladin etc. requested 6.jun ------
               sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------', after='PROCVER'
               sxaddpar, hdr, 'COMMENT',  'Further details in the Planck Legacy Archive and Explanatory Supplement '
               sxaddpar, hdr, 'COMMENT',  'http://www.cosmos.esa.int/web/planck/pla/'
               sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
               
               mwrfits, map, outdir+outfile, hdr, /create
               map=0
            endif
        end 
    end
    
    return
end


;NB to delete IDL display windows produced by mollview:  for k=32,97 do wdelete, k

;;- ;-----------------------------------------------------------------------------
;;- ; HPRs   for q=0,49 do rc4hpr, vbolos[q]
;;- ;-----------------------------------------------------------------------------
;;- pro rc4hpr, bolo, suff, stop=stop, test=test
;;- 
;;-    relnum = "_R3.00"
;;-    procver = "RD12_RC4"
;;- 
;;-    ingrp  = '/data/dmc/MISS03/DATA/PBR_JMD/'
;;-    outdir = '/redtruck/opsman/Releases/HFI_PR3_beta/'  
;;- ;   outdir=""
;;- 
;;-    stags = ['signal','bandpass','dipole','hit','phase','phi','theta','psi']
;;-    suffs = ['all_full','band_all_full','dipole_all_full','hit','ph','ptg','ptg_TUPLE_1','ptg_TUPLE_2']
;;-    if strmid(bolo,0,3) eq '353' and strlen(bolo) eq 6 then $
;;-       suffs = ['psb_full','band_psb_full','dipole_psb_full','hit','ph','ptg','ptg_TUPLE_1','ptg_TUPLE_2']
;;- 
;;-    if n_elements(bolo) eq 0 then bolo='217-1'   ; for testing
;;-    Nside = 2048LL  &   Npix = 12*Nside*Nside
;;- 
;;-    for j=0,7 do begin  ; loop over 8 suffs
;;-        objname = bolo+'_REP6_'+suffs[j]
;;-        com="begin=6639360;end=720674863"
;;-        outfile='HFI_HPR_'+bolo+'-'+stags[j]+relnum+'.fits'
;;-    
;;-        if keyword_set(test) then begin
;;-            print, objname
;;-            spawn, "dmctool getobjectinfo "+ingrp+objname+ " | grep BeginIndex"
;;-            ima = intarr(10,10)
;;-        endif else begin
;;-            print, " >> begin reading "+objname
;;-            ; reform to image of 27664 pixels (columns) x 25811 rings (rows)
;;-            ima = reform(pioread(ingrp+objname, com=com), 27664, 25811, /overwrite)
;;-            print, "   ... done; now export to fits ..."
;;-        endelse 
;;- ;stop
;;-        mkhdr, hdr, ima, /IMAGE
;;-        sxaddpar, hdr, 'BEGRING',        240,     ' first ring'    
;;-        sxaddpar, hdr, 'ENDRING',      26050,     ' last ring'    
;;-        sxaddpar, hdr, 'PIXTYPE',  'HEALPIX',     ' pixeling scheme'
;;-        sxaddpar, hdr, 'NSIDE',       Nside,      ' Healpix Nside'
;;-        sxaddpar, hdr, 'BAD_DATA',  -1.63750E+30, ' HEALPIX bad pixel value'
;;-        
;;-        sxaddpar, hdr, 'FILENAME',  outfile,    ' FITS filename'
;;-        sxaddpar, hdr, 'PROCVER',   procver,    ' Product version'
;;-        
;;-        sxaddpar, hdr, 'COMMENT', '  ' , after='PROCVER'
;;-        sxaddpar, hdr, 'COMMENT', ''
;;-        sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;-        sxaddpar, hdr, 'COMMENT', 'For further details see Planck Explanatory Supplement at:'
;;-        sxaddpar, hdr, 'COMMENT', '  http://wiki.cosmos.esa.int/planckpla2015'
;;-        sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;-        
;;-        if keyword_set(test) then print, "------------" else begin
;;-            print, ' >> Exporting '+outdir+outfile
;;-            writefits, outdir+outfile, ima, hdr
;;-            ima=1   ; recover memory
;;-        endelse        
;;- ;       print, hdr
;;-    end 
;;- 
;;-    if defined(stop) then stop
;;-    return
;;- end
;;- 
;;- 
;;- 
;;- 
;;- ;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;;- ; export single frequency map
;;- ;
;;- ;   freq:    frequency .... must fix for detset or bolo
;;- ;   suff:    nominal, halfring_n, survey_n, etc.
;;- ;    tag:    name used in map object that are of form fff{tag}_{type}
;;- ;            where type is of the form I/Q/U/H/II/IQ, .....
;;- ;  ingrp:    group name in db (with trailing '/'
;;- ; procver:   data processing version
;;- ; /temponly:  export Temp only (automatic for 545/857 maps)
;;- ;  /test:    write file with little real data (checks that data exist in db)
;;- ; /check:    only print input map names and output file names
;;- ;-----------------------------------------------------------------------------
;;- ; nov14: modify freqmap to read MAJA's round2 maps as inputs
;;- ;        add FSL correction to covariance(s) read from miss03
;;- ;
;;- ; outstanding: modif of survey maps
;;- ;              modify dsetmap in same manner
;;- ;
;;- ; ATTN: for bolo maps, 
;;- ;-----------------------------------------------------------------------------
;;- 
;;- pro freqmap, freq, relnum, suff, procver, ingrp, tag=tag, outdir=outdir, round2=r2, $
;;-              temp=temp, test=test, stop=stop, verb=verb, check=check
;;- 
;;-   on_error, 0
;;- 
;;-   if n_elements(freq) eq 0 then begin
;;-        print, 'SYNTAX: '
;;-        print, '  freqmap, freq, relnum, suff, tag, ingrp, procver, outdir=outdir,  $'
;;-        print, '           /temp, /test, /stop'    
;;-        return
;;-    endif
;;- 
;;-    if fix(freq) lt 500 then  pol=1 else pol=0      ; to separate cmb from gal channels
;;- 
;;-    if (size(freq))[1] ne 7 then freq = strtrim(freq,2)    ; to use in object names
;;- 
;;-    if n_elements(relnum)  eq 0 then relnum  = 'R0.00'
;;-    if n_elements(suff)    eq 0 then suff    = ''
;;-    if n_elements(procver) eq 0 then procver = 'DX00'
;;-    if not keyword_set(r2)      then r2 = 0
;;-    if n_elements(ingrp)   eq 0 then ingrp   = 'MAP_DX11d_noZodi_2048_GALACTIC_0240_27005/'
;;-    if keyword_set(temp)        then pol=0                ; export temp data only
;;-    if not keyword_set(tag)     then tag     = 'GHz_W_TauDeconv_odc_pol'
;;-    if not keyword_set(outdir)  then outdir  = './'
;;-    if not keyword_set(test)    then test = 0
;;-    if not keyword_set(verb)    then verb = 0
;;-    if not keyword_set(check)   then check = 0 
;;- 
;;-    db  = '/data/dmc/MISS03/DATA/'   &    grp = db+ingrp
;;-    name = 'HFI_SkyMap_'+freq
;;- 
;;-    if test then begin
;;-        Nside = 8LL 
;;-        outdir  = './'
;;-    endif else begin
;;-        Nside = 2048LL
;;-    endelse 
;;-    name = name +'_'+strtrim(Nside,2)+'_'+relnum
;;-    Npix = 12*Nside*Nside
;;-    com = 'begin=0;end='+strtrim(Npix-1,2)
;;-    outfile = name +'_'+suff+'.fits'
;;- 
;;-    ; check if SWB
;;-    if strpos(freq, '-') gt 0 then swb=1 else swb=0
;;-    if swb then pol=0
;;- 
;;-    ; get name of file (round2) or I object (DMC)
;;-    if r2 then begin
;;-        r2suff = r2tag(suff)
;;-        r2indir = '/redtruck/ashdown/repository/exchanges/dx11/maps/hfi_zodi_removed_bpleak_corrected/'
;;-        r2root = 'hfi_dx11d_'   &  r2tail = '_nozodi_bplcorrected_dust_ground.fits'
;;- 
;;-        ; use 'pix' to build name of Round2 files
;;-        if swb then pix = strmid(freq,0,3)+"_"+strmid(freq,4,1) else pix = freq
;;- 
;;-        nameI =  r2indir+r2root+pix+r2suff+r2tail  
;;-        spawn, 'ls -ld '+ nameI +' | cut -d"/" -f9', dmcname
;;- 
;;-    endif else begin
;;-        spawn, ' ls -ld '+db+ingrp+freq+tag+'_I | cut -d"/" -f6-7', dmcname
;;-    endelse
;;- 
;;-    if strlen(dmcname) eq 0 then begin   ;; check that file exists
;;-        print, " File or object not found ... quitting"
;;-        return
;;-    endif 
;;- 
;;-    if check then begin   ;; print input object name and name of output FITS file
;;-        print, form='(A-51, " <== ", A-84)', outfile, dmcname
;;-        if not test then return
;;-    endif 
;;- 
;;-    if strpos(tag, 'nominal') ge 0 then nom=1 else nom=0
;;-    fslcorr = fsl_corr(strmid(freq,0,3))
;;-    if fslcorr gt 1.00d and not nom then fslmsg = " and apply FSL correction of " else fslmsg = " " ;NO FSL corrrection "
;;- 
;;-    if r2 then begin   ;; Round-2 data from MAJA's repository
;;-        if test then begin
;;-            mapI = mrdfits(nameI, 1, h, range=[0,0])  
;;-            mapi = (mapi.(0))[0:Npix-1]
;;-            if pol then begin
;;-                mapq = mapi    &   mapu = mapi   
;;-            endif 
;;-        endif else begin
;;-            mapIn = mrdfits(nameI, 1, h)  
;;-            mapi = reorder(reform(mapIn.(0), 50331648LL), /n2r)
;;-            if pol then begin
;;-                mapq = reorder(reform(mapIn.(1), 50331648LL), /n2r)
;;-                mapu = reorder(reform(mapIn.(2), 50331648LL), /n2r)
;;-            endif 
;;-        endelse 
;;-        print, "Read I data from round-2 FITS file"
;;-        objI = dmcname 
;;-        unitI = sxpar(h, 'TUNIT1')
;;-        if pol then begin
;;-            print, "Read Q,U data from round-2 FITS file"
;;-            objq = obji  &  obju = obji
;;-            unitq = sxpar(h, 'TUNIT2') & unitu = sxpar(h, 'TUNIT3')
;;-        endif 
;;- 
;;-    endif else begin   ;;;;; NOT round-2 data .... Read from DB
;;- 
;;-        if nom then begin
;;-            print, " >> NO BP leakage of FSL correction for nominal mission maps "
;;-            bplq = 0  &  bplu = 0
;;-            fslcorr = 1.
;;-            pol=0
;;-        endif 
;;- 
;;-        print, "Read I data from MISS03"+fslmsg ; +strtrim(fslcorr,2)
;;-        nameI  = grp+freq+tag+'_I'  &    mapI = submap(nameI, 0, fslcorr, com)
;;-        xx = pioreadkeywordobject(unitI,  comm, 'Units', 'PIOSTRING', nameI)
;;-        spawn, "/bin/ls -ldh "+nameI +" | cut -c39-199 | sed s'|"+grp+"|- |'", objI
;;- 
;;-        if pol then begin
;;-            print, " >> Read Q,U etc. data from MISS03"            
;;-            if nom eq 0 then begin ;; get ground dust leakage correction and apply
;;-                print, " >> " + fslmsg; +strtrim(fslcorr,2) 
;;-                print, " >> Applying ground BP leakage correction ####"
;;-                bpldir = db+'DX11d_BP_LEAKAGE_CORRECTION/'
;;-                bpltag = bplsuffix(suff)
;;-                bplq = pioread(bpldir+freq+bpltag+"Q", com=com)   
;;-                bplu = pioread(bpldir+freq+bpltag+"U", com=com)
;;-            end 
;;- 
;;-            nameQ = grp+freq+tag+'_Q'  &  mapQ = submap(nameQ, bplq, fslcorr, com)
;;-            nameU = grp+freq+tag+'_U'  &  mapU = submap(nameU, bplu, fslcorr, com)
;;- 
;;-            xx = pioreadkeywordobject(unitQ , comm, 'Units', 'PIOSTRING', nameQ )
;;-            xx = pioreadkeywordobject(unitU , comm, 'Units', 'PIOSTRING', nameU )
;;-            spawn, "/bin/ls -ldh "+nameQ +" | cut -c39-199 | sed s'|"+grp+"|- |'", objQ 
;;-            spawn, "/bin/ls -ldh "+nameU +" | cut -c39-199 | sed s'|"+grp+"|- |'", objU 
;;-        endif 
;;-    endelse 
;;- 
;;-     ; Hits and covar are always read from MISS03
;;-     print, "Read II covar from MISS03" + fslmsg; +strtrim(fslcorr,2)
;;-     nameII = grp+freq+tag+'_II'  &  mapII = submap(nameII, 0, fslcorr * fslcorr, com)
;;-     xx = pioreadkeywordobject(unitII, comm, 'Units', 'PIOSTRING', nameII)
;;-     print, "Read Hits map from MISS03"
;;-     nameH  = grp+freq+tag+'_H'   &  mapH  = float(pioread(nameH,  com=com))
;;-     xx = pioreadkeywordobject(unitH,  comm, 'Units', 'PIOSTRING', nameH)
;;-     spawn, "/bin/ls -ldh "+nameII+" | cut -c39-199 | sed s'|"+grp+"|- |'", objII
;;-     spawn, "/bin/ls -ldh "+nameH +" | cut -c39-199 | sed s'|"+grp+"|- |'", objH
;;- 
;;-    if pol then begin
;;-        print, "Read other covars from MISS03" + fslmsg+strtrim(fslcorr,2) 
;;-        nameIQ = grp+freq+tag+'_IQ'  &  mapIQ = submap(nameIQ, 0, fslcorr * fslcorr, com)
;;-        nameIU = grp+freq+tag+'_IU'  &  mapIU = submap(nameIU, 0, fslcorr * fslcorr, com)
;;-        nameQQ = grp+freq+tag+'_QQ'  &  mapQQ = submap(nameQQ, 0, fslcorr * fslcorr, com)
;;-        nameQU = grp+freq+tag+'_QU'  &  mapQU = submap(nameQU, 0, fslcorr * fslcorr, com)
;;-        nameUU = grp+freq+tag+'_UU'  &  mapUU = submap(nameUU, 0, fslcorr * fslcorr, com)
;;-            
;;-        xx = pioreadkeywordobject(unitIQ, comm, 'Units', 'PIOSTRING', nameIQ)
;;-        xx = pioreadkeywordobject(unitIU, comm, 'Units', 'PIOSTRING', nameIU)
;;-        xx = pioreadkeywordobject(unitQQ, comm, 'Units', 'PIOSTRING', nameQQ)
;;-        xx = pioreadkeywordobject(unitQU, comm, 'Units', 'PIOSTRING', nameQU)
;;-        xx = pioreadkeywordobject(unitUU, comm, 'Units', 'PIOSTRING', nameUU)
;;- 
;;-        spawn, "/bin/ls -ldh "+nameIQ+" | cut -c39-199 | sed s'|"+grp+"|- |'", objIQ
;;-        spawn, "/bin/ls -ldh "+nameIU+" | cut -c39-199 | sed s'|"+grp+"|- |'", objIU
;;-        spawn, "/bin/ls -ldh "+nameQQ+" | cut -c39-199 | sed s'|"+grp+"|- |'", objQQ
;;-        spawn, "/bin/ls -ldh "+nameQU+" | cut -c39-199 | sed s'|"+grp+"|- |'", objQU
;;-        spawn, "/bin/ls -ldh "+nameUU+" | cut -c39-199 | sed s'|"+grp+"|- |'", objUU
;;-    endif
;;- 
;;- 
;;-    ; Fix Units kwd, i.e. remove [nuI(nu)=cst] when MJy/sr
;;-    Q = strpos( UnitI, '[' , /reverse_s)
;;-    if Q gt 0 then UnitI = strmid(UnitI, 0, q)
;;-    if strmid(UnitI,0,1) eq 'M' then begin
;;-        UnitI  = 'MJy/sr'
;;-        UnitII = '(MJy/sr)^2'
;;-    endif else begin
;;-        UnitII = 'K_CMB^2'
;;-    end 
;;- 
;;-    if not pol then $
;;-      map = replicate({ i_stokes: 0.0 , hits: 0L, ii_cov: 0.0 }, Npix) $
;;-    else $
;;-      map = replicate({i_stokes:0.0, q_stokes: 0.0, u_stokes: 0.0, hits: 0L, $
;;-            ii_cov: 0.0, iq_cov: 0.0, iu_cov: 0.0,  qq_cov: 0.0, qu_cov: 0.0, uu_cov: 0.0}, Npix)
;;-    
;;-    if defined(stop) then stop
;;-    map.i_stokes = reorder(mapI,  /r2n)
;;-    if pol then begin
;;-        map.q_stokes = reorder(mapQ,  /r2n)
;;-        map.u_stokes = reorder(mapU,  /r2n)
;;-        map.hits     = reorder(mapH,  /r2n)
;;-        map.ii_cov   = reorder(mapII, /r2n)
;;-        map.iq_cov   = reorder(mapIQ, /r2n)
;;-        map.iu_cov   = reorder(mapIU, /r2n)
;;-        map.qq_cov   = reorder(mapQQ, /r2n)
;;-        map.qu_cov   = reorder(mapQU, /r2n)
;;-        map.uu_cov   = reorder(mapUU, /r2n)
;;-    endif else begin
;;-        map.hits     = reorder(mapH,  /r2n)
;;-        map.ii_cov   = reorder(mapII, /r2n)
;;-    endelse 
;;- 
;;- 
;;-    ; write header
;;-    fxbhmake, hdr, Npix,       'FREQ-MAP',' Extension name', /init, /date, /extver
;;-    sxdelpar, hdr, 'EXTNAME'
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TFORM3'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;- 
;;-    sxaddpar, hdr, 'TUNIT1',   unitI
;;-    if pol then begin
;;-        sxaddpar, hdr, 'TUNIT2',   unitQ 
;;-        sxaddpar, hdr, 'TUNIT3',   unitU 
;;-        sxaddpar, hdr, 'TUNIT4',   unitH
;;-        sxaddpar, hdr, 'TUNIT5',   unitII
;;-        sxaddpar, hdr, 'TUNIT6',   unitII  ; IQ
;;-        sxaddpar, hdr, 'TUNIT7',   unitII  ; IU
;;-        sxaddpar, hdr, 'TUNIT8',   unitII  ; QQ
;;-        sxaddpar, hdr, 'TUNIT9',   unitII  ; QU
;;-        sxaddpar, hdr, 'TUNIT10',  unitII  ; UU
;;-        sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT10'
;;-    endif else begin
;;-        sxaddpar, hdr, 'TUNIT2',   unitH
;;-        sxaddpar, hdr, 'TUNIT3',   unitII  ;+'^2'   ; TEMPORARY FIX while OP fixes it on his side
;;-        sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT3'
;;-        ;sxaddpar, hdr, 'UNITCONV',  unitconv,     ' (MJy/sr)/K_cmb unit conv. factor'
;;-    endelse 
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, hdr, 'EXTNAME',  'FREQ-MAP',    ' Extension name'
;;-    sxaddpar, hdr, 'PIXTYPE',  'HEALPIX'
;;-    sxaddpar, hdr, 'POLCCONV', 'COSMO',       ' Polarization convention'
;;-    sxaddpar, hdr, 'COORDSYS', 'GALACTIC',    ' Coordinate system'
;;-    sxaddpar, hdr, 'ORDERING', 'NESTED',      ' Healpix ordering'
;;-    sxaddpar, hdr, 'NSIDE',     Nside,        ' Healpix Nside'
;;-    sxaddpar, hdr, 'FIRSTPIX',  0,            ' First pixel # (0 based)'
;;-    sxaddpar, hdr, 'LASTPIX',   Npix-1,       ' Last pixel # (0 based)'
;;-    sxaddpar, hdr, 'FILENAME',  outfile,      ' FITS filename'
;;-    sxaddpar, hdr, 'BAD_DATA',  -1.63750E+30, ' HEALPIX bad pixel value'
;;-    sxaddpar, hdr, 'FREQ',     freq,          ' reference frequency'
;;-    sxaddpar, hdr, 'PROCVER',  procver,       ' Product version'
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='PROCVER'
;;-    sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-    if strpos(suff, 'nominal') gt 0 then mission='nominal' else mission='full'
;;-    sxaddpar, hdr, 'COMMENT',  'Full channel sky map: '+mission+' mission'
;;-    sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT',  'For further details see Planck Explanatory Supplement at:'
;;-    sxaddpar, hdr, 'COMMENT',  '  http://wiki.cosmos.esa.int/planckpla2015'
;;-    sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT',  'HFI-DMC objects: '
;;-    sxaddpar, hdr, 'COMMENT',  'in-group: '+ingrp
;;-    sxaddpar, hdr, 'COMMENT',  'Creation date  - object name'
;;-    sxaddpar, hdr, 'COMMENT',  objI
;;-    if pol then begin
;;-        sxaddpar, hdr, 'COMMENT',  objQ 
;;-        sxaddpar, hdr, 'COMMENT',  objU 
;;-        sxaddpar, hdr, 'COMMENT',  objH 
;;-        sxaddpar, hdr, 'COMMENT',  objII
;;-        sxaddpar, hdr, 'COMMENT',  objIQ
;;-        sxaddpar, hdr, 'COMMENT',  objIU
;;-        sxaddpar, hdr, 'COMMENT',  objQQ
;;-        sxaddpar, hdr, 'COMMENT',  objQU
;;-        sxaddpar, hdr, 'COMMENT',  objUU
;;-    endif else begin
;;-        sxaddpar, hdr, 'COMMENT',  objH 
;;-        sxaddpar, hdr, 'COMMENT',  objII
;;-    endelse 
;;-    sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;- 
;;-    if test then print, ' >>> NOT EXPORTING TEST FILE: '+outdir+outfile else $
;;-      print, ' >> Exporting '+outdir+outfile
;;-    mwrfits, map, outdir+outfile, hdr, /create   ; Full chan freq (sky) maps
;;- ;   mwrfits, map, outdir+'temp.fits', hdr, /create   ; Full chan freq (sky) maps
;;- ;   addchecksum, outdir+'temp.fits', outdir+outfile
;;-    print, ""
;;- 
;;-    if verb then print, hdr
;;-    if defined(stop) then stop
;;- 
;;-    return
;;- end 
;;- 
;;- 
;;- ;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;;- ; export detset map  -  updated for PR2
;;- ;   options are full/nominal and survey 1/2/3/4 only; no halfring; other
;;- ;   params are as for freqmap
;;- ;   dset (string): detset, form= 143-DetSet1  (as used by LAL)
;;- ;                  outnames of form DetsetMap_freq-detset_{1/2}_2048_{nominal/full}
;;- ;-----------------------------------------------------------------------------
;;- 
;;- pro dsetmap, dset, relnum, suff, procver, ingrp, tag=tag, outdir=outdir,  round2=r2,  $
;;-              temponly=temp, check=check, test=test, stop=stop, verb=verb
;;- 
;;-   on_error, 0
;;- 
;;-   if n_elements(dset) eq 0 then begin
;;-        print, 'SYNTAX: '
;;-        print, '  dsetmap, dset, relnum, suff, tag, ingrp, procver, outdir=outdir,  $'
;;-        print, '           /temp, /check, /test, /stop, /verb'    
;;-        return
;;-    endif
;;- 
;;-    freq = strmid(dset,0,3)
;;-    if freq lt fix(500) then  pol=1 else pol=0      ; to separate cmb from gal channels
;;- 
;;-    if n_elements(relnum) eq 0  then relnum  = 'R0.00'
;;-    if n_elements(suff) eq 0    then suff    = ''
;;-    if n_elements(procver) eq 0 then procver = 'DX00'
;;-    if not keyword_set(r2)      then r2 = 0
;;-    if n_elements(ingrp) eq 0   then ingrp   = 'MAP_DX11d_noZodi_2048_GALACTIC_0240_27005/'
;;-    if keyword_set(temp)        then pol=0                ; export temp data only
;;-    if not keyword_set(tag)     then tag     = 'GHz_W_TauDeconv_odc_pol'
;;-    if not keyword_set(outdir)  then outdir  = './'
;;-    if not keyword_set(test)    then test = 0
;;-    if not keyword_set(verb)    then verb = 0
;;-    if not keyword_set(check)   then check = 0 
;;- 
;;-    db  = '/data/dmc/MISS03/DATA/'   &    grp = db+ingrp
;;-    name = 'HFI_SkyMap_'+strmid(dset,0,3)+'-ds'+strmid(dset,10,1)
;;- 
;;-    if test then begin
;;-        Nside = 8LL 
;;-        outdir  = './'
;;-    endif else begin
;;-        Nside = 2048LL
;;-    endelse 
;;-    name = name +'_'+strtrim(Nside,2)+'_'+relnum
;;-    Npix = 12*Nside*Nside
;;-    com = 'begin=0;end='+strtrim(Npix-1,2)
;;-    outfile = name +'_'+suff+'.fits'
;;- 
;;-    if r2 then begin
;;-        r2suff = r2tag(suff)
;;-        r2indir = '/redtruck/ashdown/repository/exchanges/dx11/maps/hfi_zodi_removed_bpleak_corrected/'
;;-        r2root = 'hfi_dx11d_'   &  r2tail = '_nozodi_bplcorrected_dust_ground.fits'
;;-        ; use 'pix' to build name of Round2 files
;;-        pix = strmid(dset,0,3)+"_ds"+strmid(dset,10,1) 
;;-        
;;-        nameI =  r2indir+r2root+pix+r2suff+r2tail  
;;-        spawn, 'ls -ld '+ nameI +' | cut -d"/" -f9', dmcname
;;- 
;;-    endif else begin
;;-        spawn, ' ls -ld '+db+ingrp+dset+tag+'_I | cut -d"/" -f6-7', dmcname
;;-    endelse
;;- 
;;-    if strlen(dmcname) lt 10 then begin
;;-        print, " File or object not found ... quitting"
;;-        return
;;-    endif 
;;- 
;;-    if check then begin
;;-        print, form='(A-51, " <== ", A-84)', outfile, dmcname
;;-        if not test then return
;;-    endif 
;;-    
;;-    fslcorr = fsl_corr( strmid(freq,0,3) )
;;-    if fslcorr gt 1.00d then fslmsg = " and apply FSL correction " else fslmsg = " "
;;- 
;;-    if r2 then begin   ;; Round-2 data from MAJA's repository
;;-        if test then begin
;;-            mapI = mrdfits(nameI, 1, h, range=[0,5])  
;;-            mapi = (reform(mapi.(0), n_elements(mapi.(0))))[0:Npix-1]
;;-            if pol then begin
;;-                mapq = mapi    &   mapu = mapi   
;;-            endif 
;;-        endif else begin
;;-            mapIn = mrdfits(nameI, 1, h)  
;;-            mapi = reorder(reform(mapIn.(0), 50331648LL), /n2r)
;;-            if pol then begin
;;-                mapq = reorder(reform(mapIn.(1), 50331648LL), /n2r)
;;-                mapu = reorder(reform(mapIn.(2), 50331648LL), /n2r)
;;-            endif 
;;-        endelse 
;;-        print, "Read I data from round-2 FITS file"
;;-        objI = dmcname 
;;-        unitI = sxpar(h, 'TUNIT1')
;;-        if pol then begin
;;-            print, "Read Q,U data from round-2 FITS file"
;;-            objq = obji  &  obju = obji
;;-            unitq = sxpar(h, 'TUNIT2') & unitu = sxpar(h, 'TUNIT3')
;;-        endif 
;;- 
;;-    endif else begin
;;-        print, "Read I data from MISS03"+fslmsg+strtrim(fslcorr,2)
;;-        nameI  = grp+dset+tag+'_I'  &    mapI  = submap(nameI, 0, fslcorr, com)
;;-        xx = pioreadkeywordobject(unitI,  comm, 'Units', 'PIOSTRING', nameI)
;;-        spawn, "/bin/ls -ldh "+nameI +" | cut -c39-199 | sed s'|"+grp+"|- |'", objI
;;-        
;;-        if pol then begin
;;-            print, "Read Q,U etc. data from MISS03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" + fslmsg ;+strtrim(fslcorr,2)
;;-            nameQ  = grp+dset+tag+'_Q'   &  mapQ  = submap(nameQ, bplq, fslcorr, com)
;;-            nameU  = grp+dset+tag+'_U'   &  mapU  = submap(nameU, bplu, fslcorr, com)
;;-            xx = pioreadkeywordobject(unitQ , comm, 'Units', 'PIOSTRING', nameQ )
;;-            xx = pioreadkeywordobject(unitU , comm, 'Units', 'PIOSTRING', nameU )
;;-            spawn, "/bin/ls -ldh "+nameQ +" | cut -c39-199 | sed s'|"+grp+"|- |'", objQ 
;;-            spawn, "/bin/ls -ldh "+nameU +" | cut -c39-199 | sed s'|"+grp+"|- |'", objU 
;;-        endif 
;;-    endelse 
;;- 
;;-     ; Hits and covar are always read from MISS03
;;-     print, "Read II covar from MISS03" + fslmsg  ;+strtrim(fslcorr,2)
;;-     nameII = grp+dset+tag+'_II'  &  mapII = submap(nameII, 0, fslcorr * fslcorr, com)
;;-     xx = pioreadkeywordobject(unitII, comm, 'Units', 'PIOSTRING', nameII)
;;-     print, "Read Hits map from MISS03"
;;-     nameH  = grp+dset+tag+'_H'   &  mapH  = float(pioread(nameH,  com=com))
;;-     xx = pioreadkeywordobject(unitH,  comm, 'Units', 'PIOSTRING', nameH)
;;-     spawn, "/bin/ls -ldh "+nameII+" | cut -c39-199 | sed s'|"+grp+"|- |'", objII
;;-     spawn, "/bin/ls -ldh "+nameH +" | cut -c39-199 | sed s'|"+grp+"|- |'", objH
;;- 
;;-    if pol then begin
;;-        print, "Read other covars from MISS03" + fslmsg;  +strtrim(fslcorr,2)
;;-        nameIQ = grp+dset+tag+'_IQ'  &  mapIQ = submap(nameIQ, 0, fslcorr * fslcorr, com)
;;-        nameIU = grp+dset+tag+'_IU'  &  mapIU = submap(nameIU, 0, fslcorr * fslcorr, com)
;;-        nameQQ = grp+dset+tag+'_QQ'  &  mapQQ = submap(nameQQ, 0, fslcorr * fslcorr, com)
;;-        nameQU = grp+dset+tag+'_QU'  &  mapQU = submap(nameQU, 0, fslcorr * fslcorr, com)
;;-        nameUU = grp+dset+tag+'_UU'  &  mapUU = submap(nameUU, 0, fslcorr * fslcorr, com)
;;-        
;;-        xx = pioreadkeywordobject(unitIQ, comm, 'Units', 'PIOSTRING', nameIQ)
;;-        xx = pioreadkeywordobject(unitIU, comm, 'Units', 'PIOSTRING', nameIU)
;;-        xx = pioreadkeywordobject(unitQQ, comm, 'Units', 'PIOSTRING', nameQQ)
;;-        xx = pioreadkeywordobject(unitQU, comm, 'Units', 'PIOSTRING', nameQU)
;;-        xx = pioreadkeywordobject(unitUU, comm, 'Units', 'PIOSTRING', nameUU)
;;- 
;;-        spawn, "/bin/ls -ldh "+nameIQ+" | cut -c39-199 | sed s'|"+grp+"|- |'", objIQ
;;-        spawn, "/bin/ls -ldh "+nameIU+" | cut -c39-199 | sed s'|"+grp+"|- |'", objIU
;;-        spawn, "/bin/ls -ldh "+nameQQ+" | cut -c39-199 | sed s'|"+grp+"|- |'", objQQ
;;-        spawn, "/bin/ls -ldh "+nameQU+" | cut -c39-199 | sed s'|"+grp+"|- |'", objQU
;;-        spawn, "/bin/ls -ldh "+nameUU+" | cut -c39-199 | sed s'|"+grp+"|- |'", objUU
;;-    endif 
;;- 
;;-    ; Fix Units kwd, i.e. remove [nuI(nu)=cst] when MJy/sr
;;-    Q = strpos( UnitI, '[' , /reverse_s )
;;-    if Q gt 0 then UnitI = strmid(UnitI, 0, q)
;;-    if strmid(UnitI,0,1) eq 'M' then begin
;;-        UnitI  = 'MJy/sr'
;;-        UnitII = '(MJy/sr)^2'
;;-    endif else begin
;;-        UnitII = 'K_CMB^2'
;;-    end 
;;- 
;;-    if not pol then $
;;-      map = replicate({ i_stokes: 0.0 , hits: 0L, ii_cov: 0.0 }, Npix) $
;;-    else $
;;-      map = replicate({i_stokes:0.0, q_stokes: 0.0, u_stokes: 0.0, hits: 0L, $
;;-            ii_cov: 0.0, iq_cov: 0.0, iu_cov: 0.0,  qq_cov: 0.0, qu_cov: 0.0, uu_cov: 0.0}, Npix)
;;-    
;;-    if defined(stop) then stop
;;-    map.i_stokes = reorder(mapI,  /r2n)
;;-    if pol then begin
;;-        map.q_stokes = reorder(mapQ,  /r2n)
;;-        map.u_stokes = reorder(mapU,  /r2n)
;;-        map.hits     = reorder(mapH,  /r2n)
;;-        map.ii_cov   = reorder(mapII, /r2n)
;;-        map.iq_cov   = reorder(mapIQ, /r2n)
;;-        map.iu_cov   = reorder(mapIU, /r2n)
;;-        map.qq_cov   = reorder(mapQQ, /r2n)
;;-        map.qu_cov   = reorder(mapQU, /r2n)
;;-        map.uu_cov   = reorder(mapUU, /r2n)
;;-    endif else begin
;;-        map.hits     = reorder(mapH,  /r2n)
;;-        map.ii_cov   = reorder(mapII, /r2n)
;;-    endelse 
;;- 
;;- 
;;-    ; write header
;;-    fxbhmake, hdr, Npix,       'FREQ-MAP',' Extension name', /init, /date, /extver
;;-    sxdelpar, hdr, 'EXTNAME'
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TFORM3'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;- 
;;-    sxaddpar, hdr, 'TUNIT1',   unitI
;;-    if pol then begin
;;-        sxaddpar, hdr, 'TUNIT2',   unitQ 
;;-        sxaddpar, hdr, 'TUNIT3',   unitU 
;;-        sxaddpar, hdr, 'TUNIT4',   unitH
;;-        sxaddpar, hdr, 'TUNIT5',   unitII
;;-        sxaddpar, hdr, 'TUNIT6',   unitII  ; IQ
;;-        sxaddpar, hdr, 'TUNIT7',   unitII  ; IU
;;-        sxaddpar, hdr, 'TUNIT8',   unitII  ; QQ
;;-        sxaddpar, hdr, 'TUNIT9',   unitII  ; QU
;;-        sxaddpar, hdr, 'TUNIT10',  unitII  ; UU
;;-        sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT10'
;;-    endif else begin
;;-        sxaddpar, hdr, 'TUNIT2',   unitH
;;-        sxaddpar, hdr, 'TUNIT3',   unitII  ;+'^2'   ; TEMPORARY FIX while OP fixes it on his side
;;-        sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT3'
;;-        ;sxaddpar, hdr, 'UNITCONV',  unitconv,     ' (MJy/sr)/K_cmb unit conv. factor'
;;-    endelse 
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, hdr, 'EXTNAME',  'FREQ-MAP',    ' Extension name'
;;-    sxaddpar, hdr, 'PIXTYPE',  'HEALPIX'
;;-    sxaddpar, hdr, 'COORDSYS', 'GALACTIC',    ' Coordinate system'
;;-    sxaddpar, hdr, 'ORDERING', 'NESTED',      ' Healpix ordering'
;;-    sxaddpar, hdr, 'NSIDE',     Nside,        ' Healpix Nside'
;;-    sxaddpar, hdr, 'FIRSTPIX',  0,            ' First pixel # (0 based)'
;;-    sxaddpar, hdr, 'LASTPIX',   Npix-1,       ' Last pixel # (0 based)'
;;-    sxaddpar, hdr, 'FILENAME',  outfile,      ' FITS filename'
;;-    sxaddpar, hdr, 'BAD_DATA',  -1.63750E+30, ' HEALPIX bad pixel value'
;;-    sxaddpar, hdr, 'FREQ',     freq,          ' reference frequency'
;;-    sxaddpar, hdr, 'PROCVER',  procver,       ' Product version'
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='PROCVER'
;;-    sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-    if strpos(suff, 'nominal') gt 0 then mission='nominal' else mission='full'
;;-    sxaddpar, hdr, 'COMMENT',  'Full channel sky map: '+mission+' mission, basic product'
;;-    sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT',  'For further details see Planck Explanatory Supplement at:'
;;-    sxaddpar, hdr, 'COMMENT',  '  http://wiki.cosmos.esa.int/planckpla2015'
;;-    sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT',  'HFI-DMC objects: '
;;-    sxaddpar, hdr, 'COMMENT',  'in-group: '+ingrp
;;-    sxaddpar, hdr, 'COMMENT',  'Creation date  - object name'
;;-    sxaddpar, hdr, 'COMMENT',  objI
;;-    if pol then begin
;;-        sxaddpar, hdr, 'COMMENT',  objQ 
;;-        sxaddpar, hdr, 'COMMENT',  objU 
;;-        sxaddpar, hdr, 'COMMENT',  objH 
;;-        sxaddpar, hdr, 'COMMENT',  objII
;;-        sxaddpar, hdr, 'COMMENT',  objIQ
;;-        sxaddpar, hdr, 'COMMENT',  objIU
;;-        sxaddpar, hdr, 'COMMENT',  objQQ
;;-        sxaddpar, hdr, 'COMMENT',  objQU
;;-        sxaddpar, hdr, 'COMMENT',  objUU
;;-    endif else begin
;;-        sxaddpar, hdr, 'COMMENT',  objH 
;;-        sxaddpar, hdr, 'COMMENT',  objII
;;-    endelse 
;;-    sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-    
;;-    if test then print, ' >>> EXPORTING TEST FILE: '+outdir+outfile else $
;;-      print, ' >> Exporting '+outdir+outfile
;;- 
;;-    mwrfits, map, outdir+outfile, hdr, /create   ; Full chan freq (sky) maps
;;- ;   mwrfits, map, outdir+'temp.fits', hdr, /create   ; Full chan freq (sky) maps
;;- ;   addchecksum, outdir+'temp.fits', outdir+outfile
;;-    print, ""
;;- 
;;-    if verb then print, hdr
;;-    if defined(stop) then stop
;;- 
;;-    return
;;- end 
;;- 
;;- 
;;- ;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;;- ; export zodi-correction maps
;;- ;   built from diff of regular - zodi-corrected maps from LAL
;;- ;-----------------------------------------------------------------------------
;;- 
;;- 
;;- pro zodicorr, freq, relnum, suff, procver, temp=temp, outdir=outdir, $
;;-               stop=stop, check=check, test=test
;;- 
;;-    on_error, 0
;;- 
;;-    if fix(freq) lt 500 then  pol=1 else pol=0      ; to separate cmb from gal channels
;;- ;   if fix(freq) eq 353 then  pol=1 else pol=0       ; to separate cmb from gal channels
;;- 
;;- ;   if (size(freq))[1] ne 7 then freq = strtrim(freq,2)    ; to use in object names
;;- 
;;-    if n_elements(relnum)  eq 0 then relnum  = 'R0.00'
;;-    if n_elements(procver) eq 0 then procver = 'DX00'
;;-    if n_elements(suff)    eq 0 then suff    = ''
;;-    if keyword_set(temp)        then pol=0                ; export temp data only
;;-    if not keyword_set(outdir)  then outdir  = './'
;;-    if not keyword_set(test)    then test = 0
;;-    if not keyword_set(verb)    then verb = 0
;;-    if not keyword_set(check)   then check = 0 
;;- 
;;-    stokes = ['_I','_Q','_U']
;;-    Nside = 2048L   &   Npix = 12*Nside*Nside
;;- 
;;-    db = "/data/dmc/MISS03/DATA/"
;;-    grpzodi  = 'MAP_DX11d_noZodi_2048_GALACTIC_0240_27005/' 
;;- 
;;-    outfile = 'HFI_CorrMap_'+freq+'-Zodi_'+strtrim(Nside,2)+'_'+relnum+'_'+suff+'.fits'
;;-    if fix(freq) le 300 then tail="GHz_W_TauDeconv_odc_pol"
;;-    if fix(freq) eq 353 then tail="GHz_W_TauDeconv_odc"
;;-    if fix(freq) ge 500 then tail="GHz_W_TauDeconv_planet"
;;-    
;;-    if strpos(suff, "survey") ge 0 then tail = tail+"_survey"+strmid(suff,7,1)
;;-    if strpos(suff, "year")   ge 0 then tail = tail+"_year"  +strmid(suff,5,1)
;;-    if strpos(suff, "halfmission") ge 0 then tail = tail+"_hm"+strmid(suff,12,1)
;;- 
;;-    if fix(freq) lt 400  then begin
;;-        grpmain = "MAP_DX11d_2048_GALACTIC_0240_27005/" 
;;-        zodi = replicate({ I_stokes:0.0, Q_stokes:0.0, U_stokes:0.0}, Npix)
;;-        Unit = 'K_CMB'
;;-    endif else begin 
;;-        grpmain = "MAP_DX11c_2048_GALACTIC_0240_27005/"
;;-        zodi = replicate({ I_stokes:0.0 }, Npix)
;;-        Unit = 'MJy/sr' 
;;-    endelse 
;;-    
;;-    print, " >> check input maps: "
;;-    spawn, "ls -d "+ db+grpmain + freq + tail + "_I"
;;-    spawn, "ls -d "+ db+grpzodi + freq + tail + "_I"
;;-    if pol then begin
;;-        spawn, "ls -d "+ db+grpmain + freq + tail + "_Q"
;;-        spawn, "ls -d "+ db+grpzodi + freq + tail + "_Q"
;;-        spawn, "ls -d "+ db+grpmain + freq + tail + "_U"
;;-        spawn, "ls -d "+ db+grpzodi + freq + tail + "_U"
;;-    endif 
;;- 
;;-    if check then begin
;;-  ;      print, "     - input:  "+grpmain + freq + tail + "_I"
;;-  ;      print, "     - output: "+outfile
;;-        print, form='(A-51, " <== ", A-84)', outfile, grpmain + freq + tail + "_I"
;;-        return
;;-    endif 
;;- 
;;-    base = pioread(db+grpmain + freq + tail + "_I")
;;-    zcor = pioread(db+grpzodi + freq + tail + "_I")
;;-    zodi.(0) = reorder( float(base - zcor), /r2n)   &  if max(zodi.(0)) eq 0 then stop
;;- 
;;-    if pol then begin
;;-        print, "  >> read "+freq+" "+suff+" Q,U maps from db"
;;-        for l=1,2 do begin
;;-            base = pioread(db+grpmain + freq + tail + stokes[l])
;;-            zcor = pioread(db+grpzodi + freq + tail + stokes[l])
;;-            zodi.(l) = reorder( float(base - zcor), /r2n)   &  if max(zodi.(l)) eq 0 then stop
;;-        endfor
;;-    endif 
;;- 
;;-    ; write header
;;-    fxbhmake, hdr, Npix,       'FREQ-MAP',' Extension name', /init, /date, /extver
;;-    sxdelpar, hdr, 'EXTNAME'
;;-    if pol then sxaddpar, hdr, 'COMMENT',  '  ' , after='TFORM3' $  
;;-      else sxaddpar, hdr, 'COMMENT',  '  ' , after='TFORM1'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, hdr, 'TUNIT1',   unit
;;-    if pol then begin
;;-        sxaddpar, hdr, 'TUNIT2',   unit
;;-        sxaddpar, hdr, 'TUNIT3',   unit
;;-    endif 
;;-    
;;-    if pol then sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT3' $  
;;-      else sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT1'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, hdr, 'EXTNAME',  'ZODI-COR',    ' Extension name'
;;-    sxaddpar, hdr, 'PIXTYPE',  'HEALPIX'
;;-    sxaddpar, hdr, 'POLCCONV', 'COSMO',       ' Polarization convention'
;;-    sxaddpar, hdr, 'COORDSYS', 'GALACTIC',    ' Coordinate system'
;;-    sxaddpar, hdr, 'ORDERING', 'NESTED',      ' Healpix ordering'
;;-    sxaddpar, hdr, 'NSIDE',     Nside,        ' Healpix Nside'
;;-    sxaddpar, hdr, 'FIRSTPIX',  0,            ' First pixel # (0 based)'
;;-    sxaddpar, hdr, 'LASTPIX',   Npix-1,       ' Last pixel # (0 based)'
;;-    sxaddpar, hdr, 'FILENAME',  outfile,      ' FITS filename'
;;-    sxaddpar, hdr, 'BAD_DATA',  -1.63750E+30, ' HEALPIX bad pixel value'
;;-    sxaddpar, hdr, 'FREQ',     freq,          ' reference frequency'
;;-    sxaddpar, hdr, 'PROCVER',  procver,       ' Product version'
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='PROCVER'
;;-    sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT',  'Full channel, '+suff+' Zodi-correction map '
;;-    sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT',  'For further details see Planck Explanatory Supplement at:'
;;-    sxaddpar, hdr, 'COMMENT',  '  http://wiki.cosmos.esa.int/planckpla2015'
;;-    sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;- 
;;-    print, "  >> Exporting "+outdir+outfile
;;-    mwrfits, zodi, outdir+outfile, hdr, /create  
;;- 
;;-    if keyword_set(stop) then stop
;;-    return
;;- end 
;;- 
;;- ;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;;- ; export leakage correction Q,U maps  -  updated for PR2
;;- ; >> MUST set tag appropriately to specify CO or dust etc.
;;- ; - file list is determined automatically
;;- ; ==> write into one file per freq. channel
;;- ; - tag:   used to build name of input object (w/o initial/final "_"
;;- ; - comp, suff:  used to build name of output file
;;- ;-----------------------------------------------------------------------------
;;- 
;;- pro leakage, relnum, procver, tag=tag, ingrp=ingrp, comp=comp, suff=suff, outdir=outdir, $
;;-              test=test, verb=verb, check=check, stop=stop
;;- 
;;-    on_error, 2
;;- 
;;-    if n_elements(procver) eq 0 then begin
;;-        print, "SYNTAX:"
;;-        print, "  leakage, procver, tag=tag, relnum=relnum, outdir=outdir, "
;;-        print, "  /verb, /test, /check "
;;-        print, "  tag: string used to build input object names, eg. 'GHz_dustleak_ground' "
;;-        return
;;-    endif 
;;- 
;;-    if keyword_set(test)       then test=1 else test=0
;;-    if not keyword_set(relnum) then relnum  = 'R0.00'
;;-    if not keyword_set(outdir) then outdir  = './'
;;-    if not keyword_set(check ) then check   = 0
;;-    if not keyword_set(suff)   then suff = ''
;;-    if keyword_set(verb)       then verb = 1 else verb = 0
;;- 
;;-    ingrp = relnum+'_BP_LEAKAGE_CORRECTION/'
;;-    db = '/data/dmc/MISS03/DATA/'
;;-    freqs = ['100','143','217','353']
;;- 
;;-    if strpos(tag, "DetSet") ge 0 then ds = 1 else ds = 0
;;-    spawn, 'dmctool getobjectlist '+db+ingrp+' | grep ' +tag+'_Q | cut -d"/" -f7 | cut -c1-3', freqs 
;;- 
;;-    freqs = ['353']    ; for PR2:
;;-    freqs = ['100','143','217']    ; other freqs (10.mar.16)
;;- 
;;-    if strlen(freqs[0]) gt 0 then begin
;;- ;       print, " >> FOUND following freqs: ", freqs, "    for "+ tag
;;- 
;;-       for f=0,n_elements(freqs)-1 do begin
;;-           freq = freqs[f]
;;-           if ds then name = 'HFI_CorrMap_'+freq+'-ds'+strmid(tag,6,1) else name = 'HFI_CorrMap_'+freq
;;-       
;;-           if test then begin
;;-               Nside = 16LL 
;;-               outdir  = './'
;;-           endif else begin
;;-               Nside = 2048LL
;;-           endelse 
;;-       
;;-           q2 = strpos(tag, "leak")   &   sub = strmid(tag, q2+5, 9)
;;-           q4 = strpos(sub, "_")      &   if q4 gt 0 then sub=strmid(sub,0,q4)
;;-    
;;-           name = name +'-'+comp+'leak-'+sub+'_'+strtrim(Nside,2)+'_R2.00'+suff
;;-           Npix = 12*Nside*Nside
;;-           com = 'begin=0;end='+strtrim(Npix-1,2)
;;-           outfile = name + '.fits'
;;-       
;;-           grp = db+ingrp
;;-           if ds then $
;;-             spawn, 'dmctool getobjectlist '+db+ingrp +' | grep '+freq+"-"+tag+'_Q ', dmcname else $
;;-             spawn, 'dmctool getobjectlist '+db+ingrp +' | grep '+freq+tag+'_Q ', dmcname
;;-       
;;-           if check then begin
;;-               print, form='("  ", A-35, " ==>   ", A-72)', dmcname, outfile
;;-           endif else begin
;;-               map = replicate({q_stokes: 0.0, u_stokes: 0.0}, Npix)
;;-               nameU = (nameQ = grp+dmcname)  &   strput, nameU, 'U', strlen(nameU)-1
;;-       
;;-               map.u_stokes = reorder(float(pioread(nameU, com=com)), /r2n)
;;-               map.q_stokes = reorder(float(pioread(nameQ, com=com)), /r2n)
;;-           
;;-               ; write header
;;-               fxbhmake, hdr, Npix,       'FREQ-MAP',' Extension name', /init, /date, /extver
;;-               sxdelpar, hdr, 'EXTNAME'
;;-       
;;-               sxaddpar, hdr, 'COMMENT',  '  ' , after='TFORM2'
;;-               sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
;;-               sxaddpar, hdr, 'COMMENT',  '  '
;;-               sxaddpar, hdr, 'TUNIT1',   'K_CMB'
;;-               sxaddpar, hdr, 'TUNIT2',   'K_CMB' 
;;-               sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT2'
;;-               sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-               sxaddpar, hdr, 'COMMENT',  '  '
;;-               
;;-               sxaddpar, hdr, 'EXTNAME',  'CORR-MAP',    ' Extension name'
;;-               sxaddpar, hdr, 'PIXTYPE',  'HEALPIX'
;;-               sxaddpar, hdr, 'POLCCONV', 'COSMO',       ' Polarization convention'
;;-               sxaddpar, hdr, 'COORDSYS', 'GALACTIC',    ' Coordinate system'
;;-               sxaddpar, hdr, 'ORDERING', 'NESTED',      ' Healpix ordering'
;;-               sxaddpar, hdr, 'NSIDE',     Nside,        ' Healpix Nside'
;;-               sxaddpar, hdr, 'FIRSTPIX',  0,            ' First pixel # (0 based)'
;;-               sxaddpar, hdr, 'LASTPIX',   Npix-1,       ' Last pixel # (0 based)'
;;-               sxaddpar, hdr, 'FILENAME',  outfile,      ' FITS filename'
;;-               sxaddpar, hdr, 'BAD_DATA',  -1.63750E+30, ' HEALPIX bad pixel value'
;;-               sxaddpar, hdr, 'FREQ',     freq,          ' reference frequency'
;;-               sxaddpar, hdr, 'PROCVER',  procver,       ' Product version'
;;-               sxaddpar, hdr, 'COMMENT',  '  ' , after='PROCVER'
;;-               sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-               sxaddpar, hdr, 'COMMENT',  'Leakage correction map'
;;-               sxaddpar, hdr, 'COMMENT',  'For further details see Planck Explanatory Supplement at:'
;;-               sxaddpar, hdr, 'COMMENT',  '  http://wiki.cosmos.esa.int/planckpla2015'
;;-               sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-               
;;-               if test then print, ' >>> EXPORTING TEST FILE: '+outdir+outfile else $
;;-                 print, ' >> Exporting '+outdir+outfile
;;-       
;;-               mwrfits, map, outdir+outfile, hdr, /create ; Full chan freq (sky) maps
;;-               ;addchecksum, outdir+'temp.fits', outdir+outfile
;;-               if verb then print, hdr
;;-               if defined(stop) then stop
;;-           endelse 
;;-       endfor 
;;-    endif  else begin
;;-        print,  ' >> NO OBJECTS with "' +tag+'" found'
;;-    endelse    
;;- 
;;-    spawn, 'rm -rf temp.fits'
;;- 
;;-    return
;;- end 
;;- 
;;- ;-----------------------------------------------------------------------------
;;- ; Generalized Global Filt leakage correction
;;- ;-----------------------------------------------------------------------------
;;- 
;;- pro ggf, relnum, procver, tag=tag, comp=comp, suff=suff, outdir=outdir, $
;;-          test=test, verb=verb, check=check, stop=stop
;;- 
;;-    on_error, 2
;;- 
;;-    if n_elements(procver) eq 0 then begin
;;-        print, "SYNTAX:"
;;-        print, "  ggf, procver, tag=tag, relnum=relnum, outdir=outdir, "
;;-        print, "  /verb, /test, /check "
;;-        print, "  tag: string used to build input object names, eg. 'GHz_dustleak_ground' "
;;-        return
;;-    endif 
;;- 
;;-    if keyword_set(test)       then test=1 else test=0
;;-    if not keyword_set(relnum) then relnum  = 'R0.00'
;;-    if not keyword_set(outdir) then outdir  = './'
;;-    if not keyword_set(check ) then check   = 0
;;-    if not keyword_set(suff)   then suff = ''
;;-    if keyword_set(verb)       then verb = 1 else verb = 0
;;- 
;;-    if test then begin
;;-        Nside = 8LL 
;;-        outdir  = './'
;;-    endif else begin
;;-        Nside = 2048LL
;;-    endelse 
;;-    Npix = 12*Nside*Nside
;;-    
;;-    grp = '/data/dmc/MISS03/DATA/MAPS_DX11d_GGF_201409/'
;;-    ;freqs = ['100','143','217','353']
;;- 
;;-    dsi = strpos(tag, "DetSet") 
;;-    if dsi ge 0 then ds = 1 else ds = 0
;;- 
;;-    ;spawn, 'dmctool getobjectlist '+grp+' | grep ' +tag+'_Q | cut -d"/" -f7 | cut -c1-3', freqs 
;;-    freqs = ['353']    ; for PR2:
;;- 
;;-    if strlen(freqs[0]) gt 0 then begin
;;-       for f=0,n_elements(freqs)-1 do begin
;;-           freq = freqs[f]
;;-           if ds then name = 'HFI_CorrMap_'+freq+'-ds'+strmid(tag,dsi+6,1) else name = 'HFI_CorrMap_'+freq
;;-    
;;-           name = name +'-leakage-global_'+strtrim(Nside,2)+'_R2.00'+suff
;;-           com = 'begin=0;end='+strtrim(Npix-1,2)
;;-           outfile = name + '.fits'
;;-           spawn, 'dmctool getobjectlist '+grp +' | grep '+freq+tag+'_Q ', dmcname
;;-       
;;-           if check then begin
;;-               print, form='("  ", A-35, " ==>   ", A-72)', dmcname, outfile
;;-               return
;;- 
;;-           endif else begin
;;-               map = replicate({q_stokes: 0.0, u_stokes: 0.0}, Npix)
;;-               nameU = (nameQ = grp+dmcname)  &   strput, nameU, 'U', strlen(nameU)-1
;;-       
;;-               map.q_stokes = reorder(float(pioread(nameQ, com=com)), /r2n)
;;-               map.u_stokes = reorder(float(pioread(nameU, com=com)), /r2n)
;;-           
;;-               ; write header
;;-               fxbhmake, hdr, Npix,       'FREQ-MAP',' Extension name', /init, /date, /extver
;;-               sxdelpar, hdr, 'EXTNAME'
;;-       
;;-               sxaddpar, hdr, 'COMMENT',  '  ' , after='TFORM2'
;;-               sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
;;-               sxaddpar, hdr, 'COMMENT',  '  '
;;-               sxaddpar, hdr, 'TUNIT1',   'K_CMB'
;;-               sxaddpar, hdr, 'TUNIT2',   'K_CMB' 
;;- 
;;-               sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT2'
;;-               sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-               sxaddpar, hdr, 'COMMENT',  '  '
;;-               
;;-               sxaddpar, hdr, 'EXTNAME',  'CORR-MAP',    ' Extension name'
;;-               sxaddpar, hdr, 'PIXTYPE',  'HEALPIX'
;;-               sxaddpar, hdr, 'POLCCONV', 'COSMO',       ' Polarization convention'
;;-               sxaddpar, hdr, 'COORDSYS', 'GALACTIC',    ' Coordinate system'
;;-               sxaddpar, hdr, 'ORDERING', 'NESTED',      ' Healpix ordering'
;;-               sxaddpar, hdr, 'NSIDE',     Nside,        ' Healpix Nside'
;;-               sxaddpar, hdr, 'FIRSTPIX',  0,            ' First pixel # (0 based)'
;;-               sxaddpar, hdr, 'LASTPIX',   Npix-1,       ' Last pixel # (0 based)'
;;-               sxaddpar, hdr, 'FILENAME',  outfile,      ' FITS filename'
;;-               sxaddpar, hdr, 'BAD_DATA',  -1.63750E+30, ' HEALPIX bad pixel value'
;;-               sxaddpar, hdr, 'FREQ',     freq,          ' reference frequency'
;;-               sxaddpar, hdr, 'PROCVER',  procver,       ' Product version'
;;-               sxaddpar, hdr, 'COMMENT',  '  ' , after='PROCVER'
;;-               sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-               sxaddpar, hdr, 'COMMENT',  'Global fit leakage correction map'
;;-               sxaddpar, hdr, 'COMMENT',  'For further details see Planck Explanatory Supplement at:'
;;-               sxaddpar, hdr, 'COMMENT',  '  http://wiki.cosmos.esa.int/planckpla2015'
;;-               sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-               
;;-               if test then print, ' >>> EXPORTING TEST FILE: '+outdir+outfile else $
;;-                 print, ' >> Exporting '+outdir+outfile
;;-       
;;-               mwrfits, map, outdir+outfile, hdr, /create ; Full chan freq (sky) maps
;;-               ;addchecksum, outdir+'temp.fits', outdir+outfile
;;-               if verb then print, hdr
;;-               if defined(stop) then stop
;;-           endelse 
;;-       endfor 
;;-    endif  else begin
;;-        print,  ' >> NO OBJECTS with "' +tag+'" found'
;;-    endelse    
;;- 
;;-    spawn, 'rm -rf temp.fits'
;;- 
;;-    return
;;- end 
;;- 
;;- ;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;;- ; export galactic plane masks
;;- ;-----------------------------------------------------------------------------
;;- 
;;- pro galmasks, relnum, procver, apo=apo, outdir=outdir, indir=indir,  $
;;-            test=test, stop=stop, verb=verb
;;- 
;;- ;  on_error, 2
;;- 
;;-   if n_elements(relnum) eq 0 then begin
;;-        print, 'SYNTAX: '
;;-        print, '  skymasks, relnum, procver, outdir=outdir, indir=indir,  $'
;;-        print, '           /test, /stop, /verb'    
;;-        return
;;-    endif
;;- 
;;-    if n_elements(relnum)  eq 0 then relnum  = 'R0.00'
;;-    if n_elements(procver) eq 0 then procver = 'test'
;;-    if not keyword_set(outdir)  then outdir  = './'
;;-    if keyword_set(test)        then test = 1 else test = 0
;;-    if keyword_set(verb)        then verb = 1 else verb = 0
;;-    if keyword_set(stop)        then stop = 1 else stop = 0
;;- 
;;-    Nside = 2048LL   &   Npix = 12*Nside*Nside
;;- 
;;-    ;indir = '/redtruck/ashdown/repository/exchanges/dx9/masks/galactic/'
;;-    indir = '/data/abenoitl/aurelien/masks_consistency/DX10/Masks/Galactic_I/mask_gal_fsky0p'
;;- 
;;-    case apo of
;;-        0: begin                 ; not apodized
;;-            tail  = '.fits.gz'          & outfile = 'HFI_Mask_GalPlane-apo0_'+strtrim(Nside,2)+'_'+relnum+'.fits'
;;-        end 
;;-        2: begin                 ; 2 deg apodization   
;;-            tail  = '_apo2_V3.fits.gz'  & outfile = 'HFI_Mask_GalPlane-apo2_'+strtrim(Nside,2)+'_'+relnum+'.fits'
;;-        end 
;;-        5: begin ; 5 deg apodization   
;;-            tail  = '_apo5_V3.fits.gz'  & outfile = 'HFI_Mask_GalPlane-apo5_'+strtrim(Nside,2)+'_'+relnum+'.fits'
;;-        end 
;;-        else: begin
;;-            print, "not available - quitting"
;;-        end 
;;-    endcase
;;-  
;;-    ; build structure, then fill it in
;;-    if apo gt 0 then $
;;-      gal = replicate({gal020: 0., gal040: 0., gal060: 0., gal070: 0., $
;;-                     gal080: 0., gal090: 0., gal097: 0., gal099: 0. }, Npix ) $
;;-      else $
;;-      gal = replicate({gal020: 0B, gal040: 0B, gal060: 0B, gal070: 0B, $
;;-                     gal080: 0B, gal090: 0B, gal097: 0B, gal099: 0B }, Npix ) 
;;- 
;;-    gal.(0)  = reorder(reform((mrdfits(indir+'20'+tail, 1, /sil, comp='gunzip')).(0), npix), /r2n)
;;-    gal.(1)  = reorder(reform((mrdfits(indir+'40'+tail, 1, /sil, comp='gunzip')).(0), npix), /r2n)
;;-    gal.(2)  = reorder(reform((mrdfits(indir+'60'+tail, 1, /sil, comp='gunzip')).(0), npix), /r2n)
;;-    gal.(3)  = reorder(reform((mrdfits(indir+'70'+tail, 1, /sil, comp='gunzip')).(0), npix), /r2n)
;;-    gal.(4)  = reorder(reform((mrdfits(indir+'80'+tail, 1, /sil, comp='gunzip')).(0), npix), /r2n)
;;-    gal.(5)  = reorder(reform((mrdfits(indir+'90'+tail, 1, /sil, comp='gunzip')).(0), npix), /r2n)
;;-    gal.(6)  = reorder(reform((mrdfits(indir+'97'+tail, 1, /sil, comp='gunzip')).(0), npix), /r2n)
;;-    gal.(7)  = reorder(reform((mrdfits(indir+'99'+tail, 1, /sil, comp='gunzip')).(0), npix), /r2n)
;;- 
;;-    ; write header for GAL masks
;;-    fxbhmake, hdr, Npix,       'GAL-MASK',' Extension name', /init, /date, /extver
;;-    sxdelpar, hdr, 'EXTNAME'
;;-    sxaddpar, hdr, 'TTYPE1',  'GAL020,', ' 20% sky coverage'
;;-    sxaddpar, hdr, 'TTYPE2',  'GAL040,', ' 40% sky coverage'
;;-    sxaddpar, hdr, 'TTYPE3',  'GAL060,', ' 60% sky coverage'
;;-    sxaddpar, hdr, 'TTYPE4',  'GAL070,', ' 70% sky coverage'
;;-    sxaddpar, hdr, 'TTYPE5',  'GAL080,', ' 80% sky coverage'
;;-    sxaddpar, hdr, 'TTYPE6',  'GAL090,', ' 90% sky coverage'
;;-    sxaddpar, hdr, 'TTYPE7',  'GAL097,', ' 97% sky coverage'
;;-    sxaddpar, hdr, 'TTYPE8',  'GAL099,', ' 99% sky coverage'
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TTYPE8'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;- 
;;-    sxaddpar, hdr, 'EXTNAME',  'GAL-MASK',    ' Extension name'
;;-    sxaddpar, hdr, 'PIXTYPE',  'HEALPIX'
;;-    sxaddpar, hdr, 'POLCCONV', 'COSMO',       ' Polarization convention'
;;-    sxaddpar, hdr, 'COORDSYS', 'GALACTIC',    ' Coordinate system'
;;-    sxaddpar, hdr, 'ORDERING', 'NESTED',      ' Healpix ordering'
;;-    sxaddpar, hdr, 'NSIDE',     Nside,        ' Healpix Nside'
;;-    sxaddpar, hdr, 'FIRSTPIX',  0,            ' First pixel # (0 based)'
;;-    sxaddpar, hdr, 'LASTPIX',   Npix-1,       ' Last pixel # (0 based)'
;;-    sxaddpar, hdr, 'FILENAME',  outfile,      ' FITS filename'
;;-    sxaddpar, hdr, 'PROCVER',  procver,       ' Product version'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------', after='PROCVER'
;;-    sxaddpar, hdr, 'COMMENT', 'Galactic emission masks, based on 353 GHz emission, apodized, and for '
;;-    sxaddpar, hdr, 'COMMENT', 'various fractions of skycoverage. For general purpose usage.'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT',  'For further details see Planck Explanatory Supplement at:'
;;-    sxaddpar, hdr, 'COMMENT',  '  http://wiki.cosmos.esa.int/planckpla2015'
;;-    sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;- 
;;-    print, ' >> Exporting '+outdir+outfile
;;-    mwrfits, gal, outdir+outfile, hdr, /create   ; Galactic Plane masks
;;- 
;;-    if verb then print, hdr
;;-    if stop then stop
;;- 
;;-    return
;;- end
;;- 
;;- ;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;;- ; export point source masks: one for Int and two for Pol in two files
;;- ; - see readme in indir for details
;;- ; updated for PR2
;;- ;-----------------------------------------------------------------------------
;;- 
;;- pro srcmasks, relnum, procver, outdir=outdir, indir=indir,  $
;;-            stop=stop, verb=verb
;;- 
;;- ;  on_error, 2
;;- 
;;-   if n_elements(relnum) eq 0 then begin
;;-        print, 'SYNTAX: '
;;-        print, '  srcmasks, relnum, procver, outdir=outdir, indir=indir,  $'
;;-        print, '           /test, /stop, /verb'    
;;-        return
;;-    endif
;;- 
;;-    if n_elements(relnum)  eq 0 then relnum  = 'R0.00'
;;-    if n_elements(procver) eq 0 then procver = 'test'
;;-    if not keyword_set(outdir)  then outdir  = './'
;;-    if keyword_set(test)        then test = 1 else test = 0
;;-    if keyword_set(verb)        then verb = 1 else verb = 0
;;-    if keyword_set(stop)        then stop = 1 else stop = 0
;;- 
;;-    indir  = '/redtruck/ashdown/repository/exchanges/dx11/masks/sources_hfi/dx11d_'
;;- 
;;-    ; Intensity: 1 extension at Nside 2048
;;-    print, "  >> intensity masks, nominal res"
;;- 
;;-    Nside = 2048LL   &   Npix = 12*Nside*Nside
;;-    ttail  = '_mask_sources_int_snr_5_radius_3sigma.fits.gz'
;;-    outfile = 'HFI_Mask_PointSrc_'+strtrim(Nside,2)+'_'+relnum+'.fits'
;;-    int  = replicate({f100: 0B, f143: 0B, f217: 0B, f353: 0B, f545: 0B, f857: 0B}, Npix )
;;- 
;;-    int.(0) = reorder(reform(byte((mrdfits(indir+'100'+ttail, 1, /sil, comp='gunzip')).(0)), npix), /r2n)
;;-    int.(1) = reorder(reform(byte((mrdfits(indir+'143'+ttail, 1, /sil, comp='gunzip')).(0)), npix), /r2n)
;;-    int.(2) = reorder(reform(byte((mrdfits(indir+'217'+ttail, 1, /sil, comp='gunzip')).(0)), npix), /r2n)
;;-    int.(3) = reorder(reform(byte((mrdfits(indir+'353'+ttail, 1, /sil, comp='gunzip')).(0)), npix), /r2n)
;;-    int.(4) = reorder(reform(byte((mrdfits(indir+'545'+ttail, 1, /sil, comp='gunzip')).(0)), npix), /r2n)
;;-    int.(5) = reorder(reform(byte((mrdfits(indir+'857'+ttail, 1, /sil, comp='gunzip')).(0)), npix), /r2n)
;;-    
;;-    fxbhmake, hdr, Npix,       'SRC-MASK',' Extension name', /init, /date, /extver
;;-    sxdelpar, hdr, 'EXTNAME'
;;- 
;;-    sxaddpar, hdr, 'TTYPE1 ',  'F100',  ' Mask for 100GHz'
;;-    sxaddpar, hdr, 'TTYPE2 ',  'F143',  ' Mask for 143GHz'
;;-    sxaddpar, hdr, 'TTYPE3 ',  'F217',  ' Mask for 217GHz'
;;-    sxaddpar, hdr, 'TTYPE4 ',  'F353',  ' Mask for 353GHz'
;;-    sxaddpar, hdr, 'TTYPE5 ',  'F545',  ' Mask for 545GHz'
;;-    sxaddpar, hdr, 'TTYPE6 ',  'F857',  ' Mask for 857GHz'
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TTYPE6'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;- 
;;-    sxaddpar, hdr, 'EXTNAME',  'SRC-INT ',    ' Extension name'
;;-    sxaddpar, hdr, 'PIXTYPE',  'HEALPIX'
;;-    sxaddpar, hdr, 'POLCCONV', 'COSMO',       ' Polarization convention'
;;-    sxaddpar, hdr, 'COORDSYS', 'GALACTIC',    ' Coordinate system'
;;-    sxaddpar, hdr, 'ORDERING', 'NESTED',      ' Healpix ordering'
;;-    sxaddpar, hdr, 'NSIDE',     Nside,        ' Healpix Nside'
;;-    sxaddpar, hdr, 'FIRSTPIX',  0,            ' First pixel # (0 based)'
;;-    sxaddpar, hdr, 'LASTPIX',   Npix-1,       ' Last pixel # (0 based)'
;;-    sxaddpar, hdr, 'FILENAME',  outfile,      ' FITS filename'
;;-    sxaddpar, hdr, 'PROCVER',  procver,       ' Product version'
;;- 
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------', after='PROCVER'
;;-    sxaddpar, hdr, 'COMMENT', 'Intensity point source masks; threshold is 5 snr'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT',  'For further details see Planck Explanatory Supplement at:'
;;-    sxaddpar, hdr, 'COMMENT',  '  http://wiki.cosmos.esa.int/planckpla2015'
;;-    sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;- 
;;-    print, ' >> Exporting '+outdir+outfile
;;-    mwrfits, int, outdir+outfile, hdr, /create, /sil   ; Point source Intensity masks
;;-    if verb then print, hdr
;;- 
;;-    ; Polar: nominal res (Nside 2048)
;;-    print, "  >> Polar masks, 3sigma radius, nominal res'"
;;- 
;;-    ptail  = '_mask_sources_pol_99.97pc_radius_3sigma.fits.gz'
;;-    pol  = replicate({f100: 0B, f143: 0B, f217: 0B, f353: 0B}, Npix )
;;- 
;;-    ; cut radius 3 beam sigma
;;-    pol.(0) = reorder(reform(byte((mrdfits(indir+'100'+ptail, 1, /sil, comp='gunzip')).(0)), npix), /r2n)
;;-    pol.(1) = reorder(reform(byte((mrdfits(indir+'143'+ptail, 1, /sil, comp='gunzip')).(0)), npix), /r2n)
;;-    pol.(2) = reorder(reform(byte((mrdfits(indir+'217'+ptail, 1, /sil, comp='gunzip')).(0)), npix), /r2n)
;;-    pol.(3) = reorder(reform(byte((mrdfits(indir+'353'+ptail, 1, /sil, comp='gunzip')).(0)), npix), /r2n)
;;-    
;;-    fxbhmake, pdr, Npix,       'SRC-MASK',' Extension name', /init, /date, /extver
;;-    sxdelpar, pdr, 'EXTNAME'
;;- 
;;-    sxaddpar, pdr, 'TTYPE1 ',  'F100',  ' Mask for 100GHz'
;;-    sxaddpar, pdr, 'TTYPE2 ',  'F143',  ' Mask for 143GHz'
;;-    sxaddpar, pdr, 'TTYPE3 ',  'F217',  ' Mask for 217GHz'
;;-    sxaddpar, pdr, 'TTYPE4 ',  'F353',  ' Mask for 353GHz'
;;- 
;;-    sxaddpar, pdr, 'COMMENT',  '  ' , after='TTYPE4'
;;-    sxaddpar, pdr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, pdr, 'COMMENT',  '  '
;;- 
;;-    sxaddpar, pdr, 'EXTNAME',  'SRC-POL ',    ' Extension name'
;;-    sxaddpar, pdr, 'PIXTYPE',  'HEALPIX'
;;-    sxaddpar, hdr, 'POLCCONV', 'COSMO',       ' Polarization convention'
;;-    sxaddpar, pdr, 'COORDSYS', 'GALACTIC',    ' Coordinate system'
;;-    sxaddpar, pdr, 'ORDERING', 'NESTED',      ' Healpix ordering'
;;-    sxaddpar, pdr, 'NSIDE',     Nside,        ' Healpix Nside'
;;-    sxaddpar, pdr, 'FIRSTPIX',  0,            ' First pixel # (0 based)'
;;-    sxaddpar, pdr, 'LASTPIX',   Npix-1,       ' Last pixel # (0 based)'
;;-    sxaddpar, pdr, 'FILENAME',  outfile,      ' FITS filename'
;;-    sxaddpar, pdr, 'PROCVER',  procver,       ' Product version'
;;- 
;;-    sxaddpar, pdr, 'COMMENT', '------------------------------------------------------------------------', after='PROCVER'
;;-    sxaddpar, pdr, 'COMMENT', 'Polarization point source masks, Nside=2048; '
;;-    sxaddpar, pdr, 'COMMENT', 'threshold is 99.97% detection significance, '
;;-    sxaddpar, pdr, 'COMMENT', 'cut radius = 3 beam sigma ~ 1.27 FWHM'
;;-    sxaddpar, pdr, 'COMMENT', '------------------------------------------------------------------------'
;;-    sxaddpar, pdr, 'COMMENT', 'For further details see Planck Explanatory Supplement at:'
;;-    sxaddpar, pdr, 'COMMENT', '  http://wiki.cosmos.esa.int/planckpla2015'
;;-    sxaddpar, pdr, 'COMMENT', '------------------------------------------------------------------------'
;;- 
;;-    print, ' >> Appending '+outdir+outfile
;;-    mwrfits, pol, outdir+outfile, pdr, /sil
;;- 
;;-    if stop then stop
;;- 
;;-    return
;;- end 
;;- 
;;- ;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;;- ; export CMB maps (Int and Pol) and associated products; updated for PR2
;;- ; Revised packaging - 14.apr.15: files with IQU in single extension; 
;;- ; separate files for different mission coverage (Eric's issue)
;;- ;-----------------------------------------------------------------------------
;;- 
;;- pro cmb, relnum, procver, outdir=outdir,  $
;;-            hires=hires, lores=lores, stop=stop, verb=verb
;;- 
;;- ;   on_error, 2
;;-    if n_elements(relnum) eq 0 or n_elements(procver) eq 0 then begin
;;-        print, 'SYNTAX:  cmb, relnum, procver, outdir=outdir, '
;;-        print, '         /test, /stop, /verb'
;;-        print, ''
;;-        return
;;-    endif 
;;- 
;;-    if not keyword_set(outdir)  then outdir  = '/redtruck/opsman/export/Releases/HFI_PR2_beta/'
;;-    if keyword_set(verb)        then verb = 1 else verb = 0
;;- 
;;-    root = 'COM_CMB_IQU-'
;;-    method = ['smica','commander','nilc','sevem' ]
;;-    indir  = '/redtruck/ashdown/repository/comparison/dx11_v2/'   ; for unfiltered Temp data
;;-    poldir  = '/redtruck/ashdown/repository/comparison/dx11_v2/highpass_20_40/' ; for high-pass flt
;;-    cover = ["full", "halfmission-1","halfmission-2", "year-1","year-2", "ringhalf-1","ringhalf-2"] 
;;-    ctag  = ["",     "_hm1","_hm2",  "_yr1","_yr2",   "_hr1","_hr2"]   ; coverage tags for MAJA's filenames
;;- 
;;-    ;-----------------------------------------------------------------------------
;;-    ; at Nside 1024, IQU maps in single header, and files by coverage
;;-    ;-----------------------------------------------------------------------------
;;- 
;;-    if keyword_set(hires) then begin
;;-        Nside = 2048LL  &  res  = '005a_2048'  &  lores=0
;;-        highpass = ""
;;-    endif 
;;-    if keyword_set(lores) then begin
;;-        Nside = 1024LL  &  res  = '010a_1024'  &  hires=0
;;-        highpass  = 'hp_20_40_'    ; for polar high-pass flt maps
;;-    endif 
;;- 
;;-    for n = 0,3 do begin    ; loop over methods
;;-        pipe = method[n]
;;-        Npix = 12*Nside*Nside
;;- 
;;-        for cc=0,6 do begin   ; loop over coverage
;;-           if cc eq 0 then full=1 else full=0
;;-           print, " ------------------------------------------------------------------"
;;-           print, " >> check that fits files exist for "+pipe+", coverage "+cover[cc]
;;-           spawn, 'ls -d1 '+indir+'dx11_v2_'+pipe+'_int_cmb'+ctag[cc]+'_' +res+'.fits', tfile
;;-           print, " >> read "+tfile
;;-           tsig  = mrdfits(tfile,  1, /sil)
;;-           if lores then begin
;;-               spawn, 'ls -d1 '+poldir+'dx11_v2_'+pipe+'_pol_case1_cmb'+ctag[cc]+'_' +highpass+res+'.fits', pfile
;;-               print, " >> read "+pfile
;;-               psig  = mrdfits(pfile,  1, /sil) 
;;-           endif 
;;- 
;;-           if lores then outfile = root+pipe+'_'+strtrim(Nside,2)+'_'+relnum+'_'+cover[cc]+'.fits' else $
;;-             outfile = root+pipe+'-field-Int_'+strtrim(Nside,2)+'_'+relnum+'_'+cover[cc]+'.fits'
;;- 
;;-           print, " >> Build I(QU) structure and fill it"
;;-           if lores then begin
;;-               if full then str = replicate({I_stokes:0.0, Q_stokes:0.0, U_stokes:0.0, tmask:0b, pmask:0b}, Npix) else $
;;-                 str = replicate({I_stokes:0.0, Q_stokes:0.0, U_stokes:0.0}, Npix)
;;-           endif else begin
;;-               if full then str = replicate({I_stokes:0.0, tmask:0b}, Npix) else $
;;-                 str = replicate({I_stokes:0.0}, Npix)
;;-           endelse               
;;- 
;;-           str.i_stokes  = reorder(reform(float(tsig.(0)), npix), /r2n)
;;-           if lores then begin
;;-               str.q_stokes  = reorder(reform(float(psig.(0)) , npix), /r2n)
;;-               str.u_stokes  = reorder(reform(float(psig.(1)) , npix), /r2n)
;;-           endif 
;;-           if full then begin
;;-               print, " ## 'full' coverage ... add confidence mask(s) to structure"
;;-               spawn, 'ls -d1 '+indir+'dx11_v2_'+pipe+'_int_mask_'+res+'.fits', tmask
;;-               print, " >> file "+tmask
;;-               tmask = mrdfits(tmask,   1, /sil) 
;;-               str.tmask  = reorder(reform(byte(tmask.(0)) , npix), /r2n)
;;-               if lores then begin   ; add also pol mask
;;-                   spawn, 'ls -d1 '+indir+'dx11_v2_'+pipe+'_pol_mask_'+res+'.fits', pmask
;;-                   print, " >> file "+pmask
;;-                   pmask = mrdfits(pmask,   1, /sil) 
;;-                   str.pmask  = reorder(reform(byte(pmask.(0)) , npix), /r2n)
;;-               endif 
;;-           endif 
;;- 
;;-           print, " >> Build primary header"
;;-           mkhdr, phdr, '', /extend
;;-           if cc eq 0 then sxaddpar, phdr, 'NUMEXT', 2, ' Number of extensions' else $
;;-             sxaddpar, phdr, 'NUMEXT', 1, ' Number of extensions', after='DATE'
;;-           sxaddpar, phdr, 'FILENAME',  outfile,    ' FITS filename', after='NUMEXT'
;;-           sxaddpar, phdr, 'COMMENT',  '  '  , after='FILENAME'
;;-           sxaddpar, phdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-           sxaddpar, phdr, 'COMMENT',  'CMB products from '+pipe+' component separation method'
;;-           sxaddpar, phdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-           sxaddpar, phdr, 'COMMENT',  'For further details see Planck Explanatory Supplement at:'
;;-           sxaddpar, phdr, 'COMMENT',  '  http://wiki.cosmos.esa.int/planckpla2015'
;;-           sxaddpar, phdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-           mwrfits, qewr, outdir+outfile, phdr, /create
;;-           
;;-           ; -----------------------------------------------------------------------------
;;-           print, " >> Build header for I(QU) maps"
;;- 
;;-           fxbhmake, hdr, Npix,       'changeme', /init, /date, /extver
;;-           sxdelpar, hdr, 'EXTNAME'
;;- 
;;-           if lores then begin
;;-               if full then sxaddpar, hdr, 'COMMENT',  '  '  , after='TFORM5' else $
;;-                 sxaddpar, hdr, 'COMMENT',  '  '  , after='TFORM3'
;;-           endif else begin 
;;-               if full then sxaddpar, hdr, 'COMMENT',  '  '  , after='TFORM2' else $
;;-                 sxaddpar, hdr, 'COMMENT', ' ', after='TFORM1'
;;-           end
;;- 
;;-           sxaddpar, hdr, 'COMMENT',  ' *** Column units ***  '
;;-           sxaddpar, hdr, 'COMMENT',  '  '  
;;-           sxaddpar, hdr, 'TUNIT1', 'K_CMB',  ' map units'
;;-           if lores then begin
;;-               sxaddpar, hdr, 'TUNIT2' , 'K_CMB ',  ' map units'
;;-               sxaddpar, hdr, 'TUNIT3' , 'K_CMB ',  ' map units'
;;-               if full then begin
;;-                   sxaddpar, hdr, 'TUNIT4' , ' ' , ' no units'
;;-                   sxaddpar, hdr, 'TUNIT5' , ' ' , ' no units'
;;-                   sxaddpar, hdr, 'COMMENT',  '  '  , after='TUNIT5'
;;-               endif else begin
;;-                   sxaddpar, hdr, 'COMMENT',  '  '  , after='TUNIT3'
;;-               endelse 
;;-           endif else begin
;;-               if full then begin
;;-                   sxaddpar, hdr, 'TUNIT2' , ' ' , ' no units'
;;-                   sxaddpar, hdr, 'COMMENT',  '  '  , after='TUNIT2'
;;-               endif else begin
;;-                   sxaddpar, hdr, 'COMMENT',  '  '  , after='TUNIT1'
;;-               end  
;;-           endelse 
;;-           
;;-           sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-           sxaddpar, hdr, 'COMMENT',  '  '
;;- 
;;-           sxaddpar, hdr, 'EXTNAME',  'COMP-MAP',    ' Extension name'
;;-           sxaddpar, hdr, 'AST-COMP', 'CMB',         ' Component'
;;-           if lores then sxaddpar, hdr, 'RESOLN', 10, ' arcmin' else $
;;-             sxaddpar, hdr, 'RESOLN', 5, 'arcmin'
;;-           sxaddpar, hdr, 'PIXTYPE',  'HEALPIX'
;;-           sxaddpar, hdr, 'POLCCONV', 'COSMO',       ' Polarization convention'
;;-           sxaddpar, hdr, 'COORDSYS', 'GALACTIC',    ' Coordinate system'
;;-           sxaddpar, hdr, 'ORDERING', 'NESTED',      ' Healpix ordering'
;;-           sxaddpar, hdr, 'NSIDE',     Nside,        ' Healpix Nside'
;;-           sxaddpar, hdr, 'FIRSTPIX',  0
;;-           sxaddpar, hdr, 'LASTPIX',   Npix-1
;;-           sxaddpar, hdr, 'FILENAME',  outfile,      ' FITS filename'
;;-           sxaddpar, hdr, 'BAD_DATA',  -1.63750E+30, ' HEALPIX bad pixel value'
;;-           sxaddpar, hdr, 'METHOD',   pipe,          ' Separation method'
;;-           sxaddpar, hdr, 'PROCVER',  procver,       ' Product version'
;;- 
;;-           sxaddpar, hdr, 'COMMENT',  '  ' , after='PROCVER'
;;-           sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-           sxaddpar, hdr, 'COMMENT',  'CMB products from '+pipe+", coverage "+cover[cc]
;;-           sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-           sxaddpar, hdr, 'COMMENT',  'For further details see Planck Explanatory Supplement at:'
;;-           sxaddpar, hdr, 'COMMENT',  '  http://wiki.cosmos.esa.int/planckpla2015'
;;-           sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-        
;;-           print, ' >> Writing I(QU) maps HDU for '+pipe+', '+cover[cc]+': '+outdir+outfile
;;-           mwrfits, str, outdir+outfile, hdr, /sil ; CMB IQU maps
;;-           if verb then print, hdr
;;- 
;;-           if full then begin
;;-              ;------------------------------------------------------------------
;;-              print, " ## Build header for beam transfer function(s)"
;;-              ;-----------------------------------------------------------------
;;-              
;;-              spawn, 'ls -d1 '+indir+'dx11_v2_'+pipe+'_int_beam_'+res+'.fits', tbeam
;;-              print, ' ## Int beam file: '+tbeam
;;-              itf = mrdfits(tbeam, 1, /sil)       ; intensity beam transfer function
;;-              if lores then begin   ; get also polar beam transfer function
;;-                  spawn, 'ls -d1 '+indir+'dx11_v2_'+pipe+'_pol_beam_'+res+'.fits', pbeam
;;-                  print, ' ## Pol beam file: '+tbeam
;;-                  ptf = mrdfits(pbeam, 1, /sil)
;;-              endif
;;-              nell = n_elements(itf)
;;- 
;;-              ; build structure and fill it in
;;-              if lores then beam = replicate({int_beam:0.0, pol_beam:0.0}, nell) $
;;-                else beam = replicate({int_beam:0.0}, nell)
;;-              beam.int_beam  = float(itf.(0))
;;-              if lores then begin
;;-                  pell = n_elements(ptf)
;;-                  beam[0:pell-1].pol_beam  = float(ptf.(0))
;;-              endif 
;;-              
;;-              ; build header
;;-              fxbhmake, bdr, nell, 'ChangeMe', /init, /date, /extver
;;-              sxdelpar, bdr, 'EXTNAME'
;;-              sxaddpar, hdr, 'TTYPE1' , 'INT_BEAM',  ' Beam Transfer function'
;;-              if lores then begin
;;-                  sxaddpar, hdr, 'TTYPE2' , 'POL_BEAM',  ' Beam Transfer function'
;;-                  sxaddpar, bdr, 'COMMENT',  '  ' , after='TFORM2'
;;-              endif else begin
;;-                  sxaddpar, bdr, 'COMMENT',  '  ' , after='TFORM1'
;;-              endelse 
;;-              sxaddpar, bdr, 'COMMENT',  ' *** Column units *** '
;;-              sxaddpar, bdr, 'COMMENT',  '  '
;;-              sxaddpar, bdr, 'TUNIT1',   'none'
;;-              if lores then begin
;;-                  sxaddpar, bdr, 'TUNIT2',   'none'
;;-                  sxaddpar, bdr, 'COMMENT',  '  ' , after='TUNIT2'
;;-              endif else begin
;;-                  sxaddpar, bdr, 'COMMENT',  '  ' , after='TUNIT1'
;;-              endelse 
;;- 
;;-              sxaddpar, bdr, 'COMMENT',  ' *** Planck params *** '
;;-              sxaddpar, bdr, 'COMMENT',  '  '
;;-              
;;-              sxaddpar, bdr, 'EXTNAME',  'BEAMTF',    ' Beam Transfer Function'
;;-              sxaddpar, bdr, 'LMIN',     0
;;-              sxaddpar, bdr, 'LMAX_I',     nell-1, ' L_max for Intnsity beam TF'
;;-              if lores then  sxaddpar, bdr, 'LMAX_P', pell-1, ' L_max for Polar beam TF'
;;-              sxaddpar, bdr, 'FILENAME',  outfile,    ' FITS filename'
;;-              sxaddpar, bdr, 'PROCVER',  procver,     ' Product version'
;;-              sxaddpar, bdr, 'METHOD',   pipe,          ' Separation method'
;;-              sxaddpar, bdr, 'COMMENT',  '  ' , after='PROCVER'
;;-              sxaddpar, bdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-              sxaddpar, bdr, 'COMMENT',  'CMB products from '+pipe+': beam transfer function'
;;-              sxaddpar, bdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-              sxaddpar, bdr, 'COMMENT',  'For further details see Planck Explanatory Supplement at:'
;;-              sxaddpar, bdr, 'COMMENT',  '  http://wiki.cosmos.esa.int/planckpla2015'
;;-              sxaddpar, bdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-              
;;-              print, ' >> Appending beam window fns HDU: '+outdir+outfile
;;-              mwrfits, beam, outdir+outfile, bdr, /sil ; - append beam window function
;;-              if verb then print, bdr
;;-          endif                   ; extension for beam
;;-      end 
;;-       if defined(stop) then stop
;;-       print, ""
;;-    end
;;- 
;;- 
;;-    return
;;-    ; to check:
;;-    mm=mrdfits('/redtruck/opsman/export/Releases/HFI_PR2/COM_CMB_IQU-smica_1024_R2.02_full.fits', 1, h)
;;-    ii=mrdfits('/redtruck/ashdown/repository/comparison/dx11_v2/dx11_v2_smica_int_cmb_010a_1024.fits', 1)
;;-    pp=mrdfits('/redtruck/ashdown/repository/comparison/dx11_v2/highpass_20_40/dx11_v2_smica_pol_case1_cmb_hp_20_40_010a_1024.fits', 1)
;;-    print, minmax(reorder( mm.(0), /n2r) - reform(ii.(0), 12*1024L*1024))
;;-    print, minmax(reorder( mm.(1), /n2r) - reform(pp.(0), 12*1024L*1024))
;;-    print, minmax(reorder( mm.(2), /n2r) - reform(pp.(1), 12*1024L*1024))
;;- end
;;-  
;;- 
;;- 
;;- ;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;;- ; export cmb kinematic quadrupole residual maps
;;- ;-----------------------------------------------------------------------------
;;- 
;;- pro kqres, relnum, procver, outdir=outdir, $
;;-             ingrp=ingrp, test=test, stop=stop, verb=verb
;;- 
;;- ;   on_error, 2
;;- 
;;-    if not defined(relnum)  then begin
;;-        print, 'SYNTAX: '
;;-        print, '  cmbmask, relnum, procver, outdir=outdir, ingrp=ingrp,  $'
;;-        print, '          /test, /stop'    
;;-        return
;;-    endif
;;- 
;;-    if not keyword_set(outdir)  then outdir  = '/redtruck/opsman/export/Releases/HFI_PR2_beta/'
;;-    if keyword_set(test)        then test  = 1 else test  = 0
;;-    if keyword_set(verb)        then verb = 1 else verb = 0
;;- 
;;-    root = 'COM_CMB_IQU-kq-resid-'
;;-    method = ['smica','commander','nilc','sevem' ]
;;-    indir  = '/redtruck/ashdown/compsep/comparison/dx11_v2/int/kq_residuals/'
;;-    Nside = 2048LL  &   Npix  = 12L*Nside*Nside    &  res  = '005a_2048' 
;;- 
;;-    for n = 0,3 do begin    ; loop over methods
;;-       pipe = method[n]
;;-       spawn, 'ls -d1 '+indir+'dx11_v2_'+pipe+'_int_kq_residuals_' +res+'.fits', file
;;-       if strlen(file[0]) lt 10 then begin
;;-           print, " No file found for "+pipe
;;-       ;    return
;;-       endif else begin
;;-          print, " working on "+file
;;-          m = mrdfits(file[0], 1, h,/sil)
;;- 
;;-          units = sxpar(h, 'TUNIT1')
;;-          str = replicate({Intensity:0.0}, Npix)
;;-          str.Intensity = reorder(reform(m.(0), npix), /r2n)
;;-          outfile = root+pipe+'-field-Int_'+strtrim(Nside,2)+'_'+relnum+'.fits'
;;-       
;;-          print, " >> Build primary header"
;;-          mkhdr, phdr, '', /extend
;;-          sxaddpar, phdr, 'NUMEXT', 1, ' Number of extensions', after='DATE'
;;-          sxaddpar, phdr, 'FILENAME',  outfile,    ' FITS filename', after='NUMEXT'
;;-          sxaddpar, phdr, 'COMMENT',  '  '  , after='FILENAME'
;;-          sxaddpar, phdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-          sxaddpar, phdr, 'COMMENT',  'CMB products from '+pipe+' component separation method'
;;-          sxaddpar, phdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-          sxaddpar, phdr, 'COMMENT',  'For further details see Planck Explanatory Supplement at:'
;;-          sxaddpar, phdr, 'COMMENT',  '  http://wiki.cosmos.esa.int/planckpla2015'
;;-          sxaddpar, phdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-          mwrfits, qewr, outdir+outfile, phdr, /create
;;-          
;;-          print, " >> Build header for kq-resid maps"
;;-          fxbhmake, hdr, Npix,       'changeme', /init, /date, /extver
;;-          sxdelpar, hdr, 'EXTNAME'
;;-          sxaddpar, hdr, 'COMMENT',  '  ' , after='TFORM1'
;;-          sxaddpar, hdr, 'COMMENT',  ' *** Column units ***  '
;;-          sxaddpar, hdr, 'COMMENT',  '  '  
;;-          sxaddpar, hdr, 'TUNIT1', 'K_CMB',  ' map units'
;;-          sxaddpar, hdr, 'EXTNAME',  'COMP-MAP',    ' Extension name'
;;-          sxaddpar, hdr, 'AST-COMP', 'KQ-RESID',    ' Component'
;;-          sxaddpar, hdr, 'RESOLN', 5, 'arcmin'
;;-          sxaddpar, hdr, 'PIXTYPE',  'HEALPIX'
;;-          ;sxaddpar, hdr, 'POLCCONV', 'COSMO',       ' Polarization convention'
;;-          sxaddpar, hdr, 'COORDSYS', 'GALACTIC',    ' Coordinate system'
;;-          sxaddpar, hdr, 'ORDERING', 'NESTED',      ' Healpix ordering'
;;-          sxaddpar, hdr, 'NSIDE',     Nside,        ' Healpix Nside'
;;-          sxaddpar, hdr, 'FIRSTPIX',  0
;;-          sxaddpar, hdr, 'LASTPIX',   Npix-1
;;-          sxaddpar, hdr, 'FILENAME',  outfile,      ' FITS filename'
;;-          sxaddpar, hdr, 'BAD_DATA',  -1.63750E+30, ' HEALPIX bad pixel value'
;;-          sxaddpar, hdr, 'METHOD',   pipe,          ' Separation method'
;;-          sxaddpar, hdr, 'PROCVER',  procver,       ' Product version'
;;-          
;;-          sxaddpar, hdr, 'COMMENT',  '  ' , after='PROCVER'
;;-          sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-          sxaddpar, hdr, 'COMMENT',  'CMB kinetic quadrupole residual map from '+pipe
;;-          sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-          sxaddpar, hdr, 'COMMENT',  'For further details see Planck Explanatory Supplement at:'
;;-          sxaddpar, hdr, 'COMMENT',  '  http://wiki.cosmos.esa.int/planckpla2015'
;;-          sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-          
;;-          print, ' >> Writing kq-resid map HDU for '+pipe+': '+outdir+outfile
;;-          mwrfits, str, outdir+outfile, hdr, /sil ; CMB IQU maps
;;-          if verb then print, hdr
;;-      endelse 
;;-    endfor 
;;-    return
;;- end
;;- 
;;- 
;;- ;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;;- ; export cmb-subtraced foreground maps
;;- ;-----------------------------------------------------------------------------
;;- 
;;- pro fgs, relnum, procver, outdir=outdir, $
;;-             ingrp=ingrp, test=test, stop=stop, verb=verb
;;- 
;;- ;   on_error, 2
;;- 
;;-    if not defined(relnum)  then begin
;;-        print, 'SYNTAX: '
;;-        print, '  fgs, relnum, procver, outdir=outdir, ingrp=ingrp,  $'
;;-        print, '          /test, /stop'    
;;-        return
;;-    endif
;;- 
;;-    if not keyword_set(outdir)  then outdir  = '/redtruck/opsman/export/Releases/HFI_PR2_beta/'
;;-    if keyword_set(test)        then test  = 1 else test  = 0
;;-    if keyword_set(verb)        then verb = 1 else verb = 0
;;- 
;;-    root = 'FI_Foregrounds-'
;;-    pipes = ['smica','commander','nilc','sevem' ]
;;-    indir = '/redtruck/ashdown/compsep/comparison/dx11_v2/int/cmb_subtracted/'
;;-    freqs = ['030','044','070','100', '143', '217', '353', '545', '857']
;;- 
;;-    Nside = 2048LL  &   Npix  = 12L*Nside*Nside
;;-    Lside = 1024LL  &   Lpix  = 12L*Lside*Lside
;;-    
;;-    Lstr = replicate({C030:0.0 , C044:0.0, C070:0.0}, Lpix)
;;-    Hstr = replicate({C100:0.0 , C143:0.0, C217:0.0, C353:0.0, C545:0.0, C857:0.0}, Npix)
;;- 
;;-    for n = 0,3 do begin    ; loop over methods
;;-       pipe = pipes[n]
;;-       for f = 0,8 do begin   ; loop over freqs
;;-           spawn, 'ls -d1 '+indir+'dx11_v2_'+pipe+'_int_' +freqs[f]+'_cmb_subtracted.fits', file
;;-           if strlen(file[0]) lt 10 then begin
;;-               print, " No file found for "+pipe
;;-               return
;;-           endif else begin
;;-               print, " reading on "+file
;;-               m = mrdfits(file[0], 1, h,/sil)
;;-               
;;-               case freqs[f] of
;;-                   '030': Lstr.c030 = reorder(reform(m.(0), lpix),/r2n)
;;-                   '044': Lstr.c044 = reorder(reform(m.(0), lpix),/r2n)
;;-                   '070': Lstr.c070 = reorder(reform(m.(0), lpix),/r2n)
;;-                   '100': Hstr.c100 = reorder(reform(m.(0), npix),/r2n)
;;-                   '143': Hstr.c143 = reorder(reform(m.(0), npix),/r2n)
;;-                   '217': Hstr.c217 = reorder(reform(m.(0), npix),/r2n)
;;-                   '353': Hstr.c353 = reorder(reform(m.(0), npix),/r2n)
;;-                   '545': Hstr.c545 = reorder(reform(m.(0), npix),/r2n)
;;-                   '857': Hstr.c857 = reorder(reform(m.(0), npix),/r2n)
;;-               endcase 
;;-           endelse 
;;-       endfor 
;;-          
;;-       outfile = "L"+root+pipe+'_'+strtrim(lside,2)+'_'+relnum+'.fits'
;;-    
;;-       ; primary header works for both output files
;;-       print, " >> Build primary header"
;;-       mkhdr, phdr, '', /extend
;;-       sxaddpar, phdr, 'NUMEXT', 1, ' Number of extensions', after='DATE'
;;-       sxaddpar, phdr, 'FILENAME',  outfile,    ' FITS filename', after='NUMEXT'
;;-       sxaddpar, phdr, 'COMMENT',  '  '  , after='FILENAME'
;;-       sxaddpar, phdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-       sxaddpar, phdr, 'COMMENT',  'For further details see Planck Explanatory Supplement at:'
;;-       sxaddpar, phdr, 'COMMENT',  '  http://wiki.cosmos.esa.int/planckpla2015'
;;-       sxaddpar, phdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-       mwrfits, qewr, outdir+outfile, phdr, /create
;;-          
;;-       print, " >> Build header for LFI maps at Nside 1024"
;;-       fxbhmake, ldr, Npix,       'changeme', /init, /date, /extver
;;-       sxdelpar, ldr, 'EXTNAME'
;;-       sxaddpar, ldr, 'TTYPE1 ',  'C030',  ' Foregrounds for 030GHz'
;;-       sxaddpar, ldr, 'TTYPE2 ',  'C044',  ' Foregrounds for 044GHz'
;;-       sxaddpar, ldr, 'TTYPE3 ',  'C070',  ' Foregrounds for 070GHz'
;;- 
;;-       sxaddpar, ldr, 'COMMENT',  '  ' , after='TTYPE3'
;;-       sxaddpar, ldr, 'COMMENT',  ' *** Column units ***  '
;;-       sxaddpar, ldr, 'COMMENT',  '  '  
;;-       sxaddpar, ldr, 'TUNIT1', 'K_CMB',  ' map units'
;;-       sxaddpar, ldr, 'TUNIT2', 'K_CMB',  ' map units'
;;-       sxaddpar, ldr, 'TUNIT3', 'K_CMB',  ' map units'
;;- 
;;-       sxaddpar, ldr, 'COMMENT',  '  ' , after='TUNIT3'
;;-       sxaddpar, ldr, 'COMMENT',  ' *** Planck params *** '
;;-       sxaddpar, ldr, 'COMMENT',  '  '
;;-       sxaddpar, ldr, 'EXTNAME',  'LFI-RESID',   ' Extension name'
;;-       sxaddpar, ldr, 'AST-COMP', 'FOREGDS',     ' Component'
;;- ;      sxaddpar, ldr, 'RESOLN', 5, 'arcmin'
;;-       sxaddpar, ldr, 'PIXTYPE',  'HEALPIX'
;;-       sxaddpar, ldr, 'COORDSYS', 'GALACTIC',    ' Coordinate system'
;;-       sxaddpar, ldr, 'ORDERING', 'NESTED',      ' Healpix ordering'
;;-       sxaddpar, ldr, 'NSIDE',     Lside,        ' Healpix Nside'
;;-       sxaddpar, ldr, 'FIRSTPIX',  0
;;-       sxaddpar, ldr, 'LASTPIX',   Lpix-1
;;-       sxaddpar, ldr, 'FILENAME',  outfile,      ' FITS filename'
;;-       sxaddpar, ldr, 'BAD_DATA',  -1.63750E+30, ' HEALPIX bad pixel value'
;;-       sxaddpar, ldr, 'METHOD',   pipe,          ' Separation method'
;;-       sxaddpar, ldr, 'PROCVER',  procver,       ' Product version'
;;-    
;;-       sxaddpar, ldr, 'COMMENT',  '  ' , after='PROCVER'
;;-       sxaddpar, ldr, 'COMMENT',  '------------------------------------------------------------------------'
;;-       sxaddpar, ldr, 'COMMENT',  'LFI foregrounds map from '+pipe
;;-       sxaddpar, ldr, 'COMMENT',  '------------------------------------------------------------------------'
;;-       sxaddpar, ldr, 'COMMENT',  'For further details see Planck Explanatory Supplement at:'
;;-       sxaddpar, ldr, 'COMMENT',  '  http://wiki.cosmos.esa.int/planckpla2015'
;;-       sxaddpar, ldr, 'COMMENT',  '------------------------------------------------------------------------'
;;-    
;;-       print, ' >> Writing LFI forgrounds for '+pipe+': '+outdir+outfile
;;-       mwrfits, Lstr, outdir+outfile, ldr, /sil ; CMB IQU maps
;;- 
;;-       outfile = "H"+root+pipe+'_'+strtrim(Nside,2)+'_'+relnum+'.fits'
;;-       mwrfits, qewr, outdir+outfile, phdr, /create
;;- 
;;-       print, " >> Build header for HFI maps at Nside 2048"
;;-       fxbhmake, hdr, Npix,       'changeme', /init, /date, /extver
;;-       sxdelpar, hdr, 'EXTNAME'
;;-       sxaddpar, hdr, 'TTYPE1 ',  'C100',  ' Foregrounds for 100GHz'
;;-       sxaddpar, hdr, 'TTYPE2 ',  'C143',  ' Foregrounds for 143GHz'
;;-       sxaddpar, hdr, 'TTYPE3 ',  'C217',  ' Foregrounds for 217GHz'
;;-       sxaddpar, hdr, 'TTYPE4 ',  'C353',  ' Foregrounds for 353GHz'
;;-       sxaddpar, hdr, 'TTYPE5 ',  'C545',  ' Foregrounds for 545GHz'
;;-       sxaddpar, hdr, 'TTYPE6 ',  'C857',  ' Foregrounds for 857GHz'
;;- 
;;-       sxaddpar, hdr, 'COMMENT',  '  ' , after='TTYPE6'
;;-       sxaddpar, hdr, 'COMMENT',  ' *** Column units ***  '
;;-       sxaddpar, hdr, 'COMMENT',  '  '  
;;-       sxaddpar, hdr, 'TUNIT1', 'K_CMB',  ' map units'
;;-       sxaddpar, hdr, 'TUNIT2', 'K_CMB',  ' map units'
;;-       sxaddpar, hdr, 'TUNIT3', 'K_CMB',  ' map units'
;;-       sxaddpar, hdr, 'TUNIT4', 'K_CMB',  ' map units'
;;-       sxaddpar, hdr, 'TUNIT5', 'MJY/SR', ' map units'
;;-       sxaddpar, hdr, 'TUNIT6', 'MJY/SR', ' map units'
;;- 
;;-       sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT6'
;;-       sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-       sxaddpar, hdr, 'COMMENT',  '  '
;;-       sxaddpar, hdr, 'EXTNAME',  'HFI-RESID',   ' Extension name'
;;-       sxaddpar, hdr, 'AST-COMP', 'FOREGDS',     ' Component'
;;- ;      sxaddpar, hdr, 'RESOLN', 5, 'arcmin'
;;-       sxaddpar, hdr, 'PIXTYPE',  'HEALPIX'
;;-       sxaddpar, hdr, 'COORDSYS', 'GALACTIC',    ' Coordinate system'
;;-       sxaddpar, hdr, 'ORDERING', 'NESTED',      ' Healpix ordering'
;;-       sxaddpar, hdr, 'NSIDE',     Nside,        ' Healpix Nside'
;;-       sxaddpar, hdr, 'FIRSTPIX',  0
;;-       sxaddpar, hdr, 'LASTPIX',   Npix-1
;;-       sxaddpar, hdr, 'FILENAME',  outfile,      ' FITS filename'
;;-       sxaddpar, hdr, 'BAD_DATA',  -1.63750E+30, ' HEALPIX bad pixel value'
;;-       sxaddpar, hdr, 'METHOD',   pipe,          ' Separation method'
;;-       sxaddpar, hdr, 'PROCVER',  procver,       ' Product version'
;;-    
;;-       sxaddpar, hdr, 'COMMENT',  '  ' , after='PROCVER'
;;-       sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-       sxaddpar, hdr, 'COMMENT',  'HFI foregrounds map from '+pipe
;;-       sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-       sxaddpar, hdr, 'COMMENT',  'For further details see Planck Explanatory Supplement at:'
;;-       sxaddpar, hdr, 'COMMENT',  '  http://wiki.cosmos.esa.int/planckpla2015'
;;-       sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-    
;;-       print, ' >> Writing HFI foregrounds for '+pipe+': '+outdir+outfile
;;-       mwrfits, Hstr, outdir+outfile, hdr, /sil ; CMB IQU maps
;;-    endfor 
;;- 
;;- 
;;-    if verb then print, hdr
;;-    return
;;- end
;;- 
;;- 
;;- ;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;;- ; export cmb common and missing pix masks: new for PR2
;;- ;-----------------------------------------------------------------------------
;;- 
;;- pro cmbmask, relnum, procver, outdir=outdir, $
;;-             ingrp=ingrp, test=test, stop=stop, verb=verb
;;- 
;;- ;   on_error, 2
;;- 
;;-    if not defined(relnum)  then begin
;;-        print, 'SYNTAX: '
;;-        print, '  cmbmask, relnum, procver, outdir=outdir, ingrp=ingrp,  $'
;;-        print, '          /test, /stop'    
;;-        return
;;-    endif
;;- 
;;-    if not keyword_set(outdir)  then outdir  = '/redtruck/opsman/export/Releases/HFI_PR2_beta/'
;;-    if keyword_set(test)        then test  = 1 else test  = 0
;;-    if keyword_set(verb)        then verb = 1 else verb = 0
;;- 
;;-    indir  = '/redtruck/ashdown/repository/comparison/dx11_v2/'   
;;-    Nside = 2048LL  &   Npix  = 12L*Nside*Nside
;;- 
;;-    ; Common Intensity mask
;;- 
;;-    res = '005a_2048'  &   Nside = 2048LL  &   Npix  = 12L*Nside*Nside
;;-    outfile = 'COM_CMB_IQU-common-field-MaskInt_'+strtrim(Nside,2)+'_'+relnum+'.fits'
;;- 
;;-    print, " >> Read Temperature common masks"
;;-    spawn, 'ls -d1 '+indir+'dx11_v2_common_int_mask_'           +res+'.fits', ut78file   ; UT78
;;-    spawn, 'ls -d1 '+indir+'dx11_v2_common_int_mask_extended_'  +res+'.fits', ut76file   ; UT76, preferred
;;-    spawn, 'ls -d1 '+indir+'dx11_v2_common_int_mask_misspix_hm_'+res+'.fits', cthmask
;;-    spawn, 'ls -d1 '+indir+'dx11_v2_common_int_mask_misspix_yr_'+res+'.fits', ctymask
;;-    print, " >> files are: "
;;-    print, "    - "+ut78file
;;-    print, "    - "+ut76file
;;-    print, "    - "+cthmask
;;-    print, "    - "+ctymask
;;- 
;;-    read_fits_map, ut78file, ut78 , 1, ordering=order 
;;-    read_fits_map, ut76file, ut76 , 1
;;-    read_fits_map, cthmask, hmis , 1
;;-    read_fits_map, ctymask, year , 1
;;- 
;;-    print, " >> Build structure and fill in "
;;-    msk = replicate({UT78:0B, UT76:0B, HMIS:0B, YEAR:0B}, npix)
;;-    
;;-    msk.ut78  = reorder(byte(ut78), /r2n)
;;-    msk.ut76  = reorder(byte(ut76), /r2n)
;;-    msk.hmis  = reorder(byte(hmis), /r2n)
;;-    msk.year  = reorder(byte(year), /r2n)
;;- 
;;-    ;-----------------------------------------------------------------------------
;;-    print, " >> Build header for Temp masks"
;;-    ;-----------------------------------------------------------------------------
;;-    fxbhmake, hdr, Npix,       'changeme', /init, /date, /extver
;;-    sxdelpar, hdr, 'EXTNAME'
;;-    
;;-    sxaddpar, hdr, 'TTYPE1' , 'UT78    ',  ' Common mask UT78'
;;-    sxaddpar, hdr, 'TTYPE2' , 'UT76    ',  ' Common mask UT76'
;;-    sxaddpar, hdr, 'TTYPE3' , 'HMIS    ',  ' Missing pixels, half-mission'
;;-    sxaddpar, hdr, 'TTYPE4' , 'YEAR    ',  ' Missing pixels, year'
;;-    
;;-    sxaddpar, hdr, 'COMMENT',  '  '  , after='TFORM4'       
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, hdr, 'TUNIT1',   '  '
;;-    sxaddpar, hdr, 'TUNIT2',   '  '
;;-    sxaddpar, hdr, 'TUNIT3',   '  '
;;-    sxaddpar, hdr, 'TUNIT4',   '  '
;;-    
;;-    sxaddpar, hdr, 'COMMENT',  '  '  , after='TUNIT4'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, hdr, 'EXTNAME',  'MASK-INT',    ' Extension name'
;;-    sxaddpar, hdr, 'AST-COMP', 'CMB',         ' Component'
;;-    sxaddpar, hdr, 'PIXTYPE',  'HEALPIX'
;;-    sxaddpar, hdr, 'POLCCONV', 'COSMO',       ' Polarization convention'
;;-    sxaddpar, hdr, 'COORDSYS', 'GALACTIC',    ' Coordinate system'
;;-    sxaddpar, hdr, 'ORDERING', 'NESTED',      ' Healpix ordering'
;;-    sxaddpar, hdr, 'NSIDE',     Nside,        ' Healpix Nside'
;;-    sxaddpar, hdr, 'FIRSTPIX',  0
;;-    sxaddpar, hdr, 'LASTPIX',   Npix-1
;;-    sxaddpar, hdr, 'FILENAME',  outfile,      ' FITS filename'
;;-    sxaddpar, hdr, 'BAD_DATA',  -1.63750E+30, ' HEALPIX bad pixel value'
;;-    sxaddpar, hdr, 'PROCVER',  procver,       ' Product version'
;;-    
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='PROCVER'
;;-    sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT',  'CMB Common masks for Intensity '
;;-    sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT',  'For further details see Planck Explanatory Supplement at:'
;;-    sxaddpar, hdr, 'COMMENT',  '  http://wiki.cosmos.esa.int/planckpla2015'
;;-    sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-    
;;-    print, ' >> Writing Temp Common Masks HDU:           '+outdir+outfile
;;-    mwrfits, msk, outdir+outfile, hdr, /create ; 
;;-    if verb then print, hdr
;;-    
;;- 
;;-    ; Common Polar mask
;;- 
;;-    res = '010a_1024'  &   Nside = 1024LL   &   Npix = 12*Nside*Nside
;;-    outfile = 'COM_CMB_IQU-common-field-MaskPol_'+strtrim(Nside,2)+'_'+relnum+'.fits'
;;- 
;;-    print, " >> Read Polar common masks"
;;-    spawn, 'ls -d1 '+indir+'dx11_v2_common_pol_mask_'           +res+'.fits', up78file
;;-    spawn, 'ls -d1 '+indir+'dx11_v2_common_pol_mask_extended1_' +res+'.fits', upa77file
;;-    spawn, 'ls -d1 '+indir+'dx11_v2_common_pol_mask_new_'       +res+'.fits', upb77file
;;-    spawn, 'ls -d1 '+indir+'dx11_v2_common_pol_mask_misspix_hm_'+res+'.fits', cphmask
;;-    spawn, 'ls -d1 '+indir+'dx11_v2_common_pol_mask_misspix_yr_'+res+'.fits', cpymask
;;-    print, " >> files are: "
;;-    print, "    - "+up78file
;;-    print, "    - "+upa77file
;;-    print, "    - "+upb77file
;;-    print, "    - "+cphmask
;;-    print, "    - "+cpymask
;;- 
;;- 
;;-    read_fits_map, up78file,  up78  , 1, ordering=order 
;;-    read_fits_map, upa77file, upa77 , 1
;;-    read_fits_map, upb77file, upb77 , 1
;;-    read_fits_map, cphmask, hmis , 1
;;-    read_fits_map, cpymask, year , 1
;;- 
;;-    print, " >> Build structure and fill in "
;;-    msk = replicate({UP78:0B, UPA77:0B, UPB77:0B, HMIS:0B, YEAR:0B}, npix)
;;-    
;;-    msk.up78   = reorder(byte(up78),  /r2n)
;;-    msk.upa77  = reorder(byte(upa77), /r2n)
;;-    msk.upb77  = reorder(byte(upb77), /r2n)
;;-    msk.hmis   = reorder(byte(hmis),  /r2n)
;;-    msk.year   = reorder(byte(year),  /r2n)
;;- 
;;-    ;-----------------------------------------------------------------------------
;;-    print, " >> Build header for Polar masks"
;;-    ;-----------------------------------------------------------------------------
;;-    fxbhmake, hdr, Npix,       'changeme', /init, /date, /extver
;;-    sxdelpar, hdr, 'EXTNAME'
;;-    
;;-    sxaddpar, hdr, 'TTYPE1' , 'UP78',   ' Common mask UP78'
;;-    sxaddpar, hdr, 'TTYPE2' , 'UPA77',  ' Common mask UPA77'
;;-    sxaddpar, hdr, 'TTYPE3' , 'UPB77',  ' Common mask UPB77'
;;-    sxaddpar, hdr, 'TTYPE4' , 'HMIS',   ' Missing pixels, half-mission'
;;-    sxaddpar, hdr, 'TTYPE5' , 'YEAR',   ' Missing pixels, year'
;;-    
;;-    sxaddpar, hdr, 'COMMENT',  '  '  , after='TFORM5'       
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, hdr, 'TUNIT1',   '  '
;;-    sxaddpar, hdr, 'TUNIT2',   '  '
;;-    sxaddpar, hdr, 'TUNIT3',   '  '
;;-    sxaddpar, hdr, 'TUNIT4',   '  '
;;-    sxaddpar, hdr, 'TUNIT5',   '  '
;;-    
;;-    sxaddpar, hdr, 'COMMENT',  '  ', after='TUNIT5'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, hdr, 'EXTNAME',  'MASK-POL',    ' Extension name'
;;-    sxaddpar, hdr, 'AST-COMP', 'CMB',         ' Component'
;;-    sxaddpar, hdr, 'PIXTYPE',  'HEALPIX'
;;-    sxaddpar, hdr, 'POLCCONV', 'COSMO',       ' Polarization convention'
;;-    sxaddpar, hdr, 'COORDSYS', 'GALACTIC',    ' Coordinate system'
;;-    sxaddpar, hdr, 'ORDERING', 'NESTED',      ' Healpix ordering'
;;-    sxaddpar, hdr, 'NSIDE',     Nside,        ' Healpix Nside'
;;-    sxaddpar, hdr, 'FIRSTPIX',  0
;;-    sxaddpar, hdr, 'LASTPIX',   Npix-1
;;-    sxaddpar, hdr, 'FILENAME',  outfile,      ' FITS filename'
;;-    sxaddpar, hdr, 'BAD_DATA',  -1.63750E+30, ' HEALPIX bad pixel value'
;;-    sxaddpar, hdr, 'PROCVER',  procver,       ' Product version'
;;-   
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='PROCVER'
;;-    sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT',  'CMB Common and Missing pixel masks for Polarisation '
;;-    sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT',  'For further details see Planck Explanatory Supplement at:'
;;-    sxaddpar, hdr, 'COMMENT',  '  http://wiki.cosmos.esa.int/planckpla2015'
;;-    sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-    
;;-    print, ' >> Writing Polar Common Masks HDU:           '+outdir+outfile
;;-    mwrfits, msk, outdir+outfile, hdr, /sil ; 
;;-    if verb then print, hdr
;;-    
;;-    if defined(stop) then stop
;;- 
;;-    return
;;- end
;;- 
;;- ; ... and the "preferred" masks used in paper
;;- 
;;- 
;;- ;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;;- ; export sevem fg-subtracted frequency maps
;;- ;-----------------------------------------------------------------------------
;;- 
;;- pro sevem, freq, relnum, procver, outdir=outdir, test=test
;;- 
;;- ;   on_error, 2
;;- 
;;-    if not defined(freq) then return
;;-    if not defined(procver) then procver = 'NotGiven'
;;-    if not defined(outdir)  then outdir  = './'
;;-    if keyword_set(test)    then test = 1 else test = 0
;;-    freq = strtrim(freq,2)
;;- 
;;-    types  = ["", "_hm1","_hm2", "_yr1","_yr2", "_hr1","_hr2"]
;;-    covers = ["full", "halfmission-1","halfmission-2", "year-1","year-2", "ringhalf-1","ringhalf-2"] 
;;- 
;;-    ; Temperature maps
;;-    if (freq eq '100' or freq eq '143' or freq eq '217' ) then begin
;;-       indir = '/redtruck/ashdown/repository/comparison/dx11_v2/freq_cleaned/'
;;-       indir = '/redtruck/ashdown/repository/comparison/dx11_v2/int/freq_cleaned/'
;;-       Nside = 2048LL  &   Npix  = 12L*Nside*Nside
;;-       for i=0,6 do begin
;;-          cover = covers[i]
;;-          outfile = 'COM_CMB_IQU-'+freq+'-fgsub-sevem-field-Int_2048_'+relnum+'_'+cover+'.fits'
;;-    
;;-          ; build structure, read file, and
;;-          print, " >> Read file, build INT structure, and fill it"
;;-          spawn, 'ls -d1 '+indir+'dx11_v2_sevem_int_'+freq+'_cmb'+types[i]+'_orig_2048.fits', tfile
;;-          print, " >> source file is ", tfile
;;-          read_fits_map, tfile,  sig , 1, hdsig, ordering=order ; 
;;-          int = replicate({I_stokes: 0.0}, Npix) 
;;-          if strmid(order,0,1) eq "R" then int.i_stokes = reorder(float(sig ), /r2n) else int.I_stokes = float(sig)
;;-       
;;-          ;-----------------------------------------------------------------------------
;;-          print, " >> Build header for Temp maps"
;;-          ;-----------------------------------------------------------------------------
;;-          fxbhmake, hdr, Npix,       'changeme', /init, /date, /extver
;;-          sxdelpar, hdr, 'EXTNAME'
;;-          
;;-          sxaddpar, hdr, 'TTYPE1' , 'I_stokes     ',  ' Stokes I'
;;-       
;;-          sxaddpar, hdr, 'COMMENT',  '  '  , after='TFORM1'       
;;-          sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
;;-          sxaddpar, hdr, 'COMMENT',  '  '
;;-       
;;-          sxaddpar, hdr, 'TUNIT1',   'K_cmb'
;;-       
;;-          sxaddpar, hdr, 'COMMENT',  '  '  , after='TUNIT1'
;;-          sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-          sxaddpar, hdr, 'COMMENT',  '  '
;;-          
;;-          sxaddpar, hdr, 'EXTNAME',  'CMB',           ' Extension name'
;;-          sxaddpar, hdr, 'FREQ',      fix(freq),    ' Frequency'
;;-          sxaddpar, hdr, 'PIXTYPE',  'HEALPIX'
;;-          sxaddpar, hdr, 'POLCCONV', 'COSMO',       ' Polarization convention'
;;-          sxaddpar, hdr, 'COORDSYS', 'GALACTIC',    ' Coordinate system'
;;-          sxaddpar, hdr, 'ORDERING', 'NESTED',      ' Healpix ordering'
;;-          sxaddpar, hdr, 'NSIDE',     Nside,        ' Healpix Nside'
;;-          sxaddpar, hdr, 'FIRSTPIX',  0
;;-          sxaddpar, hdr, 'LASTPIX',   Npix-1
;;-          sxaddpar, hdr, 'FILENAME',  outfile,      ' FITS filename'
;;-          sxaddpar, hdr, 'BAD_DATA',  -1.63750E+30, ' HEALPIX bad pixel value'
;;-          sxaddpar, hdr, 'METHOD',   'sevem',       ' Cleaning method'
;;-          sxaddpar, hdr, 'PROCVER',  procver,       ' Product version'
;;-          
;;-          sxaddpar, hdr, 'COMMENT',  '  ' , after='PROCVER'
;;-          sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-          sxaddpar, hdr, 'COMMENT',  'CMB products from foregrounds-subtracted sky map from sevem'
;;-          sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-          sxaddpar, hdr, 'COMMENT',  'For further details see Planck Explanatory Supplement at:'
;;-          sxaddpar, hdr, 'COMMENT',  '  http://wiki.cosmos.esa.int/planckpla2015'
;;-          sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-          
;;-          print, ' >> Writing sevem clean Int map HDU:    '+outdir+outfile
;;-          mwrfits, int, outdir+outfile, hdr, /create ; sevem fgsub map 
;;-          if test then print, hdr
;;-      endfor 
;;- ;     stop
;;-  end  
;;- 
;;-    ; Polar maps
;;-    if (freq eq '100' or freq eq '143' or freq eq '070' ) then begin
;;-       Nside = 1024LL  &   Npix  = 12L*Nside*Nside
;;- 
;;-       if freq eq '070' then begin
;;-           indir = '/redtruck/ashdown/repository/comparison/dx11_v2/pol_qu/freq_cleaned/'
;;-           highpass = '_orig'
;;-       endif else begin
;;-           indir = '/redtruck/ashdown/repository/comparison/dx11_v2/pol_qu_highpass_20_40/freq_cleaned/'
;;-           highpass = '_hp_20_40_010a'
;;-       end 
;;-       for i=0,6 do begin
;;-           cover = covers[i]
;;-           outfile = 'COM_CMB_IQU-'+freq+'-fgsub-sevem-field-Pol_1024_'+relnum+'_'+cover+'.fits'
;;- 
;;-           print, " >> Read file, build QU structure, and fill it"
;;-           spawn, 'ls -d1 '+indir+'dx11_v2_sevem_pol_case1_'+freq+'_cmb'+types[0]+highpass+'_1024.fits', pfile
;;-           print, " >> source file is ", pfile
;;-           read_fits_map, pfile,  sig , 1, hdsig, ordering=order ; 
;;-           qu = replicate({Q_stokes: 0.0, U_stokes: 0.0}, Npix) 
;;-           if strmid(order,0,1) eq "R" then qu.q_stokes = reorder(float(sig[*,0] ), /r2n) else qu.q_stokes = float(sig[*,0])
;;-           if strmid(order,0,1) eq "R" then qu.u_stokes = reorder(float(sig[*,1] ), /r2n) else qu.u_stokes = float(sig[*,1])
;;-       
;;-           ;-----------------------------------------------------------------------------
;;-           print, " >> Build header for QU maps and export"
;;-           ;-----------------------------------------------------------------------------
;;-           fxbhmake, qdr, Npix,       'changeme', /init, /date, /extver
;;-           sxdelpar, qdr, 'EXTNAME'
;;-       
;;-           sxaddpar, qdr, 'TTYPE1' , 'Q_stokes',         ' Stokes Q'
;;-           sxaddpar, qdr, 'TTYPE2' , 'U_stokes',         ' Stokes U'
;;-       
;;-           sxaddpar, qdr, 'COMMENT',  '  '  , after='TFORM2'       
;;-           sxaddpar, qdr, 'COMMENT',  ' *** Column units *** '
;;-           sxaddpar, qdr, 'COMMENT',  '  '
;;-           
;;-           sxaddpar, qdr, 'TUNIT1',   'K_cmb'
;;-           sxaddpar, qdr, 'TUNIT2',   'K_cmb'
;;-       
;;-           sxaddpar, qdr, 'COMMENT',  '  '  , after='TUNIT2'
;;-           sxaddpar, qdr, 'COMMENT',  ' *** Planck params *** '
;;-           sxaddpar, qdr, 'COMMENT',  '  '
;;-           
;;-           sxaddpar, qdr, 'EXTNAME',  'CMB',         ' Extension name'
;;-           sxaddpar, qdr, 'FREQ',      fix(freq),    ' Frequency'
;;-           sxaddpar, qdr, 'PIXTYPE',  'HEALPIX'
;;-           sxaddpar, qdr, 'POLCCONV', 'COSMO',       ' Polarization convention'
;;-           sxaddpar, qdr, 'COORDSYS', 'GALACTIC',    ' Coordinate system'
;;-           sxaddpar, qdr, 'ORDERING', 'NESTED',      ' Healpix ordering'
;;-           sxaddpar, qdr, 'NSIDE',     Nside,        ' Healpix Nside'
;;-           sxaddpar, qdr, 'FIRSTPIX',  0
;;-           sxaddpar, qdr, 'LASTPIX',   Npix-1
;;-           sxaddpar, qdr, 'FILENAME',  outfile,      ' FITS filename'
;;-           sxaddpar, qdr, 'BAD_DATA',  -1.63750E+30, ' HEALPIX bad pixel value'
;;-           sxaddpar, qdr, 'METHOD',   'sevem',       ' Cleaning method'
;;-           sxaddpar, qdr, 'PROCVER',  procver,       ' Product version'
;;-       
;;-           sxaddpar, qdr, 'COMMENT',  '  ' , after='PROCVER'
;;-           sxaddpar, qdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-           sxaddpar, qdr, 'COMMENT',  'CMB products from foregrounds-subtracted sky map from sevem'
;;-           sxaddpar, qdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-           sxaddpar, qdr, 'COMMENT',  'For further details see Planck Explanatory Supplement at:'
;;-           sxaddpar, qdr, 'COMMENT',  '  http://wiki.cosmos.esa.int/planckpla2015'
;;-           sxaddpar, qdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-           
;;-           print, ' >> Writing sevem clean QU maps HDU:    '+outdir+outfile
;;-           mwrfits, qu, outdir+outfile, qdr, /create ; sevem fgsub map 
;;-           if test then print, qdr
;;-       endfor 
;;- ;      stop
;;-    end 
;;- 
;;-    return
;;- end
;;- 
;;- 
;;- ;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;;- ; export CMB power spectrum - use tables produced by Silvia Galli
;;- ;-----------------------------------------------------------------------------
;;- 
;;- pro cmbspec, relnum, procver, outdir=outdir, verb=verb
;;- 
;;-    on_error, 2
;;- 
;;-    if not defined(procver) then procver = 'NotGiven'
;;-    if not defined(outdir)  then outdir  = './'
;;-    if keyword_set(verb)    then verb = 1 else verb = 0
;;-    outfile = 'COM_PowerSpect_CMB_'+relnum+'.fits'
;;- 
;;-    ; the data files
;;-    filesdir = "/mnt/gpfs3/opsman/export/Releases/HFI_PR2_beta/cmbspec/"
;;- 
;;-    ;-----------------------------------------------------------------------------
;;-    ; low-ell TT spectrum (unbinned, l=2-29)
;;-    ;-----------------------------------------------------------------------------
;;-    filecom = filesdir+"cls_commander_2014_rc_TT.txt"
;;-    readcol, filecom, form='i,x,x,f,f,f', l, ml, eupper, elower, count=nl,/sil, numline=29
;;-    spect = replicate({ell:0, d_ell:0.0, errup:0.0, errdown:0.0}, nl)
;;-    spect.ell   = l
;;-    spect.D_ell = ml
;;-    spect.errup = eupper
;;-    spect.errdown = elower
;;- 
;;-    fxbhmake, hdr, nl,       'Changeme',  ' Extension name', /init, /date, /extver
;;-    sxdelpar, hdr, 'EXTNAME'
;;- 
;;-    sxaddpar, hdr, 'TTYPE1',   'ell',      ' ell of multipole'
;;-    sxaddpar, hdr, 'TTYPE2',   'D_ell',    ' D(ell) = ell(ell+1)C_l / 2pi' 
;;-    sxaddpar, hdr, 'TTYPE3',   'errup',    ' upper error on D(ell)'
;;-    sxaddpar, hdr, 'TTYPE4',   'errdown',  ' lower error onD(ell)'
;;-    
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TFORM4'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    sxaddpar, hdr, 'TUNIT1',   'none'
;;-    sxaddpar, hdr, 'TUNIT2',   'muKcmb^2',  'uK_CMB^2'   
;;-    sxaddpar, hdr, 'TUNIT3',   'muKcmb^2',  'uK_CMB^2'   
;;-    sxaddpar, hdr, 'TUNIT4',   'muKcmb^2',  'uK_CMB^2'
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT4'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, hdr, 'EXTNAME',  'TTLOLUNB',   ' low-ell CMB TT power spectrum, unbinned'
;;-    sxaddpar, hdr, 'L_MIN',     min(l),      ' lowest ell'
;;-    sxaddpar, hdr, 'L_MAX',     max(l),      ' highest ell'
;;-    sxaddpar, hdr, 'FILENAME',  outfile,     ' FITS filename'
;;-    sxaddpar, hdr, 'PROCVER',  procver,      ' Product version'
;;- 
;;-    sxaddpar, hdr, 'COMMENT', '  ' , after='PROCVER'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT', ' Low-ell, CMB TT power sepctrum (unbinned)'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;- 
;;-    print, ' >> Exporting '+outdir+outfile
;;-    mwrfits, spect, outdir+outfile, hdr,/create         ; CMB low-ell power spectrum (commander / CAMSPEC)
;;-    if verb then print, hdr
;;- 
;;-    ;-----------------------------------------------------------------------------
;;-    ; low-ell TE spectrum (unbinned, l=2-29)
;;-    ;-----------------------------------------------------------------------------
;;-    filecom = filesdir+"bolpol_TE.txt"
;;-    readcol, filecom, form='i,x,x,f,f,f', l, ml, eupper, elower, count=nl,/sil
;;-    spect = replicate({ell:0, d_ell:0.0, errup:0.0, errdown:0.0}, nl)
;;-    spect.ell   = l
;;-    spect.D_ell = ml
;;-    spect.errup = eupper
;;-    spect.errdown = elower
;;- 
;;-    fxbhmake, hdr, nl,       'Changeme',  ' Extension name', /init, /date, /extver
;;-    sxdelpar, hdr, 'EXTNAME'
;;- 
;;-    sxaddpar, hdr, 'TTYPE1',   'ell',      ' ell of multipole'
;;-    sxaddpar, hdr, 'TTYPE2',   'D_ell',    ' D(ell) = ell(ell+1)C_l / 2pi' 
;;-    sxaddpar, hdr, 'TTYPE3',   'errup',    ' upper error on D(ell)'
;;-    sxaddpar, hdr, 'TTYPE4',   'errdown',  ' lower error onD(ell)'
;;-    
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TFORM4'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    sxaddpar, hdr, 'TUNIT1',   'none'
;;-    sxaddpar, hdr, 'TUNIT2',   'muKcmb^2',  'uK_CMB^2'   
;;-    sxaddpar, hdr, 'TUNIT3',   'muKcmb^2',  'uK_CMB^2'   
;;-    sxaddpar, hdr, 'TUNIT4',   'muKcmb^2',  'uK_CMB^2'
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT4'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, hdr, 'EXTNAME',  'TELOLUNB',   ' low-ell CMB TE power spectrum, unbinned'
;;-    sxaddpar, hdr, 'L_MIN',     min(l),      ' lowest ell'
;;-    sxaddpar, hdr, 'L_MAX',     max(l),      ' highest ell'
;;-    sxaddpar, hdr, 'FILENAME',  outfile,     ' FITS filename'
;;-    sxaddpar, hdr, 'PROCVER',  procver,      ' Product version'
;;- 
;;-    sxaddpar, hdr, 'COMMENT', '  ' , after='PROCVER'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT', ' Low-ell, CMB TE power sepctrum (unbinned)'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;- 
;;-    print, ' >> Exporting '+outdir+outfile
;;-    mwrfits, spect, outdir+outfile, hdr,/sil         ; CMB low-ell power spectrum (Natoli)
;;-    if verb then print, hdr
;;- 
;;-    ;-----------------------------------------------------------------------------
;;-    ; low-ell EE spectrum (unbinned, l=2-29)
;;-    ;-----------------------------------------------------------------------------
;;-    filecom = filesdir+"bolpol_EE.txt"
;;-    readcol, filecom, form='i,x,x,f,f,f', l, ml, eupper, elower, count=nl,/sil
;;-    spect = replicate({ell:0, d_ell:0.0, errup:0.0, errdown:0.0}, nl)
;;-    spect.ell   = l
;;-    spect.D_ell = ml
;;-    spect.errup = eupper
;;-    spect.errdown = elower
;;- 
;;-    fxbhmake, hdr, nl,       'Changeme',  ' Extension name', /init, /date, /extver
;;-    sxdelpar, hdr, 'EXTNAME'
;;- 
;;-    sxaddpar, hdr, 'TTYPE1',   'ell',      ' ell of multipole'
;;-    sxaddpar, hdr, 'TTYPE2',   'D_ell',    ' D(ell) = ell(ell+1)C_l / 2pi' 
;;-    sxaddpar, hdr, 'TTYPE3',   'errup',    ' upper error on D(ell)'
;;-    sxaddpar, hdr, 'TTYPE4',   'errdown',  ' lower error onD(ell)'
;;-    
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TFORM4'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    sxaddpar, hdr, 'TUNIT1',   'none'
;;-    sxaddpar, hdr, 'TUNIT2',   'muKcmb^2',  'uK_CMB^2'   
;;-    sxaddpar, hdr, 'TUNIT3',   'muKcmb^2',  'uK_CMB^2'   
;;-    sxaddpar, hdr, 'TUNIT4',   'muKcmb^2',  'uK_CMB^2'
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT4'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, hdr, 'EXTNAME',  'EELOLUNB',   ' low-ell CMB EE power spectrum, unbinned'
;;-    sxaddpar, hdr, 'L_MIN',     min(l),      ' lowest ell'
;;-    sxaddpar, hdr, 'L_MAX',     max(l),      ' highest ell'
;;-    sxaddpar, hdr, 'FILENAME',  outfile,     ' FITS filename'
;;-    sxaddpar, hdr, 'PROCVER',  procver,      ' Product version'
;;- 
;;-    sxaddpar, hdr, 'COMMENT', '  ' , after='PROCVER'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT', ' Low-ell, CMB EE power sepctrum (unbinned)'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;- 
;;-    print, ' >> Exporting '+outdir+outfile
;;-    mwrfits, spect, outdir+outfile, hdr,/sil         ; CMB low-ell power spectrum (Natoli)
;;-    if verb then print, hdr
;;- 
;;-    ;-----------------------------------------------------------------------------
;;-    ; low-ell TB spectrum (unbinned, l=2-29)
;;-    ;-----------------------------------------------------------------------------
;;-    filecom = filesdir+"bolpol_TB.txt"
;;-    readcol, filecom, form='i,x,x,f,f,f', l, ml, eupper, elower, count=nl,/sil
;;-    spect = replicate({ell:0, d_ell:0.0, errup:0.0, errdown:0.0}, nl)
;;-    spect.ell   = l
;;-    spect.D_ell = ml
;;-    spect.errup = eupper
;;-    spect.errdown = elower
;;- 
;;-    fxbhmake, hdr, nl,       'Changeme',  ' Extension name', /init, /date, /extver
;;-    sxdelpar, hdr, 'EXTNAME'
;;- 
;;-    sxaddpar, hdr, 'TTYPE1',   'ell',      ' ell of multipole'
;;-    sxaddpar, hdr, 'TTYPE2',   'D_ell',    ' D(ell) = ell(ell+1)C_l / 2pi' 
;;-    sxaddpar, hdr, 'TTYPE3',   'errup',    ' upper error on D(ell)'
;;-    sxaddpar, hdr, 'TTYPE4',   'errdown',  ' lower error onD(ell)'
;;-    
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TFORM4'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    sxaddpar, hdr, 'TUNIT1',   'none'
;;-    sxaddpar, hdr, 'TUNIT2',   'muKcmb^2',  'uK_CMB^2'   
;;-    sxaddpar, hdr, 'TUNIT3',   'muKcmb^2',  'uK_CMB^2'   
;;-    sxaddpar, hdr, 'TUNIT4',   'muKcmb^2',  'uK_CMB^2'
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT4'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, hdr, 'EXTNAME',  'TBLOLUNB',   ' low-ell CMB TB power spectrum, unbinned'
;;-    sxaddpar, hdr, 'L_MIN',     min(l),      ' lowest ell'
;;-    sxaddpar, hdr, 'L_MAX',     max(l),      ' highest ell'
;;-    sxaddpar, hdr, 'FILENAME',  outfile,     ' FITS filename'
;;-    sxaddpar, hdr, 'PROCVER',  procver,      ' Product version'
;;- 
;;-    sxaddpar, hdr, 'COMMENT', '  ' , after='PROCVER'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT', ' Low-ell, CMB TB power sepctrum (unbinned)'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;- 
;;-    print, ' >> Exporting '+outdir+outfile
;;-    mwrfits, spect, outdir+outfile, hdr,/sil         ; CMB low-ell power spectrum (Natoli)
;;-    if verb then print, hdr
;;- 
;;-    ;-----------------------------------------------------------------------------
;;-    ; low-ell EB spectrum (unbinned, l=2-29)
;;-    ;-----------------------------------------------------------------------------
;;-    filecom = filesdir+"bolpol_EB.txt"
;;-    readcol, filecom, form='i,x,x,f,f,f', l, ml, eupper, elower, count=nl,/sil
;;-    spect = replicate({ell:0, d_ell:0.0, errup:0.0, errdown:0.0}, nl)
;;-    spect.ell   = l
;;-    spect.D_ell = ml
;;-    spect.errup = eupper
;;-    spect.errdown = elower
;;- 
;;-    fxbhmake, hdr, nl,       'Changeme',  ' Extension name', /init, /date, /extver
;;-    sxdelpar, hdr, 'EXTNAME'
;;- 
;;-    sxaddpar, hdr, 'TTYPE1',   'ell',      ' ell of multipole'
;;-    sxaddpar, hdr, 'TTYPE2',   'D_ell',    ' D(ell) = ell(ell+1)C_l / 2pi' 
;;-    sxaddpar, hdr, 'TTYPE3',   'errup',    ' upper error on D(ell)'
;;-    sxaddpar, hdr, 'TTYPE4',   'errdown',  ' lower error onD(ell)'
;;-    
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TFORM4'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    sxaddpar, hdr, 'TUNIT1',   'none'
;;-    sxaddpar, hdr, 'TUNIT2',   'muKcmb^2',  'uK_CMB^2'   
;;-    sxaddpar, hdr, 'TUNIT3',   'muKcmb^2',  'uK_CMB^2'   
;;-    sxaddpar, hdr, 'TUNIT4',   'muKcmb^2',  'uK_CMB^2'
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT4'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, hdr, 'EXTNAME',  'EBLOLUNB',   ' low-ell CMB EB power spectrum, unbinned'
;;-    sxaddpar, hdr, 'L_MIN',     min(l),      ' lowest ell'
;;-    sxaddpar, hdr, 'L_MAX',     max(l),      ' highest ell'
;;-    sxaddpar, hdr, 'FILENAME',  outfile,     ' FITS filename'
;;-    sxaddpar, hdr, 'PROCVER',  procver,      ' Product version'
;;- 
;;-    sxaddpar, hdr, 'COMMENT', '  ' , after='PROCVER'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT', ' Low-ell, CMB EB power sepctrum (unbinned)'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;- 
;;-    print, ' >> Exporting '+outdir+outfile
;;-    mwrfits, spect, outdir+outfile, hdr,/sil         ; CMB low-ell power spectrum (Natoli)
;;-    if verb then print, hdr
;;- 
;;-    ;-----------------------------------------------------------------------------
;;-    ; low-ell BB spectrum (unbinned, l=2-29)
;;-    ;-----------------------------------------------------------------------------
;;-    filecom = filesdir+"bolpol_BB.txt"
;;-    readcol, filecom, form='i,x,x,f,f,f', l, ml, eupper, elower, count=nl,/sil
;;-    spect = replicate({ell:0, d_ell:0.0, errup:0.0, errdown:0.0}, nl)
;;-    spect.ell   = l
;;-    spect.D_ell = ml
;;-    spect.errup = eupper
;;-    spect.errdown = elower
;;- 
;;-    fxbhmake, hdr, nl,       'Changeme',  ' Extension name', /init, /date, /extver
;;-    sxdelpar, hdr, 'EXTNAME'
;;- 
;;-    sxaddpar, hdr, 'TTYPE1',   'ell',      ' ell of multipole'
;;-    sxaddpar, hdr, 'TTYPE2',   'D_ell',    ' D(ell) = ell(ell+1)C_l / 2pi' 
;;-    sxaddpar, hdr, 'TTYPE3',   'errup',    ' upper error on D(ell)'
;;-    sxaddpar, hdr, 'TTYPE4',   'errdown',  ' lower error onD(ell)'
;;-    
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TFORM4'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    sxaddpar, hdr, 'TUNIT1',   'none'
;;-    sxaddpar, hdr, 'TUNIT2',   'muKcmb^2',  'uK_CMB^2'   
;;-    sxaddpar, hdr, 'TUNIT3',   'muKcmb^2',  'uK_CMB^2'   
;;-    sxaddpar, hdr, 'TUNIT4',   'muKcmb^2',  'uK_CMB^2'
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT4'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, hdr, 'EXTNAME',  'BBLOLUNB',   ' low-ell CMB BB power spectrum, unbinned'
;;-    sxaddpar, hdr, 'L_MIN',     min(l),      ' lowest ell'
;;-    sxaddpar, hdr, 'L_MAX',     max(l),      ' highest ell'
;;-    sxaddpar, hdr, 'FILENAME',  outfile,     ' FITS filename'
;;-    sxaddpar, hdr, 'PROCVER',  procver,      ' Product version'
;;- 
;;-    sxaddpar, hdr, 'COMMENT', '  ' , after='PROCVER'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT', ' Low-ell, CMB BB power sepctrum (unbinned)'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;- 
;;-    print, ' >> Exporting '+outdir+outfile
;;-    mwrfits, spect, outdir+outfile, hdr,/sil         ; CMB low-ell power spectrum (Natoli)
;;-    if verb then print, hdr
;;- 
;;-    ;-----------------------------------------------------------------------------
;;-    ; high ell TT CMB spectrum (binned)
;;-    ;-----------------------------------------------------------------------------
;;-    filespt = filesdir+"PowerSpect_DL_TT_bin30_R2.txt"
;;-    readcol, filespt, form='f,i,i,f,f', ll,lmin,lmax,dd,err,/sil, count=nl
;;- 
;;-    spect = replicate({ell:0.0, lmin:0, lmax:0, D_ell:0.0, err:0.0}, nl )
;;-    spect.ell   = ll
;;-    spect.lmin  = lmin
;;-    spect.lmax  = lmax
;;-    spect.D_ell = dd
;;-    spect.err   = err
;;- 
;;-    fxbhmake, hdr, nl,       'Changeme',  ' Extension name', /init, /date, /extver
;;-    sxdelpar, hdr, 'EXTNAME'
;;-    sxaddpar, hdr, 'TTYPE1',   'ell',     ' mean ell of band' 
;;-    sxaddpar, hdr, 'TTYPE2',   'lmin',    ' lowest ell of band'
;;-    sxaddpar, hdr, 'TTYPE3',   'lmax',    ' highest ell of band'
;;-    sxaddpar, hdr, 'TTYPE4',   'D_ell',   ' D(ell) = ell(ell+1)C_l / 2pi'
;;-    sxaddpar, hdr, 'TTYPE5',   'err',     ' err on D(l)'
;;-    
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TFORM5'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    sxaddpar, hdr, 'TUNIT4',   'muKcmb^2',  'uK_CMB^2'   
;;-    sxaddpar, hdr, 'TUNIT5',   'muKcmb^2',  'uK_CMB^2'
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT5'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, hdr, 'EXTNAME',  'TTHILBIN',  ' high-ell CMB TT power spectrum (binned)'
;;-    sxaddpar, hdr, 'L_MIN',     min(ll),    ' lowest lmean'
;;-    sxaddpar, hdr, 'L_MAX',     max(ll),    ' highest lmean'
;;-    sxaddpar, hdr, 'FILENAME',  outfile,    ' FITS filename'
;;-    sxaddpar, hdr, 'PROCVER',  procver,     ' Product version'
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='PROCVER'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT', ' Binned CMB TT Power Sepctrum (binned)'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;- 
;;-    print, ' >> Exporting '+outdir+outfile
;;-    mwrfits, spect, outdir+outfile, hdr,/sil          ; - append (binned) high-ell TT CMB power spectrum
;;-    if verb then print, hdr
;;- 
;;-    ;-----------------------------------------------------------------------------
;;-    ; high ell TT CMB spectrum (unbinned)
;;-    ;-----------------------------------------------------------------------------
;;-    filespt = filesdir+"PowerSpect_DL_TT_unbinned_R2.txt"
;;-    readcol, filespt, form='i,f,f', ll,dd,err,/sil, count=nl
;;- 
;;-    spect = replicate({ell:0, D_ell:0.0, err:0.0}, nl )
;;-    spect.ell   = ll
;;-    spect.D_ell = dd
;;-    spect.err   = err
;;- 
;;-    fxbhmake, hdr, nl,       'Changeme',  ' Extension name', /init, /date, /extver
;;-    sxdelpar, hdr, 'EXTNAME'
;;-    sxaddpar, hdr, 'TTYPE1',   'ell',     ' ell of bin' 
;;-    sxaddpar, hdr, 'TTYPE2',   'D_ell',   ' D(ell) = ell(ell+1)C_l / 2pi'
;;-    sxaddpar, hdr, 'TTYPE3',   'err',     ' err on D(l)'
;;-    
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TFORM5'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    sxaddpar, hdr, 'TUNIT2',   'muKcmb^2',  'uK_CMB^2'   
;;-    sxaddpar, hdr, 'TUNIT3',   'muKcmb^2',  'uK_CMB^2'
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT5'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, hdr, 'EXTNAME',  'TTHILUNB',  ' high-ell CMB TT power spectrum (unbinned)'
;;-    sxaddpar, hdr, 'L_MIN',     min(ll),    ' lowest l'
;;-    sxaddpar, hdr, 'L_MAX',     max(ll),    ' highest l'
;;-    sxaddpar, hdr, 'FILENAME',  outfile,    ' FITS filename'
;;-    sxaddpar, hdr, 'PROCVER',  procver,     ' Product version'
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='PROCVER'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT', ' CMB TT Power Sepctrum (unbinned)'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;- 
;;-    print, ' >> Exporting '+outdir+outfile
;;-    mwrfits, spect, outdir+outfile, hdr,/sil          ; - append (unbinned) high-ell TT CMB power spectrum
;;-    if verb then print, hdr
;;- 
;;-    ;-----------------------------------------------------------------------------
;;-    ; high ell TE CMB spectrum (binned)
;;-    ;-----------------------------------------------------------------------------
;;-    filespt = filesdir+"PowerSpect_DL_TE_bin30_R2.txt"
;;-    readcol, filespt, form='f,i,i,f,f', ll,lmin,lmax,dd,err,/sil, count=nl
;;- 
;;-    spect = replicate({ell:0.0, lmin:0, lmax:0, D_ell:0.0, err:0.0}, nl )
;;-    spect.ell   = ll
;;-    spect.lmin  = lmin
;;-    spect.lmax  = lmax
;;-    spect.D_ell = dd
;;-    spect.err   = err
;;- 
;;-    fxbhmake, hdr, nl,       'Changeme',  ' Extension name', /init, /date, /extver
;;-    sxdelpar, hdr, 'EXTNAME'
;;-    sxaddpar, hdr, 'TTYPE1',   'ell',     ' mean ell of band' 
;;-    sxaddpar, hdr, 'TTYPE2',   'lmin',    ' lowest ell of band'
;;-    sxaddpar, hdr, 'TTYPE3',   'lmax',    ' highest ell of band'
;;-    sxaddpar, hdr, 'TTYPE4',   'D_ell',   ' D(ell) = ell(ell+1)C_l / 2pi'
;;-    sxaddpar, hdr, 'TTYPE5',   'err',     ' err on D(l)'
;;-    
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TFORM5'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    sxaddpar, hdr, 'TUNIT4',   'muKcmb^2',  'uK_CMB^2'   
;;-    sxaddpar, hdr, 'TUNIT5',   'muKcmb^2',  'uK_CMB^2'
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT5'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, hdr, 'EXTNAME',  'TEHILBIN',  ' high-ell CMB TT power spectrum (binned)'
;;-    sxaddpar, hdr, 'L_MIN',     min(ll),    ' lowest lmean'
;;-    sxaddpar, hdr, 'L_MAX',     max(ll),    ' highest lmean'
;;-    sxaddpar, hdr, 'FILENAME',  outfile,    ' FITS filename'
;;-    sxaddpar, hdr, 'PROCVER',  procver,     ' Product version'
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='PROCVER'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT', ' Binned CMB TE Power Sepctrum (binned)'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;- 
;;-    print, ' >> Exporting '+outdir+outfile
;;-    mwrfits, spect, outdir+outfile, hdr,/sil          ; - append (binned) high-ell TECMB power spectrum
;;-    if verb then print, hdr
;;- 
;;-    ;-----------------------------------------------------------------------------
;;-    ; high ell TE CMB spectrum (unbinned)
;;-    ;-----------------------------------------------------------------------------
;;-    filespt = filesdir+"PowerSpect_DL_TE_unbinned_R2.txt"
;;-    readcol, filespt, form='i,f,f', ll,dd,err,/sil, count=nl
;;- 
;;-    spect = replicate({ell:0, D_ell:0.0, err:0.0}, nl )
;;-    spect.ell   = ll
;;-    spect.D_ell = dd
;;-    spect.err   = err
;;- 
;;-    fxbhmake, hdr, nl,       'Changeme',  ' Extension name', /init, /date, /extver
;;-    sxdelpar, hdr, 'EXTNAME'
;;-    sxaddpar, hdr, 'TTYPE1',   'ell',     ' ell of bin' 
;;-    sxaddpar, hdr, 'TTYPE2',   'D_ell',   ' D(ell) = ell(ell+1)C_l / 2pi'
;;-    sxaddpar, hdr, 'TTYPE3',   'err',     ' err on D(l)'
;;-    
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TFORM5'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    sxaddpar, hdr, 'TUNIT2',   'muKcmb^2',  'uK_CMB^2'   
;;-    sxaddpar, hdr, 'TUNIT3',   'muKcmb^2',  'uK_CMB^2'
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT5'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, hdr, 'EXTNAME',  'TEHILUNB',  ' high-ell CMB TT power spectrum (unbinned)'
;;-    sxaddpar, hdr, 'L_MIN',     min(ll),    ' lowest l'
;;-    sxaddpar, hdr, 'L_MAX',     max(ll),    ' highest l'
;;-    sxaddpar, hdr, 'FILENAME',  outfile,    ' FITS filename'
;;-    sxaddpar, hdr, 'PROCVER',  procver,     ' Product version'
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='PROCVER'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT', ' CMB TE Power Sepctrum (unbinned)'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;- 
;;-    print, ' >> Exporting '+outdir+outfile
;;-    mwrfits, spect, outdir+outfile, hdr,/sil          ; - append (unbinned) high-ell TE CMB power spectrum
;;-    if verb then print, hdr
;;- 
;;-    ;-----------------------------------------------------------------------------
;;-    ; high ell EE CMB spectrum (binned)
;;-    ;-----------------------------------------------------------------------------
;;-    filespt = filesdir+"PowerSpect_DL_EE_bin30_R2.txt"
;;-    readcol, filespt, form='f,i,i,f,f', ll,lmin,lmax,dd,err,/sil, count=nl
;;- 
;;-    spect = replicate({ell:0.0, lmin:0, lmax:0, D_ell:0.0, err:0.0}, nl )
;;-    spect.ell   = ll
;;-    spect.lmin  = lmin
;;-    spect.lmax  = lmax
;;-    spect.D_ell = dd
;;-    spect.err   = err
;;- 
;;-    fxbhmake, hdr, nl,       'Changeme',  ' Extension name', /init, /date, /extver
;;-    sxdelpar, hdr, 'EXTNAME'
;;-    sxaddpar, hdr, 'TTYPE1',   'ell',     ' mean ell of band' 
;;-    sxaddpar, hdr, 'TTYPE2',   'lmin',    ' lowest ell of band'
;;-    sxaddpar, hdr, 'TTYPE3',   'lmax',    ' highest ell of band'
;;-    sxaddpar, hdr, 'TTYPE4',   'D_ell',   ' D(ell) = ell(ell+1)C_l / 2pi'
;;-    sxaddpar, hdr, 'TTYPE5',   'err',     ' err on D(l)'
;;-    
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TFORM5'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    sxaddpar, hdr, 'TUNIT4',   'muKcmb^2',  'uK_CMB^2'   
;;-    sxaddpar, hdr, 'TUNIT5',   'muKcmb^2',  'uK_CMB^2'
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT5'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, hdr, 'EXTNAME',  'EEHILBIN',  ' high-ell CMB EE power spectrum (binned)'
;;-    sxaddpar, hdr, 'L_MIN',     min(ll),    ' lowest lmean'
;;-    sxaddpar, hdr, 'L_MAX',     max(ll),    ' highest lmean'
;;-    sxaddpar, hdr, 'FILENAME',  outfile,    ' FITS filename'
;;-    sxaddpar, hdr, 'PROCVER',  procver,     ' Product version'
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='PROCVER'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT', ' Binned CMB EE Power Sepctrum (binned)'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;- 
;;-    print, ' >> Exporting '+outdir+outfile
;;-    mwrfits, spect, outdir+outfile, hdr,/sil          ; - append (binned) high-ell EE CMB power spectrum
;;-    if verb then print, hdr
;;- 
;;-    ;-----------------------------------------------------------------------------
;;-    ; high ell EE CMB spectrum (binned)
;;-    ;-----------------------------------------------------------------------------
;;-    filespt = filesdir+"PowerSpect_DL_EE_unbinned_R2.txt"
;;-    readcol, filespt, form='i,f,f', ll,dd,err,/sil, count=nl
;;- 
;;-    spect = replicate({ell:0, D_ell:0.0, err:0.0}, nl )
;;-    spect.ell   = ll
;;-    spect.D_ell = dd
;;-    spect.err   = err
;;- 
;;-    fxbhmake, hdr, nl,       'Changeme',  ' Extension name', /init, /date, /extver
;;-    sxdelpar, hdr, 'EXTNAME'
;;-    sxaddpar, hdr, 'TTYPE1',   'ell',     ' ell of bin' 
;;-    sxaddpar, hdr, 'TTYPE2',   'D_ell',   ' D(ell) = ell(ell+1)C_l / 2pi'
;;-    sxaddpar, hdr, 'TTYPE3',   'err',     ' err on D(l)'
;;-    
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TFORM5'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    sxaddpar, hdr, 'TUNIT2',   'muKcmb^2',  'uK_CMB^2'   
;;-    sxaddpar, hdr, 'TUNIT3',   'muKcmb^2',  'uK_CMB^2'
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT5'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, hdr, 'EXTNAME',  'EEHILUNB',  ' high-ell CMB EE power spectrum (unbinned)'
;;-    sxaddpar, hdr, 'L_MIN',     min(ll),    ' lowest l'
;;-    sxaddpar, hdr, 'L_MAX',     max(ll),    ' highest l'
;;-    sxaddpar, hdr, 'FILENAME',  outfile,    ' FITS filename'
;;-    sxaddpar, hdr, 'PROCVER',  procver,     ' Product version'
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='PROCVER'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT', ' CMB EE Power Sepctrum (unbinned)'
;;-    sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;- 
;;-    print, ' >> Exporting '+outdir+outfile
;;-    mwrfits, spect, outdir+outfile, hdr,/sil          ; - append (unbinned) high-ell EE CMB power spectrum
;;-    if verb then print, hdr
;;- 
;;-    return
;;- end
;;- 
;;- ; plot data in FITS file produced above
;;- pro plotcmb, file
;;- 
;;-    pp = mrdfits(file, 1)
;;-    plot, pp.ell, pp.d_ell, xr=[0,2500], ps=3, /ylog, yr=[0.01,10000]
;;- 
;;-    pp = mrdfits(file, 2)
;;-    oplot, pp.ell, pp.d_ell, col=68, ps=2
;;-    pp = mrdfits(file, 3)
;;-    oplot, pp.ell, pp.d_ell, col=222, ps=3
;;-    
;;-    pp = mrdfits(file, 4)
;;-    oplot, pp.ell, pp.d_ell, col=68, ps=1
;;- 
;;-    pp = mrdfits(file, 6)
;;-    oplot, pp.ell, pp.d_ell, col=99, ps=2
;;- stop
;;-    return
;;- end
;;- 
;;- ;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;;- ; export dust model maps to go with Boulanger et al paper
;;- ;-----------------------------------------------------------------------------
;;- 
;;- pro dustmodel, relnum, procver, outdir=outdir, $
;;-             stop=stop, verb=verb
;;- 
;;-    on_error, 2
;;- 
;;-    if not defined(relnum)  then begin
;;-        print, 'SYNTAX: '
;;-        print, '  dustopa, relnum, procver, outdir=outdir, $'
;;-        print, '          /test, /stop, /verb'    
;;-        return
;;-    endif
;;- 
;;-    if n_elements(relnum)  eq 0 then relnum  = 'R0.00'
;;-    if n_elements(procver) eq 0 then procver = 'XXX'
;;-    if not keyword_set(outdir)  then outdir  = './'
;;-    if keyword_set(test)        then test  = 1 else test  = 0
;;-    if keyword_set(verb)        then verb = 1 else verb = 0
;;- 
;;-    indir = '/wrk/fboulang/DL07/' ; 27/01/16
;;-    root = 'DL07_2048_1_11_000_111000_PR2_5am_'
;;- 
;;-    Nside = 2048LL   &   Npix = 12*Nside*Nside
;;-    com  = 'begin=0;end='+strtrim(Npix-1,2)
;;- 
;;-    ; 1. Parameters
;;-    name = 'COM_CompMap_Dust-DL07-Parameters_'+strtrim(Nside,2)+'_'+relnum
;;-    outfile = name+'.fits'
;;- 
;;-    name_mdust   = 'Parameter_M_dust.fits'
;;-    name_udust   = 'Parameter_M_dust_unc.fits'
;;-    name_qpah    = 'Parameter_q_PAH.fits'
;;-    name_upah    = 'Parameter_q_PAH_unc.fits'
;;-    name_fpdr    = 'Parameter_f_PDR.fits'
;;-    name_updr    = 'Parameter_f_PDR_unc.fits'
;;-    name_umin    = 'Parameter_U_min.fits'
;;-    name_uuin    = 'Parameter_U_min_unc.fits'
;;-    name_chi2    = 'Parameter_chi_sq_DOF.fits'
;;- 
;;-    map = replicate({sigma_mdust: 0.0,  sigma_mdust_unc: 0.0,   q_pah: 0.0,   q_pah_unc:0.0, $
;;-                           f_pdr: 0.0,   f_pdr_unc: 0.0,        u_min: 0.0,   u_min_unc: 0.0, $
;;-                        chi2_dof: 0.0}, Npix)
;;-    
;;-    read_fits_map, indir+"Parameters/"+root+name_mdust, data  &  map.sigma_mdust     = reorder(data, /r2n)
;;-    read_fits_map, indir+"Parameters/"+root+name_udust, data  &  map.sigma_mdust_unc = reorder(data, /r2n)
;;-    read_fits_map, indir+"Parameters/"+root+name_qpah,  data  &  map.q_pah     = reorder(data, /r2n)
;;-    read_fits_map, indir+"Parameters/"+root+name_upah,  data  &  map.q_pah_unc = reorder(data, /r2n)
;;-    read_fits_map, indir+"Parameters/"+root+name_fpdr,  data  &  map.f_pdr     = reorder(data, /r2n)
;;-    read_fits_map, indir+"Parameters/"+root+name_updr,  data  &  map.f_pdr_unc = reorder(data, /r2n)
;;-    read_fits_map, indir+"Parameters/"+root+name_umin,  data  &  map.u_min     = reorder(data, /r2n)
;;-    read_fits_map, indir+"Parameters/"+root+name_uuin,  data  &  map.u_min_unc = reorder(data, /r2n)
;;-    read_fits_map, indir+"Parameters/"+root+name_chi2,  data  &  map.chi2_dof  = reorder(data, /r2n)
;;- 
;;-    ; write header
;;-    fxbhmake, hdr, Npix,       'CHANGEME', ' Extension name', /init, /date, /extver
;;-    sxdelpar, hdr, 'EXTNAME'
;;-    sxaddpar, hdr, 'TTYPE1',   'sigma_mdust',     ' Dust mass surface density'
;;-    sxaddpar, hdr, 'TTYPE2',   'sigma_mdust_unc', ' Uncertainty of above'
;;-    sxaddpar, hdr, 'TTYPE3',   'q_PAH',       ' Dust mass fraction in small PAH grains'
;;-    sxaddpar, hdr, 'TTYPE4',   'q_PAH_unc',   ' Uncertainty of above'
;;-    sxaddpar, hdr, 'TTYPE5',   'f_pdr',       ' Fraction of total lumin. of hot dust'
;;-    sxaddpar, hdr, 'TTYPE6',   'f_pdr_unc',   ' Uncertainty of above'
;;-    sxaddpar, hdr, 'TTYPE7',   'u_min',       ' Intensity of starlight heating bulk of dust'
;;-    sxaddpar, hdr, 'TTYPE8',   'u_min_unc',   ' Uncertainty of above'
;;-    sxaddpar, hdr, 'TTYPE9',   'chi2_fof',    ' chi2 per deg. of freedom'
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TFORM9'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Column units *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;- 
;;-    sxaddpar, hdr, 'TUNIT1',   'Msun/kpc^2'
;;-    sxaddpar, hdr, 'TUNIT2',   'Msun/kpc^2 '
;;-    sxaddpar, hdr, 'TUNIT3',   ' ' 
;;-    sxaddpar, hdr, 'TUNIT4',   ' ' 
;;-    sxaddpar, hdr, 'TUNIT5',   ' ' 
;;-    sxaddpar, hdr, 'TUNIT6',   ' ' 
;;-    sxaddpar, hdr, 'TUNIT7',   ' ' 
;;-    sxaddpar, hdr, 'TUNIT8',   ' ' 
;;-    sxaddpar, hdr, 'TUNIT9',   ' ' 
;;-  
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='TUNIT9'
;;-    sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, hdr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, hdr, 'EXTNAME',  'COMP-MAP',    ' Extension name'
;;-    sxaddpar, hdr, 'AST-COMP', 'MODELDUST',   ' Component name'
;;-    sxaddpar, hdr, 'PIXTYPE',  'HEALPIX'
;;-    sxaddpar, hdr, 'COORDSYS', 'GALACTIC',    ' Coordinate system'
;;-    sxaddpar, hdr, 'ORDERING', 'NESTED',      ' Healpix ordering'
;;-    sxaddpar, hdr, 'NSIDE',     Nside,        ' Healpix Nside'
;;-    sxaddpar, hdr, 'FIRSTPIX',  0
;;-    sxaddpar, hdr, 'LASTPIX',   Npix-1
;;-    sxaddpar, hdr, 'BAD_DATA',  -1.63750E+30, ' HEALPIX bad pixel value'
;;-    sxaddpar, hdr, 'FILENAME',  outfile,      ' FITS filename'
;;-    sxaddpar, hdr, 'PROCVER',  procver,       ' Product version'
;;- 
;;-    sxaddpar, hdr, 'COMMENT',  '  ' , after='PROCVER'
;;-    sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-    sxaddpar, hdr, 'COMMENT',  'For further details see Planck Explanatory Supplement at:'
;;-    sxaddpar, hdr, 'COMMENT',  '  http://wiki.cosmos.esa.int/planckpla2015'
;;-    sxaddpar, hdr, 'COMMENT',  '------------------------------------------------------------------------'
;;- 
;;-    print, ' >> Exporting '+outdir+outfile
;;-    mwrfits, map, outdir+outfile, hdr, /create   ; Dust opacity map and model
;;-    if verb then print, hdr
;;- 
;;-    ; 2. Visible extinction maps
;;-    name = 'COM_CompMap_Dust-DL07-AvMaps_'+strtrim(Nside,2)+'_'+relnum
;;-    outfile = name+'.fits'
;;- 
;;-    name_av   = 'Parameter_Av_dust.fits'
;;-    name_uv   = 'Parameter_Av_dust_unc.fits'
;;-    name_rav  = 'Parameter_RAv_dust.fits'
;;-    name_urv  = 'Parameter_RAv_dust_unc.fits'
;;- 
;;-    map = replicate({AV_DL: 0.0,  AV_DL_unc: 0.0, AV_RQ: 0.0,  AV_RQ_unc: 0.0 }, Npix)
;;-    
;;-    read_fits_map, indir+"AvMaps/"+root+name_av,  data  &  map.AV_DL     = reorder(data, /r2n)
;;-    read_fits_map, indir+"AvMaps/"+root+name_uv,  data  &  map.AV_DL_unc = reorder(data, /r2n)
;;-    read_fits_map, indir+"AvMaps/"+root+name_rav, data  &  map.AV_RQ     = reorder(data, /r2n)
;;-    read_fits_map, indir+"AvMaps/"+root+name_urv, data  &  map.AV_RQ_unc = reorder(data, /r2n)
;;- 
;;-    ; write header
;;-    fxbhmake, ldr, Npix,       'CHANGEME', ' Extension name', /init, /date, /extver
;;-    sxdelpar, ldr, 'EXTNAME'
;;-    sxaddpar, ldr, 'TTYPE1',   'AV_DL',     ' Av from DLmodel'
;;-    sxaddpar, ldr, 'TTYPE2',   'AV_DL_unc', ' Uncertainty of above'
;;-    sxaddpar, ldr, 'TTYPE3',   'AV_RQ',     ' Av normalizes to SDSS QSOs'
;;-    sxaddpar, ldr, 'TTYPE4',   'AV_RQ_unc', ' Uncertainty of above'
;;- 
;;-    sxaddpar, ldr, 'COMMENT',  '  ' , after='TFORM4'
;;-    sxaddpar, ldr, 'COMMENT',  ' *** Column units *** '
;;-    sxaddpar, ldr, 'COMMENT',  '  '
;;- 
;;-    sxaddpar, ldr, 'TUNIT1',   'mag'
;;-    sxaddpar, ldr, 'TUNIT2',   'mag'
;;-    sxaddpar, ldr, 'TUNIT3',   'mag' 
;;-    sxaddpar, ldr, 'TUNIT4',   'mag' 
;;-  
;;-    sxaddpar, ldr, 'COMMENT',  '  ' , after='TUNIT4'
;;-    sxaddpar, ldr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, ldr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, ldr, 'EXTNAME',  'COMP-MAP',    ' Extension name'
;;-    sxaddpar, ldr, 'AST-COMP', 'MODELDUST',   ' Component name'
;;-    sxaddpar, ldr, 'PIXTYPE',  'HEALPIX'
;;-    sxaddpar, ldr, 'COORDSYS', 'GALACTIC',    ' Coordinate system'
;;-    sxaddpar, ldr, 'ORDERING', 'NESTED',      ' Healpix ordering'
;;-    sxaddpar, ldr, 'NSIDE',     Nside,        ' Healpix Nside'
;;-    sxaddpar, ldr, 'FIRSTPIX',  0
;;-    sxaddpar, ldr, 'LASTPIX',   Npix-1
;;-    sxaddpar, ldr, 'BAD_DATA',  -1.63750E+30, ' HEALPIX bad pixel value'
;;-    sxaddpar, ldr, 'FILENAME',  outfile,      ' FITS filename'
;;-    sxaddpar, ldr, 'PROCVER',  procver,       ' Product version'
;;- 
;;-    sxaddpar, ldr, 'COMMENT',  '  ' , after='PROCVER'
;;-    sxaddpar, ldr, 'COMMENT',  '------------------------------------------------------------------------'
;;-    sxaddpar, ldr, 'COMMENT',  'For further details see Planck Explanatory Supplement at:'
;;-    sxaddpar, ldr, 'COMMENT',  '  http://wiki.cosmos.esa.int/planckpla2015'
;;-    sxaddpar, ldr, 'COMMENT',  '------------------------------------------------------------------------'
;;- 
;;-    print, ' >> Exporting '+outdir+outfile
;;-    mwrfits, map, outdir+outfile, ldr, /create   ; Dust opacity map and model
;;-    if verb then print, ldr
;;- 
;;-    ; 3. Model fluxes
;;-    name = 'COM_CompMap_Dust-DL07-ModelFluxes_'+strtrim(Nside,2)+'_'+relnum
;;-    outfile = name+'.fits'
;;- 
;;-    name_1  = 'Model_PLANCK_857.fits'
;;-    name_2  = 'Model_PLANCK_545.fits'
;;-    name_3  = 'Model_PLANCK_353.fits'
;;-    name_4  = 'Model_WISE_12.fits'
;;-    name_5  = 'Model_IRAS_60.fits'
;;-    name_6  = 'Model_IRAS_100.fits'  
;;- 
;;-    map = replicate({Planck_857: 0.0,  Planck_545: 0.0, Planck_353: 0.0,  $
;;-                     Wise_12: 0.0, IRAS_60: 0.0, IRAS_100: 0.0  }, Npix)
;;-    
;;-    read_fits_map, indir+"ModelFluxes/"+root+name_1, data  &  map.PLANCK_857 = reorder(data, /r2n)
;;-    read_fits_map, indir+"ModelFluxes/"+root+name_2, data  &  map.PLANCK_545 = reorder(data, /r2n)
;;-    read_fits_map, indir+"ModelFluxes/"+root+name_3, data  &  map.PLANCK_353 = reorder(data, /r2n)
;;-    read_fits_map, indir+"ModelFluxes/"+root+name_4, data  &  map.WISE_12    = reorder(data, /r2n)
;;-    read_fits_map, indir+"ModelFluxes/"+root+name_5, data  &  map.IRAS_60    = reorder(data, /r2n)
;;-    read_fits_map, indir+"ModelFluxes/"+root+name_6, data  &  map.IRAS_100   = reorder(data, /r2n)
;;- 
;;-    ; write header
;;-    fxbhmake, kdr, Npix,       'CHANGEME', ' Extension name', /init, /date, /extver
;;-    sxdelpar, kdr, 'EXTNAME'
;;-    sxaddpar, kdr, 'TTYPE1',   'PLANCK_857', ' Model flux in this band'
;;-    sxaddpar, kdr, 'TTYPE2',   'PLANCK_545', ' Model flux in this band'
;;-    sxaddpar, kdr, 'TTYPE3',   'PLANCK_353', ' Model flux in this band'
;;-    sxaddpar, kdr, 'TTYPE4',   'WISE_12',    ' Model flux in this band'
;;-    sxaddpar, kdr, 'TTYPE5',   'IRAS_60',    ' Model flux in this band'
;;-    sxaddpar, kdr, 'TTYPE6',   'IRAS_100',   ' Model flux in this band'
;;- 
;;-    sxaddpar, kdr, 'COMMENT',  '  ' , after='TFORM6'
;;-    sxaddpar, kdr, 'COMMENT',  ' *** Column units *** '
;;-    sxaddpar, kdr, 'COMMENT',  '  '
;;- 
;;-    sxaddpar, kdr, 'TUNIT1',   'MJy/sr'
;;-    sxaddpar, kdr, 'TUNIT2',   'MJy/sr'
;;-    sxaddpar, kdr, 'TUNIT3',   'MJy/sr' 
;;-    sxaddpar, kdr, 'TUNIT4',   'MJy/sr' 
;;-    sxaddpar, kdr, 'TUNIT5',   'MJy/sr' 
;;-    sxaddpar, kdr, 'TUNIT6',   'MJy/sr' 
;;-  
;;-    sxaddpar, kdr, 'COMMENT',  '  ' , after='TUNIT6'
;;-    sxaddpar, kdr, 'COMMENT',  ' *** Planck params *** '
;;-    sxaddpar, kdr, 'COMMENT',  '  '
;;-    
;;-    sxaddpar, kdr, 'EXTNAME',  'COMP-MAP',    ' Extension name'
;;-    sxaddpar, kdr, 'AST-COMP', 'MODELDUST',   ' Component name'
;;-    sxaddpar, kdr, 'PIXTYPE',  'HEALPIX'
;;-    sxaddpar, kdr, 'COORDSYS', 'GALACTIC',    ' Coordinate system'
;;-    sxaddpar, kdr, 'ORDERING', 'NESTED',      ' Healpix ordering'
;;-    sxaddpar, kdr, 'NSIDE',     Nside,        ' Healpix Nside'
;;-    sxaddpar, kdr, 'FIRSTPIX',  0
;;-    sxaddpar, kdr, 'LASTPIX',   Npix-1
;;-    sxaddpar, kdr, 'BAD_DATA',  -1.63750E+30, ' HEALPIX bad pixel value'
;;-    sxaddpar, kdr, 'FILENAME',  outfile,      ' FITS filename'
;;-    sxaddpar, kdr, 'PROCVER',  procver,       ' Product version'
;;- 
;;-    sxaddpar, kdr, 'COMMENT',  '  ' , after='PROCVER'
;;-    sxaddpar, kdr, 'COMMENT',  '------------------------------------------------------------------------'
;;-    sxaddpar, kdr, 'COMMENT',  'For further details see Planck Explanatory Supplement at:'
;;-    sxaddpar, kdr, 'COMMENT',  '  http://wiki.cosmos.esa.int/planckpla2015'
;;-    sxaddpar, kdr, 'COMMENT',  '------------------------------------------------------------------------'
;;- 
;;-    print, ' >> Exporting '+outdir+outfile
;;-    mwrfits, map, outdir+outfile, kdr, /create   ; Dust opacity map and model
;;-    if verb then print, kdr
;;- 
;;-    if defined(stop) then stop
;;-    return
;;- end 
;;- 
;;- 
;;- ;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;;- ; export lensing map ... updated for PR2
;;- ; NB. gzipping has little effect on size; skip it.
;;- ;-----------------------------------------------------------------------------
;;- 
;;- pro lensing, relnum, outdir=outdir
;;- 
;;-    outdir = '/redtruck/opsman/export/Releases/HFI_PR2_beta/'
;;-    datadir= '/redtruck/dhanson/share/planck_2014_lensing_map'
;;-    spawn, "cd "+datadir+" ; tar cvf "+outdir+"COM_CompMap_Lensing_2048_"+relnum+".tar data"
;;-    spawn, "cd "+datadir+" ; tar cvf "+outdir+"COM_SimMap_Lensing_2048_"+relnum+".tar sims"
;;-    spawn, "cd "+outdir+" ; split -d -b 2048m COM_SimMap_Lensing_2048_"+relnum+".tar COM_SimMap_Lensing_2048_"+relnum+".tar."
;;-    return
;;- end
;;- 
;;- 
;;- ;-----------------------------------------------------------------------------
;;- ; export noise covariance matrices and low-res maps ... updated for PR2
;;- ;-----------------------------------------------------------------------------
;;- 
;;- pro covmat, relnum, outdir=outdir
;;- 
;;-   nside = [8, 16, 32] 
;;-   chans = ['100','143','217','353','545','857']
;;-   outdir = '/redtruck/opsman/export/Releases/HFI_PR2_beta/'
;;- 
;;-   for n = 1,1 do begin   ; nside=16 only used
;;-     ns = nside[n] 
;;-     if ns lt 10 then ss='08' else ss = strtrim(ns,2)
;;-     dir = '/data/rkeskita/dx11_lowres_CPP_scaled/nside'+strtrim(ns,2)+'/'
;;-     for k=0,5 do begin
;;-         ff = chans[k]
;;-         file = outdir+'HFI_NoiseCovMat_'+ff+'_00'+ss+'_'+relnum+'.tgz'
;;-         print,  'cd '+dir+';            ls -d1 maps/HFI_SkyMap_'+ff+'*ns00'+ss+'* ncm/dx11_ncm_'+ff+'*nside00'+ss+'.dat' 
;;-         spawn,  'cd '+dir+'; tar cvzf '+file+' maps/HFI_SkyMap_'+ff+'*ns00'+ss+'* ncm/dx11_ncm_'+ff+'*nside00'+ss+'.dat' 
;;-         print, ' >> Exporting '+file
;;-     endfor 
;;-   endfor 
;;- 
;;-   return
;;- end 
;;- 
;;- ;-----------------------------------------------------------------------------
;;- ; Scanning beams (from Gael)
;;- ;-----------------------------------------------------------------------------
;;- 
;;- pro beam, relnum, procver, det, outdir=outdir, verb=verb
;;- 
;;-    if not defined(outdir)  then outdir  = './'
;;-    if keyword_set(verb)    then verb = 1 else verb = 0
;;-    
;;-    grp = '/data/dmc/MISS03/DATA/BeamSplines_HBM_DX11v67_I5_HIGHRES/BS_HYBRID_SQUARE_RENORM_HR2_'
;;-     name = grp+toBoloID(det)
;;-     pix = (toPixName(det))[0]
;;-     outfile = "HFI_ScanBeam_"+pix+"_"+relnum+".fits"
;;- 
;;-     ima = pioread(name)
;;-     zz = pioreadkeywordobject(nx, com, 'Nx', 'PIODOUBLE', name)
;;-     zz = pioreadkeywordobject(ny, com, 'Ny', 'PIODOUBLE', name)
;;-     zz = pioreadkeywordobject(xx, com, 'Xcentre', 'PIODOUBLE', name)
;;-     zz = pioreadkeywordobject(yy, com, 'Ycentre', 'PIODOUBLE', name)
;;-     zz = pioreadkeywordobject(dx, com, 'Xdelta', 'PIODOUBLE', name)
;;-     zz = pioreadkeywordobject(dy, com, 'Ydelta', 'PIODOUBLE', name)
;;-     ima = reform(ima, nx, ny)/max(ima)
;;- 
;;-     mkhdr, hdr, ima, /IMAGE
;;- ;    sxaddpar, hdr, 'TUNIT1',   '',  ' no unit'
;;-     
;;-     sxaddpar, hdr, 'COMMENT',  '  ' , after='GCOUNT'
;;-     sxaddpar, hdr, 'COMMENT',  ' *** Planck params *** '
;;-     sxaddpar, hdr, 'COMMENT',  '  '
;;-     
;;-     sxaddpar, hdr, 'EXTNAME',  'BEAM',      ' Beam, normalized'
;;-     sxaddpar, hdr, 'BOLO',     pix,         ' Bolometer name'
;;-     sxaddpar, hdr, 'XCTR',     fix(xx),     ' X center'
;;-     sxaddpar, hdr, 'YCTR',     fix(yy),     ' Y center'
;;-     sxaddpar, hdr, 'DELTAX',   dx,          ' X step'    
;;-     sxaddpar, hdr, 'DELTAY',   dy,          ' Y step'    
;;-     sxaddpar, hdr, 'REF',      'DXX',       ' Reference system'    
;;-     
;;-     sxaddpar, hdr, 'FILENAME',  outfile,    ' FITS filename'
;;-     sxaddpar, hdr, 'PROCVER',   procver,    ' Product version'
;;-     
;;-     sxaddpar, hdr, 'COMMENT',  '  ' , after='PROCVER'
;;-     sxaddpar, hdr, 'COMMENT', ''
;;-     sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;-     sxaddpar, hdr, 'COMMENT', ' Scanning beam, normalized to 1 at maximum; no units'
;;-     sxaddpar, hdr, 'COMMENT', '------------------------------------------------------------------------'
;;- 
;;-     print, ' >> Exporting '+outdir+outfile
;;-     writefits, outdir+outfile, ima, hdr ; - beam 
;;-     if verb then print, hdr
;;-     
;;-     return
;;- end


;==================================================================================
;
;   Utilities 
; 
;==================================================================================

;;- ;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;;- ; check headers of EFDD FITS files
;;- ; - tested on RIMO, maps, and some TOI files
;;- ;
;;- ; RIMOs in /redtruck/opsman/export/RIMO  /fits
;;- ; new RIMO from Eric in /data/hivon/RIMO/HFI-RIMO-20120124+DX9beam_v1.fits
;;- ;-----------------------------------------------------------------------------
;;- 
;;- pro hdrcheck, file, verb=verb, quiet=quiet, str=str
;;- 
;;-    on_error, 2
;;- 
;;-    if not defined(file) then begin
;;-        print, ' SYNTAX:  hdrcheck, filename [, /quiet, /verb]'
;;-        return
;;-    endif 
;;- 
;;-    fits_info, file, N_ext=next, extname=extname
;;-    
;;-    if defined(verb) then begin
;;-      for k=0,next do begin
;;-        print, '  '
;;-        print, ';-----------------------------------------------------------------------------'
;;-        print, '; EXTENSION ', strtrim(k,2), ': ', extname[k]
;;-        print, '; - Header'
;;-        print, ';-----------------------------------------------------------------------------'
;;-        ;fits_read, file, data, hdr, exten_no=k, /header_only
;;-        ss = mrdfits(file, k, hdr, range=[0,8])
;;-        print, hdr
;;-        if k ge 1 and keyword_set(str) then begin
;;-            print, ';-----------------------------------------------------------------------------'
;;-            print, '; - Structure'
;;-            print, ';-----------------------------------------------------------------------------'
;;-            showstruct, ss
;;-            print, ''
;;-        end 
;;-      endfor 
;;-    endif 
;;-    print, '  '
;;- 
;;-    return
;;- end
;;- 
;;- 
;;- ;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;;- ; Add checksum parameter to headers of a FITS file and write a new file
;;- ; Will not work for FITS binary tables with variable length arrays though
;;- ; it wouldn't be hard to allow this
;;- ; - from Wayne Landsman, 13.feb.13
;;- ;-----------------------------------------------------------------------------
;;- 
;;- pro addchecksum, file, newfile
;;- 
;;-     if N_params() LT 2 then begin 
;;-         print,'Syntax - addchecksum, file, newfile'
;;-         return
;;-     endif   
;;-     fits_open,file,fcb
;;-     nextend = fcb.nextend
;;- 
;;-     ; Loop over every extension, and write each extension with checksum
;;-     ; keywords to a new file.
;;- 
;;-     for i=0,nextend do begin 
;;-         fits_read, fcb, data, hdr, exten_no=i
;;-         if i eq 0 then writefits, newfile,data, hdr, /checksum else $
;;-           writefits, newfile, data, hdr, /checksum, /append
;;-     endfor
;;- 
;;-     fits_close,fcb 
;;-     return
;;- end                

;==================================================================================
